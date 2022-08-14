MODULE ice_dynamics_module

  ! Contains all the routines needed to calculate ice-sheet geometry at the next time
  ! step, including routines to determine said time step.
  ! NOTE: routines for calculating ice velocities             have been moved to the ice_velocity_module;
  !       routines for integrating the ice thickness equation have been moved to the ice_thickness_module.

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  USE mpi
  USE parameters_module,                     ONLY: ice_density, seawater_density
  USE configuration_module,                  ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE petsc_module,                          ONLY: perr
  USE parallel_module,                       ONLY: par, sync, ierr, cerr, partition_list, &
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
  USE utilities_module,                      ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                                   check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                                   vertical_average, surface_elevation, thickness_above_floatation
  USE netcdf_module,                         ONLY: debug, write_to_debug_file
  USE data_types_module,                     ONLY: type_model_region, type_mesh, type_ice_model, type_reference_geometry, &
                                                   type_remapping_mesh_mesh, type_restart_data
  USE mesh_mapping_module,                   ONLY: remap_field_dp_2D, remap_field_dp_3D
  USE general_ice_model_data_module,         ONLY: update_general_ice_model_data, determine_masks
  USE mesh_operators_module,                 ONLY: map_a_to_c_2D, ddx_a_to_c_2D, ddy_a_to_c_2D
  USE ice_velocity_module,                   ONLY: solve_SIA, solve_SSA, solve_DIVA, initialise_velocity_solver, remap_velocities, &
                                                   map_velocities_b_to_c_2D, map_velocities_b_to_c_3D, calc_3D_vertical_velocities, &
                                                   map_velocities_b_to_a_2D, map_velocities_b_to_a_3D
  USE ice_thickness_module,                  ONLY: calc_dHi_dt
  USE basal_conditions_and_sliding_module,   ONLY: initialise_basal_conditions, remap_basal_conditions
  USE thermodynamics_module,                 ONLY: calc_ice_rheology, remap_ice_temperature
  USE calving_module,                        ONLY: run_calving_model, remove_unconnected_shelves, imposed_shelf_removal
  USE forcing_module,                        ONLY: forcing

  IMPLICIT NONE

CONTAINS

! ===== Run ice dynamics =====
! ============================

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
      CALL run_ice_dynamics_direct( region, t_end)
    ELSEIF (C%choice_timestepping == 'pc') THEN
      CALL run_ice_dynamics_pc( region, t_end)
    ELSE
      CALL crash('unknown choice_timestepping "' // TRIM( C%choice_timestepping) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_model

  SUBROUTINE run_ice_dynamics_direct( region, t_end)
    ! Ice dynamics and time-stepping with the "direct" method
    ! NOTE: this does not work for DIVA ice dynamics!

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ice_dynamics_direct'
    INTEGER                                            :: vi1, vi2
    REAL(dp)                                           :: dt_crit_SIA, dt_crit_SSA
    REAL(dp)                                           :: r_solver_acc, dt_max
    REAL(dp)                                           :: hi_memory, hi_senile, time_passed
    REAL(dp)                                           :: int_old, int_new

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Abbreviations for cleaner code
    vi1 = region%mesh%vi1
    vi2 = region%mesh%vi2

    ! Start-up phase
    ! ==============

    ! Get a more accurate velocity solution during the start-up phase to prevent initialisation "bumps"
    IF (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
      r_solver_acc = 0.01_dp * 0.99_dp * (region%time - C%start_time_of_run) / C%dt_startup_phase
    ELSE
      r_solver_acc = 1._dp
    END IF

    region%ice%DIVA_SOR_nit      = C%DIVA_SOR_nit      * CEILING( 1._dp / r_solver_acc)
    region%ice%DIVA_SOR_tol      = C%DIVA_SOR_tol      * r_solver_acc
    region%ice%DIVA_SOR_omega    = C%DIVA_SOR_omega
    region%ice%DIVA_PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
    region%ice%DIVA_PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc

    ! Reduce the time-step during the start-up phase
    IF     (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
      dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((region%time - C%start_time_of_run) / C%dt_startup_phase)**2
    ELSEIF (region%time >= C%end_time_of_run   - C%dt_startup_phase) THEN
      dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((C%end_time_of_run - region%time  ) / C%dt_startup_phase)**2
    ELSE
      dt_max = C%dt_max
    END IF

    ! Calculate ice velocities with the selected ice-dynamical approximation
    ! ======================================================================

    IF     (C%choice_ice_dynamics == 'none') THEN
      ! Fixed ice geometry

    ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
      ! Shallow ice approximation

      IF (region%time == region%t_next_SIA) THEN

        ! Calculate new ice velocities
        CALL solve_SIA( region%mesh, region%ice)

        ! Calculate critical time step
        CALL calc_critical_timestep_SIA( region%mesh, region%ice, dt_crit_SIA)

        IF (par%master) THEN

          ! Apply conditions to the time step
          dt_crit_SIA = MAX( C%dt_min, MIN( dt_max, dt_crit_SIA))

          ! Update timer
          region%dt_crit_SIA = dt_crit_SIA
          region%t_last_SIA  = region%time
          region%t_next_SIA  = region%time + region%dt_crit_SIA

        END IF
        CALL sync

      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN

    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      ! Shallow shelf approximation

      IF (region%time == region%t_next_SSA) THEN

        ! Calculate new ice velocities
        CALL solve_SSA( region%mesh, region%ice)

        ! Calculate critical time step
        CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_SSA)

        IF (par%master) THEN

          ! Apply conditions to the time step
          dt_crit_SSA = MAX( C%dt_min, MIN( dt_max, dt_crit_SSA))

          ! Update timer
          region%dt_crit_SSA = dt_crit_SSA
          region%t_last_SSA  = region%time
          region%t_next_SSA  = region%time + region%dt_crit_SSA

        END IF
        CALL sync

      END IF ! IF (ABS(region%time - region%t_next_SSA) < dt_tol) THEN

    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      ! Hybrid SIA/SSA (Bueler and Brown, 2009)

      IF (region%time == region%t_next_SIA) THEN

        ! Calculate new ice velocities
        CALL solve_SIA( region%mesh, region%ice)

        ! Calculate critical time step
        CALL calc_critical_timestep_SIA( region%mesh, region%ice, dt_crit_SIA)

        IF (par%master) THEN

          ! Apply conditions to the time step
          dt_crit_SIA = MAX( C%dt_min, MIN( dt_max, dt_crit_SIA))

          ! Update timer
          region%dt_crit_SIA = dt_crit_SIA
          region%t_last_SIA  = region%time
          region%t_next_SIA  = region%time + region%dt_crit_SIA

        END IF
        CALL sync

      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN

      IF (region%time == region%t_next_SSA) THEN

        ! Calculate new ice velocities
        CALL solve_SSA( region%mesh, region%ice)

        ! Calculate critical time step
        CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_SSA)

        IF (par%master) THEN

          ! Apply conditions to the time step
          dt_crit_SSA = MAX( C%dt_min, MIN( dt_max, dt_crit_SSA))

          ! Update timer
          region%dt_crit_SSA = dt_crit_SSA
          region%t_last_SSA  = region%time
          region%t_next_SSA  = region%time + region%dt_crit_SSA

        END IF
        CALL sync

      END IF ! IF (ABS(region%time - region%t_next_SSA) < dt_tol) THEN

    ELSE ! IF     (C%choice_ice_dynamics == 'SIA') THEN
      CALL crash('"direct" time stepping works only with SIA, SSA, or SIA/SSA ice dynamics, not with DIVA!')
    END IF ! IF     (C%choice_ice_dynamics == 'SIA') THEN

    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps_and_actions( region, t_end)

    ! Calculate new ice geometry
    ! ==========================

    IF (C%choice_ice_dynamics == 'none') THEN

      ! Fixed ice geometry
      region%ice%dHi_dt_a( vi1:vi2) = 0._dp
      region%ice%Hi_tplusdt_a( vi1:vi2) = region%ice%Hi_a( vi1:vi2)
      CALL sync

    ELSE

      ! Compute new dHi_dt and corresponding Hi_tplusdt_a
      CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt, region%mask_noice, region%refgeo_PD)

      ! Calculate ice thickness at the end of this model loop
      region%ice%Hi_tplusdt_a( vi1:vi2) = MAX( 0._dp, region%ice%Hi_a( vi1:vi2) + region%dt * region%ice%dHi_dt_a( vi1:vi2))
      CALL sync

    END IF

    ! Running dHi_dt average
    ! ======================

    ! Update the running window: drop oldest record and push the rest to the back
    region%ice%dHi_dt_window_a(vi1:vi2,2:C%dHi_dt_window_size) = region%ice%dHi_dt_window_a(vi1:vi2,1:C%dHi_dt_window_size-1)

    ! Update the running window: add new record to beginning of window
    region%ice%dHi_dt_window_a(vi1:vi2,1) = region%ice%dHi_dt_a(vi1:vi2)

    ! Compute running average
    region%ice%dHi_dt_ave_a( vi1:vi2) = SUM(region%ice%dHi_dt_window_a(vi1:vi2,:),2) / REAL(C%dHi_dt_window_size,dp)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_direct

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
    REAL(dp)                                           :: r_solver_acc, dt_max
    REAL(dp)                                           :: hi_memory, hi_senile
    REAL(dp)                                           :: time_passed, int_old, int_new!, Hi_tot_new

    ! == Initialisation
    ! =================

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Abbreviations for cleaner code
    vi1 = region%mesh%vi1
    vi2 = region%mesh%vi2

    ! Determine whether or not we need to update ice velocities
    do_update_ice_velocity = .FALSE.
    IF     (C%choice_ice_dynamics == 'none') THEN
      region%ice%dHi_dt_a(     vi1:vi2) = 0._dp
      region%ice%Hi_corr(      vi1:vi2) = region%ice%Hi_a( vi1:vi2)
      region%ice%Hi_tplusdt_a( vi1:vi2) = region%ice%Hi_a( vi1:vi2)
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
    CALL sync

    ! == Start-up / Cool-down
    ! =======================

    ! Default velocity accuracy
    r_solver_acc = 1._dp

    ! Default max time step
    dt_max = C%dt_max

    ! Start-up
    IF (C%dt_startup_phase > 0._dp) THEN
      IF (region%time < C%start_time_of_run) THEN
        ! Wind-up
        dt_max = C%dt_min
        r_solver_acc = 0.01_dp
      ELSEIF (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
        ! Reduce time maximum time step
        dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((region%time - C%start_time_of_run) / C%dt_startup_phase)**2
        ! Enforce higher accuracy
        r_solver_acc = 0.01_dp + 0.99_dp * (region%time - C%start_time_of_run) / C%dt_startup_phase
      END IF
    END IF

    ! Cool-down
    IF (C%dt_cooldown_phase > 0._dp) THEN
      IF (region%time >= C%end_time_of_run - C%dt_cooldown_phase) THEN
        ! Reduce time maximum time step
        dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((C%end_time_of_run - region%time) / C%dt_cooldown_phase)**2
        ! Enforce higher accuracy
        r_solver_acc = 0.01_dp + 0.99_dp * (C%end_time_of_run - region%time) / C%dt_cooldown_phase
      END IF
    END IF

    ! Accuracy for the SSA/DIVA solution
    region%ice%DIVA_SOR_nit      = C%DIVA_SOR_nit      * CEILING( 1._dp / r_solver_acc)
    region%ice%DIVA_SOR_tol      = C%DIVA_SOR_tol      * r_solver_acc
    region%ice%DIVA_SOR_omega    = C%DIVA_SOR_omega
    region%ice%DIVA_PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
    region%ice%DIVA_PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc

    ! == Velocity update
    ! ==================

    IF (do_update_ice_velocity) THEN

      ! Calculate time step based on the truncation error in ice thickness (Robinson et al., 2020, Eq. 33)
      CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_adv)

      IF (par%master) THEN

        ! Calculate critical time step
        region%dt_crit_ice_prev = region%dt_crit_ice
        dt_from_pc              = (C%pc_epsilon / region%ice%pc_eta)**(C%pc_k_I + C%pc_k_p) * (C%pc_epsilon / region%ice%pc_eta_prev)**(-C%pc_k_p) * region%dt
        region%dt_crit_ice      = MAX( C%dt_min, MAX( 0.5_dp * region%dt_crit_ice_prev, MINVAL([ C%dt_max, 2._dp * region%dt_crit_ice_prev, dt_crit_adv, dt_from_pc])))

        ! Apply conditions to the time step
        region%dt_crit_ice = MAX( C%dt_min, MIN( dt_max, region%dt_crit_ice))

        ! Calculate zeta
        region%ice%pc_zeta      = region%dt_crit_ice / region%dt_crit_ice_prev

      END IF
      CALL sync

      ! Predictor step
      ! ==============

      ! Calculate new ice geometry
      region%ice%pc_f2(   vi1:vi2) = region%ice%dHi_dt_a( vi1:vi2)

      CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice, region%refgeo_PD)
      region%ice%pc_f1(   vi1:vi2) = region%ice%dHi_dt_a( vi1:vi2)

      ! Robinson et al. (2020), Eq. 30)
      region%ice%Hi_pred( vi1:vi2) = MAX(0._dp, region%ice%Hi_a( vi1:vi2)   + region%dt_crit_ice * &
                                      ((1._dp + region%ice%pc_zeta / 2._dp) * region%ice%pc_f1( vi1:vi2) - &
                                               (region%ice%pc_zeta / 2._dp) * region%ice%pc_f2( vi1:vi2)))
      CALL sync

      ! Update step
      ! ===========

      ! Calculate velocities for predicted geometry
      region%ice%Hi_old( vi1:vi2) = region%ice%Hi_a(    vi1:vi2)
      region%ice%Hi_a(   vi1:vi2) = region%ice%Hi_pred( vi1:vi2)
      CALL update_general_ice_model_data( region%mesh, region%ice)

      IF     (C%choice_ice_dynamics == 'SIA') THEN

        ! Calculate velocities
        CALL solve_SIA(  region%mesh, region%ice)

        ! Update timer
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_ice
        CALL sync

      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN

        ! Calculate velocities
        CALL solve_SSA(  region%mesh, region%ice)

        ! Update timer
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_ice
        CALL sync

      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN

        ! Calculate velocities
        CALL solve_SIA(  region%mesh, region%ice)
        CALL solve_SSA(  region%mesh, region%ice)

        ! Update timer
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_ice
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_ice
        CALL sync

      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN

        ! Calculate velocities
        CALL solve_DIVA( region%mesh, region%ice)

        ! Update timer
        IF (par%master) region%t_last_DIVA = region%time
        IF (par%master) region%t_next_DIVA = region%time + region%dt_crit_ice
        CALL sync

      ELSE
        CALL crash('unknown choice_ice_dynamics "' // TRIM( C%choice_ice_dynamics) // '"!')
      END IF

      ! Corrector step
      ! ==============

      ! Calculate dHi_dt for the predicted ice thickness and updated velocity
      CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice, region%refgeo_PD)
      region%ice%pc_f3(   vi1:vi2) = region%ice%dHi_dt_a( vi1:vi2)
      region%ice%pc_f4(   vi1:vi2) = region%ice%pc_f1(    vi1:vi2)

      ! Go back to old ice thickness. Run all the other modules (climate, SMB, BMB, thermodynamics, etc.)
      ! and only go to new (corrected) ice thickness at the end of this time loop.
      region%ice%Hi_a(    vi1:vi2) = region%ice%Hi_old(   vi1:vi2)

      CALL sync

      CALL update_general_ice_model_data( region%mesh, region%ice)

      ! Calculate "corrected" ice thickness (Robinson et al. (2020), Eq. 31)
      region%ice%Hi_corr( vi1:vi2) = MAX(0._dp, region%ice%Hi_old( vi1:vi2) + 0.5_dp * region%dt_crit_ice * &
                                               (region%ice%pc_f4( vi1:vi2) + region%ice%pc_f3( vi1:vi2)))

      ! Calculate applied ice thickness rate of change
      region%ice%dHi_dt_a( vi1:vi2) = (region%ice%Hi_corr( vi1:vi2) - region%ice%Hi_a( vi1:vi2)) / region%dt_crit_ice

      CALL sync

      ! Truncation error
      ! ================

      ! Determine truncation error
      CALL calc_pc_truncation_error( region%mesh, region%ice, region%dt_crit_ice)

      ! Memory step
      ! ===========

      ! Uses the values of dHi_dt from the start of the run (or values from
      ! a previous run, in case of a restart) and computes a weighed average
      ! between the old and new values. The average starts giving full weight
      ! to the old dHi_dt (start of run), and exponentially (but slowly) shifts
      ! the weight to the new values (current time step) over time. At time
      ! C%start_time_of_run + C%get_senile_after, memory is completely lost
      ! and the newly computed dHi_dt is fully used from that time onwards.

      ! Compute how much time has passed since start of run
      time_passed = region%time - C%start_time_of_run

      ! Only if so specified and before memory fades away
      IF (C%do_use_hi_memory .AND. time_passed < C%get_senile_after) THEN

        ! Compute how senile our current run is
        IF (region%time < C%start_time_of_run) THEN
          ! For wind-up, keep memory intact
          hi_senile = 0._dp
        ELSE
          ! For actual runs, get senility level: 0 at start, 1 after C%get_senile_after years
          hi_senile = min(1._dp, max(0._dp, time_passed / C%get_senile_after))
        END IF

        ! Memory fades
        hi_memory = 1._dp - hi_senile
        ! No matter what
        hi_memory = min(hi_memory, 0.9_dp)

        IF (C%is_restart) THEN
          ! For restarts, weighed average between previous life and present
          region%ice%dHi_dt_a( vi1:vi2) = (1._dp - hi_memory) * region%ice%dHi_dt_a( vi1:vi2) &
                                                 + hi_memory  * region%ice%dHi_dt_past_a( vi1:vi2)
        ELSE
          ! Else, weighed average between birthday and present
          region%ice%dHi_dt_a( vi1:vi2) = (1._dp - hi_memory) * region%ice%dHi_dt_a( vi1:vi2) &
                                                 + hi_memory  * 0._dp
        END IF
        CALL sync

        ! Even when your neurons are no more, their atoms
        ! carry on, so make sure mass is conserved

        ! Total memory-based ice change rate
        int_new = SUM(region%ice%dHi_dt_a( vi1:vi2))
        ! Total original ice change rate
        int_old = SUM( (region%ice%Hi_corr( vi1:vi2) - region%ice%Hi_a( vi1:vi2)) / region%dt_crit_ice )
        CALL sync

        CALL MPI_ALLREDUCE( MPI_IN_PLACE, int_new, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, int_old, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

        IF (int_new > 0._dp) THEN
          ! Scale memory-based rates to conserve total mass change
          region%ice%dHi_dt_a( vi1:vi2) = region%ice%dHi_dt_a( vi1:vi2) * (int_old / int_new)
        END IF
        CALL sync

      END IF ! (C%do_use_hi_memory)

      ! Running dHi_dt average
      ! ======================

      ! Update the running window: drop oldest record and push the rest to the back
      region%ice%dHi_dt_window_a(vi1:vi2,2:C%dHi_dt_window_size) = region%ice%dHi_dt_window_a(vi1:vi2,1:C%dHi_dt_window_size-1)

      ! Update the running window: add new record to beginning of window
      region%ice%dHi_dt_window_a(vi1:vi2,1) = region%ice%dHi_dt_a(vi1:vi2)

      ! Compute running average
      region%ice%dHi_dt_ave_a( vi1:vi2) = SUM(region%ice%dHi_dt_window_a(vi1:vi2,:),2) / REAL(C%dHi_dt_window_size,dp)

    END IF ! (do_update_ice_velocity)

    ! Final quantities
    ! ================

    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps_and_actions( region, t_end)

    ! Calculate ice thickness at the end of the adjusted time step
    region%ice%Hi_tplusdt_a( vi1:vi2) = MAX( 0._dp, region%ice%Hi_a( vi1:vi2) + region%dt * region%ice%dHi_dt_a( vi1:vi2))
    CALL sync

    ! ! Apply anti-shock after mesh update
    ! IF ( (region%ice%Hi_tot_old_time /= C%start_time_of_run) .AND. &
    !      (region%time - region%ice%Hi_tot_old_time < 100._dp) ) THEN

    !   ! Get new total ice thickness
    !   Hi_tot_new = SUM(region%ice%Hi_tplusdt_a( vi1:vi2))
    !   CALL sync
    !   CALL MPI_ALLREDUCE( MPI_IN_PLACE, Hi_tot_new, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    !   IF (Hi_tot_new > 0._dp) THEN
    !     ! Scale new ice thickness to conserve total mass from before the mesh-update
    !     region%ice%Hi_tplusdt_a( vi1:vi2) = region%ice%Hi_tplusdt_a( vi1:vi2) * (region%ice%Hi_tot_old / Hi_tot_new)
    !   END IF
    !   CALL sync

    ! END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ice_dynamics_pc

! ===== Ice thickness update =====
! ================================

  SUBROUTINE update_ice_thickness( mesh, ice, mask_noice, refgeo_PD, refgeo_GIAeq)
    ! Update the ice thickness at the end of a model time loop

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                             INTENT(IN)    :: mesh
    TYPE(type_ice_model),                        INTENT(INOUT) :: ice
    INTEGER,                       DIMENSION(:), INTENT(IN)    :: mask_noice
    TYPE(type_reference_geometry),               INTENT(IN)    :: refgeo_PD
    TYPE(type_reference_geometry),               INTENT(IN)    :: refgeo_GIAeq

    ! Local variables:
    CHARACTER(LEN=256),    PARAMETER                           :: routine_name = 'update_ice_thickness'
    INTEGER                                                    :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Save the previous ice mask, for use in thermodynamics
    ice%mask_ice_a_prev( mesh%vi1:mesh%vi2) = ice%mask_ice_a( mesh%vi1:mesh%vi2)
    CALL sync

    ! Check if we want to keep the ice thickness fixed for some specific areas
    CALL fix_ice_thickness( mesh, ice)

    ! Set ice thickness to new value
    ice%Hi_a( mesh%vi1:mesh%vi2) = MAX( 0._dp, ice%Hi_tplusdt_a( mesh%vi1:mesh%vi2))
    CALL sync

    ! Apply calving
    CALL run_calving_model( mesh, ice)

    ! Remove ice following various criteria
    CALL additional_ice_removal( mesh, ice, mask_noice, refgeo_PD, refgeo_GIAeq)

    ! Update masks, surface elevation, and thickness above floatation
    CALL update_general_ice_model_data( mesh, ice)

    ! Compute ice thickness/elevation difference w.r.t PD
    ice%dHi_a( mesh%vi1:mesh%vi2) = ice%Hi_a( mesh%vi1:mesh%vi2) - refgeo_PD%Hi( mesh%vi1:mesh%vi2)
    ice%dHs_a( mesh%vi1:mesh%vi2) = ice%Hs_a( mesh%vi1:mesh%vi2) - refgeo_PD%Hs( mesh%vi1:mesh%vi2)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_ice_thickness

  SUBROUTINE fix_ice_thickness( mesh, ice)
    ! Check if we want to keep the ice thickness fixed in time for specific areas

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256),    PARAMETER                   :: routine_name = 'fix_ice_thickness'
    INTEGER                                            :: vi, ci, vc
    LOGICAL                                            :: fix_Hi_here


    ! Add routine to path
    CALL init_routine( routine_name)

    ! If so specified, keep sheet geometry fixed (incl. first adjacent shelf points)
    IF (C%fixed_sheet_geometry) THEN
      DO vi = mesh%vi1, mesh%vi2

        fix_Hi_here = .FALSE.

        IF (ice%mask_sheet_a( vi) == 1) THEN
          fix_Hi_here = .TRUE.
        ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            IF (ice%mask_sheet_a( vc) == 1) THEN
              fix_Hi_here = .TRUE.
              EXIT
            END IF
          END DO
        END IF

        IF (fix_Hi_here) THEN
          ice%Hi_tplusdt_a( vi) = ice%Hi_a( vi)
        END IF

      END DO
      CALL sync
    END IF

    ! If so specified, keep GL position fixed (both grounded and floating sides)
    IF (C%fixed_grounding_line) THEN
      DO vi = mesh%vi1, mesh%vi2

        fix_Hi_here = .FALSE.

        IF (ice%mask_sheet_a( vi) == 1) THEN
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            IF (ice%mask_shelf_a( vc) == 1) THEN
              fix_Hi_here = .TRUE.
              EXIT
            END IF
          END DO
        ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            IF (ice%mask_sheet_a( vc) == 1) THEN
              fix_Hi_here = .TRUE.
              EXIT
            END IF
          END DO
        END IF

        IF (fix_Hi_here) THEN
          ice%Hi_tplusdt_a( vi) = ice%Hi_a( vi)
        END IF

      END DO
      CALL sync
    END IF

    ! If so specified, keep shelf geometry fixed (incl. first adjacent sheet points)
    IF (C%fixed_shelf_geometry) THEN
      DO vi = mesh%vi1, mesh%vi2

        fix_Hi_here = .FALSE.

        IF (ice%mask_shelf_a( vi) == 1) THEN
          fix_Hi_here = .TRUE.
        ELSEIF (ice%mask_sheet_a( vi) == 1) THEN
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            IF (ice%mask_shelf_a( vc) == 1) THEN
              fix_Hi_here = .TRUE.
              EXIT
            END IF
          END DO
        END IF

        IF (fix_Hi_here) THEN
          ice%Hi_tplusdt_a( vi) = ice%Hi_a( vi)
        END IF

      END DO
      CALL sync
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE fix_ice_thickness

  SUBROUTINE additional_ice_removal( mesh, ice, mask_noice, refgeo_PD, refgeo_GIAeq)
    ! Remove ice following various criteria

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                             INTENT(IN)    :: mesh
    TYPE(type_ice_model),                        INTENT(INOUT) :: ice
    INTEGER,                       DIMENSION(:), INTENT(IN)    :: mask_noice
    TYPE(type_reference_geometry),               INTENT(IN)    :: refgeo_PD
    TYPE(type_reference_geometry),               INTENT(IN)    :: refgeo_GIAeq

    ! Local variables:
    CHARACTER(LEN=256),    PARAMETER                           :: routine_name = 'additional_ice_removal'
    INTEGER                                                    :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! ! Remove very thin ice
    ! DO vi = mesh%vi1, mesh%vi2
    !   IF (ice%Hi_a( vi) < 1.0_dp) THEN
    !     ice%Hi_a( vi) = 0._dp
    !   END IF
    ! END DO
    ! CALL sync

    ! Remove ice in areas where no ice is allowed (e.g. Greenland in NAM and EAS, and Ellesmere Island in GRL)
    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice( vi) == 1) THEN
        ice%Hi_a( vi) = 0._dp
      END IF
    END DO
    CALL sync

    ! If so specified, remove specific shelf areas
    CALL imposed_shelf_removal( mesh, ice, refgeo_PD, refgeo_GIAeq)

    ! Remove unconnected shelves
    CALL remove_unconnected_shelves( mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE additional_ice_removal

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

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nAc,       u_c     , wu_c     )
    CALL allocate_shared_dp_1D( mesh%nAc,       v_c     , wv_c     )
    CALL allocate_shared_dp_1D( mesh%nAc,       Hi_c    , wHi_c    )
    CALL allocate_shared_dp_1D( mesh%nAc,       dHs_dx_c, wdHs_dx_c)
    CALL allocate_shared_dp_1D( mesh%nAc,       dHs_dy_c, wdHs_dy_c)
    CALL allocate_shared_dp_2D( mesh%nAc, C%nz, u_3D_c  , wu_3D_c  )
    CALL allocate_shared_dp_2D( mesh%nAc, C%nz, v_3D_c  , wv_3D_c  )

    ! Calculate ice velocity and thickness, and surface slopes on the staggered grid
    CALL map_velocities_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, u_c, v_c)
    CALL map_a_to_c_2D( mesh, ice%Hi_a, Hi_c    )
    CALL ddx_a_to_c_2D( mesh, ice%Hs_a, dHs_dx_c)
    CALL ddy_a_to_c_2D( mesh, ice%Hs_a, dHs_dy_c)
    CALL map_velocities_b_to_c_3D( mesh, ice%u_3D_b, ice%v_3D_b, u_3D_c, v_3D_c)

    ! Initialise time step with maximum allowed value
    dt_crit_SIA = C%dt_max

    DO aci = mesh%ci1, mesh%ci2

      ! Calculate the SIA ice diffusivity
      D_SIA = 1E-9_dp
      D_SIA = MAX( D_SIA, ABS( u_c( aci) * Hi_c( aci) / dHs_dx_c( aci) ))
      D_SIA = MAX( D_SIA, ABS( v_c( aci) * Hi_c( aci) / dHs_dy_c( aci) ))

      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      dist = NORM2( mesh%V( vi,:) - mesh%V( vj,:))
      dt = dist**2 / (6._dp * D_SIA)
      dt_crit_SIA = MIN(dt_crit_SIA, dt)

      ! Also check the 3D advective time step
      DO k = 1, C%nz
        dt = dist / (ABS( u_3D_c( aci,k)) + ABS( v_3D_c( aci,k)))
        dt_crit_SIA = MIN( dt_crit_SIA, dt)
      END DO

    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_SIA, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_SIA = dt_crit_SIA * dt_correction_factor

    ! Safety
    dt_crit_SIA = MAX(dt_crit_SIA, C%dt_min)

    ! Clean up after yourself
    CALL deallocate_shared( wu_c     )
    CALL deallocate_shared( wv_c     )
    CALL deallocate_shared( wHi_c    )
    CALL deallocate_shared( wdHs_dx_c)
    CALL deallocate_shared( wdHs_dy_c)
    CALL deallocate_shared( wu_3D_c  )
    CALL deallocate_shared( wv_3D_c  )

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
    REAL(dp), DIMENSION(:    ), POINTER                ::  u_c,  v_c
    INTEGER                                            :: wu_c, wv_c
    REAL(dp)                                           :: dist, dt, uv_c
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nAc, u_c, wu_c)
    CALL allocate_shared_dp_1D( mesh%nAc, v_c, wv_c)

    ! Calculate ice velocity on the staggered grid
    CALL map_velocities_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, u_c, v_c)

    dt_crit_adv = 2._dp * C%dt_max

    DO aci = mesh%ci1, mesh%ci2

      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      dist = NORM2( mesh%V( vi,:) - mesh%V( vj,:))
      uv_c = (ABS( u_c( aci)) + ABS( v_c( aci)))
      IF (uv_c > 0._dp) THEN
        dt = dist / uv_c
      ELSE
        dt = 2._dp * C%dt_max
      END IF
      dt_crit_adv = MIN( dt_crit_adv, dt)

    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = MIN(C%dt_max, dt_crit_adv * dt_correction_factor)

    ! Safety
    dt_crit_adv = MAX(dt_crit_adv, C%dt_min)

    ! Clean up after yourself
    CALL deallocate_shared( wu_c)
    CALL deallocate_shared( wv_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_adv

  SUBROUTINE calc_pc_truncation_error( mesh, ice, dt)
    ! Calculate the truncation error in the ice thickness rate of change (Robinson et al., 2020, Eq. 32)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pc_truncation_error'
    INTEGER                                            :: vi, ci, vc
    LOGICAL                                            :: has_GL_neighbour
    REAL(dp)                                           :: eta_proc

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Save previous value of eta
    ice%pc_eta_prev  = ice%pc_eta

    ! Find maximum truncation error
    eta_proc = C%pc_eta_min

    DO vi = mesh%vi1, mesh%vi2

      ! Calculate truncation error (Robinson et al., 2020, Eq. 32)
      ice%pc_tau( vi) = ABS( ice%pc_zeta * (ice%Hi_corr( vi) - ice%Hi_pred( vi)) / &
                             ((3._dp * ice%pc_zeta + 3._dp) * dt))

      IF (ice%mask_sheet_a( vi) == 1) THEN

        has_GL_neighbour = .FALSE.
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_gl_a( vc) == 1) THEN
            has_GL_neighbour = .TRUE.
            EXIT
          END IF
        END DO

        IF (.NOT.has_GL_neighbour) THEN
          eta_proc = MAX( eta_proc, ice%pc_tau( vi))
        END IF

      END IF

    END DO
    CALL MPI_REDUCE( eta_proc, ice%pc_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pc_truncation_error

  SUBROUTINE determine_timesteps_and_actions( region, t_end)
    ! Determine how long we can run just ice dynamics before another "action" (thermodynamics,
    ! GIA, output writing, inverse routine, etc.) has to be performed.

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    REAL(dp),                            INTENT(IN)    :: t_end

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'determine_timesteps_and_actions'
    REAL(dp)                                           :: t_next

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN
      ! Determine when each model components should be updated

      ! Default time step
      t_next = MIN(t_end, region%time + C%dt_max)

      ! First the ice dynamics
      ! ======================

      IF     (C%choice_ice_dynamics == 'none') THEN
        ! Just stick to the maximum time step

      ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
        ! Use SIA time step
        t_next = MIN( t_next, region%t_next_SIA)

      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
        ! Use SSA time step
        t_next = MIN( t_next, region%t_next_SSA)

      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
        ! Use SIA/SSA time step
        t_next = MIN( t_next, region%t_next_SIA)
        t_next = MIN( t_next, region%t_next_SSA)

      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
        ! Use DIVA time step
        t_next = MIN( t_next, region%t_next_DIVA)

      ELSE
        CALL crash('unknown choice_ice_dynamics "' // TRIM( C%choice_ice_dynamics) // '"!')
      END IF

      ! Then the other model components
      ! ===============================

      region%do_thermo  = .FALSE.
      IF (region%time >= region%t_next_thermo) THEN
        region%do_thermo      = .TRUE.
        region%t_last_thermo  = region%time
        region%t_next_thermo  = region%t_last_thermo + C%dt_thermo
      END IF

      region%do_climate = .FALSE.
      IF (region%time >= region%t_next_climate) THEN
        region%do_climate     = .TRUE.
        region%t_last_climate = region%time
        region%t_next_climate = region%t_last_climate + C%dt_climate
      END IF

      region%do_ocean   = .FALSE.
      IF (region%time >= region%t_next_ocean) THEN
        region%do_ocean       = .TRUE.
        region%t_last_ocean   = region%time
        region%t_next_ocean   = region%t_last_ocean + C%dt_ocean
      END IF

      region%do_SMB     = .FALSE.
      IF (region%time >= region%t_next_SMB) THEN
        region%do_SMB         = .TRUE.
        region%t_last_SMB     = region%time
        region%t_next_SMB     = region%t_last_SMB + C%dt_SMB
      END IF

      region%do_BMB     = .FALSE.
      IF (region%time >= region%t_next_BMB) THEN
        region%do_BMB         = .TRUE.
        region%t_last_BMB     = region%time
        region%t_next_BMB     = region%t_last_BMB + C%dt_BMB
      END IF

      region%do_ELRA    = .FALSE.
      IF (C%choice_GIA_model == 'ELRA') THEN
        IF (region%time >= region%t_next_ELRA) THEN
          region%do_ELRA      = .TRUE.
          region%t_last_ELRA  = region%time
          region%t_next_ELRA  = region%t_last_ELRA + C%dt_bedrock_ELRA
        END IF
      END IF

      region%do_basal    = .FALSE.
      IF (C%do_basal_sliding_inversion) THEN
        IF (region%time >= region%t_next_basal) THEN
          region%do_basal     = .TRUE.
          region%t_last_basal = region%time
          region%t_next_basal = region%t_last_basal + C%dt_basal
        END IF
      END IF

      region%do_SMB_inv = .FALSE.
      IF (C%do_SMB_IMAUITM_inversion) THEN
        IF (region%time >= region%t_next_SMB_inv) THEN
          region%do_SMB_inv     = .TRUE.
          region%t_last_SMB_inv = region%time
          region%t_next_SMB_inv = region%t_last_SMB_inv + C%dt_SMB_inv
        END IF
      END IF

      ! Finally, the output
      ! ===================

      ! This requires a couple of additional considerations,
      ! to avoid output records at weird, decimal-laden
      ! times. We use the fact that nothing in the model
      ! really changes between one time step and the next
      ! one, i.e., the model has the same state at all times
      ! between those steps. Thus, we identify the two time
      ! steps encompasing our desired, rounded output time,
      ! and then simply ask the model to output that model
      ! state using the desired time as time-of-record.

      region%do_output  = .FALSE.
      ! Check if we will overshoot the output timestep
      IF (t_next > region%t_next_output) THEN
        ! Output before ice and time are updated
        region%do_output      = .TRUE.
        ! The output routines will record this time
        region%t_last_output  = region%t_next_output
        ! Set next output at desired time
        region%t_next_output  = region%t_last_output + C%dt_output
      END IF

      ! Set time step so that we move forward to the next action
      region%dt = t_next - region%time

    END IF ! (par%master)
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_timesteps_and_actions

! ===== Administration: allocation, initialisation, and remapping =====
! =====================================================================

  SUBROUTINE initialise_ice_model( mesh, ice, refgeo_init, refgeo_PD, restart)
    ! Allocate shared memory for all the data fields of the ice dynamical module, and
    ! initialise some of them

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    TYPE(type_restart_data),             INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_model'
    INTEGER                                            :: vi, ti, k
    REAL(dp), DIMENSION(C%nz)                          :: prof

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  Initialising ice dynamics model...'

    ! Allocate shared memory
    CALL allocate_ice_model( mesh, ice)

    IF (C%is_restart) THEN
      ! Get sea level from restart data
      ice%SL_a( mesh%vi1:mesh%vi2) = restart%SL( mesh%vi1:mesh%vi2)

    ELSE
      ! Get sea level at initial time
      IF (C%choice_sealevel_model == 'fixed') THEN
        ice%SL_a( mesh%vi1:mesh%vi2) = C%fixed_sealevel

      ELSEIF (C%choice_sealevel_model == 'prescribed') THEN
        ice%SL_a( mesh%vi1:mesh%vi2) = forcing%sealevel_obs

      ELSEIF (C%choice_sealevel_model == 'eustatic') THEN
        ice%SL_a( mesh%vi1:mesh%vi2) = C%initial_guess_sealevel

      ELSEIF (C%choice_sealevel_model == 'SELEN') THEN
       ice%SL_a( mesh%vi1:mesh%vi2) = C%initial_guess_sealevel

      ELSE
        CALL crash('unknown choice_sealevel_model "' // TRIM( C%choice_sealevel_model) // '"!')
      END IF

    END IF

    ! Initialise with data from initial file (works for both "raw" initial data and model restarts)
    DO vi = mesh%vi1, mesh%vi2
      ! Main quantities
      ice%Hi_a( vi)   = refgeo_init%Hi( vi)
      ice%Hb_a( vi)   = refgeo_init%Hb( vi)

      IF (C%is_restart) THEN
        ice%Hs_a( vi) = refgeo_init%Hs( vi)
      ELSE
        ice%Hs_a( vi) = surface_elevation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))
      END IF

      ice%TAF_a( vi)  = thickness_above_floatation( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))

      ! Differences w.r.t. present-day
      ice%dHi_a( vi)  = ice%Hi_a( vi) - refgeo_PD%Hi( vi)
      ice%dHb_a( vi)  = ice%Hb_a( vi) - refgeo_PD%Hb( vi)
      ice%dHs_a( vi)  = ice%Hs_a( vi) - refgeo_PD%Hs( vi)
    END DO
    CALL sync

    ! Determine masks
    CALL determine_masks( mesh, ice)

    IF (C%is_restart) THEN

      ice%dHb_dt_a( mesh%vi1:mesh%vi2) = restart%dHb_dt( mesh%vi1:mesh%vi2)
      ice%dHi_dt_a( mesh%vi1:mesh%vi2) = restart%dHi_dt( mesh%vi1:mesh%vi2)

      DO vi = mesh%vi1, mesh%vi2
        IF (ice%mask_land_a( vi) == 1) THEN
          ice%dHs_dt_a( vi) = ice%dHb_dt_a( vi) + ice%dHi_dt_a( vi)
        ELSE
          ice%dHs_dt_a( vi) = ice%dHi_dt_a( vi) * (1._dp - ice_density / seawater_density)
        END IF
      END DO

      ! Averaged dHi_dt from restart file
      ice%dHi_dt_past_a( mesh%vi1:mesh%vi2) = restart%dHi_dt_ave( mesh%vi1:mesh%vi2)

      ! Running window for averaged dHi_dt for current run
      DO k = 1, C%dHi_dt_window_size
        ice%dHi_dt_window_a( mesh%vi1:mesh%vi2,k) = ice%dHi_dt_past_a( mesh%vi1:mesh%vi2)
      END DO

      ! Averaged dHi_dt for current run
      ice%dHi_dt_ave_a( mesh%vi1:mesh%vi2) = ice%dHi_dt_past_a( mesh%vi1:mesh%vi2)

    ELSE

      ice%dHb_dt_a = 0._dp
      ice%dHi_dt_a = 0._dp
      ice%dHs_dt_a = 0._dp

      ice%dHi_dt_past_a   = 0._dp
      ice%dHi_dt_window_a = 0._dp
      ice%dHi_dt_ave_a    = 0._dp

    END IF

    ! Initialise the "previous ice mask", so that the first call to thermodynamics works correctly
    ice%mask_ice_a_prev( mesh%vi1:mesh%vi2) = ice%mask_ice_a( mesh%vi1:mesh%vi2)
    CALL sync

    ! Initialise some numbers for the predictor/corrector ice thickness update method
    IF (par%master) THEN
      ice%pc_zeta        = 1._dp
      ice%pc_eta         = C%pc_epsilon
      ice%pc_eta_prev    = C%pc_epsilon
    END IF
    CALL sync

    ! ! Initialise total ice thickness at start of run
    ! ice%Hi_tot_old = SUM( ice%Hi_a( mesh%vi1:mesh%vi2))
    ! CALL sync
    ! ! Communicate mid-point results among processes
    ! CALL MPI_ALLREDUCE( MPI_IN_PLACE, ice%Hi_tot_old, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    ! ! Initialise time (to make sure it is not used until a mesh update)
    ! ice%Hi_tot_old_time = C%start_time_of_run

    ! Allocate and initialise basal conditions
    CALL initialise_basal_conditions( mesh, ice, restart)

    ! Geothermal heat flux
    IF     (C%choice_geothermal_heat_flux == 'constant') THEN
      ice%GHF_a( mesh%vi1:mesh%vi2) = C%constant_geothermal_heat_flux
    ELSEIF (C%choice_geothermal_heat_flux == 'spatial') THEN
      CALL map_geothermal_heat_flux_to_mesh( mesh, ice)
    ELSE
      CALL crash('unknown choice_geothermal_heat_flux "' // TRIM( C%choice_geothermal_heat_flux) // '"!')
    END IF

    ! Initialise data and matrices for the velocity solver(s)
    ! If we're running with choice_ice_dynamics == "none", initialise the velocity
    ! solver pretending we are using the DIVA, as we will need to get reasonable
    ! velocities at least once during initialisation to get the thermodynamics right.
    IF (C%choice_ice_dynamics == 'none') THEN
      C%choice_ice_dynamics = 'DIVA'
      CALL initialise_velocity_solver( mesh, ice)
      C%choice_ice_dynamics = 'none'
    ELSE
      CALL initialise_velocity_solver( mesh, ice)
    END IF

    IF (C%is_restart) THEN

      IF (par%master) WRITE(0,*) '   Initialising ice velocities using data read from restart file...'

      ! Set 3D velocity field to restart data
      ice%u_3D_b( mesh%ti1:mesh%ti2,:) = restart%u_3D( mesh%ti1:mesh%ti2,:)
      ice%v_3D_b( mesh%ti1:mesh%ti2,:) = restart%v_3D( mesh%ti1:mesh%ti2,:)
      CALL sync

      ! Calculate 3D vertical velocity from 3D horizontal velocities and conservation of mass
      CALL calc_3D_vertical_velocities( mesh, ice)

      ! Copy surface velocity from the 3D fields
      ice%u_surf_b( mesh%ti1:mesh%ti2) = ice%u_3D_b( mesh%ti1:mesh%ti2,1)
      ice%v_surf_b( mesh%ti1:mesh%ti2) = ice%v_3D_b( mesh%ti1:mesh%ti2,1)
      CALL sync

      ! Copy basal velocity from the 3D fields
      ice%u_base_b( mesh%ti1:mesh%ti2) = ice%u_3D_b( mesh%ti1:mesh%ti2,C%nz)
      ice%v_base_b( mesh%ti1:mesh%ti2) = ice%v_3D_b( mesh%ti1:mesh%ti2,C%nz)
      CALL sync

      ! Calculate vertically averaged velocities
      DO ti = mesh%ti1, mesh%ti2
        prof = ice%u_3D_b( ti,:)
        CALL vertical_average( prof, ice%u_vav_b( ti))
        prof = ice%v_3D_b( ti,:)
        CALL vertical_average( prof, ice%v_vav_b( ti))
      END DO
      CALL sync

      ! Special cases for the little ones
      IF (C%choice_ice_dynamics == 'SIA') THEN
        ! Basal velocity is zero
        ice%u_base_b( mesh%ti1:mesh%ti2) = 0._dp
        ice%v_base_b( mesh%ti1:mesh%ti2) = 0._dp
        CALL sync

      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
        ! No vertical velocity
        ice%w_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
        CALL sync
      END IF

      ! Map velocity components to the a-grid
      CALL map_velocities_b_to_a_3D( mesh, ice%u_3D_b  , ice%v_3D_b  , ice%u_3D_a  , ice%v_3D_a  )
      CALL map_velocities_b_to_a_2D( mesh, ice%u_vav_b , ice%v_vav_b , ice%u_vav_a , ice%v_vav_a )
      CALL map_velocities_b_to_a_2D( mesh, ice%u_surf_b, ice%v_surf_b, ice%u_surf_a, ice%v_surf_a)
      CALL map_velocities_b_to_a_2D( mesh, ice%u_base_b, ice%v_base_b, ice%u_base_a, ice%v_base_a)

      ! Calculate absolute velocities on the b grid
      DO ti = mesh%ti1, mesh%ti2
        ice%uabs_vav_b(  ti) = SQRT( ice%u_vav_b(  ti)**2 + ice%v_vav_b(  ti)**2)
        ice%uabs_surf_b( ti) = SQRT( ice%u_surf_b( ti)**2 + ice%v_surf_b( ti)**2)
        ice%uabs_base_b( ti) = SQRT( ice%u_base_b( ti)**2 + ice%v_base_b( ti)**2)
      END DO

      ! Calculate absolute velocities on the a grid
      DO vi = mesh%vi1, mesh%vi2
        ice%uabs_vav_a(  vi) = SQRT( ice%u_vav_a(  vi)**2 + ice%v_vav_a(  vi)**2)
        ice%uabs_surf_a( vi) = SQRT( ice%u_surf_a( vi)**2 + ice%v_surf_a( vi)**2)
        ice%uabs_base_a( vi) = SQRT( ice%u_base_a( vi)**2 + ice%v_base_a( vi)**2)
      END DO
      CALL sync

    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = HUGE( 1))

  END SUBROUTINE initialise_ice_model

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

    ! Calculate mapping arrays
    CALL create_remapping_arrays_glob_mesh( mesh, forcing%grid_ghf, map)

    ! Map global climate data to the mesh
    CALL map_latlon2mesh_2D( mesh, map, forcing%ghf_ghf, ice%GHF_a)

    ! Deallocate mapping arrays
    CALL deallocate_remapping_arrays_glob_mesh( map)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_geothermal_heat_flux_to_mesh

  SUBROUTINE allocate_ice_model( mesh, ice)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_ice_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory
    ! ===============

    ! Basic data - ice thickness, bedrock height, surface height, mask, 3D ice velocities and temperature
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hi_a                  , ice%wHi_a                 )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hb_a                  , ice%wHb_a                 )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hs_a                  , ice%wHs_a                 )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%SL_a                  , ice%wSL_a                 )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%TAF_a                 , ice%wTAF_a                )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%Ti_a                  , ice%wTi_a                 )

    ! Ice velocities
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%u_3D_a                , ice%wu_3D_a               )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%v_3D_a                , ice%wv_3D_a               )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz       , ice%u_3D_b                , ice%wu_3D_b               )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz       , ice%v_3D_b                , ice%wv_3D_b               )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%w_3D_a                , ice%ww_3D_a               )

    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%u_vav_a               , ice%wu_vav_a              )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%v_vav_a               , ice%wv_vav_a              )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%u_vav_b               , ice%wu_vav_b              )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%v_vav_b               , ice%wv_vav_b              )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%uabs_vav_a            , ice%wuabs_vav_a           )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%uabs_vav_b            , ice%wuabs_vav_b           )

    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%u_surf_a              , ice%wu_surf_a             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%v_surf_a              , ice%wv_surf_a             )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%u_surf_b              , ice%wu_surf_b             )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%v_surf_b              , ice%wv_surf_b             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%w_surf_a              , ice%ww_surf_a             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%uabs_surf_a           , ice%wuabs_surf_a          )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%uabs_surf_b           , ice%wuabs_surf_b          )

    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%u_base_a              , ice%wu_base_a             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%v_base_a              , ice%wv_base_a             )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%u_base_b              , ice%wu_base_b             )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%v_base_b              , ice%wv_base_b             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%w_base_a              , ice%ww_base_a             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%uabs_base_a           , ice%wuabs_base_a          )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%uabs_base_b           , ice%wuabs_base_b          )

    ! Different masks
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_land_a           , ice%wmask_land_a          )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_ocean_a          , ice%wmask_ocean_a         )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_lake_a           , ice%wmask_lake_a          )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_ice_a            , ice%wmask_ice_a           )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_sheet_a          , ice%wmask_sheet_a         )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_shelf_a          , ice%wmask_shelf_a         )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_coast_a          , ice%wmask_coast_a         )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_margin_a         , ice%wmask_margin_a        )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_gl_a             , ice%wmask_gl_a            )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_cf_a             , ice%wmask_cf_a            )
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_a                , ice%wmask_a               )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%f_grnd_a              , ice%wf_grnd_a             )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%f_grnd_b              , ice%wf_grnd_b             )

    ! Ice physical properties
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%A_flow_3D_a           , ice%wA_flow_3D_a          )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%A_flow_vav_a          , ice%wA_flow_vav_a         )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%Ti_pmp_a              , ice%wTi_pmp_a             )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%Cpi_a                 , ice%wCpi_a                )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%Ki_a                  , ice%wKi_a                 )

    ! Zeta derivatives
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dzeta_dt_a            , ice%wdzeta_dt_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dzeta_dx_a            , ice%wdzeta_dx_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dzeta_dy_a            , ice%wdzeta_dy_a           )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dzeta_dz_a            , ice%wdzeta_dz_a           )

    ! Ice dynamics - ice thickness calculation
    CALL allocate_shared_dp_2D(   mesh%nV  , mesh%nC_mem,          ice%dVi_in          , ice%wdVi_in            )
    CALL allocate_shared_dp_2D(   mesh%nV  , mesh%nC_mem,          ice%dVi_out         , ice%wdVi_out           )
    CALL allocate_shared_dp_1D(   mesh%nV  ,                       ice%dHi_dt_a        , ice%wdHi_dt_a          )
    CALL allocate_shared_dp_1D(   mesh%nV  ,                       ice%dHs_dt_a        , ice%wdHs_dt_a          )
    CALL allocate_shared_dp_1D(   mesh%nV  ,                       ice%Hi_tplusdt_a    , ice%wHi_tplusdt_a      )
    CALL allocate_shared_dp_1D(   mesh%nV  ,                       ice%dHi_dt_ave_a    , ice%wdHi_dt_ave_a      )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%dHi_dt_window_size, ice%dHi_dt_window_a , ice%wdHi_dt_window_a   )
    CALL allocate_shared_dp_1D(   mesh%nV  ,                       ice%dHi_dt_past_a   , ice%wdHi_dt_past_a     )
    ! CALL allocate_shared_dp_0D(                                    ice%Hi_tot_old      , ice%wHi_tot_old        )
    ! CALL allocate_shared_dp_0D(                                    ice%Hi_tot_old_time , ice%wHi_tot_old_time   )

    ! Ice dynamics - calving
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%float_margin_frac_a   , ice%wfloat_margin_frac_a  )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hi_eff_cf_a           , ice%wHi_eff_cf_a          )

    ! Ice dynamics - predictor/corrector ice thickness update
    CALL allocate_shared_dp_0D(                           ice%pc_zeta               , ice%wpc_zeta              )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%pc_tau                , ice%wpc_tau               )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%pc_fcb                , ice%wpc_fcb               )
    CALL allocate_shared_dp_0D(                           ice%pc_eta                , ice%wpc_eta               )
    CALL allocate_shared_dp_0D(                           ice%pc_eta_prev           , ice%wpc_eta_prev          )
    CALL allocate_shared_dp_0D(                           ice%pc_beta3              , ice%wpc_beta3             )
    CALL allocate_shared_dp_0D(                           ice%pc_beta4              , ice%wpc_beta4             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%pc_f1                 , ice%wpc_f1                )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%pc_f2                 , ice%wpc_f2                )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%pc_f3                 , ice%wpc_f3                )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%pc_f4                 , ice%wpc_f4                )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hi_old                , ice%wHi_old               )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hi_pred               , ice%wHi_pred              )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hi_corr               , ice%wHi_corr              )

    ! Thermodynamics
    CALL allocate_shared_int_1D(  mesh%nV  ,              ice%mask_ice_a_prev       , ice%wmask_ice_a_prev      )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%internal_heating_a    , ice%winternal_heating_a   )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%frictional_heating_a  , ice%wfrictional_heating_a )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%GHF_a                 , ice%wGHF_a                )

    ! Mesh adaptation data
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%surf_curv             , ice%wsurf_curv            )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%log_velocity          , ice%wlog_velocity         )

    ! GIA
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dHb_a                 , ice%wdHb_a                )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dHb_dt_a              , ice%wdHb_dt_a             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dSL_dt_a              , ice%wdSL_dt_a             )

    ! Useful extra stuff
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dHi_a                 , ice%wdHi_a                )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dHs_a                 , ice%wdHs_a                )

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = HUGE( 1))

  END SUBROUTINE allocate_ice_model

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

    ! ! Save total ice thickness pre-mesh update for future reference
    ! ice%Hi_tot_old = SUM( ice%Hi_a( mesh_old%vi1:mesh_old%vi2))
    ! CALL sync
    ! ! Communicate mid-point results among processes
    ! CALL MPI_ALLREDUCE( MPI_IN_PLACE, ice%Hi_tot_old, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    ! ! Save time of this computation
    ! ice%Hi_tot_old_time = time

    ! The only fields that actually need to be mapped. The rest is either remapped in their
    ! own subroutines or are reallocated and reinitialised to 0.
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%Hi_a        ,    ice%wHi_a        ,    'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dHi_dt_a    ,    ice%wdHi_dt_a    ,    'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%Hi_tplusdt_a,    ice%wHi_tplusdt_a,    'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dHi_dt_ave_a,    ice%wdHi_dt_ave_a,    'cons_2nd_order')
    CALL remap_field_dp_3D( mesh_old, mesh_new, map, ice%dHi_dt_window_a, ice%wdHi_dt_window_a, 'cons_2nd_order')

    ! Remove very thin ice resulting from remapping errors
    DO vi = mesh_new%vi1, mesh_new%vi2
      IF (ice%Hi_a( vi) < 0.1_dp) THEN
        ice%Hi_a( vi) = 0._dp
      END IF
    END DO
    CALL sync

    ! Remove artificial ice at domain border (messes up with thermodynamics)
    DO vi = mesh_new%vi1, mesh_new%vi2
      IF ( .NOT. mesh_new%edge_index( vi) == 0) THEN
        ice%Hi_a( vi) = 0._dp
      END IF
    END DO
    CALL sync

    ! Remap englacial temperature (needs some special attention because of the discontinuity at the ice margin)
    CALL remap_ice_temperature( mesh_old, mesh_new, map, ice)

    ! Remap bedrock change and add up to PD to prevent accumulation of numerical diffusion
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dHb_a, ice%wdHb_a, 'cons_2nd_order')
    CALL reallocate_shared_dp_1D( mesh_new%nV,  ice%Hb_a, ice%wHb_a)
    ice%Hb_a( mesh_new%vi1:mesh_new%vi2) = refgeo_PD%Hb( mesh_new%vi1:mesh_new%vi2) + ice%dHb_a( mesh_new%vi1:mesh_new%vi2)

    ! Geothermal heat flux
    CALL reallocate_shared_dp_1D( mesh_new%nV, ice%GHF_a, ice%wGHF_a)
    IF     (C%choice_geothermal_heat_flux == 'constant') THEN
      ice%GHF_a( mesh_new%vi1:mesh_new%vi2) = C%constant_geothermal_heat_flux
    ELSEIF (C%choice_geothermal_heat_flux == 'spatial') THEN
      CALL map_geothermal_heat_flux_to_mesh( mesh_new, ice)
    ELSE
      CALL crash('unknown choice_geothermal_heat_flux "' // TRIM( C%choice_geothermal_heat_flux) // '"!')
    END IF

    ! Basal conditions
    CALL remap_basal_conditions( mesh_old, mesh_new, map, ice)

    ! GIA and sea level
    CALL reallocate_shared_dp_1D( mesh_new%nV, ice%dHb_a, ice%wdHb_a) ! Remapped above only to compute Hb_a. Recomputed later when needed.
    CALL remap_field_dp_2D(   mesh_old, mesh_new, map, ice%SL_a,     ice%wSL_a,     'cons_2nd_order') ! If not remapped, gets reset to 0 -> not good.
    IF (C%choice_GIA_model == 'SELEN') THEN
      CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dHb_dt_a, ice%wdHb_dt_a, 'cons_2nd_order') ! If not remapped, gets reset to and remains 0
      CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dSL_dt_a, ice%wdSL_dt_a, 'cons_2nd_order') ! until next call to SELEN -> not good.
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dHb_dt_a, ice%wdHb_dt_a, 'cons_2nd_order') ! If not remapped, gets reset to and remains 0
      CALL reallocate_shared_dp_1D( mesh_new%nV, ice%dSL_dt_a, ice%wdSL_dt_a)                         ! until next call to ELRA -> not good.
    ELSEIF (C%choice_GIA_model == 'none') THEN
      CALL reallocate_shared_dp_1D( mesh_new%nV, ice%dHb_dt_a, ice%wdHb_dt_a)
      CALL reallocate_shared_dp_1D( mesh_new%nV, ice%dSL_dt_a, ice%wdSL_dt_a)
    ELSE
      CALL crash('unknown choice_GIA_model "' // TRIM( C%choice_GIA_model) // '"!')
    END IF

    ! Simple memory reallocation for all the rest
    ! ===========================================

    ! Basic data - ice thickness, bedrock height, surface height, mask, 3D ice velocities and temperature
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%Hs_a,                 ice%wHs_a                )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%TAF_a,                ice%wTAF_a               )

    ! Different masks
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_land_a,          ice%wmask_land_a         )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_ocean_a,         ice%wmask_ocean_a        )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_lake_a,          ice%wmask_lake_a         )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_ice_a,           ice%wmask_ice_a          )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_sheet_a,         ice%wmask_sheet_a        )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_shelf_a,         ice%wmask_shelf_a        )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_coast_a,         ice%wmask_coast_a        )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_margin_a,        ice%wmask_margin_a       )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_gl_a,            ice%wmask_gl_a           )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_cf_a,            ice%wmask_cf_a           )
    CALL reallocate_shared_int_1D( mesh_new%nV,                  ice%mask_a,               ice%wmask_a              )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%f_grnd_a,             ice%wf_grnd_a            )
    CALL reallocate_shared_dp_1D(  mesh_new%nTri,                ice%f_grnd_b,             ice%wf_grnd_b            )

    ! Ice physical properties
    CALL reallocate_shared_dp_2D(  mesh_new%nV, C%nz,            ice%A_flow_3D_a,          ice%wA_flow_3D_a         )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%A_flow_vav_a,         ice%wA_flow_vav_a        )
    CALL reallocate_shared_dp_2D(  mesh_new%nV, C%nz,            ice%Ti_pmp_a,             ice%wTi_pmp_a            )
    CALL reallocate_shared_dp_2D(  mesh_new%nV, C%nz,            ice%Cpi_a,                ice%wCpi_a               )
    CALL reallocate_shared_dp_2D(  mesh_new%nV, C%nz,            ice%Ki_a,                 ice%wKi_a                )

    ! Zeta derivatives
    CALL reallocate_shared_dp_2D(  mesh_new%nV, C%nz,            ice%dzeta_dt_a,           ice%wdzeta_dt_a          )
    CALL reallocate_shared_dp_2D(  mesh_new%nV, C%nz,            ice%dzeta_dx_a,           ice%wdzeta_dx_a          )
    CALL reallocate_shared_dp_2D(  mesh_new%nV, C%nz,            ice%dzeta_dy_a,           ice%wdzeta_dy_a          )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%dzeta_dz_a,           ice%wdzeta_dz_a          )

    ! Ice dynamics - ice thickness calculation
    CALL reallocate_shared_dp_2D(  mesh_new%nV, mesh_new%nC_mem, ice%dVi_in,               ice%wdVi_in              )
    CALL reallocate_shared_dp_2D(  mesh_new%nV, mesh_new%nC_mem, ice%dVi_out,              ice%wdVi_out             )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%dHs_dt_a,             ice%wdHs_dt_a            )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%dHi_dt_past_a,        ice%wdHi_dt_past_a       )

   ! Ice dynamics - calving
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%float_margin_frac_a,  ice%wfloat_margin_frac_a )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%Hi_eff_cf_a,          ice%wHi_eff_cf_a         )

    ! Ice dynamics - predictor/corrector ice thickness update
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%pc_tau,               ice%wpc_tau              )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%pc_fcb,               ice%wpc_fcb              )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%pc_f1,                ice%wpc_f1               )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%pc_f2,                ice%wpc_f2               )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%pc_f3,                ice%wpc_f3               )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%pc_f4,                ice%wpc_f4               )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%Hi_old,               ice%wHi_old              )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%Hi_pred,              ice%wHi_pred             )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%Hi_corr,              ice%wHi_corr             )

    ! Thermodynamics
    CALL reallocate_shared_dp_2D(  mesh_new%nV, C%nz,            ice%internal_heating_a,   ice%winternal_heating_a  )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%frictional_heating_a, ice%wfrictional_heating_a)

    ! Mesh adaptation data
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%surf_curv,            ice%wsurf_curv           )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%log_velocity,         ice%wlog_velocity        )

    ! Useful extra stuff
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%dHi_a,                ice%wdHi_a               )
    CALL reallocate_shared_dp_1D(  mesh_new%nV,                  ice%dHs_a,                ice%wdHs_a               )


    ! End of memory reallocation
    ! ==========================

    ! Remap velocities
    ! If we're running with choice_ice_dynamics == "none", remap velocities
    ! pretending we are using the DIVA, just like we did  during initialisation
    ! to get reasonable velocities for the thermodynamics.
    IF (C%choice_ice_dynamics == 'none') THEN
      C%choice_ice_dynamics = 'DIVA'
      CALL remap_velocities( mesh_old, mesh_new, map, ice)
      C%choice_ice_dynamics = 'none'
    ELSE
      CALL remap_velocities( mesh_old, mesh_new, map, ice)
    END IF

    ! Recalculate ice flow factor
    CALL calc_ice_rheology( mesh_new, ice, time)

    ! Update masks
    CALL update_general_ice_model_data( mesh_new, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ice_model

END MODULE ice_dynamics_module
