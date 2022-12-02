MODULE ice_dynamics_module

  ! Contains all the routines needed to calculate ice-sheet geometry at the next time
  ! step, including routines to determine said time step.

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

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
  USE data_types_module,                   ONLY: type_model_region, type_mesh, type_ice_model, type_reference_geometry, type_restart_data
  USE utilities_module,                    ONLY: surface_elevation
  USE ice_velocity_main_module,            ONLY: initialise_velocity_solver, solve_stress_balance, remap_velocity_solver
  USE ice_thickness_module,                ONLY: calc_dHi_dt
  USE thermodynamics_module,               ONLY: calc_ice_rheology, remap_ice_temperature
  USE calving_module,                      ONLY: run_calving_model, imposed_shelf_removal, remove_unconnected_shelves
  USE basal_conditions_and_sliding_module, ONLY: initialise_basal_conditions, remap_basal_conditions
  USE netcdf_input_module,                 ONLY: read_field_from_file_2D
  USE mesh_mapping_module,                 ONLY: remap_field_dp_2D
  USE general_ice_model_data_module,       ONLY: update_general_ice_model_data

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
    INTEGER                                            :: vi1,vi2
    REAL(dp)                                           :: dt_crit_SIA, dt_crit_SSA, r_solver_acc, dt_max

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Abbreviations for cleaner code
    vi1 = region%mesh%vi1
    vi2 = region%mesh%vi2

  ! == Start-up phase

    ! Get a more accurate velocity solution during the start-up phase to prevent initialisation "bumps"
    IF (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
      r_solver_acc = 0.01_dp * 0.99_dp * (region%time - C%start_time_of_run) / C%dt_startup_phase
    ELSE
      r_solver_acc = 1._dp
    END IF

    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! Velocities are set to zero (geometry only subjected to mass balance;
      ! to keep geometry entirely fixed, set choice_ice_integration_method to 'none')
    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      ! No iterative solver involved
    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      region%ice%SSA%PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
      region%ice%SSA%PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc
    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      region%ice%SSA%PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
      region%ice%SSA%PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc
    ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
      region%ice%DIVA%PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
      region%ice%DIVA%PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc
    ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
      region%ice%BPA%PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
      region%ice%BPA%PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc
    ELSE
      CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    END IF

    ! Reduce the time-step during the start-up phase
    IF     (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
      dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((region%time - C%start_time_of_run) / C%dt_startup_phase)**2
    ELSEIF (region%time >= C%end_time_of_run   - C%dt_startup_phase) THEN
      dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((C%end_time_of_run - region%time  ) / C%dt_startup_phase)**2
    ELSE
      dt_max = C%dt_max
    END IF

  ! == Calculate ice velocities with the selected ice-dynamical approximation

    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! Velocities are set to zero (geometry only subjected to mass balance;
      ! to keep geometry entirely fixed, set choice_ice_integration_method to 'none')

    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      ! Shallow ice approximation

      IF (region%time == region%t_next_SIA) THEN

        ! Calculate new ice velocities
        CALL solve_stress_balance( region%mesh, region%ice)

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

    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      ! Shallow shelf approximation

      IF (region%time == region%t_next_SSA) THEN

        ! Calculate new ice velocities
        CALL solve_stress_balance( region%mesh, region%ice)

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

    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      ! Hybrid SIA/SSA (Bueler and Brown, 2009)

      IF (region%time == region%t_next_SIA) THEN

        ! Calculate new ice velocities
        CALL solve_stress_balance( region%mesh, region%ice)

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
        CALL solve_stress_balance( region%mesh, region%ice)

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

      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN

    ELSE   ! IF     (C%choice_stress_balance_approximation == 'SIA') THEN
      CALL crash('"direct" time stepping works only with SIA, SSA, or SIA/SSA ice dynamics, not with DIVA or BPA!')
    END IF ! IF     (C%choice_stress_balance_approximation == 'SIA') THEN

    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps_and_actions( region, t_end)

    !IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_SIA = ', dt_crit_SIA, ', dt_crit_SSA = ', dt_crit_SSA, ', dt = ', region%dt

  ! == Calculate new ice geometry

    ! Calculate dH/dt based on the new velocity solution
    CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt, region%mask_noice, region%refgeo_PD)

    ! Calculate ice thickness at the end of this model loop
    region%ice%Hi_tplusdt_a( vi1:vi2) = MAX( 0._dp, region%ice%Hi_a( vi1:vi2) + region%dt * region%ice%dHi_dt_a( vi1:vi2))
    CALL sync

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
    REAL(dp)                                           :: dt_from_pc, dt_crit_adv, r_solver_acc, dt_max

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Abbreviations for cleaner code
    vi1 = region%mesh%vi1
    vi2 = region%mesh%vi2

  ! == Start-up phase

    ! Get a more accurate velocity solution during the start-up phase to prevent initialisation "bumps"
    IF (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
      r_solver_acc = 0.01_dp * 0.99_dp * (region%time - C%start_time_of_run) / C%dt_startup_phase
    ELSE
      r_solver_acc = 1._dp
    END IF

    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! Velocities are set to zero (geometry only subjected to mass balance;
      ! to keep geometry entirely fixed, set choice_ice_integration_method to 'none')
    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      ! No iterative solver involved
    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      region%ice%SSA%PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
      region%ice%SSA%PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc
    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      region%ice%SSA%PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
      region%ice%SSA%PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc
    ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
      region%ice%DIVA%PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
      region%ice%DIVA%PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc
    ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
      region%ice%BPA%PETSc_rtol   = C%DIVA_PETSc_rtol   * r_solver_acc
      region%ice%BPA%PETSc_abstol = C%DIVA_PETSc_abstol * r_solver_acc
    ELSE
      CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    END IF

    ! Reduce the time-step during the start-up phase
    IF     (region%time <= C%start_time_of_run + C%dt_startup_phase) THEN
      dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((region%time - C%start_time_of_run) / C%dt_startup_phase)**2
    ELSEIF (region%time >= C%end_time_of_run   - C%dt_startup_phase) THEN
      dt_max = C%dt_min + (C%dt_max - C%dt_min) * ((C%end_time_of_run - region%time  ) / C%dt_startup_phase)**2
    ELSE
      dt_max = C%dt_max
    END IF

  ! == Update ice velocities

    ! Determine whether or not we need to update ice velocities
    do_update_ice_velocity = .FALSE.
    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! Velocities are set to zero (geometry only subjected to mass balance;
      ! to keep geometry entirely fixed, set choice_ice_integration_method to 'none')
    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      IF (region%time == region%t_next_SIA ) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      IF (region%time == region%t_next_SSA ) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      IF (region%time == region%t_next_SIA ) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
      IF (region%time == region%t_next_DIVA) do_update_ice_velocity = .TRUE.
    ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
      IF (region%time == region%t_next_BPA ) do_update_ice_velocity = .TRUE.
    ELSE
      CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    END IF

    IF (do_update_ice_velocity) THEN

      ! Calculate time step based on the truncation error in ice thickness (Robinson et al., 2020, Eq. 33)
      CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_adv)
      IF (par%master) THEN
        region%dt_crit_ice_prev = region%dt_crit_ice
        region%ice%pc_eta_prev  = region%ice%pc_eta
        dt_from_pc              = (C%pc_epsilon / region%ice%pc_eta)**(C%pc_k_I + C%pc_k_p) * (C%pc_epsilon / region%ice%pc_eta_prev)**(-C%pc_k_p) * region%dt
        region%dt_crit_ice      = MAX(C%dt_min, MINVAL([ dt_max, 2._dp * region%dt_crit_ice_prev, dt_crit_adv, dt_from_pc]))
        region%ice%pc_zeta      = region%dt_crit_ice / region%dt_crit_ice_prev
        region%ice%pc_beta1     = 1._dp + region%ice%pc_zeta / 2._dp
        region%ice%pc_beta2     =       - region%ice%pc_zeta / 2._dp
      END IF
      CALL sync

      ! Predictor step
      ! ==============

      ! Calculate new ice geometry
      region%ice%pc_f2(   vi1:vi2) = region%ice%dHi_dt_a( vi1:vi2)
      CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice, region%refgeo_PD)
      region%ice%pc_f1(   vi1:vi2) = region%ice%dHi_dt_a( vi1:vi2)
      region%ice%Hi_pred( vi1:vi2) = MAX(0._dp, region%ice%Hi_a(     vi1:vi2) + region%dt_crit_ice * region%ice%dHi_dt_a( vi1:vi2))
      CALL sync

      ! Update step
      ! ===========

      ! Solve the stress balance and calculate ice velocities for predicted geometry
      region%ice%Hi_old( vi1:vi2) = region%ice%Hi_a(    vi1:vi2)
      region%ice%Hi_a(   vi1:vi2) = region%ice%Hi_pred( vi1:vi2)
      CALL update_general_ice_model_data( region%mesh, region%ice)
      CALL solve_stress_balance( region%mesh, region%ice)

      ! Update timer

      IF     (C%choice_stress_balance_approximation == 'none') THEN
        ! Velocities are set to zero (geometry only subjected to mass balance;
        ! to keep geometry entirely fixed, set choice_ice_integration_method to 'none')
      ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_ice
        CALL sync
      ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_ice
        CALL sync
      ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_ice
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_ice
        CALL sync
      ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
        IF (par%master) region%t_last_DIVA = region%time
        IF (par%master) region%t_next_DIVA = region%time + region%dt_crit_ice
        CALL sync
      ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
        IF (par%master) region%t_last_BPA = region%time
        IF (par%master) region%t_next_BPA = region%time + region%dt_crit_ice
        CALL sync
      ELSE
        CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
      END IF

      ! Corrector step
      ! ==============

      ! Go back to old ice thickness. Run all the other modules (climate, SMB, BMB, thermodynamics, etc.)
      ! and only go to new (corrected) ice thickness at the end of this time loop.
      region%ice%Hi_a(    vi1:vi2) = region%ice%Hi_old(   vi1:vi2)
      CALL update_general_ice_model_data( region%mesh, region%ice)

      ! Calculate "corrected" ice thickness based on new velocities
      CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt_crit_ice, region%mask_noice, region%refgeo_PD)
      region%ice%pc_f3(   vi1:vi2) = region%ice%dHi_dt_a( vi1:vi2)
      region%ice%pc_f4(   vi1:vi2) = region%ice%pc_f1(    vi1:vi2)
      region%ice%Hi_corr( vi1:vi2) = MAX(0._dp, region%ice%Hi_a( vi1:vi2) + 0.5_dp * region%dt_crit_ice * (region%ice%pc_f3( vi1:vi2) + region%ice%pc_f4( vi1:vi2)))
      CALL sync

      ! Determine truncation error
      CALL calc_pc_truncation_error( region%mesh, region%ice, region%dt_crit_ice, region%dt_prev)

    END IF ! IF (do_update_ice_velocity) THEN

    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps_and_actions( region, t_end)

    ! Calculate ice thickness at the end of this model loop
    region%ice%Hi_tplusdt_a( vi1:vi2) = MAX( 0._dp, region%ice%Hi_a( vi1:vi2) + region%dt * region%ice%dHi_dt_a( vi1:vi2))
    CALL sync

    !IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_adv = ', dt_crit_adv, ', dt_from_pc = ', dt_from_pc, ', dt = ', region%dt

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
    INTEGER                                            :: ti, via, vib, vic
    REAL(dp)                                           :: d_ab, d_bc, d_ca, d_min, Hi, D_SIA, dt
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise time step with maximum allowed value
    dt_crit_SIA = C%dt_max

    DO ti = mesh%ti1, mesh%ti2

      ! Calculate shortest triangle side
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      d_ab = NORM2( mesh%V( vib,:) - mesh%V( via,:))
      d_bc = NORM2( mesh%V( vic,:) - mesh%V( vib,:))
      d_ca = NORM2( mesh%V( via,:) - mesh%V( vic,:))

      d_min = MINVAL([ d_ab, d_bc, d_ca])

      ! Find maximum diffusivity in the vertical column
      D_SIA = MAXVAL( ABS( ice%SIA%D_3D_b( ti,:)))

      ! Calculate critical timestep
      Hi = MAXVAL( [ice%Hi_a( via), ice%Hi_a( vib), ice%Hi_a( vic)])
      dt = d_min**2 / (6._dp * Hi * D_SIA)
      dt_crit_SIA = MIN( dt_crit_SIA, dt)

    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_SIA, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_SIA = dt_crit_SIA * dt_correction_factor

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
    REAL(dp)                                           :: dist, dt
    REAL(dp), PARAMETER                                :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.

    ! Add routine to path
    CALL init_routine( routine_name)

    dt_crit_adv = 2._dp * C%dt_max

    DO aci = mesh%ci1, mesh%ci2

      ! Only check at ice-covered vertices
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)
      IF (ice%Hi_a( vi) == 0._dp .OR. ice%Hi_a( vj) == 0._dp) CYCLE

      dist = NORM2( mesh%V( vi,:) - mesh%V( vj,:))
      dt = dist / (ABS( ice%u_vav_c( aci)) + ABS( ice%v_vav_c( aci)))
      dt_crit_adv = MIN( dt_crit_adv, dt)

    END DO

    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = MIN(C%dt_max, dt_crit_adv * dt_correction_factor)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_critical_timestep_adv

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

    ! Ratio of time steps
    zeta = dt / dt_prev

    ! Find maximum truncation error
    eta_proc = C%pc_eta_min

    DO vi = mesh%vi1, mesh%vi2

      ! Calculate truncation error (Robinson et al., 2020, Eq. 32)
      ice%pc_tau( vi) = ABS( zeta * (ice%Hi_corr( vi) - ice%Hi_pred( vi)) / ((3._dp * zeta + 3._dp) * dt))

      IF (ice%mask_sheet_a( vi) == 1) THEN

        has_GL_neighbour = .FALSE.
        DO ci = 1, mesh%nC( vi)
          vc = mesh%C( vi,ci)
          IF (ice%mask_gl_a( vc) == 1) THEN
            has_GL_neighbour = .TRUE.
            EXIT
          END IF
        END DO

        IF (.NOT.has_GL_neighbour) eta_proc = MAX( eta_proc, ice%pc_tau( vi))

      END IF

    END DO
    CALL MPI_REDUCE( eta_proc, ice%pc_eta, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pc_truncation_error

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

    IF (par%master) THEN

      ! Determine when each model components should be updated

      t_next = MIN(t_end, region%time + C%dt_max)

      ! First the ice dynamics
      ! ======================

      IF     (C%choice_stress_balance_approximation == 'none') THEN
        ! Just stick to the maximum time step
      ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
        t_next = MIN( t_next, region%t_next_SIA)
      ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
        t_next = MIN( t_next, region%t_next_SSA)
      ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
        t_next = MIN( t_next, region%t_next_SIA)
        t_next = MIN( t_next, region%t_next_SSA)
      ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
        t_next = MIN( t_next, region%t_next_DIVA)
      ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
        t_next = MIN( t_next, region%t_next_BPA)
      ELSE
        CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
      END IF ! IF (C%choice_stress_balance_approximation == 'SIA') THEN

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
      IF (C%choice_GIA_model == 'ELRA') THEN
        IF (region%time == region%t_next_ELRA) THEN
          region%do_ELRA      = .TRUE.
          region%t_last_ELRA  = region%time
          region%t_next_ELRA  = region%t_last_ELRA + C%dt_bedrock_ELRA
        END IF
        t_next = MIN( t_next, region%t_next_ELRA)
      END IF

      region%do_basal    = .FALSE.
      IF (region%time == region%t_next_basal) THEN
        region%do_basal       = .TRUE.
        region%t_last_basal   = region%time
        region%t_next_basal   = region%t_last_basal + C%BIVgeo_dt
      END IF
      t_next = MIN( t_next, region%t_next_basal)

      region%do_SMB_inv = .FALSE.
      IF (region%time == region%t_next_SMB_inv) THEN
        region%do_SMB_inv       = .TRUE.
        region%t_last_SMB_inv   = region%time
        region%t_next_SMB_inv   = region%t_last_SMB_inv + C%dt_SMB_inv
      END IF
      t_next = MIN( t_next, region%t_next_SMB_inv)

      region%do_output  = .FALSE.
      IF (region%time == region%t_next_output) THEN
        region%do_output      = .TRUE.
        region%t_last_output  = region%time
        region%t_next_output  = region%t_last_output + C%dt_output
      END IF
      t_next = MIN( t_next, region%t_next_output)

      ! Set time step so that we move forward to the next action
      region%dt = t_next - region%time

    END IF ! IF (par%master) THEN
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
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    TYPE(type_restart_data),             INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_model'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '  Initialising ice dynamics model...'

    ! Allocate shared memory
    CALL allocate_ice_model( mesh, ice)

    ! Initialise with data from initial file (works for both "raw" initial data and model restarts)
    DO vi = mesh%vi1, mesh%vi2
      ! Main quantities
      ice%Hi_a( vi) = refgeo_init%Hi( vi)
      ice%Hb_a( vi) = refgeo_init%Hb( vi)
      ice%Hs_a( vi) = surface_elevation( ice%Hi_a( vi), ice%Hb_a( vi), 0._dp)

      ! Differences w.r.t. present-day
      ice%dHi_a( vi) = ice%Hi_a( vi) - refgeo_PD%Hi( vi)
      ice%dHb_a( vi) = ice%Hb_a( vi) - refgeo_PD%Hb( vi)
      ice%dHs_a( vi) = ice%Hs_a( vi) - refgeo_PD%Hs( vi)
    END DO
    CALL sync

    ! Initialise masks and slopes
    CALL update_general_ice_model_data( mesh, ice)

    ! Allocate and initialise basal conditions
    CALL initialise_basal_conditions( mesh, ice, restart)

    ! Geothermal heat flux
    CALL initialise_geothermal_heat_flux( mesh, ice)

    ! Initialise data and matrices for the velocity solver(s)
    CALL initialise_velocity_solver( mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = HUGE( 1))

  END SUBROUTINE initialise_ice_model

  SUBROUTINE initialise_geothermal_heat_flux( mesh, ice)
    ! Initialise the 2-D geothermal heat flux, either with a uniform value
    ! or from an external file.
    !
    ! Assumes memory has already been allocated.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_geothermal_heat_flux'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_geothermal_heat_flux == 'constant' .OR. C%choice_geothermal_heat_flux == 'uniform') THEN
      ! Use a spatially uniform value
      ice%GHF_a( mesh%vi1:mesh%vi2) = C%constant_geothermal_heat_flux
    ELSEIF (C%choice_geothermal_heat_flux == 'spatial') THEN
      ! Read a 2-D field from an external file
      CALL read_field_from_file_2D( C%filename_geothermal_heat_flux, 'hflux', mesh, ice%GHF_a, 'ANT')
    ELSE
      CALL crash('unknown choice_geothermal_heat_flux "' // TRIM( C%choice_geothermal_heat_flux) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_geothermal_heat_flux

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
    CALL allocate_shared_dp_1D(   mesh%nAc ,              ice%u_vav_c               , ice%wu_vav_c              )
    CALL allocate_shared_dp_1D(   mesh%nAc ,              ice%v_vav_c               , ice%wv_vav_c              )
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

    ! Strain rates
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%du_dx_3D_a            , ice%wdu_dx_3D_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%du_dy_3D_a            , ice%wdu_dy_3D_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%du_dz_3D_a            , ice%wdu_dz_3D_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dv_dx_3D_a            , ice%wdv_dx_3D_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dv_dy_3D_a            , ice%wdv_dy_3D_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dv_dz_3D_a            , ice%wdv_dz_3D_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dw_dx_3D_a            , ice%wdw_dx_3D_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dw_dy_3D_a            , ice%wdw_dy_3D_a           )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dw_dz_3D_a            , ice%wdw_dz_3D_a           )

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

    ! Basal friction coefficient
    CALL allocate_shared_dp_1D(   mesh%nV               , ice%beta_b_a              , ice%wbeta_b_a             )

    ! Zeta gradients
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dzeta_dt_ak           , ice%wdzeta_dt_ak          )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dzeta_dx_ak           , ice%wdzeta_dx_ak          )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dzeta_dy_ak           , ice%wdzeta_dy_ak          )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%dzeta_dz_ak           , ice%wdzeta_dz_ak          )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%d2zeta_dx2_ak         , ice%wd2zeta_dx2_ak        )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%d2zeta_dxdy_ak        , ice%wd2zeta_dxdy_ak       )
    CALL allocate_shared_dp_2D(   mesh%nV  , C%nz       , ice%d2zeta_dy2_ak         , ice%wd2zeta_dy2_ak        )

    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz       , ice%dzeta_dx_bk           , ice%wdzeta_dx_bk          )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz       , ice%dzeta_dy_bk           , ice%wdzeta_dy_bk          )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz       , ice%dzeta_dz_bk           , ice%wdzeta_dz_bk          )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz       , ice%d2zeta_dx2_bk         , ice%wd2zeta_dx2_bk        )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz       , ice%d2zeta_dxdy_bk        , ice%wd2zeta_dxdy_bk       )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz       , ice%d2zeta_dy2_bk         , ice%wd2zeta_dy2_bk        )

    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz-1     , ice%dzeta_dx_bks          , ice%wdzeta_dx_bks         )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz-1     , ice%dzeta_dy_bks          , ice%wdzeta_dy_bks         )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz-1     , ice%dzeta_dz_bks          , ice%wdzeta_dz_bks         )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz-1     , ice%d2zeta_dx2_bks        , ice%wd2zeta_dx2_bks       )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz-1     , ice%d2zeta_dxdy_bks       , ice%wd2zeta_dxdy_bks      )
    CALL allocate_shared_dp_2D(   mesh%nTri, C%nz-1     , ice%d2zeta_dy2_bks        , ice%wd2zeta_dy2_bks       )

    ! Ice dynamics - ice thickness calculation
    CALL allocate_shared_dp_2D(   mesh%nV  , mesh%nC_mem, ice%dVi_in                , ice%wdVi_in               )
    CALL allocate_shared_dp_2D(   mesh%nV  , mesh%nC_mem, ice%dVi_out               , ice%wdVi_out              )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dHi_dt_a              , ice%wdHi_dt_a             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dHs_dt_a              , ice%wdHs_dt_a             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hi_tplusdt_a          , ice%wHi_tplusdt_a         )

    ! Ice dynamics - calving
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%float_margin_frac_a   , ice%wfloat_margin_frac_a  )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hi_eff_cf_a           , ice%wHi_eff_cf_a          )

    ! Ice dynamics - predictor/corrector ice thickness update
    CALL allocate_shared_dp_0D(                           ice%pc_zeta               , ice%wpc_zeta              )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%pc_tau                , ice%wpc_tau               )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%pc_fcb                , ice%wpc_fcb               )
    CALL allocate_shared_dp_0D(                           ice%pc_eta                , ice%wpc_eta               )
    CALL allocate_shared_dp_0D(                           ice%pc_eta_prev           , ice%wpc_eta_prev          )
    CALL allocate_shared_dp_0D(                           ice%pc_beta1              , ice%wpc_beta1             )
    CALL allocate_shared_dp_0D(                           ice%pc_beta2              , ice%wpc_beta2             )
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

  SUBROUTINE remap_ice_model( mesh_old, mesh_new, ice, refgeo_PD, time)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_ice_model'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! The only fields that actually need to be mapped. The rest only needs memory reallocation.
    CALL remap_field_dp_2D( mesh_old, mesh_new, ice%Hi_a        , ice%wHi_a        , 'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, ice%dHi_dt_a    , ice%wdHi_dt_a    , 'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, ice%Hi_tplusdt_a, ice%wHi_tplusdt_a, 'cons_2nd_order')
   !CALL remap_field_dp_3D( mesh_old, mesh_new, ice%Ti_a        , ice%wTi_a        , 'cons_2nd_order')

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
    CALL remap_ice_temperature( mesh_old, mesh_new, ice)

    ! Remap bedrock change and add up to PD to prevent accumulation of numerical diffusion
    CALL remap_field_dp_2D( mesh_old, mesh_new, ice%dHb_a,               ice%wdHb_a,               'cons_2nd_order')
    CALL reallocate_shared_dp_1D(    mesh_new%nV,  ice%Hb_a,                   ice%wHb_a                     )
    ice%Hb_a( mesh_new%vi1:mesh_new%vi2) = refgeo_PD%Hb( mesh_new%vi1:mesh_new%vi2) + ice%dHb_a( mesh_new%vi1:mesh_new%vi2)

    ! Geothermal heat flux
    CALL reallocate_shared_dp_1D( mesh_new%nV, ice%GHF_a, ice%wGHF_a)
    CALL initialise_geothermal_heat_flux( mesh_new, ice)

    ! Basal conditions
    CALL remap_basal_conditions( mesh_old, mesh_new, ice)

    ! GIA and sea level
    CALL reallocate_shared_dp_1D( mesh_new%nV, ice%dHb_a, ice%wdHb_a) ! Remapped above only to compute Hb_a. Recomputed later when needed.
    CALL remap_field_dp_2D(   mesh_old, mesh_new, ice%SL_a,     ice%wSL_a,     'cons_2nd_order') ! If not remapped, gets reset to 0 -> not good.
    IF (C%choice_GIA_model == 'SELEN') THEN
      CALL remap_field_dp_2D( mesh_old, mesh_new, ice%dHb_dt_a, ice%wdHb_dt_a, 'cons_2nd_order') ! If not remapped, gets reset to and remains 0
      CALL remap_field_dp_2D( mesh_old, mesh_new, ice%dSL_dt_a, ice%wdSL_dt_a, 'cons_2nd_order') ! until next call to SELEN -> not good.
    ELSEIF (C%choice_GIA_model == 'ELRA') THEN
      CALL remap_field_dp_2D( mesh_old, mesh_new, ice%dHb_dt_a, ice%wdHb_dt_a, 'cons_2nd_order') ! If not remapped, gets reset to and remains 0
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
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hi_a                  , ice%wHi_a                 )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hb_a                  , ice%wHb_a                 )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hs_a                  , ice%wHs_a                 )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%SL_a                  , ice%wSL_a                 )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%TAF_a                 , ice%wTAF_a                )
   !CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%Ti_a                  , ice%wTi_a                 )

    ! Ice velocities
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%u_3D_a                , ice%wu_3D_a               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%v_3D_a                , ice%wv_3D_a               )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%u_3D_b                , ice%wu_3D_b               )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%v_3D_b                , ice%wv_3D_b               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%w_3D_a                , ice%ww_3D_a               )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_vav_a               , ice%wu_vav_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_vav_a               , ice%wv_vav_a              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_vav_b               , ice%wu_vav_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_vav_b               , ice%wv_vav_b              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_vav_a            , ice%wuabs_vav_a           )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_vav_b            , ice%wuabs_vav_b           )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_surf_a              , ice%wu_surf_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_surf_a              , ice%wv_surf_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_surf_b              , ice%wu_surf_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_surf_b              , ice%wv_surf_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%w_surf_a              , ice%ww_surf_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_surf_a           , ice%wuabs_surf_a          )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_surf_b           , ice%wuabs_surf_b          )

    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_base_a              , ice%wu_base_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_base_a              , ice%wv_base_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_base_b              , ice%wu_base_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_base_b              , ice%wv_base_b             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%w_base_a              , ice%ww_base_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_base_a           , ice%wuabs_base_a          )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_base_b           , ice%wuabs_base_b          )

    ! Strain rates
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%du_dx_3D_a            , ice%wdu_dx_3D_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%du_dy_3D_a            , ice%wdu_dy_3D_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%du_dz_3D_a            , ice%wdu_dz_3D_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dv_dx_3D_a            , ice%wdv_dx_3D_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dv_dy_3D_a            , ice%wdv_dy_3D_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dv_dz_3D_a            , ice%wdv_dz_3D_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dw_dx_3D_a            , ice%wdw_dx_3D_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dw_dy_3D_a            , ice%wdw_dy_3D_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dw_dz_3D_a            , ice%wdw_dz_3D_a           )

    ! Different masks
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_land_a           , ice%wmask_land_a          )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_ocean_a          , ice%wmask_ocean_a         )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_lake_a           , ice%wmask_lake_a          )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_ice_a            , ice%wmask_ice_a           )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_sheet_a          , ice%wmask_sheet_a         )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_shelf_a          , ice%wmask_shelf_a         )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_coast_a          , ice%wmask_coast_a         )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_margin_a         , ice%wmask_margin_a        )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_gl_a             , ice%wmask_gl_a            )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_cf_a             , ice%wmask_cf_a            )
    CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_a                , ice%wmask_a               )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%f_grnd_a              , ice%wf_grnd_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%f_grnd_b              , ice%wf_grnd_b             )

    ! Ice physical properties
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%A_flow_3D_a           , ice%wA_flow_3D_a          )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%A_flow_vav_a          , ice%wA_flow_vav_a         )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%Ti_pmp_a              , ice%wTi_pmp_a             )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%Cpi_a                 , ice%wCpi_a                )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%Ki_a                  , ice%wKi_a                 )

    ! Basal friction coefficient
    CALL reallocate_shared_dp_1D(   mesh_new%nV                   , ice%beta_b_a              , ice%wbeta_b_a             )

    ! Zeta gradients
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dzeta_dt_ak           , ice%wdzeta_dt_ak          )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dzeta_dx_ak           , ice%wdzeta_dx_ak          )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dzeta_dy_ak           , ice%wdzeta_dy_ak          )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dzeta_dz_ak           , ice%wdzeta_dz_ak          )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%d2zeta_dx2_ak         , ice%wd2zeta_dx2_ak        )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%d2zeta_dxdy_ak        , ice%wd2zeta_dxdy_ak       )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%d2zeta_dy2_ak         , ice%wd2zeta_dy2_ak        )

    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%dzeta_dx_bk           , ice%wdzeta_dx_bk          )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%dzeta_dy_bk           , ice%wdzeta_dy_bk          )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%dzeta_dz_bk           , ice%wdzeta_dz_bk          )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%d2zeta_dx2_bk         , ice%wd2zeta_dx2_bk        )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%d2zeta_dxdy_bk        , ice%wd2zeta_dxdy_bk       )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%d2zeta_dy2_bk         , ice%wd2zeta_dy2_bk        )

    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz-1         , ice%dzeta_dx_bks          , ice%wdzeta_dx_bks         )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz-1         , ice%dzeta_dy_bks          , ice%wdzeta_dy_bks         )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz-1         , ice%dzeta_dz_bks          , ice%wdzeta_dz_bks         )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz-1         , ice%d2zeta_dx2_bks        , ice%wd2zeta_dx2_bks       )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz-1         , ice%d2zeta_dxdy_bks       , ice%wd2zeta_dxdy_bks      )
    CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz-1         , ice%d2zeta_dy2_bks        , ice%wd2zeta_dy2_bks       )

    ! Ice dynamics - ice thickness calculation
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , mesh_new%nC_mem, ice%dVi_in                , ice%wdVi_in               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , mesh_new%nC_mem, ice%dVi_out               , ice%wdVi_out              )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dHi_dt_a              , ice%wdHi_dt_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dHs_dt_a              , ice%wdHs_dt_a             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hi_tplusdt_a          , ice%wHi_tplusdt_a         )

   ! Ice dynamics - calving
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%float_margin_frac_a   , ice%wfloat_margin_frac_a  )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hi_eff_cf_a           , ice%wHi_eff_cf_a          )

    ! Ice dynamics - predictor/corrector ice thickness update
   !CALL reallocate_shared_dp_0D(                                   ice%pc_zeta               , ice%wpc_zeta              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%pc_tau                , ice%wpc_tau               )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%pc_fcb                , ice%wpc_fcb               )
   !CALL reallocate_shared_dp_0D(                                   ice%pc_eta                , ice%wpc_eta               )
   !CALL reallocate_shared_dp_0D(                                   ice%pc_eta_prev           , ice%wpc_eta_prev          )
   !CALL reallocate_shared_dp_0D(                                   ice%pc_beta1              , ice%wpc_beta1             )
   !CALL reallocate_shared_dp_0D(                                   ice%pc_beta2              , ice%wpc_beta2             )
   !CALL reallocate_shared_dp_0D(                                   ice%pc_beta3              , ice%wpc_beta3             )
   !CALL reallocate_shared_dp_0D(                                   ice%pc_beta4              , ice%wpc_beta4             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%pc_f1                 , ice%wpc_f1                )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%pc_f2                 , ice%wpc_f2                )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%pc_f3                 , ice%wpc_f3                )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%pc_f4                 , ice%wpc_f4                )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hi_old                , ice%wHi_old               )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hi_pred               , ice%wHi_pred              )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hi_corr               , ice%wHi_corr              )

    ! Thermodynamics
   !CALL reallocate_shared_int_1D(  mesh_new%nV  ,                  ice%mask_ice_a_prev       , ice%wmask_ice_a_prev      )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%internal_heating_a    , ice%winternal_heating_a   )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%frictional_heating_a  , ice%wfrictional_heating_a )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%GHF_a                 , ice%wGHF_a                )

    ! Mesh adaptation data
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%surf_curv             , ice%wsurf_curv            )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%log_velocity          , ice%wlog_velocity         )

    ! Useful extra stuff
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dHi_a                 , ice%wdHi_a                )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dHs_a                 , ice%wdHs_a                )


    ! End of memory reallocation
    ! ==========================

    ! Remap velocities
    CALL remap_velocity_solver( mesh_old, mesh_new, ice)

    ! Recalculate ice flow factor
    CALL calc_ice_rheology( mesh_new, ice, time)

    ! Update masks
    CALL update_general_ice_model_data( mesh_new, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ice_model

END MODULE ice_dynamics_module
