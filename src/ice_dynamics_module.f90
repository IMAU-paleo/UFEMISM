MODULE ice_dynamics_module

  ! Contains all the routines needed to calculate ice-sheet geometry at the next time
  ! step, including routines to determine said time step.
  ! NOTE: routines for calculating ice velocities             have been moved to the ice_velocity_module;
  !       routines for integrating the ice thickness equation have been moved to the ice_thickness_module.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
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
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  
  ! Import specific functionality                     
  USE data_types_module,               ONLY: type_model_region, type_mesh, type_ice_model, type_reference_geometry, &
                                             type_remapping_mesh_mesh
  USE utilities_module,                ONLY: vertical_average, surface_elevation                 
  USE mesh_mapping_module,             ONLY: remap_field_dp_2D, remap_field_dp_3D
  USE general_ice_model_data_module,   ONLY: update_general_ice_model_data
  USE mesh_operators_module,           ONLY: map_a_to_c_2D, ddx_a_to_c_2D, ddy_a_to_c_2D
  USE ice_velocity_module,             ONLY: solve_SIA, solve_SSA, solve_DIVA, initialise_velocity_solver, remap_velocities, &
                                             map_velocities_b_to_c_2D, map_velocities_b_to_c_3D
  USE ice_thickness_module,            ONLY: calc_dHi_dt
  USE basal_conditions_and_sliding_module, ONLY: initialise_basal_conditions, remap_basal_conditions
  USE thermodynamics_module,           ONLY: calc_ice_rheology, remap_ice_temperature

  IMPLICIT NONE
  
CONTAINS

! == The main ice dynamics routine
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
    REAL(dp)                                           :: dt_crit_SIA, dt_crit_SSA
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
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
        IF (par%master) region%dt_crit_SIA = dt_crit_SIA
        
        ! Update timer
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_SIA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN
      
    ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
      ! Shallow shelf approximation
    
      IF (region%time == region%t_next_SSA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SSA( region%mesh, region%ice)
    
        ! Calculate critical time step
        CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_SSA)
        IF (par%master) region%dt_crit_SSA = dt_crit_SSA
        
        ! Update timer
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_SSA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SSA) < dt_tol) THEN
      
    ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
      ! Hybrid SIA/SSA (Bueler and Brown, 2009)
    
      IF (region%time == region%t_next_SIA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SIA( region%mesh, region%ice)
    
        ! Calculate critical time step
        CALL calc_critical_timestep_SIA( region%mesh, region%ice, dt_crit_SIA)
        IF (par%master) region%dt_crit_SIA = dt_crit_SIA
        
        ! Update timer
        IF (par%master) region%t_last_SIA = region%time
        IF (par%master) region%t_next_SIA = region%time + region%dt_crit_SIA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN
      
      IF (region%time == region%t_next_SSA) THEN
    
        ! Calculate new ice velocities
        CALL solve_SSA( region%mesh, region%ice)
    
        ! Calculate critical time step
        CALL calc_critical_timestep_adv( region%mesh, region%ice, dt_crit_SSA)
        IF (par%master) region%dt_crit_SSA = dt_crit_SSA
        
        ! Update timer
        IF (par%master) region%t_last_SSA = region%time
        IF (par%master) region%t_next_SSA = region%time + region%dt_crit_SSA
        CALL sync
        
      END IF ! IF (ABS(region%time - region%t_next_SIA) < dt_tol) THEN
      
    ELSE ! IF     (C%choice_ice_dynamics == 'SIA') THEN
      CALL crash('"direct" time stepping works only with SIA, SSA, or SIA/SSA ice dynamics, not with DIVA!')
    END IF ! IF     (C%choice_ice_dynamics == 'SIA') THEN
    
    ! Adjust the time step to prevent overshooting other model components (thermodynamics, SMB, output, etc.)
    CALL determine_timesteps_and_actions( region, t_end)
      
    !IF (par%master) WRITE(0,'(A,F7.4,A,F7.4,A,F7.4)') 'dt_crit_SIA = ', dt_crit_SIA, ', dt_crit_SSA = ', dt_crit_SSA, ', dt = ', region%dt
    
    ! Calculate new ice geometry
    ! ==========================
    
    IF (C%choice_ice_dynamics == 'none') THEN
      ! Fixed ice geometry
      region%ice%dHi_dt_a( region%mesh%vi1 : region%mesh%vi2) = 0._dp
      CALL sync
    ELSE
      CALL calc_dHi_dt( region%mesh, region%ice, region%SMB, region%BMB, region%dt, region%mask_noice, region%refgeo_PD)
    END IF
    
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
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Abbreviations for cleaner code
    vi1 = region%mesh%vi1
    vi2 = region%mesh%vi2
    
    ! Determine whether or not we need to update ice velocities
    do_update_ice_velocity = .FALSE.
    IF     (C%choice_ice_dynamics == 'none') THEN
      region%ice%dHi_dt_a( vi1:vi2) = 0._dp
      CALL sync
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
      IF (par%master) THEN
        region%dt_crit_ice_prev = region%dt_crit_ice
        region%ice%pc_eta_prev  = region%ice%pc_eta
        dt_from_pc              = (C%pc_epsilon / region%ice%pc_eta)**(C%pc_k_I + C%pc_k_p) * (C%pc_epsilon / region%ice%pc_eta_prev)**(-C%pc_k_p) * region%dt
        region%dt_crit_ice      = MAX(C%dt_min, MINVAL([ C%dt_max, 2._dp * region%dt_crit_ice_prev, dt_crit_adv, dt_from_pc]))
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
  
! == Update the ice thickness at the end of a model time loop
  SUBROUTINE update_ice_thickness( mesh, ice)
    ! Update the ice thickness at the end of a model time loop
    
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_ice_thickness'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Save the previous ice mask, for use in thermodynamics
    ice%mask_ice_a_prev( mesh%vi1:mesh%vi2) = ice%mask_ice_a( mesh%vi1:mesh%vi2)
    CALL sync
    
    ! Set ice thickness to new value
    ice%Hi_a( mesh%vi1:mesh%vi2) = MAX( 0._dp, ice%Hi_tplusdt_a( mesh%vi1:mesh%vi2))
    CALL sync
    
    ! Apply calving law
    ! 
!    ! NOTE: done twice, so that the calving front is also allowed to retreat
!    IF (.NOT. C%choice_calving_law == 'none') THEN
!      CALL determine_masks_ice(                mesh, ice)
!      CALL determine_masks_transitions(        mesh, ice)
!      CALL determine_floating_margin_fraction( mesh, ice)
!      CALL apply_calving_law(                  mesh, ice)
!      
!      CALL determine_masks_ice(                mesh, ice)
!      CALL determine_masks_transitions(        mesh, ice)
!      CALL determine_floating_margin_fraction( mesh, ice)
!      CALL apply_calving_law(                  mesh, ice)
!      
!      ! Remove unconnected shelves
!      CALL determine_masks_ice(                mesh, ice)
!      CALL determine_masks_transitions(        mesh, ice)
!      CALL remove_unconnected_shelves(         mesh, ice)
!    END IF
    
    CALL update_general_ice_model_data(      mesh, ice)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE update_ice_thickness
  
! == Time stepping
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
    REAL(dp)                                           :: dist, dt
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
      dt = dist / (ABS( u_c( aci)) + ABS( v_c( aci)))
      dt_crit_adv = MIN( dt_crit_adv, dt)
       
    END DO
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = MIN(C%dt_max, dt_crit_adv * dt_correction_factor)
    
    ! Clean up after yourself
    CALL deallocate_shared( wu_c)
    CALL deallocate_shared( wv_c)
    
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
      
      IF     (C%choice_ice_dynamics == 'none') THEN
        ! Just stick to the maximum time step
      ELSEIF (C%choice_ice_dynamics == 'SIA') THEN
        t_next = MIN( t_next, region%t_next_SIA)
      ELSEIF (C%choice_ice_dynamics == 'SSA') THEN
        t_next = MIN( t_next, region%t_next_SSA)
      ELSEIF (C%choice_ice_dynamics == 'SIA/SSA') THEN
        t_next = MIN( t_next, region%t_next_SIA)
        t_next = MIN( t_next, region%t_next_SSA)
      ELSEIF (C%choice_ice_dynamics == 'DIVA') THEN
        t_next = MIN( t_next, region%t_next_DIVA)
      ELSE
        CALL crash('unknown choice_ice_dynamics "' // TRIM( C%choice_ice_dynamics) // '"!')
      END IF ! IF (C%choice_ice_dynamics == 'SIA') THEN
      
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
!      IF (C%do_asynchronous_BMB) THEN
        IF (region%time == region%t_next_BMB) THEN
          region%do_BMB         = .TRUE.
          region%t_last_BMB     = region%time
          region%t_next_BMB     = region%t_last_BMB + C%dt_BMB
        END IF
        t_next = MIN( t_next, region%t_next_BMB)
!      ELSE
!        ! Don't use separate timestepping for the BMB; just run it in every ice dynamics time step
!        region%do_BMB = .TRUE.
!      END IF
      
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
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE determine_timesteps_and_actions
  
! == Administration: allocation, initialisation, and remapping
  SUBROUTINE initialise_ice_model( mesh, ice, refgeo_init)
    ! Allocate shared memory for all the data fields of the ice dynamical module, and
    ! initialise some of them    
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_init
    
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
      ice%Hi_a( vi) = refgeo_init%Hi( vi)
      ice%Hb_a( vi) = refgeo_init%Hb( vi)
      ice%Hs_a( vi) = surface_elevation( ice%Hi_a( vi), ice%Hb_a( vi), 0._dp)
    END DO
    CALL sync
    
    ! Initialise masks and slopes
    CALL update_general_ice_model_data( mesh, ice)
    
    ! Allocate and initialise basal conditions
    CALL initialise_basal_conditions( mesh, ice)
    
    ! Geothermal heat flux
    IF     (C%choice_geothermal_heat_flux == 'constant') THEN
      ice%GHF_a( mesh%vi1:mesh%vi2) = C%constant_geothermal_heat_flux
    ELSEIF (C%choice_geothermal_heat_flux == 'spatial') THEN
      CALL map_geothermal_heat_flux_to_mesh( mesh, ice)
    ELSE
      CALL crash('unknown choice_geothermal_heat_flux "' // TRIM( C%choice_geothermal_heat_flux) // '"!')
    END IF
    
    ! Initialise data and matrices for the velocity solver(s)
    CALL initialise_velocity_solver( mesh, ice)
    
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
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%uabs_surf_a           , ice%wuabs_surf_a          )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%uabs_surf_b           , ice%wuabs_surf_b          )
    
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%u_base_a              , ice%wu_base_a             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%v_base_a              , ice%wv_base_a             )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%u_base_b              , ice%wu_base_b             )
    CALL allocate_shared_dp_1D(   mesh%nTri,              ice%v_base_b              , ice%wv_base_b             )
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
    CALL allocate_shared_dp_2D(   mesh%nV  , mesh%nC_mem, ice%dVi_in                , ice%wdVi_in               )
    CALL allocate_shared_dp_2D(   mesh%nV  , mesh%nC_mem, ice%dVi_out               , ice%wdVi_out              ) 
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%dHi_dt_a              , ice%wdHi_dt_a             )
    CALL allocate_shared_dp_1D(   mesh%nV  ,              ice%Hi_tplusdt_a          , ice%wHi_tplusdt_a         )
    
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
        
    ! The only fields that actually need to be mapped. The rest only needs memory reallocation.
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%Hi_a        , ice%wHi_a        , 'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dHi_dt_a    , ice%wdHi_dt_a    , 'cons_2nd_order')
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%Hi_tplusdt_a, ice%wHi_tplusdt_a, 'cons_2nd_order')
   !CALL remap_field_dp_3D( mesh_old, mesh_new, map, ice%Ti_a        , ice%wTi_a        , 'cons_2nd_order')
    
    ! Remove very thin ice resulting from remapping errors
    DO vi = mesh_new%vi1, mesh_new%vi2
      IF (ice%Hi_a( vi) < 0.1_dp) THEN
        ice%Hi_a( vi) = 0._dp
      END IF
    END DO
    CALL sync
    
    ! Remap englacial temperature (needs some special attention because of the discontinuity at the ice margin)
    CALL remap_ice_temperature( mesh_old, mesh_new, map, ice)
    
    ! Remap bedrock change and add up to PD to prevent accumulation of numerical diffusion
    CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%dHb_a,               ice%wdHb_a,               'cons_2nd_order')
    CALL reallocate_shared_dp_1D(    mesh_new%nV,  ice%Hb_a,                   ice%wHb_a                     )
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
    
    ! Simple memory reallocation for all the rest
    ! ===========================================
    
    ! Basic data - ice thickness, bedrock height, surface height, mask, 3D ice velocities and temperature
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hi_a                  , ice%wHi_a                 )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hb_a                  , ice%wHb_a                 )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hs_a                  , ice%wHs_a                 )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%SL_a                  , ice%wSL_a                 )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%TAF_a                 , ice%wTAF_a                )
   !CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%Ti_a                  , ice%wTi_a                 ) 
    
    ! Ice velocities (these are remapped/reallocated in remap_velocities)
   !CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%u_3D_a                , ice%wu_3D_a               )
   !CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%v_3D_a                , ice%wv_3D_a               )
   !CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%u_3D_b                , ice%wu_3D_b               )
   !CALL reallocate_shared_dp_2D(   mesh_new%nTri, C%nz           , ice%v_3D_b                , ice%wv_3D_b               )
   !CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%w_3D_a                , ice%ww_3D_a               )
    
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_vav_a               , ice%wu_vav_a              )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_vav_a               , ice%wv_vav_a              )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_vav_b               , ice%wu_vav_b              )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_vav_b               , ice%wv_vav_b              )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_vav_a            , ice%wuabs_vav_a           )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_vav_b            , ice%wuabs_vav_b           )
    
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_surf_a              , ice%wu_surf_a             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_surf_a              , ice%wv_surf_a             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_surf_b              , ice%wu_surf_b             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_surf_b              , ice%wv_surf_b             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_surf_a           , ice%wuabs_surf_a          )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_surf_b           , ice%wuabs_surf_b          )
    
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%u_base_a              , ice%wu_base_a             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%v_base_a              , ice%wv_base_a             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%u_base_b              , ice%wu_base_b             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%v_base_b              , ice%wv_base_b             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%uabs_base_a           , ice%wuabs_base_a          )
   !CALL reallocate_shared_dp_1D(   mesh_new%nTri,                  ice%uabs_base_b           , ice%wuabs_base_b          )
    
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
    
    ! Zeta derivatives
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dzeta_dt_a            , ice%wdzeta_dt_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dzeta_dx_a            , ice%wdzeta_dx_a           )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , C%nz           , ice%dzeta_dy_a            , ice%wdzeta_dy_a           )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dzeta_dz_a            , ice%wdzeta_dz_a           )
    
    ! Ice dynamics - ice thickness calculation
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , mesh_new%nC_mem, ice%dVi_in                , ice%wdVi_in               )
    CALL reallocate_shared_dp_2D(   mesh_new%nV  , mesh_new%nC_mem, ice%dVi_out               , ice%wdVi_out              ) 
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dHi_dt_a              , ice%wdHi_dt_a             )
   !CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%Hi_tplusdt_a          , ice%wHi_tplusdt_a         )
    
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
    
    ! GIA
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dHb_a                 , ice%wdHb_a                )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dHb_dt_a              , ice%wdHb_dt_a             )
    CALL reallocate_shared_dp_1D(   mesh_new%nV  ,                  ice%dSL_dt_a              , ice%wdSL_dt_a             )
    
    ! End of memory reallocation
    ! ==========================
    
    ! Remap velocities
    CALL remap_velocities( mesh_old, mesh_new, map, ice)
    
    ! Recalculate ice flow factor
    CALL calc_ice_rheology( mesh_new, ice, time)
    
    ! Update masks
    CALL update_general_ice_model_data( mesh_new, ice)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE remap_ice_model
  
END MODULE ice_dynamics_module
