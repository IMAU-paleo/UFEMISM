MODULE basal_conditions_and_sliding_module

  ! Contains all the routines for calculating the basal conditions underneath the ice.

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
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_1D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             SSA_Schoof2006_analytical_solution
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_remapping_mesh_mesh, &
                                             type_reference_geometry, type_grid, type_restart_data
  USE mesh_mapping_module,             ONLY: remap_field_dp_2D, smooth_Gaussian_2D

  IMPLICIT NONE

CONTAINS

  ! The main routine, to be called from the ice_velocity_module
  SUBROUTINE calc_basal_conditions( mesh, ice)
    ! Determine the basal conditions underneath the ice

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_conditions'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Basal hydrology
    CALL calc_basal_hydrology( mesh, ice)

    ! Bed roughness
    CALL calc_bed_roughness( mesh, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_conditions
  SUBROUTINE initialise_basal_conditions( mesh, ice, restart)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_restart_data),             INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_conditions'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Basal hydrology
    CALL initialise_basal_hydrology( mesh, ice)

    ! Bed roughness
    CALL initialise_bed_roughness( mesh, ice, restart)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = HUGE( 1))

  END SUBROUTINE initialise_basal_conditions

! == Basal hydrology
! ==================

  SUBROUTINE calc_basal_hydrology( mesh, ice)
    ! Calculate the pore water pressure and effective basal pressure

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_basal_hydrology'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate pore water pressure using the chosen basal hydrology model
    ! ====================================================================

    IF     (C%choice_basal_hydrology == 'saturated') THEN
      ! Assume all marine till is saturated (i.e. pore water pressure is equal to water pressure at depth everywhere)
      CALL calc_pore_water_pressure_saturated( mesh, ice)
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      ! The Martin et al. (2011) parameterisation of pore water pressure
      CALL calc_pore_water_pressure_Martin2011( mesh, ice)
    ELSE
      CALL crash('unknown choice_basal_hydrology "' // TRIM( C%choice_basal_hydrology) // '"!')
    END IF

    ! Calculate overburden and effective pressure
    ! ===========================================

    DO vi = mesh%vi1, mesh%vi2
      ice%overburden_pressure_a( vi) = ice_density * grav * ice%Hi_a( vi)
      ice%Neff_a(                vi) = MAX(0._dp, ice%overburden_pressure_a( vi) - ice%pore_water_pressure_a( vi))
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_hydrology
  SUBROUTINE initialise_basal_hydrology( mesh, ice)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_hydrology'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    IF     (C%choice_basal_hydrology == 'saturated') THEN
      CALL allocate_shared_dp_1D( mesh%nV, ice%pore_water_pressure_a, ice%wpore_water_pressure_a)
      CALL allocate_shared_dp_1D( mesh%nV, ice%overburden_pressure_a, ice%woverburden_pressure_a)
      CALL allocate_shared_dp_1D( mesh%nV, ice%Neff_a               , ice%wNeff_a               )
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      CALL allocate_shared_dp_1D( mesh%nV, ice%pore_water_pressure_a, ice%wpore_water_pressure_a)
      CALL allocate_shared_dp_1D( mesh%nV, ice%overburden_pressure_a, ice%woverburden_pressure_a)
      CALL allocate_shared_dp_1D( mesh%nV, ice%Neff_a               , ice%wNeff_a               )
    ELSE
      CALL crash('unknown choice_basal_hydrology "' // TRIM( C%choice_basal_hydrology) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 3)

  END SUBROUTINE initialise_basal_hydrology

  SUBROUTINE calc_pore_water_pressure_saturated( mesh, ice)
    ! Calculate the pore water pressure
    !
    ! Assume all till is saturated, i.e. pore water pressure = -rho_w * g * Hb

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_saturated'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      ice%pore_water_pressure_a( vi) = -seawater_density * grav * ice%Hb_a( vi)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_saturated
  SUBROUTINE calc_pore_water_pressure_Martin2011( mesh, ice)
    ! Calculate the pore water pressure
    !
    ! Use the parameterisation from Martin et al. (2011)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_pore_water_pressure_Martin2011'
    INTEGER                                            :: vi
    REAL(dp)                                           :: lambda_p

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Pore water pressure scaling factor (Martin et al., 2011, Eq. 12)
      lambda_p = MIN( 1._dp, MAX( 0._dp, 1._dp - (ice%Hb_a( vi) - ice%SL_a( vi) - C%Martin2011_hydro_Hb_min) / (C%Martin2011_hydro_Hb_max - C%Martin2011_hydro_Hb_min) ))

      ! Pore water pressure (Martin et al., 2011, Eq. 11)
      ice%pore_water_pressure_a( vi) = 0.96_dp * ice_density * grav * ice%Hi_a( vi) * lambda_p

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_pore_water_pressure_Martin2011

! == Bed roughness
! ================

  SUBROUTINE calc_bed_roughness( mesh, ice)
    ! Calculate the bed roughness

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! In case of no sliding or "idealised" sliding (e.g. ISMIP-HOM experiments), no bed roughness is required
    IF (C%choice_sliding_law == 'no_sliding' .OR. &
        C%choice_sliding_law == 'idealised') THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    IF (C%choice_basal_roughness == 'uniform') THEN
      ! Basal roughness values are constant and were already initialised; no need to do anything

    ELSEIF (C%choice_basal_roughness == 'parameterised') THEN
      ! Apply the chosen parameterisation of bed roughness

      IF     (C%choice_param_basal_roughness == 'none') THEN
        ! Nothing - apparently we're using an idealised sliding law where basal roughness is already included
      ELSEIF (C%choice_param_basal_roughness == 'Martin2011') THEN
        ! The Martin et al. (2011) parameterisation of basal roughness (specifically the till friction angle and till yield stress)
        CALL calc_bed_roughness_Martin2011( mesh, ice)
      ELSEIF (C%choice_param_basal_roughness == 'SSA_icestream') THEN
        ! The basal roughness parameterisation in the SSA_icestream idealised-geometry experiment
        CALL calc_bed_roughness_SSA_icestream( mesh, ice)
      ELSEIF (C%choice_param_basal_roughness == 'MISMIPplus') THEN
        ! The basal roughness parameterisation in the MISMIP+ idealised-geometry experiment
        CALL calc_bed_roughness_MISMIPplus( mesh, ice)
      ELSE
        CALL crash('unknown choice_param_basal_roughness "' // TRIM( C%choice_param_basal_roughness) // '"!')
      END IF

    ELSEIF (C%choice_basal_roughness == 'prescribed') THEN
      ! Basal roughness has been initialised from an external file; no need to do anything

    ELSEIF (C%choice_basal_roughness == 'restart') THEN
      ! Basal roughness has been initialised from a restart file; no need to do anything

    ELSE
      CALL crash('unknown choice_basal_roughness "' // TRIM( C%choice_basal_roughness) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness
  SUBROUTINE initialise_bed_roughness( mesh, ice, restart)
    ! Allocation and initialisation

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_restart_data),             INTENT(IN)    :: restart

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed

    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments

    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL allocate_shared_dp_1D( mesh%nV, ice%beta_sq_a , ice%wbeta_sq_a )

    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      CALL allocate_shared_dp_1D( mesh%nV, ice%phi_fric_a, ice%wphi_fric_a)
      CALL allocate_shared_dp_1D( mesh%nV, ice%tauc_a    , ice%wtauc_a    )

    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL allocate_shared_dp_1D( mesh%nV, ice%alpha_sq_a, ice%walpha_sq_a)
      CALL allocate_shared_dp_1D( mesh%nV, ice%beta_sq_a , ice%wbeta_sq_a )

    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL allocate_shared_dp_1D( mesh%nV, ice%alpha_sq_a, ice%walpha_sq_a)
      CALL allocate_shared_dp_1D( mesh%nV, ice%beta_sq_a , ice%wbeta_sq_a )

    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL allocate_shared_dp_1D( mesh%nV, ice%phi_fric_a, ice%wphi_fric_a)
      CALL allocate_shared_dp_1D( mesh%nV, ice%tauc_a    , ice%wtauc_a    )

    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF
    CALL sync

    ! Initialise values
    IF (C%choice_basal_roughness == 'uniform') THEN
      ! Apply a uniform bed roughness

      IF     (C%choice_sliding_law == 'Weertman') THEN
        ! Weertman sliding law; bed roughness is described by beta_sq
        ice%beta_sq_a( mesh%vi1:mesh%vi2) = C%slid_Weertman_beta_sq_uniform

      ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
        ! Coulomb sliding law; bed roughness is described by phi_fric
        ice%phi_fric_a( mesh%vi1:mesh%vi2) = C%slid_Coulomb_phi_fric_uniform

      ELSEIF (C%choice_sliding_law == 'Coulomb_regularised') THEN
        ! Regularised Coulomb sliding law; bed roughness is described by phi_fric
        ice%phi_fric_a( mesh%vi1:mesh%vi2) = C%slid_Coulomb_phi_fric_uniform

      ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
        ! Tsai2015 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part
        ice%alpha_sq_a( mesh%vi1:mesh%vi2) = C%slid_Tsai2015_alpha_sq_uniform
        ice%beta_sq_a(  mesh%vi1:mesh%vi2) = C%slid_Tsai2015_beta_sq_uniform

      ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
        ! Schoof2005 sliding law; bed roughness is described by alpha_sq for the Coulomb part, and beta_sq for the Weertman part
        ice%alpha_sq_a( mesh%vi1:mesh%vi2) = C%slid_Schoof2005_alpha_sq_uniform
        ice%beta_sq_a(  mesh%vi1:mesh%vi2) = C%slid_Schoof2005_beta_sq_uniform

      ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
        ! Zoet-Iverson sliding law; bed roughness is described by phi_fric
        ice%phi_fric_a( mesh%vi1:mesh%vi2) = C%slid_ZI_phi_fric_uniform

      ELSE
        CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
      END IF

    ELSEIF (C%choice_basal_roughness == 'parameterised') THEN
      ! Apply the chosen parameterisation of bed roughness
      CALL calc_bed_roughness( mesh, ice)

    ELSEIF (C%choice_basal_roughness == 'prescribed') THEN
      ! If bed roughness is prescribed, read it from the provided NetCDF file
      CALL initialise_bed_roughness_from_file( mesh, ice)

    ELSEIF (C%choice_basal_roughness == 'restart') THEN
      ! Assign the values that have been already read from a restart file

      IF (par%master) WRITE(0,*) '   Initialising bed roughness using data read from restart file...'

      IF     (C%choice_sliding_law == 'Weertman') THEN
        ice%beta_sq_a( mesh%vi1:mesh%vi2) = restart%beta_sq( mesh%vi1:mesh%vi2)

      ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
        ice%phi_fric_a( mesh%vi1:mesh%vi2) = restart%phi_fric( mesh%vi1:mesh%vi2)

      ELSEIF (C%choice_sliding_law == 'Coulomb_regularised') THEN
        ice%phi_fric_a( mesh%vi1:mesh%vi2) = restart%phi_fric( mesh%vi1:mesh%vi2)

      ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
        ice%alpha_sq_a( mesh%vi1:mesh%vi2) = C%slid_Tsai2015_alpha_sq_uniform
        ice%beta_sq_a(  mesh%vi1:mesh%vi2) = restart%beta_sq( mesh%vi1:mesh%vi2)

      ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
        ice%alpha_sq_a( mesh%vi1:mesh%vi2) = C%slid_Schoof2005_alpha_sq_uniform
        ice%beta_sq_a(  mesh%vi1:mesh%vi2) = restart%beta_sq( mesh%vi1:mesh%vi2)

      ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
        ice%phi_fric_a( mesh%vi1:mesh%vi2) = restart%phi_fric( mesh%vi1:mesh%vi2)

      ELSE
        CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
      END IF

    ELSE
      CALL crash('unknown choice_basal_roughness "' // TRIM( C%choice_basal_roughness) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = HUGE( 1))

  END SUBROUTINE initialise_bed_roughness

  ! The Martin et al. (2011) till parameterisation
  SUBROUTINE calc_bed_roughness_Martin2011( mesh, ice)
    ! Calculate the till friction angle phi_fric and till yield stress tauc,
    ! using the till model by Martin et al. (2011).
    !
    ! Only applicable when choice_sliding_law = "Coulomb" or "Coulomb_regularised"

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_Martin2011'
    INTEGER                                            :: vi
    REAL(dp)                                           :: w_Hb

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (.NOT. (C%choice_sliding_law == 'Coulomb' .OR. C%choice_sliding_law == 'Coulomb_regularised' .OR. C%choice_sliding_law == 'Zoet-Iverson')) THEN
      CALL crash('only applicable when choice_sliding_law = "Coulomb", "Coulomb_regularised", or "Zoet-Iverson"!')
    END IF

    DO vi = mesh%vi1, mesh%vi2

      ! Martin et al. (2011) Eq. 10
      w_Hb = MIN( 1._dp, MAX( 0._dp, (ice%Hb_a( vi) - C%Martin2011till_phi_Hb_min) / (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))
      ice%phi_fric_a( vi) = (1._dp - w_Hb) * C%Martin2011till_phi_min + w_Hb * C%Martin2011till_phi_max

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%phi_fric_a, 'ice%phi_fric_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_Martin2011

  ! Idealised cases
  SUBROUTINE calc_bed_roughness_SSA_icestream( mesh, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: SSA_icestream (i.e. the Schoof 2006 analytical solution)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_SSA_icestream'
    INTEGER                                            :: vi
    REAL(dp)                                           :: y, dummy1

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      y = mesh%V( vi,2)
      CALL SSA_Schoof2006_analytical_solution( 0.001_dp, ice%Hi_a( vi), ice%A_flow_vav_a( vi), y, dummy1, ice%tauc_a( vi))
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_SSA_icestream
  SUBROUTINE calc_bed_roughness_MISMIPplus( mesh, ice)
    ! Determine the basal conditions underneath the ice
    !
    ! Idealised case: MISMIP+ (see Asay-Davis et al., 2016)

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_bed_roughness_MISMIPplus'
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: vi
    REAL(dp), PARAMETER                                :: MISMIPplus_alpha_sq = 0.5_dp   ! Coulomb-law friction coefficient [unitless];         see Asay-Davis et al., 2016
    REAL(dp), PARAMETER                                :: MISMIPplus_beta_sq  = 1.0E4_dp ! Power-law friction coefficient   [Pa m^−1/3 yr^1/3]; idem dito

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_sliding_law == 'Weertman') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the first (Weertman) sliding law option

      DO vi = mesh%vi1, mesh%vi2
        ice%beta_sq_a( vi) = MISMIPplus_beta_sq
      END DO
      CALL sync

    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the second (Tsai et al., 2015) sliding law option

      DO vi = mesh%vi1, mesh%vi2
        ice%alpha_sq_a( vi) = MISMIPplus_alpha_sq
        ice%beta_sq_a(  vi) = MISMIPplus_beta_sq
      END DO
      CALL sync

    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Uniform sliding factor for the MISMIP+ configuration, using the third (Schoof, 2005) sliding law option

      DO vi = mesh%vi1, mesh%vi2
        ice%alpha_sq_a( vi) = MISMIPplus_alpha_sq
        ice%beta_sq_a(  vi) = MISMIPplus_beta_sq
      END DO
      CALL sync

    ELSE
      CALL crash('only defined when choice_sliding_law = "Weertman", "Tsai2015", or "Schoof2005"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_bed_roughness_MISMIPplus

  ! Initialise bed roughness from a file
  SUBROUTINE initialise_bed_roughness_from_file( mesh, ice)
    ! Initialise bed roughness with data from an external NetCDF file

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_sliding_law == 'no_sliding' .OR. &
            C%choice_sliding_law == 'idealised') THEN
      ! No sliding allowed / sliding laws for some idealised experiments
      CALL crash('not defined for choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL initialise_bed_roughness_from_file_Weertman( mesh, ice)
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Coulomb-type sliding law
      CALL initialise_bed_roughness_from_file_Coulomb( mesh, ice)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL initialise_bed_roughness_from_file_Tsai2015( mesh, ice)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL initialise_bed_roughness_from_file_Schoof2005( mesh, ice)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL initialise_bed_roughness_from_file_ZoetIverson( mesh, ice)
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file
  SUBROUTINE initialise_bed_roughness_from_file_Weertman( mesh, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Weertman-type sliding law: bed roughness described by beta_sq

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_Weertman'

    REAL(dp) :: dummy_dp
    dummy_dp = mesh%V( 1,1)
    dummy_dp = ice%Hi_a( 1)

    CALL crash('FIXME!')

!    ! Local variables:
!    TYPE(type_BIV_bed_roughness)                       :: BIV
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Determine filename
!    BIV%netcdf%filename = C%basal_roughness_filename
!
!    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
!
!    ! Inquire mesh data from the NetCDF file
!    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
!    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
!
!    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Allocate memory - mesh
!    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
!    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
!    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%beta_sq , BIV%wbeta_sq )
!
!    ! Read mesh & bed roughness data from file
!    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( BIV%beta_sq,  'BIV%beta_sq')
!
!    ! Since we want data represented as [j,i] internally, transpose the data we just read.
!    CALL transpose_dp_2D( BIV%beta_sq,  BIV%wbeta_sq )
!
!    ! Map (transposed) raw data to the model mesh
!    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%beta_sq , ice%beta_sq_a )
!
!    ! Deallocate raw data
!    CALL deallocate_shared( BIV%wnx      )
!    CALL deallocate_shared( BIV%wny      )
!    CALL deallocate_shared( BIV%wx       )
!    CALL deallocate_shared( BIV%wy       )
!    CALL deallocate_shared( BIV%wbeta_sq )
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_Weertman
  SUBROUTINE initialise_bed_roughness_from_file_Coulomb( mesh, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Coulomb-type sliding law: bed roughness described by phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_Coulomb'

    REAL(dp) :: dummy_dp
    dummy_dp = mesh%V( 1,1)
    dummy_dp = ice%Hi_a( 1)

    CALL crash('FIXME!')

!    ! Local variables:
!    TYPE(type_BIV_bed_roughness)                       :: BIV
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Determine filename
!    BIV%netcdf%filename = C%basal_roughness_filename
!
!    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
!
!    ! Inquire mesh data from the NetCDF file
!    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
!    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
!
!    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Allocate memory - mesh
!    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
!    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
!    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%phi_fric, BIV%wphi_fric)
!
!    ! Read mesh & bed roughness data from file
!    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( BIV%phi_fric, 'BIV%phi_fric')
!
!    ! Since we want data represented as [j,i] internally, transpose the data we just read.
!    CALL transpose_dp_2D( BIV%phi_fric, BIV%wphi_fric)
!
!    ! Map (transposed) raw data to the model mesh
!    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%phi_fric, ice%phi_fric_a)
!
!    ! Deallocate raw data
!    CALL deallocate_shared( BIV%wnx      )
!    CALL deallocate_shared( BIV%wny      )
!    CALL deallocate_shared( BIV%wx       )
!    CALL deallocate_shared( BIV%wy       )
!    CALL deallocate_shared( BIV%wphi_fric)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_Coulomb
  SUBROUTINE initialise_bed_roughness_from_file_Tsai2015( mesh, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Tsai 2015 sliding law: bed roughness described by alpha_sq & beta_sq

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_Tsai2015'

    REAL(dp) :: dummy_dp
    dummy_dp = mesh%V( 1,1)
    dummy_dp = ice%Hi_a( 1)

    CALL crash('FIXME!')

!    ! Local variables:
!    TYPE(type_BIV_bed_roughness)                       :: BIV
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Determine filename
!    BIV%netcdf%filename = C%basal_roughness_filename
!
!    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
!
!    ! Inquire mesh data from the NetCDF file
!    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
!    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
!
!    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Allocate memory - mesh
!    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
!    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
!    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%alpha_sq, BIV%walpha_sq)
!    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%beta_sq , BIV%wbeta_sq )
!
!    ! Read mesh & bed roughness data from file
!    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( BIV%alpha_sq, 'BIV%alpha_sq')
!    CALL check_for_NaN_dp_1D( BIV%beta_sq,  'BIV%beta_sq' )
!
!    ! Since we want data represented as [j,i] internally, transpose the data we just read.
!    CALL transpose_dp_2D( BIV%alpha_sq, BIV%walpha_sq)
!    CALL transpose_dp_2D( BIV%beta_sq,  BIV%wbeta_sq )
!
!    ! Map (transposed) raw data to the model mesh
!    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%alpha_sq, ice%alpha_sq_a)
!    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%beta_sq , ice%beta_sq_a )
!
!    ! Deallocate raw data
!    CALL deallocate_shared( BIV%wnx      )
!    CALL deallocate_shared( BIV%wny      )
!    CALL deallocate_shared( BIV%wx       )
!    CALL deallocate_shared( BIV%wy       )
!    CALL deallocate_shared( BIV%walpha_sq)
!    CALL deallocate_shared( BIV%wbeta_sq )
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_Tsai2015
  SUBROUTINE initialise_bed_roughness_from_file_Schoof2005( mesh, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Schoof 2005 sliding law: bed roughness described by alpha_sq & beta_sq

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_Schoof2005'

    REAL(dp) :: dummy_dp
    dummy_dp = mesh%V( 1,1)
    dummy_dp = ice%Hi_a( 1)

    CALL crash('FIXME!')

!    ! Local variables:
!    TYPE(type_BIV_bed_roughness)                       :: BIV
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Determine filename
!    BIV%netcdf%filename = C%basal_roughness_filename
!
!    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
!
!    ! Inquire mesh data from the NetCDF file
!    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
!    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
!
!    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Allocate memory - mesh
!    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
!    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
!    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%alpha_sq, BIV%walpha_sq)
!    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%beta_sq , BIV%wbeta_sq )
!
!    ! Read mesh & bed roughness data from file
!    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( BIV%alpha_sq, 'BIV%alpha_sq')
!    CALL check_for_NaN_dp_1D( BIV%beta_sq,  'BIV%beta_sq' )
!
!    ! Since we want data represented as [j,i] internally, transpose the data we just read.
!    CALL transpose_dp_2D( BIV%alpha_sq, BIV%walpha_sq)
!    CALL transpose_dp_2D( BIV%beta_sq,  BIV%wbeta_sq )
!
!    ! Map (transposed) raw data to the model mesh
!    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%alpha_sq, ice%alpha_sq_a)
!    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%beta_sq , ice%beta_sq_a )
!
!    ! Deallocate raw data
!    CALL deallocate_shared( BIV%wnx      )
!    CALL deallocate_shared( BIV%wny      )
!    CALL deallocate_shared( BIV%wx       )
!    CALL deallocate_shared( BIV%wy       )
!    CALL deallocate_shared( BIV%walpha_sq)
!    CALL deallocate_shared( BIV%wbeta_sq )
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_Schoof2005
  SUBROUTINE initialise_bed_roughness_from_file_ZoetIverson( mesh, ice)
    ! Initialise bed roughness with data from an external NetCDF file
    !
    ! Zoet-Iverson sliding law: bed roughness described by phi_fric

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_bed_roughness_from_file_ZoetIverson'

    REAL(dp) :: dummy_dp
    dummy_dp = mesh%V( 1,1)
    dummy_dp = ice%Hi_a( 1)

    CALL crash('FIXME!')

!    ! Local variables:
!    TYPE(type_BIV_bed_roughness)                       :: BIV
!
!    ! Add routine to path
!    CALL init_routine( routine_name)
!
!    ! Determine filename
!    BIV%netcdf%filename = C%basal_roughness_filename
!
!    IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'
!
!    ! Inquire mesh data from the NetCDF file
!    CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
!    CALL allocate_shared_int_0D( BIV%ny, BIV%wny)
!
!    IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Allocate memory - mesh
!    CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
!    CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
!    CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%phi_fric, BIV%wphi_fric)
!
!    ! Read mesh & bed roughness data from file
!    IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
!    CALL sync
!
!    ! Safety
!    CALL check_for_NaN_dp_1D( BIV%phi_fric, 'BIV%phi_fric')
!
!    ! Since we want data represented as [j,i] internally, transpose the data we just read.
!    CALL transpose_dp_2D( BIV%phi_fric, BIV%wphi_fric)
!
!    ! Map (transposed) raw data to the model mesh
!    CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%phi_fric, ice%phi_fric_a)
!
!    ! Deallocate raw data
!    CALL deallocate_shared( BIV%wnx      )
!    CALL deallocate_shared( BIV%wny      )
!    CALL deallocate_shared( BIV%wx       )
!    CALL deallocate_shared( BIV%wy       )
!    CALL deallocate_shared( BIV%wphi_fric)
!
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_ZoetIverson

! == Sliding laws
! ===============

  SUBROUTINE calc_sliding_law( mesh, ice, u_a, v_a, beta_a)
    ! Calculate the sliding term beta in the SSA/DIVA using the specified sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed (choice of beta is trivial)
      beta_a( mesh%vi1:mesh%vi2) = 0._dp
      CALL sync
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
      CALL calc_sliding_law_idealised(           mesh, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL calc_sliding_law_Weertman(            mesh, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
      ! Coulomb-type sliding law
      CALL calc_sliding_law_Coulomb(             mesh, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      CALL calc_sliding_law_Coulomb_regularised( mesh, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL calc_sliding_law_Tsai2015(            mesh, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL calc_sliding_law_Schoof2005(          mesh, ice, u_a, v_a, beta_a)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL calc_sliding_law_ZoetIverson(         mesh, ice, u_a, v_a, beta_a)
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law

  SUBROUTINE calc_sliding_law_Weertman( mesh, ice, u_a, v_a, beta_a)
    ! Weertman-type ("power law") sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Weertman'
    INTEGER                                            :: vi
    REAL(dp)                                           :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Asay-Davis et al. (2016), Eq. 6
      beta_a( vi) = ice%beta_sq_a( vi) * uabs ** (1._dp / C%slid_Weertman_m - 1._dp)

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Weertman
  SUBROUTINE calc_sliding_law_Coulomb( mesh, ice, u_a, v_a, beta_a)
    ! Coulomb-type sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Coulomb'
    INTEGER                                            :: vi
    REAL(dp)                                           :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the till yield stress from the till friction angle and the effective pressure
    DO vi = mesh%vi1, mesh%vi2
      ice%tauc_a( vi) = TAN((pi / 180._dp) * ice%phi_fric_a( vi)) * ice%Neff_a( vi)
    END DO
    CALL sync

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      beta_a( vi) = ice%tauc_a( vi) / uabs

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Coulomb
  SUBROUTINE calc_sliding_law_Coulomb_regularised( mesh, ice, u_a, v_a, beta_a)
    ! Regularised Coulomb-type sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Coulomb_regularised'
    INTEGER                                            :: vi
    REAL(dp)                                           :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the till yield stress from the till friction angle and the effective pressure
    DO vi = mesh%vi1, mesh%vi2
      ice%tauc_a( vi) = TAN((pi / 180._dp) * ice%phi_fric_a( vi)) * ice%Neff_a( vi)
    END DO
    CALL sync

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      beta_a( vi) = ice%tauc_a( vi) * uabs ** (C%slid_Coulomb_reg_q_plastic - 1._dp) / (C%slid_Coulomb_reg_u_threshold ** C%slid_Coulomb_reg_q_plastic)

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Coulomb_regularised
  SUBROUTINE calc_sliding_law_Tsai2015(  mesh, ice, u_a, v_a, beta_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Dvis et al., 2016)
    !
    ! Asay-Dvis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    !
    ! Tsai et al.: Marine ice-sheet profiles and stability under Coulomb basal conditions,
    ! Journal of Glaciology 61, 205–215, doi:10.3189/2015JoG14J221, 2015.

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Tsai2015'
    INTEGER                                            :: vi
    REAL(dp)                                           :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Asay-Dvis et al. (2016), Eq. 7
      beta_a( vi) = MIN( ice%alpha_sq_a( vi) * ice%Neff_a( vi), ice%beta_sq_a( vi) * uabs ** (1._dp / C%slid_Weertman_m)) * uabs**(-1._dp)

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Tsai2015
  SUBROUTINE calc_sliding_law_Schoof2005(  mesh, ice, u_a, v_a, beta_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Dvis et al., 2016)
    !
    ! Asay-Dvis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    !
    ! Schoof: The effect of cvitation on glacier sliding, P. Roy. Soc. A-Math. Phy., 461, 609–627, doi:10.1098/rspa.2004.1350, 2005

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_Schoof2005'
    INTEGER                                            :: vi
    REAL(dp)                                           :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Asay-Dvis et al. (2016), Eq. 11
      beta_a( vi) = ((ice%beta_sq_a( vi) * uabs**(1._dp / C%slid_Weertman_m) * ice%alpha_sq_a( vi) * ice%Neff_a( vi)) / &
        ((ice%beta_sq_a( vi)**C%slid_Weertman_m * uabs + (ice%alpha_sq_a( vi) * ice%Neff_a( vi))**C%slid_Weertman_m)**(1._dp / C%slid_Weertman_m))) * uabs**(-1._dp)

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Schoof2005
  SUBROUTINE calc_sliding_law_ZoetIverson( mesh, ice, u_a, v_a, beta_a)
    ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_ZoetIverson'
    INTEGER                                            :: vi
    REAL(dp)                                           :: uabs

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate the till yield stress from the till friction angle and the effective pressure
    DO vi = mesh%vi1, mesh%vi2
      ice%tauc_a( vi) = TAN((pi / 180._dp) * ice%phi_fric_a( vi)) * ice%Neff_a( vi)
    END DO
    CALL sync

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Zoet & Iverson (2020), Eq. (3) (divided by u to give beta = tau_b / u)
      beta_a( vi) = ice%tauc_a( vi) * (uabs**(1._dp / C%slid_ZI_p - 1._dp)) * ((uabs + C%slid_ZI_ut)**(-1._dp / C%slid_ZI_p))

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_ZoetIverson

  SUBROUTINE calc_sliding_law_idealised(  mesh, ice, u_a, v_a, beta_a)
    ! Sliding laws for some idealised experiments

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised'
    REAL(dp) :: dummy_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings...
    dummy_dp = u_a( 1)
    dummy_dp = v_a( 1)

    IF     (C%choice_idealised_sliding_law == 'ISMIP_HOM_C') THEN
      ! ISMIP-HOM experiment C

      CALL calc_sliding_law_idealised_ISMIP_HOM_C( mesh, beta_a)

    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_D') THEN
      ! ISMIP-HOM experiment D

      CALL calc_sliding_law_idealised_ISMIP_HOM_D( mesh, beta_a)

    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_E') THEN
      ! ISMIP-HOM experiment E

      CALL crash('the Glacier Arolla experiment is not implemented in UFEMISM!')

    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_F') THEN
      ! ISMIP-HOM experiment F

      CALL calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice, beta_a)

    ELSE
      CALL crash('unknown choice_idealised_sliding_law "' // TRIM( C%choice_idealised_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C( mesh, beta_a)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment C

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_C'
    INTEGER                                            :: vi
    REAL(dp)                                           :: x,y

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      beta_a( vi) = 1000._dp + 1000._dp * SIN( 2._dp * pi * x / C%ISMIP_HOM_L) * SIN( 2._dp * pi * y / C%ISMIP_HOM_L)
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D( mesh, beta_a)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment D

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_D'
    INTEGER                                            :: vi
    REAL(dp)                                           :: x

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      beta_a( vi) = 1000._dp + 1000._dp * SIN( 2._dp * pi * x / C%ISMIP_HOM_L)
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D
  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice, beta_a)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment F

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_F'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      beta_a( vi) = (ice%A_flow_vav_a( vi) * 1000._dp)**(-1._dp)
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F

! == Remapping
! ============

  SUBROUTINE remap_basal_conditions( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_basal_conditions'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Basal hydrology
    CALL remap_basal_hydrology( mesh_old, mesh_new, map, ice)

    ! Bed roughness
    CALL remap_bed_roughness( mesh_old, mesh_new, map, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_basal_conditions
  SUBROUTINE remap_basal_hydrology( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_basal_hydrology'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV
    int_dummy = map%int_dummy

    ! Allocate shared memory
    IF     (C%choice_basal_hydrology == 'saturated') THEN
      CALL reallocate_shared_dp_1D( mesh_new%nV, ice%pore_water_pressure_a, ice%wpore_water_pressure_a)
      CALL reallocate_shared_dp_1D( mesh_new%nV, ice%overburden_pressure_a, ice%woverburden_pressure_a)
      CALL reallocate_shared_dp_1D( mesh_new%nV, ice%Neff_a               , ice%wNeff_a               )
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      CALL reallocate_shared_dp_1D( mesh_new%nV, ice%pore_water_pressure_a, ice%wpore_water_pressure_a)
      CALL reallocate_shared_dp_1D( mesh_new%nV, ice%overburden_pressure_a, ice%woverburden_pressure_a)
      CALL reallocate_shared_dp_1D( mesh_new%nV, ice%Neff_a               , ice%wNeff_a               )
    ELSE
      CALL crash('unknown choice_basal_hydrology "' // TRIM( C%choice_basal_hydrology) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_basal_hydrology
  SUBROUTINE remap_bed_roughness( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_bed_roughness'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV
    int_dummy = map%int_dummy

    ! Allocate shared memory
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      IF (C%do_basal_sliding_inversion) THEN
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%beta_sq_a, ice%wbeta_sq_a, 'cons_2nd_order')
      ELSE
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%beta_sq_a , ice%wbeta_sq_a )
      END IF
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      IF (C%do_basal_sliding_inversion) THEN
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%phi_fric_a, ice%wphi_fric_a, 'cons_2nd_order')
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%tauc_a,     ice%wtauc_a,     'cons_2nd_order')
      ELSE
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%phi_fric_a, ice%wphi_fric_a)
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%tauc_a    , ice%wtauc_a    )
      END IF
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      IF (C%do_basal_sliding_inversion) THEN
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%alpha_sq_a, ice%walpha_sq_a, 'cons_2nd_order')
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%beta_sq_a,  ice%wbeta_sq_a,  'cons_2nd_order')
      ELSE
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%alpha_sq_a, ice%walpha_sq_a)
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%beta_sq_a , ice%wbeta_sq_a )
      END IF
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      IF (C%do_basal_sliding_inversion) THEN
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%alpha_sq_a, ice%walpha_sq_a, 'cons_2nd_order')
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%beta_sq_a,  ice%wbeta_sq_a,  'cons_2nd_order')
      ELSE
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%alpha_sq_a, ice%walpha_sq_a)
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%beta_sq_a , ice%wbeta_sq_a )
      END IF
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      IF (C%do_basal_sliding_inversion) THEN
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%phi_fric_a, ice%wphi_fric_a, 'cons_2nd_order')
        CALL remap_field_dp_2D( mesh_old, mesh_new, map, ice%tauc_a,     ice%wtauc_a,     'cons_2nd_order')
      ELSE
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%phi_fric_a, ice%wphi_fric_a)
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%tauc_a    , ice%wtauc_a    )
      END IF
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! If bed roughness is prescribed (and we are not doing
    ! an inversion), read it from the provided NetCDF file
    IF (C%choice_basal_roughness == 'prescribed' .AND. (.NOT. C%do_basal_sliding_inversion)) THEN
      CALL initialise_bed_roughness_from_file( mesh_new, ice)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_bed_roughness

! == Inversion
! ============

  SUBROUTINE basal_sliding_inversion( mesh, grid, ice, refgeo)
    ! Iteratively invert for basal friction conditions under the grounded ice sheet

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'basal_sliding_inversion'
    INTEGER                                            :: vi
    REAL(dp)                                           :: h_scale, h_delta, new_val, w_smooth
    REAL(dp),           DIMENSION(SIZE(ice%Hi_a))      :: rough_smoothed

    ! Add routine to path
    CALL init_routine( routine_name)

    h_scale = 1.0_dp/C%basal_sliding_inv_scale

    DO vi = mesh%vi1, mesh%vi2

      h_delta = ice%Hi_a( vi) - refgeo%Hi( vi)

      ! Invert only where the reference/model is grounded ice
      IF ( ice%mask_sheet_a( vi) == 1 ) THEN

        h_delta = MAX(-1.5_dp, MIN(1.5_dp, h_delta * h_scale))

        ! Further adjust only where the previous value is not improving the result
        IF ( (h_delta > 0._dp .AND. ice%dHi_dt_a( vi) >= 0._dp) .OR. &
             (h_delta < 0._dp .AND. ice%dHi_dt_a( vi) <= 0._dp) ) THEN

          IF (C%choice_sliding_law == 'Weertman' .OR. &
              C%choice_sliding_law == 'Tsai2015' .OR. &
              C%choice_sliding_law == 'Schoof2005') THEN

            new_val = ice%beta_sq_a( vi)
            new_val = new_val * (10._dp ** h_delta)
            new_val = MIN(MAX(new_val, 1._dp), 10._dp)

            ice%beta_sq_a( vi) = new_val

          ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
                  C%choice_sliding_law == 'Coulomb_regularised' .OR. &
                  C%choice_sliding_law == 'Zoet-Iverson') THEN

            new_val = ice%phi_fric_a( vi)
            new_val = new_val * (10._dp ** (-h_delta))
            new_val = MIN(MAX(new_val, 0.01_dp), 45._dp)

            ice%phi_fric_a( vi) = new_val

          ELSE
            CALL crash('choice_sliding_law "' // TRIM( C%choice_sliding_law) // '" not compatible with basal sliding inversion!')
          END IF

        END IF ! else phi_fric_a does not change from previous time step

      END IF ! else the reference is not grounded ice sheet, so leave it alone

    END DO
    CALL sync

    ! Smooth the resulting field
    ! ==========================

    IF (C%do_basal_sliding_smoothing) THEN

      IF (C%choice_sliding_law == 'Weertman' .OR. &
          C%choice_sliding_law == 'Tsai2015' .OR. &
          C%choice_sliding_law == 'Schoof2005') THEN

        ! CALL smooth_Gaussian_2D( mesh, grid, ice%beta_sq_a, C%basal_sliding_inv_rsmooth)

      ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
              C%choice_sliding_law == 'Coulomb_regularised' .OR. &
              C%choice_sliding_law == 'Zoet-Iverson') THEN

        rough_smoothed( 1:mesh%nV) = ice%phi_fric_a( 1:mesh%nV)
        CALL sync

        CALL smooth_Gaussian_2D( mesh, grid, rough_smoothed, C%basal_sliding_inv_rsmooth)

        DO vi = mesh%vi1, mesh%vi2
            ice%phi_fric_a( vi) = (1._dp - C%basal_sliding_inv_wsmooth) * ice%phi_fric_a( vi) + C%basal_sliding_inv_wsmooth * rough_smoothed( vi)
        END DO
        CALL sync

      ELSE
        CALL crash('choice_sliding_law "' // TRIM( C%choice_sliding_law) // '" not compatible with basal sliding inversion!')
      END IF

    END IF
    CALL sync

    ! Extrapolate the resulting field
    ! ===============================

    ! WIP. Not needed yet really, as the smoothing leaks over the entire domain. Kinda extrapolating.
    ! CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE basal_sliding_inversion

END MODULE basal_conditions_and_sliding_module
