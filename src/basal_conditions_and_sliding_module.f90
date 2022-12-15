MODULE basal_conditions_and_sliding_module

  ! Contains all the routines for calculating the basal conditions underneath the ice.

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

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
                                             checksum_dp_1D, checksum_dp_2D, checksum_dp_3D, &
                                             SSA_Schoof2006_analytical_solution, extrapolate_Gaussian_floodfill_mesh
  USE netcdf_debug_module,             ONLY: debug
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_reference_geometry, type_grid, type_restart_data
  USE data_types_netcdf_module,        ONLY: type_netcdf_BIV_target_velocity
  USE mesh_mapping_module,             ONLY: smooth_Gaussian_2D, map_from_xy_grid_to_mesh_2D, remap_field_dp_2D
  USE mesh_help_functions_module,      ONLY: mesh_bilinear_dp, find_containing_vertex
  USE netcdf_basic_module,             ONLY: open_existing_netcdf_file_for_reading, close_netcdf_file
  USE netcdf_input_module,             ONLY: read_field_from_xy_file_2D
  USE netcdf_output_module,            ONLY: create_new_netcdf_file_for_writing, setup_mesh_in_netcdf_file, setup_xy_grid_in_netcdf_file, &
                                             add_field_mesh_dp_2D_notime, add_field_grid_dp_2D_notime, write_to_field_multiple_options_mesh_dp_2D_notime, &
                                             write_to_field_multiple_options_grid_dp_2D_notime
  USE mesh_operators_module,           ONLY: map_b_to_a_2D, map_a_to_b_2D

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

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
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
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

    ! Basal inversion
    IF (C%do_BIVgeo) THEN
      CALL initialise_basal_inversion( mesh, ice)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = HUGE( 1))

  END SUBROUTINE initialise_basal_conditions

! ===== Basal hydrology =====
! ===========================

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

! ===== Bed roughness =====
! =========================

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

    ! Shared memory allocation
    ! ========================

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

    ! Inversion-stuff allocation
    ! ==========================

    IF (C%do_BIVgeo) THEN
      IF (C%choice_sliding_law == 'Weertman' .OR. &
          C%choice_sliding_law == 'Tsai2015' .OR. &
          C%choice_sliding_law == 'Schoof2005') THEN

        CALL allocate_shared_dp_1D( mesh%nV, ice%beta_sq_inv_a , ice%wbeta_sq_inv_a )

      ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
              C%choice_sliding_law == 'Coulomb_regularised' .OR. &
              C%choice_sliding_law == 'Zoet-Iverson') THEN

        CALL allocate_shared_dp_1D( mesh%nV, ice%phi_fric_inv_a, ice%wphi_fric_inv_a)

      ELSE
        CALL crash('choice_sliding_law "' // TRIM( C%choice_sliding_law) // '" not compatible with basal sliding inversion!')
      END IF
    END IF

    ! Initialisation
    ! ==============

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

    ! Inversion initialisation
    ! ========================

    IF (C%do_BIVgeo) THEN
      IF (C%choice_sliding_law == 'Weertman' .OR. &
          C%choice_sliding_law == 'Tsai2015' .OR. &
          C%choice_sliding_law == 'Schoof2005') THEN

        ice%beta_sq_inv_a(  mesh%vi1:mesh%vi2) = ice%beta_sq_a(  mesh%vi1:mesh%vi2)

      ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
              C%choice_sliding_law == 'Coulomb_regularised' .OR. &
              C%choice_sliding_law == 'Zoet-Iverson') THEN

        ice%phi_fric_inv_a( mesh%vi1:mesh%vi2) = ice%phi_fric_a( mesh%vi1:mesh%vi2)

      ELSE
        CALL crash('choice_sliding_law "' // TRIM( C%choice_sliding_law) // '" not compatible with basal sliding inversion!')
      END IF
    END IF

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

   ! ! Local variables:
   ! TYPE(type_BIV_bed_roughness)                       :: BIV

   ! ! Add routine to path
   ! CALL init_routine( routine_name)

   ! ! Determine filename
   ! BIV%netcdf%filename = C%basal_roughness_filename

   ! IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'

   ! ! Inquire mesh data from the NetCDF file
   ! CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
   ! CALL allocate_shared_int_0D( BIV%ny, BIV%wny)

   ! IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Allocate memory - mesh
   ! CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
   ! CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
   ! CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%beta_sq , BIV%wbeta_sq )

   ! ! Read mesh & bed roughness data from file
   ! IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Safety
   ! CALL check_for_NaN_dp_1D( BIV%beta_sq,  'BIV%beta_sq')

   ! ! Since we want data represented as [j,i] internally, transpose the data we just read.
   ! CALL transpose_dp_2D( BIV%beta_sq,  BIV%wbeta_sq )

   ! ! Map (transposed) raw data to the model mesh
   ! CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%beta_sq , ice%beta_sq_a )

   ! ! Deallocate raw data
   ! CALL deallocate_shared( BIV%wnx      )
   ! CALL deallocate_shared( BIV%wny      )
   ! CALL deallocate_shared( BIV%wx       )
   ! CALL deallocate_shared( BIV%wy       )
   ! CALL deallocate_shared( BIV%wbeta_sq )

   ! ! Finalise routine path
   ! CALL finalise_routine( routine_name)

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

   ! ! Local variables:
   ! TYPE(type_BIV_bed_roughness)                       :: BIV

   ! ! Add routine to path
   ! CALL init_routine( routine_name)

   ! ! Determine filename
   ! BIV%netcdf%filename = C%basal_roughness_filename

   ! IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'

   ! ! Inquire mesh data from the NetCDF file
   ! CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
   ! CALL allocate_shared_int_0D( BIV%ny, BIV%wny)

   ! IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Allocate memory - mesh
   ! CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
   ! CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
   ! CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%phi_fric, BIV%wphi_fric)

   ! ! Read mesh & bed roughness data from file
   ! IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Safety
   ! CALL check_for_NaN_dp_1D( BIV%phi_fric, 'BIV%phi_fric')

   ! ! Since we want data represented as [j,i] internally, transpose the data we just read.
   ! CALL transpose_dp_2D( BIV%phi_fric, BIV%wphi_fric)

   ! ! Map (transposed) raw data to the model mesh
   ! CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%phi_fric, ice%phi_fric_a)

   ! ! Deallocate raw data
   ! CALL deallocate_shared( BIV%wnx      )
   ! CALL deallocate_shared( BIV%wny      )
   ! CALL deallocate_shared( BIV%wx       )
   ! CALL deallocate_shared( BIV%wy       )
   ! CALL deallocate_shared( BIV%wphi_fric)

   ! ! Finalise routine path
   ! CALL finalise_routine( routine_name)

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

   ! ! Local variables:
   ! TYPE(type_BIV_bed_roughness)                       :: BIV

   ! ! Add routine to path
   ! CALL init_routine( routine_name)

   ! ! Determine filename
   ! BIV%netcdf%filename = C%basal_roughness_filename

   ! IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'

   ! ! Inquire mesh data from the NetCDF file
   ! CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
   ! CALL allocate_shared_int_0D( BIV%ny, BIV%wny)

   ! IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Allocate memory - mesh
   ! CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
   ! CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
   ! CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%alpha_sq, BIV%walpha_sq)
   ! CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%beta_sq , BIV%wbeta_sq )

   ! ! Read mesh & bed roughness data from file
   ! IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Safety
   ! CALL check_for_NaN_dp_1D( BIV%alpha_sq, 'BIV%alpha_sq')
   ! CALL check_for_NaN_dp_1D( BIV%beta_sq,  'BIV%beta_sq' )

   ! ! Since we want data represented as [j,i] internally, transpose the data we just read.
   ! CALL transpose_dp_2D( BIV%alpha_sq, BIV%walpha_sq)
   ! CALL transpose_dp_2D( BIV%beta_sq,  BIV%wbeta_sq )

   ! ! Map (transposed) raw data to the model mesh
   ! CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%alpha_sq, ice%alpha_sq_a)
   ! CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%beta_sq , ice%beta_sq_a )

   ! ! Deallocate raw data
   ! CALL deallocate_shared( BIV%wnx      )
   ! CALL deallocate_shared( BIV%wny      )
   ! CALL deallocate_shared( BIV%wx       )
   ! CALL deallocate_shared( BIV%wy       )
   ! CALL deallocate_shared( BIV%walpha_sq)
   ! CALL deallocate_shared( BIV%wbeta_sq )

   ! ! Finalise routine path
   ! CALL finalise_routine( routine_name)

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

   ! ! Local variables:
   ! TYPE(type_BIV_bed_roughness)                       :: BIV

   ! ! Add routine to path
   ! CALL init_routine( routine_name)

   ! ! Determine filename
   ! BIV%netcdf%filename = C%basal_roughness_filename

   ! IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'

   ! ! Inquire mesh data from the NetCDF file
   ! CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
   ! CALL allocate_shared_int_0D( BIV%ny, BIV%wny)

   ! IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Allocate memory - mesh
   ! CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
   ! CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
   ! CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%alpha_sq, BIV%walpha_sq)
   ! CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%beta_sq , BIV%wbeta_sq )

   ! ! Read mesh & bed roughness data from file
   ! IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Safety
   ! CALL check_for_NaN_dp_1D( BIV%alpha_sq, 'BIV%alpha_sq')
   ! CALL check_for_NaN_dp_1D( BIV%beta_sq,  'BIV%beta_sq' )

   ! ! Since we want data represented as [j,i] internally, transpose the data we just read.
   ! CALL transpose_dp_2D( BIV%alpha_sq, BIV%walpha_sq)
   ! CALL transpose_dp_2D( BIV%beta_sq,  BIV%wbeta_sq )

   ! ! Map (transposed) raw data to the model mesh
   ! CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%alpha_sq, ice%alpha_sq_a)
   ! CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%beta_sq , ice%beta_sq_a )

   ! ! Deallocate raw data
   ! CALL deallocate_shared( BIV%wnx      )
   ! CALL deallocate_shared( BIV%wny      )
   ! CALL deallocate_shared( BIV%wx       )
   ! CALL deallocate_shared( BIV%wy       )
   ! CALL deallocate_shared( BIV%walpha_sq)
   ! CALL deallocate_shared( BIV%wbeta_sq )

   ! ! Finalise routine path
   ! CALL finalise_routine( routine_name)

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

   ! ! Local variables:
   ! TYPE(type_BIV_bed_roughness)                       :: BIV

   ! ! Add routine to path
   ! CALL init_routine( routine_name)

   ! ! Determine filename
   ! BIV%netcdf%filename = C%basal_roughness_filename

   ! IF (par%master) WRITE(0,*) '  Initialising basal roughness from file ', TRIM( BIV%netcdf%filename), '...'

   ! ! Inquire mesh data from the NetCDF file
   ! CALL allocate_shared_int_0D( BIV%nx, BIV%wnx)
   ! CALL allocate_shared_int_0D( BIV%ny, BIV%wny)

   ! IF (par%master) CALL inquire_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Allocate memory - mesh
   ! CALL allocate_shared_dp_1D( BIV%nx,         BIV%x       , BIV%wx       )
   ! CALL allocate_shared_dp_1D(         BIV%ny, BIV%y       , BIV%wy       )
   ! CALL allocate_shared_dp_2D( BIV%ny, BIV%nx, BIV%phi_fric, BIV%wphi_fric)

   ! ! Read mesh & bed roughness data from file
   ! IF (par%master) CALL read_BIV_bed_roughness_file( BIV)
   ! CALL sync

   ! ! Safety
   ! CALL check_for_NaN_dp_1D( BIV%phi_fric, 'BIV%phi_fric')

   ! ! Since we want data represented as [j,i] internally, transpose the data we just read.
   ! CALL transpose_dp_2D( BIV%phi_fric, BIV%wphi_fric)

   ! ! Map (transposed) raw data to the model mesh
   ! CALL map_square_to_square_cons_2nd_order_2D( BIV%nx, BIV%ny, BIV%x, BIV%y, mesh%nx, mesh%ny, mesh%x, mesh%y, BIV%phi_fric, ice%phi_fric_a)

   ! ! Deallocate raw data
   ! CALL deallocate_shared( BIV%wnx      )
   ! CALL deallocate_shared( BIV%wny      )
   ! CALL deallocate_shared( BIV%wx       )
   ! CALL deallocate_shared( BIV%wy       )
   ! CALL deallocate_shared( BIV%wphi_fric)

   ! ! Finalise routine path
   ! CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_bed_roughness_from_file_ZoetIverson

! ===== Sliding laws =====
! ========================

  SUBROUTINE calc_basal_friction_coefficient( mesh, ice, u_b, v_b)
    ! Calculate the basal friction coefficient betab using the specified sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law'
    REAL(dp), DIMENSION(:    ), POINTER                :: u_a
    REAL(dp), DIMENSION(:    ), POINTER                :: v_a
    INTEGER                                            :: wu_a, wv_a

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, u_a    , wu_a    )
    CALL allocate_shared_dp_1D( mesh%nV, v_a    , wv_a    )

    ! Map velocities to the a-grid
    CALL map_b_to_a_2D( mesh, u_b, u_a)
    CALL map_b_to_a_2D( mesh, v_b, v_a)

    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed (choice of beta is trivial)
      ice%beta_b_a( mesh%vi1:mesh%vi2) = 0._dp
      CALL sync
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
      CALL calc_sliding_law_idealised(           mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Weertman-type ("power law") sliding law
      CALL calc_sliding_law_Weertman(            mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb') THEN
      ! Coulomb-type sliding law
      CALL calc_sliding_law_Coulomb(             mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      CALL calc_sliding_law_Coulomb_regularised( mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
      ! Modified power-law relation according to Tsai et al. (2015)
      CALL calc_sliding_law_Tsai2015(            mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation according to Schoof (2005)
      CALL calc_sliding_law_Schoof2005(          mesh, ice, u_a, v_a)
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL calc_sliding_law_ZoetIverson(         mesh, ice, u_a, v_a)
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wu_a    )
    CALL deallocate_shared( wv_a    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_basal_friction_coefficient

  SUBROUTINE calc_sliding_law_Weertman( mesh, ice, u_a, v_a)
    ! Weertman-type ("power law") sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a

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
      ice%beta_b_a( vi) = ice%beta_sq_a( vi) * uabs ** (1._dp / C%slid_Weertman_m - 1._dp)

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Weertman

  SUBROUTINE calc_sliding_law_Coulomb( mesh, ice, u_a, v_a)
    ! Coulomb-type sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a

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

      ice%beta_b_a( vi) = ice%tauc_a( vi) / uabs

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Coulomb

  SUBROUTINE calc_sliding_law_Coulomb_regularised( mesh, ice, u_a, v_a)
    ! Regularised Coulomb-type sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a

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

      ice%beta_b_a( vi) = ice%tauc_a( vi) * uabs ** (C%slid_Coulomb_reg_q_plastic - 1._dp) / (C%slid_Coulomb_reg_u_threshold ** C%slid_Coulomb_reg_q_plastic)

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Coulomb_regularised

  SUBROUTINE calc_sliding_law_Tsai2015(  mesh, ice, u_a, v_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
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
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a

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
      ice%beta_b_a( vi) = MIN( ice%alpha_sq_a( vi) * ice%Neff_a( vi), ice%beta_sq_a( vi) * uabs ** (1._dp / C%slid_Weertman_m)) * uabs**(-1._dp)

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Tsai2015

  SUBROUTINE calc_sliding_law_Schoof2005(  mesh, ice, u_a, v_a)
    ! Modified power-law relation according to Tsai et al. (2015)
    ! (implementation based on equations provided by Asay-Davis et al., 2016)
    !
    ! Asay-Dvis et al.: Experimental design for three interrelated marine ice sheet and ocean model
    ! intercomparison projects: MISMIP v. 3 (MISMIP+), ISOMIP v. 2 (ISOMIP+) and MISOMIP v. 1 (MISOMIP1),
    ! Geoscientific Model Development 9, 2471-2497, 2016
    !
    ! Schoof: The effect of cvitation on glacier sliding, P. Roy. Soc. A-Math. Phy., 461, 609–627, doi:10.1098/rspa.2004.1350, 2005

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a

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
      ice%beta_b_a( vi) = ((ice%beta_sq_a( vi) * uabs**(1._dp / C%slid_Weertman_m) * ice%alpha_sq_a( vi) * ice%Neff_a( vi)) / &
        ((ice%beta_sq_a( vi)**C%slid_Weertman_m * uabs + (ice%alpha_sq_a( vi) * ice%Neff_a( vi))**C%slid_Weertman_m)**(1._dp / C%slid_Weertman_m))) * uabs**(-1._dp)

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Schoof2005

  SUBROUTINE calc_sliding_law_ZoetIverson( mesh, ice, u_a, v_a)
    ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a

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
      ice%beta_b_a( vi) = ice%tauc_a( vi) * (uabs**(1._dp / C%slid_ZI_p - 1._dp)) * ((uabs + C%slid_ZI_ut)**(-1._dp / C%slid_ZI_p))

    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_ZoetIverson

  SUBROUTINE calc_sliding_law_idealised(  mesh, ice, u_a, v_a)
    ! Sliding laws for some idealised experiments

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_a
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_a

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

      CALL calc_sliding_law_idealised_ISMIP_HOM_C( mesh, ice)

    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_D') THEN
      ! ISMIP-HOM experiment D

      CALL calc_sliding_law_idealised_ISMIP_HOM_D( mesh, ice)

    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_E') THEN
      ! ISMIP-HOM experiment E

      CALL crash('the Glacier Arolla experiment is not implemented in UFEMISM!')

    ELSEIF (C%choice_idealised_sliding_law == 'ISMIP_HOM_F') THEN
      ! ISMIP-HOM experiment F

      CALL calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice)

    ELSE
      CALL crash('unknown choice_idealised_sliding_law "' // TRIM( C%choice_idealised_sliding_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised

  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C( mesh, ice)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment C

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_C'
    INTEGER                                            :: vi
    REAL(dp)                                           :: x,y

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      ice%beta_b_a( vi) = 1000._dp + 1000._dp * SIN( 2._dp * pi * x / C%ISMIP_HOM_L) * SIN( 2._dp * pi * y / C%ISMIP_HOM_L)
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_C

  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D( mesh, ice)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment D

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_D'
    INTEGER                                            :: vi
    REAL(dp)                                           :: x

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      ice%beta_b_a( vi) = 1000._dp + 1000._dp * SIN( 2._dp * pi * x / C%ISMIP_HOM_L)
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_D

  SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F( mesh, ice)
    ! Sliding laws for some idealised experiments
    !
    ! ISMIP-HOM experiment F

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised_ISMIP_HOM_F'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2
      ice%beta_b_a( vi) = (ice%A_flow_vav_a( vi) * 1000._dp)**(-1._dp)
    END DO
    CALL sync

    ! Safety
    CALL check_for_NaN_dp_1D( ice%beta_b_a, 'beta_b_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised_ISMIP_HOM_F

! ===== Remapping =====
! =====================

  SUBROUTINE remap_basal_conditions( mesh_old, mesh_new, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_basal_conditions'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Basal hydrology
    CALL remap_basal_hydrology( mesh_old, mesh_new, ice)

    ! Bed roughness
    CALL remap_bed_roughness( mesh_old, mesh_new, ice)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_basal_conditions

  SUBROUTINE remap_basal_hydrology( mesh_old, mesh_new, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_basal_hydrology'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV

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

  SUBROUTINE remap_bed_roughness( mesh_old, mesh_new, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_bed_roughness'
    INTEGER                                            :: int_dummy

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV

    ! == Reallocate/Remap shared memory
    ! =================================

    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed

    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments

    ELSEIF (C%choice_sliding_law == 'Weertman') THEN
      ! Power-law sliding law
      IF (C%do_BIVgeo .OR. C%choice_basal_roughness == 'restart') THEN
        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%beta_sq_a, ice%wbeta_sq_a, 'cons_2nd_order')
        IF (C%do_BIVgeo) THEN
          CALL remap_field_dp_2D( mesh_old, mesh_new, ice%beta_sq_inv_a, ice%wbeta_sq_inv_a, 'cons_2nd_order')
        END IF
      ELSE
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%beta_sq_a , ice%wbeta_sq_a )
      END IF

    ELSEIF (C%choice_sliding_law == 'Tsai2015' .OR. &
            C%choice_sliding_law == 'Schoof2005') THEN
      ! Modified power-law relation
      IF (C%do_BIVgeo .OR. C%choice_basal_roughness == 'restart') THEN
        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%alpha_sq_a, ice%walpha_sq_a, 'cons_2nd_order')
        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%beta_sq_a,  ice%wbeta_sq_a,  'cons_2nd_order')
        IF (C%do_BIVgeo) THEN
          CALL remap_field_dp_2D( mesh_old, mesh_new, ice%beta_sq_inv_a, ice%wbeta_sq_inv_a, 'cons_2nd_order')
        END IF
      ELSE
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%alpha_sq_a, ice%walpha_sq_a)
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%beta_sq_a , ice%wbeta_sq_a )
      END IF

    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised' .OR. &
            C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Yield-stress sliding law
      IF (C%do_BIVgeo .OR. C%choice_basal_roughness == 'restart') THEN
        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%phi_fric_a, ice%wphi_fric_a, 'cons_2nd_order')
        CALL remap_field_dp_2D( mesh_old, mesh_new, ice%tauc_a,     ice%wtauc_a,     'cons_2nd_order')
        IF (C%do_BIVgeo) THEN
          CALL remap_field_dp_2D( mesh_old, mesh_new, ice%phi_fric_inv_a, ice%wphi_fric_inv_a, 'cons_2nd_order')
        END IF
      ELSE
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%phi_fric_a, ice%wphi_fric_a)
        CALL reallocate_shared_dp_1D( mesh_new%nV, ice%tauc_a    , ice%wtauc_a    )
      END IF

    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! == Reinitialise values
    ! ======================

    IF (.NOT. C%do_BIVgeo) THEN
      ! Do not reset the values if we are doing an inversion (values from previous
      ! time-step and mesh remapped in the previous step above)

      IF (C%choice_basal_roughness == 'uniform') THEN
        ! Apply a uniform bed roughness

        IF     (C%choice_sliding_law == 'Weertman') THEN
          ice%beta_sq_a( mesh_new%vi1:mesh_new%vi2) = C%slid_Weertman_beta_sq_uniform
        ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
                C%choice_sliding_law == 'Coulomb_regularised') THEN
          ice%phi_fric_a( mesh_new%vi1:mesh_new%vi2) = C%slid_Coulomb_phi_fric_uniform
        ELSEIF (C%choice_sliding_law == 'Tsai2015') THEN
          ice%alpha_sq_a( mesh_new%vi1:mesh_new%vi2) = C%slid_Tsai2015_alpha_sq_uniform
          ice%beta_sq_a(  mesh_new%vi1:mesh_new%vi2) = C%slid_Tsai2015_beta_sq_uniform
        ELSEIF (C%choice_sliding_law == 'Schoof2005') THEN
          ice%alpha_sq_a( mesh_new%vi1:mesh_new%vi2) = C%slid_Schoof2005_alpha_sq_uniform
          ice%beta_sq_a(  mesh_new%vi1:mesh_new%vi2) = C%slid_Schoof2005_beta_sq_uniform
        ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
          ice%phi_fric_a( mesh_new%vi1:mesh_new%vi2) = C%slid_ZI_phi_fric_uniform
        ELSE
          CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
        END IF ! (C%choice_sliding_law)

      ELSEIF (C%choice_basal_roughness == 'parameterised') THEN
        ! Apply the chosen parameterisation of bed roughness
        CALL calc_bed_roughness( mesh_new, ice)

      ELSEIF (C%choice_basal_roughness == 'prescribed') THEN
        ! If bed roughness is prescribed, read it from the provided NetCDF file
        CALL initialise_bed_roughness_from_file( mesh_new, ice)

      ELSEIF (C%choice_basal_roughness == 'restart') THEN
        ! Do nothing, as these values were already remapped in the previous step above

      ELSE
        CALL crash('unknown choice_basal_roughness "' // TRIM( C%choice_basal_roughness) // '"!')

      END IF ! (C%choice_basal_roughness)
      CALL sync

    END IF ! (.NOT. C%do_BIVgeo)

    ! Basal inversion target velocity
    ! ===============================

    IF (C%do_BIVgeo) THEN
      IF (C%choice_BIVgeo_method == 'Berends2022') THEN
        CALL deallocate_shared( ice%wBIV_uabs_surf_target)
        CALL initialise_basal_inversion_target_velocity( mesh_new, ice)
      END IF
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_bed_roughness

! ===== Inversion =====
! =====================

  SUBROUTINE basal_sliding_inversion( mesh, grid, ice, refgeo, dt)
    ! Iteratively invert for basal friction conditions under the grounded ice sheet,
    ! and extrapolate the resulting field over the rest of the domain

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'basal_sliding_inversion'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Apply the selected inversion scheme
    IF     (C%choice_BIVgeo_method == 'Bernales2017') THEN
      CALL basal_sliding_inversion_Bernales2017( mesh, grid, ice, refgeo)
    ELSEIF (C%choice_BIVgeo_method == 'Berends2022') THEN
      CALL basal_sliding_inversion_Berends2022( mesh, grid, ice, refgeo, dt)
    ELSE
      CALL crash('unknown choice_BIVgeo_method "' // TRIM( C%choice_BIVgeo_method) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE basal_sliding_inversion

  SUBROUTINE basal_sliding_inversion_Bernales2017( mesh, grid, ice, refgeo)
    ! Iteratively invert for basal friction conditions under the grounded ice sheet,
    ! and extrapolate the resulting field over the rest of the domain

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'basal_sliding_inversion_Bernales2017'
    INTEGER                                            :: vi
    REAL(dp)                                           :: h_scale, h_delta, h_dfrac, new_val, w_smooth
    REAL(dp), DIMENSION(SIZE(ice%Hi_a))                :: rough_smoothed
    INTEGER,  DIMENSION(:    ), POINTER                ::  mask,  mask_filled
    INTEGER                                            :: wmask, wmask_filled

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate masks for extrapolation
    CALL allocate_shared_int_1D( mesh%nV, mask,        wmask       )
    CALL allocate_shared_int_1D( mesh%nV, mask_filled, wmask_filled)

    ! === Inversion scaling ===
    ! =========================

    ! Define the ice thickness factor for scaling of inversion
    h_scale = 1.0_dp/C%BIVgeo_Bernales_scale

    ! === Inversion ===
    ! =================

    DO vi = mesh%vi1, mesh%vi2

      ! Ice thickness difference w.r.t. reference thickness
      h_delta = ice%Hi_a( vi) - refgeo%Hi( vi)
      ! Ratio between this difference and the reference ice thickness
      h_dfrac = h_delta / MAX(refgeo%Hi( vi), 1._dp)

      ! Invert only where the model has grounded ice
      IF (ice%mask_sheet_a( vi) == 1) THEN

        ! Check if grounded vertex is marginal or interior ice.
        ! If marginal, override its value during extrapolation.
        IF (ice%mask_margin_a( vi) == 1 .OR. ice%mask_gl_a( vi) == 1) THEN
          ! Mark this vertex as ice margin/grounding-line
          mask( vi) = 1
        ELSE
          ! Mark this vertex as grounded ice
          mask( vi) = 2
        END IF

        ! If the difference/fraction is outside the specified tolerance
        IF (ABS(h_delta) >= C%BIVgeo_Bernales_tol_diff .OR. &
            ABS(h_dfrac) >= C%BIVgeo_Bernales_tol_frac) THEN

          ! Scale the difference and restrict it to the [-1.5 1.5] range
          h_delta = MAX(-1.5_dp, MIN(1.5_dp, h_delta * h_scale))

          ! Further adjust only where the previous value is not significantly improving the result
          IF ( (h_delta > 0._dp .AND. ice%dHi_dt_a( vi) >= -0.1_dp) .OR. &
               (h_delta < 0._dp .AND. ice%dHi_dt_a( vi) <=  0.1_dp) ) THEN

            ! Power-law-style sliding laws
            IF (C%choice_sliding_law == 'Weertman' .OR. &
                C%choice_sliding_law == 'Tsai2015' .OR. &
                C%choice_sliding_law == 'Schoof2005') THEN

              ! Save bed roughness at this vertex in aux variable
              new_val = ice%beta_sq_inv_a( vi)
              ! Adjust based on scaled ice thickness difference
              new_val = new_val * (10._dp ** h_delta)
              ! Constrain adjusted value to roughness limits
              new_val = MIN(MAX(new_val, 1._dp), 10._dp)
              ! Replace old bed roughness value with the adjusted one
              ice%beta_sq_inv_a( vi) = new_val

            ! Coulomb-style sliding laws
            ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
                    C%choice_sliding_law == 'Coulomb_regularised' .OR. &
                    C%choice_sliding_law == 'Zoet-Iverson') THEN

              ! Save bed roughness at this vertex in aux variable
              new_val = ice%phi_fric_inv_a( vi)
              ! Adjust based on scaled ice thickness difference
              new_val = new_val * (10._dp ** (-h_delta))
              ! Constrain adjusted value to roughness limits
              new_val = MIN(MAX(new_val, C%BIVgeo_Bernales_phi_min), C%BIVgeo_Bernales_phi_max)
              ! Replace old bed roughness value with the adjusted one
              ice%phi_fric_inv_a( vi) = new_val

            ELSE
              CALL crash('choice_sliding_law "' // TRIM( C%choice_sliding_law) // '" not compatible with basal sliding inversion!')
            END IF

          END IF ! else the fit is already improving for some other reason, so leave it alone

        END IF ! else the difference is within the specified tolerance, so leave it alone

      ELSE

        ! This vertex is not grounded ice sheet, so mark it for later extrapolation
        mask( vi) = 1

      END IF ! (ice%mask_sheet_a( vi) == 1)

    END DO
    CALL sync

    ! === Smoothing ===
    ! =================

    IF (C%BIVgeo_Bernales_do_smooth) THEN
      ! Smooth the resulting field

      IF (C%choice_sliding_law == 'Weertman' .OR. &
          C%choice_sliding_law == 'Tsai2015' .OR. &
          C%choice_sliding_law == 'Schoof2005') THEN

        ! Store the inverted parameters in a local variable
        rough_smoothed( 1:mesh%nV) = ice%beta_sq_inv_a( 1:mesh%nV)
        CALL sync

        ! Smooth the local variable
        CALL smooth_Gaussian_2D( mesh, grid, rough_smoothed, C%BIVgeo_Bernales_rsmooth)

        ! Combined the smoothed and raw inverted parameter through a weighed average
        DO vi = mesh%vi1, mesh%vi2
            ice%beta_sq_a( vi) = (1._dp - C%BIVgeo_Bernales_wsmooth) * ice%beta_sq_inv_a( vi) + C%BIVgeo_Bernales_wsmooth * rough_smoothed( vi)
            ! Make sure the variable stays within the prescribed limits
            ice%beta_sq_a( vi) = MIN(MAX(ice%beta_sq_a( vi), 1._dp), 10._dp)
        END DO
        CALL sync

      ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
              C%choice_sliding_law == 'Coulomb_regularised' .OR. &
              C%choice_sliding_law == 'Zoet-Iverson') THEN

        ! Store the inverted parameters in a local variable
        rough_smoothed( 1:mesh%nV) = ice%phi_fric_inv_a( 1:mesh%nV)
        CALL sync

        ! Smooth the local variable
        CALL smooth_Gaussian_2D( mesh, grid, rough_smoothed, C%BIVgeo_Bernales_rsmooth)

        ! Combined the smoothed and raw inverted parameter through a weighed average
        DO vi = mesh%vi1, mesh%vi2
            ice%phi_fric_a( vi) = (1._dp - C%BIVgeo_Bernales_wsmooth) * ice%phi_fric_inv_a( vi) + C%BIVgeo_Bernales_wsmooth * rough_smoothed( vi)
            ! Make sure the variable stays within the prescribed limits
            ice%phi_fric_a( vi) = MIN(MAX(ice%phi_fric_a( vi), C%BIVgeo_Bernales_phi_min), C%BIVgeo_Bernales_phi_max)
        END DO
        CALL sync

      ELSE
        CALL crash('choice_sliding_law "' // TRIM( C%choice_sliding_law) // '" not compatible with basal sliding inversion!')
      END IF

    ELSE
      ! Don't smooth the resulting field; simply copy it into the main variable

      IF (C%choice_sliding_law == 'Weertman' .OR. &
          C%choice_sliding_law == 'Tsai2015' .OR. &
          C%choice_sliding_law == 'Schoof2005') THEN

        ! Replace old bed roughness field with the adjusted one
        ice%beta_sq_a( mesh%vi1:mesh%vi2) = ice%beta_sq_inv_a( mesh%vi1:mesh%vi2)

      ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
              C%choice_sliding_law == 'Coulomb_regularised' .OR. &
              C%choice_sliding_law == 'Zoet-Iverson') THEN

        ! Replace old bed roughness field with the adjusted one
        ice%phi_fric_a( mesh%vi1:mesh%vi2) = ice%phi_fric_inv_a( mesh%vi1:mesh%vi2)

      ELSE
        CALL crash('choice_sliding_law "' // TRIM( C%choice_sliding_law) // '" not compatible with basal sliding inversion!')
      END IF

    END IF ! (C%BIVgeo_Bernales_do_smooth)
    CALL sync

    ! === Extrapolation ===
    ! =====================

    ! Perform the extrapolation
    IF (par%master) THEN
      CALL extrapolate_Gaussian_floodfill_mesh( mesh, mask, ice%phi_fric_a, 40000._dp, mask_filled)
    END IF
    CALL sync

    ! === Finalisation ===

    IF (C%choice_sliding_law == 'Weertman' .OR. &
        C%choice_sliding_law == 'Tsai2015' .OR. &
        C%choice_sliding_law == 'Schoof2005') THEN

      ! Safety
      CALL check_for_NaN_dp_1D( ice%beta_sq_a, 'beta_sq_a')

    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised' .OR. &
            C%choice_sliding_law == 'Zoet-Iverson') THEN

      ! Safety
      CALL check_for_NaN_dp_1D( ice%phi_fric_a, 'phi_fric_a')

    END IF

    ! Clean up after yourself
    CALL deallocate_shared( wmask_filled)
    CALL deallocate_shared( wmask)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE basal_sliding_inversion_Bernales2017

  SUBROUTINE basal_sliding_inversion_Berends2022( mesh, grid, ice, refgeo, dt)
    ! Iteratively invert for basal friction conditions under the grounded ice sheet,
    ! and extrapolate the resulting field over the rest of the domain

    IMPLICIT NONE

    ! In/output variables
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'basal_sliding_inversion_Berends2022'
    INTEGER,  DIMENSION(:    ), POINTER                ::  mask,  mask_filled
    INTEGER                                            :: wmask, wmask_filled
    REAL(dp), DIMENSION(:    ), POINTER                ::  dphi_dt
    INTEGER                                            :: wdphi_dt
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE            :: trace_up, trace_down
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: w_up, w_down
    INTEGER                                            :: vi,ti,n_up,n_down,k
    REAL(dp), DIMENSION(2)                             :: p,pt
    REAL(dp)                                           :: Hs_mod, Hs_target, u_mod, u_target
    REAL(dp)                                           :: I1, I2, I3, I_tot
    REAL(dp)                                           :: R
    REAL(dp)                                           :: h_delta, h_dfrac
    REAL(dp)                                           :: sigma

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_int_1D( mesh%nV, mask       , wmask       )
    CALL allocate_shared_int_1D( mesh%nV, mask_filled, wmask_filled)
    CALL allocate_shared_dp_1D(  mesh%nV, dphi_dt    , wdphi_dt    )

    ! Allocate memory for the up- and downstream traces
    ALLOCATE( trace_up(   mesh%nV, 2))
    ALLOCATE( trace_down( mesh%nV, 2))

    ! Allocate memory for the linear scaling functions
    ALLOCATE( w_up(   mesh%nV))
    ALLOCATE( w_down( mesh%nV))

    DO vi = mesh%vi1, mesh%vi2

      ! Obviously can only be done where there's (grounded) ice
      IF (ice%mask_sheet_a( vi) == 0 .OR. ice%mask_margin_a( vi) == 1 .OR. ice%Hi_a( vi) < 1._dp .OR. &
          refgeo%Hi( vi) < 1._dp .OR. ice%BIV_uabs_surf_target( vi) == 0._dp .OR. ice%mask_gl_a( vi) == 1) THEN
        mask( vi) = 1
        CYCLE
      ELSE
        mask( vi) = 2
      END IF

      ! The point p
      p = [mesh%V( vi,1), mesh%V( vi,2)]

      ! Trace the flowline upstream and downstream (Berends et al., 2022, Eqs. 2)
      ti = mesh%iTri( vi,1)
      CALL trace_flowline_upstream(   mesh, ice, p, trace_up  , n_up          , ti)
      CALL trace_flowline_downstream( mesh, ice, p, trace_down, n_down, refgeo, ti)

      ! Calculate linear scaling functions (Berends et al., Eqs. 5)
      w_up = 0._dp
      DO k = 1, n_up
        w_up( k) = REAL( n_up + 1 - k, dp)
      END DO
      w_up = w_up / SUM( w_up)

      w_down = 0._dp
      DO k = 1, n_down
        w_down( k) = REAL( n_down + 1 - k, dp)
      END DO
      w_down = w_down / SUM( w_down)

      ! Calculate upstream line integrals
      I1 = 0._dp
      I3 = 0._dp
      DO k = 1, n_up

        pt = trace_up( k,:)
        CALL mesh_bilinear_dp(  mesh, ice%uabs_surf_a         , pt, ti, u_mod)
        CALL mesh_bilinear_dp(  mesh, ice%BIV_uabs_surf_target, pt, ti, u_target)
        CALL mesh_bilinear_dp(  mesh, ice%Hs_a                , pt, ti, Hs_mod)
        CALL mesh_bilinear_dp(  mesh, refgeo%Hs               , pt, ti, Hs_target)

        ! If no target velocity data is available, assume zero difference
        IF (u_target /= u_target) u_target = u_mod

        I1 = I1 - ( u_mod -  u_target) * w_up( k) / C%BIVgeo_Berends2022_u0    ! Berends et al., (2022), Eq. 4a
        I3 = I3 + (Hs_mod - Hs_target) * w_up( k) / C%BIVgeo_Berends2022_H0    ! Berends et al., (2022), Eq. 4c

      END DO

      ! Calculate downstream line integral
      I2 = 0._dp
      DO k = 1, n_down

        pt = trace_down( k,:)
        CALL mesh_bilinear_dp(  mesh, ice%uabs_surf_a         , pt, ti, u_mod)
        CALL mesh_bilinear_dp(  mesh, ice%BIV_uabs_surf_target, pt, ti, u_target)

        ! If no target velocity data is available, assume zero difference
        IF (u_target /= u_target) u_target = u_mod

        I2 = I2 - ( u_mod -  u_target) * w_down( k) / C%BIVgeo_Berends2022_u0  ! Berends et al., (2022), Eq. 4b

      END DO

      ! Scale weights with local ice thickness * velocity
      ! (thinner and/or ice experiences less basal friction, so the solution is less sensitive to changes in bed roughness there)

      ! Berends et al. (2022), Eq. 7
      R = MAX( 0._dp, MIN( 1._dp, ((ice%uabs_surf_a( vi) * ice%Hi_a( vi)) / (C%BIVgeo_Berends2022_u_scale * C%BIVgeo_Berends2022_Hi_scale)) ))
      ! Berends et al. (2022), Eq. 6
      I_tot = (I1 + I2 + I3) * R

      ! Ice thickness difference w.r.t. reference thickness
      h_delta = ice%Hi_a( vi) - refgeo%Hi( vi)
      ! Ratio between this difference and the reference ice thickness
      h_dfrac = h_delta / MAX( refgeo%Hi( vi), 1._dp)

      ! If the difference/fraction is outside the specified tolerance
      IF (ABS( h_delta) >= C%BIVgeo_Bernales_tol_diff .OR. &
          ABS( h_dfrac) >= C%BIVgeo_Bernales_tol_frac) THEN

        ! Further adjust only where the previous value is not improving the result
        IF ( (h_delta > 0._dp .AND. ice%dHi_dt_a( vi) >= 0.0_dp) .OR. &
             (h_delta < 0._dp .AND. ice%dHi_dt_a( vi) <= 0.0_dp) ) THEN

          ! Calculate rate of change of bed roughness (Berends et al. (2022), Eq. 8)
          dphi_dt( vi) = -ice%phi_fric_a( vi) * I_tot / C%BIVgeo_Berends2022_tauc

        END IF
      END IF

    END DO
    CALL sync

    ! Extrapolate new values outside the ice sheet
    CALL extrapolate_Gaussian_floodfill_mesh( mesh, mask, dphi_dt, 40000._dp, mask_filled)

    ! Update bed roughness (Berends et al. (2022), Eq. 9)

    ! First regularisation step (function F1 in Berends et al. (2022), Eq. 9)
    sigma = mesh%resolution_min / 1.5_dp
    CALL smooth_Gaussian_2D( mesh, grid, dphi_dt, sigma)

    DO vi = mesh%vi1, mesh%vi2
      ice%phi_fric_a( vi) = MAX( C%BIVgeo_Berends2022_phimin, MIN( C%BIVgeo_Berends2022_phimax, ice%phi_fric_a( vi) + dphi_dt( vi) * dt ))
    END DO
    CALL sync

    ! Second regularisation step (function F2 in Berends et al. (2022), Eq. 9)
    sigma = mesh%resolution_min / 4._dp
    CALL smooth_Gaussian_2D( mesh, grid, ice%phi_fric_a, sigma)

    ! Clean up after yourself
    DEALLOCATE( trace_up  )
    DEALLOCATE( trace_down)
    DEALLOCATE( w_up      )
    DEALLOCATE( w_down    )
    CALL deallocate_shared( wmask       )
    CALL deallocate_shared( wmask_filled)
    CALL deallocate_shared( wdphi_dt    )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE basal_sliding_inversion_Berends2022

  SUBROUTINE trace_flowline_upstream( mesh, ice, p, T, n, ti)
    ! Trace the flowline passing through point p upstream.
    ! Returns a list T of n points on the flowline spaced dx_trace apart.
    !
    ! Stop the trace when it encounters the ice divide (defined as an ice velocity lower than 1 m/yr)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: T
    INTEGER,                             INTENT(OUT)   :: n
    INTEGER,                             INTENT(INOUT) :: ti

    ! Local variables:
    INTEGER                                            :: ii,jj,it,nmax
    REAL(dp), DIMENSION(2)                             :: pt, u
    REAL(dp)                                           :: dx_trace

    nmax = SIZE( T,1)

    dx_trace = mesh%resolution_min

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p
    it = 0

    DO WHILE (n < nmax)

      ! Interpolate surface velocity to the tracer location
      CALL mesh_bilinear_dp(  mesh, ice%u_surf_a, pt, ti, u( 1))
      CALL mesh_bilinear_dp(  mesh, ice%v_surf_a, pt, ti, u( 2))

      ! If we've reached the ice divide, end the trace
      IF (NORM2( u) < 1._dp) EXIT

      ! Add current position to the traces
      n = n + 1
      T( n,:) = pt

      ! Normalise velocity vector
      u = u / NORM2( u)

      ! Move the tracer upstream
      pt = pt - u * dx_trace

      ! If the new tracer location is outside the domain, end the trace
      IF (pt( 1) <= mesh%xmin .OR. pt( 2) >= mesh%xmax .OR. &
          pt( 2) <= mesh%ymin .OR. pt( 2) >= mesh%ymax) EXIT

      ! Safety
      it = it + 1
      IF (it > nmax) EXIT

    END DO ! DO WHILE (n_up < MAX(grid%nx, grid%ny))

    ! Safety
    IF (n == 0) THEN
      n = 1
      T( 1,:) = p
    END IF

  END SUBROUTINE trace_flowline_upstream

  SUBROUTINE trace_flowline_downstream( mesh, ice, p, T, n, refgeo, ti)
    ! Trace the flowline passing through point p downstream.
    ! Returns a list T of n points on the flowline spaced dx_trace_rel * grid%dx apart.
    !
    ! Stop the trace when it encounters the ice margin (either in the modelled or the target geometry)

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(2),              INTENT(IN)    :: p
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: T
    INTEGER,                             INTENT(OUT)   :: n
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo
    INTEGER,                             INTENT(INOUT) :: ti

    ! Local variables:
    INTEGER                                            :: vi,it,nmax
    REAL(dp), DIMENSION(2)                             :: pt, u
    REAL(dp)                                           :: dx_trace

    nmax = SIZE( T,1)

    dx_trace = mesh%resolution_min

    ! Initialise
    T  = 0._dp
    n  = 0
    pt = p
    it = 0

    DO WHILE (n < nmax)

      ! Find current tracer location grid indices
      vi = mesh%Tri( ti,1)
      CALL find_containing_vertex( mesh, pt, vi)

      ! If we've reached the ice margin, end the trace
      IF (ice%mask_ice_a( vi) == 0 .OR. ice%mask_margin_a( vi) == 1 .OR. ice%Hi_a( vi) < 1._dp .OR. &
          refgeo%Hi( vi) < 1._dp .OR. ice%BIV_uabs_surf_target( vi) == 0._dp) EXIT

      ! Interpolate surface velocity to the tracer location
      CALL mesh_bilinear_dp(  mesh, ice%u_surf_a, pt, ti, u( 1))
      CALL mesh_bilinear_dp(  mesh, ice%v_surf_a, pt, ti, u( 2))

      ! Safety
      IF (u( 1) == 0._dp .AND. u( 2) == 0._dp) EXIT

      ! Add current position to the traces
      n = n + 1
      T( n,:) = pt

      ! Normalise velocity vector
      u = u / NORM2( u)

      ! Move the tracer downstream
      pt = pt + u * dx_trace

      ! If the new tracer location is outside the domain, end the trace
      IF (pt( 1) <= mesh%xmin .OR. pt( 2) >= mesh%xmax .OR. &
          pt( 2) <= mesh%ymin .OR. pt( 2) >= mesh%ymax) EXIT

      ! Safety
      it = it + 1
      IF (it > nmax) EXIT

    END DO ! DO WHILE (n_down < MAX(grid%nx, grid%ny))

    ! Safety
    IF (n == 0) THEN
      n = 1
      T( 1,:) = p
    END IF

  END SUBROUTINE trace_flowline_downstream

  SUBROUTINE write_inverted_bed_roughness_to_file( mesh, grid, ice, region_name)
    ! Create a new NetCDF file and write the inverted bed roughness to it

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'write_inverted_bed_roughness_to_file'
    CHARACTER(LEN=256)                                 :: filename_mesh, filename_grid
    INTEGER                                            :: i
    INTEGER                                            :: ncid_mesh, ncid_grid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! File names
    filename_mesh = TRIM(C%output_dir) // TRIM( C%BIVgeo_filename_output)
    filename_grid = TRIM(C%output_dir) // TRIM( C%BIVgeo_filename_output)
    i = INDEX( filename_grid,'.nc')
    filename_grid = filename_grid( 1:i-1) // '_grid.nc'

    ! Create files
    CALL create_new_netcdf_file_for_writing( filename_mesh, ncid_mesh)
    CALL create_new_netcdf_file_for_writing( filename_grid, ncid_grid)

    ! Set up mesh and grid
    CALL setup_mesh_in_netcdf_file(    filename_mesh, ncid_mesh, mesh)
    CALL setup_xy_grid_in_netcdf_file( filename_grid, ncid_grid, grid)

    ! Add variables and write data
    IF (C%choice_sliding_law == 'Weertman' .OR. &
        C%choice_sliding_law == 'Tsai2015' .OR. &
        C%choice_sliding_law == 'Schoof2005') THEN

      ! Add variables
      CALL add_field_mesh_dp_2D_notime( filename_mesh, ncid_mesh, 'beta_sq', long_name = 'beta_sq', units = 'Pa m^−1/3 yr^1/3')
      CALL add_field_grid_dp_2D_notime( filename_grid, ncid_grid, 'beta_sq', long_name = 'beta_sq', units = 'Pa m^−1/3 yr^1/3')
      ! Write data
      CALL write_to_field_multiple_options_mesh_dp_2D_notime( filename_mesh, ncid_mesh,                    'beta_sq', ice%beta_sq_a)
      CALL write_to_field_multiple_options_grid_dp_2D_notime( filename_grid, ncid_grid, mesh, region_name, 'beta_sq', ice%beta_sq_a)

    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised' .OR. &
            C%choice_sliding_law == 'Zoet-Iverson') THEN

      ! Add variables
      CALL add_field_mesh_dp_2D_notime( filename_mesh, ncid_mesh, 'phi_fric', long_name = 'phi_fric', units = 'degrees')
      CALL add_field_grid_dp_2D_notime( filename_grid, ncid_grid, 'phi_fric', long_name = 'phi_fric', units = 'degrees')
      ! Write data
      CALL write_to_field_multiple_options_mesh_dp_2D_notime( filename_mesh, ncid_mesh,                    'phi_fric', ice%phi_fric_a)
      CALL write_to_field_multiple_options_grid_dp_2D_notime( filename_grid, ncid_grid, mesh, region_name, 'phi_fric', ice%phi_fric_a)

    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! Close the files
    CALL close_netcdf_file( ncid_mesh)
    CALL close_netcdf_file( ncid_grid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_inverted_bed_roughness_to_file

  SUBROUTINE initialise_basal_inversion( mesh, ice)
    ! Fill in the initial guess for the bed roughness

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_inversion'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If needed, initialise target velocity fields
    ! ============================================

    IF     (C%choice_BIVgeo_method == 'PDC2012' .OR. &
            C%choice_BIVgeo_method == 'Lipscomb2021' .OR. &
            C%choice_BIVgeo_method == 'Bernales2017') THEN
      ! Not needed in these methods
    ELSEIF (C%choice_BIVgeo_method == 'CISM+' .OR. &
            C%choice_BIVgeo_method == 'Berends2022') THEN
      ! Needed in these methods

      CALL initialise_basal_inversion_target_velocity( mesh, ice)

    ELSE
      CALL crash('unknown choice_BIVgeo_method "' // TRIM(C%choice_BIVgeo_method) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_inversion

  SUBROUTINE initialise_basal_inversion_target_velocity( mesh, ice)
    ! Initialise the target velocity fields used in a velocity-based basal inversion routine

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_basal_inversion_target_velocity'
    CHARACTER(LEN=256)                                 :: filename
    TYPE(type_grid)                                    :: grid_raw
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  u_surf_raw,  v_surf_raw,  uabs_surf_raw
    INTEGER                                            :: wu_surf_raw, wv_surf_raw, wuabs_surf_raw
    INTEGER                                            :: i,j
    REAL(dp)                                           :: NaN

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename
    IF     (C%choice_BIVgeo_method == 'CISM+' .OR. &
            C%choice_BIVgeo_method == 'Berends2022') THEN
      filename = C%BIVgeo_target_velocity_filename
    ELSE
      CALL crash('unknown choice_BIVgeo_method "' // TRIM(C%choice_BIVgeo_method) // '"!')
    END IF

    IF (par%master) WRITE(0,*) '  Initialising basal inversion target velocity from file ', TRIM( filename), '...'

    ! Read x and y velocities from the NetCDF file
    CALL read_field_from_xy_file_2D( filename, 'u_surf', 'ANT', grid_raw, u_surf_raw, wu_surf_raw)
    CALL read_field_from_xy_file_2D( filename, 'v_surf', 'ANT', grid_raw, v_surf_raw, wv_surf_raw)

    ! NOTE: do not check for NaNs in the target velocity field. The Rignot 2011 Antarctica velocity product
    !       has a lot of missing data points, indicated by NaN values. This is acceptable, the basal inversion
    !       routine can handle that.

!    ! Safety
!    CALL check_for_NaN_dp_2D( u_surf_raw, 'u_surf_raw')
!    CALL check_for_NaN_dp_2D( v_surf_raw, 'v_surf_raw')

    ! Set missing values to NaN
    NaN = 0._dp
    NaN = 0._dp / NaN
    DO i = grid_raw%i1, grid_raw%i2
    DO j = 1, grid_raw%ny
      IF (u_surf_raw( i,j) == 0._dp .AND. v_surf_raw( i,j) == 0._dp) THEN
        u_surf_raw( i,j) = NaN
        v_surf_raw( i,j) = NaN
      END IF
    END DO
    END DO
    CALL sync

    ! Get absolute velocity
    DO i = grid_raw%i1, grid_raw%i2
    DO j = 1, grid_raw%ny
      uabs_surf_raw( i,j) = SQRT( u_surf_raw( i,j)**2 + v_surf_raw( i,j)**2)
    END DO
    END DO
    CALL sync

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, ice%BIV_uabs_surf_target, ice%wBIV_uabs_surf_target)

    ! Map raw data to the model mesh
    CALL map_from_xy_grid_to_mesh_2D( grid_raw, mesh, uabs_surf_raw, ice%BIV_uabs_surf_target)

    ! Deallocate raw data
    CALL deallocate_shared( grid_raw%wnx  )
    CALL deallocate_shared( grid_raw%wny  )
    CALL deallocate_shared( grid_raw%wdx  )
    CALL deallocate_shared( grid_raw%wx   )
    CALL deallocate_shared( grid_raw%wy   )
    CALL deallocate_shared( wu_surf_raw   )
    CALL deallocate_shared( wv_surf_raw   )
    CALL deallocate_shared( wuabs_surf_raw)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_basal_inversion_target_velocity


END MODULE basal_conditions_and_sliding_module
