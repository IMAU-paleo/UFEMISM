MODULE general_sea_level_module

! ===== Preamble =====
! ====================

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module,               ONLY: ice_density, ocean_area, seawater_density
  USE parallel_module,                 ONLY: par, sync, ierr
  USE data_types_module,               ONLY: type_model_region, type_global_scalar_data
  USE forcing_module,                  ONLY: forcing, update_sealevel_record_at_model_time
  USE utilities_module,                ONLY: thickness_above_floatation

  IMPLICIT NONE

CONTAINS

! ===== Regional sea level =====
! ==============================

  ! Update regional sea level
  SUBROUTINE update_regional_sea_level( NAM, EAS, GRL, ANT, global_data, time)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),       INTENT(INOUT)  :: NAM, EAS, GRL, ANT
    TYPE(type_global_scalar_data), INTENT(INOUT)  :: global_data
    REAL(dp),                      INTENT(IN)     :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'update_regional_sea_level'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_sealevel_model == 'fixed') THEN

      IF (C%do_NAM) NAM%ice%SL_a( NAM%mesh%vi1:NAM%mesh%vi2) = C%fixed_sealevel
      IF (C%do_EAS) EAS%ice%SL_a( EAS%mesh%vi1:EAS%mesh%vi2) = C%fixed_sealevel
      IF (C%do_GRL) GRL%ice%SL_a( GRL%mesh%vi1:GRL%mesh%vi2) = C%fixed_sealevel
      IF (C%do_ANT) ANT%ice%SL_a( ANT%mesh%vi1:ANT%mesh%vi2) = C%fixed_sealevel

    ELSEIF (C%choice_sealevel_model == 'eustatic') THEN

      IF (C%do_NAM) NAM%ice%SL_a( NAM%mesh%vi1:NAM%mesh%vi2) = global_data%GMSL
      IF (C%do_EAS) EAS%ice%SL_a( EAS%mesh%vi1:EAS%mesh%vi2) = global_data%GMSL
      IF (C%do_GRL) GRL%ice%SL_a( GRL%mesh%vi1:GRL%mesh%vi2) = global_data%GMSL
      IF (C%do_ANT) ANT%ice%SL_a( ANT%mesh%vi1:ANT%mesh%vi2) = global_data%GMSL

    ELSEIF (C%choice_sealevel_model == 'prescribed') THEN

      CALL update_sealevel_record_at_model_time( time)

      IF (C%do_NAM) NAM%ice%SL_a( NAM%mesh%vi1:NAM%mesh%vi2) = forcing%sealevel_obs
      IF (C%do_EAS) EAS%ice%SL_a( EAS%mesh%vi1:EAS%mesh%vi2) = forcing%sealevel_obs
      IF (C%do_GRL) GRL%ice%SL_a( GRL%mesh%vi1:GRL%mesh%vi2) = forcing%sealevel_obs
      IF (C%do_ANT) ANT%ice%SL_a( ANT%mesh%vi1:ANT%mesh%vi2) = forcing%sealevel_obs

    ELSEIF (C%choice_sealevel_model == 'SELEN') THEN

      ! Sea level fields are filled in the SELEN routines

    ELSE
      CALL crash('unknown choice_sealevel_model "' // TRIM(C%choice_sealevel_model) // '"!')
    END IF

    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_regional_sea_level

! ===== Sea-level contributions =====
! ===================================

  ! Calculate regional ice volume and area
  SUBROUTINE calculate_icesheet_volume_and_area( region, do_ddt)
    ! Calculate this region's ice sheet's volume and area

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region), INTENT(INOUT) :: region
    LOGICAL,                 INTENT(IN)    :: do_ddt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER          :: routine_name = 'calculate_icesheet_volume_and_area'
    INTEGER                                :: vi
    REAL(dp)                               :: ice_area, ice_volume, ice_volume_above_flotation
    REAL(dp)                               :: ice_area_prev, ice_volume_prev, ice_volume_above_flotation_prev


    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise
    ice_area                   = 0._dp
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp

    ! Save previous value
    ice_area_prev                   = region%ice_area
    ice_volume_prev                 = region%ice_volume
    ice_volume_above_flotation_prev = region%ice_volume_above_flotation
    CALL sync

    ! Calculate ice area and volume for process domain
    DO vi = region%mesh%vi1, region%mesh%vi2

      IF (region%ice%mask_ice_a( vi) == 1) THEN
        ice_volume                 = ice_volume                 + (region%ice%Hi_a( vi) * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area))
        ice_area                   = ice_area                   + region%mesh%A( vi) * 1.0E-06_dp ! [km^2]
        ice_volume_above_flotation = ice_volume_above_flotation + region%ice%TAF_a( vi) * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)
      END IF

    END DO
    CALL sync

    CALL MPI_REDUCE( ice_area,                   region%ice_area,                   1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume,                 region%ice_volume,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    IF (par%master) THEN

      ! Calculate GMSL contribution
      region%GMSL_contribution = -1._dp * (region%ice_volume_above_flotation - region%ice_volume_above_flotation_PD)

      IF (do_ddt) THEN
        ! Calculate variation in area
        region%dice_area_dt = (region%ice_area - ice_area_prev) / region%dt
        ! Calculate variation in volume
        region%dice_volume_dt = (region%ice_volume - ice_volume_prev) / region%dt
        ! Calculate variation in volume above floatation
        region%dice_volume_above_flotation_dt = (region%ice_volume_above_flotation - ice_volume_above_flotation_prev) / region%dt
      END IF

    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_icesheet_volume_and_area

  ! Calculate present-day sea level contribution
  SUBROUTINE calculate_PD_sealevel_contribution( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calculate_PD_sealevel_contribution'
    INTEGER                                       :: vi
    REAL(dp)                                      :: ice_volume, thickness_above_flotation, ice_volume_above_flotation

    ! Add routine to path
    CALL init_routine( routine_name)

    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp

    DO vi = region%mesh%vi1, region%mesh%vi2

      ! Thickness above flotation
      IF (region%refgeo_PD%Hi( vi) > 0._dp) THEN
        thickness_above_flotation = MAX(0._dp, region%refgeo_PD%Hi( vi) - MAX(0._dp, (0._dp - region%refgeo_PD%Hb( vi)) * (seawater_density / ice_density)))
      ELSE
        thickness_above_flotation = 0._dp
      END IF

      ! Ice volume (above flotation) in m.s.l.e
      ice_volume                 = ice_volume                 + region%refgeo_PD%Hi( vi)   * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)
      ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation  * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)

    END DO
    CALL sync

    CALL MPI_REDUCE( ice_volume                , region%ice_volume_PD,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation_PD, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_PD_sealevel_contribution

  ! Determine current GMSL contributions from all simulated ice sheets
  SUBROUTINE determine_GMSL_contributions( NAM, EAS, GRL, ANT, global_data, time)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),       INTENT(IN)     :: NAM, EAS, GRL, ANT
    TYPE(type_global_scalar_data), INTENT(INOUT)  :: global_data
    REAL(dp),                      INTENT(IN)     :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'determine_GMSL_contributions'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      ! Set GMSL contributions of all simulated ice sheets (NetCDF version)
      global_data%GMSL_NAM = 0._dp
      global_data%GMSL_EAS = 0._dp
      global_data%GMSL_GRL = 0._dp
      global_data%GMSL_ANT = 0._dp

      IF (C%do_NAM) global_data%GMSL_NAM = NAM%GMSL_contribution
      IF (C%do_EAS) global_data%GMSL_EAS = EAS%GMSL_contribution
      IF (C%do_GRL) global_data%GMSL_GRL = GRL%GMSL_contribution
      IF (C%do_ANT) global_data%GMSL_ANT = ANT%GMSL_contribution

    END IF
    CALL sync

    ! Determine global mean sea level (NetCDF version)
    IF     (C%choice_sealevel_model == 'fixed') THEN
      global_data%GMSL = C%fixed_sealevel

    ELSEIF (C%choice_sealevel_model == 'eustatic' .OR. C%choice_sealevel_model == 'SELEN') THEN
      global_data%GMSL = global_data%GMSL_NAM + global_data%GMSL_EAS + global_data%GMSL_GRL + global_data%GMSL_ANT

    ELSEIF (C%choice_sealevel_model == 'prescribed') THEN
      CALL update_sealevel_record_at_model_time( time)
      global_data%GMSL = forcing%sealevel_obs

    ELSE
      CALL crash('unknown choice_sealevel_model "' // TRIM(C%choice_sealevel_model) // '"!')
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_GMSL_contributions

END MODULE general_sea_level_module
