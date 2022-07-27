MODULE general_sea_level_module

! ===== Preamble =====
! ====================

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module,               ONLY: ice_density, ocean_area, seawater_density
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, allocate_shared_dp_0D, allocate_shared_dp_1D
  USE data_types_module,               ONLY: type_model_region, type_global_scalar_data
  USE forcing_module,                  ONLY: forcing, update_sealevel_record_at_model_time

  IMPLICIT NONE

CONTAINS

! ===== Regional sea level =====
! ==============================

  ! Update regional sea level
  SUBROUTINE update_regional_sea_level( NAM, EAS, GRL, ANT, GMSL_glob, time)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: NAM, EAS, GRL, ANT
    REAL(dp),                   INTENT(IN)        :: GMSL_glob
    REAL(dp),                   INTENT(IN)        :: time

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

      IF (C%do_NAM) NAM%ice%SL_a( NAM%mesh%vi1:NAM%mesh%vi2) = GMSL_glob
      IF (C%do_EAS) EAS%ice%SL_a( EAS%mesh%vi1:EAS%mesh%vi2) = GMSL_glob
      IF (C%do_GRL) GRL%ice%SL_a( GRL%mesh%vi1:GRL%mesh%vi2) = GMSL_glob
      IF (C%do_ANT) ANT%ice%SL_a( ANT%mesh%vi1:ANT%mesh%vi2) = GMSL_glob

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

  ! Calculate present-day sea level contribution
  SUBROUTINE calculate_PD_sealevel_contribution( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),    INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calculate_PD_sealevel_contribution'
    INTEGER                                       :: i, j
    REAL(dp)                                      :: ice_volume, thickness_above_flotation, ice_volume_above_flotation

    ! Add routine to path
    CALL init_routine( routine_name)

    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp

    DO i = region%refgeo_PD%grid%i1, region%refgeo_PD%grid%i2
    DO j = 1, region%refgeo_PD%grid%ny

      ! Thickness above flotation
      IF (region%refgeo_PD%Hi_grid( i,j) > 0._dp) THEN
        thickness_above_flotation = MAX(0._dp, region%refgeo_PD%Hi_grid( i,j) - MAX(0._dp, (0._dp - region%refgeo_PD%Hb_grid( i,j)) * (seawater_density / ice_density)))
      ELSE
        thickness_above_flotation = 0._dp
      END IF

      ! Ice volume (above flotation) in m.s.l.e
      ice_volume                 = ice_volume                 + region%refgeo_PD%Hi_grid( i,j)   * region%refgeo_PD%grid%dx * region%refgeo_PD%grid%dx * ice_density / (seawater_density * ocean_area)
      ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation        * region%refgeo_PD%grid%dx * region%refgeo_PD%grid%dx * ice_density / (seawater_density * ocean_area)

    END DO
    END DO

    CALL MPI_REDUCE( ice_volume                , region%ice_volume_PD,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_REDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation_PD, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calculate_PD_sealevel_contribution

  ! Determine current GMSL contributions from all simulated ice sheets
  SUBROUTINE determine_GMSL_contributions( GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, GMSL_glob, NAM, EAS, GRL, ANT, global_data)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),       INTENT(IN)     :: NAM, EAS, GRL, ANT
    REAL(dp),                      INTENT(INOUT)  :: GMSL_NAM, GMSL_EAS, GMSL_GRL, GMSL_ANT, GMSL_glob
    TYPE(type_global_scalar_data), INTENT(INOUT)  :: global_data

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'determine_GMSL_contributions'

    IF (.NOT. par%master) RETURN

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set GMSL contributions of all simulated ice sheets (text version)
    GMSL_NAM = 0._dp
    GMSL_EAS = 0._dp
    GMSL_GRL = 0._dp
    GMSL_ANT = 0._dp

    IF (C%do_NAM) GMSL_NAM = NAM%GMSL_contribution
    IF (C%do_EAS) GMSL_EAS = EAS%GMSL_contribution
    IF (C%do_GRL) GMSL_GRL = GRL%GMSL_contribution
    IF (C%do_ANT) GMSL_ANT = ANT%GMSL_contribution

    ! Determine global mean sea level (text version)
    GMSL_glob = GMSL_NAM + GMSL_EAS + GMSL_GRL + GMSL_ANT

    ! Set GMSL contributions of all simulated ice sheets (NetCDF version)
    global_data%GMSL_NAM = 0._dp
    global_data%GMSL_EAS = 0._dp
    global_data%GMSL_GRL = 0._dp
    global_data%GMSL_ANT = 0._dp

    IF (C%do_NAM) global_data%GMSL_NAM = NAM%GMSL_contribution
    IF (C%do_EAS) global_data%GMSL_EAS = EAS%GMSL_contribution
    IF (C%do_GRL) global_data%GMSL_GRL = GRL%GMSL_contribution
    IF (C%do_ANT) global_data%GMSL_ANT = ANT%GMSL_contribution

    ! Determine global mean sea level (NetCDF version)
    IF     (C%choice_sealevel_model == 'fixed') THEN
      global_data%GMSL = C%fixed_sealevel
    ELSEIF (C%choice_sealevel_model == 'eustatic' .OR. C%choice_sealevel_model == 'SELEN') THEN
      global_data%GMSL = global_data%GMSL_NAM + global_data%GMSL_EAS + global_data%GMSL_GRL + global_data%GMSL_ANT
    ELSEIF (C%choice_sealevel_model == 'prescribed') THEN
      global_data%GMSL = forcing%sealevel_obs
    ELSE
      CALL crash('unknown choice_sealevel_model "' // TRIM(C%choice_sealevel_model) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE determine_GMSL_contributions

END MODULE general_sea_level_module
