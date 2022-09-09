module general_sea_level_module

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module, only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module,    only : ice_density, ocean_area, seawater_density
  use parallel_module,      only : par, sync, ierr
  use data_types_module,    only : type_model_region, type_global_scalar_data
  use forcing_module,       only : forcing, update_sealevel_at_model_time
  implicit none

contains

! ===== Regional sea level =====
! ==============================

  subroutine update_regional_sea_level( NAM, EAS, GRL, ANT, global_data, time)
    ! Update regional sea level

    implicit none

    ! In/output variables:
    type(type_model_region),       intent(inout) :: NAM, EAS, GRL, ANT
    type(type_global_scalar_data), intent(inout) :: global_data
    real(dp),                      intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'update_regional_sea_level'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_sealevel_model)

      case ('fixed')
        ! Fixed value

        if (C%do_NAM) NAM%ice%SL_a( NAM%mesh%vi1:NAM%mesh%vi2) = C%fixed_sealevel
        if (C%do_EAS) EAS%ice%SL_a( EAS%mesh%vi1:EAS%mesh%vi2) = C%fixed_sealevel
        if (C%do_GRL) GRL%ice%SL_a( GRL%mesh%vi1:GRL%mesh%vi2) = C%fixed_sealevel
        if (C%do_ANT) ANT%ice%SL_a( ANT%mesh%vi1:ANT%mesh%vi2) = C%fixed_sealevel

      case ('eustatic')
        ! Use eustatic sea level

        if (C%do_NAM) NAM%ice%SL_a( NAM%mesh%vi1:NAM%mesh%vi2) = global_data%GMSL
        if (C%do_EAS) EAS%ice%SL_a( EAS%mesh%vi1:EAS%mesh%vi2) = global_data%GMSL
        if (C%do_GRL) GRL%ice%SL_a( GRL%mesh%vi1:GRL%mesh%vi2) = global_data%GMSL
        if (C%do_ANT) ANT%ice%SL_a( ANT%mesh%vi1:ANT%mesh%vi2) = global_data%GMSL

      case ('prescribed')
        ! Sea level from record

        call update_sealevel_at_model_time( time)

        if (C%do_NAM) NAM%ice%SL_a( NAM%mesh%vi1:NAM%mesh%vi2) = forcing%sealevel_obs
        if (C%do_EAS) EAS%ice%SL_a( EAS%mesh%vi1:EAS%mesh%vi2) = forcing%sealevel_obs
        if (C%do_GRL) GRL%ice%SL_a( GRL%mesh%vi1:GRL%mesh%vi2) = forcing%sealevel_obs
        if (C%do_ANT) ANT%ice%SL_a( ANT%mesh%vi1:ANT%mesh%vi2) = forcing%sealevel_obs

      case ('SELEN')
        ! Sea level fields are filled in the SELEN routines

      case default
        ! Unknown case

        call crash('unknown choice_sealevel_model "' // &
                    trim(C%choice_sealevel_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_regional_sea_level

! ===== Sea-level contributions =====
! ===================================

  subroutine calculate_icesheet_volume_and_area( region)
    ! Calculate this region's ice sheet's volume and area

    implicit none

    ! In/output variables:
    type(type_model_region), intent(inout)  :: region

    ! Local variables:
    character(len=256), parameter           :: routine_name = 'calculate_icesheet_volume_and_area'
    integer                                 :: vi
    real(dp)                                :: ice_area, ice_volume, thickness_above_flotation, ice_volume_above_flotation

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ice_area                   = 0._dp
    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp

    ! === Area and volume ===
    ! =======================

    ! Calculate ice area and volume for this process' domain
    do vi = region%mesh%vi1, region%mesh%vi2

      if (region%ice%mask_ice_a( vi) == 1) then

        ! Thickness above flotation
        thickness_above_flotation = max(0._dp, region%ice%Hi_a( vi) - max(0._dp, (region%ice%SL_a( vi) - region%ice%Hb_a( vi)) * (seawater_density / ice_density)))

        ! Total ice area in km^3 for this process' domain
        ice_area   = ice_area   + region%mesh%A( vi) * 1.0E-06_dp ! [km^3]

        ! Total ice volume in m.s.l.e for this process' domain
        ice_volume = ice_volume + region%ice%Hi_a( vi) * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)

        ! Ice volume (above flotation) in m.s.l.e for this process' domain
        ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)

      end if

    end do

    ! Get total area, volume, and volume above floatation adding up all process domains
    call MPI_ALLREDUCE( ice_area,                   region%ice_area,                   1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( ice_volume,                 region%ice_volume,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! === GMSL contribution ===
    ! =========================

    ! Calculate GMSL contribution
    region%GMSL_contribution = -1._dp * (region%ice_volume_above_flotation - region%ice_volume_above_flotation_PD)

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calculate_icesheet_volume_and_area

  subroutine calculate_PD_sealevel_contribution( region)
    ! Calculate present-day sea level contribution

    implicit none

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'calculate_PD_sealevel_contribution'
    integer                                :: vi
    real(dp)                               :: ice_volume, thickness_above_flotation, ice_volume_above_flotation

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ice_volume                 = 0._dp
    ice_volume_above_flotation = 0._dp

    ! === Volume ===
    ! ==============

    do vi = region%mesh%vi1, region%mesh%vi2

      if (region%ice%mask_ice_a( vi) == 1) then

        ! Thickness above flotation
        thickness_above_flotation = max(0._dp, region%refgeo_PD%Hi( vi) - max(0._dp, (0._dp - region%refgeo_PD%Hb( vi)) * (seawater_density / ice_density)))

        ! Total ice volume in m.s.l.e for this process' domain
        ice_volume = ice_volume + region%refgeo_PD%Hi( vi) * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)

        ! Ice volume (above flotation) in m.s.l.e for this process' domain
        ice_volume_above_flotation = ice_volume_above_flotation + thickness_above_flotation * region%mesh%A( vi) * ice_density / (seawater_density * ocean_area)

      end if

    end do

    ! Get total volume and volume above floatation adding up all process domains
    call MPI_ALLREDUCE( ice_volume                , region%ice_volume_PD,                 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( ice_volume_above_flotation, region%ice_volume_above_flotation_PD, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calculate_PD_sealevel_contribution

  subroutine determine_GMSL_contributions( NAM, EAS, GRL, ANT, global_data, time)
    ! Determine current GMSL contributions from all simulated ice sheets

    implicit none

    ! In/output variables:
    type(type_model_region),       intent(in)    :: NAM, EAS, GRL, ANT
    type(type_global_scalar_data), intent(inout) :: global_data
    real(dp),                      intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'determine_GMSL_contributions'

    ! Add routine to path
    call init_routine( routine_name)

    ! Set GMSL contributions of all simulated ice sheets (NetCDF version)
    global_data%GMSL_NAM = 0._dp
    global_data%GMSL_EAS = 0._dp
    global_data%GMSL_GRL = 0._dp
    global_data%GMSL_ANT = 0._dp

    if (C%do_NAM) global_data%GMSL_NAM = NAM%GMSL_contribution
    if (C%do_EAS) global_data%GMSL_EAS = EAS%GMSL_contribution
    if (C%do_GRL) global_data%GMSL_GRL = GRL%GMSL_contribution
    if (C%do_ANT) global_data%GMSL_ANT = ANT%GMSL_contribution

    ! Determine global mean sea level (NetCDF version)
    select case (C%choice_sealevel_model)

      case ('fixed')
        ! Fixed value
        global_data%GMSL = C%fixed_sealevel

      case ('eustatic','SELEN')
        ! Eustatic sea level or SELEN
        global_data%GMSL = global_data%GMSL_NAM + global_data%GMSL_EAS + global_data%GMSL_GRL + global_data%GMSL_ANT

      case ('prescribed')
        ! Sea level from record
        CALL update_sealevel_at_model_time( time)
        global_data%GMSL = forcing%sealevel_obs

      case default
        ! Unknown case
        call crash('unknown choice_sealevel_model "' // &
                    trim(C%choice_sealevel_model) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_GMSL_contributions

end module general_sea_level_module
