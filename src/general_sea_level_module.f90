module general_sea_level_module

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module, only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module,    only : ice_density, ocean_area, seawater_density
  use parallel_module,      only : par, sync, ierr
  use data_types_module,    only : type_model_region
  use forcing_module,       only : forcing
  implicit none

contains

! ===== Sea-level contributions =====
! ===================================

  ! Calculate regional ice volume and area
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

  ! Calculate present-day sea level contribution
  subroutine calculate_PD_sealevel_contribution( region)

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

end module general_sea_level_module
