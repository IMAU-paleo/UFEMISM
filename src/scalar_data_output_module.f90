module scalar_data_output_module
  ! Contains all the routines for writing to the scalar output NetCDF files.

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module, only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parallel_module,      only : par, sync, ierr
  use data_types_module,    only : type_global_scalar_data, type_model_region, type_forcing_data
  use netcdf_module,        only : create_global_scalar_output_file, write_to_global_scalar_output_file, &
                                   create_regional_scalar_output_file, write_to_regional_scalar_output_file
  implicit none

contains

! ===== Global =====
! ==================

  subroutine write_global_scalar_data( NAM, EAS, GRL, ANT, forcing, global_data, time)
    ! Collect some global scalar data values and write to the NetCDF output file

    implicit none

    ! Input variables:
    type(type_model_region),       intent(in)    :: NAM, EAS, GRL, ANT
    type(type_forcing_data),       intent(in)    :: forcing
    type(type_global_scalar_data), intent(inout) :: global_data
    real(dp),                      intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'write_global_scalar_data'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then

      ! Sea level: variables have already been set in the UFEMISM main loop

      ! CO2
      if     (C%choice_forcing_method == 'none') then
        ! Do nothing
      elseif (C%choice_forcing_method == 'CO2_direct') then
        global_data%CO2_obs = forcing%CO2_obs
      elseif (C%choice_forcing_method == 'd18O_inverse_dT_glob') then
        ! Do nothing
      elseif (C%choice_forcing_method == 'd18O_inverse_CO2') then
        global_data%CO2_obs = forcing%CO2_obs
        global_data%CO2_mod = forcing%CO2_mod
      else
        call crash('unknown choice_forcing_method "' // trim( C%choice_forcing_method) // '"!')
      end if

      ! d18O
      if     (C%do_calculate_benthic_d18O) then
        global_data%dT_glob  = forcing%dT_glob
        global_data%dT_dw    = forcing%dT_deepwater
        global_data%d18O_obs = forcing%d18O_obs
        global_data%d18O_mod = forcing%d18O_mod
        global_data%d18O_ice = forcing%d18O_from_ice_volume_mod
        global_data%d18O_Tdw = forcing%d18O_from_temperature_mod
        global_data%d18O_NAM = forcing%d18O_NAM
        global_data%d18O_EAS = forcing%d18O_EAS
        global_data%d18O_GRL = forcing%d18O_GRL
        global_data%d18O_ANT = forcing%d18O_ANT
      end if

      ! Computation times
      global_data%tcomp_total   = 0._dp
      global_data%tcomp_ice     = 0._dp
      global_data%tcomp_thermo  = 0._dp
      global_data%tcomp_climate = 0._dp
      global_data%tcomp_GIA     = 0._dp

      if (C%do_NAM) then
        global_data%tcomp_total   = global_data%tcomp_total   + NAM%tcomp_total
        global_data%tcomp_ice     = global_data%tcomp_ice     + NAM%tcomp_ice
        global_data%tcomp_thermo  = global_data%tcomp_thermo  + NAM%tcomp_thermo
        global_data%tcomp_climate = global_data%tcomp_climate + NAM%tcomp_climate
        global_data%tcomp_GIA     = global_data%tcomp_GIA     + NAM%tcomp_GIA
      end if
      if (C%do_EAS) then
        global_data%tcomp_total   = global_data%tcomp_total   + EAS%tcomp_total
        global_data%tcomp_ice     = global_data%tcomp_ice     + EAS%tcomp_ice
        global_data%tcomp_thermo  = global_data%tcomp_thermo  + EAS%tcomp_thermo
        global_data%tcomp_climate = global_data%tcomp_climate + EAS%tcomp_climate
        global_data%tcomp_GIA     = global_data%tcomp_GIA     + EAS%tcomp_GIA
      end if
      if (C%do_GRL) then
        global_data%tcomp_total   = global_data%tcomp_total   + GRL%tcomp_total
        global_data%tcomp_ice     = global_data%tcomp_ice     + GRL%tcomp_ice
        global_data%tcomp_thermo  = global_data%tcomp_thermo  + GRL%tcomp_thermo
        global_data%tcomp_climate = global_data%tcomp_climate + GRL%tcomp_climate
        global_data%tcomp_GIA     = global_data%tcomp_GIA     + GRL%tcomp_GIA
      end if
      if (C%do_ANT) then
        global_data%tcomp_total   = global_data%tcomp_total   + ANT%tcomp_total
        global_data%tcomp_ice     = global_data%tcomp_ice     + ANT%tcomp_ice
        global_data%tcomp_thermo  = global_data%tcomp_thermo  + ANT%tcomp_thermo
        global_data%tcomp_climate = global_data%tcomp_climate + ANT%tcomp_climate
        global_data%tcomp_GIA     = global_data%tcomp_GIA     + ANT%tcomp_GIA
      end if

      ! Write to output file
      call write_to_global_scalar_output_file( global_data, time)

    end if ! (par%master)
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_global_scalar_data

  subroutine initialise_global_scalar_data( global_data)

    implicit none

    ! Input variables:
    type(type_global_scalar_data), intent(inout) :: global_data

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'initialise_global_scalar_data'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      ! Create the netcdf file
      call create_global_scalar_output_file( global_data%netcdf)
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_global_scalar_data

! ===== Regional =====
! ====================

  subroutine write_regional_scalar_data( region, time)
    ! Write some regionally integrated scalar values to the NetCDF output file

    implicit none

    ! Input variables:
    type(type_model_region), intent(inout) :: region
    real(dp),                intent(in)    :: time

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'write_regional_scalar_data'
    integer                                :: vi,m
    real(dp)                               :: T2m_mean
    real(dp)                               :: total_snowfall
    real(dp)                               :: total_rainfall
    real(dp)                               :: total_melt
    real(dp)                               :: total_refreezing
    real(dp)                               :: total_runoff
    real(dp)                               :: total_SMB
    real(dp)                               :: total_BMB
    real(dp)                               :: total_MB

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Annual mean surface temperature ===
    ! =======================================

    if (C%choice_climate_model == 'none') then
      ! In this case, no surface temperature is calculated at all
    else

      T2m_mean = 0._dp

      do vi = region%mesh%vi1, region%mesh%vi2
        T2m_mean = T2m_mean + sum( region%climate_matrix%applied%T2m( vi,:)) * region%mesh%A( vi) &
                            / (12._dp * (region%mesh%xmax - region%mesh%xmin) * (region%mesh%ymax - region%mesh%ymin))
      end do

      call MPI_ALLREDUCE( MPI_IN_PLACE, T2m_mean , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      region%int_T2m = T2m_mean

    end if

    ! === Total mass balance ===
    ! ==========================

    total_SMB = 0._dp
    total_BMB = 0._dp

    do vi = region%mesh%vi1, region%mesh%vi2
      if (region%ice%mask_ice_a( vi) == 1) then
        total_SMB = total_SMB + (region%SMB%SMB_year( vi) * region%mesh%A( vi) / 1E9_dp)
        total_BMB = total_BMB + (region%BMB%BMB(      vi) * region%mesh%A( vi) / 1E9_dp)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, T2m_mean , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, total_SMB, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, total_BMB, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    total_MB = total_SMB + total_BMB

    region%int_SMB = total_SMB
    region%int_BMB = total_BMB
    region%int_MB  = total_MB

    ! === SMB ===
    ! ===========

    ! Individual SMB components
    if     (C%choice_SMB_model == 'uniform') then
      ! Do nothing
    elseif (C%choice_SMB_model == 'IMAU-ITM') then

      total_snowfall             = 0._dp
      total_rainfall             = 0._dp
      total_melt                 = 0._dp
      total_refreezing           = 0._dp
      total_runoff               = 0._dp

      DO vi = region%mesh%vi1, region%mesh%vi2

        if (region%ice%Hi_a( vi) > 0._dp) then

          DO m = 1, 12
            total_snowfall   = total_snowfall   + (region%SMB%Snowfall(   vi,m) * region%mesh%A( vi) / 1E9_dp)
            total_rainfall   = total_rainfall   + (region%SMB%Rainfall(   vi,m) * region%mesh%A( vi) / 1E9_dp)
            total_melt       = total_melt       + (region%SMB%Melt(       vi,m) * region%mesh%A( vi) / 1E9_dp)
            total_refreezing = total_refreezing + (region%SMB%Refreezing( vi,m) * region%mesh%A( vi) / 1E9_dp)
            total_runoff     = total_runoff     + (region%SMB%Runoff(     vi,m) * region%mesh%A( vi) / 1E9_dp)
          END DO

        end if

      END DO

      call MPI_ALLREDUCE( MPI_IN_PLACE, total_snowfall  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, total_rainfall  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, total_melt      , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, total_refreezing, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, total_runoff    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)


      region%int_snowfall   = total_snowfall
      region%int_rainfall   = total_rainfall
      region%int_melt       = total_melt
      region%int_refreezing = total_refreezing
      region%int_runoff     = total_runoff

    else
      call crash('unknown choice_SMB_model "' // trim( C%choice_SMB_model) // '"!')
    end if

    ! === Write data ===
    ! ==================

    if (par%master) then
      ! Write to NetCDF file
      call write_to_regional_scalar_output_file( region, time)
    end if
    call sync

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_regional_scalar_data

  subroutine initialise_regional_scalar_data( region)

    implicit none

    ! Input variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'initialise_regional_scalar_data'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      ! Create the netcdf file
      call create_regional_scalar_output_file( region)
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_regional_scalar_data

end module scalar_data_output_module
