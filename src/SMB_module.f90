module SMB_module
  ! All the routines for calculating the surface mass balance from the climate forcing

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module, only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module,    only : ice_density, L_fusion, sec_per_year, T0, pi
  use parallel_module,      only : par, sync, ierr, cerr, partition_list
  use utilities_module,     only : check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                   check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  use netcdf_module,        only : debug, inquire_restart_file_SMB, &
                                   read_restart_file_SMB
  use data_types_module,    only : type_mesh, type_ice_model, &
                                   type_SMB_model, type_remapping_mesh_mesh, &
                                   type_climate_matrix_regional, &
                                   type_climate_snapshot_regional, &
                                   type_restart_data
  use forcing_module,       only : forcing
  use mesh_mapping_module,  only : remap_field_dp_2D, remap_field_dp_3D, &
                                   calc_remapping_operator_grid2mesh, map_grid2mesh_2D, map_grid2mesh_3D, &
                                   deallocate_remapping_operators_grid2mesh
  use reallocate_mod,       only : reallocate_bounds

  implicit none

  real(dp), parameter :: albedo_water        = 0.1_dp
  real(dp), parameter :: albedo_soil         = 0.2_dp
  real(dp), parameter :: albedo_ice          = 0.5_dp
  real(dp), parameter :: albedo_snow         = 0.85_dp
  real(dp), parameter :: initial_snow_depth  = 0.1_dp

contains

! ===== Main =====
! ================

  subroutine run_SMB_model( mesh, ice, climate_matrix, time, SMB, mask_noice)
    ! Run the selected regional SMB model

    implicit none

    ! In/output variables
    type(type_mesh),                       intent(in)    :: mesh
    type(type_ice_model),                  intent(in)    :: ice
    type(type_climate_matrix_regional),    intent(in)    :: climate_matrix
    real(dp),                              intent(in)    :: time
    type(type_SMB_model),                  intent(inout) :: SMB
    integer, dimension(mesh%vi1:mesh%vi2), intent(in)    :: mask_noice

    ! Local variables:
    character(len=256), parameter                        :: routine_name = 'run_SMB_model'
    integer                                              :: vi

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Choice of model ===
    ! =======================

    select case (C%choice_SMB_model)

      case ('uniform')
        ! Apply a simple uniform SMB

        do vi = mesh%vi1, mesh%vi2
          if (mask_noice( vi) == 0) then
            SMB%SMB_year( vi) = C%SMB_uniform
          else
            SMB%SMB_year( vi) = 0._dp
          end if
        end do

      case ('IMAU-ITM')
        ! Run the IMAU-ITM SMB model
        call run_SMB_model_IMAUITM( mesh, ice, climate_matrix%applied, SMB, mask_noice)

      case default
        !Unknown case
        call crash('unknown choice_SMB_model "' // &
                    trim(C%choice_SMB_model) // '"!')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model

  subroutine initialise_SMB_model( mesh, ice, SMB, region_name)
    ! Allocate memory for the data fields of the SMB model.

    implicit none

    ! In/output variables
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in)    :: ice
    type(type_SMB_model), intent(inout) :: SMB
    character(len=3),     intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_SMB_model'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write (*,"(3A)") '  Initialising regional SMB model "', &
                          trim(C%choice_SMB_model), '"...'
    end if
    call sync

    ! === Choice of model ===
    ! =======================

    select case (C%choice_SMB_model)

      case ('uniform')
        ! Only need yearly total SMB
        allocate( SMB%SMB_year(mesh%vi1:mesh%vi2))

      case ('IMAU-ITM')
        ! Initialise the IMAU-ITM SMB model
        call initialise_SMB_model_IMAU_ITM( mesh, ice, SMB, region_name)

      case default
        !Unknown case
        call crash('unknown choice_SMB_model "' // &
                    trim(C%choice_SMB_model) // '"!')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model

! ===== The IMAU-ITM SMB model =====
! ==================================

  subroutine run_SMB_model_IMAUITM( mesh, ice, climate, SMB, mask_noice)
    ! Run the IMAU-ITM SMB model.

    ! NOTE: all the SMB components are in meters of water equivalent;
    !       the end result (SMB and SMB_year) are in meters of ice equivalent.

    implicit none

    ! In/output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_climate_snapshot_regional),   intent(in)    :: climate
    type(type_SMB_model),                   intent(inout) :: SMB
    integer,  dimension(mesh%vi1:mesh%vi2), intent(in)    :: mask_noice

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'run_SMB_model_IMAUITM'
    integer                                               :: vi, m
    integer                                               :: mprev
    real(dp)                                              :: snowfrac, liquid_water, sup_imp_wat

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! === Background albedo ===
      ! =========================

      ! Default value
      SMB%AlbedoSurf( vi) = albedo_soil

      ! Open water value
      if (ice%mask_ocean_a( vi) == 1 .and. &
         (ice%mask_shelf_a( vi) == 0 .or. mask_noice( vi) == 1)) then
        SMB%AlbedoSurf( vi) = albedo_water
      end if

      ! Ice value
      if ( ice%mask_ice_a( vi) == 1) then
        SMB%AlbedoSurf( vi) = albedo_ice
      end if

      ! === Month loop ===
      ! ==================

      do m = 1, 12  ! Month loop

        ! Set previous month index
        ! ========================

        mprev = m - 1

        ! If January, then December
        if (mprev==0) then
          mprev = 12
        end if

        ! Monthly albedo
        ! ==============

        ! Compute value
        SMB%Albedo( vi,m) = min(albedo_snow, max( SMB%AlbedoSurf( vi), &
                                albedo_snow - (albedo_snow - SMB%AlbedoSurf( vi))  * &
                                exp(-15._dp * SMB%FirnDepth( vi,mprev)) - 0.015_dp * SMB%MeltPreviousYear( vi)))

        ! Reset for open water
        if (ice%mask_ocean_a( vi) == 1 .and. &
           (ice%mask_shelf_a( vi) == 0 .or. mask_noice( vi) == 1)) then
          SMB%Albedo( vi,m) = albedo_water
        end if

        ! Monthly ablation
        ! ================

        ! Determine ablation as function af surface temperature and
        ! albedo/insolation according to Bintanja et al. (2002)

        SMB%Melt( vi,m) = max(0._dp, ( SMB%C_abl_Ts        * max(0._dp, climate%T2m( vi,m) - T0) + &
                                       SMB%C_abl_Q         * (1.0_dp - SMB%Albedo( vi,m)) * climate%Q_TOA( vi,m) - &
                                       SMB%C_abl_constant) * sec_per_year / (L_fusion * 1000._dp * 12._dp))

        ! Monthly accumulation
        ! ====================

        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999),
        ! liquid water content (rain and melt water) and snowdepth
        snowfrac = MAX(0._dp, MIN(1._dp, 0.725_dp * (1 - ATAN((climate%T2m( vi,m) - T0) / 5.95_dp) / 1.8566_dp)))

        SMB%Snowfall( vi,m) = climate%Precip( vi,m) *          snowfrac
        SMB%Rainfall( vi,m) = climate%Precip( vi,m) * (1._dp - snowfrac)

        ! Add this month's snow accumulation to next month's initial snow depth.
        SMB%AddedFirn( vi,m) = SMB%Snowfall( vi,m) - SMB%Melt( vi,m)
        SMB%FirnDepth( vi,m) = MIN(10._dp, MAX(0._dp, SMB%FirnDepth( vi,mprev) + SMB%AddedFirn( vi,m) ))

      end do ! m = 1, 12

      ! === Annual Refreezing ===
      ! =========================

      ! According to Janssens & Huybrechts (1999, 2000)
      ! The refreezing (=effective retention) is the minimum value of the amount
      ! of super imposed water and the available liquid water. Calculated for the
      ! whole year, then divided equally over the 12 months. This is done to account
      ! for the fact that liquid water is mostly available in summer, but "refreezing
      ! potential" is mostly available in winter. Not capped at total precipitation
      ! to account for cases where there is more melt than precipitation (e.g. ANT).

      ! Total refreezing estimation
      sup_imp_wat  = SMB%C_refr * max(0._dp, T0 - sum(climate%T2m( vi,:))/12._dp)

      ! Total amount of available liquid water
      liquid_water = sum(SMB%Rainfall( vi,:)) + sum(SMB%Melt( vi,:))

      ! Effective refreezing
      SMB%Refreezing_year( vi) = min( sup_imp_wat, liquid_water)

      ! Limit it to ice areas only
      if (ice%mask_ice_a( vi)==0) then
        SMB%Refreezing_year( vi) = 0._dp
      end if

      ! === Monthly refreezing, runoff, and SMB ===
      ! ===========================================

      do m = 1, 12
        SMB%Refreezing( vi,m) = SMB%Refreezing_year( vi) / 12._dp
        SMB%Runoff(     vi,m) = SMB%Melt( vi,m) + SMB%Rainfall( vi,m) - SMB%Refreezing( vi,m)
        SMB%SMB(        vi,m) = SMB%Snowfall( vi,m) + SMB%Refreezing( vi,m) - SMB%Melt( vi,m)
      end do

      ! === Annual SMB ===
      ! ==================

      ! Sum of monthly fields
      SMB%SMB_year( vi) = sum(SMB%SMB( vi,:))

      ! === Annual ablation ===
      ! =======================

      ! Calculate total melt over this year, to be used for determining next year's albedo
      SMB%MeltPreviousYear( vi) = sum(SMB%Melt( vi,:))

    end do ! vi = mesh%vi1, mesh%vi2

    ! === Final quantities ===
    ! ========================

    ! Convert final SMB from water to ice equivalent
    SMB%SMB(      mesh%vi1:mesh%vi2,:) = SMB%SMB(      mesh%vi1:mesh%vi2,:) * 1000._dp / ice_density
    SMB%SMB_year( mesh%vi1:mesh%vi2  ) = SMB%SMB_year( mesh%vi1:mesh%vi2  ) * 1000._dp / ice_density

    ! === Finalisation ===
    ! ====================

    ! Safety
    call check_for_NaN_dp_1D( SMB%AlbedoSurf      , 'SMB%AlbedoSurf'      )
    call check_for_NaN_dp_2D( SMB%Albedo          , 'SMB%Albedo'          )
    call check_for_NaN_dp_2D( SMB%Melt            , 'SMB%Melt'            )
    call check_for_NaN_dp_2D( SMB%Snowfall        , 'SMB%Snowfall'        )
    call check_for_NaN_dp_2D( SMB%Rainfall        , 'SMB%Rainfall'        )
    call check_for_NaN_dp_2D( SMB%Refreezing      , 'SMB%Refreezing'      )
    call check_for_NaN_dp_2D( SMB%Runoff          , 'SMB%Runoff'          )
    call check_for_NaN_dp_2D( SMB%SMB             , 'SMB%SMB'             )
    call check_for_NaN_dp_2D( SMB%AddedFirn       , 'SMB%AddedFirn'       )
    call check_for_NaN_dp_2D( SMB%FirnDepth       , 'SMB%FirnDepth'       )
    call check_for_NaN_dp_1D( SMB%SMB_year        , 'SMB%SMB_year'        )
    call check_for_NaN_dp_1D( SMB%MeltPreviousYear, 'SMB%MeltPreviousYear')
    !call check_for_NaN_dp_1D( SMB%Albedo_year     , 'SMB%Albedo_year'     ) ! Unused, undefined

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_IMAUITM

  subroutine initialise_SMB_model_IMAU_ITM( mesh, ice, SMB, region_name)
    ! Allocate memory for the data fields of the SMB model.

    implicit none

    ! In/output variables
    type(type_mesh),      intent(inout) :: mesh
    type(type_ice_model), intent(in)    :: ice
    type(type_SMB_model), intent(inout) :: SMB
    character(len=3),     intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_SMB_model_IMAU_ITM'
    integer                             :: vi
    character(len=256)                  :: SMB_IMAUITM_choice_init_firn

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Data fields
    allocate(  SMB%AlbedoSurf      (mesh%vi1:mesh%vi2    ))
    allocate(  SMB%MeltPreviousYear(mesh%vi1:mesh%vi2    ))
    allocate(  SMB%FirnDepth       (mesh%vi1:mesh%vi2, 12))
    allocate(  SMB%Rainfall        (mesh%vi1:mesh%vi2, 12))
    allocate(  SMB%Snowfall        (mesh%vi1:mesh%vi2, 12))
    allocate(  SMB%AddedFirn       (mesh%vi1:mesh%vi2, 12))
    allocate(  SMB%Melt            (mesh%vi1:mesh%vi2, 12))
    allocate(  SMB%Refreezing      (mesh%vi1:mesh%vi2, 12))
    allocate(  SMB%Refreezing_year (mesh%vi1:mesh%vi2    ))
    allocate(  SMB%Runoff          (mesh%vi1:mesh%vi2, 12))
    allocate(  SMB%Albedo          (mesh%vi1:mesh%vi2, 12))
    !allocate(  SMB%Albedo_year     (mesh%vi1:mesh%vi2    ))
    allocate(  SMB%SMB             (mesh%vi1:mesh%vi2, 12))
    allocate(  SMB%SMB_year        (mesh%vi1:mesh%vi2    ))

    ! === Parameters ===
    ! ==================

    ! Initialise regional scalar parameters to specified values
    if     (region_name == 'NAM') then
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_NAM
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_NAM
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_NAM
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_NAM
    elseif (region_name == 'EAS') then
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_EAS
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_EAS
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_EAS
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_EAS
    elseif (region_name == 'GRL') then
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_GRL
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_GRL
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_GRL
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_GRL
    elseif (region_name == 'ANT') then
      SMB%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_ANT
      SMB%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_ANT
      SMB%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_ANT
      SMB%C_refr                   = C%SMB_IMAUITM_C_refr_ANT
    end if

    ! === Firn ===
    ! ============

    ! Initialisation choice
    if     (region_name == 'NAM') then
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_NAM
    elseif (region_name == 'EAS') then
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_EAS
    elseif (region_name == 'GRL') then
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_GRL
    elseif (region_name == 'ANT') then
      SMB_IMAUITM_choice_init_firn = C%SMB_IMAUITM_choice_init_firn_ANT
    end if

    select case (SMB_IMAUITM_choice_init_firn)

      case ('uniform')
        ! Initialise with a uniform firn layer over the ice sheet

        do vi = mesh%vi1, mesh%vi2
          if (ice%Hi_a( vi) > 0._dp) then
            SMB%FirnDepth(        vi,:) = C%SMB_IMAUITM_initial_firn_thickness
            SMB%MeltPreviousYear( vi  ) = 0._dp
          else
            SMB%FirnDepth(        vi,:) = 0._dp
            SMB%MeltPreviousYear( vi  ) = 0._dp
          end if
        end do

      case ('restart')
        ! Initialise with the firn layer of a previous run
        call crash('firn from restart not implemented yet!')

      case default
        ! Unknown case
        call crash('unknown SMB_IMAUITM_choice_init_firn "' // &
                    trim(SMB_IMAUITM_choice_init_firn) // '"!')

    end select

    ! === Albedo ===
    ! ==============

    ! Initialise albedo
    do vi = mesh%vi1, mesh%vi2

      ! Background albedo
      if (ice%Hb_a( vi) < 0._dp) then
        SMB%AlbedoSurf( vi) = albedo_water
      else
        SMB%AlbedoSurf( vi) = albedo_soil
      end if

      if (ice%Hi_a( vi) > 0._dp) then
        SMB%AlbedoSurf(  vi) = albedo_snow
      end if

      SMB%Albedo( vi,:) = SMB%AlbedoSurf( vi)

    end do

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_IMAU_ITM

! ===== Remapping =====
! =====================

  subroutine remap_SMB_model( mesh_old, mesh_new, map, SMB)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    type(type_mesh),                     intent(in)    :: mesh_old
    type(type_mesh),                     intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh),      intent(in)    :: map
    type(type_SMB_model),                intent(inout) :: SMB

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'remap_SMB_model'
    integer                                            :: int_dummy

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV
    int_dummy = int_dummy

    ! === Annual SMB ===
    ! ==================

    ! Reallocate annual SMB, needed by all methods
    call reallocate_bounds( SMB%SMB_year, mesh_new%vi1,mesh_new%vi2)

    ! === IMAU-ITM ===
    ! ================

    ! Reallocate IMAU-ITM stuff
    if (C%choice_SMB_model == 'IMAU-ITM') then
      ! Firn depth and melt-during-previous-year must be remapped
      call remap_field_dp_2D( mesh_old, mesh_new, map, SMB%MeltPreviousYear, 'trilin')
      call remap_field_dp_3D( mesh_old, mesh_new, map, SMB%FirnDepth,        'trilin')

      ! Reallocate rather than remap; after a mesh update we'll immediately run the BMB model anyway
      call reallocate_bounds( SMB%Q_TOA           ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%AlbedoSurf      ,mesh_new%vi1,mesh_new%vi2    )
      call reallocate_bounds( SMB%Rainfall        ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Snowfall        ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%AddedFirn       ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Melt            ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Refreezing      ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Refreezing_year ,mesh_new%vi1,mesh_new%vi2    )
      call reallocate_bounds( SMB%Runoff          ,mesh_new%vi1,mesh_new%vi2, 12)
      call reallocate_bounds( SMB%Albedo          ,mesh_new%vi1,mesh_new%vi2, 12)
      !call reallocate_bounds( SMB%Albedo_year     ,mesh_new%vi1,mesh_new%vi2    )
      call reallocate_bounds( SMB%SMB             ,mesh_new%vi1,mesh_new%vi2, 12)
    end if

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model

end module SMB_module
