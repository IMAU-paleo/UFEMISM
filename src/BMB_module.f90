module BMB_module
  ! Contains all the routines for calculating the basal mass balance.

  use mpi
  use configuration_module, only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module,    only : ice_density, seawater_density, cp_ocean, L_fusion, sec_per_year
  use parallel_module,      only : par, sync, ierr, cerr, partition_list
  use utilities_module,     only : check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                   check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                   interpolate_ocean_depth
  use netcdf_module,        only : debug
  use data_types_module,    only : type_mesh, type_ice_model, type_BMB_model, type_remapping_mesh_mesh, &
                                   type_climate_snapshot_regional, type_ocean_snapshot_regional, &
                                   type_remapping_mesh_mesh
  use forcing_module,       only : forcing, get_insolation_at_time_month_and_lat
  use reallocate_mod,       only : reallocate_bounds

  implicit none

contains

! ===== Main =====
! ================

  subroutine run_BMB_model( mesh, ice, ocean, BMB, region_name, time)
    ! Run the selected BMB model

    implicit none

    ! In/output variables
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ice_model),               intent(in)    :: ice
    type(type_ocean_snapshot_regional), intent(inout) :: ocean
    type(type_BMB_model),               intent(inout) :: BMB
    character(len=3),                   intent(in)    :: region_name
    real(dp),                           intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'run_BMB_model'
    integer                                           :: vi
    real(dp)                                          :: amplification_factor

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise at zero
    BMB%BMB(       mesh%vi1:mesh%vi2) = 0._dp
    BMB%BMB_sheet( mesh%vi1:mesh%vi2) = 0._dp
    BMB%BMB_shelf( mesh%vi1:mesh%vi2) = 0._dp

    ! === Ice sheet ===
    ! =================

    select case (C%choice_BMB_sheet_model)

      case ('uniform')
        ! Uniform BMB over whole domain
        BMB%BMB_sheet( mesh%vi1:mesh%vi2) = C%BMB_sheet_uniform

      case default
        ! Unknown case
        call crash('unknown choice_BMB_sheet_model "' // &
                    trim(C%choice_BMB_sheet_model) // '"')

    end select

    ! === Ice shelves ===
    ! ===================

    select case (C%choice_BMB_shelf_model)

      case ('uniform')
        ! Uniform BMB over whole domain
        BMB%BMB_shelf( mesh%vi1:mesh%vi2) = C%BMB_shelf_uniform

      case ('Favier2019_quad')
        ! Favier et al., 2019 quadratic parameterisation
        call run_BMB_model_Favier2019_quadratic( mesh, ice, ocean, BMB)

      case default
        ! Unknown case
        call crash('choice_BMB_shelf_model "' // &
                    TRIM(C%choice_BMB_shelf_model) // &
                    '" not implemented in run_BMB_model!')

    end select

    ! === Extrapolation ===
    ! =====================

    ! Extrapolate melt field from the regular (FCMP) mask to the (extended) PMP mask
    if (C%choice_BMB_subgrid == 'PMP') then
      call extrapolate_melt_from_FCMP_to_PMP( mesh, ice, BMB)
    end if

    ! === Final BMB ===
    ! =================

    ! Add sheet and shelf melt rates together, applying the selected scheme for sub-grid shelf melt
    ! (see Leguy et al. 2021 for explanations of the three schemes)
    do vi = mesh%vi1, mesh%vi2

      ! No sub-grid scaling for sub-sheet melt yet
      BMB%BMB( vi) = 0._dp
      if (ice%mask_sheet_a( vi) == 1._dp) then
        BMB%BMB( vi) = BMB%BMB_sheet( vi)
      end if

      ! Different sub-grid schemes for sub-shelf melt
      select case (C%choice_BMB_subgrid)

      case ('FCMP')
        ! Combination based on floating condition
        if (ice%mask_shelf_a( vi) == 1) then
          BMB%BMB( vi) = BMB%BMB( vi) + BMB%BMB_shelf( vi)
        end if

      case ('PMP')
        ! Combination based on floating area fractions
        BMB%BMB( vi) = BMB%BMB( vi) + (1._dp - ice%f_grnd_a( vi)) * BMB%BMB_shelf( vi)

      case ('NMP')
        ! Combination based on grounded area fractions
        if (ice%f_grnd_a( vi) == 0._dp) then
          BMB%BMB( vi) = BMB%BMB( vi) + BMB%BMB_shelf( vi)
        end if

      case default
        ! Unknown case
        call crash('unknown choice_BMB_subgrid "' // trim(C%choice_BMB_subgrid) // '"')

      end select

    end do

    ! === Finalisation ===
    ! ====================

    ! Safety
    call check_for_NaN_dp_1D( BMB%BMB_sheet, 'BMB%BMB_sheet')
    call check_for_NaN_dp_1D( BMB%BMB_shelf, 'BMB%BMB_shelf')
    call check_for_NaN_dp_1D( BMB%BMB,       'BMB%BMB'      )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_BMB_model

  subroutine initialise_BMB_model( mesh, ice, BMB, region_name)
    ! Allocate memory for the data fields of the SMB model.

    implicit none

    ! In/output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(in)    :: ice
    type(type_BMB_model), intent(inout) :: BMB
    character(len=3),     intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_BMB_model'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    if (par%master) then
      write (*,"(5A)") '  Initialising BMB model: sheet = "', &
                          trim(C%choice_BMB_sheet_model), &
                          '", shelf = "', TRIM(C%choice_BMB_shelf_model), '"...'
    end if

    ! === General ===
    ! ===============

    allocate(BMB%BMB      (mesh%vi1:mesh%vi2))
    allocate(BMB%BMB_shelf(mesh%vi1:mesh%vi2))
    allocate(BMB%BMB_sheet(mesh%vi1:mesh%vi2))

    ! === Ice sheet ===
    ! =================

    select case (C%choice_BMB_sheet_model)

      case ('uniform')
        ! Nothing else needs to be done

      case default
        ! Unknown case
        call crash('unknown choice_BMB_sheet_model "' // &
                    trim(C%choice_BMB_sheet_model) // '"')

    end select

    ! === Ice shelves ===
    ! ===================

    select case (C%choice_BMB_shelf_model)

      case ('uniform')
        ! Nothing else needs to be done

      case ('Favier2019_quad')
        ! Favier et al., 2019 parameterisations
        call initialise_BMB_model_Favier2019( mesh, BMB)

      case default
        ! Unknown case
        call crash('unknown choice_BMB_shelf_model "' // &
                    trim(C%choice_BMB_shelf_model) // '"')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_BMB_model

! ===== Favier et al. 2019 =====
! ==============================

  subroutine run_BMB_model_Favier2019_quadratic( mesh, ice, ocean, BMB)
    ! Calculate sub-shelf melt with Favier et al. (2019) quadratic parameterisation

    implicit none

    ! In/output variables
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ice_model),               intent(in)    :: ice
    type(type_ocean_snapshot_regional), intent(in)    :: ocean
    type(type_BMB_model),               intent(inout) :: BMB

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'run_BMB_model_Favier2019_quadratic'
    integer                                           :: vi
    real(dp)                                          :: dT

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Basal temperatures ===
    ! ==========================

    ! Calculate ocean temperature and freezing point at the base of the shelf
    call calc_ocean_temperature_at_shelf_base(    mesh, ice, ocean, BMB)
    call calc_ocean_freezing_point_at_shelf_base( mesh, ice, ocean, BMB)

    ! === Basal melt ===
    ! ==================

    do vi = mesh%vi1, mesh%vi2

      ! Initialise
      BMB%BMB_shelf( vi) = 0._dp

      if (ice%mask_shelf_a( vi) == 1) then

        ! Temperature forcing
        dT = BMB%T_ocean_base( vi) - BMB%T_ocean_freeze_base( vi)

        ! Favier et al. (2019), Eq. 4
        ! Altered to allow for negative basal melt (i.e. refreezing) when dT < 0
        BMB%BMB_shelf( vi) = -sec_per_year * C%BMB_Favier2019_quad_GammaT * sign(dT,1._dp) * (seawater_density * cp_ocean * dT / (ice_density * L_fusion))**2._dp

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_BMB_model_Favier2019_quadratic

  subroutine initialise_BMB_model_Favier2019( mesh, BMB)
    ! Allocate memory for the data fields of the Favier et al. (2019) shelf BMB parameterisations.

    implicit none

    ! In/output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_BMB_model), intent(inout) :: BMB

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_BMB_model_Favier2019'

    ! Add routine to path
    call init_routine( routine_name)

    ! Variables
    allocate( BMB%T_ocean_base       (mesh%vi1:mesh%vi2))
    allocate( BMB%T_ocean_freeze_base(mesh%vi1:mesh%vi2))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_BMB_model_Favier2019

! ===== Tools =====
! =================

  subroutine calc_ocean_temperature_at_shelf_base( mesh, ice, ocean, BMB)
    ! Calculate ocean temperature at the base of the shelf by interpolating
    ! the 3-D ocean temperature field in the vertical column

    implicit none

    ! In/output variables
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ice_model),               intent(in)    :: ice
    type(type_ocean_snapshot_regional), intent(in)    :: ocean
    type(type_BMB_model),               intent(inout) :: BMB

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'calc_ocean_temperature_at_shelf_base'
    integer                                           :: vi
    real(dp)                                          :: depth

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Initialise at zero
      BMB%T_ocean_base( vi) = 0._dp

      if (ice%mask_shelf_a( vi) == 1) then

        ! Calculate depth
        depth = max( 0.1_dp, ice%Hi_a( vi) * ice_density / seawater_density)

        ! Find ocean temperature at this depth
        call interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%T_ocean_corr_ext( vi,:), depth, BMB%T_ocean_base( vi))

      end if ! IF (ice%mask_shelf_a( vi) == 1) THEN

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ocean_temperature_at_shelf_base

  subroutine calc_ocean_freezing_point_at_shelf_base( mesh, ice, ocean, BMB)
    ! Calculate the ocean freezing point at the base of the shelf, needed to calculate
    ! basal melt in the different parameterisations from Favier et al. (2019)

    implicit none

    ! In/output variables
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ice_model),               intent(in)    :: ice
    type(type_ocean_snapshot_regional), intent(in)    :: ocean
    type(type_BMB_model),               intent(inout) :: BMB

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'calc_ocean_freezing_point_at_shelf_base'
    integer                                           :: vi
    real(dp)                                          :: depth
    real(dp)                                          :: S0                   ! Practical salinity [PSU]
    real(dp), parameter                               :: lambda1 = -0.0575_dp ! Liquidus slope                [degC PSU^-1] (Favier et al. (2019), Table 2)
    real(dp), parameter                               :: lambda2 = 0.0832_dp  ! Liquidus intercept            [degC]        (Favier et al. (2019), Table 2)
    real(dp), parameter                               :: lambda3 = 7.59E-4_dp ! Liquidus pressure coefficient [degC m^-1]   (Favier et al. (2019), Table 2)

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Initialise at zero
      BMB%T_ocean_freeze_base( vi) = 0._dp

      if (ice%mask_shelf_a( vi) == 1) then

        ! Calculate depth
        depth = max( 0.1_dp, ice%Hi_a( vi) * ice_density / seawater_density)

        ! Find salinity at this depth
        call interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S_ocean_corr_ext( vi,:), depth, S0)

        ! Calculate ocean freezing temperature (Favier et al. (2019), Eq. 3) in degrees Celsius
        BMB%T_ocean_freeze_base( vi) = lambda1 * S0 + lambda2 - lambda3 * depth

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ocean_freezing_point_at_shelf_base

  subroutine extrapolate_melt_from_FCMP_to_PMP( mesh, ice, BMB)
    ! All the BMB parameterisations are implicitly run using the FCMP sub-grid scheme
    ! (i.e. they are only applied to grid cells whose centre is floating).
    ! Calculating melt rates for partially-floating-but-grounded-at-the-centre cells (which
    ! is needed in the PMP scheme) is not straightforward; instead, just extrapolate values into
    ! the grid cells where this is the case.

    implicit none

    ! In/output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(in)    :: ice
    type(type_BMB_model), intent(inout) :: BMB

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'extrapolate_melt_from_FCMP_to_PMP'
    integer                             :: vi,i1,i2,j1,j2,ii,jj,n,n_ext
    integer,  dimension(:), allocatable :: mask_FCMP, mask_PMP
    real(dp), dimension(:), allocatable :: BMB_shelf_extra
    real(dp), dimension(:), allocatable :: dBMBdx, dBMBdy
    integer                             :: wmask_FCMP, wmask_PMP, wBMB_shelf_extra, wdBMBdx, wdBMBdy
    real(dp)                            :: BMB_av

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === WIP ===
    ! ===========

    call crash('this extrapolation has not been implemented yet!')

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extrapolate_melt_from_FCMP_to_PMP

! ===== Remapping =====
! =====================

  subroutine remap_BMB_model( mesh_old, mesh_new, map, BMB)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_BMB_model),           intent(inout) :: BMB

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_BMB_model'
    integer                                       :: int_dummy

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV
    int_dummy = map%int_dummy

    ! === General ===
    ! ===============

    ! Reallocate rather than remap; after a mesh update we'll immediately run the BMB model anyway
    call reallocate_bounds( BMB%BMB      , mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( BMB%BMB_sheet, mesh_new%vi1, mesh_new%vi2 )
    call reallocate_bounds( BMB%BMB_shelf, mesh_new%vi1, mesh_new%vi2 )

    ! === Ice sheet ===
    ! =================

    select case (C%choice_BMB_sheet_model)

      case ('uniform')
        ! Nothing else needs to be done

      case default
        ! Unknown case
        call crash('unknown choice_BMB_sheet_model "' // &
                    trim(C%choice_BMB_sheet_model) // '"')

    end select

    ! === Ice shelves ===
    ! ===================

    select case (C%choice_BMB_shelf_model)

      case ('uniform')
        ! Nothing else needs to be done

      case ('Favier2019_quad')
        ! Favier et al., 2019 parameterisations
        call reallocate_bounds( BMB%T_ocean_base, mesh_new%vi1, mesh_new%vi2 )
        call reallocate_bounds( BMB%T_ocean_freeze_base, mesh_new%vi1, mesh_new%vi2)

      case default
        ! Unknown case
        call crash('unknown choice_BMB_shelf_model "' // &
                    trim(C%choice_BMB_shelf_model) // '"')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_BMB_model

end module BMB_module
