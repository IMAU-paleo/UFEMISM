module basal_conditions_and_sliding_module
  ! Contains all the routines for calculating the basal conditions underneath the ice.

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module, only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module,    only : pi, grav, seawater_density, ice_density
  use parallel_module,      only : par, sync, ierr, cerr, partition_list
  use utilities_module,     only : check_for_NaN_dp_1D, extrapolate_Gaussian_floodfill_mesh
  USE mesh_mapping_module,  only : smooth_Gaussian_2D, remap_field_dp_2D
  use data_types_module,    only : type_mesh, type_ice_model, type_remapping_mesh_mesh, type_grid, type_reference_geometry
  use reallocate_mod,       only : reallocate_bounds
  use mpi_module,           only : allgather_array

  implicit none

contains

! ===== Main =====
! ================

  subroutine calc_sliding_law( mesh, ice, u_a, v_a, beta_a)
    ! Calculate the sliding term beta in the SSA/DIVA using the specified sliding law

    implicit none

    ! In- and output variables:
    type(type_mesh),        intent(in)    :: mesh
    type(type_ice_model),   intent(inout) :: ice
    real(dp), dimension(:), intent(in)    :: u_a
    real(dp), dimension(:), intent(in)    :: v_a
    real(dp), dimension(:), intent(out)   :: beta_a

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'calc_sliding_law'

    ! Add routine to path
    call init_routine( routine_name)

    select case(C%choice_sliding_law)

    case ('no_sliding')
      ! No sliding allowed (choice of beta is trivial)
      beta_a( mesh%vi1:mesh%vi2) = 0._dp

    case ('idealised')
      ! Sliding laws for some idealised experiments
      call calc_sliding_law_idealised( mesh, ice, beta_a)

    case ('Coulomb_regularised')
      ! Regularised Coulomb-type sliding law
      call calc_sliding_law_Coulomb_regularised( mesh, ice, u_a, v_a, beta_a)

    case ('Zoet-Iverson')
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      call calc_sliding_law_ZoetIverson( mesh, ice, u_a, v_a, beta_a)

    case default
      ! Unkown case
      call crash('unknown choice_sliding_law "' // &
                  trim( C%choice_sliding_law) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law

  subroutine calc_basal_conditions( mesh, ice)
    ! Determine the basal conditions underneath the ice

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_basal_conditions'

    ! Add routine to path
    call init_routine( routine_name)

    ! Basal hydrology
    call calc_basal_hydrology( mesh, ice)

    ! Bed roughness
    call calc_bed_roughness( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_conditions

  subroutine initialise_basal_conditions( mesh, ice)
    ! Allocation and initialisation

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_basal_conditions'

    ! Add routine to path
    call init_routine( routine_name)

    ! Basal hydrology
    call initialise_basal_hydrology( mesh, ice)

    ! Bed roughness
    call initialise_bed_roughness( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = huge( 1))

  end subroutine initialise_basal_conditions

! ===== Basal hydrology =====
! ===========================

  subroutine calc_basal_hydrology( mesh, ice)
    ! Calculate the pore water pressure and effective basal pressure

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_basal_hydrology'
    integer                             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! === Pore water pressure ===
    ! ===========================

    select case (C%choice_basal_hydrology)

      case ('saturated')
        ! Assume all marine till is saturated (i.e. pore water
        ! pressure is equal to water pressure at depth everywhere)
        call calc_pore_water_pressure_saturated( mesh, ice)

      case ('Martin2011')
        ! Martin et al. (2011) parameterisation
        call calc_pore_water_pressure_Martin2011( mesh, ice)

      case default
        ! Unknown case
        call crash('unknown choice_basal_hydrology "' // &
                    trim( C%choice_basal_hydrology) // '"!')

    end select

    ! === Overburden and effective pressure ===
    ! =========================================

    do vi = mesh%vi1, mesh%vi2
      ice%overburden_pressure_a( vi) = ice_density * grav * ice%Hi_a( vi)
      ice%Neff_a(                vi) = max( 0._dp, ice%overburden_pressure_a( vi) &
                                                 - ice%pore_water_pressure_a( vi))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_hydrology

  subroutine initialise_basal_hydrology( mesh, ice)
    ! Allocation and initialisation

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_basal_hydrology'

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    select case (C%choice_basal_hydrology)

      case ('saturated')
        ! Perma-saturated basal conditions
        allocate( ice%pore_water_pressure_a( mesh%vi1:mesh%vi2))
        allocate( ice%overburden_pressure_a( mesh%vi1:mesh%vi2))
        allocate( ice%Neff_a               ( mesh%vi1:mesh%vi2))

      case ('Martin2011')
        ! Parameterisation from Martin et al., 2011
        allocate( ice%pore_water_pressure_a( mesh%vi1:mesh%vi2))
        allocate( ice%overburden_pressure_a( mesh%vi1:mesh%vi2))
        allocate( ice%Neff_a               ( mesh%vi1:mesh%vi2))

      case default
        ! Unknown case
        call crash('unknown choice_basal_hydrology "' // &
                    trim( C%choice_basal_hydrology) // '"!')

    end select

    ! Initialise values
    call calc_basal_hydrology( mesh, ice)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_windows_expected = 3)

  end subroutine initialise_basal_hydrology

  subroutine calc_pore_water_pressure_saturated( mesh, ice)
    ! Calculate the pore water pressure
    !
    ! Assume all till is saturated, i.e. pore water pressure = -rho_w * g * Hb

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_pore_water_pressure_saturated'
    integer                             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      ice%pore_water_pressure_a( vi) = -seawater_density * grav * ice%Hb_a( vi)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_pore_water_pressure_saturated

  subroutine calc_pore_water_pressure_Martin2011( mesh, ice)
    ! Calculate the pore water pressure
    !
    ! Use the parameterisation from Martin et al. (2011)

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_pore_water_pressure_Martin2011'
    integer                             :: vi
    real(dp)                            :: lambda_p

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Pore water pressure scaling factor (Martin et al., 2011, Eq. 12)
      lambda_p = min( 1._dp, max( 0._dp, &
                                  1._dp - (ice%Hb_a( vi) - ice%SL_a( vi) - &
                                  C%Martin2011_hydro_Hb_min) / (C%Martin2011_hydro_Hb_max - &
                                                                C%Martin2011_hydro_Hb_min) ))

      ! Pore water pressure (Martin et al., 2011, Eq. 11)
      ice%pore_water_pressure_a( vi) = 0.96_dp * ice_density * grav * ice%Hi_a( vi) * lambda_p

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_pore_water_pressure_Martin2011

! ===== Bed roughness =====
! =========================

  subroutine calc_bed_roughness( mesh, ice)
    ! Calculate the bed roughness

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_bed_roughness'

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_sliding_law == 'no_sliding' .or. &
        C%choice_sliding_law == 'idealised') then
      call finalise_routine( routine_name)
      return
    end if

    select case (C%choice_basal_roughness)

      case ('uniform')
        ! No need to recalculate. Do nothing.

      case ('parameterised')
        ! Apply the chosen parameterisation of bed roughness

        if (C%choice_param_basal_roughness == 'Martin2011') THEN
          ! The Martin et al. (2011) parameterisation of basal roughness &
          ! (specifically the till friction angle and till yield stress)
          call calc_bed_roughness_Martin2011( mesh, ice)
        else
          call crash('unknown choice_param_basal_roughness "' // &
                      trim( C%choice_param_basal_roughness) // '"!')
        end if

      case ('prescribed')
        ! Values were read from an external file; do nothing

      case ('restart')
        ! Values were read from a restart file; do nothing

      case default
        ! Unknown case
        call crash('unknown choice_basal_roughness "' // &
                    trim( C%choice_basal_roughness) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bed_roughness

  subroutine initialise_bed_roughness( mesh, ice)
    ! Allocation and initialisation

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_bed_roughness'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Allocation ===
    ! ==================

    select case (C%choice_sliding_law)

      case ('no_sliding')
        ! No sliding allowed - do nothing

      case ('idealised')
        ! Idealised sliding law

      case ('Coulomb_regularised')
        ! Regularised Coulomb-type sliding law
        allocate( ice%phi_fric_a (mesh%vi1:mesh%vi2))
        allocate( ice%tauc_a(mesh%vi1:mesh%vi2))

        ! Initialise with uniform value
        ice%phi_fric_a( mesh%vi1:mesh%vi2) = C%slid_Coulomb_phi_fric_uniform

      case ('Zoet-Iverson')
        ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
        allocate( ice%phi_fric_a (mesh%vi1:mesh%vi2))
        allocate( ice%tauc_a(mesh%vi1:mesh%vi2))

        ! Initialise with uniform value
        ice%phi_fric_a( mesh%vi1:mesh%vi2) = C%slid_Coulomb_phi_fric_uniform

      case default
        ! Unknown case
        call crash('unknown choice_sliding_law "' // &
                    trim( C%choice_sliding_law) // '"!')

    end select

    ! === Initial values ===
    ! ======================

    select case (C%choice_basal_roughness)

      case ('uniform')
        ! Uniform values already assigned

      case ('parameterised')
        ! Apply the chosen parameterisation of bed roughness
        call calc_bed_roughness( mesh, ice)

      case ('prescribed')
        ! If bed roughness is prescribed, read it from the provided NetCDF file
        call initialise_bed_roughness_from_file( mesh, ice)

      case ('restart')
        ! Assign the values that have been already read from a restart file
        call crash('bed roughness from restart not yet implemented!')
        ! call initialise_bed_roughness_from_restart_data( mesh, ice)

      case default
        ! Unknown case
        call crash('unknown choice_basal_roughness "' // &
                    trim( C%choice_basal_roughness) // '"!')

    end select

    ! === Basal sliding inversion ===
    ! ===============================

    ! Basal inversion
    if (C%do_slid_inv) then
      call initialise_basal_sliding_inversion( mesh, ice)
    end if

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
      call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness

  subroutine calc_bed_roughness_Martin2011( mesh, ice)
    ! Calculate the till friction angle phi_fric and till yield stress
    ! tauc, using the till model by Martin et al. (2011), Eq. 10.

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_bed_roughness_Martin2011'
    integer                             :: vi
    real(dp)                            :: w_Hb

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. (C%choice_sliding_law == 'Coulomb_regularised' .or. &
               C%choice_sliding_law == 'Zoet-Iverson')) then
      call crash('only applicable when choice_sliding_law = ' // &
                 'Coulomb_regularised", or "Zoet-Iverson"!')
    end if

    ! === Bed roughness ===
    ! =====================

    do vi = mesh%vi1, mesh%vi2

      ! Martin et al. (2011)
      w_Hb = MIN( 1._dp, MAX( 0._dp, &
                              (ice%Hb_a( vi) - C%Martin2011till_phi_Hb_min) / &
                              (C%Martin2011till_phi_Hb_max - C%Martin2011till_phi_Hb_min) ))

      ice%phi_fric_a( vi) = (1._dp - w_Hb) * C%Martin2011till_phi_min + &
                                     w_Hb  * C%Martin2011till_phi_max

    end do

    ! === Finalisation ===
    ! ====================

    ! Safety
    call check_for_NaN_dp_1D( ice%phi_fric_a, 'ice%phi_fric_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_bed_roughness_Martin2011

  subroutine initialise_bed_roughness_from_file( mesh, ice)
    ! Initialise bed roughness with data from an external NetCDF file

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_bed_roughness_from_file'

    ! Add routine to path
    call init_routine( routine_name)

    ! WIP
    call crash('bed roughness from file not yet implemented!')

    select case (C%choice_sliding_law)

      case ('no_sliding')
        ! No sliding - do nothing

      case ('idealised')
        ! Idealised bed roughness - do nothing

      case ('Coulomb_regularised')
        ! Coulomb-type sliding law
        ! call initialise_bed_roughness_from_file_Coulomb( mesh, ice)

      case ('Zoet-Iverson')
        ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
        ! call initialise_bed_roughness_from_file_ZoetIverson( mesh, ice)

      case default
        call crash('unknown choice_sliding_law "' // &
                    trim( C%choice_sliding_law) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_from_file

! ===== Sliding laws =====
! ========================

  subroutine calc_sliding_law_Coulomb_regularised( mesh, ice, u_a, v_a, beta_a)
    ! Regularised Coulomb-type sliding law

    implicit none

    ! In- and output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)    :: u_a
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)    :: v_a
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(out)   :: beta_a

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_sliding_law_Coulomb_regularised'
    integer                                               :: vi
    real(dp)                                              :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the till yield stress from the till friction angle and the effective pressure
    do vi = mesh%vi1, mesh%vi2
      ice%tauc_a( vi) = tan((pi / 180._dp) * ice%phi_fric_a( vi)) * ice%Neff_a( vi)
    end do

    ! Calculate beta
    do vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = sqrt( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      beta_a( vi) = ice%tauc_a( vi) * uabs ** (C%slid_Coulomb_reg_q_plastic - 1._dp) / (C%slid_Coulomb_reg_u_threshold ** C%slid_Coulomb_reg_q_plastic)

    end do

    ! Safety
    call check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_Coulomb_regularised

  subroutine calc_sliding_law_ZoetIverson( mesh, ice, u_a, v_a, beta_a)
    ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)

    implicit none

    ! In- and output variables:
    type(type_mesh),        intent(in)    :: mesh
    type(type_ice_model),   intent(inout) :: ice
    real(dp), dimension(:), intent(in)    :: u_a
    real(dp), dimension(:), intent(in)    :: v_a
    real(dp), dimension(:), intent(out)   :: beta_a

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'calc_sliding_law_ZoetIverson'
    integer                               :: vi
    real(dp)                              :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the till yield stress from the till friction angle and the effective pressure
    do vi = mesh%vi1, mesh%vi2
      ice%tauc_a( vi) = tan((pi / 180._dp) * ice%phi_fric_a( vi)) * ice%Neff_a( vi)
    end do

    ! Calculate beta
    do vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = sqrt( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Zoet & Iverson (2020), Eq. (3) (divided by u to give beta = tau_b / u)
      beta_a( vi) = ice%tauc_a( vi) * (uabs**(1._dp / C%slid_ZI_p - 1._dp)) * ((uabs + C%slid_ZI_ut)**(-1._dp / C%slid_ZI_p))

    end do

    ! Safety
    call check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_ZoetIverson

  subroutine calc_sliding_law_idealised(  mesh, ice, beta_a)
    ! Sliding laws for some idealised experiments

    implicit none

    ! In- and output variables:
    type(type_mesh),        intent(in)  :: mesh
    type(type_ice_model),   intent(in)  :: ice
    real(dp), dimension(:), intent(out) :: beta_a

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'calc_sliding_law_idealised'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('idealised sliding stuff not yet implemented!')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_sliding_law_idealised

! ===== Remapping =====
! =====================

  subroutine remap_basal_conditions( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    implicit none

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_ice_model),           intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_basal_conditions'

    ! Add routine to path
    call init_routine( routine_name)

    ! Basal hydrology
    call remap_basal_hydrology( mesh_old, mesh_new, map, ice)

    ! Bed roughness
    call remap_bed_roughness( mesh_old, mesh_new, map, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_basal_conditions

  subroutine remap_basal_hydrology( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    implicit none

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_ice_model),           intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_basal_hydrology'
    integer                                       :: int_dummy

    ! Add routine to path
    call init_routine( routine_name)

    ! To prevent compiler warnings for unused variables
    int_dummy = mesh_old%nV
    int_dummy = mesh_new%nV
    int_dummy = map%int_dummy

    ! Allocate shared memory
    select case (C%choice_basal_hydrology)

      case ('saturated')
        call reallocate_bounds( ice%pore_water_pressure_a, mesh_new%vi1, mesh_new%vi2 )
        call reallocate_bounds( ice%overburden_pressure_a, mesh_new%vi1, mesh_new%vi2 )
        call reallocate_bounds( ice%Neff_a               , mesh_new%vi1, mesh_new%vi2 )

      case ('Martin2011')
        call reallocate_bounds( ice%pore_water_pressure_a, mesh_new%vi1, mesh_new%vi2 )
        call reallocate_bounds( ice%overburden_pressure_a, mesh_new%vi1, mesh_new%vi2 )
        call reallocate_bounds( ice%Neff_a               , mesh_new%vi1, mesh_new%vi2 )

      case default
        call crash('unknown choice_basal_hydrology "' // trim( C%choice_basal_hydrology) // '"!')

    end select

    ! Reinitialise values
    call calc_basal_hydrology( mesh_new, ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_basal_hydrology

  subroutine remap_bed_roughness( mesh_old, mesh_new, map, ice)
    ! Remap or reallocate all the data fields

    implicit none

    ! In/output variables:
    type(type_mesh),                intent(in)    :: mesh_old
    type(type_mesh),                intent(in)    :: mesh_new
    type(type_remapping_mesh_mesh), intent(in)    :: map
    type(type_ice_model),           intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter                 :: routine_name = 'remap_bed_roughness'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Reallocation ===
    ! ====================

    select case (C%choice_sliding_law)

      case ('no_sliding')
        ! No sliding allowed. Nothing to do.

      case ('idealised')
        ! Sliding laws for some idealised experiments. Nothing to do.

      case ('Coulomb_regularised')
        ! Regularised Coulomb-type sliding law

        if (C%do_slid_inv) then
          ! Remap for inversion
          call remap_field_dp_2D( mesh_old, mesh_new, map, ice%phi_fric_a, 'nearest_neighbour')
          call remap_field_dp_2D( mesh_old, mesh_new, map, ice%phi_fric_inv_a, 'nearest_neighbour')
        else
          ! Reallocate
          call reallocate_bounds( ice%phi_fric_a , mesh_new%vi1, mesh_new%vi2 )
          ! Reinitialise
          ice%phi_fric_a( mesh_new%vi1:mesh_new%vi2) = C%slid_Coulomb_phi_fric_uniform
        end if

        call reallocate_bounds( ice%tauc_a , mesh_new%vi1, mesh_new%vi2 )

      case ('Zoet-Iverson')
        ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
        if (C%do_slid_inv) then
          ! Remap for inversion
          call remap_field_dp_2D( mesh_old, mesh_new, map, ice%phi_fric_a, 'nearest_neighbour')
          call remap_field_dp_2D( mesh_old, mesh_new, map, ice%phi_fric_inv_a, 'nearest_neighbour')
        else
          ! Reallocate
          call reallocate_bounds( ice%phi_fric_a , mesh_new%vi1, mesh_new%vi2 )
          ! Reinitialise
          ice%phi_fric_a( mesh_new%vi1:mesh_new%vi2) = C%slid_ZI_phi_fric_uniform
        end if

        call reallocate_bounds( ice%tauc_a , mesh_new%vi1, mesh_new%vi2 )

      case default
        call crash('unknown choice_sliding_law "' // trim( C%choice_sliding_law) // '"!')

    end select

    ! === Reinitialisation ===
    ! ========================

    select case (C%choice_basal_roughness)

      case ('uniform')
        ! Uniform values already assigned

      case ('parameterised')
        ! Apply the chosen parameterisation of bed roughness
        call calc_bed_roughness( mesh_new, ice)

      case ('prescribed')
        ! If bed roughness is prescribed, read it from the provided NetCDF file
        call initialise_bed_roughness_from_file( mesh_new, ice)

      case ('restart')
        ! Assign the values that have been already read from a restart file
        call crash('bed roughness from restart not yet implemented!')
        ! call initialise_bed_roughness_from_restart_data( mesh, ice)

      case default
        ! Unknown case
        call crash('unknown choice_basal_roughness "' // &
                    trim( C%choice_basal_roughness) // '"!')

    end select

    ! === Target velocities ===
    ! =========================

    ! For basal sliding inversion
    if (C%do_slid_inv .and. C%choice_slid_inv_method == 'Berends2022') then
      ! deallocate( ice%BIV_uabs_surf_target)
      ! call initialise_basal_inversion_target_velocity( mesh_new, ice)
      call crash('Berends2022 inversion not implemented yet!')
    end if

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_bed_roughness

! ===== Inversion =====
! =====================

  subroutine basal_sliding_inversion( mesh, grid, ice, refgeo, time)
    ! Iteratively invert for basal friction conditions under the grounded ice sheet,
    ! and extrapolate the resulting field over the rest of the domain

    implicit none

    ! In/output variables
    type(type_mesh),               intent(inout) :: mesh
    type(type_grid),               intent(in)    :: grid
    type(type_ice_model),          intent(inout) :: ice
    type(type_reference_geometry), intent(in)    :: refgeo
    real(dp),                      intent(in)    :: time

    ! Local variables
    character(len=256), parameter                :: routine_name = 'basal_sliding_inversion'

    ! Add routine to path
    call init_routine( routine_name)

    ! Apply the selected inversion scheme
    select case (C%choice_slid_inv_method)

      case ('Bernales2017')
        ! Ice thickness-based inversion
        call basal_sliding_inversion_Bernales2017( mesh, grid, ice, refgeo, time)

      case ('Berends2022')
        ! Ice thickness + velocity-based inversion
        call crash('Berends2022 inversion not implemented yet!')

      case default
        ! Unknown case
        call crash('unknown choice_slid_inv_method "' // trim( C%choice_slid_inv_method) // '"!')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine basal_sliding_inversion

  subroutine basal_sliding_inversion_Bernales2017( mesh, grid, ice, refgeo, time)
    ! Iteratively invert for basal friction conditions under the grounded ice sheet,
    ! and extrapolate the resulting field over the rest of the domain

    implicit none

    ! In/output variables
    type(type_mesh),               intent(inout) :: mesh
    type(type_grid),               intent(in)    :: grid
    type(type_ice_model),          intent(inout) :: ice
    type(type_reference_geometry), intent(in)    :: refgeo
    real(dp),                      intent(in)    :: time

    ! Local variables
    character(len=256), parameter                :: routine_name = 'basal_sliding_inversion_Bernales2017'
    integer                                      :: vi
    real(dp)                                     :: h_delta, t_scale, a_scale, m_scale, w_smooth
    integer,  dimension(:), allocatable          :: mask, mask_filled
    real(dp), dimension(:), allocatable          :: phi_inv

    ! Initialisation
    ! ==============

    ! Add routine to path
    call init_routine( routine_name)

    ! Initial checks
    if (time < C%slid_inv_t_start .or. &
        time > C%slid_inv_t_end) then
      ! Nothing yet/else to do for now. Just return.
      call finalise_routine( routine_name)
      return
    end if

    ! Allocate arrays for extrapolation
    allocate( mask (mesh%nV))
    allocate( mask_filled (mesh%nV))
    allocate( phi_inv (mesh%nV))
    mask = 0
    mask_filled = 0
    phi_inv = 0._dp

    ! Adjustment magnitude
    ! ====================

    ! Default values
    a_scale = C%slid_inv_Bernales2017_scale_start
    m_scale = .0_dp
    t_scale = 0._dp

    ! Time scale
    if (C%slid_inv_t_start < C%slid_inv_t_end) then
      ! Compute how much time has passed since start of inversion
      t_scale = (time - C%slid_inv_t_start) / (C%slid_inv_t_end - C%slid_inv_t_start)
      ! Limit t_scale to [0 1]
      t_scale = max( 0._dp, min( t_scale, 1._dp))
    end if

    ! Magnitude decay
    if (C%do_slid_inv_Bernales2017_decay) then
      ! Reduce adjustment amount as time goes on
      a_scale = C%slid_inv_Bernales2017_scale_start * (1._dp - t_scale) + C%slid_inv_Bernales2017_scale_end * t_scale
    end if

    ! Do the inversion
    ! ================

    do vi = mesh%vi1, mesh%vi2

      ! Ice thickness difference w.r.t. reference thickness
      h_delta = ice%Hi_a( vi) - refgeo%Hi( vi)

      ! Invert only where the model has grounded ice
      if (ice%mask_sheet_a( vi) == 1) then

        ! Mark this vertex as grounded ice
        mask( vi) = 2

          if ( h_delta >= 0._dp .and. ice%dHi_dt_a( vi) >= .0_dp ) then

            if (t_scale < .9_dp) then
              ice%phi_fric_inv_a( vi) = ice%phi_fric_inv_a( vi) - a_scale * (1._dp - exp( -abs( h_delta*ice%dHi_dt_a( vi))))
            else
              ice%phi_fric_inv_a( vi) = ice%phi_fric_inv_a( vi) - a_scale * (1._dp - exp( -abs( ice%dHi_dt_a( vi))))
            end if

          elseif ( h_delta <= 0._dp .and. ice%dHi_dt_a( vi) <= .0_dp ) then

            if (t_scale < .9_dp) then
              ice%phi_fric_inv_a( vi) = ice%phi_fric_inv_a( vi) + a_scale/3._dp * (1._dp - exp( -abs( h_delta*ice%dHi_dt_a( vi))))
            else
              ice%phi_fric_inv_a( vi) = ice%phi_fric_inv_a( vi) + a_scale/3._dp * (1._dp - exp( -abs( ice%dHi_dt_a( vi))))
            end if

          end if

          ! Constrain adjusted value to roughness limits
          ice%phi_fric_inv_a( vi) = min( max( ice%phi_fric_inv_a( vi), C%slid_inv_phi_min), C%slid_inv_phi_max)

      else

        ! This vertex is not grounded ice sheet, so mark it for later extrapolation
        mask( vi) = 1

      end if ! ice%mask_sheet_a( vi) == 1

    end do

    ! Communicate results
    ! ===================

    ! Gather mask info
    call allgather_array(mask)
    call allgather_array(mask_filled)

    ! Gather inverted bed roughness
    phi_inv( mesh%vi1:mesh%vi2) = ice%phi_fric_inv_a
    call allgather_array(phi_inv)

    ! Extrapolate the resulting field
    ! ===============================

    if (C%do_slid_inv_Bernales2017_extrap) then
      ! Perform the extrapolation
      call extrapolate_Gaussian_floodfill_mesh( mesh, mask, phi_inv, 40000._dp, mask_filled)
      ! Copy results to main variable
      ice%phi_fric_inv_a = phi_inv( mesh%vi1:mesh%vi2)
    end if

    ! Smoothing
    ! =========

    if (C%do_slid_inv_Bernales2017_smooth) then

      ! Smooth the resulting field
      call smooth_Gaussian_2D( mesh, grid, phi_inv, C%slid_inv_Bernales2017_smooth_r)

      do vi = mesh%vi1, mesh%vi2
        ! Combine the smoothed and raw inverted parameter through a weighed average
        ice%phi_fric_a( vi) = (1._dp - C%slid_inv_Bernales2017_smooth_w) * ice%phi_fric_inv_a( vi) + C%slid_inv_Bernales2017_smooth_w * phi_inv( vi)
        ! Make sure the variable stays within the prescribed limits
        ice%phi_fric_a( vi) = min( max( ice%phi_fric_a( vi), C%slid_inv_phi_min), C%slid_inv_phi_max)
      end do

    else
      ! Don't smooth the resulting field; simply copy it into the main variable
      ice%phi_fric_a = ice%phi_fric_inv_a

    end if

    ! Finalisation
    ! ============

    ! Clean up after yourself
    deallocate( mask_filled)
    deallocate( mask)
    deallocate( phi_inv)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine basal_sliding_inversion_Bernales2017

  subroutine initialise_basal_sliding_inversion( mesh, ice)
    ! Fill in the initial guess for the bed roughness

    implicit none

    ! Input variables:
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'initialise_basal_sliding_inversion'
    integer                             :: k

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Allocation ===
    ! ==================

    allocate( ice%phi_fric_inv_a (mesh%vi1:mesh%vi2))

    ! === Initial value ===
    ! =====================

    ice%phi_fric_inv_a = ice%phi_fric_a

    ! === Target velocity fields ===
    ! ==============================

    select case (C%choice_slid_inv_method)

      case ('Bernales2017')
        ! Not needed in this method

      case ('Berends2022')
        ! Needed in this method
        call crash('Berends2022 inversion not implemented yet!')

      case default
        call crash('unknown choice_slid_inv_method "' // trim(C%choice_slid_inv_method) // '"!')

    end select

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_basal_sliding_inversion

end module basal_conditions_and_sliding_module
