module basal_conditions_and_sliding_module
  ! Contains all the routines for calculating the basal conditions underneath the ice.

! ===== Preamble =====
! ====================

  use mpi
  use configuration_module, only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module,    only : pi, grav, seawater_density, ice_density
  use parallel_module,      only : par, sync, ierr, cerr, partition_list
  use utilities_module,     only : check_for_NaN_dp_1D, SSA_Schoof2006_analytical_solution
  use data_types_module,    only : type_mesh, type_ice_model, type_remapping_mesh_mesh
  use reallocate_mod,       only : reallocate_bounds

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

    ! Initial values
    call calc_basal_conditions( mesh, ice)

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
        ! Already initialised; do nothing

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
        ! Values read from an external file; do nothing

      case ('inversion')
        ! Updated by the inversion routines; do nothing

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

  SUBROUTINE calc_sliding_law_Coulomb_regularised( mesh, ice, u_a, v_a, beta_a)
    ! Regularised Coulomb-type sliding law

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN) :: u_a
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(IN) :: v_a
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2), INTENT(OUT):: beta_a

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

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_Coulomb_regularised

  SUBROUTINE calc_sliding_law_ZoetIverson( mesh, ice, u_a, v_a, beta_a)
    ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
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

    ! Calculate beta
    DO vi = mesh%vi1, mesh%vi2

      ! Include a normalisation term following Bueler & Brown (2009) to prevent divide-by-zero errors.
      uabs = SQRT( C%slid_delta_v**2 + u_a( vi)**2 + v_a( vi)**2)

      ! Zoet & Iverson (2020), Eq. (3) (divided by u to give beta = tau_b / u)
      beta_a( vi) = ice%tauc_a( vi) * (uabs**(1._dp / C%slid_ZI_p - 1._dp)) * ((uabs + C%slid_ZI_ut)**(-1._dp / C%slid_ZI_p))

    END DO

    ! Safety
    CALL check_for_NaN_dp_1D( beta_a, 'beta_a')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_ZoetIverson

  SUBROUTINE calc_sliding_law_idealised(  mesh, ice, beta_a)
    ! Sliding laws for some idealised experiments

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(IN)    :: ice
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: beta_a

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_sliding_law_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    call crash('idealised sliding stuff not yet implemented!')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_sliding_law_idealised

! ===== Remapping =====
! =====================

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
    CALL remap_bed_roughness( mesh_new,  ice)

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
      CALL reallocate_bounds( ice%pore_water_pressure_a, mesh_new%vi1, mesh_new%vi2 )
      CALL reallocate_bounds( ice%overburden_pressure_a, mesh_new%vi1, mesh_new%vi2 )
      CALL reallocate_bounds( ice%Neff_a               , mesh_new%vi1, mesh_new%vi2 )
    ELSEIF (C%choice_basal_hydrology == 'Martin2011') THEN
      CALL reallocate_bounds( ice%pore_water_pressure_a, mesh_new%vi1, mesh_new%vi2 )
      CALL reallocate_bounds( ice%overburden_pressure_a, mesh_new%vi1, mesh_new%vi2 )
      CALL reallocate_bounds( ice%Neff_a               , mesh_new%vi1, mesh_new%vi2 )
    ELSE
      CALL crash('unknown choice_basal_hydrology "' // TRIM( C%choice_basal_hydrology) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_basal_hydrology

  SUBROUTINE remap_bed_roughness( mesh_new, ice)
    ! Remap or reallocate all the data fields

    IMPLICIT NONE

    ! In/output variables:
    !TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    !TYPE(type_remapping_mesh_mesh),      INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_bed_roughness'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! To prevent compiler warnings for unused variables

    ! Allocate shared memory
    IF     (C%choice_sliding_law == 'no_sliding') THEN
      ! No sliding allowed
    ELSEIF (C%choice_sliding_law == 'idealised') THEN
      ! Sliding laws for some idealised experiments
    ELSEIF (C%choice_sliding_law == 'Coulomb' .OR. &
            C%choice_sliding_law == 'Coulomb_regularised') THEN
      ! Regularised Coulomb-type sliding law
      CALL reallocate_bounds( ice%phi_fric_a , mesh_new%vi1, mesh_new%vi2 )
      CALL reallocate_bounds( ice%tauc_a , mesh_new%vi1, mesh_new%vi2 )
    ELSEIF (C%choice_sliding_law == 'Zoet-Iverson') THEN
      ! Zoet-Iverson sliding law (Zoet & Iverson, 2020)
      CALL reallocate_bounds( ice%phi_fric_a , mesh_new%vi1, mesh_new%vi2 )
      CALL reallocate_bounds( ice%tauc_a , mesh_new%vi1, mesh_new%vi2 )
    ELSE
      CALL crash('unknown choice_sliding_law "' // TRIM( C%choice_sliding_law) // '"!')
    END IF

    ! ! If bed roughness is prescribed, read it from the provided NetCDF file
    ! IF (C%choice_basal_roughness == 'prescribed') THEN
    !   CALL initialise_bed_roughness_from_file( mesh_new, ice)
    ! END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_bed_roughness

END MODULE basal_conditions_and_sliding_module
