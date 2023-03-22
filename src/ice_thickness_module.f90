MODULE ice_thickness_module
  !
  ! Contains the routines for solving the ice thickness equation

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module,               ONLY: pi, ice_density, seawater_density
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
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                             is_floating
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_SMB_model, type_BMB_model, &
                                             type_reference_geometry
  USE mesh_help_functions_module,      ONLY: rotate_xy_to_po_stag, find_containing_vertex
  USE ice_velocity_module,             ONLY: map_velocities_b_to_c_2D

  IMPLICIT NONE

CONTAINS

! ===== Compute new ice thickness at t+dt =====
! =============================================

  ! The main routine that is called from "run_ice_model" in the ice_dynamics_module
  SUBROUTINE calc_dHi_dt( mesh, ice, SMB, BMB, dt, mask_noice, refgeo_PD, time, do_dt_lim)
    ! Use the total ice velocities to update the ice thickness

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),               INTENT(IN)    :: mesh
    TYPE(type_ice_model),          INTENT(INOUT) :: ice
    TYPE(type_SMB_model),          INTENT(IN)    :: SMB
    TYPE(type_BMB_model),          INTENT(IN)    :: BMB
    REAL(dp),                      INTENT(INOUT) :: dt
    INTEGER,  DIMENSION(:    ),    INTENT(IN)    :: mask_noice
    TYPE(type_reference_geometry), INTENT(IN)    :: refgeo_PD
    REAL(dp),                      INTENT(IN)    :: time
    LOGICAL,                       INTENT(IN)    :: do_dt_lim

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                :: routine_name = 'calc_dHi_dt'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Use the specified time integration method to calculate the ice thickness at t+dt
    IF     (C%choice_ice_integration_method == 'none') THEN
      ice%dHi_dt_a( mesh%vi1:mesh%vi2) = 0._dp
      CALL sync
    ELSEIF (C%choice_ice_integration_method == 'explicit') THEN
      CALL calc_dHi_dt_explicit( mesh, ice, SMB, BMB, dt)
    ELSEIF (C%choice_ice_integration_method == 'dynamic') THEN
      CALL calc_dHi_dt_dynamic(  mesh, ice, SMB, BMB, dt, do_dt_lim, time)
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
      CALL crash('calc_dHi_dt_semiimplicit: FIXME!')
      !CALL calc_dHi_dt_semiimplicit( mesh, ice, SMB, BMB, dt)
    ELSE
      CALL crash('unknown choice_ice_integration_method "' // TRIM( C%choice_ice_integration_method) // '"')
    END IF

    ! Apply boundary conditions
    CALL apply_ice_thickness_BC( mesh, ice, dt, refgeo_PD)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt

  ! Different solvers for the ice thickness equation (explicit & semi-implicit)
  SUBROUTINE calc_dHi_dt_explicit( mesh, ice, SMB, BMB, dt)
    ! The explicit solver for the ice thickness equation

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_dHi_dt_explicit'
    REAL(dp), DIMENSION(:    ), POINTER                ::  u_c,  v_c,  up_c,  uo_c
    INTEGER                                            :: wu_c, wv_c, wup_c, wuo_c
    INTEGER                                            :: aci, vi, vj, cii, ci, cji, cj, n_ocn
    REAL(dp)                                           :: dVi, Vi_out, Vi_in, Vi_available, rescale_factor, m_ocn
    REAL(dp), DIMENSION(mesh%nV)                       :: Vi_MB

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise at zero
    ice%dVi_in(  mesh%vi1:mesh%vi2, :) = 0._dp
    ice%dVi_out( mesh%vi1:mesh%vi2, :) = 0._dp
    CALL sync

    Vi_in          = 0._dp
    Vi_out         = 0._dp
    Vi_available   = 0._dp
    rescale_factor = 0._dp
    Vi_MB          = 0._dp

    ! Calculate vertically averaged ice velocities along vertex connections
    CALL allocate_shared_dp_1D( mesh%nAc, u_c , wu_c )
    CALL allocate_shared_dp_1D( mesh%nAc, v_c , wv_c )
    CALL allocate_shared_dp_1D( mesh%nAc, up_c, wup_c)
    CALL allocate_shared_dp_1D( mesh%nAc, uo_c, wuo_c)

    CALL map_velocities_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, u_c, v_c)
    CALL rotate_xy_to_po_stag( mesh, u_c, v_c, up_c, uo_c)

    ! Calculate ice fluxes across all Aa vertex connections
    ! based on ice velocities calculated on Ac mesh
    ! =============================================

    DO aci = mesh%ci1, mesh%ci2

      ! The two Aa vertices connected by the Ac vertex
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)

      ! Find their own respective connectivity, for storing the ice flux
      ci = 0
      cj = 0
      DO cii = 1, mesh%nC( vi)
        IF (mesh%C( vi,cii)==vj) THEN
          ci = cii
          EXIT
        END IF
      END DO
      DO cji = 1, mesh%nC( vj)
        IF (mesh%C( vj,cji)==vi) THEN
          cj = cji
          EXIT
        END IF
      END DO

      ! Calculate ice volume per year moving along connection from vi to vj as the product of:
      ! - width          (m   - determined by distance between adjacent triangle circumcenters)
      ! - ice thickness  (m   - at flux origin (upwind scheme, because this is really advection)
      ! - ice velocity   (m/y - calculated at midpoint, using surface slope along connection)
      ! - time step      (y)
      IF (up_c( aci) > 0._dp) THEN
        dVi = ice%Hi_a( vi) * up_c( aci) * mesh%Cw( vi,ci) * dt ! m3
      ELSE
        dVi = ice%Hi_a( vj) * up_c( aci) * mesh%Cw( vi,ci) * dt ! m3
      END IF

      ! Keep track of ice fluxes across individual connections, to correct for
      ! negative ice thicknesses if necessary
      ice%dVi_in( vi, ci) = -dVi ! m3
      ice%dVi_in( vj, cj) =  dVi ! m3

    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync

    ! Correct outfluxes for possible resulting negative ice thicknesses
    ! =================================================================

    DO vi = mesh%vi1, mesh%vi2

      ! Ice volume added to each grid cell through the (surface + basal) mass balance
      ! => With an exception for the calving front, where we only apply
      !    the mass balance to the floating fraction
      ! => Plus a term accounting for the lateral melt at the calving front
      ! => And for ice-free ocean, where no accumulation is allowed
      IF (ice%mask_cf_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 1) THEN

        ! Basal melt propotional to the floating area fraction
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi)) * mesh%A( vi) * dt * MAX( .1_dp, ice%float_margin_frac_a( vi))

        ! Count number of ice-free ocean neighbours
        n_ocn = 0
        DO ci = 1, mesh%nC( vi)
          IF (ice%mask_ocean_a( ci) == 1 .AND. ice%mask_shelf_a( ci) == 0) THEN
            n_ocn = n_ocn + 1
          END IF
        END DO

        ! IF (mesh%lat( vi) > -80._dp .AND. mesh%lon( vi) > 240._dp .AND. mesh%lon( vi) < 270._dp) THEN

          ! Additional lateral mass balance proportional to the length of the ice column underwater and
          ! an approximation to the fraction of the Voronoi cell perimeter in contact with ocean neighbours
          Vi_MB( vi) = Vi_MB( vi) + BMB%BMB_shelf( vi) * dt * &
                                    (ice%Hi_eff_cf_a( vi) * ice_density / seawater_density) * &
                                    (2._dp * pi * mesh%R( vi) * REAL(n_ocn,dp) / REAL(mesh%nC( vi),dp))

        ! END IF

      ELSEIF (ice%mask_cf_a( vi) == 1 .AND. ice%mask_sheet_a( vi) == 1) THEN

        ! Basal melt propotional to ... nothing special for now
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi)) * mesh%A( vi) * dt

        ! Count number of ice-free ocean neighbours and get the sum of their basal melt rates
        n_ocn = 0
        m_ocn = 0._dp
        DO ci = 1, mesh%nC( vi)
          IF (ice%mask_ocean_a( ci) == 1 .AND. ice%mask_shelf_a( ci) == 0) THEN
            n_ocn = n_ocn + 1
            m_ocn = m_ocn + BMB%BMB_shelf( ci)
          END IF
        END DO

        ! IF (mesh%lat( vi) > -80._dp .AND. mesh%lon( vi) > 240._dp .AND. mesh%lon( vi) < 270._dp) THEN

          ! Additional lateral mass balance proportional to the length of the ice column underwater and
          ! an approximation to the fraction of the Voronoi cell perimeter in contact with ocean neighbours
          Vi_MB( vi) = Vi_MB( vi) + m_ocn / REAL(n_ocn,dp) * dt * &
                                    ( MAX( 0._dp, ice%SL_a( vi) - ice%Hb_a(vi)) * ice_density / seawater_density) * &
                                    (2._dp * pi * mesh%R( vi) * REAL(n_ocn,dp) / REAL(mesh%nC( vi),dp))

        ! END IF

      ELSEIF (ice%mask_ocean_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 0) THEN
        ! No mass balance for ice-free ocean vertices
        Vi_MB( vi) = 0._dp
      ELSE
        ! Any other point
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi)) * mesh%A( vi) * dt
      END IF

      ! Check how much ice is available for melting or removing (in m^3)
      Vi_available = mesh%A( vi) * ice%Hi_a( vi)

      dVi = SUM( ice%dVi_in( vi,:))

      Vi_in  = 0._dp
      Vi_out = 0._dp
      DO ci = 1, mesh%nC( vi)
        IF (ice%dVi_in( vi,ci) > 0._dp) THEN
          Vi_in  = Vi_in  + ice%dVi_in( vi,ci)
        ELSE
          Vi_out = Vi_out - ice%dVi_in( vi,ci)
        END IF
      END DO

      ! Default value
      rescale_factor = 1._dp

      ! If all the ice already present melts away, there can be no outflux.
      IF (-Vi_MB( vi) >= Vi_available) THEN
        ! All ice in this vertex melts, nothing remains to move around. Rescale outfluxes to zero.
        rescale_factor = 0._dp
        ! Limit the (negative, we know from the IF condition) magnitude of the mass balance to be
        ! at most equal to the mass already there plus the incoming flux
        Vi_MB( vi) = MAX( Vi_MB( vi), -(Vi_available + Vi_in))
      END IF

      ! If the total outflux exceeds the available ice plus SMB plus total influx, rescale outfluxes
      IF (Vi_out > Vi_available + Vi_MB( vi)) THEN
        ! Total outflux plus melt exceeds available ice volume. Rescale outfluxes to correct for this.
        rescale_factor = MIN( 1._dp, MAX( 0._dp, (Vi_available + Vi_MB( vi)) / Vi_out))
      END IF

      ! If this vertex is a sub-grid calving front, don't allow for outfluxes until it is full
      IF (ice%mask_cf_a( vi) == 1 .AND. &
          ice%mask_shelf_a( vi) == 1 .AND. &
          ice%float_margin_frac_a( vi) < 1._dp) THEN
        ! Rescale outfluxes to zero
        rescale_factor = 0._dp
      END IF

      ! Rescale ice outfluxes out of vi and into vi's neighbours
      IF (rescale_factor < 1._dp) THEN
        DO ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)

          IF (ice%dVi_in( vi,ci) < 0._dp) THEN
            ice%dVi_in( vi,ci) = ice%dVi_in( vi,ci) * rescale_factor

            DO cji = 1, mesh%nC( vj)
              IF (mesh%C( vj,cji) == vi) THEN
                ice%dVi_in( vj,cji) = -ice%dVi_in( vi,ci)
                EXIT
              END IF
            END DO

          END IF ! (ice%dVi_in( vi,ci) < 0._dp)
        END DO ! ci = 1, mesh%nC( vi)
      END IF ! (rescale_factor < 1._dp)

    END DO ! vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================

    DO vi = mesh%vi1, mesh%vi2
      dVi  = SUM( ice%dVi_in( vi,:))
      ice%dHi_dt_a(     vi) = (dVi + Vi_MB( vi)) / (mesh%A( vi) * dt)
      ice%Hi_tplusdt_a( vi) = MAX( 0._dp, ice%Hi_a( vi) + ice%dHi_dt_a( vi) * dt)
      ice%dHi_dt_a(     vi) = (ice%Hi_tplusdt_a( vi) - ice%Hi_a( vi)) / dt
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wu_c )
    CALL deallocate_shared( wv_c )
    CALL deallocate_shared( wup_c)
    CALL deallocate_shared( wuo_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt_explicit

  SUBROUTINE calc_dHi_dt_dynamic( mesh, ice, SMB, BMB, dt, do_dt_lim, time)
    ! Another explicit solver for the ice thickness equation

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),      INTENT(IN)    :: mesh
    TYPE(type_ice_model), INTENT(INOUT) :: ice
    TYPE(type_SMB_model), INTENT(IN)    :: SMB
    TYPE(type_BMB_model), INTENT(IN)    :: BMB
    REAL(dp),             INTENT(INOUT) :: dt
    LOGICAL,              INTENT(IN)    :: do_dt_lim
    REAL(dp),             INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER       :: routine_name = 'calc_dHi_dt_dynamic'
    REAL(dp), DIMENSION(:), POINTER     ::  u_c,  v_c,  up_c,  uo_c,  Vi_MB
    INTEGER                             :: wu_c, wv_c, wup_c, wuo_c, wVi_MB
    INTEGER                             :: aci, vi, vj, cii, ci, cji, cj, n_ocn
    REAL(dp)                            :: dVi, Vi_available, m_ocn, dt_lim, dt_lim_min

    ! Initialisation
    ! ==============

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Calculate surface + basal + lateral mass balance (in m^3/yr)
    ! ============================================================

    CALL allocate_shared_dp_1D( mesh%nV, Vi_MB , wVi_MB)

    DO vi = mesh%vi1, mesh%vi2

      IF (ice%mask_cf_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 1) THEN
        ! Floating calving front

        ! Basal melt propotional to the floating area fraction
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi)) * mesh%A( vi) * MAX( .1_dp, ice%float_margin_frac_a( vi))

        ! Count number of ice-free ocean neighbours
        n_ocn = 0
        DO ci = 1, mesh%nC( vi)
          IF (ice%mask_ocean_a( ci) == 1 .AND. ice%mask_shelf_a( ci) == 0) THEN
            n_ocn = n_ocn + 1
          END IF
        END DO

        ! Additional lateral mass balance proportional to the length of the ice column underwater and
        ! an approximation to the fraction of the Voronoi cell perimeter in contact with ocean neighbours
        Vi_MB( vi) = Vi_MB( vi) + BMB%BMB_shelf( vi) * &
                                  (ice%Hi_eff_cf_a( vi) * ice_density / seawater_density) * &
                                  (2._dp * pi * mesh%R( vi) * REAL(n_ocn,dp) / REAL(mesh%nC( vi),dp))

      ELSEIF (ice%mask_cf_a( vi) == 1 .AND. ice%mask_sheet_a( vi) == 1) THEN
        ! Grounded calving front

        ! Basal melt propotional to ... nothing special for now
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi)) * mesh%A( vi)

        ! Count number of ice-free ocean neighbours and get the sum of their basal melt rates
        n_ocn = 0
        m_ocn = 0._dp
        DO ci = 1, mesh%nC( vi)
          IF (ice%mask_ocean_a( ci) == 1 .AND. ice%mask_shelf_a( ci) == 0) THEN
            n_ocn = n_ocn + 1
            m_ocn = m_ocn + BMB%BMB_shelf( ci)
          END IF
        END DO

        ! Additional lateral mass balance proportional to the length of the ice column underwater and
        ! an approximation to the fraction of the Voronoi cell perimeter in contact with ocean neighbours
        Vi_MB( vi) = Vi_MB( vi) + m_ocn / REAL(n_ocn,dp) * &
                                  ( MAX( 0._dp, ice%SL_a( vi) - ice%Hb_a(vi)) * ice_density / seawater_density) * &
                                  (2._dp * pi * mesh%R( vi) * REAL(n_ocn,dp) / REAL(mesh%nC( vi),dp))

      ELSEIF (ice%mask_ocean_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 0) THEN
        ! Ice-free ocean

        ! No accumulation allowed
        Vi_MB( vi) = 0._dp

      ELSE
        ! Any other point

        ! Full mass balance applied
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi)) * mesh%A( vi)

      END IF

    END DO
    CALL sync

    ! Calculate vertically averaged ice
    ! velocities along vertex connections
    ! ===================================

    CALL allocate_shared_dp_1D( mesh%nAc, u_c , wu_c )
    CALL allocate_shared_dp_1D( mesh%nAc, v_c , wv_c )
    CALL allocate_shared_dp_1D( mesh%nAc, up_c, wup_c)
    CALL allocate_shared_dp_1D( mesh%nAc, uo_c, wuo_c)

    CALL map_velocities_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, u_c, v_c)
    CALL rotate_xy_to_po_stag( mesh, u_c, v_c, up_c, uo_c)

    ! Calculate ice fluxes across all Aa vertex connections
    ! based on ice velocities calculated on Ac mesh
    ! =============================================

    ! Initialise at zero
    ice%dVi_in(  mesh%vi1:mesh%vi2, :) = 0._dp
    ice%dVi_out( mesh%vi1:mesh%vi2, :) = 0._dp
    CALL sync

    DO aci = mesh%ci1, mesh%ci2

      ! The two Aa vertices connected by the Ac vertex
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)

      ! Find their own respective connectivity, for storing the ice flux
      ci = 0
      cj = 0
      DO cii = 1, mesh%nC( vi)
        IF (mesh%C( vi,cii)==vj) THEN
          ci = cii
          EXIT
        END IF
      END DO
      DO cji = 1, mesh%nC( vj)
        IF (mesh%C( vj,cji)==vi) THEN
          cj = cji
          EXIT
        END IF
      END DO

      ! Calculate ice volume per year moving along connection from vi to vj as the product of:
      ! - width          (m   - determined by distance between adjacent triangle circumcenters)
      ! - ice thickness  (m   - at flux origin (upwind scheme, because this is really advection)
      ! - ice velocity   (m/y - calculated at midpoint, using surface slope along connection)
      IF (up_c( aci) > 0._dp) THEN
        dVi = ice%Hi_a( vi) * up_c( aci) * mesh%Cw( vi,ci) ! m3/yr
      ELSE
        dVi = ice%Hi_a( vj) * up_c( aci) * mesh%Cw( vi,ci) ! m3/yr
      END IF

      ! Keep track of ice fluxes across individual connections, to correct for
      ! negative ice thicknesses if necessary
      ice%dVi_in( vi, ci) = -dVi ! m3/yr
      ice%dVi_in( vj, cj) =  dVi ! m3/yr

    END DO ! DO aci = mesh%ci1, mesh%ci2
    CALL sync

    ! Compute correction for sub-grid floating calving fronts
    ! =======================================================

    DO vi = mesh%vi1, mesh%vi2

      ! If this vertex is a sub-grid floating calving front
      IF (ice%mask_cf_a( vi) == 1 .AND. &
          ice%mask_shelf_a( vi) == 1 .AND. &
          ice%float_margin_frac_a( vi) < 1._dp) THEN

        ! Don't allow for outfluxes until it is full
        DO ci = 1, mesh%nC( vi)

          ! Connection to one of the neighbours
          vj = mesh%C( vi,ci)

          ! If there is outflux from current vertex to that neighbour
          IF (ice%dVi_in( vi,ci) < 0._dp) THEN
            ! Set that outflux to zero
            ice%dVi_in( vi,ci) = 0._dp

            ! Find the connection from the neighbour's POV
            DO cji = 1, mesh%nC( vj)
              IF (mesh%C( vj,cji) == vi) THEN
                ! Set the corresponding influx to zero
                ice%dVi_in( vj,cji) = 0._dp
                ! And stop the search
                EXIT
              END IF
            END DO

          END IF ! ice%dVi_in( vi,ci) < 0._dp
        END DO ! ci = 1, mesh%nC( vi)

      END IF ! sub-grid floating calving front

    END DO ! vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Compute maximum time step that
    ! avoids a negative ice thickness
    ! ===============================

    ! Use the current vertex ice volume and the rate of ice outflux
    ! to determine the maximum dt that prevents moving more ice than
    ! available between neighbouring vertices (thus violating mass
    ! conservation). Since the mass balance term Vi_MB is considered a
    ! source or sink, it is not used in the actual computation of dt.

    IF (do_dt_lim) THEN

      ! Initialise time step limit
      dt_lim_min = 1e20

      DO vi = mesh%vi1, mesh%vi2

        ! Ice available for removing (in m^3)
        Vi_available = mesh%A( vi) * ice%Hi_a( vi)
        ! Influx of mass through vertex (in m^3/yr)
        dVi  = SUM( ice%dVi_in( vi,:))

        ! If there is ice, and there is a mass loss
        IF ( Vi_available + dVi*dt < 0._dp) THEN

          ! Compute maximum dt (in yr) based on
          ! available ice and the ice outflux
          dt_lim = -Vi_available / dVi
          ! Limit maximum dt to positive numbers greater than the min time step
          dt_lim = MAX ( C%dt_min, dt_lim)

          ! Save smallest time step found so far
          dt_lim_min = MIN( dt_lim, dt_lim_min)

        END IF

      END DO
      CALL sync

      ! Limit model time step accordingly, if necessary
      dt = MIN( dt, dt_lim_min)

    END IF ! do_dt_lim

    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================

    DO vi = mesh%vi1, mesh%vi2
      ! Influx of mass through vertex (in m^3/yr)
      dVi  = SUM( ice%dVi_in( vi,:))
      ! Change in ice thickness (in m/yr)
      IF (C%do_dHdt_target .AND. time < C%dHdt_target_turnoff) THEN
        ice%dHi_dt_a( vi) = (dVi + Vi_MB( vi)) / mesh%A( vi) - ice%dHi_dt_target_a( vi)
      ELSE
        ice%dHi_dt_a( vi) = (dVi + Vi_MB( vi)) / mesh%A( vi)
      END IF
      ! Predicted new ice thickness (in m)
      ice%Hi_tplusdt_a( vi) = MAX( 0._dp, ice%Hi_a( vi) + ice%dHi_dt_a( vi) * dt)
      ! Diagnostic ice thickness change (in m/yr)
      ice%dHi_dt_a(     vi) = (ice%Hi_tplusdt_a( vi) - ice%Hi_a( vi)) / dt
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wu_c  )
    CALL deallocate_shared( wv_c  )
    CALL deallocate_shared( wup_c )
    CALL deallocate_shared( wuo_c )
    CALL deallocate_shared( wVi_MB)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt_dynamic

! ===== Boundary conditions =====
! ===============================

  ! Apply ice thickness boundary conditions
  SUBROUTINE apply_ice_thickness_BC( mesh, ice, dt, refgeo_PD)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: dt
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_ice_thickness_BC'
    INTEGER                                            :: vi, vvi, vj, n

    ! Add routine to path
    CALL init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! West
      IF (mesh%edge_index( vi) == 6 .OR. mesh%edge_index( vi) == 7 .OR. mesh%edge_index( vi) == 8) THEN
        IF     (C%ice_thickness_west_BC == 'zero') THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        ELSEIF (C%ice_thickness_west_BC == 'infinite') THEN
          CALL crash('ice_thickness_west_BC = infinite - FIXME!')
        ELSEIF (C%ice_thickness_west_BC == 'periodic') THEN
          CALL crash('ice_thickness_west_BC = periodic - FIXME!')
        ELSEIF (C%ice_thickness_west_BC == 'ISMIP_HOM_F') THEN
          ice%dHi_dt_a(     vi) = (1000._dp - ice%Hi_a( vi)) / dt
          ice%Hi_tplusdt_a( vi) = 1000._dp
        ELSEIF (C%ice_thickness_west_BC == 'fixed') THEN
          ice%dHi_dt_a(     vi) = 0._dp
          ice%Hi_tplusdt_a( vi) = ice%Hi_a( vi)
        ELSE
          CALL crash('unknown ice_thickness_west_BC "' // TRIM( C%ice_thickness_west_BC) // '"!')
        END IF
      END IF

      ! East
      IF (mesh%edge_index( vi) == 2 .OR. mesh%edge_index( vi) == 3 .OR. mesh%edge_index( vi) == 4) THEN
        IF     (C%ice_thickness_east_BC == 'zero') THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        ELSEIF (C%ice_thickness_east_BC == 'infinite') THEN
          CALL crash('ice_thickness_east_BC = infinite - FIXME!')
        ELSEIF (C%ice_thickness_east_BC == 'periodic') THEN
          CALL crash('ice_thickness_east_BC = periodic - FIXME!')
        ELSEIF (C%ice_thickness_east_BC == 'ISMIP_HOM_F') THEN
          ice%dHi_dt_a(     vi) = (1000._dp - ice%Hi_a( vi)) / dt
          ice%Hi_tplusdt_a( vi) = 1000._dp
        ELSEIF (C%ice_thickness_east_BC == 'fixed') THEN
          ice%dHi_dt_a(     vi) = 0._dp
          ice%Hi_tplusdt_a( vi) = ice%Hi_a( vi)
        ELSE
          CALL crash('unknown ice_thickness_east_BC "' // TRIM( C%ice_thickness_east_BC) // '"!')
        END IF
      END IF

      ! South
      IF (mesh%edge_index( vi) == 4 .OR. mesh%edge_index( vi) == 5 .OR. mesh%edge_index( vi) == 6) THEN
        IF     (C%ice_thickness_south_BC == 'zero') THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        ELSEIF (C%ice_thickness_south_BC == 'infinite') THEN
          CALL crash('ice_thickness_south_BC = infinite - FIXME!')
        ELSEIF (C%ice_thickness_south_BC == 'periodic') THEN
          CALL crash('ice_thickness_south_BC = periodic - FIXME!')
        ELSEIF (C%ice_thickness_south_BC == 'ISMIP_HOM_F') THEN
          ice%dHi_dt_a(     vi) = (1000._dp - ice%Hi_a( vi)) / dt
          ice%Hi_tplusdt_a( vi) = 1000._dp
        ELSEIF (C%ice_thickness_south_BC == 'fixed') THEN
          ice%dHi_dt_a(     vi) = 0._dp
          ice%Hi_tplusdt_a( vi) = ice%Hi_a( vi)
        ELSE
          CALL crash('unknown ice_thickness_south_BC "' // TRIM( C%ice_thickness_south_BC) // '"!')
        END IF
      END IF

      ! North
      IF (mesh%edge_index( vi) == 8 .OR. mesh%edge_index( vi) == 1 .OR. mesh%edge_index( vi) == 2) THEN
        IF     (C%ice_thickness_north_BC == 'zero') THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        ELSEIF (C%ice_thickness_north_BC == 'infinite') THEN
          CALL crash('ice_thickness_north_BC = infinite - FIXME!')
        ELSEIF (C%ice_thickness_north_BC == 'periodic') THEN
          CALL crash('ice_thickness_north_BC = periodic - FIXME!')
        ELSEIF (C%ice_thickness_north_BC == 'ISMIP_HOM_F') THEN
          ice%dHi_dt_a(     vi) = (1000._dp - ice%Hi_a( vi)) / dt
          ice%Hi_tplusdt_a( vi) = 1000._dp
        ELSEIF (C%ice_thickness_north_BC == 'fixed') THEN
          ice%dHi_dt_a(     vi) = 0._dp
          ice%Hi_tplusdt_a( vi) = ice%Hi_a( vi)
        ELSE
          CALL crash('unknown ice_thickness_north_BC "' // TRIM( C%ice_thickness_north_BC) // '"!')
        END IF
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_ice_thickness_BC

END MODULE ice_thickness_module
