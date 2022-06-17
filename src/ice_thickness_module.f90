MODULE ice_thickness_module
  !
  ! Contains the routines for solving the ice thickness equation

#include <petsc/finclude/petscksp.h>

! ===== Preamble =====
! ====================

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
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
  SUBROUTINE calc_dHi_dt( mesh, ice, SMB, BMB, dt, mask_noice, refgeo_PD)
    ! Use the total ice velocities to update the ice thickness

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_dHi_dt'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Use the specified time integration method to calculate the ice thickness at t+dt
    IF     (C%choice_ice_integration_method == 'none') THEN
      ice%dHi_dt_a( mesh%vi1:mesh%vi2) = 0._dp
      CALL sync
    ELSEIF (C%choice_ice_integration_method == 'explicit') THEN
      CALL calc_dHi_dt_explicit(     mesh, ice, SMB, BMB, dt)
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
    INTEGER                                            :: aci, vi, vj, cii, ci, cji, cj
    REAL(dp)                                           :: dVi, Vi_out, Vi_in, Vi_available, rescale_factor
    REAL(dp), DIMENSION(mesh%nV)                       :: Vi_MB

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise at zero
    ice%dVi_in(  mesh%vi1:mesh%vi2, :) = 0._dp
    ice%dVi_out( mesh%vi1:mesh%vi2, :) = 0._dp
    CALL sync

    Vi_in          = 0._dp
    Vi_available   = 0._dp
    rescale_factor = 0._dp
    Vi_MB         = 0._dp

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
      ! => And for ice-free ocean, where no accumulation is allowed
      IF (ice%mask_cf_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 1) THEN
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi))  * mesh%A( vi) * dt * ice%float_margin_frac_a( vi)
      ELSEIF (ice%mask_ocean_a( vi) == 1 .AND. ice%mask_shelf_a( vi) == 0) THEN
        Vi_MB( vi) = 0._dp
      ELSE
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi))  * mesh%A( vi) * dt
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

      rescale_factor = 1._dp

      ! If all the ice already present melts away, there can be no outflux.
      IF (-Vi_MB( vi) >= Vi_available) THEN
        ! All ice in this vertex melts, nothing remains to move around. Rescale outfluxes to zero.
        Vi_MB( vi) = -Vi_available
        rescale_factor = 0._dp
      END IF

      ! If the total outflux exceeds the available ice plus SMB plus total influx, rescale outfluxes
      IF (Vi_out > Vi_available + Vi_MB( vi)) THEN
        ! Total outflux plus melt exceeds available ice volume. Rescale outfluxes to correct for this.
        rescale_factor = (Vi_available + Vi_MB( vi)) / Vi_out
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

          END IF ! IF (ice%dVi_in( vi,ci) < 0._dp) THEN
        END DO ! DO ci = 1, mesh%nC( vi)
      END IF ! IF (rescale_factor < 1._dp) THEN

    END DO ! DO vi = mesh%vi1, mesh%vi2
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
