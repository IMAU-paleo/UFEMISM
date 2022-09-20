MODULE ice_thickness_module

  ! Contains the routines for solving the ice thickness equation

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_SMB_model, type_BMB_model, &
                                             type_reference_geometry
  USE utilities_module,                ONLY: is_floating
  USE mesh_help_functions_module,      ONLY: rotate_xy_to_po_stag, find_containing_vertex
  USE ice_velocity_module,             ONLY: map_velocities_b_to_c_2D
  use mpi_module,                      only: allgather_array

  IMPLICIT NONE

CONTAINS

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
    INTEGER, DIMENSION(mesh%vi1:mesh%vi2), INTENT(in)  :: mask_noice
    TYPE(type_reference_geometry),       INTENT(IN)    :: refgeo_PD

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_dHi_dt'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Use the specified time integration method to calculate the ice thickness at t+dt
    IF     (C%choice_ice_integration_method == 'none') THEN
      ice%dHi_dt_a( mesh%vi1:mesh%vi2) = 0._dp
    ELSEIF (C%choice_ice_integration_method == 'explicit') THEN
       call calc_dHi_dt_explicit(     mesh, ice, SMB, BMB, dt)
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
       call crash("not implemented")
   !   CALL crash('calc_dHi_dt_semiimplicit: FIXME!')
      !CALL calc_dHi_dt_semiimplicit( mesh, ice, SMB, BMB, dt)
    ELSE
      CALL crash('unknown choice_ice_integration_method "' // TRIM( C%choice_ice_integration_method) // '"')
    END IF

    ! Apply boundary conditions
    CALL apply_ice_thickness_BC( mesh, ice, dt, mask_noice, refgeo_PD)

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
    REAL(dp), DIMENSION(:    ), allocatable            ::  u_c,  v_c,  up_c,  uo_c
    real(dp), dimension(:    ), allocatable            :: v_vav_b, u_vav_b, Hi_a
    real(dp), dimension(:,:  ), allocatable            :: Cw
    INTEGER                                            :: aci, vi, vj, cii, ci, cji, cj
    REAL(dp)                                           :: dVi, Vi_out, Vi_in, Vi_available, rescale_factor
    REAL(dp), DIMENSION(mesh%nV)                       :: Vi_SMB

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Initialise at zero
    ice%dVi_in     = 0._dp

    Vi_in          = 0._dp
    Vi_available   = 0._dp
    rescale_factor = 0._dp
    Vi_SMB         = 0._dp

    ! Calculate vertically averaged ice velocities along vertex connections
    allocate(u_c ( 1:mesh%nAc))
    allocate(v_c ( 1:mesh%nAc))
    allocate(up_c( 1:mesh%nAc))
    allocate(uo_c( 1:mesh%nAc))

    allocate(u_vav_b(mesh%nTri))
    allocate(v_vav_b(mesh%nTri))
    u_vav_b(mesh%ti1:mesh%ti2) = ice%u_vav_b
    v_vav_b(mesh%ti1:mesh%ti2) = ice%v_vav_b
    call allgather_array(u_vav_b)
    call allgather_array(v_vav_b)

    allocate(Hi_a( 1:mesh%nV))
    allocate(Cw  ( 1:mesh%nV, mesh%nC_mem))
    Hi_a(mesh%vi1:mesh%vi2) = ice%Hi_a
    Cw  (mesh%vi1:mesh%vi2,:) = mesh%Cw
    call allgather_array(Hi_a)
    call allgather_array(Cw)


    CALL map_velocities_b_to_c_2D( mesh, u_vav_b, v_vav_b, u_c, v_c)
    CALL rotate_xy_to_po_stag( mesh, u_c, v_c, up_c, uo_c)

    ! Calculate ice fluxes across all Aa vertex connections
    ! based on ice velocities calculated on Ac mesh
    ! =============================================

    DO aci = 1, mesh%nAc

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
        dVi = Hi_a( vi) * up_c( aci) * Cw( vi,ci) * dt ! m3
      ELSE
        dVi = Hi_a( vj) * up_c( aci) * Cw( vi,ci) * dt ! m3
      END IF

      ! Keep track of ice fluxes across individual connections, to correct for
      ! negative ice thicknesses if necessary
      !TODO because of this indexing, its not straightforward to parallelize (and therefore we don't parallelize, currently)
      ice%dVi_in( vi, ci) = -dVi ! m3
      ice%dVi_in( vj, cj) =  dVi ! m3

    END DO ! DO aci = mesh%ci1, mesh%ci2

    ! Correct outfluxes for possible resulting negative ice thicknesses
    ! =================================================================
    Vi_SMB( mesh%vi1:mesh%vi2) = (SMB%SMB_year( mesh%vi1:mesh%vi2) + BMB%BMB( mesh%vi1:mesh%vi2))  * mesh%A( mesh%vi1:mesh%vi2) * dt

    DO vi = mesh%vi1, mesh%vi2

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
      IF (-Vi_SMB( vi) >= Vi_available) THEN
        ! All ice in this vertex melts, nothing remains to move around. Rescale outfluxes to zero.
        Vi_SMB( vi) = -Vi_available
        rescale_factor = 0._dp
      END IF

      ! If the total outflux exceeds the available ice plus SMB plus total influx, rescale outfluxes
      IF (Vi_out > Vi_available + Vi_SMB( vi)) THEN
        ! Total outflux plus melt exceeds available ice volume. Rescale outfluxes to correct for this.
        rescale_factor = (Vi_available + Vi_SMB( vi)) / Vi_out
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

    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================

    DO vi = mesh%vi1, mesh%vi2
      dVi  = SUM( ice%dVi_in( vi,:))
      ice%dHi_dt_a(     vi) = (dVi + Vi_SMB( vi)) / (mesh%A( vi) * dt)
      ice%Hi_tplusdt_a( vi) = MAX( 0._dp, ice%Hi_a( vi) + ice%dHi_dt_a( vi) * dt)
      ice%dHi_dt_a(     vi) = (ice%Hi_tplusdt_a( vi) - ice%Hi_a( vi)) / dt
    END DO

    ! Clean up after yourself
    deallocate( u_c )
    deallocate( v_c )
    deallocate( up_c)
    deallocate( uo_c)
    deallocate( Hi_a)
    deallocate( Cw  )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_dHi_dt_explicit
  ! Some useful tools
  SUBROUTINE apply_ice_thickness_BC( mesh, ice, dt, mask_noice, refgeo_PD)
    ! Apply ice thickness boundary conditions (at the domain boundary, and through the mask_noice)

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: dt
    INTEGER, DIMENSION(mesh%vi1:mesh%vi2), INTENT(in)  :: mask_noice
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

    ! Remove ice in areas where no ice is allowed (e.g. Greenland in NAM and EAS, and Ellesmere Island in GRL)
    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice(     vi) == 1) THEN
        ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
        ice%Hi_tplusdt_a( vi) = 0._dp
      END IF
    END DO

    ! If so specified, remove all floating ice
    IF (C%do_remove_shelves) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (is_floating( ice%Hi_tplusdt_a( vi), ice%Hb_a( vi), ice%SL_a( vi))) THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        END IF
      END DO
    END IF ! IF (C%do_remove_shelves) THEN

    ! If so specified, remove all floating ice beyond the present-day calving front
    IF (C%remove_shelves_larger_than_PD) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (refgeo_PD%Hi( vi) == 0._dp .AND. refgeo_PD%Hb( vi) < 0._dp) THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        END IF
      END DO
    END IF ! IF (C%remove_shelves_larger_than_PD) THEN

    ! If so specified, remove all floating ice crossing the continental shelf edge
    IF (C%continental_shelf_calving) THEN
      CALL crash('continental_shelf_calving: FIXME!')
!      DO i = grid%i1, grid%i2
!      DO j = 1, grid%ny
!        IF (refgeo_GIAeq%Hi( j,i) == 0._dp .AND. refgeo_GIAeq%Hb( j,i) < C%continental_shelf_min_height) THEN
!          ice%dHi_dt_a(     j,i) = -ice%Hi_a( j,i) / dt
!          ice%Hi_tplusdt_a( j,i) = 0._dp
!        END IF
!      END DO
!      END DO
!      CALL sync
    END IF ! IF (C%continental_shelf_calving) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_ice_thickness_BC

END MODULE ice_thickness_module
