module ice_thickness_module
  ! Contains the routines for solving the ice thickness equation

! === Preamble ===
! ================

  use mpi
  use configuration_module,            only : dp, C, routine_path, init_routine, finalise_routine, crash
  use data_types_module,               only : type_mesh, type_ice_model, type_SMB_model, type_BMB_model, &
                                              type_reference_geometry
  use mesh_help_functions_module,      only : rotate_xy_to_po_stag
  use ice_velocity_module,             only : map_velocities_b_to_c_2D
  use mpi_module,                      only : allgather_array

  implicit none

contains

! ===== Compute new ice thickness at t+dt =====
! =============================================

  subroutine calc_dHi_dt( mesh, ice, SMB, BMB, dt, mask_noice, refgeo_PD)
    ! Use the total ice velocities to update the ice thickness

    implicit none

    ! In- and output variables
    type(type_mesh),                       intent(in)    :: mesh
    type(type_ice_model),                  intent(inout) :: ice
    type(type_SMB_model),                  intent(in)    :: SMB
    type(type_BMB_model),                  intent(in)    :: BMB
    real(dp),                              intent(in)    :: dt
    integer, dimension(mesh%vi1:mesh%vi2), intent(in)    :: mask_noice
    type(type_reference_geometry),         intent(in)    :: refgeo_PD

    ! Local variables:
    character(len=256), parameter                        :: routine_name = 'calc_dHi_dt'

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! === Time integration ===
    ! ========================

    select case (C%choice_ice_integration_method)

    case ('none')
      ice%dHi_dt_a( mesh%vi1:mesh%vi2) = 0._dp

    case ('explicit')
       call calc_dHi_dt_explicit( mesh, ice, SMB, BMB, dt)

    case ('semi-implicit')
       call crash("semi-implicit method not implemented yet!")

    case default
      call crash('unknown choice_ice_integration_method "' // trim( C%choice_ice_integration_method) // '"')

    end select

    ! === Boundary conditions ===
    ! ===========================

    ! Apply boundary conditions
    call apply_ice_thickness_BC( mesh, ice, dt, mask_noice, refgeo_PD)

    ! === Finalisation ===
    ! ====================

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dHi_dt

  subroutine calc_dHi_dt_explicit( mesh, ice, SMB, BMB, dt)
    ! The explicit solver for the ice thickness equation

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)      :: mesh
    type(type_ice_model), intent(inout)   :: ice
    type(type_SMB_model), intent(in)      :: SMB
    type(type_BMB_model), intent(in)      :: BMB
    real(dp),             intent(in)      :: dt

    ! Local variables:
    character(len=256), parameter         :: routine_name = 'calc_dHi_dt_explicit'
    real(dp), dimension(:  ), allocatable ::  u_c,  v_c,  up_c,  uo_c
    real(dp), dimension(:  ), allocatable :: v_vav_b, u_vav_b, Hi_a
    real(dp), dimension(:,:), allocatable :: Cw
    integer                               :: aci, vi, vj, cii, ci, cji, cj
    real(dp)                              :: dVi, Vi_out, Vi_in, Vi_available, rescale_factor
    real(dp), dimension(mesh%nV)          :: Vi_MB

    ! === Initialisation ===
    ! ======================

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise at zero
    ice%dVi_in     = 0._dp

    Vi_in          = 0._dp
    Vi_available   = 0._dp
    rescale_factor = 0._dp
    Vi_MB          = 0._dp

    ! Calculate vertically averaged ice velocities along vertex connections
    allocate(u_c ( mesh%ci1:mesh%ci2))
    allocate(v_c ( mesh%ci1:mesh%ci2))
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
    CALL rotate_xy_to_po_stag( mesh, u_c, v_c, up_c(mesh%ci1:mesh%ci2), uo_c(mesh%ci1:mesh%ci2))

    call allgather_array(up_c)
    call allgather_array(uo_c)

    ! Calculate ice fluxes across all Aa vertex connections
    ! based on ice velocities calculated on Ac mesh
    ! =============================================

    do aci = 1, mesh%nAc

      ! The two Aa vertices connected by the Ac vertex
      vi = mesh%Aci( aci,1)
      vj = mesh%Aci( aci,2)

      ! Find their own respective connectivity, for storing the ice flux
      ci = 0
      cj = 0
      do cii = 1, mesh%nC( vi)
        if (mesh%C( vi,cii)==vj) then
          ci = cii
          exit
        end if
      end do
      do cji = 1, mesh%nC( vj)
        if (mesh%C( vj,cji)==vi) then
          cj = cji
          exit
        end if
      end do

      ! Calculate ice volume per year moving along connection from vi to vj as the product of:
      ! - width          (m   - determined by distance between adjacent triangle circumcenters)
      ! - ice thickness  (m   - at flux origin (upwind scheme, because this is really advection)
      ! - ice velocity   (m/y - calculated at midpoint, using surface slope along connection)
      ! - time step      (y)
      if (up_c( aci) > 0._dp) then
        dVi = Hi_a( vi) * up_c( aci) * Cw( vi,ci) * dt ! m3
      else
        dVi = Hi_a( vj) * up_c( aci) * Cw( vi,ci) * dt ! m3
      end if

      ! Keep track of ice fluxes across individual connections, to correct for
      ! negative ice thicknesses if necessary
      !TODO because of this indexing, its not straightforward to parallelize (and therefore we don't parallelize, currently)
      ice%dVi_in( vi, ci) = -dVi ! m3
      ice%dVi_in( vj, cj) =  dVi ! m3

    end do ! aci = 1, mesh%nAc

    ! Correct outfluxes for possible resulting negative ice thicknesses
    ! =================================================================

    do vi = mesh%vi1, mesh%vi2

      ! Ice volume added to each grid cell through the (surface + basal) mass balance
      ! => With an exception for the calving front, where we only apply
      !    the mass balance to the floating fraction
      ! => And for ice-free ocean, where no accumulation is allowed
      if (ice%mask_cf_a( vi) == 1 .and. ice%mask_shelf_a( vi) == 1) then
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi))  * mesh%A( vi) * dt * ice%float_margin_frac_a( vi)
      elseif (ice%mask_ocean_a( vi) == 1 .and. ice%mask_shelf_a( vi) == 0) then
        Vi_MB( vi) = 0._dp
      else
        Vi_MB( vi) = (SMB%SMB_year( vi) + BMB%BMB( vi))  * mesh%A( vi) * dt
      end if

      ! Check how much ice is available for melting or removing (in m^3)
      Vi_available = mesh%A( vi) * ice%Hi_a( vi)

      dVi = sum( ice%dVi_in( vi,:))

      Vi_in  = 0._dp
      Vi_out = 0._dp
      do ci = 1, mesh%nC( vi)
        if (ice%dVi_in( vi,ci) > 0._dp) then
          Vi_in  = Vi_in  + ice%dVi_in( vi,ci)
        else
          Vi_out = Vi_out - ice%dVi_in( vi,ci)
        end if
      end do

      rescale_factor = 1._dp

      ! If all the ice already present melts away, there can be no outflux.
      if (-Vi_MB( vi) >= Vi_available) then
        ! All ice in this vertex melts, nothing remains to move around. Rescale outfluxes to zero.
        Vi_MB( vi) = -Vi_available
        rescale_factor = 0._dp
      end if

      ! If the total outflux exceeds the available ice plus SMB plus total influx, rescale outfluxes
      if (Vi_out > Vi_available + Vi_MB( vi)) then
        ! Total outflux plus melt exceeds available ice volume. Rescale outfluxes to correct for this.
        rescale_factor = (Vi_available + Vi_MB( vi)) / Vi_out
      end if

      ! Rescale ice outfluxes out of vi and into vi's neighbours
      if (rescale_factor < 1._dp) then
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)

          if (ice%dVi_in( vi,ci) < 0._dp) then
            ice%dVi_in( vi,ci) = ice%dVi_in( vi,ci) * rescale_factor

            do cji = 1, mesh%nC( vj)
              if (mesh%C( vj,cji) == vi) then
                ice%dVi_in( vj,cji) = -ice%dVi_in( vi,ci)
                exit
              end if
            end do

          end if ! IF (ice%dVi_in( vi,ci) < 0._dp) THEN
        end do ! DO ci = 1, mesh%nC( vi)
      end if ! IF (rescale_factor < 1._dp) THEN

    end do ! vi = mesh%vi1, mesh%vi2

    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================

    do vi = mesh%vi1, mesh%vi2
      dVi  = sum( ice%dVi_in( vi,:))
      ice%dHi_dt_a(     vi) = (dVi + Vi_MB( vi)) / (mesh%A( vi) * dt)
      ice%Hi_tplusdt_a( vi) = max( 0._dp, ice%Hi_a( vi) + ice%dHi_dt_a( vi) * dt)
      ice%dHi_dt_a(     vi) = (ice%Hi_tplusdt_a( vi) - ice%Hi_a( vi)) / dt
    end do

    ! === Finalisation ===
    ! ====================

    ! Clean up after yourself
    deallocate( u_c )
    deallocate( v_c )
    deallocate( up_c)
    deallocate( uo_c)
    deallocate( Hi_a)
    deallocate( Cw  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dHi_dt_explicit

! ===== Boundary conditions =====
! ===============================

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_ice_thickness_BC

END MODULE ice_thickness_module
