MODULE ice_thickness_module

  ! Contains the routines for solving the ice thickness equation

  ! Import basic functionality
  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, write_to_memory_log, &
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
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  
  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_ice_model, type_SMB_model, type_BMB_model, type_PD_data_fields
  USE utilities_module,                ONLY: is_floating
  USE mesh_help_functions_module,      ONLY: rotate_xy_to_po_stag

  IMPLICIT NONE
  
CONTAINS
  
  ! The main routine that is called from "run_ice_model" in the ice_dynamics_module
  SUBROUTINE calc_dHi_dt( mesh, ice, SMB, BMB, dt, mask_noice, PD)
    ! Use the total ice velocities to update the ice thickness
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_SMB_model),                INTENT(IN)    :: SMB
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB
    REAL(dp),                            INTENT(IN)    :: dt
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice
    TYPE(type_PD_data_fields),           INTENT(IN)    :: PD
    
    ! Use the specified time integration method to calculate the ice thickness at t+dt
    IF     (C%choice_ice_integration_method == 'none') THEN
      ice%dHi_dt_a( mesh%vi1:mesh%vi2) = 0._dp
      CALL sync
    ELSEIF (C%choice_ice_integration_method == 'explicit') THEN
      CALL calc_dHi_dt_explicit(     mesh, ice, SMB, BMB, dt)
    ELSEIF (C%choice_ice_integration_method == 'semi-implicit') THEN
      WRITE(0,*) 'calc_dHi_dt - calc_dHi_dt_semiimplicit: FIXME!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      !CALL calc_dHi_dt_semiimplicit( mesh, ice, SMB, BMB, dt)
    ELSE
      IF (par%master) WRITE(0,*) 'calc_dHi_dt - ERROR: unknown choice_ice_integration_method "', TRIM(C%choice_ice_integration_method), '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Apply boundary conditions
    CALL apply_ice_thickness_BC( mesh, ice, dt, mask_noice, PD)
    
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
    
    ! Local variables
    REAL(dp), DIMENSION(:    ), POINTER                ::  up_c,  uo_c
    INTEGER                                            :: wup_c, wuo_c
    INTEGER                                            :: aci, vi, vj, cii, ci, cji, cj
    REAL(dp)                                           :: dVi, Vi_out, Vi_in, Vi_available, rescale_factor
    REAL(dp), DIMENSION(mesh%nV)                       :: Vi_SMB
            
    ! Initialise at zero
    ice%dVi_in(  mesh%vi1:mesh%vi2, :) = 0._dp
    ice%dVi_out( mesh%vi1:mesh%vi2, :) = 0._dp
    CALL sync
    
    Vi_in          = 0._dp
    Vi_available   = 0._dp
    rescale_factor = 0._dp
    Vi_SMB         = 0._dp
    
    ! Calculate vertically averaged ice velocities along vertex connections
    CALL allocate_shared_dp_1D( mesh%nAc, up_c, wup_c)
    CALL allocate_shared_dp_1D( mesh%nAc, uo_c, wuo_c)
    CALL rotate_xy_to_po_stag( mesh, ice%u_vav_c, ice%v_vav_c, up_c, uo_c)
        
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
    
    Vi_SMB( mesh%vi1:mesh%vi2) = (SMB%SMB_year( mesh%vi1:mesh%vi2) + BMB%BMB( mesh%vi1:mesh%vi2))  * mesh%A( mesh%vi1:mesh%vi2) * dt! * ice_density / 1000._dp     ! m3 ice equivalent
    CALL sync
    
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
    CALL sync
    
    ! Calculate change in ice thickness over time at every vertex
    ! ===========================================================
      
    DO vi = mesh%vi1, mesh%vi2
      dVi  = SUM( ice%dVi_in( vi,:))
      ice%dHi_dt_a(     vi) = (dVi + Vi_SMB( vi)) / (mesh%A( vi) * dt)  
      ice%Hi_tplusdt_a( vi) = MAX( 0._dp, ice%Hi_a( vi) + ice%dHi_dt_a( vi))
      ice%dHi_dt_a(     vi) = (ice%Hi_tplusdt_a( vi) - ice%Hi_a( vi)) / dt
    END DO    
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wup_c)
    CALL deallocate_shared( wuo_c)
    
  END SUBROUTINE calc_dHi_dt_explicit
    
  ! Some useful tools
  SUBROUTINE apply_ice_thickness_BC( mesh, ice, dt, mask_noice, PD)
    ! Apply ice thickness boundary conditions (at the domain boundary, and through the mask_noice)
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    REAL(dp),                            INTENT(IN)    :: dt
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_noice
    TYPE(type_PD_data_fields),           INTENT(IN)    :: PD
    
    ! Local variables:
    INTEGER                                            :: vi
    
    DO vi = mesh%vi1, mesh%vi2
      
      ! West
      IF (mesh%edge_index( vi) == 6 .OR. mesh%edge_index( vi) == 7 .OR. mesh%edge_index( vi) == 8) THEN
        IF     (C%ice_thickness_west_BC == 'zero') THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        ELSEIF (C%ice_thickness_west_BC == 'infinite') THEN
          WRITE(0,*) 'apply_ice_thickness_BC - ice_thickness_west_BC: FIXME!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSEIF (C%ice_thickness_west_BC == 'periodic') THEN
          WRITE(0,*) 'apply_ice_thickness_BC - ice_thickness_west_BC: FIXME!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSEIF (C%ice_thickness_west_BC == 'ISMIP_HOM_F') THEN
          ice%dHi_dt_a(     vi) = (1000._dp - ice%Hi_a( vi)) / dt
          ice%Hi_tplusdt_a( vi) = 1000._dp
        ELSE
          IF (par%master) WRITE(0,*) 'apply_ice_thickness_BC - ERROR: unknown ice_thickness_west_BC "', TRIM(C%ice_thickness_west_BC), '"!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END IF
        
      ! East
      IF (mesh%edge_index( vi) == 2 .OR. mesh%edge_index( vi) == 3 .OR. mesh%edge_index( vi) == 4) THEN
        IF     (C%ice_thickness_east_BC == 'zero') THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        ELSEIF (C%ice_thickness_east_BC == 'infinite') THEN
          WRITE(0,*) 'apply_ice_thickness_BC - ice_thickness_east_BC: FIXME!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSEIF (C%ice_thickness_east_BC == 'periodic') THEN
          WRITE(0,*) 'apply_ice_thickness_BC - ice_thickness_east_BC: FIXME!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSEIF (C%ice_thickness_east_BC == 'ISMIP_HOM_F') THEN
          ice%dHi_dt_a(     vi) = (1000._dp - ice%Hi_a( vi)) / dt
          ice%Hi_tplusdt_a( vi) = 1000._dp
        ELSE
          IF (par%master) WRITE(0,*) 'apply_ice_thickness_BC - ERROR: unknown ice_thickness_east_BC "', TRIM(C%ice_thickness_east_BC), '"!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END IF
        
      ! South
      IF (mesh%edge_index( vi) == 4 .OR. mesh%edge_index( vi) == 5 .OR. mesh%edge_index( vi) == 6) THEN
        IF     (C%ice_thickness_south_BC == 'zero') THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        ELSEIF (C%ice_thickness_south_BC == 'infinite') THEN
          WRITE(0,*) 'apply_ice_thickness_BC - ice_thickness_south_BC: FIXME!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSEIF (C%ice_thickness_south_BC == 'periodic') THEN
          WRITE(0,*) 'apply_ice_thickness_BC - ice_thickness_south_BC: FIXME!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSEIF (C%ice_thickness_south_BC == 'ISMIP_HOM_F') THEN
          ice%dHi_dt_a(     vi) = (1000._dp - ice%Hi_a( vi)) / dt
          ice%Hi_tplusdt_a( vi) = 1000._dp
        ELSE
          IF (par%master) WRITE(0,*) 'apply_ice_thickness_BC - ERROR: unknown ice_thickness_south_BC "', TRIM(C%ice_thickness_south_BC), '"!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END IF
        
      ! North
      IF (mesh%edge_index( vi) == 8 .OR. mesh%edge_index( vi) == 1 .OR. mesh%edge_index( vi) == 2) THEN
        IF     (C%ice_thickness_north_BC == 'zero') THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        ELSEIF (C%ice_thickness_north_BC == 'infinite') THEN
          WRITE(0,*) 'apply_ice_thickness_BC - ice_thickness_north_BC: FIXME!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSEIF (C%ice_thickness_north_BC == 'periodic') THEN
          WRITE(0,*) 'apply_ice_thickness_BC - ice_thickness_north_BC: FIXME!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        ELSEIF (C%ice_thickness_north_BC == 'ISMIP_HOM_F') THEN
          ice%dHi_dt_a(     vi) = (1000._dp - ice%Hi_a( vi)) / dt
          ice%Hi_tplusdt_a( vi) = 1000._dp
        ELSE
          IF (par%master) WRITE(0,*) 'apply_ice_thickness_BC - ERROR: unknown ice_thickness_north_BC "', TRIM(C%ice_thickness_north_BC), '"!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END IF
      
    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync
    
    ! Remove ice in areas where no ice is allowed (e.g. Greenland in NAM and EAS, and Ellesmere Island in GRL)
    DO vi = mesh%vi1, mesh%vi2
      IF (mask_noice(     vi) == 1) THEN
        ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
        ice%Hi_tplusdt_a( vi) = 0._dp
      END IF
    END DO
    CALL sync
    
    ! If so specified, remove all floating ice
    IF (C%do_remove_shelves) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (is_floating( ice%Hi_tplusdt_a( vi), ice%Hb_a( vi), ice%SL_a( vi))) THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        END IF
      END DO
      CALL sync
    END IF ! IF (C%do_remove_shelves) THEN
    
    ! If so specified, remove all floating ice beyond the present-day calving front
    IF (C%remove_shelves_larger_than_PD) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (PD%Hi( vi) == 0._dp .AND. PD%Hb( vi) < 0._dp) THEN
          ice%dHi_dt_a(     vi) = -ice%Hi_a( vi) / dt
          ice%Hi_tplusdt_a( vi) = 0._dp
        END IF
      END DO
      CALL sync
    END IF ! IF (C%remove_shelves_larger_than_PD) THEN
    
    ! If so specified, remove all floating ice crossing the continental shelf edge
    IF (C%continental_shelf_calving) THEN
      WRITE(0,*) 'apply_ice_thickness_BC - continental_shelf_calving: FIXME!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
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
    
  END SUBROUTINE apply_ice_thickness_BC
  
END MODULE ice_thickness_module
