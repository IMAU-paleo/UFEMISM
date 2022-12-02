MODULE ice_velocity_main_module

! == Contains all the routines needed to calculate instantaneous ice velocities for the current modelled ice-sheet geometry.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE petscksp
  USE configuration_module,                ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                        ONLY: perr, petscmat_checksum
  USE parallel_module,                     ONLY: par, sync, ierr, cerr, partition_list, &
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
  USE utilities_module,                    ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                                 check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D, &
                                                 checksum_dp_1D, checksum_dp_2D, checksum_dp_3D
  USE netcdf_debug_module,                 ONLY: debug, write_to_debug_file, &
                                                 save_variable_as_netcdf_int_1D, save_variable_as_netcdf_dp_1D, &
                                                 save_variable_as_netcdf_int_2D, save_variable_as_netcdf_dp_2D, &
                                                 save_variable_as_netcdf_int_3D, save_variable_as_netcdf_dp_3D, &
                                                 write_PETSc_matrix_to_NetCDF, write_CSR_matrix_to_NetCDF

  ! Import specific functionality
  USE data_types_module,                   ONLY: type_mesh, type_ice_model, type_velocity_solver_SIA, type_velocity_solver_SSA, &
                                                 type_velocity_solver_DIVA, type_velocity_solver_BPA, type_BMB_model
  USE ice_velocity_SIA_module,             ONLY: initialise_SIA_solver , solve_SIA , remap_SIA_solver
  USE ice_velocity_SSA_module,             ONLY: initialise_SSA_solver , solve_SSA , remap_SSA_solver
  USE ice_velocity_DIVA_module,            ONLY: initialise_DIVA_solver, solve_DIVA, remap_DIVA_solver
  USE ice_velocity_BPA_module,             ONLY: initialise_BPA_solver , solve_BPA , remap_BPA_solver
  USE utilities_module,                    ONLY: vertical_average
  USE mesh_operators_module,               ONLY: map_b_to_a_2D, map_b_to_a_3D, ddx_a_to_a_2D, ddy_a_to_a_2D, &
                                                 field2vec_ak, field2vec_bk, vec2field_ak, vec2field_bk
  USE petsc_module,                        ONLY: multiply_PETSc_matrix_with_vector_1D
  USE mesh_help_functions_module,          ONLY: rotate_xy_to_po_stag_3D

  IMPLICIT NONE

CONTAINS

! == The main routines, to be called from the ice dynamics module

  SUBROUTINE initialise_velocity_solver( mesh, ice)
    ! Initialise the velocity solver for the chosen stress balance approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_velocity_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      CALL initialise_SIA_solver(  mesh, ice%SIA)
    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      CALL initialise_SSA_solver(  mesh, ice%SSA)
    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      CALL initialise_SIA_solver(  mesh, ice%SIA)
      CALL initialise_SSA_solver(  mesh, ice%SSA)
    ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
      CALL initialise_DIVA_solver( mesh, ice%DIVA)
    ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
      CALL initialise_BPA_solver(  mesh, ice%BPA)
    ELSE
      CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_velocity_solver

  SUBROUTINE solve_stress_balance( mesh, ice)
    ! Calculate all ice velocities based on the chosen stress balance approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'solve_stress_balance'
    INTEGER                                            :: vi,ti,k

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      ! Calculate velocities according to the Shallow Ice Approximation

      CALL solve_SIA( mesh, ice, ice%SIA)
      CALL set_ice_velocities_to_SIA_results( mesh, ice, ice%SIA)

    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      ! Calculate velocities according to the Shallow Shelf Approximation

      CALL solve_SSA( mesh, ice, ice%SSA)
      CALL set_ice_velocities_to_SSA_results( mesh, ice, ice%SSA)

    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      ! Calculate velocities according to the hybrid SIA/SSA

      CALL solve_SIA( mesh, ice, ice%SIA)
      CALL solve_SSA( mesh, ice, ice%SSA)
      CALL set_ice_velocities_to_SIASSA_results( mesh, ice, ice%SIA, ice%SSA)

    ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
      ! Calculate velocities according to the Depth-Integrated Viscosity Approximation

      CALL solve_DIVA( mesh, ice, ice%DIVA)
      CALL set_ice_velocities_to_DIVA_results( mesh, ice, ice%DIVA)

    ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
      ! Calculate velocities according to the Depth-Integrated Viscosity Approximation

      CALL solve_BPA( mesh, ice, ice%BPA)
      CALL set_ice_velocities_to_BPA_results( mesh, ice, ice%BPA)

    ELSE
      CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    END IF

  ! == Fill in derived velocity fields (surface, base, vertical average)

    DO ti = mesh%ti1, mesh%ti2

      ! Surface
      ice%u_surf_b(    ti) = ice%u_3D_b( ti,1)
      ice%v_surf_b(    ti) = ice%v_3D_b( ti,1)
      ice%uabs_surf_b( ti) = SQRT( ice%u_surf_b( ti)**2 + ice%v_surf_b( ti)**2)

      ! Base
      ice%u_base_b(    ti) = ice%u_3D_b( ti,C%nz)
      ice%v_base_b(    ti) = ice%v_3D_b( ti,C%nz)
      ice%uabs_base_b( ti) = SQRT( ice%u_base_b( ti)**2 + ice%v_base_b( ti)**2)

      ! Vertical average
      CALL vertical_average( C%zeta, ice%u_3D_b( ti,:), ice%u_vav_b( ti))
      CALL vertical_average( C%zeta, ice%v_3D_b( ti,:), ice%v_vav_b( ti))
      ice%uabs_vav_b( ti) = SQRT( ice%u_vav_b( ti)**2 + ice%v_vav_b( ti)**2)

    END DO
    CALL sync

  ! == Calculate velocities on the a-grid (needed to calculate the vertical velocity w, and for writing to output)

    ! 3-D
    CALL map_b_to_a_3D( mesh, ice%u_3D_b  , ice%u_3D_a  )
    CALL map_b_to_a_3D( mesh, ice%v_3D_b  , ice%v_3D_a  )

    ! Surface
    CALL map_b_to_a_2D( mesh, ice%u_surf_b, ice%u_surf_a)
    CALL map_b_to_a_2D( mesh, ice%v_surf_b, ice%v_surf_a)

    ! Base
    CALL map_b_to_a_2D( mesh, ice%u_base_b, ice%u_base_a)
    CALL map_b_to_a_2D( mesh, ice%v_base_b, ice%v_base_a)

    ! Vertical average
    CALL map_b_to_a_2D( mesh, ice%u_vav_b , ice%u_vav_a )
    CALL map_b_to_a_2D( mesh, ice%v_vav_b , ice%v_vav_a )

    ! Absolute
    DO vi = mesh%vi1, mesh%vi2
      ice%uabs_surf_a( vi) = SQRT( ice%u_surf_a( vi)**2 + ice%v_surf_a( vi)**2)
      ice%uabs_base_a( vi) = SQRT( ice%u_base_a( vi)**2 + ice%v_base_a( vi)**2)
      ice%uabs_vav_a(  vi) = SQRT( ice%u_vav_a(  vi)**2 + ice%v_vav_a(  vi)**2)
    END DO

  ! == Calculate velocities on the c-grid (needed to solve the ice thickness equation)

    CALL map_velocities_from_b_to_c_2D( mesh, ice%u_vav_b, ice%v_vav_b, ice%u_vav_c, ice%v_vav_c)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE solve_stress_balance

  SUBROUTINE calc_vertical_velocities( mesh, ice, BMB)
    ! Calculate vertical velocity w from conservation of mass
    !
    ! NOTE: since the vertical velocities for floating ice depend on
    !       the thinning rate dH/dt, this routine must be called
    !       after having calculated dHi_dt!

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_vertical_velocities'
    INTEGER                                            :: vi,k,ci,vj,aci
    REAL(dp), DIMENSION(:    ), POINTER                :: H_base_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dt_a
    INTEGER                                            :: wH_base_a, wdH_base_dx_a, wdH_base_dy_a, wdH_base_dt_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: u_3D_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: v_3D_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: p_3D_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: o_3D_c
    INTEGER                                            :: wu_3D_c, wv_3D_c, wp_3D_c, wo_3D_c
    REAL(dp)                                           :: dzeta, up, Q_hor, w_base, Q_base, Q_top

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV  ,         H_base_a       , wH_base_a       )
    CALL allocate_shared_dp_1D( mesh%nV  ,         dH_base_dx_a   , wdH_base_dx_a   )
    CALL allocate_shared_dp_1D( mesh%nV  ,         dH_base_dy_a   , wdH_base_dy_a   )
    CALL allocate_shared_dp_1D( mesh%nV  ,         dH_base_dt_a   , wdH_base_dt_a   )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , u_3D_c         , wu_3D_c         )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , v_3D_c         , wv_3D_c         )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , p_3D_c         , wp_3D_c         )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , o_3D_c         , wo_3D_c         )

    DO vi = mesh%vi1, mesh%vi2

      ! Calculate elevation of the ice base
      H_base_a( vi) = ice%Hs_a( vi) - ice%Hi_a( vi)

      ! Calculate rate of change of ice base elevation
      IF     (ice%mask_sheet_a( vi) == 1) THEN
        ! For grounded ice, the ice base simply moves with the bedrock
        dH_base_dt_a( vi) =  ice%dHb_dt_a( vi)
      ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
        ! For floating ice, the ice base moves according to the thinning rate times the density fraction
        dH_base_dt_a( vi) = -ice%dHi_dt_a( vi) * ice_density / seawater_density
      ELSE
        ! No ice, so no vertical velocity
        dH_base_dt_a( vi) = 0._dp
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Calculate slopes of the ice base
    CALL ddx_a_to_a_2D( mesh, H_base_a, dH_base_dx_a)
    CALL ddy_a_to_a_2D( mesh, H_base_a, dH_base_dy_a)

    ! Calculate 3-D ice velocities on the c-grid
    CALL map_velocities_from_b_to_c_3D( mesh, ice%u_3D_b, ice%v_3D_b, u_3D_c, v_3D_c)

    ! Calculate 3-D velocity components parallel to mesh edges on the c-grid
    CALL rotate_xy_to_po_stag_3D( mesh, u_3D_c, v_3D_c, p_3D_c, o_3D_c)

    ! Calculate vertical velocities by solving conservation of mass in each 3-D cell
    DO vi = mesh%vi1, mesh%vi2

      ! No ice means no velocity
      IF (ice%mask_ice_a( vi) == 0) THEN
        ice%w_3D_a( vi,:) = 0._dp
        CYCLE
      END IF

      ! Calculate the vertical velocity at the ice base, which is equal to
      ! the horizontal motion along the sloping ice base, plus the vertical
      ! motion of the ice base itself, plus the vertical motion of an ice
      ! particle with respect to the ice base (i.e. the basal melt rate).
      !
      ! NOTE: BMB is defined so that a positive number means accumulation of ice;
      !       at the ice base, that means that a positive BMB means a positive
      !       value of w

      ice%w_3D_a( vi,C%nz) = (ice%u_3D_a( vi,C%nz) * dH_base_dx_a( vi)) + &
                             (ice%v_3D_a( vi,C%nz) * dH_base_dy_a( vi)) + &
                              dH_base_dt_a( vi) + BMB%BMB( vi)

      DO k = C%nz-1, 1, -1

        dzeta = C%zeta( k) - C%zeta( k+1)

        ! Calculate net horizontal ice flux into this 3-D grid box
        Q_hor = 0._dp
        DO ci = 1, mesh%nC( vi)
          vj  = mesh%C(    vi,ci)
          aci = mesh%iAci( vi,ci)
          up  = 0.5_dp * (p_3D_c( aci,k) + p_3D_c( aci,k+1))
          IF (p_3D_c( aci,k) > 0._dp) THEN
            ! Ice moves from vi to vj
            Q_hor = Q_hor - up * dzeta * ice%Hi_a( vi) * mesh%Cw( vi,ci)   ! m^3/yr
          ELSE
            ! Ice moves from vj to vi
            Q_hor = Q_hor + up * dzeta * ice%Hi_a( vj) * mesh%Cw( vi,ci)   ! m^3/yr
          END IF
        END DO

        ! Calculate vertical ice flux at the base of this 3-D grid box
        Q_base = ice%w_3D_a( vi,k+1) * mesh%A( vi)

        ! Calculate vertical ice flux through the top of this 3-D grid box,
        ! such that conservation of mass is satisfied
        Q_top = -(Q_hor + Q_base)

        ! Calculate vertical ice velocity at the top of this 3-D grid box
        ice%w_3D_a( vi,k) = Q_top / mesh%A( vi)

      END DO ! DO k = C%nz-1, 1, -1

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Clean up after yourself
    CALL deallocate_shared( wH_base_a       )
    CALL deallocate_shared( wdH_base_dx_a   )
    CALL deallocate_shared( wdH_base_dy_a   )
    CALL deallocate_shared( wdH_base_dt_a   )
    CALL deallocate_shared( wu_3D_c         )
    CALL deallocate_shared( wv_3D_c         )
    CALL deallocate_shared( wp_3D_c         )
    CALL deallocate_shared( wo_3D_c         )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_velocities

  SUBROUTINE calc_vertical_velocities2( mesh, ice, BMB)
    ! Calculate vertical velocity w from conservation of mass
    !
    ! NOTE: since the vertical velocities for floating ice depend on
    !       the thinning rate dH/dt, this routine must be called
    !       after having calculated dHi_dt!

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_vertical_velocities'
    INTEGER                                            :: vi,k,ci,vj,aci
    REAL(dp), DIMENSION(:    ), POINTER                :: H_base_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dt_a
    INTEGER                                            :: wH_base_a, wdH_base_dx_a, wdH_base_dy_a, wdH_base_dt_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: u_3D_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: v_3D_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: p_3D_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: o_3D_c
    INTEGER                                            :: wu_3D_c, wv_3D_c, wp_3D_c, wo_3D_c
    REAL(dp), DIMENSION(:    ), POINTER                :: w_base_a
    REAL(dp), DIMENSION(:    ), POINTER                :: w_surf_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: w_3D_aks
    INTEGER                                            :: ww_base_a, ww_surf_a, ww_3D_aks
    REAL(dp)                                           :: zeta_top, zeta_bottom, dzeta
    REAL(dp)                                           :: Q_hor, w_base, Q_base, Q_top, w_top

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV  ,         H_base_a       , wH_base_a       )
    CALL allocate_shared_dp_1D( mesh%nV  ,         dH_base_dx_a   , wdH_base_dx_a   )
    CALL allocate_shared_dp_1D( mesh%nV  ,         dH_base_dy_a   , wdH_base_dy_a   )
    CALL allocate_shared_dp_1D( mesh%nV  ,         dH_base_dt_a   , wdH_base_dt_a   )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , u_3D_c         , wu_3D_c         )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , v_3D_c         , wv_3D_c         )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , p_3D_c         , wp_3D_c         )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , o_3D_c         , wo_3D_c         )
    CALL allocate_shared_dp_1D( mesh%nV  ,         w_base_a       , ww_base_a       )
    CALL allocate_shared_dp_1D( mesh%nV  ,         w_surf_a       , ww_surf_a       )
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz-1, w_3D_aks       , ww_3D_aks       )

    DO vi = mesh%vi1, mesh%vi2

      ! Calculate elevation of the ice base
      H_base_a( vi) = ice%Hs_a( vi) - ice%Hi_a( vi)

      ! Calculate rate of change of ice base elevation
      IF     (ice%mask_sheet_a( vi) == 1) THEN
        ! For grounded ice, the ice base simply moves with the bedrock
        dH_base_dt_a( vi) =  ice%dHb_dt_a( vi)
      ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
        ! For floating ice, the ice base moves according to the thinning rate times the density fraction
        dH_base_dt_a( vi) = -ice%dHi_dt_a( vi) * ice_density / seawater_density
      ELSE
        ! No ice, so no vertical velocity
        dH_base_dt_a( vi) = 0._dp
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Calculate slopes of the ice base
    CALL ddx_a_to_a_2D( mesh, H_base_a, dH_base_dx_a)
    CALL ddy_a_to_a_2D( mesh, H_base_a, dH_base_dy_a)

    ! Calculate 3-D ice velocities on the c-grid
    CALL map_velocities_from_b_to_c_3D( mesh, ice%u_3D_b, ice%v_3D_b, u_3D_c, v_3D_c)

    ! Calculate 3-D velocity components parallel to mesh edges on the c-grid
    CALL rotate_xy_to_po_stag_3D( mesh, u_3D_c, v_3D_c, p_3D_c, o_3D_c)

    ! Calculate vertical velocities by solving conservation of mass in each 3-D cell
    DO vi = mesh%vi1, mesh%vi2

      ! No ice means no velocity
      IF (ice%mask_ice_a( vi) == 0) THEN
        ice%w_3D_a( vi,:) = 0._dp
        CYCLE
      END IF

      ! Calculate the vertical velocity at the ice base, which is equal to
      ! the horizontal motion along the sloping ice base, plus the vertical
      ! motion of the ice base itself, plus the vertical motion of an ice
      ! particle with respect to the ice base (i.e. the basal melt rate).
      !
      ! NOTE: BMB is defined so that a positive number means accumulation of ice;
      !       at the ice base, that means that a positive BMB means a positive
      !       value of w

      w_base_a( vi) = (ice%u_3D_a( vi,C%nz) * dH_base_dx_a( vi)) + &
                      (ice%v_3D_a( vi,C%nz) * dH_base_dy_a( vi)) + &
                       dH_base_dt_a( vi) + BMB%BMB( vi)

      DO k = C%nz-1, 1, -1

        ! Height of this grid cell
        IF     (k == 1) THEN
          ! Ice surface
          zeta_top    = 0._dp
          zeta_bottom = (C%zeta( k+1) + C%zeta( k)) / 2._dp
        ELSEIF (k == C%nz) THEN
          ! Ice base
          zeta_top    = (C%zeta( k-1) + C%zeta( k)) / 2._dp
          zeta_bottom = 1._dp
        ELSE
          ! Ice column
          zeta_top    = (C%zeta( k-1) + C%zeta( k)) / 2._dp
          zeta_bottom = (C%zeta( k+1) + C%zeta( k)) / 2._dp
        END IF

        dzeta = zeta_top - zeta_bottom

        ! Calculate net horizontal ice flux into this 3-D grid cell
        Q_hor = 0._dp
        DO ci = 1, mesh%nC( vi)
          vj  = mesh%C(    vi,ci)
          aci = mesh%iAci( vi,ci)
          IF (p_3D_c( aci,k) > 0._dp) THEN
            ! Ice moves from vi to vj
            Q_hor = Q_hor - p_3D_c( aci,k) * dzeta * ice%Hi_a( vi) * mesh%Cw( vi,ci)   ! m^3/yr
          ELSE
            ! Ice moves from vj to vi
            Q_hor = Q_hor + p_3D_c( aci,k) * dzeta * ice%Hi_a( vj) * mesh%Cw( vi,ci)   ! m^3/yr
          END IF
        END DO

        ! Calculate vertical ice flux at the base of this 3-D grid cell
        IF (k == C%nz) THEN
          ! Ice base
          w_base = w_base_a( vi)
        ELSE
          ! Ice column
          w_base = w_3D_aks( vi,k)
        END IF
        Q_base = w_base * mesh%A( vi)

        ! Calculate vertical ice flux through the top of this 3-D grid cell,
        ! such that conservation of mass is satisfied
        Q_top = -(Q_hor + Q_base)

        ! Calculate vertical ice velocity at the top of this 3-D grid cell
        w_top = Q_top / mesh%A( vi)

        ! Fill into arrays
        IF (k == 1) THEN
          ! Ice surface
          w_surf_a( vi) = w_top
        ELSE
          ! Ice column
          w_3D_aks( vi,k-1) = w_top
        END IF

      END DO ! DO k = C%nz-1, 1, -1

      ! Un-stagger vertical velocities
      ice%w_3D_a( vi,1   ) = w_surf_a( vi)
      ice%w_3D_a( vi,C%nz) = w_base_a( vi)
      DO k = 2, C%nz-1
        ice%w_3D_a( vi,k) = (w_3D_aks( vi,k-1) + w_3D_aks( vi,k)) / 2._dp
      END DO

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Clean up after yourself
    CALL deallocate_shared( wH_base_a       )
    CALL deallocate_shared( wdH_base_dx_a   )
    CALL deallocate_shared( wdH_base_dy_a   )
    CALL deallocate_shared( wdH_base_dt_a   )
    CALL deallocate_shared( wu_3D_c         )
    CALL deallocate_shared( wv_3D_c         )
    CALL deallocate_shared( wp_3D_c         )
    CALL deallocate_shared( wo_3D_c         )
    CALL deallocate_shared( ww_base_a       )
    CALL deallocate_shared( ww_surf_a       )
    CALL deallocate_shared( ww_3D_aks       )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_velocities2

  SUBROUTINE calc_vertical_velocities_old( mesh, ice, BMB)
    ! Calculate vertical velocity w from conservation of mass
    !
    ! NOTE: since the vertical velocities for floating ice depend on
    !       the thinning rate dH/dt, this routine must be called
    !       after having calculated dHi_dt!

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_vertical_velocities'
    INTEGER                                            :: vi,ti,k,n
    REAL(dp), DIMENSION(:    ), POINTER                :: H_base_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dt_a
    REAL(dp), DIMENSION(:    ), POINTER                :: u_bk_vec
    REAL(dp), DIMENSION(:    ), POINTER                :: v_bk_vec
    REAL(dp), DIMENSION(:    ), POINTER                :: du_dxp_ak_vec
    REAL(dp), DIMENSION(:    ), POINTER                :: dv_dyp_ak_vec
    REAL(dp), DIMENSION(:    ), POINTER                :: du_dzeta_ak_vec
    REAL(dp), DIMENSION(:    ), POINTER                :: dv_dzeta_ak_vec
    REAL(dp), DIMENSION(:,:  ), POINTER                :: du_dxp_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: dv_dyp_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: du_dzeta_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: dv_dzeta_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: du_dx_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: dv_dy_3D_a
    REAL(dp), DIMENSION(:,:  ), POINTER                :: dw_dz_3D_a
    INTEGER                                            :: wH_base_a, wdH_base_dx_a, wdH_base_dy_a, wdH_base_dt_a
    INTEGER                                            :: wu_bk_vec, wv_bk_vec, wdu_dxp_ak_vec, wdv_dyp_ak_vec, wdu_dzeta_ak_vec, wdv_dzeta_ak_vec
    INTEGER                                            :: wdu_dxp_3D_a, wdu_dzeta_3D_a, wdv_dyp_3D_a, wdv_dzeta_3D_a, wdu_dx_3D_a, wdv_dy_3D_a, wdw_dz_3D_a
    REAL(dp)                                           :: dw_dz, dz

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV  ,       H_base_a       , wH_base_a       )
    CALL allocate_shared_dp_1D( mesh%nV  ,       dH_base_dx_a   , wdH_base_dx_a   )
    CALL allocate_shared_dp_1D( mesh%nV  ,       dH_base_dy_a   , wdH_base_dy_a   )
    CALL allocate_shared_dp_1D( mesh%nV  ,       dH_base_dt_a   , wdH_base_dt_a   )
    CALL allocate_shared_dp_1D( mesh%nnbk,       u_bk_vec       , wu_bk_vec       )
    CALL allocate_shared_dp_1D( mesh%nnbk,       v_bk_vec       , wv_bk_vec       )
    CALL allocate_shared_dp_1D( mesh%nnak,       du_dxp_ak_vec  , wdu_dxp_ak_vec  )
    CALL allocate_shared_dp_1D( mesh%nnak,       du_dzeta_ak_vec, wdu_dzeta_ak_vec)
    CALL allocate_shared_dp_1D( mesh%nnak,       dv_dyp_ak_vec  , wdv_dyp_ak_vec  )
    CALL allocate_shared_dp_1D( mesh%nnak,       dv_dzeta_ak_vec, wdv_dzeta_ak_vec)
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz, du_dxp_3D_a    , wdu_dxp_3D_a    )
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz, du_dzeta_3D_a  , wdu_dzeta_3D_a  )
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz, dv_dyp_3D_a    , wdv_dyp_3D_a    )
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz, dv_dzeta_3D_a  , wdv_dzeta_3D_a  )
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz, du_dx_3D_a     , wdu_dx_3D_a     )
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz, dv_dy_3D_a     , wdv_dy_3D_a     )
    CALL allocate_shared_dp_2D( mesh%nV  , C%nz, dw_dz_3D_a     , wdw_dz_3D_a     )

    DO vi = mesh%vi1, mesh%vi2

      ! Calculate elevation of the ice base
      H_base_a( vi) = ice%Hs_a( vi) - ice%Hi_a( vi)

      ! Calculate rate of change of ice base elevation
      IF     (ice%mask_sheet_a( vi) == 1) THEN
        ! For grounded ice, the ice base simply moves with the bedrock
        dH_base_dt_a( vi) =  ice%dHb_dt_a( vi)
      ELSEIF (ice%mask_shelf_a( vi) == 1) THEN
        ! For floating ice, the ice base moves according to the thinning rate times the density fraction
        dH_base_dt_a( vi) = -ice%dHi_dt_a( vi) * ice_density / seawater_density
      ELSE
        ! No ice, so no vertical velocity
        dH_base_dt_a( vi) = 0._dp
      END IF

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Calculate slopes of the ice base
    CALL ddx_a_to_a_2D( mesh, H_base_a, dH_base_dx_a)
    CALL ddy_a_to_a_2D( mesh, H_base_a, dH_base_dy_a)

    ! Transform 3-D horizontal velocities from field form to vector form
    DO ti = mesh%ti1, mesh%ti2
    DO k  = 1, C%nz
      n = mesh%tik2n( ti,k)
      u_bk_vec( n) = ice%u_3D_b( ti,k)
      v_bk_vec( n) = ice%v_3D_b( ti,k)
    END DO
    END DO

    ! Calculate du/dxp, du/dzeta, dv/dyp, and dv/dzeta in vector form
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddxp_bk_ak  , u_bk_vec, du_dxp_ak_vec  )
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddzeta_bk_ak, u_bk_vec, du_dzeta_ak_vec)
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddyp_bk_ak  , v_bk_vec, dv_dyp_ak_vec  )
    CALL multiply_PETSc_matrix_with_vector_1D( mesh%M_ddzeta_bk_ak, v_bk_vec, dv_dzeta_ak_vec)

    ! Convert du/dxp, du/dzeta, dv/dyp, and dv/dzeta from vector form to field form
    DO vi = mesh%vi1, mesh%vi2
    DO k  = 1, C%nz
      n = mesh%vik2n( vi,k)
      du_dxp_3D_a(   vi,k) = du_dxp_ak_vec(   n)
      du_dzeta_3D_a( vi,k) = du_dzeta_ak_vec( n)
      dv_dyp_3D_a(   vi,k) = dv_dyp_ak_vec(   n)
      dv_dzeta_3D_a( vi,k) = dv_dzeta_ak_vec( n)
    END DO
    END DO

    ! Calculate du/dx, dv/dy, dw/dz
    DO vi = mesh%vi1, mesh%vi2
    DO k  = 1, C%nz

      ! Coordinate transformation
      du_dx_3D_a( vi,k) = du_dxp_3D_a( vi,k) + ice%dzeta_dx_ak( vi,k) * du_dzeta_3D_a( vi,k)
      dv_dy_3D_a( vi,k) = dv_dyp_3D_a( vi,k) + ice%dzeta_dy_ak( vi,k) * dv_dzeta_3D_a( vi,k)

      ! Conservation of mass: du/dx + dv/dy + dw/dz = 0
      dw_dz_3D_a( vi,k) = -1._dp * (du_dx_3D_a( vi,k) + dv_dy_3D_a( vi,k))

    END DO
    END DO

    ! Integrate over the incompressibility condition to find the vertical velocity profiles
    DO vi = mesh%vi1, mesh%vi2

      ! Exception for ice-free grid cells
      IF (ice%mask_ice_a( vi) == 0) THEN
        ice%w_3D_a( vi,:) = 0._dp
        CYCLE
      END IF

      ! Calculate the vertical velocity at the ice base, which is equal to
      ! the horizontal motion along the sloping ice base, plus the vertical
      ! motion of the ice base itself, plus the vertical motion of an ice
      ! particle with respect to the ice base (i.e. the basal melt rate).
      !
      ! NOTE: BMB is defined so that a positive number means accumulation of ice;
      !       at the ice base, that means that a positive BMB means a positive
      !       value of w

      ice%w_3D_a( vi,C%nz) = (ice%u_3D_a( vi,C%nz) * dH_base_dx_a( vi)) + &
                             (ice%v_3D_a( vi,C%nz) * dH_base_dy_a( vi)) + &
                             dH_base_dt_a( vi) + BMB%BMB( vi)

      ! Conservation of mass: du/dx + dv/dy + dw/dz = 0
      DO k = C%nz-1, 1, -1

        ! w in the middle of the current layer (= w( k))
        ! equals:
        ! w in the middle of the layer below this one ( = w( k+1))
        ! plus
        ! dw/dz halfway between these layers (= dw/dz( k+1/2)) times the distance dz between these layers

        ! Let dw/dz( k+1/2) = (dw/dz( k+1) + dw/dz( k)) / 2
        dw_dz = (dw_dz_3D_a( vi,k+1) + dw_dz_3D_a( vi,k)) / 2._dp

        ! Distance dz between these layers
        dz = ice%Hi_a( vi) * (C%zeta( k+1) - C%zeta( k))

        ! Calculate w in the current layer
        ice%w_3D_a( vi,k) = ice%w_3D_a( vi,k+1) + dw_dz * dz

      END DO ! DO k = C%nz-1, -1, 1

    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync

    ! DENK DROM
    IF (MAXVAL( ice%Hi_a) > 0.24900E04_dp) THEN
      CALL save_variable_as_netcdf_int_1D( mesh%vi2n, 'vi2n')
      CALL save_variable_as_netcdf_int_1D( mesh%ti2n, 'ti2n')
      CALL save_variable_as_netcdf_int_2D( mesh%n2vi, 'n2vi')
      CALL save_variable_as_netcdf_int_2D( mesh%vik2n, 'vik2n')
      CALL save_variable_as_netcdf_int_2D( mesh%tik2n, 'tik2n')

      CALL write_PETSc_matrix_to_NetCDF( mesh%M_map_b_a, 'M_map_b_a')
      CALL write_PETSc_matrix_to_NetCDF( mesh%M_ddx_b_a, 'M_ddx_b_a')
      CALL write_PETSc_matrix_to_NetCDF( mesh%M_ddy_b_a, 'M_ddy_b_a')
      CALL write_PETSc_matrix_to_NetCDF( mesh%M_map_bk_ak, 'M_map_bk_ak')
      CALL write_PETSc_matrix_to_NetCDF( mesh%M_ddxp_bk_ak, 'M_ddxp_bk_ak')
      CALL write_PETSc_matrix_to_NetCDF( mesh%M_ddyp_bk_ak, 'M_ddyp_bk_ak')
      CALL write_PETSc_matrix_to_NetCDF( mesh%M_ddzeta_bk_ak, 'M_ddzeta_bk_ak')

      CALL save_variable_as_netcdf_dp_1D( C%zeta, 'zeta')
      CALL save_variable_as_netcdf_dp_2D( ice%dzeta_dx_ak, 'dzeta_dx_ak')
      CALL save_variable_as_netcdf_dp_2D( ice%dzeta_dy_ak, 'dzeta_dy_ak')
      CALL save_variable_as_netcdf_dp_1D( ice%Hi_a, 'Hi_a')
      CALL save_variable_as_netcdf_dp_1D( ice%Hb_a, 'Hb_a')
      CALL save_variable_as_netcdf_dp_1D( ice%Hs_a, 'Hs_a')
      CALL save_variable_as_netcdf_dp_1D( H_base_a, 'H_base_a')
      CALL save_variable_as_netcdf_dp_1D( dH_base_dx_a, 'dH_base_dx_a')
      CALL save_variable_as_netcdf_dp_1D( dH_base_dy_a, 'dH_base_dy_a')
      CALL save_variable_as_netcdf_dp_1D( dH_base_dt_a, 'dH_base_dt_a')
      CALL save_variable_as_netcdf_dp_2D( ice%u_3D_b, 'u_3D_b')
      CALL save_variable_as_netcdf_dp_2D( ice%v_3D_b, 'v_3D_b')
      CALL save_variable_as_netcdf_dp_2D( ice%u_3D_a, 'u_3D_a')
      CALL save_variable_as_netcdf_dp_2D( ice%v_3D_a, 'v_3D_a')
      CALL save_variable_as_netcdf_dp_1D( ice%u_vav_a, 'u_vav_a')
      CALL save_variable_as_netcdf_dp_1D( ice%v_vav_a, 'v_vav_a')
    END IF

    ! Extract vertical velocities at the ice surface and base
    DO vi = mesh%vi1, mesh%vi2
      ice%w_surf_a( vi) = ice%w_3D_a( vi,1   )
      ice%w_base_a( vi) = ice%w_3D_a( vi,C%nz)
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wH_base_a       )
    CALL deallocate_shared( wdH_base_dx_a   )
    CALL deallocate_shared( wdH_base_dy_a   )
    CALL deallocate_shared( wdH_base_dt_a   )
    CALL deallocate_shared( wu_bk_vec       )
    CALL deallocate_shared( wv_bk_vec       )
    CALL deallocate_shared( wdu_dxp_ak_vec  )
    CALL deallocate_shared( wdu_dzeta_ak_vec)
    CALL deallocate_shared( wdv_dyp_ak_vec  )
    CALL deallocate_shared( wdv_dzeta_ak_vec)
    CALL deallocate_shared( wdu_dxp_3D_a    )
    CALL deallocate_shared( wdu_dzeta_3D_a  )
    CALL deallocate_shared( wdv_dyp_3D_a    )
    CALL deallocate_shared( wdv_dzeta_3D_a  )
    CALL deallocate_shared( wdu_dx_3D_a     )
    CALL deallocate_shared( wdv_dy_3D_a     )
    CALL deallocate_shared( wdw_dz_3D_a     )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_velocities_old

  SUBROUTINE remap_velocity_solver( mesh_old, mesh_new, ice)
    ! Remap the velocity solver for the chosen stress balance approximation

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_old
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh_new
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remap_velocity_solver'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_stress_balance_approximation == 'none') THEN
      ! No need to do anything
    ELSEIF (C%choice_stress_balance_approximation == 'SIA') THEN
      CALL remap_SIA_solver(  mesh_old, mesh_new, ice, ice%SIA)
    ELSEIF (C%choice_stress_balance_approximation == 'SSA') THEN
      CALL remap_SSA_solver(  mesh_old, mesh_new, ice, ice%SSA)
    ELSEIF (C%choice_stress_balance_approximation == 'SIA/SSA') THEN
      CALL remap_SIA_solver(  mesh_old, mesh_new, ice, ice%SIA)
      CALL remap_SSA_solver(  mesh_old, mesh_new, ice, ice%SSA)
    ELSEIF (C%choice_stress_balance_approximation == 'DIVA') THEN
      CALL remap_DIVA_solver( mesh_old, mesh_new, ice, ice%DIVA)
    ELSEIF (C%choice_stress_balance_approximation == 'BPA') THEN
      CALL remap_BPA_solver(  mesh_old, mesh_new, ice, ice%BPA)
    ELSE
      CALL crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_velocity_solver

! == Set applied ice model velocities to stress balance results

  SUBROUTINE set_ice_velocities_to_SIA_results( mesh, ice, SIA)
    ! Set applied ice model velocities and strain rates to SIA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_velocity_solver_SIA),      INTENT(IN)    :: SIA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_SIA_results'
    INTEGER                                            :: vi,ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Velocities
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = SIA%u_3D_b( ti,:)
      ice%v_3D_b( ti,:) = SIA%v_3D_b( ti,:)
    END DO

    ! Strain rates
    DO vi = mesh%vi1, mesh%vi2
      ice%du_dz_3D_a( vi,:) = SIA%du_dz_3D_a( vi,:)
      ice%dv_dz_3D_a( vi,:) = SIA%dv_dz_3D_a( vi,:)
    END DO
    ! In the SIA, horizontal gradients of u,v, and all gradients of w, are neglected
    ice%du_dx_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%du_dy_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dv_dx_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dv_dy_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dx_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dy_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dz_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_SIA_results

  SUBROUTINE set_ice_velocities_to_SSA_results( mesh, ice, SSA)
    ! Set applied ice model velocities to SSA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_velocity_solver_SSA),      INTENT(IN)    :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_SSA_results'
    INTEGER                                            :: ti,vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Velocities
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = SSA%u_b( ti)
      ice%v_3D_b( ti,:) = SSA%v_b( ti)
    END DO

    ! Strain rates
    DO vi = mesh%vi1, mesh%vi2
      ice%du_dx_3D_a( vi,:) = SSA%du_dx_a( vi)
      ice%du_dy_3D_a( vi,:) = SSA%du_dy_a( vi)
      ice%dv_dx_3D_a( vi,:) = SSA%dv_dx_a( vi)
      ice%dv_dy_3D_a( vi,:) = SSA%dv_dy_a( vi)
    END DO
    ! In the SSA, vertical gradients of u,v, and all gradients of w, are neglected
    ice%du_dz_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dv_dz_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dx_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dy_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dz_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_SSA_results

  SUBROUTINE set_ice_velocities_to_SIASSA_results( mesh, ice, SIA, SSA)
    ! Set applied ice model velocities to hybrid SIA/SSA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_velocity_solver_SIA),      INTENT(IN)    :: SIA
    TYPE(type_velocity_solver_SSA),      INTENT(IN)    :: SSA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_SIASSA_results'
    INTEGER                                            :: ti,vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (C%choice_hybrid_SIASSA_scheme == 'add') THEN
      ! u = u_SIA + u_SSA

      ! Velocities
      DO ti = mesh%ti1, mesh%ti2
        ice%u_3D_b( ti,:) = SIA%u_3D_b( ti,:) + SSA%u_b( ti)
        ice%v_3D_b( ti,:) = SIA%v_3D_b( ti,:) + SSA%v_b( ti)
      END DO

      ! Strain rates
      DO vi = mesh%vi1, mesh%vi2
        ice%du_dz_3D_a( vi,:) = SIA%du_dz_3D_a( vi,:)
        ice%dv_dz_3D_a( vi,:) = SIA%dv_dz_3D_a( vi,:)
        ice%du_dx_3D_a( vi,:) = SSA%du_dx_a(    vi  )
        ice%du_dy_3D_a( vi,:) = SSA%du_dy_a(    vi  )
        ice%dv_dx_3D_a( vi,:) = SSA%dv_dx_a(    vi  )
        ice%dv_dy_3D_a( vi,:) = SSA%dv_dy_a(    vi  )
      END DO
      ! In the hybrid SIA/SSA, gradients of w are neglected
      ice%dw_dx_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
      ice%dw_dy_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
      ice%dw_dz_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
      CALL sync

    ELSE
      CALL crash('unknown choice_hybrid_SIASSA_scheme_config "' // TRIM( C%choice_hybrid_SIASSA_scheme) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_SIASSA_results

  SUBROUTINE set_ice_velocities_to_DIVA_results( mesh, ice, DIVA)
    ! Set applied ice model velocities to DIVA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_velocity_solver_DIVA),     INTENT(IN)    :: DIVA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_DIVA_results'
    INTEGER                                            :: ti,vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Velocities
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = DIVA%u_3D_b( ti,:)
      ice%v_3D_b( ti,:) = DIVA%v_3D_b( ti,:)
    END DO

    ! Strain rates
    DO vi = mesh%vi1, mesh%vi2
      ice%du_dx_3D_a( vi,:) = DIVA%du_dx_a(    vi  )
      ice%du_dy_3D_a( vi,:) = DIVA%du_dy_a(    vi  )
      ice%du_dz_3D_a( vi,:) = DIVA%du_dz_3D_a( vi,:)
      ice%dv_dx_3D_a( vi,:) = DIVA%dv_dx_a(    vi  )
      ice%dv_dy_3D_a( vi,:) = DIVA%dv_dy_a(    vi  )
      ice%dv_dz_3D_a( vi,:) = DIVA%dv_dz_3D_a( vi,:)
    END DO
    ! In the DIVA, gradients of w are neglected
    ice%dw_dx_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dy_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dz_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_DIVA_results

  SUBROUTINE set_ice_velocities_to_BPA_results( mesh, ice, BPA)
    ! Set applied ice model velocities and strain rates to BPA results

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(INOUT) :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_velocity_solver_BPA),      INTENT(IN)    :: BPA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'set_ice_velocities_to_BPA_results'
    INTEGER                                            :: ti

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Velocities
    DO ti = mesh%ti1, mesh%ti2
      ice%u_3D_b( ti,:) = BPA%u_3D_b( ti,:)
      ice%v_3D_b( ti,:) = BPA%v_3D_b( ti,:)
    END DO

    ! Strain rates
    CALL vec2field_ak( mesh, BPA%du_dx_ak_vec, ice%du_dx_3D_a)
    CALL vec2field_ak( mesh, BPA%du_dy_ak_vec, ice%du_dy_3D_a)
    CALL vec2field_ak( mesh, BPA%du_dz_ak_vec, ice%du_dz_3D_a)
    CALL vec2field_ak( mesh, BPA%dv_dx_ak_vec, ice%dv_dx_3D_a)
    CALL vec2field_ak( mesh, BPA%dv_dy_ak_vec, ice%dv_dy_3D_a)
    CALL vec2field_ak( mesh, BPA%dv_dz_ak_vec, ice%dv_dz_3D_a)
    ! In the BPA, gradients of w are neglected
    ice%dw_dx_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dy_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    ice%dw_dz_3D_a( mesh%vi1:mesh%vi2,:) = 0._dp
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE set_ice_velocities_to_BPA_results

! == Calculate velocities on the c-grid for solving the ice thickness equation

  SUBROUTINE map_velocities_from_b_to_c_2D( mesh, u_b, v_b, u_c, v_c)
    ! Calculate velocities on the c-grid for solving the ice thickness equation
    !
    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: u_b
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: v_b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: u_c
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: v_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_from_b_to_c_2D'
    INTEGER                                            :: aci, til, tir

    ! Add routine to path
    CALL init_routine( routine_name)

    DO aci = mesh%ci1, mesh%ci2

      til = mesh%Aci( aci,5)
      tir = mesh%Aci( aci,6)

      IF     (til == 0 .AND. tir > 0) THEN
        u_c( aci) = u_b( tir)
        v_c( aci) = v_b( tir)
      ELSEIF (tir == 0 .AND. til > 0) THEN
        u_c( aci) = u_b( til)
        v_c( aci) = v_b( til)
      ELSEIF (til >  0 .AND. tir > 0) THEN
        u_c( aci) = (u_b( til) + u_b( tir)) / 2._dp
        v_c( aci) = (v_b( til) + v_b( tir)) / 2._dp
      ELSE
        CALL crash('something is seriously wrong with the Aci array of this mesh!')
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_from_b_to_c_2D

  SUBROUTINE map_velocities_from_b_to_c_3D( mesh, u_b, v_b, u_c, v_c)
    ! Calculate velocities on the c-grid for solving the ice thickness equation
    !
    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: u_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: v_b
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: u_c
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: v_c

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'map_velocities_from_b_to_c_2D'
    INTEGER                                            :: aci, til, tir

    ! Add routine to path
    CALL init_routine( routine_name)

    DO aci = mesh%ci1, mesh%ci2

      til = mesh%Aci( aci,5)
      tir = mesh%Aci( aci,6)

      IF     (til == 0 .AND. tir > 0) THEN
        u_c( aci,:) = u_b( tir,:)
        v_c( aci,:) = v_b( tir,:)
      ELSEIF (tir == 0 .AND. til > 0) THEN
        u_c( aci,:) = u_b( til,:)
        v_c( aci,:) = v_b( til,:)
      ELSEIF (til >  0 .AND. tir > 0) THEN
        u_c( aci,:) = (u_b( til,:) + u_b( tir,:)) / 2._dp
        v_c( aci,:) = (v_b( til,:) + v_b( tir,:)) / 2._dp
      ELSE
        CALL crash('something is seriously wrong with the Aci array of this mesh!')
      END IF

    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE map_velocities_from_b_to_c_3D

END MODULE ice_velocity_main_module
