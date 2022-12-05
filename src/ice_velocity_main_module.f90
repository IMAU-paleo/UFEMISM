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
                                                 field2vec_ak, field2vec_bk, vec2field_ak, vec2field_bk, ddx_b_to_a_3D, ddy_b_to_a_3D
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
    !
    ! Derivation:
    !
    ! Conservation of mass, combined with the incompressibility
    ! condition (i.e. constant density) of ice, is described by:
    !
    !   du/dx + dv/dy + dw/dz = 0
    !
    ! Applying the zeta coordinate transformation yields:
    !
    !   du/dxp + dzeta/dx du/dzeta + dv/dxp + dzeta/dy dv/dzeta + dzeta/dz dw/dzeta = 0
    !
    ! The terms du/dxp + dv/dyp describe the two-dimensional divergence in scaled coordinates:
    !
    !   grad uv = du/dxp + dv/dyp
    !
    ! The average value over a single grid cell (Voronoi cell) of this divergence is:
    !
    !   grad uv = intint_Voronoi (grad uv) dA / intint dA = 1/A intint_Voronoi (grad uv) dA
    !
    ! By applying the divergence theorem, the surface integral over the Voronoi cell
    ! can be transformed into a loop integral over the boundary of that Voronoi cell:
    !
    !   grad uv = 1/A cint (uv * n_hat) dS
    !
    ! Here, n_hat is the outward unit normal to the Voronoi cell boundary. Substituting
    ! this into the equation for conservation of mass yields:
    !
    !   dw/dzeta = =1 / dzeta/dz [ 1/A cint (uv * n_hat) dS + dzeta/dx du/zeta + dzeta/dy dv/dzeta]
    !
    ! The vertical velocity w at the ice base is equal to the horizontal motion along
    ! the sloping ice base, plus the vertical motion of the ice base itself, plus the
    ! vertical motion of an ice particle with respect to the ice base (i.e. the basal melt rate):
    !
    !   w( z=b) = u( z=b) * dH_base/dx + v( z=b) * dH_base/dy + dH_base/dt + M_base
    !
    ! With this boundary condition, dw/dzeta can be integrated over zeta to yield w( z).

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_BMB_model),                INTENT(IN)    :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'calc_vertical_velocities'
    INTEGER                                            :: vi,k,ks,ci,vj,aci
    REAL(dp), DIMENSION(:    ), POINTER                :: H_base_a
    INTEGER                                            :: wH_base_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dx_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dy_a
    REAL(dp), DIMENSION(:    ), POINTER                :: dH_base_dt_a
    INTEGER                                            :: wdH_base_dx_a, wdH_base_dy_a, wdH_base_dt_a
    REAL(dp)                                           :: dzeta
    REAL(dp), DIMENSION(:,:  ), POINTER                :: u_3D_c
    REAL(dp), DIMENSION(:,:  ), POINTER                :: v_3D_c
    INTEGER                                            :: wu_3D_c, wv_3D_c
    REAL(dp)                                           :: cint_un_dS, dS, u_ks, v_ks, un_dS, grad_uv_ks
    REAL(dp), DIMENSION(2)                             :: n_hat
    REAL(dp)                                           :: du_dzeta_ks, dv_dzeta_ks
    REAL(dp)                                           :: dzeta_dx_ks, dzeta_dy_ks, dzeta_dz_ks
    REAL(dp)                                           :: dw_dzeta_ks

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV  ,         H_base_a       , wH_base_a       )
    CALL allocate_shared_dp_1D( mesh%nV  ,         dH_base_dx_a   , wdH_base_dx_a   )
    CALL allocate_shared_dp_1D( mesh%nV  ,         dH_base_dy_a   , wdH_base_dy_a   )
    CALL allocate_shared_dp_1D( mesh%nV  ,         dH_base_dt_a   , wdH_base_dt_a   )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , u_3D_c         , wu_3D_c         )
    CALL allocate_shared_dp_2D( mesh%nAc , C%nz  , v_3D_c         , wv_3D_c         )

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

    ! Calculate u,v on the c-grid (edges)
    CALL map_velocities_from_b_to_c_3D( mesh, ice%u_3D_b, ice%v_3D_b, u_3D_c, v_3D_c)

    ! Calculate vertical velocities by solving conservation of mass in each 3-D cell
    DO vi = mesh%vi1, mesh%vi2

      ! No ice means no velocity
      IF (ice%mask_ice_a( vi) == 0) THEN
        ice%w_3D_a( vi,:) = 0._dp
        CYCLE
      END IF

      ! Calculate the vertical velocity at the ice base
      !
      ! NOTE: BMB is defined so that a positive number means accumulation of ice;
      !       at the ice base, that means that a positive BMB means a positive
      !       value of w

      ice%w_3D_a( vi,C%nz) = (ice%u_3D_a( vi,C%nz) * dH_base_dx_a( vi)) + &
                             (ice%v_3D_a( vi,C%nz) * dH_base_dy_a( vi)) + &
                              dH_base_dt_a( vi) + BMB%BMB( vi)

      ! Exception for very thin ice / ice margin: assume horizontal stretching
      ! is negligible, so that w( z) = w( z = b)
      IF (ice%Hi_a( vi) < 10._dp) THEN
        ice%w_3D_a( vi,:) = ice%w_3D_a( vi,C%nz)
        CYCLE
      END IF ! IF (ice%mask_margin_a( vi) == 1 .OR. ice%Hi_a( vi) < 10._dp) THEN

      ! Calculate vertical velocities by integrating dw/dz over the vertical column

      DO ks = C%nz-1, 1, -1

        dzeta = C%zeta( ks+1) - C%zeta( ks)

        ! Integrate u*n_hat around the Voronoi cell boundary
        cint_un_dS = 0._dp
        DO ci = 1, mesh%nC( vi)
          vj  = mesh%C(    vi,ci)
          aci = mesh%iAci( vi,ci)
          ! Velocities at this section of the boundary
          u_ks = 0.5_dp * (u_3D_c( aci,ks) + u_3D_c( aci,ks+1))
          v_ks = 0.5_dp * (v_3D_c( aci,ks) + v_3D_c( aci,ks+1))
          ! Length of this section of the boundary
          dS = mesh%Cw( vi,ci)
          ! Outward normal vector to this section of the boundary
          n_hat = mesh%V( vj,:) - mesh%V( vi,:)
          n_hat = n_hat / NORM2( n_hat)
          ! Line integral over this section of the boundary
          un_dS = (u_ks * n_hat( 1) + v_ks * n_hat( 2)) * dS
          ! Add to loop integral
          cint_un_dS = cint_un_dS + un_dS
        END DO

        ! Calculate grad uv from the divergence theorem
        grad_uv_ks = cint_un_dS / mesh%A( vi)

        ! Calculate du/dzeta, dv/dzeta
        du_dzeta_ks = (ice%u_3D_a( vi,ks+1) - ice%u_3D_a( vi,ks)) / dzeta
        dv_dzeta_ks = (ice%v_3D_a( vi,ks+1) - ice%v_3D_a( vi,ks)) / dzeta

        ! Calculate dzeta/dx, dzeta/dy, dzeta/dz
        dzeta_dx_ks = 0.5_dp * (ice%dzeta_dx_ak( vi,ks) + ice%dzeta_dx_ak( vi,ks+1))
        dzeta_dy_ks = 0.5_dp * (ice%dzeta_dy_ak( vi,ks) + ice%dzeta_dy_ak( vi,ks+1))
        dzeta_dz_ks = 0.5_dp * (ice%dzeta_dz_ak( vi,ks) + ice%dzeta_dz_ak( vi,ks+1))

        ! Calculate dw/dzeta
        dw_dzeta_ks = -1._dp / dzeta_dz_ks * (grad_uv_ks + dzeta_dx_ks * du_dzeta_ks + dzeta_dy_ks * dv_dzeta_ks)

        ! Calculate w
        ice%w_3D_a( vi,ks) = ice%w_3D_a( vi,ks+1) - dzeta * dw_dzeta_ks

      END DO ! DO k = C%nz-1, 1, -1

    END DO ! DO vi = mesh%vi1, mesh%vi2
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wH_base_a    )
    CALL deallocate_shared( wdH_base_dx_a)
    CALL deallocate_shared( wdH_base_dy_a)
    CALL deallocate_shared( wdH_base_dt_a)
    CALL deallocate_shared( wu_3D_c      )
    CALL deallocate_shared( wv_3D_c      )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE calc_vertical_velocities

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
