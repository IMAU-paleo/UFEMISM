MODULE zeta_module

  ! Contains all the routines for setting up the vertical zeta grid (regular and staggered),
  ! with different options for the grid spacing.
  !
  ! NOTE: the zeta discretisation stuff is calculated in the mesh_operators module!

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list

  ! Import specific functionality

  IMPLICIT NONE

CONTAINS

  SUBROUTINE setup_vertical_grid
    ! Set up the vertical zeta grid

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_vertical_grid'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%nz < 3) THEN
      CALL crash('need at least 3 vertical layers!')
    ELSEIF (C%nz > 1000) THEN
      CALL warning('running the model with nz = {int_01}; are you quite sure about this?', int_01 = C%nz)
    END IF

    ! Allocate memory
    ALLOCATE( C%zeta(      C%nz  ))
    ALLOCATE( C%zeta_stag( C%nz-1))

    ! Calculate zeta
    IF     (C%choice_zeta_grid == 'regular') THEN
      CALL setup_vertical_grid_regular
    ELSEIF (C%choice_zeta_grid == 'irregular_log') THEN
      CALL setup_vertical_grid_irregular_log
    ELSEIF (C%choice_zeta_grid == 'old_15_layer_zeta') THEN
      CALL setup_vertical_grid_old_15_layer_zeta
    ELSE
      CALL crash('unknown choice_zeta_grid "' // TRIM( C%choice_zeta_grid) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_vertical_grid

  SUBROUTINE setup_vertical_grid_regular
    ! Set up the vertical zeta grid
    !
    ! Just use a simple regular grid

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_vertical_grid_regular'
    INTEGER                                            :: k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Fill zeta values
    DO k = 1, C%nz
      C%zeta( k) = REAL( k-1,dp) / REAL( C%nz-1,dp)
    END DO

    ! Calculate zeta_stag
    C%zeta_stag = (C%zeta( 1:C%nz-1) + C%zeta( 2:C%nz)) / 2._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_vertical_grid_regular

  SUBROUTINE setup_vertical_grid_irregular_log
    ! Set up the vertical zeta grid
    !
    ! This scheme ensures that the ratio between subsequent grid spacings is
    ! constant, and that the ratio between the first (surface) and last (basal)
    ! layer thickness is (approximately) equal to R

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_vertical_grid_irregular_log'
    INTEGER                                            :: k, ks
    REAL(dp)                                           :: sigma, sigma_stag

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (C%zeta_irregular_log_R <= 0._dp) CALL crash('zeta_irregular_log_R should be positive!')

    ! Exception: R = 1 implies a regular grid, but the equation becomes 0/0
    IF (C%zeta_irregular_log_R == 1._dp) THEN
      CALL setup_vertical_grid_regular
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    DO k = 1, C%nz
      ! Regular grid
      sigma = REAL( k-1,dp) / REAL( C%nz-1,dp)
      C%zeta( C%nz + 1 - k) = 1._dp - (C%zeta_irregular_log_R**sigma - 1._dp) / (C%zeta_irregular_log_R - 1._dp)
      ! Staggered grid
      IF (k < C%nz) THEN
        sigma_stag = sigma + 0.5_dp / REAL( C%nz-1,dp)
        C%zeta_stag( C%nz - k) = 1._dp - (C%zeta_irregular_log_R**sigma_stag - 1._dp) / (C%zeta_irregular_log_R - 1._dp)
      END IF
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_vertical_grid_irregular_log

  SUBROUTINE setup_vertical_grid_old_15_layer_zeta
    ! Set up the vertical zeta grid
    !
    ! Use the old irregularly-spaced 15 layers from ANICE

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_vertical_grid_old_15_layer_zeta'
    INTEGER                                            :: k
    REAL(dp)                                           :: dzeta_prev, dzeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (C%nz /= 15) CALL crash('only works when nz = 15!')

    C%zeta = (/ 0.00_dp, 0.10_dp, 0.20_dp, 0.30_dp, 0.40_dp, 0.50_dp, 0.60_dp, 0.70_dp, 0.80_dp, 0.90_dp, 0.925_dp, 0.95_dp, 0.975_dp, 0.99_dp, 1.00_dp /)

    ! Calculate zeta_stag
    C%zeta_stag = (C%zeta( 1:C%nz-1) + C%zeta( 2:C%nz)) / 2._dp

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_vertical_grid_old_15_layer_zeta

END MODULE zeta_module