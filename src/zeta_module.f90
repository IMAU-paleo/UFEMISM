MODULE zeta_module

  ! Contains all the routines for setting up the vertical zeta grid and its discretisation

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning

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

    ! Calculate zeta_stag
    C%zeta_stag = (C%zeta( 1:C%nz-1) + C%zeta( 2:C%nz)) / 2._dp

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

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_vertical_grid_regular

  SUBROUTINE setup_vertical_grid_irregular_log
    ! Set up the vertical zeta grid
    !
    ! Let the distance dzeta between each consecutive pair of layers decrease by R (=0.95 by default)

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'setup_vertical_grid_irregular_log'
    INTEGER                                            :: k
    REAL(dp)                                           :: dzeta_prev, dzeta

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Fill zeta values
    C%zeta( 1) = 0._dp
    C%zeta( 2) = 1._dp

    DO k = 3, C%nz
      dzeta_prev = C%zeta( k-1) - C%zeta( k-2)
      dzeta = dzeta_prev * C%zeta_irregular_log_R
      C%zeta( k) = C%zeta( k-1) + dzeta
    END DO

    C%zeta = C%zeta / C%zeta( C%nz)

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
    IF (C%nz /= 15) CALL crash('onyl works when nz = 15!')

    C%zeta = (/ 0.00_dp, 0.10_dp, 0.20_dp, 0.30_dp, 0.40_dp, 0.50_dp, 0.60_dp, 0.70_dp, 0.80_dp, 0.90_dp, 0.925_dp, 0.95_dp, 0.975_dp, 0.99_dp, 1.00_dp /)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE setup_vertical_grid_old_15_layer_zeta

END MODULE zeta_module