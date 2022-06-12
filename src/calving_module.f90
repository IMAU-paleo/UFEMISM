MODULE calving_module

  ! Contains all the routines for calving.

  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parallel_module,                 ONLY: par, sync
  USE data_types_module,               ONLY: type_mesh, type_ice_model

  IMPLICIT NONE

CONTAINS

! == The main routines that should be called from the main ice model/program
! ==========================================================================

  SUBROUTINE apply_calving_law( mesh, ice, calving_event)
    ! Apply the selected calving law

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    LOGICAL, DIMENSION(par%n),           INTENT(INOUT) :: calving_event

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_calving_law'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Set the default state: if no calving, exit loop after this round
    calving_event(par%i+1) = .FALSE.

    ! Apply the selected calving law
    IF     (C%choice_calving_law == 'none') THEN
      ! No calving at all; just return and exit loop
    ELSEIF (C%choice_calving_law == 'threshold_thickness') THEN
      CALL threshold_thickness_calving( mesh, ice, calving_event)
    ELSE
      CALL crash('unknown choice_calving_law"' // TRIM(C%choice_calving_law) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_calving_law

! == Routines for different calving laws
! ===================================

  SUBROUTINE threshold_thickness_calving( mesh, ice, calving_event)
    ! A simple threshold thickness calving law

    IMPLICIT NONE

    ! Input variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    LOGICAL, DIMENSION(par%n),           INTENT(INOUT) :: calving_event

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'threshold_thickness_calving'
    INTEGER                                            :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Apply calving law
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_cf_a( vi) == 1 .AND. ice%Hi_eff_cf_a( vi) < C%calving_threshold_thickness) THEN
        ice%Hi_a( vi) = 0._dp
        calving_event(par%i+1) = .TRUE.
      END IF
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE threshold_thickness_calving

  SUBROUTINE remove_unconnected_shelves( mesh, ice)
    ! Use a flood-fill algorithm to find all shelves connected to sheets.
    ! Remove all other shelves.

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remove_unconnected_shelves'
    INTEGER                                            :: vi, ci, vc
    INTEGER, DIMENSION(:    ), ALLOCATABLE             :: map
    INTEGER, DIMENSION(:    ), ALLOCATABLE             :: stack
    INTEGER                                            :: stackN

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      ! Allocate memory for map and stack
      ALLOCATE( map(   mesh%nV ))
      ALLOCATE( stack( mesh%nV ))

      map    = 0
      stack  = 0
      stackN = 0

      ! Fill the stack with all shelf-next-to-sheet points
      DO vi = 1, mesh%nV
        IF (ice%mask_shelf_a( vi) == 1) THEN
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            IF (ice%mask_sheet_a( vc) == 1) THEN
              stackN = stackN + 1                  ! Increase stack total
              stack( stackN) = vi                  ! Add vertex to the stack
              map( vi) = 1                         ! Mark it on the map
              EXIT                                 ! No need to check further
            END IF
          END DO
        END IF
      END DO

      ! Mark all connected shelf points on the map using a flood fill
      DO WHILE (stackN > 0)

        ! Remove the last element from the stack
        vi = stack( stackN)
        stackN = stackN - 1

        ! Check if this element is shelf. If so, mark it on the map,
        ! and add its neighbours to the stack (if they're not in it yet).
        ! If not, do nothing
        IF (ice%mask_shelf_a( vi) == 1) THEN
          ! This element is shelf. Mark it on the map
          map( vi) = 2
          ! Add its neighbours to the stack
          DO ci = 1, mesh%nC( vi)
            vc = mesh%C( vi,ci)
            IF (map( vc)==0) THEN
              map( vc) = 1
              stackN = stackN + 1
              stack( stackN) = vc
            END IF
          END DO
        END IF

      END DO ! DO WHILE (stackN > 0)

      ! Remove ice for all unconnected shelves
      DO vi = 1, mesh%nV
        IF (ice%mask_shelf_a( vi) == 1 .AND. map( vi) == 0) THEN
          print*, 'bang!'
          ice%Hi_a( vi) = 0._dp
        END IF
      END DO

      ! Clean up after yourself
      DEALLOCATE( map)
      DEALLOCATE( stack)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remove_unconnected_shelves

END MODULE calving_module
