module calving_module

  ! Contains all the routines for calving.

  use mpi
  use configuration_module,          only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parallel_module,               only : par, sync
  use data_types_module,             only : type_mesh, type_ice_model, type_reference_geometry
  use utilities_module,              only : is_floating
  use general_ice_model_data_module, only : determine_masks

  implicit none

contains

! ===== The main routines =====
! =============================

  ! Apply calving law in chain mode: hunt until extintion
  ! SUBROUTINE run_calving_model( mesh, ice)

  !   IMPLICIT NONE

  !   ! In- and output variables:
  !   TYPE(type_mesh),                     INTENT(IN)    :: mesh
  !   TYPE(type_ice_model),                INTENT(INOUT) :: ice

  !   ! Local variables:
  !   CHARACTER(LEN=256),    PARAMETER                   :: routine_name = 'run_calving_model'
  !   LOGICAL, DIMENSION(:), POINTER                     :: calving_event
  !   INTEGER                                            :: wcalving_event, calving_round

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! Allocate the calving event flag
  !   CALL allocate_shared_bool_1D(par%n, calving_event, wcalving_event)

  !   calving_round = 0                                       ! Initialise loop counter
  !   calving_event = .TRUE.                                  ! Let all processes enter the loop
  !   DO WHILE ( ANY(calving_event) .AND. &                   ! Exit loop if no more calving occurs
  !              calving_round < C%max_calving_rounds)        ! Exit loop if it exceeds the max rounds allowed
  !     CALL determine_masks_ice( mesh, ice)                  ! Update the ice sheet and ice shelf masks
  !     CALL determine_masks_transitions( mesh, ice)          ! Update the calving front mask
  !     CALL determine_floating_margin_fraction( mesh, ice)   ! Update the fractions for new calving fronts
  !     CALL apply_calving_law( mesh, ice, calving_event)     ! Apply calving law and update calving event flag
  !     calving_round = calving_round + 1                     ! Increase the counter
  !     CALL sync                                             ! Let all processes reach this point before next loop
  !   END DO
  !   ! IF (par%master .AND. calving_round == C%max_calving_rounds .AND. C%max_calving_rounds > 5) THEN
  !   !   CALL warning('max_calving_rounds reached! Thin ice potentially still floating around...')
  !   ! END IF

  !   ! Clean up after yourself
  !   CALL deallocate_shared(wcalving_event)

  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE run_calving_model

  ! ! Apply the selected calving law
  ! SUBROUTINE apply_calving_law( mesh, ice, calving_event)

  !   IMPLICIT NONE

  !   ! Input variables:
  !   TYPE(type_mesh),                     INTENT(IN)    :: mesh
  !   TYPE(type_ice_model),                INTENT(INOUT) :: ice
  !   LOGICAL, DIMENSION(par%n),           INTENT(INOUT) :: calving_event

  !   ! Local variables:
  !   CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'apply_calving_law'

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! Reset the flag to false and check if any more calving occurs
  !   calving_event(par%i+1) = .FALSE.

  !   ! Apply the selected calving law
  !   IF     (C%choice_calving_law == 'none') THEN
  !     ! No calving at all; just return and exit loop
  !   ELSEIF (C%choice_calving_law == 'threshold_thickness') THEN
  !     CALL threshold_thickness_calving( mesh, ice, calving_event)
  !   ELSE
  !     CALL crash('unknown choice_calving_law"' // TRIM(C%choice_calving_law) // '"!')
  !   END IF

  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE apply_calving_law

! ===== Routines for different calving laws =====
! ===============================================

  ! ! A simple threshold thickness calving law
  ! SUBROUTINE threshold_thickness_calving( mesh, ice, calving_event)

  !   IMPLICIT NONE

  !   ! Input variables:
  !   TYPE(type_mesh),                     INTENT(IN)    :: mesh
  !   TYPE(type_ice_model),                INTENT(INOUT) :: ice
  !   LOGICAL, DIMENSION(par%n),           INTENT(INOUT) :: calving_event

  !   ! Local variables:
  !   CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'threshold_thickness_calving'
  !   INTEGER                                            :: vi

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! Apply calving law
  !   DO vi = mesh%vi1, mesh%vi2
  !     ! If vertex is calving front
  !     IF (ice%mask_cf_a( vi) == 1) THEN
  !       ! If its effective (shelf) or modelled (sheet) thickness is below the threshold
  !       IF ( (ice%mask_shelf_a( vi) == 1 .AND. ice%Hi_eff_cf_a( vi) < C%calving_threshold_thickness_shelf) .OR. &
  !            (ice%mask_sheet_a( vi) == 1 .AND. ice%Hi_a( vi)        < C%calving_threshold_thickness_sheet) ) THEN
  !         ! Remove ice from this vertex
  !         ice%Hi_a( vi) = 0._dp
  !         ! Calving event occurred. This will cause the calving loop to
  !         ! do a whole another check over the entire (new) calving front
  !         ! after this iteration.
  !         calving_event(par%i+1) = .TRUE.
  !       END IF
  !     END IF
  !   END DO
  !   CALL sync

  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE threshold_thickness_calving

! ===== Extra tools for ice shelf removal =====
! =============================================

  ! Force the removal of specific ice shelf areas
  subroutine imposed_shelf_removal(mesh, ice, refgeo_PD, refgeo_GIAeq)

    implicit none

    ! In- and output variables:
    type(type_mesh),               intent(in)    :: mesh
    type(type_ice_model),          intent(inout) :: ice
    type(type_reference_geometry), intent(in)    :: refgeo_PD
    type(type_reference_geometry), intent(in)    :: refgeo_GIAeq

    ! Local variables:
    character(len=256), parameter                :: routine_name = 'imposed_shelf_removal'
    integer                                      :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! If so specified, remove all floating ice
    if (C%do_remove_shelves) then
      do vi = mesh%vi1, mesh%vi2
        if (is_floating( ice%Hi_a( vi), ice%Hb_a( vi), ice%SL_a( vi))) then
          ice%Hi_a( vi) = 0._dp
        end if
      end do
    end if

    ! If so specified, remove all floating ice beyond the present-day calving front
    if (C%remove_shelves_larger_than_PD) then
      do vi = mesh%vi1, mesh%vi2
        if (refgeo_PD%Hi( vi) == 0._dp .and. refgeo_PD%Hb( vi) < 0._dp) then
          ice%Hi_a( vi) = 0._dp
        end if
      end do
    end if

    ! If so specified, remove all floating ice crossing the continental shelf edge
    if (C%continental_shelf_calving) then
      do vi = mesh%vi1, mesh%vi2
        if (refgeo_GIAeq%Hi( vi) == 0._dp .and. refgeo_GIAeq%Hb( vi) < C%continental_shelf_min_height) then
          ice%Hi_a( vi) = 0._dp
        end if
      end do
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine imposed_shelf_removal

  ! Use a flood-fill algorithm to remove all shelves not connected to sheets
  subroutine remove_unconnected_shelves( mesh, ice)

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=256), parameter       :: routine_name = 'remove_unconnected_shelves'
    integer                             :: vi, ci, vc
    integer, dimension(:), allocatable  :: map
    integer, dimension(:), allocatable  :: stack
    integer                             :: stackN

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine current calving front
    call determine_masks( mesh, ice)

    ! Allocate memory for map and stack
    allocate( map(   mesh%nV ))
    allocate( stack( mesh%nV ))

    ! Initialise stack
    map    = 0
    stack  = 0
    stackN = 0

    ! Fill the stack with all shelf-next-to-sheet points
    do vi = 1, mesh%nV
      if (ice%mask_shelf_a( vi) == 1) then       ! If this element is shelf
        do ci = 1, mesh%nC( vi)                  ! Loop over all neighbours
          vc = mesh%C( vi,ci)                    ! Get index of neighbour
          if (ice%mask_sheet_a( vc) == 1) then   ! If neighbour is sheet
            stackN = stackN + 1                  ! Increase stack total
            stack( stackN) = vi                  ! Add vertex to the stack
            map( vi) = 1                         ! Mark it on the map
            exit                                 ! No need to check further
          end if
        end do
      end if
    end do

    ! Mark all connected shelf points on the map using a flood fill
    do while (stackN > 0)
      vi = stack( stackN)                        ! Get index of last element
      stackN = stackN - 1                        ! Remove element from stack
      if (ice%mask_shelf_a( vi) == 1) then       ! If this element is shelf
        map( vi) = 2                             ! Mark it on the map
        do ci = 1, mesh%nC( vi)                  ! Loop over all neighbours
          vc = mesh%C( vi,ci)                    ! Get index of neighbour
          if (map( vc)==0) then                  ! If not already in the map
            map( vc) = 1                         ! Add this neighbour to the map
            stackN = stackN + 1                  ! Increase the stack counter
            stack( stackN) = vc                  ! Add this neighbour to the stack
          end if
        end do
      end if                                     ! Else, do nothing (removal from stack is enough)
    end do

    ! Remove ice for all unconnected shelves
    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_shelf_a( vi) == 1 .and. map( vi) == 0) then
        ice%Hi_a( vi) = 0._dp
      end if
    end do

    ! Clean up after yourself
    deallocate( map)
    deallocate( stack)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_unconnected_shelves

end module calving_module
