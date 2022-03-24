! Created by Victor Azizi
! Module to expose all mpi functionality to UFEMISM
! created to introduced distributed memory paradigm
! A lot of thing are copied from the parallel_module.f90 file which uses shared memory paradigm
module mpi_module
  use mpi_f08
  use parallel_module, only: par ! port it here finally
  
  implicit none

  integer :: ierr ! Error flag for MPI routine

!  type parallel_info
!    integer  :: i        ! ID of this process (0 = master, >0 = slave)
!    integer  :: n        ! Total number of processes (0 = single-core, >0 = master+slaves)
!    logical  :: master   ! Whether or not the current process is the master process
!  end type

!  type(parallel_info), save :: par

contains
  ! Initialise the MPI parallelisation
  SUBROUTINE initialise_parallelisation
    implicit none
    
    ! MPI Initialisation
    ! ==================
    
    ! Use MPI to create copies of the program on all the processors, so the model can run in parallel.
    CALL MPI_INIT(ierr)
    
    ! Get rank of current process and total number of processes
    CALL MPI_COMM_RANK(       MPI_COMM_WORLD, par%i, ierr)
    CALL MPI_COMM_SIZE(       MPI_COMM_WORLD, par%n, ierr)
    par%master = (par%i == 0)
    
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)

    ! Memory use tracker !TODO REMOVE!
    IF (par%master) THEN
      par%mem%n = 0
      ALLOCATE( par%mem%h( 500000))
    END IF
  END SUBROUTINE initialise_parallelisation

  ! shortcut for global barrier
  SUBROUTINE SYNC
    implicit none
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)    
  END SUBROUTINE SYNC


end module
