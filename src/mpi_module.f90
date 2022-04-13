! Created by Victor Azizi
! Module to expose all mpi functionality to UFEMISM
! created to introduced distributed memory paradigm
! A lot of thing are copied from the parallel_module.f90 file which uses shared memory paradigm
module mpi_module
  use mpi_f08
  use configuration_module, only: dp
  use parallel_module, only: par, partition_list ! port it here finally
  
  implicit none

  integer :: ierr ! Error flag for MPI routine

!  type parallel_info
!    integer  :: i        ! ID of this process (0 = master, >0 = slave)
!    integer  :: n        ! Total number of processes (0 = single-core, >0 = master+slaves)
!    logical  :: master   ! Whether or not the current process is the master process
!  end type

!  type(parallel_info), save :: par

interface allgather_array
  procedure :: allgather_array_dp
  procedure :: allgather_array_int
end interface

contains
  subroutine allgather_array_dp(array, i1_, i2_)
    implicit none
    real(dp),              dimension(:), intent(inout) :: array
    integer, optional,                   intent(in)    :: i1_,i2_

    integer                                            :: i1,i2, err, n
    integer, dimension(1:par%n)                        :: counts, displs

    ! Gather sizes that will be sent
    if (present(i1_) .and. present(i2_)) then
      i1 = i1_
      i2 = i2_
    else
      call partition_list(size(array), par%i, par%n, i1, i2)
    endif
    call mpi_allgather( i2-i1+1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, err)

    ! Calculate offsets through the sizes
    displs(1) = 0
    do n=2,size(displs)
      displs(n) = displs(n-1) + counts(n-1)
    end do
      
    ! Send everything to master
    call mpi_allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL &
                       , array, counts, displs, MPI_REAL8, MPI_COMM_WORLD, err)
  end subroutine allgather_array_dp
  subroutine allgather_array_int(array,i1_,i2_)
    implicit none
    integer,               dimension(:), intent(inout) :: array
    integer, optional,                   intent(in)    :: i1_,i2_

    integer                                            :: i1,i2, err, n
    integer, dimension(1:par%n)                        :: counts, displs

    ! Gather sizes that will be sent
    if (present(i1_) .and. present(i2_)) then
      i1 = i1_
      i2 = i2_
    else
      call partition_list(size(array), par%i, par%n, i1, i2)
    endif
    call mpi_allgather( i2-i1+1, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, err)

    ! Calculate offsets through the sizes
    displs(1) = 0
    do n=2,size(displs)
      displs(n) = displs(n-1) + counts(n-1)
    end do
      
    ! Send everything to master
    call mpi_allgatherv( MPI_IN_PLACE, 0, MPI_DATATYPE_NULL &
                       , array, counts, displs, MPI_INTEGER, MPI_COMM_WORLD, err)
  end subroutine allgather_array_int


end module
