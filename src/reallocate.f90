module reallocate_mod
  use configuration_module, only: dp

implicit none
interface reallocate
  procedure :: reallocate_dp_1d
  procedure :: reallocate_dp_2d
  procedure :: reallocate_int_1d
  procedure :: reallocate_int_2d
end interface

contains
subroutine reallocate_dp_1d(array,newx)
  implicit none
  real(dp), allocatable, dimension(:), intent(inout)   :: array
  integer                            , intent(in)      :: newx
  real(dp), allocatable, dimension(:)                  :: newarray

  ! allocate, move, swap pointer(bonus implicit deallocate)
  allocate(newarray(newx))
  newarray(1:min(newx,size(array,1))) = array
  call move_alloc(newarray, array)
end subroutine
subroutine reallocate_dp_2d(array,newx, newy)
  implicit none
  real(dp), allocatable, dimension(:,:), intent(inout) :: array
  integer                              , intent(in)    :: newx
  integer                              , intent(in)    :: newy
  real(dp), allocatable, dimension(:,:)                :: newarray
  
  ! allocate, move, swap pointer(bonus implicit deallocate)
  allocate(newarray(newx,newy))
  newarray(1:min(newx,size(array,1)),1:min(newy,size(array,2))) = array
  call move_alloc(newarray, array)
end subroutine
subroutine reallocate_int_1d(array,newx)
  implicit none
  integer,  allocatable, dimension(:), intent(inout)   :: array
  integer                            , intent(in)      :: newx
  integer,  allocatable, dimension(:)                  :: newarray

  ! allocate, move, swap pointer(bonus implicit deallocate)
  allocate(newarray(newx))
  newarray(1:min(newx,size(array,1))) = array
  call move_alloc(newarray, array)
end subroutine
subroutine reallocate_int_2d(array,newx, newy)
  implicit none
  integer,  allocatable, dimension(:,:), intent(inout) :: array
  integer                              , intent(in)    :: newx
  integer                              , intent(in)    :: newy
  integer,  allocatable, dimension(:,:)                :: newarray
  
  ! allocate, move, swap pointer(bonus implicit deallocate)
  allocate(newarray(newx,newy))
  newarray(1:min(newx,size(array,1)),1:min(newy,size(array,2))) = array
  call move_alloc(newarray, array)
end subroutine
end module
