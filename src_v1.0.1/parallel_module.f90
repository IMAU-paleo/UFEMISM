MODULE parallel_module
  ! A collection of different routines that make parallel programming in UFEMISM a lot easier.

  USE mpi
  USE configuration_module,        ONLY: dp
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_LOC
  
  IMPLICIT NONE
    
  TYPE parallel_info
    
    INTEGER :: i        ! ID of this process (0 = master, >0 = slave)
    INTEGER :: n        ! Total number of processes (0 = single-core, >0 = master+slaves)
    LOGICAL :: master   ! Whether or not the current process is the master process
    
  END TYPE parallel_info
    
  TYPE(parallel_info), SAVE :: par

CONTAINS

  ! Synchronise the different processes
  SUBROUTINE sync
    ! Use MPI_BARRIER to synchronise all the processes         
    
    IMPLICIT NONE
    
    INTEGER                                           :: ierr
    
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr)
  
  END SUBROUTINE sync
    
  ! Allocate some shared memory space on all the processes (cannot be done on only one process).
  ! The slave process allocate zero space. The master process allocates actual space,
  ! and sends the window to this space to all the slaves. All processes then associate
  ! their own pointers with the master's memory space, so that a subroutine that accepts a pointer as
  ! an input argument works for any process.
  ! All processes save the "window" to their own memory space, which can later be used to deallocate that memory space.
  SUBROUTINE allocate_shared_memory_int_0D(              p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                    POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = 4_MPI_ADDRESS_KIND
      disp_unit   = 4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN  
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p)
  
  END SUBROUTINE allocate_shared_memory_int_0D  
  SUBROUTINE allocate_shared_memory_int_1D(  n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    INTEGER,  DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*4_MPI_ADDRESS_KIND
      disp_unit   = 4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN  
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
  
  END SUBROUTINE allocate_shared_memory_int_1D  
  SUBROUTINE allocate_shared_memory_int_2D(  n1, n2,     p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2     ! Dimension(s) of memory to be allocated
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*4_MPI_ADDRESS_KIND
      disp_unit   = n1*4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN 
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2])
  
  END SUBROUTINE allocate_shared_memory_int_2D  
  SUBROUTINE allocate_shared_memory_int_3D(  n1, n2, n3, p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2, n3 ! Dimension(s) of memory to be allocated
    INTEGER,  DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*n3*4_MPI_ADDRESS_KIND
      disp_unit   = n1*n2*4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN   
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2, n3])
  
  END SUBROUTINE allocate_shared_memory_int_3D  
  SUBROUTINE allocate_shared_memory_dp_0D(               p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    REAL(dp),                   POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = 8_MPI_ADDRESS_KIND
      disp_unit   = 8
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p)
  
  END SUBROUTINE allocate_shared_memory_dp_0D  
  SUBROUTINE allocate_shared_memory_dp_1D(   n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*8_MPI_ADDRESS_KIND
      disp_unit   = 8
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN  
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
  
  END SUBROUTINE allocate_shared_memory_dp_1D  
  SUBROUTINE allocate_shared_memory_dp_2D(   n1, n2,     p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2     ! Dimension(s) of memory to be allocated
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*8_MPI_ADDRESS_KIND
      disp_unit   = n1*8
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN    
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2])
  
  END SUBROUTINE allocate_shared_memory_dp_2D  
  SUBROUTINE allocate_shared_memory_dp_3D(   n1, n2, n3, p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2, n3 ! Dimension(s) of memory to be allocated
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*n2*n3*8_MPI_ADDRESS_KIND
      disp_unit   = n1*n2*8
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN   
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2, n3])
  
  END SUBROUTINE allocate_shared_memory_dp_3D  
  SUBROUTINE allocate_shared_memory_bool_1D( n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    LOGICAL,  DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    IF (par%master) THEN
      windowsize  = n1*4_MPI_ADDRESS_KIND
      disp_unit   = 4
    ELSE
      windowsize  = 0_MPI_ADDRESS_KIND
      disp_unit   = 1
    END IF
    
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    IF (.NOT. par%master) THEN   
      ! Get the baseptr, size and disp_unit values of the master's memory space.
      CALL MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    END IF
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
  
  END SUBROUTINE allocate_shared_memory_bool_1D
    
  ! Use the "window" to the allocated memory space (which consists of zero bytes for the slave)
  ! to deallocate that memory.
  SUBROUTINE deallocate_shared_memory( win)
    ! Use MPI_WIN_FREE to deallocate shared memory space for an array.
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(INOUT) :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    
    CALL MPI_WIN_FREE( win, ierr)  
  
  END SUBROUTINE deallocate_shared_memory
  
  ! Change the amount of shared memory that is used for a data array. First, save the data
  ! in a temporary array (done by the master process), then deallocate the shared memory
  ! (including the zero-bytes memory allocated by the slaves), then call one of the
  ! allocate_shared_memory subroutines to allocate new memory space. Copy the data from
  ! the temporary array back to the new space, then deallocate the temporary array.
  SUBROUTINE adapt_shared_memory_int_1D(  nin, n1,         p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1
    INTEGER,  DIMENSION(:    ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: temp_mem
        
    ! Only the master process allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    IF (par%master) THEN
      ALLOCATE(temp_mem(nin))
      temp_mem = p(1:nin)
    ELSE
      ALLOCATE(temp_mem(1))
      temp_mem = 0
    END IF
    
    ! Deallocate old memory (through MPI window), nullify pointer, allocate new memory
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_int_1D( n1, p, win)
    
    ! Only the master moves the data from the temporary memory to the new memory, deallocates temporary memory.
    IF (par%master) THEN
      p = 0
      p(1:nin) = temp_mem
    END IF
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_int_1D
  SUBROUTINE adapt_shared_memory_int_2D(  nin, n1, n2,     p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1,n2
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE       :: temp_mem
            
    ! Only the master process allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    IF (par%master) THEN
      ALLOCATE(temp_mem(nin,n2))
      temp_mem = p(1:nin,:)
    ELSE
      ALLOCATE(temp_mem(1,1))
      temp_mem = 0
    END IF
    
    ! Deallocate old memory (through MPI window), nullify pointer, allocate new memory
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_int_2D( n1, n2, p, win)
    
    ! Only the master moves the data from the temporary memory to the new memory, deallocates temporary memory.
    IF (par%master) THEN
      p = 0
      p(1:nin,:) = temp_mem
    END IF
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_int_2D
  SUBROUTINE adapt_shared_memory_int_3D(  nin, n1, n2, n3, p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1,n2,n3
    INTEGER,  DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE       :: temp_mem
            
    ! Only the master process allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    IF (par%master) THEN
      ALLOCATE(temp_mem(nin,n2,n3))
      temp_mem = p(1:nin,:,:)
    ELSE
      ALLOCATE(temp_mem(1,1,1))
      temp_mem = 0
    END IF
    
    ! Deallocate old memory (through MPI window), nullify pointer, allocate new memory
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_int_3D( n1, n2, n3, p, win)
    
    ! Only the master moves the data from the temporary memory to the new memory, deallocates temporary memory.
    IF (par%master) THEN
      p = 0
      p(1:nin,:,:) = temp_mem
    END IF
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_int_3D
  SUBROUTINE adapt_shared_memory_dp_1D(   nin, n1,         p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: temp_mem
        
    ! Only the master process allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    IF (par%master) THEN
      ALLOCATE(temp_mem(nin))
      temp_mem = p(1:nin)
    ELSE
      ALLOCATE(temp_mem(1))
      temp_mem = 0._dp
    END IF
    
    ! Deallocate old memory (through MPI window), nullify pointer, allocate new memory
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dp_1D( n1, p, win)
    
    ! Only the master moves the data from the temporary memory to the new memory, deallocates temporary memory.
    IF (par%master) THEN
      p = 0._dp
      p(1:nin) = temp_mem
    END IF
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dp_1D
  SUBROUTINE adapt_shared_memory_dp_2D(   nin, n1, n2,     p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1,n2
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: temp_mem
            
    ! Only the master process allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    IF (par%master) THEN
      ALLOCATE(temp_mem(nin,n2))
      temp_mem = p(1:nin,:)
    ELSE
      ALLOCATE(temp_mem(1,1))
      temp_mem = 0._dp
    END IF
    
    ! Deallocate old memory (through MPI window), nullify pointer, allocate new memory
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dp_2D( n1, n2, p, win)
    
    ! Only the master moves the data from the temporary memory to the new memory, deallocates temporary memory.
    IF (par%master) THEN
      p = 0._dp
      p(1:nin,:) = temp_mem
    END IF
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dp_2D
  SUBROUTINE adapt_shared_memory_dp_3D(   nin, n1, n2, n3, p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1,n2,n3
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: temp_mem
            
    ! Only the master process allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    IF (par%master) THEN
      ALLOCATE(temp_mem(nin,n2,n3))
      temp_mem = p(1:nin,:,:)
    ELSE
      ALLOCATE(temp_mem(1,1,1))
      temp_mem = 0._dp
    END IF
    
    ! Deallocate old memory (through MPI window), nullify pointer, allocate new memory
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dp_3D( n1, n2, n3, p, win)
    
    ! Only the master moves the data from the temporary memory to the new memory, deallocates temporary memory.
    IF (par%master) THEN
      p = 0._dp
      p(1:nin,:,:) = temp_mem
    END IF
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dp_3D
  SUBROUTINE adapt_shared_memory_bool_1D( nin, n1,         p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1
    LOGICAL,  DIMENSION(:    ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE       :: temp_mem
        
    ! Only the master process allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    IF (par%master) THEN
      ALLOCATE(temp_mem(nin))
      temp_mem = p(1:nin)
    ELSE
      ALLOCATE(temp_mem(1))
      temp_mem = .FALSE.
    END IF
    
    ! Deallocate old memory (through MPI window), nullify pointer, allocate new memory
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_bool_1D( n1, p, win)
    
    ! Only the master moves the data from the temporary memory to the new memory, deallocates temporary memory.
    IF (par%master) THEN
      p = .FALSE.
      p(1:nin) = temp_mem
    END IF
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_bool_1D
  
  ! Properly "distributed" memory, right now used only in mesh generation
  ! =====================================================================
  SUBROUTINE allocate_shared_memory_dist_int_0D(              p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                    POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    windowsize  = 4_MPI_ADDRESS_KIND
    disp_unit   = 4
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
         
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, par%i, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p)
  
  END SUBROUTINE allocate_shared_memory_dist_int_0D  
  SUBROUTINE allocate_shared_memory_dist_int_1D(  n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space.    
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    INTEGER,  DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! Allocate MPI-shared memory for data array, with an associated window object    
    windowsize  = n1*4_MPI_ADDRESS_KIND
    disp_unit   = 4
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
    
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, par%i, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
    
  END SUBROUTINE allocate_shared_memory_dist_int_1D
  SUBROUTINE allocate_shared_memory_dist_int_2D(  n1, n2,     p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2     ! Dimension(s) of memory to be allocated
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    windowsize  = n1*n2*4_MPI_ADDRESS_KIND
    disp_unit   = n1*4
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
      
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, par%i, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2])
  
  END SUBROUTINE allocate_shared_memory_dist_int_2D  
  SUBROUTINE allocate_shared_memory_dist_int_3D(  n1, n2, n3, p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2, n3 ! Dimension(s) of memory to be allocated
    INTEGER,  DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    windowsize  = n1*n2*n3*4_MPI_ADDRESS_KIND
    disp_unit   = n1*n2*4
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
          
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, par%i, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2, n3])
  
  END SUBROUTINE allocate_shared_memory_dist_int_3D 
  SUBROUTINE allocate_shared_memory_dist_dp_0D(               p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    REAL(dp),                   POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    windowsize  = 8_MPI_ADDRESS_KIND
    disp_unit   = 8
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
          
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, par%i, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p)
  
  END SUBROUTINE allocate_shared_memory_dist_dp_0D  
  SUBROUTINE allocate_shared_memory_dist_dp_1D(   n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    windowsize  = n1*8_MPI_ADDRESS_KIND
    disp_unit   = 8
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
       
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, par%i, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
  
  END SUBROUTINE allocate_shared_memory_dist_dp_1D  
  SUBROUTINE allocate_shared_memory_dist_dp_2D(   n1, n2,     p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2     ! Dimension(s) of memory to be allocated
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    windowsize  = n1*n2*8_MPI_ADDRESS_KIND
    disp_unit   = n1*8
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
      
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, par%i, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2])
  
  END SUBROUTINE allocate_shared_memory_dist_dp_2D  
  SUBROUTINE allocate_shared_memory_dist_dp_3D(   n1, n2, n3, p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1, n2, n3 ! Dimension(s) of memory to be allocated
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    windowsize  = n1*n2*n3*8_MPI_ADDRESS_KIND
    disp_unit   = n1*n2*8
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
         
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, par%i, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1, n2, n3])
  
  END SUBROUTINE allocate_shared_memory_dist_dp_3D  
  SUBROUTINE allocate_shared_memory_dist_bool_1D( n1,         p, win)
    ! Use MPI_WIN_ALLOCATE_SHARED to allocate shared memory space for an array.
    ! Return a pointer associated with that memory space. Makes it so that all processes
    ! can access the same memory directly.      
    
    IMPLICIT NONE
  
    INTEGER,                             INTENT(IN)    :: n1         ! Dimension(s) of memory to be allocated
    LOGICAL,  DIMENSION(:    ), POINTER, INTENT(OUT)   :: p          ! Pointer to memory space
    INTEGER,                             INTENT(OUT)   :: win        ! MPI window to the allocated memory space
    
    INTEGER                                            :: ierr
    INTEGER(KIND=MPI_ADDRESS_KIND)                     :: windowsize
    INTEGER                                            :: disp_unit
    TYPE(C_PTR)                                        :: baseptr    
    
    ! ==========
    ! Allocate MPI-shared memory for data array, with an associated window object
    ! (Needs to be called in Master and Slaves, but only Master actually allocates any space)
    
    windowsize  = n1*4_MPI_ADDRESS_KIND
    disp_unit   = 4
    CALL MPI_WIN_ALLOCATE_SHARED(windowsize, disp_unit, MPI_INFO_NULL, MPI_COMM_WORLD, baseptr, win, ierr)
       
    ! Get the baseptr, size and disp_unit values of this memory space.
    CALL MPI_WIN_SHARED_QUERY( win, par%i, windowsize, disp_unit, baseptr, ierr)
    
    ! Associate a pointer with this memory space.
    CALL C_F_POINTER(baseptr, p, [n1])
  
  END SUBROUTINE allocate_shared_memory_dist_bool_1D
  
  SUBROUTINE adapt_shared_memory_dist_int_1D(  nin, n1,         p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1
    INTEGER,  DIMENSION(:    ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    INTEGER,  DIMENSION(:    ), ALLOCATABLE       :: temp_mem
        
    ! Allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    
    ALLOCATE(temp_mem(nin))
    temp_mem = p(1:nin)    
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dist_int_1D( n1, p, win)
    p = 0
    p(1:nin) = temp_mem
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dist_int_1D
  SUBROUTINE adapt_shared_memory_dist_int_2D(  nin, n1, n2,     p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1,n2
    INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    INTEGER,  DIMENSION(:,:  ), ALLOCATABLE       :: temp_mem
        
    ! Allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    
    ALLOCATE(temp_mem(nin,n2))
    temp_mem = p(1:nin,:)
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dist_int_2D( n1, n2, p, win)
    p = 0
    p(1:nin,:) = temp_mem
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dist_int_2D
  SUBROUTINE adapt_shared_memory_dist_int_3D(  nin, n1, n2, n3, p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1,n2,n3
    INTEGER,  DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    INTEGER,  DIMENSION(:,:,:), ALLOCATABLE       :: temp_mem
        
    ! Allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    
    ALLOCATE(temp_mem(nin,n2,n3))
    temp_mem = p(1:nin,:,:)
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dist_int_3D( n1, n2, n3, p, win)
    p = 0
    p(1:nin,:,:) = temp_mem
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dist_int_3D
  SUBROUTINE adapt_shared_memory_dist_dp_1D(   nin, n1,         p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    REAL(dp), DIMENSION(:    ), ALLOCATABLE       :: temp_mem
        
    ! Allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    
    ALLOCATE(temp_mem(nin))
    temp_mem = p(1:nin)    
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dist_dp_1D( n1, p, win)
    p = 0._dp
    p(1:nin) = temp_mem
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dist_dp_1D
  SUBROUTINE adapt_shared_memory_dist_dp_2D(   nin, n1, n2,     p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1,n2
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE       :: temp_mem
        
    ! Allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    
    ALLOCATE(temp_mem(nin,n2))
    temp_mem = p(1:nin,:)
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dist_dp_2D( n1, n2, p, win)
    p = 0._dp
    p(1:nin,:) = temp_mem
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dist_dp_2D
  SUBROUTINE adapt_shared_memory_dist_dp_3D(   nin, n1, n2, n3, p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1,n2,n3
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE       :: temp_mem
        
    ! Allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    
    ALLOCATE(temp_mem(nin,n2,n3))
    temp_mem = p(1:nin,:,:)
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dist_dp_3D( n1, n2, n3, p, win)
    p = 0._dp
    p(1:nin,:,:) = temp_mem
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dist_dp_3D
  SUBROUTINE adapt_shared_memory_dist_bool_1D( nin, n1,         p, win)    
    
    IMPLICIT NONE
 
    INTEGER,                             INTENT(IN)    :: nin,n1
    LOGICAL,  DIMENSION(:    ), POINTER, INTENT(INOUT) :: p
    INTEGER,                             INTENT(INOUT) :: win
    
    LOGICAL,  DIMENSION(:    ), ALLOCATABLE       :: temp_mem
        
    ! Allocates temporary memory where the data is stored so the old memory can be deallocated and extended.
    
    ALLOCATE(temp_mem(nin))
    temp_mem = p(1:nin)    
    CALL deallocate_shared_memory( win)
    NULLIFY( p)
    CALL allocate_shared_memory_dist_bool_1D( n1, p, win)
    p = .FALSE.
    p(1:nin) = temp_mem
    DEALLOCATE(temp_mem)
   
  END SUBROUTINE adapt_shared_memory_dist_bool_1D
  
  SUBROUTINE ShareMemoryAccess_int_0D( p_left, p_right, p_other, w_self, w_other)
    ! Give process p_left access to the a variable of submesh memory of p_right
    ! Used in submesh merging    
    
    IMPLICIT NONE
     
    INTEGER,                            INTENT(IN)        :: p_left   ! Left process ID
    INTEGER,                            INTENT(IN)        :: p_right  ! Right process ID
    INTEGER,                   POINTER, INTENT(INOUT)     :: p_other  ! Pointer to data from other process
    INTEGER,                            INTENT(IN)        :: w_self   ! Window  to data from own process
    INTEGER,                            INTENT(INOUT)     :: w_other  ! Window  to data from other process
    
    INTEGER                                               :: ierr, status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND)                        :: windowsize
    INTEGER                                               :: disp_unit
    TYPE(C_PTR)                                           :: baseptr  
      
    IF (par%i == p_left) THEN
      ! Receive MPI window from p_right
      CALL MPI_RECV( w_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      ! Query this window to get a C baseptr
      CALL MPI_WIN_SHARED_QUERY( w_other, p_right, windowsize, disp_unit, baseptr, ierr)
      ! Convert to fortran pointer
      CALL C_F_POINTER(baseptr, p_other)
    ELSEIF (par%i == p_right) THEN
      ! Send MPI window to p_left
      CALL MPI_SEND( w_self,  1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)      
    END IF ! IF (par%i == p_left) THEN 
    
  END SUBROUTINE ShareMemoryAccess_int_0D
  SUBROUTINE ShareMemoryAccess_int_1D( p_left, p_right, p_other, w_self, w_other, n1)
    ! Give process p_left access to the a variable of submesh memory of p_right
    ! Used in submesh merging    
    
    IMPLICIT NONE
     
    INTEGER,                            INTENT(IN)        :: p_left   ! Left process ID
    INTEGER,                            INTENT(IN)        :: p_right  ! Right process ID
    INTEGER, DIMENSION(:    ), POINTER, INTENT(INOUT)     :: p_other  ! Pointer to data from other process
    INTEGER,                            INTENT(IN)        :: w_self   ! Window  to data from own process
    INTEGER,                            INTENT(INOUT)     :: w_other  ! Window  to data from other process
    INTEGER,                            INTENT(IN)        :: n1
    
    INTEGER                                               :: ierr, status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND)                        :: windowsize
    INTEGER                                               :: disp_unit
    TYPE(C_PTR)                                           :: baseptr  
      
    IF (par%i == p_left) THEN
      ! Receive MPI window from p_right
      CALL MPI_RECV( w_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      ! Query this window to get a C baseptr
      CALL MPI_WIN_SHARED_QUERY( w_other, p_right, windowsize, disp_unit, baseptr, ierr)
      ! Convert to fortran pointer
      CALL C_F_POINTER(baseptr, p_other, [n1])
    ELSEIF (par%i == p_right) THEN
      ! Send MPI window to p_left
      CALL MPI_SEND( w_self,  1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)      
    END IF ! IF (par%i == p_left) THEN 
    
  END SUBROUTINE ShareMemoryAccess_int_1D
  SUBROUTINE ShareMemoryAccess_int_2D( p_left, p_right, p_other, w_self, w_other, n1, n2)
    ! Give process p_left access to the a variable of submesh memory of p_right
    ! Used in submesh merging    
    
    IMPLICIT NONE
     
    INTEGER,                            INTENT(IN)        :: p_left   ! Left process ID
    INTEGER,                            INTENT(IN)        :: p_right  ! Right process ID
    INTEGER, DIMENSION(:,:  ), POINTER, INTENT(INOUT)     :: p_other  ! Pointer to data from other process
    INTEGER,                            INTENT(IN)        :: w_self   ! Window  to data from own process
    INTEGER,                            INTENT(INOUT)     :: w_other  ! Window  to data from other process
    INTEGER,                            INTENT(IN)        :: n1, n2
    
    INTEGER                                               :: ierr, status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND)                        :: windowsize
    INTEGER                                               :: disp_unit
    TYPE(C_PTR)                                           :: baseptr  
      
    IF (par%i == p_left) THEN
      ! Receive MPI window from p_right
      CALL MPI_RECV( w_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      ! Query this window to get a C baseptr
      CALL MPI_WIN_SHARED_QUERY( w_other, p_right, windowsize, disp_unit, baseptr, ierr)
      ! Convert to fortran pointer
      CALL C_F_POINTER(baseptr, p_other, [n1, n2])
    ELSEIF (par%i == p_right) THEN
      ! Send MPI window to p_left
      CALL MPI_SEND( w_self,  1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)      
    END IF ! IF (par%i == p_left) THEN 
    
  END SUBROUTINE ShareMemoryAccess_int_2D
  SUBROUTINE ShareMemoryAccess_int_3D( p_left, p_right, p_other, w_self, w_other, n1, n2, n3)
    ! Give process p_left access to the a variable of submesh memory of p_right
    ! Used in submesh merging    
    
    IMPLICIT NONE
     
    INTEGER,                            INTENT(IN)        :: p_left   ! Left process ID
    INTEGER,                            INTENT(IN)        :: p_right  ! Right process ID
    INTEGER, DIMENSION(:,:,:), POINTER, INTENT(INOUT)     :: p_other  ! Pointer to data from other process
    INTEGER,                            INTENT(IN)        :: w_self   ! Window  to data from own process
    INTEGER,                            INTENT(INOUT)     :: w_other  ! Window  to data from other process
    INTEGER,                            INTENT(IN)        :: n1, n2, n3
    
    INTEGER                                               :: ierr, status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND)                        :: windowsize
    INTEGER                                               :: disp_unit
    TYPE(C_PTR)                                           :: baseptr  
      
    IF (par%i == p_left) THEN
      ! Receive MPI window from p_right
      CALL MPI_RECV( w_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      ! Query this window to get a C baseptr
      CALL MPI_WIN_SHARED_QUERY( w_other, p_right, windowsize, disp_unit, baseptr, ierr)
      ! Convert to fortran pointer
      CALL C_F_POINTER(baseptr, p_other, [n1, n2, n3])
    ELSEIF (par%i == p_right) THEN
      ! Send MPI window to p_left
      CALL MPI_SEND( w_self,  1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)      
    END IF ! IF (par%i == p_left) THEN 
    
  END SUBROUTINE ShareMemoryAccess_int_3D
  SUBROUTINE ShareMemoryAccess_dp_0D(  p_left, p_right, p_other, w_self, w_other)
    ! Give process p_left access to the a variable of submesh memory of p_right
    ! Used in submesh merging    
    
    IMPLICIT NONE
     
    INTEGER,                            INTENT(IN)        :: p_left   ! Left process ID
    INTEGER,                            INTENT(IN)        :: p_right  ! Right process ID
    REAL(dp),                  POINTER, INTENT(INOUT)     :: p_other  ! Pointer to data from other process
    INTEGER,                            INTENT(IN)        :: w_self   ! Window  to data from own process
    INTEGER,                            INTENT(INOUT)     :: w_other  ! Window  to data from other process
    
    INTEGER                                               :: ierr, status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND)                        :: windowsize
    INTEGER                                               :: disp_unit
    TYPE(C_PTR)                                           :: baseptr  
      
    IF (par%i == p_left) THEN
      ! Receive MPI window from p_right
      CALL MPI_RECV( w_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      ! Query this window to get a C baseptr
      CALL MPI_WIN_SHARED_QUERY( w_other, p_right, windowsize, disp_unit, baseptr, ierr)
      ! Convert to fortran pointer
      CALL C_F_POINTER(baseptr, p_other)
    ELSEIF (par%i == p_right) THEN
      ! Send MPI window to p_left
      CALL MPI_SEND( w_self,  1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)      
    END IF ! IF (par%i == p_left) THEN 
    
  END SUBROUTINE ShareMemoryAccess_dp_0D
  SUBROUTINE ShareMemoryAccess_dp_1D(  p_left, p_right, p_other, w_self, w_other, n1)
    ! Give process p_left access to the a variable of submesh memory of p_right
    ! Used in submesh merging    
    
    IMPLICIT NONE
     
    INTEGER,                             INTENT(IN)        :: p_left   ! Left process ID
    INTEGER,                             INTENT(IN)        :: p_right  ! Right process ID
    REAL(dp), DIMENSION(:    ), POINTER, INTENT(INOUT)     :: p_other  ! Pointer to data from other process
    INTEGER,                             INTENT(IN)        :: w_self   ! Window  to data from own process
    INTEGER,                             INTENT(INOUT)     :: w_other  ! Window  to data from other process
    INTEGER,                             INTENT(IN)        :: n1
    
    INTEGER                                                :: ierr, status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND)                         :: windowsize
    INTEGER                                                :: disp_unit
    TYPE(C_PTR)                                            :: baseptr  
      
    IF (par%i == p_left) THEN
      ! Receive MPI window from p_right
      CALL MPI_RECV( w_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      ! Query this window to get a C baseptr
      CALL MPI_WIN_SHARED_QUERY( w_other, p_right, windowsize, disp_unit, baseptr, ierr)
      ! Convert to fortran pointer
      CALL C_F_POINTER(baseptr, p_other, [n1])
    ELSEIF (par%i == p_right) THEN
      ! Send MPI window to p_left
      CALL MPI_SEND( w_self,  1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)      
    END IF ! IF (par%i == p_left) THEN 
    
  END SUBROUTINE ShareMemoryAccess_dp_1D
  SUBROUTINE ShareMemoryAccess_dp_2D(  p_left, p_right, p_other, w_self, w_other, n1, n2)
    ! Give process p_left access to the a variable of submesh memory of p_right
    ! Used in submesh merging    
    
    IMPLICIT NONE
     
    INTEGER,                             INTENT(IN)        :: p_left   ! Left process ID
    INTEGER,                             INTENT(IN)        :: p_right  ! Right process ID
    REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(INOUT)     :: p_other  ! Pointer to data from other process
    INTEGER,                             INTENT(IN)        :: w_self   ! Window  to data from own process
    INTEGER,                             INTENT(INOUT)     :: w_other  ! Window  to data from other process
    INTEGER,                             INTENT(IN)        :: n1, n2
    
    INTEGER                                                :: ierr, status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND)                         :: windowsize
    INTEGER                                                :: disp_unit
    TYPE(C_PTR)                                            :: baseptr  
      
    IF (par%i == p_left) THEN
      ! Receive MPI window from p_right
      CALL MPI_RECV( w_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      ! Query this window to get a C baseptr
      CALL MPI_WIN_SHARED_QUERY( w_other, p_right, windowsize, disp_unit, baseptr, ierr)
      ! Convert to fortran pointer
      CALL C_F_POINTER(baseptr, p_other, [n1, n2])
    ELSEIF (par%i == p_right) THEN
      ! Send MPI window to p_left
      CALL MPI_SEND( w_self,  1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)      
    END IF ! IF (par%i == p_left) THEN 
    
  END SUBROUTINE ShareMemoryAccess_dp_2D
  SUBROUTINE ShareMemoryAccess_dp_3D(  p_left, p_right, p_other, w_self, w_other, n1, n2, n3)
    ! Give process p_left access to the a variable of submesh memory of p_right
    ! Used in submesh merging    
    
    IMPLICIT NONE
     
    INTEGER,                             INTENT(IN)        :: p_left   ! Left process ID
    INTEGER,                             INTENT(IN)        :: p_right  ! Right process ID
    REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(INOUT)     :: p_other  ! Pointer to data from other process
    INTEGER,                             INTENT(IN)        :: w_self   ! Window  to data from own process
    INTEGER,                             INTENT(INOUT)     :: w_other  ! Window  to data from other process
    INTEGER,                             INTENT(IN)        :: n1, n2, n3
    
    INTEGER                                                :: ierr, status(MPI_STATUS_SIZE)
    INTEGER(KIND=MPI_ADDRESS_KIND)                         :: windowsize
    INTEGER                                                :: disp_unit
    TYPE(C_PTR)                                            :: baseptr  
      
    IF (par%i == p_left) THEN
      ! Receive MPI window from p_right
      CALL MPI_RECV( w_other, 1, MPI_INTEGER, p_right, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
      ! Query this window to get a C baseptr
      CALL MPI_WIN_SHARED_QUERY( w_other, p_right, windowsize, disp_unit, baseptr, ierr)
      ! Convert to fortran pointer
      CALL C_F_POINTER(baseptr, p_other, [n1, n2, n3])
    ELSEIF (par%i == p_right) THEN
      ! Send MPI window to p_left
      CALL MPI_SEND( w_self,  1, MPI_INTEGER, p_left, 0, MPI_COMM_WORLD, ierr)      
    END IF ! IF (par%i == p_left) THEN 
    
  END SUBROUTINE ShareMemoryAccess_dp_3D
  
END MODULE parallel_module
