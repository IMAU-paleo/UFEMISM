MODULE zeta_module

  ! Import basic functionality
  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parameters_module
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, write_to_memory_log, &
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
  
  ! Import specific functionality 
  USE data_types_module,               ONLY: type_mesh, type_ice_model
  USE mesh_operators_module,           ONLY: ddx_a_to_a_2D, ddy_a_to_a_2D
  
  IMPLICIT NONE
  
  TYPE type_zeta_parameters
    ! Some useful constants for the scaled vertical coordinate transformation
    
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: c_zeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: c_zetazeta
    REAL(dp), DIMENSION(:), ALLOCATABLE :: a_zeta_minus
    REAL(dp), DIMENSION(:), ALLOCATABLE :: b_zeta_minus
    REAL(dp), DIMENSION(:), ALLOCATABLE :: z_zeta_minus
    
    REAL(dp) :: a_N
    REAL(dp) :: b_1
    REAL(dp) :: b_t
    REAL(dp) :: s_t
  
  END TYPE type_zeta_parameters
  
  TYPE( type_zeta_parameters), SAVE :: p_zeta
  
CONTAINS

  SUBROUTINE calculate_zeta_derivatives( mesh, ice)
    ! This subroutine calculates the global struct which contains the derivatives of 
    ! zeta, which are used in the transformation to the t,x,y,zeta coordinates. 
    ! This zeta derivatives are the Jacobians.    
    
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_ice_model),                INTENT(INOUT) :: ice

    ! Local variables:
    INTEGER                                            :: vi, k
    REAL(dp), DIMENSION(:    ), POINTER                ::  dHs_dt_a,  dHs_dx_a,  dHs_dy_a,  dHi_dx_a,  dHi_dy_a
    INTEGER                                            :: wdHs_dt_a, wdHs_dx_a, wdHs_dy_a, wdHi_dx_a, wdHi_dy_a
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, dHs_dt_a, wdHs_dt_a)
    CALL allocate_shared_dp_1D( mesh%nV, dHs_dx_a, wdHs_dx_a)
    CALL allocate_shared_dp_1D( mesh%nV, dHs_dy_a, wdHs_dy_a)
    CALL allocate_shared_dp_1D( mesh%nV, dHi_dx_a, wdHi_dx_a)
    CALL allocate_shared_dp_1D( mesh%nV, dHi_dy_a, wdHi_dy_a)
    
    ! Calculate dHs_dt
    DO vi = mesh%vi1, mesh%vi2
      IF (ice%mask_land_a( vi) == 1) THEN
        dHs_dt_a( vi) = ice%dHb_dt_a( vi) + ice%dHi_dt_a( vi)
      ELSE
        dHs_dt_a( vi) = (1._dp - ice_density / seawater_density) * ice%dHi_dt_a( vi)
      END IF
    END DO
    CALL sync
    
    ! Calculate spatial derivatives of Hi and Hs
    CALL ddx_a_to_a_2D( mesh, ice%Hs_a, dHs_dx_a)
    CALL ddy_a_to_a_2D( mesh, ice%Hs_a, dHs_dy_a)
    CALL ddx_a_to_a_2D( mesh, ice%Hi_a, dHi_dx_a)
    CALL ddy_a_to_a_2D( mesh, ice%Hi_a, dHi_dy_a)
    
    ! Calculate derivatives of zeta
    DO vi = mesh%vi1, mesh%vi2
    
      ice%dzeta_dz_a(vi) = -1._dp / ice%Hi_a( vi)   
                                                         
      DO k = 1, C%NZ
        ice%dzeta_dt_a( vi,k)  =  (dHs_dt_a( vi)  - C%zeta(k) * ice%dHi_dt_a( vi)) / ice%Hi_a( vi)
        ice%dzeta_dx_a( vi,k)  =  (dHs_dx_a( vi)  - C%zeta(k) *     dHi_dx_a( vi)) / ice%Hi_a( vi)
        ice%dzeta_dy_a( vi,k)  =  (dHs_dy_a( vi)  - C%zeta(k) *     dHi_dy_a( vi)) / ice%Hi_a( vi)
      END DO
      
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wdHs_dt_a)
    CALL deallocate_shared( wdHs_dx_a)
    CALL deallocate_shared( wdHs_dy_a)
    CALL deallocate_shared( wdHi_dx_a)
    CALL deallocate_shared( wdHi_dy_a)
    
  END SUBROUTINE calculate_zeta_derivatives
  SUBROUTINE initialize_zeta_discretization
    ! See table 3 and table 1.
    ! To be called in anice, not only usefull for uvw-mumps, but also for temperature
    ! Calculating the discretization coefficents for a Arakawa A-grid discretization    
    
    IMPLICIT NONE
    
    ! Local variables:
    INTEGER                       :: k
    REAL(dp), DIMENSION(2:C%NZ)   :: a_k
    REAL(dp), DIMENSION(1:C%NZ-1) :: b_k
    REAL(dp), DIMENSION(3:C%NZ)   :: c_k
    REAL(dp), DIMENSION(1:C%NZ-2) :: d_k

    ALLOCATE(p_zeta%a_zeta(2:C%NZ-1))
    ALLOCATE(p_zeta%b_zeta(2:C%NZ-1))
    ALLOCATE(p_zeta%c_zeta(2:C%NZ-1))
    ALLOCATE(p_zeta%a_zetazeta(2:C%NZ-1))
    ALLOCATE(p_zeta%b_zetazeta(2:C%NZ-1))
    ALLOCATE(p_zeta%c_zetazeta(2:C%NZ-1))
    
    ALLOCATE(p_zeta%z_zeta_minus(3:C%NZ))
    ALLOCATE(p_zeta%a_zeta_minus(3:C%NZ))
    ALLOCATE(p_zeta%b_zeta_minus(3:C%NZ))
    
    DO k = 2, C%NZ
     a_k(k) = C%zeta(k)   - C%zeta(k-1)
    END DO
    DO k = 1, C%NZ-1
     b_k(k) = C%zeta(k+1) - C%zeta(k)
    END DO
    DO k = 3, C%NZ
     c_k(k) = C%zeta(k)   - C%zeta(k-2)
    END DO
    DO k = 1, C%NZ-2
     d_k(k) = C%zeta(k+2) - C%zeta(k)
    END DO
    p_zeta%a_N = a_k(C%NZ)
    p_zeta%b_1 = b_k(1)
    
    DO k = 2, C%NZ-1
     p_zeta%a_zeta(k)     =         - b_k(k)  / (a_k(k) * (a_k(k) + b_k(k))) 
     p_zeta%b_zeta(k)     = (b_k(k) - a_k(k)) / (a_k(k) *  b_k(k))  
     p_zeta%c_zeta(k)     =           a_k(k)  / (b_k(k) * (a_k(k) + b_k(k)))
     p_zeta%a_zetazeta(k) =           2._dp   / (a_k(k) * (a_k(k) + b_k(k)))
     p_zeta%b_zetazeta(k) =         - 2._dp   / (a_k(k) *  b_k(k))       
     p_zeta%c_zetazeta(k) =           2._dp   / (b_k(k) * (a_k(k) + b_k(k)))
    END DO

    ! Not all of these are in use:
    DO k = 3, C%NZ
      p_zeta%z_zeta_minus(k) =           a_k(k)  / (c_k(k) * (c_k(k) - a_k(k)))   
      p_zeta%a_zeta_minus(k) =           c_k(k)  / (a_k(k) * (a_k(k) - c_k(k)))   
      p_zeta%b_zeta_minus(k) = (a_k(k) + c_k(k)) / (a_k(k) * c_k(k))          
    END DO
    
    ! Discretization coefficients in time
    p_zeta%b_t = -1._dp / C%dt_thermo                      ! See equations 14.38
    p_zeta%s_t =  1._dp / C%dt_thermo                      ! See equations 14.39
    
  END SUBROUTINE initialize_zeta_discretization

END MODULE zeta_module
