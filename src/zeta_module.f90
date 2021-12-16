MODULE zeta_module

  USE configuration_module,        ONLY: dp, C
  USE data_types_module,           ONLY: type_mesh, type_ice_model
  USE parallel_module,             ONLY: par, sync, ierr, cerr
  
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
    REAL(dp)                                           :: inverse_Hi                 ! Contains the inverse of Hi
    
    DO vi = mesh%v1, mesh%v2
      inverse_Hi = 1._dp / ice%Hi_a( vi)
      ice%dzeta_dz_a(vi) = -inverse_Hi                                                      
      DO k = 1, C%NZ
        ice%dzeta_dt_a( vi,k)  =  inverse_Hi * (ice%dHs_dt_a( vi)  - C%zeta(k) * ice%dHi_dt_a( vi))
        ice%dzeta_dx_a( vi,k)  =  inverse_Hi * (ice%dHs_dx_a( vi)  - C%zeta(k) * ice%dHi_dx_a( vi))
        ice%dzeta_dy_a( vi,k)  =  inverse_Hi * (ice%dHs_dy_a( vi)  - C%zeta(k) * ice%dHi_dy_a( vi))
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
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
