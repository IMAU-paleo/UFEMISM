module zeta_module

  use mpi
  use configuration_module,  only : dp, C, routine_path, init_routine, finalise_routine, crash, warning
  use parameters_module,     only : ice_density, seawater_density
  use data_types_module,     only : type_mesh, type_ice_model
  use mesh_operators_module, only : ddx_a_to_a_2D, ddy_a_to_a_2D

  implicit none

  type type_zeta_parameters
    ! Some useful constants for the scaled vertical coordinate transformation

    real(dp), dimension(:), allocatable :: a_zeta
    real(dp), dimension(:), allocatable :: b_zeta
    real(dp), dimension(:), allocatable :: c_zeta
    real(dp), dimension(:), allocatable :: a_zetazeta
    real(dp), dimension(:), allocatable :: b_zetazeta
    real(dp), dimension(:), allocatable :: c_zetazeta
    real(dp), dimension(:), allocatable :: a_zeta_minus
    real(dp), dimension(:), allocatable :: b_zeta_minus
    real(dp), dimension(:), allocatable :: z_zeta_minus

    real(dp) :: a_N
    real(dp) :: b_1
    real(dp) :: b_t
    real(dp) :: s_t

  end type type_zeta_parameters

  type( type_zeta_parameters), save :: p_zeta

contains

  subroutine calculate_zeta_derivatives( mesh, ice)
    ! Derivatives of zeta, which are used in the transformation to the
    ! t,x,y,zeta coordinates. These zeta derivatives are the Jacobians.

    implicit none

    ! In- and output variables
    type(type_mesh),      intent(in)    :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    integer                             :: vi, k
    real(dp), dimension(:), allocatable :: dHs_dt_a, dHs_dx_a, dHs_dy_a, dHi_dx_a, dHi_dy_a
    real(dp)                            :: Hi_a

    ! Allocate shared memory
    allocate( dHs_dt_a (mesh%vi1:mesh%vi2) )
    allocate( dHs_dx_a (mesh%vi1:mesh%vi2) )
    allocate( dHs_dy_a (mesh%vi1:mesh%vi2) )
    allocate( dHi_dx_a (mesh%vi1:mesh%vi2) )
    allocate( dHi_dy_a (mesh%vi1:mesh%vi2) )

    ! Calculate dHs_dt
    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_land_a( vi) == 1) then
        dHs_dt_a( vi) = ice%dHb_dt_a( vi) + ice%dHi_dt_a( vi)
      else
        dHs_dt_a( vi) = (1._dp - ice_density / seawater_density) * ice%dHi_dt_a( vi)
      end if
    end do

    ! Calculate spatial derivatives of Hi and Hs
    call ddx_a_to_a_2D( mesh, ice%Hs_a, dHs_dx_a)
    call ddy_a_to_a_2D( mesh, ice%Hs_a, dHs_dy_a)
    call ddx_a_to_a_2D( mesh, ice%Hi_a, dHi_dx_a)
    call ddy_a_to_a_2D( mesh, ice%Hi_a, dHi_dy_a)

    ! Calculate derivatives of zeta
    do vi = mesh%vi1, mesh%vi2
      Hi_a = max(ice%Hi_a(vi),0.00001) ! Cannot divide by zero 

      ice%dzeta_dz_a(vi) = -1._dp / Hi_a

      do k = 1, C%nz
        ice%dzeta_dt_a( vi,k)  =  (dHs_dt_a( vi)  - C%zeta(k) * ice%dHi_dt_a( vi)) / Hi_a
        ice%dzeta_dx_a( vi,k)  =  (dHs_dx_a( vi)  - C%zeta(k) *     dHi_dx_a( vi)) / Hi_a
        ice%dzeta_dy_a( vi,k)  =  (dHs_dy_a( vi)  - C%zeta(k) *     dHi_dy_a( vi)) / Hi_a
      end do

    end do

    ! Clean up after yourself
    deallocate( dHs_dt_a)
    deallocate( dHs_dx_a)
    deallocate( dHs_dy_a)
    deallocate( dHi_dx_a)
    deallocate( dHi_dy_a)

  end subroutine calculate_zeta_derivatives

  subroutine initialise_zeta_discretisation
    ! Discretization coefficents for an Arakawa A-grid

    implicit none

    ! Local variables:
    integer                       :: k
    real(dp), dimension(2:C%NZ)   :: a_k
    real(dp), dimension(1:C%NZ-1) :: b_k
    real(dp), dimension(3:C%NZ)   :: c_k
    real(dp), dimension(1:C%NZ-2) :: d_k

    allocate(p_zeta%a_zeta(2:C%NZ-1))
    allocate(p_zeta%b_zeta(2:C%NZ-1))
    allocate(p_zeta%c_zeta(2:C%NZ-1))
    allocate(p_zeta%a_zetazeta(2:C%NZ-1))
    allocate(p_zeta%b_zetazeta(2:C%NZ-1))
    allocate(p_zeta%c_zetazeta(2:C%NZ-1))

    allocate(p_zeta%z_zeta_minus(3:C%NZ))
    allocate(p_zeta%a_zeta_minus(3:C%NZ))
    allocate(p_zeta%b_zeta_minus(3:C%NZ))

    do k = 2, C%NZ
     a_k(k) = C%zeta(k)   - C%zeta(k-1)
    end do
    do k = 1, C%NZ-1
     b_k(k) = C%zeta(k+1) - C%zeta(k)
    end do
    do k = 3, C%NZ
     c_k(k) = C%zeta(k)   - C%zeta(k-2)
    end do
    do k = 1, C%NZ-2
     d_k(k) = C%zeta(k+2) - C%zeta(k)
    end do
    p_zeta%a_N = a_k(C%NZ)
    p_zeta%b_1 = b_k(1)

    do k = 2, C%NZ-1
     p_zeta%a_zeta(k)     =         - b_k(k)  / (a_k(k) * (a_k(k) + b_k(k)))
     p_zeta%b_zeta(k)     = (b_k(k) - a_k(k)) / (a_k(k) *  b_k(k))
     p_zeta%c_zeta(k)     =           a_k(k)  / (b_k(k) * (a_k(k) + b_k(k)))
     p_zeta%a_zetazeta(k) =           2._dp   / (a_k(k) * (a_k(k) + b_k(k)))
     p_zeta%b_zetazeta(k) =         - 2._dp   / (a_k(k) *  b_k(k))
     p_zeta%c_zetazeta(k) =           2._dp   / (b_k(k) * (a_k(k) + b_k(k)))
    end do

    ! Not all of these are in use:
    do k = 3, C%NZ
      p_zeta%z_zeta_minus(k) =           a_k(k)  / (c_k(k) * (c_k(k) - a_k(k)))
      p_zeta%a_zeta_minus(k) =           c_k(k)  / (a_k(k) * (a_k(k) - c_k(k)))
      p_zeta%b_zeta_minus(k) = (a_k(k) + c_k(k)) / (a_k(k) * c_k(k))
    end do

    ! Discretization coefficients in time
    p_zeta%b_t = -1._dp / C%dt_thermo
    p_zeta%s_t =  1._dp / C%dt_thermo

  end subroutine initialise_zeta_discretisation

end module zeta_module
