MODULE utilities_module

  ! Some generally useful tools

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>

  USE mpi
  USE configuration_module,          ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                  ONLY: perr
  USE parallel_module,               ONLY: par, sync, ierr, cerr, partition_list
  USE data_types_module,             ONLY: type_mesh, type_grid, type_model_region

  implicit none
  ! Interfaces to LAPACK, which are otherwise implicitly generated (taken from
  ! LAPACK source)
  !  *
  !  *  -- LAPACK routine (version 3.1) --
  !  *     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
  !  *     November 2006
  interface
    SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
      INTEGER            INFO, LDA, M, N
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
    END SUBROUTINE
    SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
      INTEGER            INFO, LDA, LWORK, N
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( * )
    END SUBROUTINE
    SUBROUTINE dgtsv( N, NRHS, DL, D, DU, B, LDB, INFO )
      INTEGER            INFO, LDB, N, NRHS
      DOUBLE PRECISION   B( LDB, * ), D( * ), DL( * ), DU( * )
    END SUBROUTINE
  end interface

  interface check_for_nan
    procedure :: check_for_NaN_dp_1D
    procedure :: check_for_NaN_dp_2D
    procedure :: check_for_NaN_dp_3D
    procedure :: check_for_NaN_int_1D
    procedure :: check_for_NaN_int_2D
    procedure :: check_for_NaN_int_3D
  end interface
CONTAINS

! == Some operations on the scaled vertical coordinate
  SUBROUTINE vertical_integration_from_bottom_to_zeta( f, integral_f)
    ! This subroutine integrates f from the bottom level at C%zeta(k=C%nz) = 1 up to the level C%zeta(k):
    !  See Eq. (12.1)
    ! If the integrand f is positive (our case) the integral is negative because the integration is in
    ! the opposite zeta direction. A 1D array which contains for each k-layer the integrated value from
    ! the bottom up to that k-layer is returned. The value of the integrand f at some integration step k
    ! is the average of f(k+1) and f(k):
    !  integral_f(k) = integral_f(k+1) + 0.5 * (f(k+1) + f(k)) * (-dzeta)
    ! with dzeta = C%zeta(k+1) - C%zeta(k). So for f > 0  integral_f < 0.

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%nz), INTENT(IN)  :: f
    REAL(dp), DIMENSION(C%nz), INTENT(OUT) :: integral_f

    ! Local variables:
    INTEGER                                :: k

    integral_f(C%nz) = 0._dp
    DO k = C%nz-1, 1, -1
      integral_f(k) = integral_f(k+1) - 0.5_dp * (f(k+1) + f(k)) * (C%zeta(k+1) - C%zeta(k))
    END DO

  END SUBROUTINE vertical_integration_from_bottom_to_zeta
  SUBROUTINE vertical_integration_from_top_to_zeta(    f, integral_f)
    ! This subroutine integrates f from the top level at C%zeta(k=1) = 0 down to the level C%zeta(k): Eq. (12.2)
    ! Similar to Eq. (12.1) but in the other direction.
    ! If the integrand f is positive (our case) the integral is positive because the integration is in
    ! the zeta direction. A 1D array which contains for each k-layer the integrated value from
    ! the top down to that k-layer is returned. The value of the integrand f at some integration step k
    ! is the average of f(k) and f(k-1):
    ! integral_f(k) = integral_f(k-1) + 0.5 * (f(k) + f(k-1)) * (dzeta); with dzeta = C%zeta(k+1) - C%zeta(k).
    ! Heiko Goelzer (h.goelzer@uu.nl) Jan 2016

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%nz), INTENT(IN)  :: f
    REAL(dp), DIMENSION(C%nz), INTENT(OUT) :: integral_f

    ! Local variables:
    INTEGER                                :: k

    integral_f(1) = 0._dp
    DO k = 2, C%nz, 1
      integral_f(k) = integral_f(k-1) + 0.5_dp * (f(k) + f(k-1)) * (C%zeta(k) - C%zeta(k-1))
    END DO

  END SUBROUTINE vertical_integration_from_top_to_zeta
  SUBROUTINE vertical_integrate(                       f, integral_f)
    ! Integrate f over the ice column (from the base to the surface)

    IMPLICIT NONE

    ! Input variable:
    REAL(dp), DIMENSION(C%nz), INTENT(IN) :: f
    REAL(dp),                  INTENT(OUT):: integral_f

    ! Local variable:
    INTEGER                               :: k

    ! Initial value is zero
    integral_f = 0.0_dp

    ! Intermediate values include sum of all previous values
    ! Take current value as average between points
    DO k = 2, C%nz
       integral_f = integral_f + 0.5_dp*(f(k)+f(k-1))*(C%zeta(k) - C%zeta(k-1))
    END DO

  END SUBROUTINE vertical_integrate
  SUBROUTINE vertical_average(                         f, average_f)
    ! Calculate the vertical average of any given function f defined at the vertical zeta grid.
    !  See Eq. (11.3) in DOCUMENTATION/icedyn-documentation/icedyn-documentation.pdf.
    ! The integration is in the direction of the positive zeta-axis from C%zeta(k=1) = 0 up to C%zeta(k=C%nz) = 1.
    ! Numerically: de average between layer k and k+1 is calculated and multiplied by the distance between those
    ! layers k and k+1, which is imediately the weight factor for this contribution because de total layer distance
    ! is scaled to 1. The sum of all the weighted contribution gives average_f the vertical average of f.

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(C%nz), INTENT(IN) :: f
    REAL(dp),                  INTENT(OUT):: average_f

    ! Local variables:
    INTEGER                               :: k

    !  See Eq. (11.4) in DOCUMENTATION/icedyn-documentation/icedyn-documentation.pdf
    average_f = 0._dp
    DO k = 1, C%nz-1
      average_f = average_f + 0.5_dp * (f(k+1) + f(k)) * (C%zeta(k+1) - C%zeta(k))
    END DO

  END SUBROUTINE vertical_average

! == Floatation criterion, surface elevation, and thickness above floatation
  FUNCTION is_floating( Hi, Hb, SL) RESULT( isso)
    ! The flotation criterion

    IMPLICIT NONE

    REAL(dp),                            INTENT(IN)    :: Hi, Hb, SL
    LOGICAL                                            :: isso

    isso = .FALSE.
    IF (Hi < (SL - Hb) * seawater_density/ice_density) isso = .TRUE.

  END FUNCTION is_floating
  FUNCTION surface_elevation( Hi, Hb, SL) RESULT( Hs)
    ! The surface elevation equation

    IMPLICIT NONE

    REAL(dp),                            INTENT(IN)    :: Hi, Hb, SL
    REAL(dp)                                           :: Hs

    Hs = Hi + MAX( SL - ice_density / seawater_density * Hi, Hb)

  END FUNCTION surface_elevation
  FUNCTION thickness_above_floatation( Hi, Hb, SL) RESULT( TAF)
    ! The thickness-above-floatation equation

    IMPLICIT NONE

    REAL(dp),                            INTENT(IN)    :: Hi, Hb, SL
    REAL(dp)                                           :: TAF

    TAF = Hi - MAX(0._dp, (SL - Hb) * (seawater_density / ice_density))

  END FUNCTION thickness_above_floatation

! == The error function (used in the Roe&Lindzen precipitation model)
  subroutine error_function(X, ERR)
    ! Purpose: Compute error function erf(x)
    ! Input:   x   --- Argument of erf(x)
    ! Output:  ERR --- erf(x)

    implicit none

    ! Input/Output variables:
    real(dp), intent(in)  :: X
    real(dp), intent(out) :: ERR

    ! Local variables:
    real(dp)              :: EPS
    real(dp)              :: X2
    real(dp)              :: ER
    real(dp)              :: R
    real(dp)              :: C0
    integer               :: k

    EPS = 1.0E-15_dp
    X2  = X * X
    if ( abs(X) < 3.5_dp) then
     ER = 1.0_dp
     R  = 1.0_dp
     do k = 1, 50
       R  = R * X2 / (real(k, dp) + 0.5_dp)
       ER = ER+R
       if ( abs(R) < abs(ER) * EPS) then
        C0  = 2.0_dp / sqrt(pi) * X * exp(-X2)
        ERR = C0 * ER
        exit
       end if
     end do
    else
     ER = 1.0_dp
     R  = 1.0_dp
     do k = 1, 12
       R  = -R * (real(k, dp) - 0.5_dp) / X2
       ER = ER + R
       C0  = exp(-X2) / (abs(X) * sqrt(pi))
       ERR = 1.0_dp - C0 * ER
       if ( X < 0.0_dp) ERR = -ERR
     end do
    end if

  end subroutine error_function

! ===== The oblique stereographic projection =====
! ================================================

  SUBROUTINE get_lat_lon_coordinates( mesh)
   ! Use the inverse stereographic projection for the mesh model region to calculate
   ! lat/lon coordinates for all the vertices

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),          INTENT(INOUT)       :: mesh

    ! Local variables
    INTEGER                                       :: vi

    ! Calculate lat and lon directly from X and Y using inverse projection
    DO vi = mesh%vi1, mesh%vi2
      CALL inverse_oblique_sg_projection(mesh%V(vi,1), mesh%V(vi,2), mesh%lambda_M, &
                                         mesh%phi_M, mesh%alpha_stereo, mesh%lon(vi), &
                                         mesh%lat(vi))
    END DO

  END SUBROUTINE get_lat_lon_coordinates

  SUBROUTINE oblique_sg_projection(lambda, phi, lambda_M_deg, phi_M_deg, alpha_deg, x_IM_P_prime, y_IM_P_prime, k_P)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates to a rectangular coordinate system, with coordinates (x,y).
    !
    ! For more information about M, alpha_deg, the center of projection and the used
    ! projection method see: Reerink et al. (2010), Mapping technique of climate fields
    ! between GCM's and ice models, GMD

    ! For North and South Pole: lambda_M_deg = 0._dp, to generate the correct coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)            :: lambda        ! in degrees
    REAL(dp), INTENT(IN)            :: phi           ! in degrees

    ! Polar stereographic projection parameters
    REAL(dp), INTENT(IN)            :: lambda_M_deg  ! in degrees
    REAL(dp), INTENT(IN)            :: phi_M_deg     ! in degrees
    REAL(dp), INTENT(IN)            :: alpha_deg     ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT)           :: x_IM_P_prime  ! in metres
    REAL(dp), INTENT(OUT)           :: y_IM_P_prime  ! in metres
    REAL(dp), INTENT(OUT), OPTIONAL :: k_P           ! Length scale factor [-],  k in Snyder (1987)

    ! Local variables:
    REAL(dp)                        :: phi_P         ! in radians
    REAL(dp)                        :: lambda_P      ! in radians
    REAL(dp)                        :: t_P_prime
    REAL(dp)                        :: lambda_M, phi_M, alpha

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = (pi / 180._dp) * phi
    lambda_P = (pi / 180._dp) * lambda

    ! Convert projection parameters to radians:
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = (1._dp + COS(alpha)) / (1._dp + COS(phi_P) * COS(phi_M) * COS(lambda_P - lambda_M) + SIN(phi_P) * SIN(phi_M))

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x_IM_P_prime =  earth_radius * (COS(phi_P) * SIN(lambda_P - lambda_M)) * t_P_prime
    y_IM_P_prime =  earth_radius * (SIN(phi_P) * COS(phi_M) - (COS(phi_P) * SIN(phi_M)) * COS(lambda_P - lambda_M)) * t_P_prime

    ! See equation (21-4) on page 157 in Snyder (1987):
    IF(PRESENT(k_P)) k_P = (1._dp + COS(alpha)) / (1._dp + SIN(phi_M) * SIN(phi_P) + COS(phi_M) * COS(phi_P) * COS(lambda_P - lambda_M))

  END SUBROUTINE oblique_sg_projection

  SUBROUTINE inverse_oblique_sg_projection(x_IM_P_prime, y_IM_P_prime, lambda_M_deg, phi_M_deg, alpha_deg, lambda_P, phi_P)
    ! This subroutine projects with an inverse oblique stereographic projection the
    ! (x,y) coordinates to a longitude-latitude coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, alpha_deg, the center of projection and the used
    ! projection method see: Reerink et al. (2010), Mapping technique of climate fields
    ! between GCM's and ice models, GMD

    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in metres
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in metres

    ! Polar stereographic projection parameters
    REAL(dp), INTENT(IN)  :: lambda_M_deg  ! in degrees
    REAL(dp), INTENT(IN)  :: phi_M_deg     ! in degrees
    REAL(dp), INTENT(IN)  :: alpha_deg     ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: x_3D_P_prime  ! in metres
    REAL(dp)              :: y_3D_P_prime  ! in metres
    REAL(dp)              :: z_3D_P_prime  ! in metres
    REAL(dp)              :: a
    REAL(dp)              :: t_P
    REAL(dp)              :: x_3D_P        ! in metres
    REAL(dp)              :: y_3D_P        ! in metres
    REAL(dp)              :: z_3D_P        ! in metres
    REAL(dp)              :: lambda_M, phi_M, alpha

    ! Convert projection parameters to radians:
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = earth_radius * COS(alpha) * COS(lambda_M) * COS(phi_M) - SIN(lambda_M) * x_IM_P_prime - COS(lambda_M) * SIN(phi_M) * y_IM_P_prime
    y_3D_P_prime = earth_radius * COS(alpha) * SIN(lambda_M) * COS(phi_M) + COS(lambda_M) * x_IM_P_prime - SIN(lambda_M) * SIN(phi_M) * y_IM_P_prime
    z_3D_P_prime = earth_radius * COS(alpha) *                 SIN(phi_M)                                +                 COS(phi_M) * y_IM_P_prime

    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = COS(lambda_M) * COS(phi_M) * x_3D_P_prime  +  SIN(lambda_M) * COS(phi_M) * y_3D_P_prime  +  SIN(phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * earth_radius**2 + 2._dp * earth_radius * a) / (earth_radius**2 + 2._dp * earth_radius * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  earth_radius * COS(lambda_M) * COS(phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  earth_radius * SIN(lambda_M) * COS(phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  earth_radius *                 SIN(phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

    ! See equation (2.7) or equation (B.24) in Reerink et al. (2010):
    IF(x_3D_P <  0._dp                      ) THEN
      lambda_P = 180._dp + (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P >= 0._dp) THEN
      lambda_P =           (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P >  0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 360._dp + (180._dp / pi) * ATAN(y_3D_P / x_3D_P)
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P >  0._dp) THEN
      lambda_P =  90._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P <  0._dp) THEN
      lambda_P = 270._dp
    ELSE IF(x_3D_P == 0._dp .AND. y_3D_P == 0._dp) THEN
      lambda_P =   0._dp
    END IF

    ! See equation (2.8) or equation (B.25) in Reerink et al. (2010):
    IF(x_3D_P /= 0._dp .OR. y_3D_P /= 0._dp) THEN
      phi_P = (180._dp / pi) * ATAN(z_3D_P / sqrt(x_3D_P**2 + y_3D_P**2))
    ELSE IF(z_3D_P >  0._dp) THEN
      phi_P =   90._dp
    ELSE IF(z_3D_P <  0._dp) THEN
      phi_P =  -90._dp
    END IF

  END SUBROUTINE inverse_oblique_sg_projection

! ! == Map data between two square grids using 2nd-order conservative remapping
!   SUBROUTINE map_square_to_square_cons_2nd_order_2D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src, d_dst)
!     ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)

!     IMPLICIT NONE

!     ! Input and output variables
!     INTEGER,                            INTENT(IN)    :: nx_src
!     INTEGER,                            INTENT(IN)    :: ny_src
!     REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_src
!     REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_src
!     INTEGER,                            INTENT(IN)    :: nx_dst
!     INTEGER,                            INTENT(IN)    :: ny_dst
!     REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_dst
!     REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_dst
!     REAL(dp), DIMENSION(:,:  ),         INTENT(IN)    :: d_src
!     REAL(dp), DIMENSION(:,:  ),         INTENT(OUT)   :: d_dst

!     ! Local variables:
!     CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_square_to_square_cons_2nd_order_2D'
!     INTEGER                                           :: i,j,i_src,j_src,i1,i2,igmin,igmax,jgmin,jgmax,j1,j2
!     REAL(dp)                                          :: dx_src, dy_src, dx_dst, dy_dst, xcmin, xcmax, ycmin, ycmax
!     INTEGER,  DIMENSION(nx_dst,2)                     :: ir_src
!     INTEGER,  DIMENSION(ny_dst,2)                     :: jr_src
!     REAL(dp)                                          :: xomin, xomax, yomin, yomax, w0, w1x, w1y
!     REAL(dp)                                          :: Ad, Asd, Asum
!     REAL(dp), DIMENSION(:,:  ), POINTER               ::  ddx_src,  ddy_src
!     INTEGER                                           :: wddx_src, wddy_src
!     INTEGER,  DIMENSION(:,:  ), POINTER               ::  mask_dst_outside_src
!     INTEGER                                           :: wmask_dst_outside_src
!     REAL(dp)                                          :: LI_mxydx1, LI_mxydx2, LI_mxydx3, LI_mxydx4
!     REAL(dp)                                          :: LI_xydy1, LI_xydy2, LI_xydy3, LI_xydy4

!     ! Allocate shared memory
!     CALL allocate_shared_dp_2D(  ny_src, nx_src, ddx_src,              wddx_src             )
!     CALL allocate_shared_dp_2D(  ny_src, nx_src, ddy_src,              wddy_src             )
!     CALL allocate_shared_int_2D( ny_dst, nx_dst, mask_dst_outside_src, wmask_dst_outside_src)

!     ! Find grid spacings
!     dx_src = x_src(2) - x_src(1)
!     dy_src = y_src(2) - y_src(1)
!     dx_dst = x_dst(2) - x_dst(1)
!     dy_dst = y_dst(2) - y_dst(1)
!     Ad = dx_dst * dy_dst

!     ! If the grids are equal, the solution is trivial; just copy the data
!     IF (dx_src == dx_dst .AND. dy_src == dy_dst .AND. nx_src == nx_dst .AND. ny_src == ny_dst) THEN
!       CALL partition_list( nx_dst, par%i, par%n, i1, i2)
!       d_dst( :,i1:i2) = d_src( :,i1:i2)
!       CALL sync
!       RETURN
!     END IF

!     ! Find overlaps between grids
!     DO i = 1, nx_dst
!       ! Dst cell i overlaps with src cells ir_src( i,1) to ir_src( i,2)
!       xcmin = x_dst( i) - dx_dst/2._dp
!       xcmax = x_dst( i) + dx_dst/2._dp
!       ir_src( i,:) = MAX( 1, MIN( nx_src, [CEILING(-1.5_dp + FLOOR(nx_src/2._dp) + xcmin / dx_src), &
!                                            CEILING( 1.5_dp + FLOOR(nx_src/2._dp) + xcmax / dx_src)] ))
!     END DO ! DO i = 1, nx_dst
!     DO j = 1, ny_dst
!       ! Dst cell j overlaps with src cells jr_src( j,1) to ir_src( j,2)
!       ycmin = y_dst( j) - dy_dst/2._dp
!       ycmax = y_dst( j) + dy_dst/2._dp
!       jr_src( j,:) = MAX( 1, MIN( ny_src, [CEILING(-1.5_dp + FLOOR(ny_src/2._dp) + ycmin / dy_src), &
!                                            CEILING( 1.5_dp + FLOOR(ny_src/2._dp) + ycmax / dy_src)] ))
!     END DO ! DO j = 1, ny_dst

!     ! Get derivatives of d_src
!     CALL partition_list( nx_src, par%i, par%n, i1, i2)
!     DO i = MAX(2,i1), MIN(nx_src-1,i2)
!     DO j = 2, ny_src-1
!       ddx_src( j,i) = (d_src( j,i+1) - d_src( j,i-1)) / (2._dp * dx_src)
!       ddy_src( j,i) = (d_src( j+1,i) - d_src( j-1,i)) / (2._dp * dy_src)
!     END DO
!     END DO
!     CALL sync

!     ! Find parallelisation domains
!     CALL partition_list( nx_dst, par%i, par%n, i1, i2)
!     CALL partition_list( ny_dst, par%i, par%n, j1, j2)

!     DO i = i1, i2
!     DO j = 1, ny_dst

!       d_dst(                j,i) = 0._dp
!       mask_dst_outside_src( j,i) = 0
!       Asum                       = 0._dp

!       DO i_src = ir_src( i,1), ir_src( i,2)
!       DO j_src = jr_src( j,1), jr_src( j,2)

!         xomin = MAX( x_dst( i) - dx_dst/2._dp, x_src( i_src) - dx_src/2._dp)
!         xomax = MIN( x_dst( i) + dx_dst/2._dp, x_src( i_src) + dx_src/2._dp)
!         yomin = MAX( y_dst( j) - dy_dst/2._dp, y_src( j_src) - dy_src/2._dp)
!         yomax = MIN( y_dst( j) + dy_dst/2._dp, y_src( j_src) + dy_src/2._dp)

!         IF (xomax <= xomin .OR. yomax <= yomin) CYCLE

!         Asd  = (xomax - xomin) * (yomax - yomin)
!         Asum = Asum + Asd

!         w0  = Asd / Ad

!         CALL line_integral_mxydx( [xomin,yomin], [xomax,yomin], 1E-9_dp, LI_mxydx1)
!         CALL line_integral_mxydx( [xomax,yomin], [xomax,yomax], 1E-9_dp, LI_mxydx2)
!         CALL line_integral_mxydx( [xomax,yomax], [xomin,yomax], 1E-9_dp, LI_mxydx3)
!         CALL line_integral_mxydx( [xomin,yomax], [xomin,yomin], 1E-9_dp, LI_mxydx4)

!         w1x = 1._dp / Ad * (LI_mxydx1 + LI_mxydx2 + LI_mxydx3 + LI_mxydx4) - w0 * x_src( i_src)

!         CALL line_integral_xydy(  [xomin,yomin], [xomax,yomin], 1E-9_dp, LI_xydy1)
!         CALL line_integral_xydy(  [xomax,yomin], [xomax,yomax], 1E-9_dp, LI_xydy2)
!         CALL line_integral_xydy(  [xomax,yomax], [xomin,yomax], 1E-9_dp, LI_xydy3)
!         CALL line_integral_xydy(  [xomin,yomax], [xomin,yomin], 1E-9_dp, LI_xydy4)

!         w1y = 1._dp / Ad * (LI_xydy1  + LI_xydy2  + LI_xydy3  + LI_xydy4 ) - w0 * y_src( j_src)

!         d_dst( j,i) = d_dst( j,i) + w0  * d_src(   j_src,i_src) + &
!                                     w1x * ddx_src( j_src,i_src) + &
!                                     w1y * ddy_src( j_src,i_src)

!       END DO ! DO j_src = jr_src( j,1), jr_src( j,2)
!       END DO ! DO i_src = ir_src( i,1), ir_src( i,2)

!       IF (Asum < Ad) mask_dst_outside_src( j,i) = 1

!     END DO ! DO j = 1, ny_dst
!     END DO ! DO i = i1, i2
!     CALL sync

!     ! Use nearest-neighbour extrapolation for dst cells outside of the src grid
!     ! =========================================================================

!     ! Find the range of grid cells that were mapped correctly
!     igmin = 0
!     igmax = 0
!     jgmin = 0
!     jgmax = 0

!     j = INT( REAL(ny_dst,dp)/2._dp)
!     DO i = 1, nx_dst
!       IF (mask_dst_outside_src( j,i) == 0) THEN
!         igmin = i
!         EXIT
!       END IF
!     END DO
!     DO i = nx_dst, 1, -1
!       IF (mask_dst_outside_src( j,i) == 0) THEN
!         igmax = i
!         EXIT
!       END IF
!     END DO

!     i = INT( REAL(nx_dst,dp)/2._dp)
!     DO j = 1, ny_dst
!       IF (mask_dst_outside_src( j,i) == 0) THEN
!         jgmin = j
!         EXIT
!       END IF
!     END DO
!     DO j = ny_dst, 1, -1
!       IF (mask_dst_outside_src( j,i) == 0) THEN
!         jgmax = j
!         EXIT
!       END IF
!     END DO

!     ! Corners
!     IF (par%master) THEN
!       ! Southwest
!       d_dst( 1      :jgmin-1 ,1      :igmin-1) = d_dst( jgmin,igmin)
!       ! Southeast
!       d_dst( 1      :jgmin-1 ,igmax+1:nx_dst ) = d_dst( jgmin,igmax)
!       ! Northwest
!       d_dst( jgmax+1:ny_dst  ,1      :igmin-1) = d_dst( jgmax,igmin)
!       ! Northeast
!       d_dst( jgmax+1:ny_dst  ,igmax+1:nx_dst ) = d_dst( jgmax,igmax)
!     END IF ! IF (par%master) THEN
!     CALL sync

!     ! Borders
!     DO i = MAX(i1,igmin), MIN(i2,igmax)
!       ! South
!       d_dst( 1      :jgmin-1,i) = d_dst( jgmin,i)
!       ! North
!       d_dst( jgmax+1:ny_dst ,i) = d_dst( jgmax,i)
!     END DO
!     DO j = MAX(j1,jgmin), MIN(j2,jgmax)
!       ! West
!       d_dst( j,1      :igmin-1) = d_dst( j,igmin)
!       ! East
!       d_dst( j,igmax+1:nx_dst ) = d_dst( j,igmax)
!     END DO
!     CALL sync

!     ! Clean up after yourself
!     CALL deallocate_shared( wddx_src             )
!     CALL deallocate_shared( wddy_src             )
!     CALL deallocate_shared( wmask_dst_outside_src)

!   END SUBROUTINE map_square_to_square_cons_2nd_order_2D
!   SUBROUTINE map_square_to_square_cons_2nd_order_3D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src, d_dst)
!     ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)

!     IMPLICIT NONE

!     ! Input and output variables
!     INTEGER,                            INTENT(IN)    :: nx_src
!     INTEGER,                            INTENT(IN)    :: ny_src
!     REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_src
!     REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_src
!     INTEGER,                            INTENT(IN)    :: nx_dst
!     INTEGER,                            INTENT(IN)    :: ny_dst
!     REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_dst
!     REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_dst
!     REAL(dp), DIMENSION(:,:,:),         INTENT(IN)    :: d_src
!     REAL(dp), DIMENSION(:,:,:),         INTENT(OUT)   :: d_dst

!     ! Local variables:
!     CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_square_to_square_cons_2nd_order_3D'
!     INTEGER                                           :: nz, k, i1_src, i2_src, i1_dst, i2_dst
!     REAL(dp), DIMENSION(:,:  ), POINTER               ::  d_src_2D,  d_dst_2D
!     INTEGER                                           :: wd_src_2D, wd_dst_2D

!     nz = SIZE( d_src,1)

!     ! Safety
!     IF (SIZE(d_dst,1) /= nz) THEN
!       IF (par%master) WRITE(0,*) 'map_square_to_square_cons_2nd_order_3D - ERROR: nz_src /= nz_dst!'
!       CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!     END IF

!     CALL partition_list( nx_src, par%i, par%n, i1_src, i2_src)
!     CALL partition_list( nx_dst, par%i, par%n, i1_dst, i2_dst)

!     ! Allocate shared memory
!     CALL allocate_shared_dp_2D( ny_src, nx_src, d_src_2D, wd_src_2D)
!     CALL allocate_shared_dp_2D( ny_dst, nx_dst, d_dst_2D, wd_dst_2D)

!     ! Remap the fields one layer at a time
!     DO k = 1, nz
!       d_src_2D(   :,i1_src:i2_src) = d_src(    k,:,i1_src:i2_src)
!       CALL map_square_to_square_cons_2nd_order_2D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src_2D, d_dst_2D)
!       d_dst(    k,:,i1_dst:i2_dst) = d_dst_2D(   :,i1_dst:i2_dst)
!     END DO

!     ! Clean up after yourself
!     CALL deallocate_shared( wd_src_2D)
!     CALL deallocate_shared( wd_dst_2D)

!   END SUBROUTINE map_square_to_square_cons_2nd_order_3D

! == Line integrals used in conservative remapping
  SUBROUTINE line_integral_xdy(   p, q, tol_dist, I_pq)
    ! Calculate the line integral x dy from p to q

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
    REAL(dp),                                INTENT(IN)    :: tol_dist
    REAL(dp),                                INTENT(OUT)   :: I_pq

    ! Local variables:
    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy

    xp = p(1)
    yp = p(2)
    xq = q(1)
    yq = q(2)

    IF (ABS(yp-yq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF

    dx = q(1)-p(1)
    dy = q(2)-p(2)

    I_pq = xp*dy - yp*dx + (dx / (2._dp*dy)) * (yq**2 - yp**2)

  END SUBROUTINE line_integral_xdy
  SUBROUTINE line_integral_mxydx( p, q, tol_dist, I_pq)
    ! Calculate the line integral -xy dx from p to q

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
    REAL(dp),                                INTENT(IN)    :: tol_dist
    REAL(dp),                                INTENT(OUT)   :: I_pq

    ! Local variables:
    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy

    xp = p(1)
    yp = p(2)
    xq = q(1)
    yq = q(2)

    IF (ABS(xp-xq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF

    dx = q(1)-p(1)
    dy = q(2)-p(2)

    I_pq = (1._dp/2._dp * (xp*dy/dx - yp) * (xq**2-xp**2)) - (1._dp/3._dp * dy/dx * (xq**3-xp**3))

  END SUBROUTINE line_integral_mxydx
  SUBROUTINE line_integral_xydy(  p, q, tol_dist, I_pq)
    ! Calculate the line integral xy dy from p to q

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(2),                  INTENT(IN)    :: p, q
    REAL(dp),                                INTENT(IN)    :: tol_dist
    REAL(dp),                                INTENT(OUT)   :: I_pq

    ! Local variables:
    REAL(dp)                                               :: xp, yp, xq, yq, dx, dy

    xp = p(1)
    yp = p(2)
    xq = q(1)
    yq = q(2)

    IF (ABS(yp-yq) < tol_dist) THEN
      I_pq = 0._dp
      RETURN
    END IF

    dx = q(1)-p(1)
    dy = q(2)-p(2)

    I_pq = (1._dp/2._dp * (xp - yp*dx/dy) * (yq**2-yp**2)) + (1._dp/3._dp * dx/dy * (yq**3-yp**3))

  END SUBROUTINE line_integral_xydy

! == Smoothing operations
  SUBROUTINE smooth_Gaussian_2D_grid( grid, d, r)
    ! Apply a Gaussian smoothing filter of with sigma = n*dx to the 2D data field d

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in m

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Gaussian_2D_grid'
    INTEGER                                            :: i,j,ii,jj,n
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_ext,  d_ext_smooth
    INTEGER                                            :: wd_ext, wd_ext_smooth
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: f

    n = CEILING( r / grid%dx) * 3  ! Number of cells to extend the data by (3 standard deviations is enough to capture the tails of the normal distribution)

    ! Fill in the smoothing filters
    ALLOCATE( f( -n:n))
    f = 0._dp
    DO i = -n, n
      f(i) = EXP( -0.5_dp * (REAL(i,dp) * grid%dx/r)**2)
    END DO
    f = f / SUM(f)

    ! Allocate temporary shared memory for the extended and smoothed data fields
    allocate(d_ext(        grid%nx + 2*n, grid%ny + 2*n ))
    allocate(d_ext_smooth( grid%nx + 2*n, grid%ny + 2*n ))

    ! Copy data to the extended array and fill in the margins
    d(1+n:grid%nx+n,1+n:grid%ny+n) = d
    ! West
    d_ext( n+1:n+grid%nx, 1            ) = d( :      ,1      )
    ! East
    d_ext( n+1:n+grid%nx, grid%ny+2*n  ) = d( :      ,grid%ny)
    ! South
    d_ext( 1            , n+1:n+grid%ny) = d( 1      ,:      )
    ! North
    d_ext( grid%nx+2*n  , n+1:n+grid%ny) = d( grid%nx,:      )
    ! Corners
    d_ext( 1:n,                     1:n                    ) = d( 1      ,1      )
    d_ext( 1:n,                     grid%ny+n+1:grid%ny+2*n) = d( 1      ,grid%ny)
    d_ext( grid%nx+n+1:grid%nx+2*n, 1:n                    ) = d( grid%nx,1      )
    d_ext( grid%nx+n+1:grid%nx+2*n, grid%ny+n+1:grid%ny+2*n) = d( grid%nx,grid%ny)

    ! Convolute extended data with the smoothing filter
    d_ext_smooth( 1+n:grid%nx+n,:) = 0._dp

    DO j = 1, grid%ny
    DO i = 1, grid%nx
      DO ii = -n, n
        d_ext_smooth( i+n,j+n) = d_ext_smooth( i+n,j+n) + d_ext( i+n+ii,j+n) * f(ii)
      END DO
    END DO
    END DO

    d_ext( 1+n:grid%nx+n,:) = d_ext_smooth( 1+n:grid%nx+n,:)

    DO j = 1, grid%ny
      d_ext(           1:          n,j) = d( 1      ,j)
      d_ext( grid%nx+n+1:grid%nx+2*n,j) = d( grid%nx,j)
    END DO

    d_ext_smooth( 1+n:grid%nx+n,:) = 0._dp

    DO j = 1, grid%ny
    DO i = 1, grid%nx
      DO jj = -n, n
        d_ext_smooth( i+n,j+n) = d_ext_smooth( i+n,j+n) + d_ext( i+n,j+n+jj) * f(jj)
      END DO
    END DO
    END DO

    ! Copy data back
    d = d_ext_smooth( 1+n:grid%nx+n, 1+n:grid%ny+n)

    ! Clean up after yourself
    DEALLOCATE( f)
    deallocate( d_ext)
    deallocate( d_ext_smooth)

  END SUBROUTINE smooth_Gaussian_2D_grid
  SUBROUTINE smooth_Gaussian_3D_grid( grid, d, r)
    ! Apply a Gaussian smoothing filter of with sigma = n*dx to the 3D data field d

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in km

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Gaussian_3D_grid'
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_2D

    ! Allocate temporary shared memory for the extended and smoothed data fields
    allocate(d_2D( grid%nx, grid%ny))

    DO k = 1, SIZE( d,3)
      d_2D = d( :,:,k)
      CALL smooth_Gaussian_2D_grid( grid, d_2D, r)
      d( :,:,k) = d_2D
    END DO

    ! Clean up after yourself
    deallocate( d_2D)

  END SUBROUTINE smooth_Gaussian_3D_grid
!   SUBROUTINE smooth_Shepard_2D_grid( grid, d, r)
!     ! Apply a Shepard smoothing filter of with sigma = n*dx to the 2D data field d

!     IMPLICIT NONE

!     ! In/output variables:
!     TYPE(type_grid),                     INTENT(IN)    :: grid
!     REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
!     REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in m

!     ! Local variables:
!     CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Shepard_2D_grid'
!     INTEGER                                            :: i,j,k,l,n
!     REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_ext,  d_ext_smooth
!     INTEGER                                            :: wd_ext, wd_ext_smooth
!     REAL(dp)                                           :: ShepNumSum     ! The sum of the numerators in the Shepard weighting
!     REAL(dp)                                           :: ShepDenSum     ! The sum of the denumerators in the Shepard weighting
!     REAL(dp)                                           :: distance       ! in gridsize units
!     REAL(dp)                                           :: smooth_radius  ! in gridsize units
!     REAL(dp)                                           :: exponent       ! The distance weighting exponent in the Shepard weighting

!     n = CEILING( r / grid%dx) ! Number of cells to extend the data by (3 standard deviations is enough to capture the tails of the normal distribution)

!     smooth_radius = REAL(n,dp)
!     exponent      = 2._dp

!     ! Allocate temporary shared memory for the extended and smoothed data fields
!     CALL allocate_shared_dp_2D( grid%ny + 2*n, grid%nx + 2*n, d_ext,        wd_ext       )
!     CALL allocate_shared_dp_2D( grid%ny + 2*n, grid%nx + 2*n, d_ext_smooth, wd_ext_smooth)

!     ! Copy data to the extended array and fill in the margins
!     DO i = grid%i1, grid%i2
!     DO j = 1, grid%ny
!       d_ext( j+n,i+n) = d( j,i)
!     END DO
!     END DO
!     IF (par%master) THEN
!       ! West
!       d_ext( n+1:n+grid%ny, 1            ) = d( :      ,1      )
!       ! East
!       d_ext( n+1:n+grid%ny, grid%nx+2*n  ) = d( :      ,grid%nx)
!       ! South
!       d_ext( 1            , n+1:n+grid%nx) = d( 1      ,:      )
!       ! North
!       d_ext( grid%ny+2*n  , n+1:n+grid%nx) = d( grid%ny,:      )
!       ! Corners
!       d_ext( 1:n,                     1:n                    ) = d( 1      ,1      )
!       d_ext( 1:n,                     grid%nx+n+1:grid%nx+2*n) = d( 1      ,grid%nx)
!       d_ext( grid%ny+n+1:grid%ny+2*n, 1:n                    ) = d( grid%ny,1      )
!       d_ext( grid%ny+n+1:grid%ny+2*n, grid%nx+n+1:grid%nx+2*n) = d( grid%ny,grid%nx)
!     END IF
!     CALL sync

!     ! Convolute extended data with the smoothing filter
!     d_ext_smooth( :,grid%i1+n:grid%i2+n) = 0._dp
!     CALL sync

!     DO i = grid%i1, grid%i2
!     DO j = 1,       grid%ny
!       ShepNumSum = 0._dp
!       ShepDenSum = 0._dp
!       DO k = -n, n
!       DO l = -n, n
!         distance = SQRT(REAL(k,dp)**2 + REAL(l,dp)**2)
!         IF (distance <= smooth_radius .AND. distance > 0._dp) THEN
!           ShepNumSum = ShepNumSum + (d_ext( j+n+l,i+n+k)/(distance**exponent))
!           ShepDenSum = ShepDenSum + (1._dp/(distance**exponent))
!         END IF
!       END DO
!       END DO
!       d_ext_smooth( j+n,i+n) = ShepNumSum/ShepDenSum
!     END DO
!     END DO
!     CALL sync

!     ! Copy data back
!     DO i = grid%i1, grid%i2
!     DO j = 1, grid%ny
!       d( j,i) = d_ext_smooth( j+n, i+n)
!     END DO
!     END DO

!     ! Clean up after yourself
!     CALL deallocate_shared( wd_ext)
!     CALL deallocate_shared( wd_ext_smooth)

!   END SUBROUTINE smooth_Shepard_2D_grid
!   SUBROUTINE smooth_Shepard_3D_grid( grid, d, r)
!     ! Apply a Shepard smoothing filter of with sigma = n*dx to the 3D data field d

!     IMPLICIT NONE

!     ! In/output variables:
!     TYPE(type_grid),                     INTENT(IN)    :: grid
!     REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
!     REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in km

!     ! Local variables:
!     CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Shepard_3D_grid'
!     INTEGER                                            :: k
!     REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_2D
!     INTEGER                                            :: wd_2D

!     ! Allocate temporary shared memory for the extended and smoothed data fields
!     CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_2D, wd_2D)

!     DO k = 1, SIZE( d,2)
!       d_2D( :,grid%i1:grid%i2) = d( k,:,grid%i1:grid%i2)
!       CALL smooth_Shepard_2D_grid( grid, d_2D, r)
!       d( k,:,grid%i1:grid%i2) = d_2D( :,grid%i1:grid%i2)
!     END DO

!     ! Clean up after yourself
!     CALL deallocate_shared( wd_2D)

!   END SUBROUTINE smooth_Shepard_3D_grid

! == Remove Lake Vostok from Antarctic input geometry data
  SUBROUTINE remove_Lake_Vostok( x, y, Hi, Hb, Hs)
    ! Remove Lake Vostok from Antarctic input geometry data
    ! by manually increasing ice thickness so that Hi = Hs - Hb
    !
    ! NOTE: since UFEMISM doesn't consider subglacial lakes, Vostok simply shows
    !       up as a "dip" in the initial geometry. The model will run fine, the dip
    !       fills up in a few centuries, but it slows down the model for a while and
    !       it looks ugly, so we just remove it right away.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: x,y
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: Hi
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Hb
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: Hs

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'remove_Lake_Vostok'
    INTEGER                                       :: i,j,nx,ny
    REAL(dp), PARAMETER                           :: lake_Vostok_xmin = 1164250.0
    REAL(dp), PARAMETER                           :: lake_Vostok_xmax = 1514250.0
    REAL(dp), PARAMETER                           :: lake_Vostok_ymin = -470750.0
    REAL(dp), PARAMETER                           :: lake_Vostok_ymax = -220750.0
    INTEGER                                       :: il,iu,jl,ju


    nx = SIZE( Hi,1)
    ny = SIZE( Hi,2)

    il = 1
    DO WHILE (x( il) < lake_Vostok_xmin)
      il = il+1
    END DO
    iu = nx
    DO WHILE (x( iu) > lake_Vostok_xmax)
      iu = iu-1
    END DO
    jl = 1
    DO WHILE (y( jl) < lake_Vostok_ymin)
      jl = jl+1
    END DO
    ju = ny
    DO WHILE (y( ju) > lake_Vostok_ymax)
      ju = ju-1
    END DO

    DO i = il, iu
    DO j = jl, ju
      Hi( i,j) = Hs( i,j) - Hb( i,j)
    END DO
    END DO

  END SUBROUTINE remove_Lake_Vostok

! == Analytical solution by Schoof 2006 for the "SSA_icestream" benchmark experiment
  SUBROUTINE SSA_Schoof2006_analytical_solution( tantheta, h0, A_flow, y, U, tauc)

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: tantheta   ! Surface slope in the x-direction
    REAL(dp),                            INTENT(IN)    :: h0         ! Ice thickness
    REAL(dp),                            INTENT(IN)    :: A_flow     ! Ice flow factor
    REAL(dp),                            INTENT(IN)    :: y          ! y-coordinate
    REAL(dp),                            INTENT(OUT)   :: U          ! Ice velocity in the x-direction
    REAL(dp),                            INTENT(OUT)   :: tauc       ! Till yield stress

    ! Local variables:
    REAL(dp)                                           :: m, B, f, W, ua, ub, uc, ud, ue
    REAL(dp)                                           :: L = 40000._dp     ! Ice-stream width (m)

    m = C%SSA_icestream_m

    ! Calculate the gravitational driving stress f
    f = ice_density * grav * h0 * tantheta

    ! Calculate the ice hardness factor B
    B = A_flow**(-1._dp/C%n_flow)

    ! Calculate the "ice stream half-width" W
    W = L * (m+1._dp)**(1._dp/m)

    ! Calculate the till yield stress across the stream
    tauc = f * ABS(y/L)**m

    ! Calculate the analytical solution for u
    ua = -2._dp * f**3 * L**4 / (B**3 * h0**3)
    ub = ( 1._dp / 4._dp                           ) * (   (y/L)**     4._dp  - (m+1._dp)**(       4._dp/m) )
    uc = (-3._dp / ((m+1._dp)    * (      m+4._dp))) * (ABS(y/L)**(  m+4._dp) - (m+1._dp)**(1._dp+(4._dp/m)))
    ud = ( 3._dp / ((m+1._dp)**2 * (2._dp*m+4._dp))) * (ABS(y/L)**(2*m+4._dp) - (m+1._dp)**(2._dp+(4._dp/m)))
    ue = (-1._dp / ((m+1._dp)**3 * (3._dp*m+4._dp))) * (ABS(y/L)**(3*m+4._dp) - (m+1._dp)**(3._dp+(4._dp/m)))
    u = ua * (ub + uc + ud + ue)

    ! Outside the ice-stream, velocity is zero
    IF (ABS(y) > w) U = 0._dp

  END SUBROUTINE SSA_Schoof2006_analytical_solution

! == Some wrappers for LAPACK matrix functionality
  function tridiagonal_solve( ldiag, diag, udiag, rhs) result(x)
    ! Lapack tridiagnal solver (in double precision):
    ! Matrix system solver for tridiagonal matrices.
    ! Used e.g. in solving the ADI scheme.
    ! ldiag = lower diagonal elements (j,j-1) of the matrix
    ! diag  = diagonal elements (j,j) of the matrix
    ! udiag = upper diagonal elements (j,j+1) of the matrix
    ! rhs   = right hand side of the matrix equation in the ADI scheme

    implicit none

    ! Input variables:
    real(dp), dimension(:),            intent(in) :: diag
    real(dp), dimension(size(diag)-1), intent(in) :: udiag, ldiag
    real(dp), dimension(size(diag)),   intent(in) :: rhs

    ! Result variable:
    real(dp), dimension(size(diag))               :: x

    ! Local variables:
    integer                                       :: info
    real(dp), dimension(size(diag))               :: diag_copy
    real(dp), dimension(size(udiag))              :: udiag_copy, ldiag_copy

    ! The LAPACK solver will overwrite the rhs with the solution x. Therefore we
    ! first copy the rhs in the solution vector x:
    x = rhs

    ! The LAPACK solver will change the elements in the matrix, therefore we copy them:
    diag_copy  =  diag
    udiag_copy = udiag
    ldiag_copy = ldiag

    call DGTSV(size(diag), 1, ldiag_copy, diag_copy, udiag_copy, x, size(diag), info)

    if (info /= 0) then
      write(0,*) 'tridiagonal_solve - ERROR: LAPACK solver DGTSV returned error message info = ', info
      call MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    end if

  end function tridiagonal_solve

! == Finding inverses of some very small matrices
  SUBROUTINE calc_matrix_inverse_2_by_2( A, Ainv)
    ! Direct inversion of a 2-by-2 matrix

    IMPLICIT NONE

    ! In- and output variables:
    REAL(dp), DIMENSION(2,2  ),          INTENT(IN)    :: A
    REAL(dp), DIMENSION(2,2  ),          INTENT(INOUT) :: Ainv

    ! Local variables
    REAL(dp)                                           :: detA

    ! Calculate the determinant of A
    detA = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

    ! Safety
    IF (detA == 0.0d0 ) THEN
      PRINT *, A(1,1), ',', A(1,2)
      PRINT *, A(2,1), ',', A(2,2)
      write(0,*) 'determinant:', detA 
      WRITE(0,*) 'calc_matrix_inverse_2_by_2 - ERROR: matrix is numerically singular!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Find the inverse of A
    Ainv( 1,1) =  A( 2,2) / detA
    Ainv( 1,2) = -A( 1,2) / detA
    Ainv( 2,1) = -A( 2,1) / detA
    Ainv( 2,2) =  A( 1,1) / detA

  END SUBROUTINE calc_matrix_inverse_2_by_2
  SUBROUTINE calc_matrix_inverse_3_by_3( A, Ainv)
    ! Direct inversion of a 3-by-3 matrix
    !
    ! See: https://metric.ma.ic.ac.uk/metric_public/matrices/inverses/inverses2.html

    IMPLICIT NONE

    ! In- and output variables:
    REAL(dp), DIMENSION(3,3  ),          INTENT(IN)    :: A
    REAL(dp), DIMENSION(3,3  ),          INTENT(INOUT) :: Ainv

    ! Local variables
    REAL(dp)                                           :: detA

    ! Calculate the minors of A
    Ainv( 1,1) = A( 2,2) * A( 3,3) - A( 2,3) * A( 3,2)
    Ainv( 1,2) = A( 2,1) * A( 3,3) - A( 2,3) * A( 3,1)
    Ainv( 1,3) = A( 2,1) * A( 3,2) - A( 2,2) * A( 3,1)
    Ainv( 2,1) = A( 1,2) * A( 3,3) - A( 1,3) * A( 3,2)
    Ainv( 2,2) = A( 1,1) * A( 3,3) - A( 1,3) * A( 3,1)
    Ainv( 2,3) = A( 1,1) * A( 3,2) - A( 1,2) * A( 3,1)
    Ainv( 3,1) = A( 1,2) * A( 2,3) - A( 1,3) * A( 2,2)
    Ainv( 3,2) = A( 1,1) * A( 2,3) - A( 1,3) * A( 2,1)
    Ainv( 3,3) = A( 1,1) * A( 2,2) - A( 1,2) * A( 2,1)

    ! Calculate the determinant of A
    detA = A( 1,1) * Ainv( 1,1) - A( 1,2) * Ainv( 1,2) + A( 1,3) * Ainv( 1,3)

    ! Safety
    IF (detA == 0.0d0 ) THEN
      PRINT *, A(1,1), ',', A(1,2), ',', A(1,3)
      PRINT *, A(2,1), ',', A(2,2), ',', A(2,3)
      PRINT *, A(3,1), ',', A(3,2), ',', A(3,3)
      write(0,*) 'determinant:', detA 
      WRITE(0,*) 'calc_matrix_inverse_3_by_3 - ERROR: matrix is numerically singular!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF

    ! Change matrix of minors to get the matrix of cofactors
    Ainv( 1,2) = -Ainv( 1,2)
    Ainv( 2,1) = -Ainv( 2,1)
    Ainv( 2,3) = -Ainv( 2,3)
    Ainv( 3,2) = -Ainv( 3,2)

    ! Transpose matrix of cofactors
    Ainv( 1,2) = Ainv( 1,2) + Ainv( 2,1)
    Ainv( 2,1) = Ainv( 1,2) - Ainv( 2,1)
    Ainv( 1,2) = Ainv( 1,2) - Ainv( 2,1)

    Ainv( 1,3) = Ainv( 1,3) + Ainv( 3,1)
    Ainv( 3,1) = Ainv( 1,3) - Ainv( 3,1)
    Ainv( 1,3) = Ainv( 1,3) - Ainv( 3,1)

    Ainv( 2,3) = Ainv( 2,3) + Ainv( 3,2)
    Ainv( 3,2) = Ainv( 2,3) - Ainv( 3,2)
    Ainv( 2,3) = Ainv( 2,3) - Ainv( 3,2)

    ! Divide by det(A)
    Ainv = Ainv / detA

  END SUBROUTINE calc_matrix_inverse_3_by_3
  SUBROUTINE calc_matrix_inverse_5_by_5( A, Ainv)
    ! Calculate the inverse Ainv of an n-by-n matrix A using LAPACK

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(5,5  )            , INTENT(IN)    :: A
    REAL(dp), DIMENSION(5,5  )            , INTENT(OUT)   :: Ainv

    ! Local variables:
    integer                                            :: n
    logical                                            :: ok_flag
    real(dp)                                           :: detA

    call m55inv(A,Ainv,detA,ok_flag)

    ! Safety
    IF (.not. ok_flag) THEN
      do n = 1, 5
        write(0,*) A(1,n), ',', A(2,n), ',', A(3,n), ',', A(4,n), ',', A(5,n)
      end do
      write(0,*) 'determinant:', detA 
      write(0,*) 'calc_matrix_inverse_5_by_5 - error: matrix inversion failed (singular matrix)!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
  contains
    include 'm55inv.f90'

  END SUBROUTINE calc_matrix_inverse_5_by_5

! == Debugging
  SUBROUTINE check_for_NaN_dp_1D(  d, d_name)
    ! Check if NaN values occur in the 1-D dp data field d
    ! NOTE: parallelised!

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name

    ! Local variables:
    INTEGER                                                :: nx,i,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc

    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN

    ! Field size
    nx = SIZE(d,1)

    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = i1, i2

      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...

      IF     (d( i) /= d( i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
      ELSEIF (d( i) > HUGE( d( i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
      ELSEIF (d( i)  < -HUGE( d( i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]'), int_01 = i)
      END IF

    END DO
    CALL sync

  END SUBROUTINE check_for_NaN_dp_1D
  SUBROUTINE check_for_NaN_dp_2D(  d, d_name)
    ! Check if NaN values occur in the 2-D dp data field d
    ! NOTE: parallelised!

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name

    ! Local variables:
    INTEGER                                                :: nx,ny,i,j,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc

    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN

    ! Field size
    nx = SIZE(d,2)
    ny = SIZE(d,1)

    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny

      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...

      IF     (d( j,i) /= d( j,i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      ELSEIF (d( j,i) > HUGE( d( j,i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      ELSEIF (d( j,i) < -HUGE( d( j,i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      END IF

    END DO
    END DO
    CALL sync

  END SUBROUTINE check_for_NaN_dp_2D
  SUBROUTINE check_for_NaN_dp_3D(  d, d_name)
    ! Check if NaN values occur in the 3-D dp data field d
    ! NOTE: parallelised!

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name

    ! Local variables:
    INTEGER                                                :: nx,ny,nz,i,j,k,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc

    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN

    ! Field size
    nx = SIZE(d,3)
    ny = SIZE(d,2)
    nz = SIZE(d,1)

    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    DO k = 1, nz

      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...

      IF     (d( k,j,i) /= d( k,j,i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (d( k,j,i) > HUGE( d( k,j,i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (d( k,j,i) < -HUGE( d( k,j,i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      END IF

    END DO
    END DO
    END DO
    CALL sync

  END SUBROUTINE check_for_NaN_dp_3D
  SUBROUTINE check_for_NaN_int_1D( d, d_name)
    ! Check if NaN values occur in the 1-D int data field d
    ! NOTE: parallelised!

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name

    ! Local variables:
    INTEGER                                                :: nx,i,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc

    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN

    ! Field size
    nx = SIZE(d,1)

    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = i1, i2

      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...

      IF     (d( i) /= d( i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
      ELSEIF (d( i) > HUGE( d( i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]' ), int_01 = i)
      ELSEIF (d( i) < -HUGE( d( i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01}]'), int_01 = i)
      END IF

    END DO
    CALL sync

  END SUBROUTINE check_for_NaN_int_1D
  SUBROUTINE check_for_NaN_int_2D( d, d_name)
    ! Check if NaN values occur in the 2-D int data field d
    ! NOTE: parallelised!

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name

    ! Local variables:
    INTEGER                                                :: nx,ny,i,j,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc

    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN

    ! Field size
    nx = SIZE(d,2)
    ny = SIZE(d,1)

    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny

      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...

      IF     (d( j,i) /= d( j,i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      ELSEIF (d( j,i) > HUGE( d( j,i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      ELSEIF (d( j,i) < -HUGE( d( j,i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02}]'), int_01 = j, int_02 = i)
      END IF

    END DO
    END DO
    CALL sync

  END SUBROUTINE check_for_NaN_int_2D
  SUBROUTINE check_for_NaN_int_3D( d, d_name)
    ! Check if NaN values occur in the 3-D int data field d
    ! NOTE: parallelised!

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,  DIMENSION(:,:,:),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name

    ! Local variables:
    INTEGER                                                :: nx,ny,nz,i,j,k,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc

    ! Only do this when so specified in the config
    IF (.NOT. C%do_check_for_NaN) RETURN

    ! Field size
    nx = SIZE(d,3)
    ny = SIZE(d,2)
    nz = SIZE(d,1)

    ! Parallelisation range
    CALL partition_list( nx, par%i, par%n, i1, i2)

    ! Variable name and routine name
    IF (PRESENT( d_name)) THEN
      d_name_loc = TRIM(d_name)
    ELSE
      d_name_loc = '?'
    END IF

    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    DO k = 1, nz

      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...

      IF     (d( k,j,i) /= d( k,j,i)) THEN
        CALL crash(  ('detected NaN in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (d( k,j,i) > HUGE( d( k,j,i))) THEN
        CALL crash(  ('detected Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      ELSEIF (d( k,j,i) < -HUGE( d( k,j,i))) THEN
        CALL crash( ('detected -Inf in variable "' // TRIM( d_name_loc) // '" at [{int_01},{int_02},{int_03}]'), int_01 = k, int_02 = j, int_03 = i)
      END IF

    END DO
    END DO
    END DO
    CALL sync

  END SUBROUTINE check_for_NaN_int_3D

!   ! == Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)
!   SUBROUTINE transpose_dp_2D( d, wd)
!     ! Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)

!     IMPLICIT NONE

!     ! In/output variables:
!     REAL(dp), DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: d
!     INTEGER,                             INTENT(INOUT) :: wd

!     ! Local variables:
!     INTEGER                                      :: i,j,nx,ny,i1,i2
!     REAL(dp), DIMENSION(:,:  ), POINTER          ::  d_temp
!     INTEGER                                      :: wd_temp

!     nx = SIZE( d,1)
!     ny = SIZE( d,2)
!     CALL partition_list( nx, par%i, par%n, i1, i2)

!     ! Allocate temporary memory
!     CALL allocate_shared_dp_2D( nx, ny, d_temp, wd_temp)

!     ! Copy data to temporary memory
!     DO i = i1,i2
!     DO j = 1, ny
!       d_temp( i,j) = d( i,j)
!     END DO
!     END DO
!     CALL sync

!     ! Deallocate memory
!     CALL deallocate_shared( wd)

!     ! Reallocate transposed memory
!     CALL allocate_shared_dp_2D( ny, nx, d, wd)

!     ! Copy and transpose data from temporary memory
!     DO i = i1, i2
!     DO j = 1, ny
!       d( j,i) = d_temp( i,j)
!     END DO
!     END DO
!     CALL sync

!     ! Deallocate temporary memory
!     CALL deallocate_shared( wd_temp)

!   END SUBROUTINE transpose_dp_2D
!   SUBROUTINE transpose_dp_3D( d, wd)
!     ! Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)

!     IMPLICIT NONE

!     ! In/output variables:
!     REAL(dp), DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: d
!     INTEGER,                             INTENT(INOUT) :: wd

!     ! Local variables:
!     INTEGER                                      :: i,j,k,nx,ny,nz,i1,i2
!     REAL(dp), DIMENSION(:,:,:), POINTER          ::  d_temp
!     INTEGER                                      :: wd_temp

!     nx = SIZE( d,1)
!     ny = SIZE( d,2)
!     nz = SIZE( d,3)
!     CALL partition_list( nx, par%i, par%n, i1, i2)

!     ! Allocate temporary memory
!     CALL allocate_shared_dp_3D( nx, ny, nz, d_temp, wd_temp)

!     ! Copy data to temporary memory
!     DO i = i1,i2
!     DO j = 1, ny
!     DO k = 1, nz
!       d_temp( i,j,k) = d( i,j,k)
!     END DO
!     END DO
!     END DO
!     CALL sync

!     ! Deallocate memory
!     CALL deallocate_shared( wd)

!     ! Reallocate transposed memory
!     CALL allocate_shared_dp_3D( nz, ny, nx, d, wd)

!     ! Copy and transpose data from temporary memory
!     DO i = i1, i2
!     DO j = 1, ny
!     DO k = 1, nz
!       d( k,j,i) = d_temp( i,j,k)
!     END DO
!     END DO
!     END DO
!     CALL sync

!     ! Deallocate temporary memory
!     CALL deallocate_shared( wd_temp)

!   END SUBROUTINE transpose_dp_3D
!   SUBROUTINE transpose_int_2D( d, wd)
!     ! Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)

!     IMPLICIT NONE

!     ! In/output variables:
!     INTEGER,  DIMENSION(:,:  ), POINTER, INTENT(INOUT) :: d
!     INTEGER,                             INTENT(INOUT) :: wd

!     ! Local variables:
!     INTEGER                                      :: i,j,nx,ny,i1,i2
!     INTEGER,  DIMENSION(:,:  ), POINTER          ::  d_temp
!     INTEGER                                      :: wd_temp

!     nx = SIZE( d,1)
!     ny = SIZE( d,2)
!     CALL partition_list( nx, par%i, par%n, i1, i2)

!     ! Allocate temporary memory
!     CALL allocate_shared_int_2D( nx, ny, d_temp, wd_temp)

!     ! Copy data to temporary memory
!     DO i = i1,i2
!     DO j = 1, ny
!       d_temp( i,j) = d( i,j)
!     END DO
!     END DO
!     CALL sync

!     ! Deallocate memory
!     CALL deallocate_shared( wd)

!     ! Reallocate transposed memory
!     CALL allocate_shared_int_2D( ny, nx, d, wd)

!     ! Copy and transpose data from temporary memory
!     DO i = i1, i2
!     DO j = 1, ny
!       d( j,i) = d_temp( i,j)
!     END DO
!     END DO
!     CALL sync

!     ! Deallocate temporary memory
!     CALL deallocate_shared( wd_temp)

!   END SUBROUTINE transpose_int_2D
!   SUBROUTINE transpose_int_3D( d, wd)
!     ! Transpose a data field (i.e. go from [i,j] to [j,i] indexing or the other way round)

!     IMPLICIT NONE

!     ! In/output variables:
!     INTEGER,  DIMENSION(:,:,:), POINTER, INTENT(INOUT) :: d
!     INTEGER,                             INTENT(INOUT) :: wd

!     ! Local variables:
!     INTEGER                                      :: i,j,k,nx,ny,nz,i1,i2
!     INTEGER,  DIMENSION(:,:,:), POINTER          ::  d_temp
!     INTEGER                                      :: wd_temp

!     nx = SIZE( d,1)
!     ny = SIZE( d,2)
!     nz = SIZE( d,3)
!     CALL partition_list( nx, par%i, par%n, i1, i2)

!     ! Allocate temporary memory
!     CALL allocate_shared_int_3D( nx, ny, nz, d_temp, wd_temp)

!     ! Copy data to temporary memory
!     DO i = i1,i2
!     DO j = 1, ny
!     DO k = 1, nz
!       d_temp( i,j,k) = d( i,j,k)
!     END DO
!     END DO
!     END DO
!     CALL sync

!     ! Deallocate memory
!     CALL deallocate_shared( wd)

!     ! Reallocate transposed memory
!     CALL allocate_shared_int_3D( nz, ny, nx, d, wd)

!     ! Copy and transpose data from temporary memory
!     DO i = i1, i2
!     DO j = 1, ny
!     DO k = 1, nz
!       d( k,j,i) = d_temp( i,j,k)
!     END DO
!     END DO
!     END DO
!     CALL sync

!     ! Deallocate temporary memory
!     CALL deallocate_shared( wd_temp)

!   END SUBROUTINE transpose_int_3D

  ! == 2nd-order conservative remapping of a 1-D variable
  SUBROUTINE remap_cons_2nd_order_1D( z_src, mask_src, d_src, z_dst, mask_dst, d_dst)
    ! 2nd-order conservative remapping of a 1-D variable
    !
    ! Used to remap ocean data from the provided vertical grid to the UFEMISM ocean vertical grid
    !
    ! Both z_src and z_dst can be irregular.
    !
    ! Both the src and dst data have a mask, with 0 indicating grid points where no data is defined.
    !
    ! This subroutine is serial, as it will be applied to single grid cells when remapping 3-D data fields,
    !   with the parallelisation being done by distributing the 2-D grid cells over the processes.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_src
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_src
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: d_src
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_dst
    INTEGER,  DIMENSION(:    ),          INTENT(IN)    :: mask_dst
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: d_dst

    ! Local variables:
    LOGICAL                                            :: all_are_masked
    INTEGER                                            :: nz_src, nz_dst
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: ddz_src
    INTEGER                                            :: k_src, k_dst
    REAL(dp)                                           :: zl_src, zu_src, zl_dst, zu_dst, z_lo, z_hi, z, d
    REAL(dp)                                           :: dz_overlap, dz_overlap_tot, d_int, d_int_tot
    REAL(dp)                                           :: dist_to_dst, dist_to_dst_min, max_dist
    INTEGER                                            :: k_src_nearest_to_dst

    ! Initialise
    d_dst = 0._dp

    ! Sizes
    nz_src = SIZE( z_src,1)
    nz_dst = SIZE( z_dst,1)

    ! Maximum distance on combined grids
    max_dist = MAXVAL([ ABS( z_src( nz_src) - z_src( 1)), &
                        ABS( z_dst( nz_dst) - z_dst( 1)), &
                        ABS( z_src( nz_src) - z_dst( 1)), &
                        ABS( z_dst( nz_dst) - z_src( 1))])

    ! Exception for when the entire src field is masked
    all_are_masked = .TRUE.
    DO k = 1, nz_src
      IF (mask_src( k) == 1) all_are_masked = .FALSE.
    END DO
    IF (all_are_masked) RETURN

    ! Exception for when the entire dst field is masked
    all_are_masked = .TRUE.
    DO k = 1, nz_dst
      IF (mask_dst( k) == 1) all_are_masked = .FALSE.
    END DO
    IF (all_are_masked) RETURN

    ! Calculate derivative d_src/dz (one-sided differencing at the boundary, central differencing everywhere else)
    ALLOCATE( ddz_src( nz_src))
    DO k = 2, nz_src-1
      ddz_src( k    ) = (d_src( k+1   ) - d_src( k-1     )) / (z_src( k+1   ) - z_src( k-1     ))
    END DO
    ddz_src(  1     ) = (d_src( 2     ) - d_src( 1       )) / (z_src( 2     ) - z_src( 1       ))
    ddz_src(  nz_src) = (d_src( nz_src) - d_src( nz_src-1)) / (z_src( nz_src) - z_src( nz_src-1))

    ! Perform conservative remapping by finding regions of overlap
    ! between source and destination grid cells

    DO k_dst = 1, nz_dst

      ! Skip masked grid cells
      IF (mask_dst( k_dst) == 0) THEN
        d_dst( k_dst) = 0._dp
        CYCLE
      END IF

      ! Find z range covered by this dst grid cell
      IF (k_dst > 1) THEN
        zl_dst = 0.5_dp * (z_dst( k_dst - 1) + z_dst( k_dst))
      ELSE
        zl_dst = z_dst( 1) - 0.5_dp * (z_dst( 2) - z_dst( 1))
      END IF
      IF (k_dst < nz_dst) THEN
        zu_dst = 0.5_dp * (z_dst( k_dst + 1) + z_dst( k_dst))
      ELSE
        zu_dst = z_dst( nz_dst) + 0.5_dp * (z_dst( nz_dst) - z_dst( nz_dst-1))
      END IF

      ! Find all overlapping src grid cells
      d_int_tot      = 0._dp
      dz_overlap_tot = 0._dp
      DO k_src = 1, nz_src

        ! Skip masked grid cells
        IF (mask_src( k_src) == 0) CYCLE

        ! Find z range covered by this src grid cell
        IF (k_src > 1) THEN
          zl_src = 0.5_dp * (z_src( k_src - 1) + z_src( k_src))
        ELSE
          zl_src = z_src( 1) - 0.5_dp * (z_src( 2) - z_src( 1))
        END IF
        IF (k_src < nz_src) THEN
          zu_src = 0.5_dp * (z_src( k_src + 1) + z_src( k_src))
        ELSE
          zu_src = z_src( nz_src) + 0.5_dp * (z_src( nz_src) - z_src( nz_src-1))
        END IF

        ! Find region of overlap
        z_lo = MAX( zl_src, zl_dst)
        z_hi = MIN( zu_src, zu_dst)
        dz_overlap = MAX( 0._dp, z_hi - z_lo)

        ! Calculate integral over region of overlap and add to sum
        IF (dz_overlap > 0._dp) THEN
          z = 0.5_dp * (z_lo + z_hi)
          d = d_src( k_src) + ddz_src( k_src) * (z - z_src( k_src))
          d_int = d * dz_overlap

          d_int_tot      = d_int_tot      + d_int
          dz_overlap_tot = dz_overlap_tot + dz_overlap
        END IF

      END DO ! DO k_src = 1, nz_src

      IF (dz_overlap_tot > 0._dp) THEN
        ! Calculate dst value
        d_dst( k_dst) = d_int_tot / dz_overlap_tot
      ELSE
        ! Exception for when no overlapping src grid cells were found; use nearest-neighbour extrapolation

        k_src_nearest_to_dst = 0._dp
        dist_to_dst_min      = max_dist
        DO k_src = 1, nz_src
          IF (mask_src( k_src) == 1) THEN
            dist_to_dst = ABS( z_src( k_src) - z_dst( k_dst))
            IF (dist_to_dst < dist_to_dst_min) THEN
              dist_to_dst_min      = dist_to_dst
              k_src_nearest_to_dst = k_src
            END IF
          END IF
        END DO

        ! Safety
        IF (k_src_nearest_to_dst == 0) THEN
          WRITE(0,*) '  remap_cons_2nd_order_1D - ERROR: couldnt find nearest neighbour on source grid!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF

        d_dst( k_dst) = d_src( k_src_nearest_to_dst)

      END IF ! IF (dz_overlap_tot > 0._dp) THEN

    END DO ! DO k_dst = 1, nz_dst

    ! Clean up after yourself
    DEALLOCATE( ddz_src)

  END SUBROUTINE remap_cons_2nd_order_1D

!   ! == Map data from a global lon/lat-grid to the model y/x-grid
!   SUBROUTINE map_glob_to_grid_2D( nlat, nlon, lat, lon, grid, d_glob, d_grid)
!     ! Map a data field from a global lat-lon grid to the regional square grid

!     IMPLICIT NONE

!     ! In/output variables:
!     INTEGER,                         INTENT(IN)  :: nlat
!     INTEGER,                         INTENT(IN)  :: nlon
!     REAL(dp), DIMENSION(nlat),       INTENT(IN)  :: lat
!     REAL(dp), DIMENSION(nlon),       INTENT(IN)  :: lon
!     TYPE(type_grid),                 INTENT(IN)  :: grid
!     REAL(dp), DIMENSION(:,:  ),      INTENT(IN)  :: d_glob
!     REAL(dp), DIMENSION(:,:  ),      INTENT(OUT) :: d_grid

!     ! Local variables:
!     INTEGER                                                :: i, j, il, iu, jl, ju
!     REAL(dp)                                               :: wil, wiu, wjl, wju

!     DO j = 1, grid%ny
!     DO i = grid%i1, grid%i2

!       ! Find enveloping lat-lon indices
!       il  = MAX(1,MIN(nlon-1, 1 + FLOOR((grid%lon(i,j)-MINVAL(lon)) / (lon(2)-lon(1)))))
!       iu  = il+1
!       wil = (lon(iu) - grid%lon(i,j))/(lon(2)-lon(1))
!       wiu = 1-wil

!       ! Exception for pixels near the zero meridian
!       IF (grid%lon(i,j) < MINVAL(lon)) THEN
!         il = nlon
!         iu = 1
!         wil = (lon(iu) - grid%lon(i,j))/(lon(2)-lon(1))
!         wiu = 1-wil
!       ELSEIF (grid%lon(i,j) > MAXVAL(lon)) THEN
!         il = nlon
!         iu = 1
!         wiu = (grid%lon(i,j) - lon(il))/(lon(2)-lon(1))
!         wil = 1-wiu
!       END IF

!       jl  = MAX(1,MIN(nlat-1, 1 + FLOOR((grid%lat(i,j)-MINVAL(lat)) / (lat(2)-lat(1)))))
!       ju  = jl+1
!       wjl = (lat(ju) - grid%lat(i,j))/(lat(2)-lat(1))
!       wju = 1-wjl

!       ! Interpolate data
!       d_grid( i,j) = (d_glob( il,jl) * wil * wjl) + &
!                      (d_glob( il,ju) * wil * wju) + &
!                      (d_glob( iu,jl) * wiu * wjl) + &
!                      (d_glob( iu,ju) * wiu * wju)

!     END DO
!     END DO
!     CALL sync

!   END SUBROUTINE map_glob_to_grid_2D

  subroutine map_glob_to_grid_3D( nlat, nlon, lat, lon, grid, d_glob, d_grid)
    ! Map a data field from a global lat-lon grid to the regional square grid

    implicit none

    ! In/output variables:
    integer,                    intent(in)  :: nlat
    integer,                    intent(in)  :: nlon
    real(dp), dimension(nlat),  intent(in)  :: lat
    real(dp), dimension(nlon),  intent(in)  :: lon
    type(type_grid),            intent(in)  :: grid
    real(dp), dimension(:,:,:), intent(in)  :: d_glob
    real(dp), dimension(:,:,:), intent(out) :: d_grid

    ! Local variables:
    integer                                 :: i, j, il, iu, jl, ju, k, nz
    real(dp)                                :: wil, wiu, wjl, wju

    nz = size( d_glob,3)

    ! Safety
    if (size(d_grid,3) /= nz) then
      call crash('Z-dimensions of global data and grid data fields do not match!')
    end if

    do j = 1, grid%ny
    do i = grid%i1, grid%i2

      ! Find enveloping lat-lon indices
      il  = max(1,min(nlon-1, 1 + floor((grid%lon(i,j)-minval(lon)) / (lon(2)-lon(1)))))
      iu  = il+1
      wil = (lon(iu) - grid%lon(i,j))/(lon(2)-lon(1))
      wiu = 1-wil

      ! Exception for pixels near the zero meridian
      if (grid%lon(i,j) < minval(lon)) then
        il = nlon
        iu = 1
        wil = (lon(iu) - grid%lon(i,j))/(lon(2)-lon(1))
        wiu = 1-wil
      elseif (grid%lon(i,j) > maxval(lon)) then
        il = nlon
        iu = 1
        wiu = (grid%lon(i,j) - lon(il))/(lon(2)-lon(1))
        wil = 1-wiu
      end if

      jl  = max(1,min(nlat-1, 1 + floor((grid%lat(i,j)-minval(lat)) / (lat(2)-lat(1)))))
      ju  = jl+1
      wjl = (lat(ju) - grid%lat(i,j))/(lat(2)-lat(1))
      wju = 1-wjl

      ! Interpolate data
      do k = 1, nz
        d_grid( i,j,k) = (d_glob( il,jl,k) * wil * wjl) + &
                         (d_glob( il,ju,k) * wil * wju) + &
                         (d_glob( iu,jl,k) * wiu * wjl) + &
                         (d_glob( iu,ju,k) * wiu * wju)
      end do

    end do
    end do

  end subroutine map_glob_to_grid_3D

  ! == Gaussian extrapolation (used for ocean data)
  subroutine extrapolate_Gaussian_floodfill( grid, mask, d, sigma, mask_filled)
    ! Extrapolate the data field d into the area designated by the mask,
    ! using Gaussian extrapolation of sigma
    !
    ! NOTE: not parallelised! This is done instead by dividing vertical
    !       ocean layers over the processes.
    !
    ! Note about the mask:
    !    2 = data provided
    !    1 = no data provided, fill allowed
    !    0 = no fill allowed
    ! (so basically this routine extrapolates data from the area
    !  where mask == 2 into the area where mask == 1)

    implicit none

    ! In/output variables:
    type(type_grid),            intent(in)    :: grid
    integer,  dimension(:,:  ), intent(in)    :: mask
    real(dp), dimension(:,:  ), intent(inout) :: d
    real(dp),                   intent(in)    :: sigma
    integer,  dimension(:,:  ), intent(out)   :: mask_filled   ! 1 = successfully filled, 2 = failed to fill (region of to-be-filled pixels not connected to source data)

    ! Local variables:
    integer                                   :: i,j,k,ii,jj,it
    integer                                   :: stackN1, stackN2
    integer,  dimension(:,:  ), allocatable   :: stack1, stack2
    integer,  dimension(:,:  ), allocatable   :: map
    integer                                   :: n_search
    logical                                   :: has_filled_neighbours
    integer                                   :: n
    real(dp)                                  :: sum_d, w, sum_w

    n_search = 1 + ceiling( 2._dp * sigma / grid%dx)

    ! Allocate map and stacks. Amount for memory is an estimation; if there
    ! are out-of-bounds errors, increase the memory for the stacks.
    allocate( map(       grid%nx,  grid%ny))
    allocate( stack1( 8*(grid%nx + grid%ny),2))
    allocate( stack2( 8*(grid%nx + grid%ny),2))

    map         = 0
    stack1      = 0
    stack2      = 0
    stackN1     = 0
    stackN2     = 0
    mask_filled = 0

    ! Initialise the map from the mask
    do j = 1, grid%ny
    do i = 1, grid%nx
      if (mask( i,j) == 2) then
        map( i,j) = 2
      end if
    end do
    end do

    ! Initialise the stack with all empty-next-to-filled grid cells
    do j = 1, grid%ny
    do i = 1, grid%nx

      if (mask( i,j) == 1) then
        ! This grid cell is empty and should be filled

        has_filled_neighbours = .false.
        do jj = max(1 ,j-1), min(grid%ny,j+1)
        do ii = max(1 ,i-1), min(grid%nx,i+1)
          if (mask( ii,jj) == 2) then
            has_filled_neighbours = .true.
            exit
          end if
        end do
        if (has_filled_neighbours) then
          exit
        end if
        end do

        if (has_filled_neighbours) then
          ! Add this empty-with-filled-neighbours grid cell to the stack,
          ! and mark it as stacked on the map
          map( i,j) = 1
          stackN2 = stackN2 + 1
          stack2( stackN2,:) = [i,j]
        end if

      end if

    end do
    end do

    ! Perform the flood-fill
    it = 0
    do while (stackN2 > 0)

      it = it + 1

      ! Go over all the stacked empty-next-to-filled grid cells, perform the
      ! Gaussian-kernel extrapolation to fill, and mark them as as such on the map
      do k = 1, stackN2

        ! Get grid cell indices
        i = stack2( k,1)
        j = stack2( k,2)

        ! Find Gaussian-weighted average value over nearby filled pixels within the basin
        n     = 0
        sum_d = 0._dp
        sum_w = 0._dp

        do jj = max( 1 ,j - n_search), min( grid%ny,j + n_search)
        do ii = max( 1 ,i - n_search), min( grid%nx,i + n_search)

          if (map( ii,jj) == 2) then
            n     = n + 1
            w     = exp( -0.5_dp * (sqrt(real(ii-i,dp)**2 + real(jj-j,dp)**2) / sigma)**2)
            sum_w = sum_w + w
            sum_d = sum_d + w * d( ii,jj)
          end if

        end do
        end do

        ! Fill in averaged value
        d( i,j) = sum_d / sum_w

        ! Mark grid cell as filled
        map( i,j) = 2
        mask_filled( i,j) = 1

      end do ! k = 1, stackN2

      ! Cycle the stacks
      stack1  = stack2
      stackN1 = stackN2
      stack2  = 0
      stackN2 = 0

      ! List new empty-next-to-filled grid cells
      do k = 1, stackN1

        ! Get grid cell indices
        i = stack1( k,1)
        j = stack1( k,2)

        ! Find empty neighbours; if unlisted, list them and mark them on the map
        do jj = max( 1, j-1), min( grid%ny, j+1)
        do ii = max( 1, i-1), min( grid%nx, i+1)

          if (map( ii,jj) == 0 .and. mask( ii,jj) == 1) then
            map( ii,jj) = 1
            stackN2 = stackN2 + 1
            stack2( stackN2,:) = [ii,jj]
          end if

        end do
        end do

      end do ! k = 1, stackN1

      ! Safety
      if (it > 2 * max( grid%ny, grid%nx)) then
        call crash('flood-fill got stuck!')
      end if

    end do ! while (stackN2 > 0)

    ! Mark grid cells that could not be filled
    do j = 1, grid%ny
    do i = 1, grid%nx
      if (mask_filled( i,j) == 0 .and. mask( i,j) == 1) then
        mask_filled( i,j) = 2
      end if
    end do
    end do

    ! Clean up after yourself
    deallocate( map   )
    deallocate( stack1)
    deallocate( stack2)

  end subroutine extrapolate_Gaussian_floodfill

!   ! == Interpolate ocean column data to a queried depth
  SUBROUTINE interpolate_ocean_depth( nz_ocean, z_ocean, f_ocean, z_query, f_query)
    ! Interpolate ocean column data to a queried depth using a simple bisection method.

    IMPLICIT NONE

    ! In/output variables:
    INTEGER,                             INTENT(IN)    :: nz_ocean    ! Number of vertical layers
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: z_ocean     ! Depth of layers (assumed to be monotonically increasing, does not need to be regularly spaced)
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: f_ocean     ! Value of whatever function we want to interpolate
    REAL(dp),                            INTENT(IN)    :: z_query     ! Depth at which we want to know the function
    REAL(dp),                            INTENT(OUT)   :: f_query     ! Interpolated function value

    ! Local variables:
    INTEGER                                            :: k_lo,k_hi,k_mid
    LOGICAL                                            :: foundit
    REAL(dp)                                           :: w

    ! Safety
    IF (z_query < 0._dp) THEN
      WRITE(0,*) '  interpolate_ocean_depth - ERROR: z_query = ', z_query, '< 0; cannot extrapolate above the sea surface, obviously!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSEIF (z_query > 12000._dp) THEN
      WRITE(0,*) '  interpolate_ocean_depth - ERROR: z_query = ', z_query, '> 12 km; the ocean is not that deep!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSEIF (SIZE(z_ocean,1) /= nz_ocean) THEN
      WRITE(0,*) '  interpolate_ocean_depth - ERROR: SIZE(z_ocean,1) = ', SIZE(z_ocean,1), ' /= nz_ocean = ', nz_ocean, '!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSEIF (SIZE(f_ocean,1) /= nz_ocean) THEN
      WRITE(0,*) '  interpolate_ocean_depth - ERROR: SIZE(f_ocean,1) = ', SIZE(f_ocean,1), ' /= nz_ocean = ', nz_ocean, '!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    ELSEIF (z_query > MAXVAL(z_ocean)) THEN
      !WRITE(0,*) '  interpolate_ocean_depth - ERROR: z_query = ', z_query, '> MAXVAL(z_ocean) = ', MAXVAL(z_ocean), '!'
      !CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

      ! Nearest-neighbour extrapolation when querying data beneath the end of the ocean data column
      f_query = f_ocean( nz_ocean)
      RETURN
    END IF

    ! Exception for when z_query = 0 (the World Ocean Atlas depth starts at 1.25...)
    IF (z_query < MINVAL(z_ocean)) THEN
      f_query = f_ocean(1)
      RETURN
    END IF

    ! Bisection method
    k_lo  = 1
    k_hi  = nz_ocean
    k_mid = INT( REAL(k_lo + k_hi,dp) / 2._dp)

    ! Exceptions
    IF     (ABS(z_query - z_ocean( k_lo )) < 1E-4_dp) THEN
      f_query = f_ocean( k_lo)
      RETURN
    ELSEIF (ABS(z_query - z_ocean( k_hi )) < 1E-4_dp) THEN
      f_query = f_ocean( k_hi)
      RETURN
    ELSEIF (ABS(z_query - z_ocean( k_mid)) < 1E-4_dp) THEN
      f_query = f_ocean( k_mid)
      RETURN
    END IF

    ! Bisection method
    foundit = .FALSE.
    DO WHILE (.NOT. foundit)

      IF (ABS(z_query - z_ocean( k_mid)) < 1E-4_dp) THEN
        ! Exception for when the queried depth is exactly at the midpoint index depth
        f_query = f_ocean( k_mid)
        RETURN
      ELSEIF (z_query > z_ocean( k_mid)) THEN
        ! Queried depth lies to the right of the midpoint
        k_lo = k_mid
        k_mid = INT( REAL(k_lo + k_hi,dp) / 2._dp)
      ELSE
        ! Queried depth lies to the left of the midpoint
        k_hi = k_mid
        k_mid = INT( REAL(k_lo + k_hi,dp) / 2._dp)
      END IF

      ! Stop iterating when endpoints lie next to each other; then just do linear interpolation between those two.
      IF (k_hi == k_lo+1) foundit = .TRUE.

    END DO ! DO WHILE (.NOT. foundit)

    ! Linear interpolation between nearest layers
    w = (z_query - z_ocean( k_lo)) / (z_ocean( k_hi) - z_ocean( k_lo))
    f_query = w * f_ocean( k_hi) + (1._dp - w) * f_ocean( k_lo)

  END SUBROUTINE interpolate_ocean_depth

  subroutine time_display( region, t_end, dt_ave, it)
    ! Little time display for the screen

    implicit none

    ! Input/Ouput variables
    type(type_model_region), intent(in)  :: region
    real(dp),                intent(in)  :: t_end
    real(dp),                intent(in)  :: dt_ave
    integer,                 intent(in)  :: it

    ! Local variables
    character(len=9)                     :: r_time, r_step, r_adv, r_ave

    if (region%time + region%dt < t_end) then
      r_adv = "no"
      write(r_time,"(F8.3)") min(region%time,t_end) / 1000._dp
      write(r_step,"(F6.3)") max(region%dt,0.001_dp)
      write(*,"(A)",advance=trim(r_adv)) "\r"// &
              "  t = " // trim(r_time) // " kyr - dt = " // trim(r_step) // " yr"
    else
      r_adv = "yes"
      write(r_time,"(F8.3)") min(region%time,t_end) / 1000._dp
      write(r_step, "(F6.3)") dt_ave / real(it-1,dp)
      write(*,"(A)",advance=trim(r_adv)) '\r'// &
            "  t = " // trim(r_time) // " kyr - dt_ave = " // trim(r_step) // " yr"
    end if
    if (region%do_output) then
      r_adv = "no"
      write(*,"(A)",advance=trim(r_adv)) "\r"
    end if

  end subroutine time_display

END MODULE utilities_module
