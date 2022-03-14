MODULE utilities_module

  ! Some generally useful tools

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: perr
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
  USE parallel_module,                 ONLY: allocate_shared_dist_int_0D, allocate_shared_dist_dp_0D, &
                                             allocate_shared_dist_int_1D, allocate_shared_dist_dp_1D, &
                                             allocate_shared_dist_int_2D, allocate_shared_dist_dp_2D, &
                                             allocate_shared_dist_int_3D, allocate_shared_dist_dp_3D, &
                                             allocate_shared_dist_bool_1D, &
                                             adapt_shared_dist_int_1D,    adapt_shared_dist_dp_1D, &
                                             adapt_shared_dist_int_2D,    adapt_shared_dist_dp_2D, &
                                             adapt_shared_dist_int_3D,    adapt_shared_dist_dp_3D, &
                                             adapt_shared_dist_bool_1D
  USE data_types_module,               ONLY: type_mesh, type_grid

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
  SUBROUTINE error_function(X, ERR)
    ! Purpose: Compute error function erf(x)
    ! Input:   x   --- Argument of erf(x)
    ! Output:  ERR --- erf(x)
    
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: X
    
    ! Output variables:
    REAL(dp), INTENT(OUT) :: ERR
    
    ! Local variables:
    REAL(dp)              :: EPS
    REAL(dp)              :: X2
    REAL(dp)              :: ER
    REAL(dp)              :: R
    REAL(dp)              :: C0
    INTEGER               :: k
    
    EPS = 1.0E-15_dp
    X2  = X * X
    IF(ABS(X) < 3.5_dp) THEN
     ER = 1.0_dp
     R  = 1.0_dp
     DO k = 1, 50
       R  = R * X2 / (REAL(k, dp) + 0.5_dp)
       ER = ER+R
       IF(ABS(R) < ABS(ER) * EPS) THEN
        C0  = 2.0_dp / SQRT(pi) * X * EXP(-X2)
        ERR = C0 * ER
        EXIT
       END IF
     END DO
    ELSE
     ER = 1.0_dp
     R  = 1.0_dp
     DO k = 1, 12
       R  = -R * (REAL(k, dp) - 0.5_dp) / X2
       ER = ER + R
       C0  = EXP(-X2) / (ABS(X) * SQRT(pi))
       ERR = 1.0_dp - C0 * ER
       IF(X < 0.0_dp) ERR = -ERR
     END DO
    ENDIF

    RETURN
  END SUBROUTINE error_function
  
! == The oblique stereographic projection 
  SUBROUTINE oblique_sg_projection( lambda, phi, lambda_M_deg, phi_M_deg, alpha_deg, x_IM_P_prime, y_IM_P_prime)
    ! This subroutine projects with an oblique stereographic projection the longitude-latitude
    ! coordinates which coincide with the GCM grid points to the rectangular IM coordinate
    ! system, with coordinates (x,y).
    !
    ! For more information about M, C%alpha_stereographic, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)            :: lambda        ! lon in degrees
    REAL(dp), INTENT(IN)            :: phi           ! lat in degrees
    
    ! Polar stereographic projection parameters
    REAL(dp), INTENT(IN)            :: lambda_M_deg  ! in degrees
    REAL(dp), INTENT(IN)            :: phi_M_deg     ! in degrees
    REAL(dp), INTENT(IN)            :: alpha_deg     ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT)           :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(OUT)           :: y_IM_P_prime  ! in meter

    ! Local variables:
    REAL(dp)                        :: phi_P         ! in radians
    REAL(dp)                        :: lambda_P      ! in radians
    REAL(dp)                        :: t_P_prime
    
    REAL(dp)                        :: lambda_M, phi_M, alpha
    
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! For North and South Pole: C%lambda_M = 0._dp, to generate the correct IM coordinate
    ! system, see equation (2.3) or equation (A.53) in Reerink et al. (2010).

    ! Convert longitude-latitude coordinates to radians:
    phi_P    = (pi / 180._dp) * phi
    lambda_P = (pi / 180._dp) * lambda

    ! See equation (2.6) or equation (A.56) in Reerink et al. (2010):
    t_P_prime = ((1._dp + COS(alpha)) / (1._dp + COS(phi_P) * COS(phi_M) * COS(lambda_P - lambda_M) + SIN(phi_P) * SIN(phi_M))) / (pi / 180._dp)

    ! See equations (2.4-2.5) or equations (A.54-A.55) in Reerink et al. (2010):
    x_IM_P_prime =  earth_radius * (COS(phi_P) * SIN(lambda_P - lambda_M)) * t_P_prime
    y_IM_P_prime =  earth_radius * (SIN(phi_P) * COS(phi_M) - (COS(phi_P) * SIN(phi_M)) * COS(lambda_P - lambda_M)) * t_P_prime   
    
    
  END SUBROUTINE oblique_sg_projection
  SUBROUTINE inverse_oblique_sg_projection( x_IM_P_prime, y_IM_P_prime, lambda_M_deg, phi_M_deg, alpha_deg, lambda_P, phi_P)
    ! This subroutine projects with an inverse oblique stereographic projection the
    ! (x,y) coordinates which coincide with the IM grid points to the longitude-latitude
    ! coordinate system, with coordinates (lambda, phi) in degrees.
    !
    ! For more information about M, alpha, the center of projection and the used
    ! projection method see:
    !  Reerink et al. (2010), Mapping technique of climate fields between GCM's and ice models, GMD
    
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), INTENT(IN)  :: x_IM_P_prime  ! in meter
    REAL(dp), INTENT(IN)  :: y_IM_P_prime  ! in meter
    
    ! Polar stereographic projection parameters
    REAL(dp), INTENT(IN)  :: lambda_M_deg  ! in degrees
    REAL(dp), INTENT(IN)  :: phi_M_deg     ! in degrees
    REAL(dp), INTENT(IN)  :: alpha_deg     ! in degrees

    ! Output variables:
    REAL(dp), INTENT(OUT) :: lambda_P      ! in degrees
    REAL(dp), INTENT(OUT) :: phi_P         ! in degrees

    ! Local variables:
    REAL(dp)              :: x_3D_P_prime  ! in meter
    REAL(dp)              :: y_3D_P_prime  ! in meter
    REAL(dp)              :: z_3D_P_prime  ! in meter
    REAL(dp)              :: a
    REAL(dp)              :: t_P
    REAL(dp)              :: x_3D_P        ! in meter
    REAL(dp)              :: y_3D_P        ! in meter
    REAL(dp)              :: z_3D_P        ! in meter
    
    REAL(dp)              :: lambda_M, phi_M, alpha
    
    lambda_M = (pi / 180._dp) * lambda_M_deg
    phi_M    = (pi / 180._dp) * phi_M_deg
    alpha    = (pi / 180._dp) * alpha_deg

    ! See equations (2.14-2.16) or equations (B.21-B.23) in Reerink et al. (2010):
    x_3D_P_prime = earth_radius * COS(alpha) * COS(lambda_M) * COS(phi_M) - SIN(lambda_M) * (x_IM_P_prime*(pi / 180._dp)) - COS(lambda_M) * SIN(phi_M) * (y_IM_P_prime*(pi / 180._dp))
    y_3D_P_prime = earth_radius * COS(alpha) * SIN(lambda_M) * COS(phi_M) + COS(lambda_M) * (x_IM_P_prime*(pi / 180._dp)) - SIN(lambda_M) * SIN(phi_M) * (y_IM_P_prime*(pi / 180._dp))
    z_3D_P_prime = earth_radius * COS(alpha) *                 SIN(phi_M)                                          +                   COS(phi_M) * (y_IM_P_prime*(pi / 180._dp))

    ! See equation (2.13) or equation (B.20) in Reerink et al. (2010):
    a = COS(lambda_M) * COS(phi_M) * x_3D_P_prime  +  SIN(lambda_M) * COS(phi_M) * y_3D_P_prime  +  SIN(phi_M) * z_3D_P_prime

    ! See equation (2.12) or equation (B.19) in Reerink et al. (2010):
    t_P = (2._dp * earth_radius**2 + 2._dp * earth_radius * a) / (earth_radius**2 + 2._dp * earth_radius * a + x_3D_P_prime**2 + y_3D_P_prime**2 + z_3D_P_prime**2)

    ! See equations (2.9-2.11) or equations (B.16-B.18) in Reerink et al. (2010):
    x_3D_P =  earth_radius * COS(lambda_M) * COS(phi_M) * (t_P - 1._dp) + x_3D_P_prime * t_P
    y_3D_P =  earth_radius * SIN(lambda_M) * COS(phi_M) * (t_P - 1._dp) + y_3D_P_prime * t_P
    z_3D_P =  earth_radius *                   SIN(phi_M) * (t_P - 1._dp) + z_3D_P_prime * t_P

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
  
! == Map data between two square grids using 2nd-order conservative remapping
  SUBROUTINE map_square_to_square_cons_2nd_order_2D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src, d_dst)
    ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)
    
    IMPLICIT NONE
  
    ! Input and output variables
    INTEGER,                            INTENT(IN)    :: nx_src
    INTEGER,                            INTENT(IN)    :: ny_src
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_src
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_src
    INTEGER,                            INTENT(IN)    :: nx_dst
    INTEGER,                            INTENT(IN)    :: ny_dst
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_dst
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_dst
    REAL(dp), DIMENSION(:,:  ),         INTENT(IN)    :: d_src
    REAL(dp), DIMENSION(:,:  ),         INTENT(OUT)   :: d_dst
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_square_to_square_cons_2nd_order_2D'
    INTEGER                                           :: i,j,i_src,j_src,i1,i2,igmin,igmax,jgmin,jgmax,j1,j2
    REAL(dp)                                          :: dx_src, dy_src, dx_dst, dy_dst, xcmin, xcmax, ycmin, ycmax
    INTEGER,  DIMENSION(nx_dst,2)                     :: ir_src
    INTEGER,  DIMENSION(ny_dst,2)                     :: jr_src
    REAL(dp)                                          :: xomin, xomax, yomin, yomax, w0, w1x, w1y
    REAL(dp)                                          :: Ad, Asd, Asum
    REAL(dp), DIMENSION(:,:  ), POINTER               ::  ddx_src,  ddy_src
    INTEGER                                           :: wddx_src, wddy_src
    INTEGER,  DIMENSION(:,:  ), POINTER               ::  mask_dst_outside_src
    INTEGER                                           :: wmask_dst_outside_src
    REAL(dp)                                          :: LI_mxydx1, LI_mxydx2, LI_mxydx3, LI_mxydx4
    REAL(dp)                                          :: LI_xydy1, LI_xydy2, LI_xydy3, LI_xydy4
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D(  ny_src, nx_src, ddx_src,              wddx_src             )
    CALL allocate_shared_dp_2D(  ny_src, nx_src, ddy_src,              wddy_src             )
    CALL allocate_shared_int_2D( ny_dst, nx_dst, mask_dst_outside_src, wmask_dst_outside_src)
    
    ! Find grid spacings
    dx_src = x_src(2) - x_src(1)
    dy_src = y_src(2) - y_src(1)
    dx_dst = x_dst(2) - x_dst(1)
    dy_dst = y_dst(2) - y_dst(1)
    Ad = dx_dst * dy_dst
    
    ! If the grids are equal, the solution is trivial; just copy the data
    IF (dx_src == dx_dst .AND. dy_src == dy_dst .AND. nx_src == nx_dst .AND. ny_src == ny_dst) THEN
      CALL partition_list( nx_dst, par%i, par%n, i1, i2)
      d_dst( :,i1:i2) = d_src( :,i1:i2)
      CALL sync
      RETURN
    END IF
    
    ! Find overlaps between grids
    DO i = 1, nx_dst
      ! Dst cell i overlaps with src cells ir_src( i,1) to ir_src( i,2)
      xcmin = x_dst( i) - dx_dst/2._dp
      xcmax = x_dst( i) + dx_dst/2._dp
      ir_src( i,:) = MAX( 1, MIN( nx_src, [CEILING(-1.5_dp + FLOOR(nx_src/2._dp) + xcmin / dx_src), &
                                           CEILING( 1.5_dp + FLOOR(nx_src/2._dp) + xcmax / dx_src)] ))
    END DO ! DO i = 1, nx_dst
    DO j = 1, ny_dst
      ! Dst cell j overlaps with src cells jr_src( j,1) to ir_src( j,2)
      ycmin = y_dst( j) - dy_dst/2._dp
      ycmax = y_dst( j) + dy_dst/2._dp
      jr_src( j,:) = MAX( 1, MIN( ny_src, [CEILING(-1.5_dp + FLOOR(ny_src/2._dp) + ycmin / dy_src), &
                                           CEILING( 1.5_dp + FLOOR(ny_src/2._dp) + ycmax / dy_src)] ))
    END DO ! DO j = 1, ny_dst
    
    ! Get derivatives of d_src
    CALL partition_list( nx_src, par%i, par%n, i1, i2)
    DO i = MAX(2,i1), MIN(nx_src-1,i2)
    DO j = 2, ny_src-1
      ddx_src( j,i) = (d_src( j,i+1) - d_src( j,i-1)) / (2._dp * dx_src)
      ddy_src( j,i) = (d_src( j+1,i) - d_src( j-1,i)) / (2._dp * dy_src)
    END DO
    END DO
    CALL sync
    
    ! Find parallelisation domains
    CALL partition_list( nx_dst, par%i, par%n, i1, i2)
    CALL partition_list( ny_dst, par%i, par%n, j1, j2)
    
    DO i = i1, i2
    DO j = 1, ny_dst
      
      d_dst(                j,i) = 0._dp
      mask_dst_outside_src( j,i) = 0
      Asum                       = 0._dp
      
      DO i_src = ir_src( i,1), ir_src( i,2)
      DO j_src = jr_src( j,1), jr_src( j,2)
        
        xomin = MAX( x_dst( i) - dx_dst/2._dp, x_src( i_src) - dx_src/2._dp)
        xomax = MIN( x_dst( i) + dx_dst/2._dp, x_src( i_src) + dx_src/2._dp)
        yomin = MAX( y_dst( j) - dy_dst/2._dp, y_src( j_src) - dy_src/2._dp)
        yomax = MIN( y_dst( j) + dy_dst/2._dp, y_src( j_src) + dy_src/2._dp)
        
        IF (xomax <= xomin .OR. yomax <= yomin) CYCLE
        
        Asd  = (xomax - xomin) * (yomax - yomin)
        Asum = Asum + Asd
        
        w0  = Asd / Ad
        
        CALL line_integral_mxydx( [xomin,yomin], [xomax,yomin], 1E-9_dp, LI_mxydx1)
        CALL line_integral_mxydx( [xomax,yomin], [xomax,yomax], 1E-9_dp, LI_mxydx2)
        CALL line_integral_mxydx( [xomax,yomax], [xomin,yomax], 1E-9_dp, LI_mxydx3)
        CALL line_integral_mxydx( [xomin,yomax], [xomin,yomin], 1E-9_dp, LI_mxydx4)
        
        w1x = 1._dp / Ad * (LI_mxydx1 + LI_mxydx2 + LI_mxydx3 + LI_mxydx4) - w0 * x_src( i_src)
        
        CALL line_integral_xydy(  [xomin,yomin], [xomax,yomin], 1E-9_dp, LI_xydy1)
        CALL line_integral_xydy(  [xomax,yomin], [xomax,yomax], 1E-9_dp, LI_xydy2)
        CALL line_integral_xydy(  [xomax,yomax], [xomin,yomax], 1E-9_dp, LI_xydy3)
        CALL line_integral_xydy(  [xomin,yomax], [xomin,yomin], 1E-9_dp, LI_xydy4)
        
        w1y = 1._dp / Ad * (LI_xydy1  + LI_xydy2  + LI_xydy3  + LI_xydy4 ) - w0 * y_src( j_src)
        
        d_dst( j,i) = d_dst( j,i) + w0  * d_src(   j_src,i_src) + &
                                    w1x * ddx_src( j_src,i_src) + &
                                    w1y * ddy_src( j_src,i_src)
        
      END DO ! DO j_src = jr_src( j,1), jr_src( j,2)
      END DO ! DO i_src = ir_src( i,1), ir_src( i,2)
      
      IF (Asum < Ad) mask_dst_outside_src( j,i) = 1
      
    END DO ! DO j = 1, ny_dst
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Use nearest-neighbour extrapolation for dst cells outside of the src grid
    ! =========================================================================
    
    ! Find the range of grid cells that were mapped correctly
    igmin = 0
    igmax = 0
    jgmin = 0
    jgmax = 0
    
    j = INT( REAL(ny_dst,dp)/2._dp)
    DO i = 1, nx_dst
      IF (mask_dst_outside_src( j,i) == 0) THEN
        igmin = i
        EXIT
      END IF
    END DO
    DO i = nx_dst, 1, -1
      IF (mask_dst_outside_src( j,i) == 0) THEN
        igmax = i
        EXIT
      END IF
    END DO
    
    i = INT( REAL(nx_dst,dp)/2._dp)
    DO j = 1, ny_dst
      IF (mask_dst_outside_src( j,i) == 0) THEN
        jgmin = j
        EXIT
      END IF
    END DO
    DO j = ny_dst, 1, -1
      IF (mask_dst_outside_src( j,i) == 0) THEN
        jgmax = j
        EXIT
      END IF
    END DO
    
    ! Corners
    IF (par%master) THEN
      ! Southwest
      d_dst( 1      :jgmin-1 ,1      :igmin-1) = d_dst( jgmin,igmin)
      ! Southeast
      d_dst( 1      :jgmin-1 ,igmax+1:nx_dst ) = d_dst( jgmin,igmax)
      ! Northwest
      d_dst( jgmax+1:ny_dst  ,1      :igmin-1) = d_dst( jgmax,igmin)
      ! Northeast
      d_dst( jgmax+1:ny_dst  ,igmax+1:nx_dst ) = d_dst( jgmax,igmax)
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Borders
    DO i = MAX(i1,igmin), MIN(i2,igmax)
      ! South
      d_dst( 1      :jgmin-1,i) = d_dst( jgmin,i)
      ! North
      d_dst( jgmax+1:ny_dst ,i) = d_dst( jgmax,i)
    END DO
    DO j = MAX(j1,jgmin), MIN(j2,jgmax)
      ! West
      d_dst( j,1      :igmin-1) = d_dst( j,igmin)
      ! East
      d_dst( j,igmax+1:nx_dst ) = d_dst( j,igmax)
    END DO
    CALL sync
    
    ! Clean up after yourself
    CALL deallocate_shared( wddx_src             )
    CALL deallocate_shared( wddy_src             )
    CALL deallocate_shared( wmask_dst_outside_src)
  
  END SUBROUTINE map_square_to_square_cons_2nd_order_2D
  SUBROUTINE map_square_to_square_cons_2nd_order_3D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src, d_dst)
    ! Map data from one square grid to another (e.g. PD ice thickness from the square grid in the input file to the model square grid)
    
    IMPLICIT NONE
  
    ! Input and output variables
    INTEGER,                            INTENT(IN)    :: nx_src
    INTEGER,                            INTENT(IN)    :: ny_src
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_src
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_src
    INTEGER,                            INTENT(IN)    :: nx_dst
    INTEGER,                            INTENT(IN)    :: ny_dst
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: x_dst
    REAL(dp), DIMENSION(:    ),         INTENT(IN)    :: y_dst
    REAL(dp), DIMENSION(:,:,:),         INTENT(IN)    :: d_src
    REAL(dp), DIMENSION(:,:,:),         INTENT(OUT)   :: d_dst
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'map_square_to_square_cons_2nd_order_3D'
    INTEGER                                           :: nz, k, i1_src, i2_src, i1_dst, i2_dst
    REAL(dp), DIMENSION(:,:  ), POINTER               ::  d_src_2D,  d_dst_2D
    INTEGER                                           :: wd_src_2D, wd_dst_2D
    
    nz = SIZE( d_src,1)
    
    ! Safety
    IF (SIZE(d_dst,1) /= nz) THEN
      IF (par%master) WRITE(0,*) 'map_square_to_square_cons_2nd_order_3D - ERROR: nz_src /= nz_dst!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    CALL partition_list( nx_src, par%i, par%n, i1_src, i2_src)
    CALL partition_list( nx_dst, par%i, par%n, i1_dst, i2_dst)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( ny_src, nx_src, d_src_2D, wd_src_2D)
    CALL allocate_shared_dp_2D( ny_dst, nx_dst, d_dst_2D, wd_dst_2D)
    
    ! Remap the fields one layer at a time
    DO k = 1, nz
      d_src_2D(   :,i1_src:i2_src) = d_src(    k,:,i1_src:i2_src)
      CALL map_square_to_square_cons_2nd_order_2D( nx_src, ny_src, x_src, y_src, nx_dst, ny_dst, x_dst, y_dst, d_src_2D, d_dst_2D)
      d_dst(    k,:,i1_dst:i2_dst) = d_dst_2D(   :,i1_dst:i2_dst)
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_src_2D)
    CALL deallocate_shared( wd_dst_2D)
  
  END SUBROUTINE map_square_to_square_cons_2nd_order_3D
  
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
    CALL allocate_shared_dp_2D( grid%nx + 2*n, grid%ny + 2*n, d_ext,        wd_ext       )
    CALL allocate_shared_dp_2D( grid%nx + 2*n, grid%ny + 2*n, d_ext_smooth, wd_ext_smooth)

    ! Copy data to the extended array and fill in the margins
    DO j = 1, grid%ny
    DO i = grid%i1, grid%i2
      d_ext( i+n,j+n) = d( i,j)
    END DO
    END DO
    IF (par%master) THEN
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
    END IF
    CALL sync

    ! Convolute extended data with the smoothing filter
    d_ext_smooth( grid%i1+n:grid%i2+n,:) = 0._dp
    CALL sync

    DO j = 1,       grid%ny
    DO i = grid%i1, grid%i2
      DO ii = -n, n
        d_ext_smooth( i+n,j+n) = d_ext_smooth( i+n,j+n) + d_ext( i+n+ii,j+n) * f(ii)
      END DO
    END DO
    END DO
    CALL sync

    d_ext( grid%i1+n:grid%i2+n,:) = d_ext_smooth( grid%i1+n:grid%i2+n,:)
    CALL sync

    DO j = grid%j1, grid%j2
      d_ext(           1:          n,j) = d( 1      ,j)
      d_ext( grid%nx+n+1:grid%nx+2*n,j) = d( grid%nx,j)
    END DO
    CALL sync

    d_ext_smooth( grid%i1+n:grid%i2+n,:) = 0._dp
    CALL sync

    DO j = grid%j1, grid%j2
    DO i = 1,       grid%nx
      DO jj = -n, n
        d_ext_smooth( i+n,j+n) = d_ext_smooth( i+n,j+n) + d_ext( i+n,j+n+jj) * f(jj)
      END DO
    END DO
    END DO
    CALL sync

    ! Copy data back
    DO j = 1, grid%ny
    DO i = grid%i1, grid%i2
      d( i,j) = d_ext_smooth( i+n, j+n)
    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    DEALLOCATE( f)
    CALL deallocate_shared( wd_ext)
    CALL deallocate_shared( wd_ext_smooth)

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
    INTEGER                                            :: wd_2D
    
    ! Allocate temporary shared memory for the extended and smoothed data fields
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_2D, wd_2D)
    
    DO k = 1, SIZE( d,3)
      d_2D( grid%i1:grid%i2,:) = d( grid%i1:grid%i2,:,k)
      CALL smooth_Gaussian_2D_grid( grid, d_2D, r)
      d( grid%i1:grid%i2,:,k) = d_2D( grid%i1:grid%i2,:)
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_2D)
    
  END SUBROUTINE smooth_Gaussian_3D_grid
  SUBROUTINE smooth_Shepard_2D_grid( grid, d, r)
    ! Apply a Shepard smoothing filter of with sigma = n*dx to the 2D data field d
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in m
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Shepard_2D_grid'
    INTEGER                                            :: i,j,k,l,n
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_ext,  d_ext_smooth
    INTEGER                                            :: wd_ext, wd_ext_smooth
    REAL(dp)                                           :: ShepNumSum     ! The sum of the numerators in the Shepard weighting
    REAL(dp)                                           :: ShepDenSum     ! The sum of the denumerators in the Shepard weighting
    REAL(dp)                                           :: distance       ! in gridsize units
    REAL(dp)                                           :: smooth_radius  ! in gridsize units
    REAL(dp)                                           :: exponent       ! The distance weighting exponent in the Shepard weighting
    
    n = CEILING( r / grid%dx) ! Number of cells to extend the data by (3 standard deviations is enough to capture the tails of the normal distribution)
    
    smooth_radius = REAL(n,dp)
    exponent      = 2._dp
    
    ! Allocate temporary shared memory for the extended and smoothed data fields
    CALL allocate_shared_dp_2D( grid%ny + 2*n, grid%nx + 2*n, d_ext,        wd_ext       )
    CALL allocate_shared_dp_2D( grid%ny + 2*n, grid%nx + 2*n, d_ext_smooth, wd_ext_smooth)
    
    ! Copy data to the extended array and fill in the margins
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      d_ext( j+n,i+n) = d( j,i)
    END DO
    END DO
    IF (par%master) THEN
      ! West
      d_ext( n+1:n+grid%ny, 1            ) = d( :      ,1      )
      ! East
      d_ext( n+1:n+grid%ny, grid%nx+2*n  ) = d( :      ,grid%nx)
      ! South
      d_ext( 1            , n+1:n+grid%nx) = d( 1      ,:      )
      ! North
      d_ext( grid%ny+2*n  , n+1:n+grid%nx) = d( grid%ny,:      )
      ! Corners
      d_ext( 1:n,                     1:n                    ) = d( 1      ,1      )
      d_ext( 1:n,                     grid%nx+n+1:grid%nx+2*n) = d( 1      ,grid%nx)
      d_ext( grid%ny+n+1:grid%ny+2*n, 1:n                    ) = d( grid%ny,1      )
      d_ext( grid%ny+n+1:grid%ny+2*n, grid%nx+n+1:grid%nx+2*n) = d( grid%ny,grid%nx)
    END IF
    CALL sync
    
    ! Convolute extended data with the smoothing filter
    d_ext_smooth( :,grid%i1+n:grid%i2+n) = 0._dp
    CALL sync
    
    DO i = grid%i1, grid%i2
    DO j = 1,       grid%ny
      ShepNumSum = 0._dp
      ShepDenSum = 0._dp
      DO k = -n, n
      DO l = -n, n
        distance = SQRT(REAL(k,dp)**2 + REAL(l,dp)**2)
        IF (distance <= smooth_radius .AND. distance > 0._dp) THEN
          ShepNumSum = ShepNumSum + (d_ext( j+n+l,i+n+k)/(distance**exponent))
          ShepDenSum = ShepDenSum + (1._dp/(distance**exponent))
        END IF
      END DO
      END DO
      d_ext_smooth( j+n,i+n) = ShepNumSum/ShepDenSum
    END DO
    END DO
    CALL sync
    
    ! Copy data back
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      d( j,i) = d_ext_smooth( j+n, i+n)
    END DO
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_ext)
    CALL deallocate_shared( wd_ext_smooth)
    
  END SUBROUTINE smooth_Shepard_2D_grid
  SUBROUTINE smooth_Shepard_3D_grid( grid, d, r)
    ! Apply a Shepard smoothing filter of with sigma = n*dx to the 3D data field d
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in km
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'smooth_Shepard_3D_grid'
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_2D
    INTEGER                                            :: wd_2D
    
    ! Allocate temporary shared memory for the extended and smoothed data fields
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_2D, wd_2D)
    
    DO k = 1, SIZE( d,2)
      d_2D( :,grid%i1:grid%i2) = d( k,:,grid%i1:grid%i2)
      CALL smooth_Shepard_2D_grid( grid, d_2D, r)
      d( k,:,grid%i1:grid%i2) = d_2D( :,grid%i1:grid%i2)
    END DO
    
    ! Clean up after yourself
    CALL deallocate_shared( wd_2D)
    
  END SUBROUTINE smooth_Shepard_3D_grid
  
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
    
    IF (par%master) THEN
      
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
      
    END IF ! IF (par%master) THEN
    CALL sync
    
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
  FUNCTION tridiagonal_solve( ldiag, diag, udiag, rhs) RESULT(x)
    ! Lapack tridiagnal solver (in double precision):
    ! Matrix system solver for tridiagonal matrices. 
    ! Used e.g. in solving the ADI scheme. 
    ! ldiag = lower diagonal elements (j,j-1) of the matrix
    ! diag  = diagonal elements (j,j) of the matrix
    ! udiag = upper diagonal elements (j,j+1) of the matrix
    ! rhs   = right hand side of the matrix equation in the ADI scheme
    
    IMPLICIT NONE

    ! Input variables:
    REAL(dp), DIMENSION(:),            INTENT(IN) :: diag
    REAL(dp), DIMENSION(SIZE(diag)-1), INTENT(IN) :: udiag, ldiag
    REAL(dp), DIMENSION(SIZE(diag)),   INTENT(IN) :: rhs

    ! Result variable:
    REAL(dp), DIMENSION(SIZE(diag))               :: x
    
    ! Local variables:     
    INTEGER                                       :: info
    REAL(dp), DIMENSION(SIZE(diag))               :: diag_copy
    REAL(dp), DIMENSION(SIZE(udiag))              :: udiag_copy, ldiag_copy

    ! The LAPACK solver will overwrite the rhs with the solution x. Therefore we 
    ! first copy the rhs in the solution vector x:
    x = rhs

    ! The LAPACK solver will change the elements in the matrix, therefore we copy them:
    diag_copy  =  diag
    udiag_copy = udiag
    ldiag_copy = ldiag

    CALL DGTSV(SIZE(diag), 1, ldiag_copy, diag_copy, udiag_copy, x, SIZE(diag), info)
    
    IF (info /= 0) THEN
      WRITE(0,*) 'tridiagonal_solve - ERROR: LAPACK solver DGTSV returned error message info = ', info
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END FUNCTION tridiagonal_solve
  
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
    IF (ABS( detA) < TINY( detA)) THEN
      !PRINT *, A(1,1), A(1,2)
      !PRINT *, A(2,1), A(2,2)
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
    IF (ABS( detA) < TINY( detA)) THEN
      ! PRINT *, A(1,1), A(1,2), A(1,3)
      ! PRINT *, A(2,1), A(2,2), A(2,3)
      ! PRINT *, A(3,1), A(3,2), A(3,3)
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
  SUBROUTINE calc_matrix_inverse_general( A, Ainv)
    ! Calculate the inverse Ainv of an n-by-n matrix A using LAPACK
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: A
    REAL(dp), DIMENSION(:,:  ),          INTENT(INOUT) :: Ainv
    
    ! Local variables:
    REAL(dp), DIMENSION( SIZE( A,1))                   :: work     ! work array for LAPACK
    INTEGER,  DIMENSION( SIZE( A,1))                   :: ipiv     ! pivot indices
    INTEGER                                            :: n,info

    n = size( A,1)
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    
    ! DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
    CALL DGETRF( n, n, Ainv, n, ipiv, info)
    
    ! Safety
    IF (info /= 0) THEN
      ! IF (n == 5) THEN
        ! PRINT *, A(1,1), A(1,2), A(1,3), A(1,4), A(1,5)
        ! PRINT *, A(2,1), A(2,2), A(2,3), A(2,4), A(2,5)
        ! PRINT *, A(3,1), A(3,2), A(3,3), A(3,4), A(3,5)
        ! PRINT *, A(4,1), A(4,2), A(4,3), A(4,4), A(4,5)
        ! PRINT *, A(5,1), A(5,2), A(5,3), A(5,4), A(5,5)
      ! END IF
      WRITE(0,*) 'calc_matrix_inverse_general - DGETRF error: matrix is numerically singular!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
    CALL DGETRI( n, Ainv, n, ipiv, work, n, info)
    
    ! Safety
    IF (info /= 0) THEN
      WRITE(0,*) 'calc_matrix_inverse_general - DGETRI error: matrix inversion failed!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE calc_matrix_inverse_general
  
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

END MODULE utilities_module
