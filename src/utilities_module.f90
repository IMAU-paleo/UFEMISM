MODULE utilities_module

  ! Some generally useful tools

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
  USE parallel_module,                 ONLY: allocate_shared_dist_int_0D, allocate_shared_dist_dp_0D, &
                                             allocate_shared_dist_int_1D, allocate_shared_dist_dp_1D, &
                                             allocate_shared_dist_int_2D, allocate_shared_dist_dp_2D, &
                                             allocate_shared_dist_int_3D, allocate_shared_dist_dp_3D, &
                                             allocate_shared_dist_bool_1D, &
                                             adapt_shared_dist_int_1D,    adapt_shared_dist_dp_1D, &
                                             adapt_shared_dist_int_2D,    adapt_shared_dist_dp_2D, &
                                             adapt_shared_dist_int_3D,    adapt_shared_dist_dp_3D, &
                                             adapt_shared_dist_bool_1D
  USE data_types_module,               ONLY: type_mesh, type_grid, type_sparse_matrix_CSR

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
    
    ! Local variables
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
    
    ! Local variables
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
      DO jj = -n, n
        d_ext_smooth( j+n,i+n) = d_ext_smooth( j+n,i+n) + d_ext( j+n+jj,i+n) * f(jj)
      END DO
    END DO
    END DO
    CALL sync
    
    d_ext( :,grid%i1+n:grid%i2+n) = d_ext_smooth( :,grid%i1+n:grid%i2+n)
    CALL sync
    
    DO j = grid%j1, grid%j2
      d_ext( j,           1:          n) = d( j,1      )
      d_ext( j, grid%nx+n+1:grid%nx+2*n) = d( j,grid%nx)
    END DO
    CALL sync
    
    d_ext_smooth( :,grid%i1+n:grid%i2+n) = 0._dp
    CALL sync
    
    DO j = grid%j1, grid%j2
    DO i = 1,       grid%nx
      DO ii = -n, n
        d_ext_smooth( j+n,i+n) = d_ext_smooth( j+n,i+n) + d_ext( j+n,i+n+ii) * f(ii)
      END DO
    END DO
    END DO
    CALL sync
    
    ! Copy data back
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      d( j,i) = d_ext_smooth( j+n, i+n)
    END DO
    END DO
    CALL sync
    
    ! Clean up after yourself
    DEALLOCATE( f)
    CALL deallocate_shared( wd_ext)
    CALL deallocate_shared( wd_ext_smooth)
    
  END SUBROUTINE smooth_Gaussian_2D_grid
  SUBROUTINE smooth_Gaussian_3D_grid( grid, d, r, nz)
    ! Apply a Gaussian smoothing filter of with sigma = n*dx to the 3D data field d
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in km
    INTEGER,                             INTENT(IN)    :: nz
    
    ! Local variables:
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_2D
    INTEGER                                            :: wd_2D
    
    ! Allocate temporary shared memory for the extended and smoothed data fields
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_2D, wd_2D)
    
    DO k = 1, nz
      d_2D( :,grid%i1:grid%i2) = d( k,:,grid%i1:grid%i2)
      CALL smooth_Gaussian_2D_grid( grid, d_2D, r)
      d( k,:,grid%i1:grid%i2) = d_2D( :,grid%i1:grid%i2)
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
  SUBROUTINE smooth_Shepard_3D_grid( grid, d, r, nz)
    ! Apply a Shepard smoothing filter of with sigma = n*dx to the 3D data field d
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                     INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:,:),          INTENT(INOUT) :: d
    REAL(dp),                            INTENT(IN)    :: r      ! Smoothing radius in km
    INTEGER,                             INTENT(IN)    :: nz
    
    ! Local variables:
    INTEGER                                            :: k
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  d_2D
    INTEGER                                            :: wd_2D
    
    ! Allocate temporary shared memory for the extended and smoothed data fields
    CALL allocate_shared_dp_2D( grid%ny, grid%nx, d_2D, wd_2D)
    
    DO k = 1, nz
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
    ! NOTE: since IMAU-ICE doesn't consider subglacial lakes, Vostok simply shows
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
  
! == Some operations on sparse matrices in CSR format
  SUBROUTINE multiply_matrix_matrix_CSR( AA, BB, CC)
    ! Perform the matrix multiplication C = A*B, where all three
    ! matrices are provided in CSR format (parallelised)
    !
    ! NOTE: probably not very efficient, better to use PETSc in the future!
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA, BB
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: CC
    
    ! Local variables:
    INTEGER                                            :: nnz_max
    TYPE(type_sparse_matrix_CSR)                       :: BBT
    INTEGER                                            :: i1, i2, ic, jc
    INTEGER                                            :: ia,  ka1,  ka2,  ka,  ja1,  ja2,  ja
    INTEGER                                            :: ibt, kbt1, kbt2, kbt, jbt1, jbt2, jbt
    LOGICAL                                            :: is_nonzero
    REAL(dp)                                           :: Cij
    
    ! Safety
    IF (AA%n /= BB%m) THEN
      IF (par%master) WRITE(0,*) 'multiply_matrix_matrix_CSR - ERROR: sizes of A and B dont match!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared distributed memory for CC
    nnz_max = C%nconmax * (AA%nnz + BB%nnz)
    CALL allocate_matrix_CSR_dist( CC, AA%m, BB%n, nnz_max)
    
    ! Calculate the transpose BT of B
    CALL transpose_matrix_CSR( BB, BBT)
    
    ! Partition rows of over the processors
    CALL partition_list( CC%m, par%i, par%n, i1, i2)
    
    CC%nnz = 0
    CC%ptr = 1
    DO ic = i1, i2
    
      DO jc = 1, CC%n
        
        ! Calculate dot product between the i-th row of A and the j-th row of BT
        ia   = ic
        ka1  = AA%ptr( ia)
        ka2  = AA%ptr( ia+1) - 1
        IF (ka2 < ka1) CYCLE ! This row of A has no entries
        ja1  = AA%index( ka1)
        ja2  = AA%index( ka2)
        
        ibt  = jc
        kbt1 = BBT%ptr( ibt)
        kbt2 = BBT%ptr( ibt+1) - 1
        IF (kbt2 < kbt1) CYCLE ! This row of BT has no entries
        jbt1 = BBT%index( kbt1)
        jbt2 = BBT%index( kbt2)
        
        ! Compute the dot product
        Cij = 0._dp
        is_nonzero = .FALSE.
        DO ka = ka1, ka2
          ja = AA%index( ka)
          DO kbt = kbt1, kbt2
            jbt = BBT%index( kbt)
            IF (jbt == ja) THEN
              is_nonzero = .TRUE.
              Cij = Cij + (AA%val( ka) * BBT%val( kbt))
              EXIT
            END IF
          END DO
        END DO
        
        ! If the dot product is non-zero, list it in C
        IF (is_nonzero) THEN
          CC%nnz  = CC%nnz + 1
          CC%index( CC%nnz) = jc
          CC%val(   CC%nnz) = Cij
        END IF
        
      END DO ! DO jc = 1, CC%n
      
      ! Finalise this row of C
      CC%ptr( ic+1) = CC%nnz + 1
      
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( CC, i1, i2)
    
    ! Deallocate BT
    CALL deallocate_matrix_CSR( BBT)
    
  END SUBROUTINE multiply_matrix_matrix_CSR
  SUBROUTINE multiply_matrix_vector_CSR( AA, BB, CC)
    ! Perform the matrix multiplication C = A*B, where 
    ! A is a CSR-format matrix, and B and C are regular vectors.
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: BB
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: CC
    
    ! Local variables:
    INTEGER                                            :: mb, mc, i1, i2, i, k, j
    
    mb = SIZE( BB,1)
    mc = SIZE( CC,1)
    
    ! Safety
    IF (mb /= AA%n .OR. mc /= AA%m) THEN
      IF (par%master) WRITE(0,*) 'multiply_matrix_vector_CSR - ERROR: sizes of A, B, and C dont match!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Partition rows over the processors
    CALL partition_list( mc, par%i, par%n, i1, i2)
    
    DO i = i1, i2
      CC( i) = 0._dp
      DO k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%index( k)
        CC( i) = CC( i) + AA%val( k) * BB( j)
      END DO
    END DO
    CALL sync
    
  END SUBROUTINE multiply_matrix_vector_CSR
  SUBROUTINE multiply_matrix_vector_2D_CSR( AA, BB, CC)
    ! Perform the matrix multiplication C = A*B, where 
    ! A is a CSR-format matrix, and B and C are regular vectors.
    ! 
    ! NOTE: A = [m-by-n], B,C = [n-by-nz]
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: BB
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: CC
    
    ! Local variables:
    INTEGER                                            :: mb, mc, i1, i2, i, k, j
    
    mb = SIZE( BB,1)
    mc = SIZE( CC,1)
    
    ! Safety
    IF (mb /= AA%n .OR. mc /= AA%m) THEN
      IF (par%master) WRITE(0,*) 'multiply_matrix_vector_CSR - ERROR: sizes of A, B, and C dont match!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Partition rows over the processors
    CALL partition_list( mc, par%i, par%n, i1, i2)
    
    DO i = i1, i2
      CC( i,:) = 0._dp
      DO k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%index( k)
        CC( i,:) = CC( i,:) + AA%val( k) * BB( j,:)
      END DO
    END DO
    CALL sync
    
  END SUBROUTINE multiply_matrix_vector_2D_CSR
  SUBROUTINE add_matrix_matrix_CSR( AA, BB, CC, alpha, beta)
    ! Perform the matrix multiplication C = (alpha*A)+(beta*B), where all three
    ! matrices are provided in CSR format (parallelised)
    ! 
    ! NOTE: assumes column entries and A and B are sorted ascending!
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA, BB
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: CC
    REAL(dp), OPTIONAL,                  INTENT(IN)    :: alpha, beta
    
    ! Local variables:
    REAL(dp)                                           :: alphap, betap
    INTEGER                                            :: i1, i2, ic
    INTEGER                                            :: ka1, ka2, nnz_row_a
    INTEGER                                            :: kb1, kb2, nnz_row_b
    INTEGER                                            :: ka, kb, ja, jb
    LOGICAL                                            :: finished_a, finished_b
    
    ! Safety
    IF (AA%m /= BB%m .OR. AA%n /= BB%n) THEN
      IF (par%master) WRITE(0,*) 'add_matrix_matrix_CSR - ERROR: A and B are not of the same size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! If no coefficients are provided, assume unity
    IF (PRESENT( alpha)) THEN
      alphap = alpha
    ELSE
      alphap = 1._dp
    END IF
    IF (PRESENT( beta )) THEN
      betap  = beta
    ELSE
      betap  = 1._dp
    END IF
    
    ! Partition rows over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    ! Allocate distributed shared memory for C
    CALL allocate_matrix_CSR_dist( CC, AA%m, AA%n, AA%nnz + BB%nnz)
    
    ! Initialise
    CC%ptr = 1
    CC%nnz = 0
    
    DO ic = i1, i2
      
      ka1 = AA%ptr( ic)
      ka2 = AA%ptr( ic+1) - 1
      nnz_row_a = ka2 + 1 - ka1
      
      kb1 = BB%ptr( ic)
      kb2 = BB%ptr( ic+1) - 1
      nnz_row_b = kb2 + 1 - kb1
      
      IF (nnz_row_a == 0 .AND. nnz_row_b == 0) THEN
        ! Neither A nor B has entries in this row
        CYCLE
      ELSEIF (nnz_row_a == 0 .AND. nnz_row_b > 0) THEN
        ! A has no entries in this row, but B does; copy data from B
        CC%index( CC%nnz+1: CC%nnz+nnz_row_b) = BB%index( kb1:kb2)
        CC%val(   CC%nnz+1: CC%nnz+nnz_row_b) = BB%val(   kb1:kb2)
        CC%nnz = CC%nnz + nnz_row_b
      ELSEIF (nnz_row_a > 0 .AND. nnz_row_b == 0) THEN
        ! B has no entries in this row, but A does; copy data from A
        CC%index( CC%nnz+1: CC%nnz+nnz_row_a) = AA%index( ka1:ka2)
        CC%val(   CC%nnz+1: CC%nnz+nnz_row_a) = AA%val(   ka1:ka2)
        CC%nnz = CC%nnz + nnz_row_a
      ELSE
        ! Both A and B have entries in this row
        
        ka = ka1
        kb = kb1
        finished_a = .FALSE.
        finished_b = .FALSE.
        
        DO WHILE ((.NOT. finished_a) .OR. (.NOT. finished_b))
          
          IF ((.NOT. finished_a) .AND. (.NOT. finished_b)) THEN
          
            ja = AA%index( ka)
            jb = BB%index( kb)
            
            IF (ja < jb) THEN
              CC%nnz =  CC%nnz + 1
              CC%index( CC%nnz) = AA%index( ka)
              CC%val(   CC%nnz) = alphap * AA%val( ka)
              ka = ka + 1
            ELSEIF (jb < ja) THEN
              CC%nnz =  CC%nnz + 1
              CC%index( CC%nnz) = BB%index( kb)
              CC%val(   CC%nnz) = betap  * BB%val( kb)
              kb = kb + 1
            ELSEIF (jb == ja) THEN
              CC%nnz =  CC%nnz + 1
              CC%index( CC%nnz) = AA%index( ka)
              CC%val(   CC%nnz) = alphap * AA%val( ka) + betap * BB%val( kb)
              ka = ka + 1
              kb = kb + 1
            END IF
          
          ELSEIF (.NOT. finished_a) THEN
          
            CC%nnz =  CC%nnz + 1
            CC%index( CC%nnz) = AA%index( ka)
            CC%val(   CC%nnz) = alphap * AA%val( ka)
            ka = ka + 1
            
          ELSEIF (.NOT. finished_b) THEN
          
            CC%nnz =  CC%nnz + 1
            CC%index( CC%nnz) = BB%index( kb)
            CC%val(   CC%nnz) = betap * BB%val( kb)
            kb = kb + 1
            
          END IF ! IF ((.NOT. finished_a) .AND. (.NOT. finished_b)) THEN
          
          IF (ka == ka2 + 1) finished_a = .TRUE.
          IF (kb == kb2 + 1) finished_b = .TRUE.
          
        END DO ! DO WHILE ((.NOT. finished_a) .OR. (.NOT. finished_b))
        
      END IF ! IF (nnz_row_a == 0 .AND. nnz_row_b == 0) THEN
      
      ! Finalise this row of C
      CC%ptr( ic+1:CC%m+1) = CC%nnz+1
      
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( CC, i1, i2)
    
  END SUBROUTINE add_matrix_matrix_CSR
  SUBROUTINE overwrite_rows_CSR( AA, BB)
    ! Given sparse matrices B and A (both represented in CSR format),
    ! whenever a row in B has a non-zero entry, overwrite that row in A with that row in B
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: AA
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: BB
    
    ! Local variables:
    TYPE(type_sparse_matrix_CSR)                       :: AA_old
    INTEGER                                            :: i1, i2, k1, k2, i
    INTEGER                                            :: ka1, ka2, ka, ja
    INTEGER                                            :: kb1, kb2, kb, jb
    
    ! Safety
    IF (AA%m /= BB%m .OR. AA%n /= BB%n) THEN
      IF (par%master) WRITE(0,*) 'overwrite_rows_CSR - ERROR: A and B are not of the same size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Partition rows over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    ! Allocate temporary shared memory for old matrix A
    CALL allocate_matrix_CSR_shared( AA_old, AA%m, AA%n, AA%nnz)
    
    ! Copy old A data to temporary memory
    CALL partition_list( AA%nnz, par%i, par%n, k1, k2)
    IF (par%master) THEN
      AA_old%m       = AA%m
      AA_old%n       = AA%n
      AA_old%nnz     = AA%nnz
      AA_old%nnz_max = AA%nnz_max
      AA_old%ptr( AA%m+1) = AA%ptr( AA%m+1)
    END IF
    CALL sync
    AA_old%ptr(   i1:i2) = AA%ptr(   i1:i2)
    AA_old%index( k1:k2) = AA%index( k1:k2)
    AA_old%val(   k1:k2) = AA%val(   k1:k2)
    CALL sync
    
    ! Deallocate A
    CALL deallocate_matrix_CSR( AA)
    
    ! Allocate distributed shared memory for new A
    CALL allocate_matrix_CSR_dist( AA, AA_old%m, AA_old%n, AA_old%nnz + BB%nnz)
    
    ! Initialise
    AA%ptr = 1
    AA%nnz = 0
    
    DO i = i1, i2
    
      ! Check if this row in B has any non-zero entries
      kb1 = BB%ptr( i)
      kb2 = BB%ptr( i+1) - 1
      IF (kb2 >= kb1) THEN
        ! This row in B has non-zero entries; write it to A
        
        DO kb = kb1, kb2
          jb = BB%index( kb)
          AA%nnz  = AA%nnz + 1
          AA%index( AA%nnz) = jb
          AA%val(   AA%nnz) = BB%val( kb)
        END DO
        
      ELSE ! IF (kb1 > kb2) THEN
        ! This row in B has no non-zero entries, use data from A_old
        
        ka1 = AA_old%ptr( i)
        ka2 = AA_old%ptr( i+1) - 1
        DO ka = ka1, ka2
          ja = AA_old%index( ka)
          AA%nnz  = AA%nnz + 1
          AA%index( AA%nnz) = ja
          AA%val(   AA%nnz) = AA_old%val( ka)
        END DO
        
      END IF
      
      ! Finalise this row of A
      AA%ptr( i+1) = AA%nnz+1
      
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Deallocate A_old
    CALL deallocate_matrix_CSR( AA_old)
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( AA, i1, i2)
    
  END SUBROUTINE overwrite_rows_CSR
  SUBROUTINE transpose_matrix_CSR( AA, AAT)
    ! Calculate the transpose AT of the matrix A, where both
    ! matrices are provided in CSR format (parallelised)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: AAT
    
    ! Local variables:
    INTEGER                                            :: i1, i2, i, ii, kk1, kk2, kk, jj
    
    ! Allocate shared distributed memory for AAT
    CALL allocate_matrix_CSR_dist( AAT, AA%n, AA%m, AA%nnz)
    
    ! Partition rows of AT (so columns of A) over the processors
    CALL partition_list( AAT%m, par%i, par%n, i1, i2)
    
    AAT%nnz = 0
    AAT%ptr = 1
    DO i = i1, i2
      
      ! Find all entries in column i of A
      DO ii = 1, AA%m ! ii loops over all the rows of A
        kk1 = AA%ptr( ii)
        kk2 = AA%ptr( ii+1) - 1
        DO kk = kk1, kk2 ! kk loops over all the entries of A(ii,:)
          jj = AA%index( kk)
          IF (jj == i) THEN
            ! Found an entry in the i-th column, ii-th row of A; add it in the i-th row, ii-th column of AT
            AAT%nnz = AAT%nnz + 1
            AAT%index( AAT%nnz) = ii
            AAT%val(   AAT%nnz) = AA%val( kk)
          END IF
        END DO
      END DO
      
      ! Finalise this row of AT
      AAT%ptr( i+1) = AAT%nnz + 1
      
    END DO ! DO i = i1, i2
    CALL sync
    
    ! Combine results from the different processes
    CALL finalise_matrix_CSR_dist( AAT, i1, i2)
    
  END SUBROUTINE transpose_matrix_CSR
  SUBROUTINE multiply_matrix_rows_with_vector( AA, BB, CC)
    ! Multiply the rows of the m-by-n matrix A with the elements
    ! of the m-by-1 vector B, call the result C
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: BB
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: CC
    
    ! Local variables:
    INTEGER                                            :: i1,i2,i,k
    
    ! Safety
    IF (SIZE(BB,1) /= AA%m) THEN
      IF (par%master) WRITE(0,*) 'multiply_matrix_rows_with_vector - ERROR: sizes of A and B dont match!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Allocate shared memory for C
    CALL allocate_matrix_CSR_shared( CC, AA%m, AA%n, AA%nnz)
    
    ! Partition rows over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    IF (par%master) THEN
      CC%m       = AA%m
      CC%n       = AA%n
      CC%nnz     = AA%nnz
      CC%nnz_max = AA%nnz_max
      CC%ptr( CC%m+1) = AA%ptr( AA%m+1)
    END IF
    CALL sync
    
    DO i = i1, i2
      CC%ptr( i) = AA%ptr( i)
      DO k = AA%ptr( i), AA%ptr( i+1)-1
        CC%index( k) = AA%index( k)
        CC%val(   k) = AA%val(   k) * BB( i)
      END DO
    END DO
    CALL sync
    
  END SUBROUTINE multiply_matrix_rows_with_vector
  SUBROUTINE convert_vector_to_diag_CSR( A_vec, AA)
    ! Convert the vector A_vec to the diagonal matrix A
      
    IMPLICIT NONE
    
    ! In- and output variables:
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: A_vec
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: AA
    
    ! Local variables:
    INTEGER                                            :: m, i, i1, i2
    
    m = SIZE( A_vec, 1)
    
    ! Allocate shared memory for A
    CALL allocate_matrix_CSR_shared( AA, m, m, m)
    IF (par%master) THEN
      AA%m       = m
      AA%n       = m
      AA%nnz_max = m
      AA%nnz     = m
    END IF
    CALL sync
    
    ! Partition rows over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    ! Fill in values
    DO i = i1, i2
      AA%ptr(   i) = i
      AA%index( i) = i
      AA%val(   i) = A_vec( i)
    END DO
    IF (par%master) AA%ptr( m+1) = m+1
    CALL sync
    
  END SUBROUTINE convert_vector_to_diag_CSR
  SUBROUTINE are_identical_matrices_CSR( AA, BB, isso)
    ! Check if CSR-formatted matrices A and B are identical
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA, BB
    LOGICAL,                             INTENT(OUT)   :: isso
    
    ! Local variables:
    INTEGER                                            :: i1, i2, k1, k2, i, k
    
    isso = .TRUE.
    
    ! Simple dimension check
    IF (AA%m /= BB%m .OR. AA%n /= BB%n .OR. AA%nnz /= BB%nnz) THEN
      isso = .FALSE.
      RETURN
    END IF
    
    ! ptr
    CALL partition_list( AA%m  , par%i, par%n, i1, i2)
    DO i = i1, i2
      IF (AA%ptr( i) /= BB%ptr( i)) THEN
        isso = .FALSE.
        EXIT
      END IF
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, isso, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. isso) RETURN
    
    ! index, val
    CALL partition_list( AA%nnz, par%i, par%n, k1, k2)
    DO k = k1, k2
      IF (AA%index( k) /= BB%index( k) .OR. AA%val( k) /= BB%val( k)) THEN
        isso = .FALSE.
        EXIT
      END IF
    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, isso, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    IF (.NOT. isso) RETURN
    
  END SUBROUTINE are_identical_matrices_CSR
  
  SUBROUTINE solve_matrix_equation_CSR( AA, b, x, choice_matrix_solver, SOR_nit, SOR_tol, SOR_omega, PETSc_rtol, PETSc_abstol, colour_v1, colour_v2, colour_vi)
    ! Solve the matrix equation Ax = b
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: x
    CHARACTER(LEN=256),                  INTENT(IN)    :: choice_matrix_solver
    INTEGER,                             INTENT(IN)    :: SOR_nit
    REAL(dp),                            INTENT(IN)    :: SOR_tol
    REAL(dp),                            INTENT(IN)    :: SOR_omega
    REAL(dp),                            INTENT(IN)    :: PETSc_rtol
    REAL(dp),                            INTENT(IN)    :: PETSc_abstol
    INTEGER,  DIMENSION(5)    , OPTIONAL, INTENT(IN)   :: colour_v1, colour_v2
    INTEGER,  DIMENSION(:,:  ), OPTIONAL, INTENT(IN)   :: colour_vi
    
    ! Safety
    IF (AA%m /= AA%n .OR. AA%m /= SIZE( b,1) .OR. AA%m /= SIZE( x,1)) THEN
      IF (par%master) WRITE(0,*) 'solve_matrix_equation_CSR - ERROR: matrix sizes dont match!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    CALL check_for_NaN_dp_1D(  AA%val, 'AA%val', 'solve_matrix_equation_CSR')
    CALL check_for_NaN_dp_1D(  b,      'b',      'solve_matrix_equation_CSR')
    CALL check_for_NaN_dp_1D(  x,      'x',      'solve_matrix_equation_CSR')
    
    IF (choice_matrix_solver == 'SOR') THEN
      ! Use the old simple SOR solver
      
      CALL solve_matrix_equation_CSR_SOR( AA, b, x, SOR_nit, SOR_tol, SOR_omega, colour_v1, colour_v2, colour_vi)
      
    ELSEIF (choice_matrix_solver == 'PETSc') THEN
      ! Use the PETSc solver (much preferred, this is way faster and more stable!)
      
      IF (par%master) WRITE(0,*) 'solve_matrix_equation_CSR - ERROR: PETSc is not yet available in UFEMISM!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    
      !CALL solve_matrix_equation_CSR_PETSc( AA, b, x, PETSc_rtol, PETSc_abstol)
    
    ELSE
      IF (par%master) WRITE(0,*) 'solve_matrix_equation_CSR - ERROR: unknown choice_matrix_solver "', choice_matrix_solver, '"!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
  END SUBROUTINE solve_matrix_equation_CSR
  SUBROUTINE solve_matrix_equation_CSR_SOR( AA, b, x, nit, tol, omega, colour_v1, colour_v2, colour_vi)
    ! Solve the matrix equation Ax = b using successive over-relaxation (SOR)
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: x
    INTEGER,                             INTENT(IN)    :: nit
    REAL(dp),                            INTENT(IN)    :: tol
    REAL(dp),                            INTENT(IN)    :: omega
    INTEGER,  DIMENSION(5)    , OPTIONAL, INTENT(IN)   :: colour_v1, colour_v2
    INTEGER,  DIMENSION(:,:  ), OPTIONAL, INTENT(IN)   :: colour_vi
    
    IF (PRESENT( colour_v1)) THEN
      ! Safety
      IF ((.NOT. PRESENT( colour_v2)) .OR. (.NOT. PRESENT( colour_vi))) THEN
        IF (par%master) WRITE(0,*) 'solve_matrix_equation_CSR_SOR - ERROR: needs all three colour arguments!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      CALL solve_matrix_equation_CSR_SOR_coloured( AA, b, x, nit, tol, omega, colour_v1, colour_v2, colour_vi)
    ELSE
      CALL solve_matrix_equation_CSR_SOR_regular(  AA, b, x, nit, tol, omega)
    END IF
    
  END SUBROUTINE solve_matrix_equation_CSR_SOR
  SUBROUTINE solve_matrix_equation_CSR_SOR_regular( AA, b, x, nit, tol, omega)
    ! Solve the matrix equation Ax = b using successive over-relaxation (SOR)
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: x
    INTEGER,                             INTENT(IN)    :: nit
    REAL(dp),                            INTENT(IN)    :: tol
    REAL(dp),                            INTENT(IN)    :: omega
    
    ! Local variables:
    INTEGER                                            :: i,j,k,it,i1,i2
    REAL(dp)                                           :: lhs, res, res_max, omega_dyn
    REAL(dp), DIMENSION(:    ), POINTER                ::  diagA
    INTEGER                                            :: wdiagA
    LOGICAL                                            :: found_it
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( AA%m, diagA, wdiagA)
    
    ! Partition equations over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    ! Store the central diagonal separately for faster access
    DO i = i1, i2
      found_it = .FALSE.
      DO k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%index( k)
        IF (j == i) THEN
          diagA( i) = AA%val( k)
          found_it = .TRUE.
        END IF
      END DO
      IF (.NOT. found_it) THEN
        WRITE(0,*) 'solve_matrix_equation_CSR_SOR - ERROR: matrix is missing a diagonal element!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END DO
    CALL sync
    
    ! Perform the successive over-relaxation
    omega_dyn = omega
    
    res_max = tol * 2._dp
    it = 0
    SOR_iterate: DO WHILE (res_max > tol .AND. it < nit)
      it = it+1
      
      res_max = 0._dp

      DO i = i1, i2
      
        lhs = 0._dp
        DO k = AA%ptr( i), AA%ptr( i+1)-1
          j = AA%index( k)
          lhs = lhs + AA%val( k) * x( j)
        END DO
        
        res = (lhs - b( i)) / diagA( i)
        res_max = MAX( res_max, res)
        
        x( i) = x( i) - omega_dyn * res
        
      END DO ! DO i = i1, i2
      CALL sync
      
      ! Check if we've reached a stable solution
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, res_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      
      !IF (par%master) WRITE(0,*) '      SOR iteration ', it, ': res_max = ', res_max
      
      IF (it > 100 .AND. res_max > 1E3_dp ) THEN
        
        ! Divergence detected - decrease omega, reset solution to zero, restart SOR.
        IF (par%master) WRITE(0,*) '  solve_matrix_equation_CSR_SOR - divergence detected; decrease omega, reset solution to zero, restart SOR'
        omega_dyn = omega_dyn - 0.1_dp
        it = 0
        x( i1:i2) = 0._dp
        CALL sync
        
        IF (omega_dyn <= 0.1_dp) THEN
          IF (par%master) WRITE(0,*) '  solve_matrix_equation_CSR_SOR - ERROR: divergence detected even with extremely low relaxation parameter!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END IF
      
    END DO SOR_iterate
    
    ! Clean up after yourself
    CALL deallocate_shared( wdiagA)
    
  END SUBROUTINE solve_matrix_equation_CSR_SOR_regular
  SUBROUTINE solve_matrix_equation_CSR_SOR_coloured( AA, b, x, nit, tol, omega, colour_v1, colour_v2, colour_vi)
    ! Solve the matrix equation Ax = b using successive over-relaxation (SOR)
    ! The matrix A is provided in Compressed Sparse Row (CSR) format
    ! 
    ! Use a five-colouring map for perfect parallelisation
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: b
    REAL(dp), DIMENSION(:    ),          INTENT(OUT)   :: x
    INTEGER,                             INTENT(IN)    :: nit
    REAL(dp),                            INTENT(IN)    :: tol
    REAL(dp),                            INTENT(IN)    :: omega
    INTEGER,  DIMENSION(5)    ,          INTENT(IN)   :: colour_v1, colour_v2
    INTEGER,  DIMENSION(:,:  ),          INTENT(IN)   :: colour_vi
    
    ! Local variables:
    INTEGER                                            :: fci,fcvi,i,j,k,it,i1,i2
    REAL(dp)                                           :: lhs, res, res_max, omega_dyn
    REAL(dp), DIMENSION(:    ), POINTER                ::  diagA
    INTEGER                                            :: wdiagA
    LOGICAL                                            :: found_it
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( AA%m, diagA, wdiagA)
    
    ! Partition equations over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    ! Store the central diagonal separately for faster access
    DO i = i1, i2
      found_it = .FALSE.
      DO k = AA%ptr( i), AA%ptr( i+1)-1
        j = AA%index( k)
        IF (j == i) THEN
          diagA( i) = AA%val( k)
          found_it = .TRUE.
        END IF
      END DO
      IF (.NOT. found_it) THEN
        WRITE(0,*) 'solve_matrix_equation_CSR_SOR - ERROR: matrix is missing a diagonal element!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END DO
    CALL sync
    
    ! Perform the successive over-relaxation
    omega_dyn = omega
    
    res_max = tol * 2._dp
    it = 0
    SOR_iterate: DO WHILE (res_max > tol .AND. it < nit)
      it = it+1
      
      res_max = 0._dp

      DO fci = 1, 5
        DO fcvi = colour_v1( fci), colour_v2( fci)
        
          i = colour_vi( fcvi, fci)
      
          lhs = 0._dp
          DO k = AA%ptr( i), AA%ptr( i+1)-1
            j = AA%index( k)
            lhs = lhs + AA%val( k) * x( j)
          END DO
          
          res = (lhs - b( i)) / diagA( i)
          res_max = MAX( res_max, res)
          
          x( i) = x( i) - omega_dyn * res
        
        END DO ! DO fcvi = fcvi1, fcvi2
        CALL sync
      END DO ! DO fci = 1, 5
      
      ! Check if we've reached a stable solution
      CALL MPI_ALLREDUCE( MPI_IN_PLACE, res_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      
      !IF (par%master) WRITE(0,*) '      SOR iteration ', it, ': res_max = ', res_max
      
      IF (it > 100 .AND. res_max > 1E3_dp ) THEN
        
        ! Divergence detected - decrease omega, reset solution to zero, restart SOR.
        IF (par%master) WRITE(0,*) '  solve_matrix_equation_CSR_SOR - divergence detected; decrease omega, reset solution to zero, restart SOR'
        omega_dyn = omega_dyn - 0.1_dp
        it = 0
        x( i1:i2) = 0._dp
        CALL sync
        
        IF (omega_dyn <= 0.1_dp) THEN
          IF (par%master) WRITE(0,*) '  solve_matrix_equation_CSR_SOR - ERROR: divergence detected even with extremely low relaxation parameter!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      END IF
      
    END DO SOR_iterate
    
    ! Clean up after yourself
    CALL deallocate_shared( wdiagA)
    
  END SUBROUTINE solve_matrix_equation_CSR_SOR_coloured
  
  SUBROUTINE allocate_matrix_CSR_shared( A, m, n, nnz_max)
    ! Allocate shared memory for a CSR-format sparse m-by-n matrix A
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    INTEGER,                             INTENT(IN)    :: m, n, nnz_max
    
    CALL allocate_shared_int_0D( A%m,       A%wm      )
    CALL allocate_shared_int_0D( A%n,       A%wn      )
    CALL allocate_shared_int_0D( A%nnz_max, A%wnnz_max)
    CALL allocate_shared_int_0D( A%nnz,     A%wnnz    )
    
    IF (par%master) THEN
      A%m       = m
      A%n       = n
      A%nnz_max = nnz_max
      A%nnz     = 0
    END IF
    CALL sync
    
    CALL allocate_shared_int_1D( A%m+1,     A%ptr,   A%wptr  )
    CALL allocate_shared_int_1D( A%nnz_max, A%index, A%windex)
    CALL allocate_shared_dp_1D(  A%nnz_max, A%val,   A%wval  )
    
  END SUBROUTINE allocate_matrix_CSR_shared
  SUBROUTINE allocate_matrix_CSR_dist( A, m, n, nnz_max_loc)
    ! Allocate shared memory for a CSR-format sparse m-by-n matrix A
    ! 
    ! NOTE: uses distributed shared memory, so the processes can create their own
    !       CSR lists in parallel; use "finalise_matrix_CSR" after this is done!
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    INTEGER,                             INTENT(IN)    :: m, n, nnz_max_loc
    
    CALL allocate_shared_dist_int_0D( A%m,       A%wm      )
    CALL allocate_shared_dist_int_0D( A%n,       A%wn      )
    CALL allocate_shared_dist_int_0D( A%nnz_max, A%wnnz_max)
    CALL allocate_shared_dist_int_0D( A%nnz,     A%wnnz    )
    
    A%m       = m
    A%n       = n
    A%nnz_max = nnz_max_loc
    A%nnz     = 0
    
    CALL allocate_shared_dist_int_1D( A%m+1,     A%ptr,   A%wptr  )
    CALL allocate_shared_dist_int_1D( A%nnz_max, A%index, A%windex)
    CALL allocate_shared_dist_dp_1D(  A%nnz_max, A%val,   A%wval  )
    
  END SUBROUTINE allocate_matrix_CSR_dist
  SUBROUTINE extend_matrix_CSR_dist( A, nnz_max_loc_new)
    ! Extend shared memory for a CSR-format sparse m-by-n matrix A
    ! 
    ! NOTE: uses distributed shared memory, so the processes can create their own
    !       CSR lists in parallel; use "finalise_matrix_CSR_dist" after this is done!
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    INTEGER,                             INTENT(IN)    :: nnz_max_loc_new
    
    A%nnz_max = nnz_max_loc_new
    
    CALL adapt_shared_dist_int_1D( A%nnz, A%nnz_max, A%index, A%windex)
    CALL adapt_shared_dist_dp_1D(  A%nnz, A%nnz_max, A%val,   A%wval  )
  
  END SUBROUTINE extend_matrix_CSR_dist
  SUBROUTINE sort_columns_in_CSR_dist( A)
    ! Sort the columns in each row of CSR-formatted matrix A in ascending order
    ! 
    ! NOTE: A is still stored distributed over the processes
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    
    ! Local variables:
    INTEGER                                            :: i,k1,k2,nnz_row,k,kk,jk,jkk
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: index_row
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: val_row
    
    DO i = 1, A%m
      
      k1 = A%ptr( i)
      k2 = A%ptr( i+1) - 1
      nnz_row = k2 + 1 - k1
      
      IF (nnz_row <= 0) CYCLE
      
      ! Allocate temporary memory for this row
      ALLOCATE( index_row( nnz_row))
      ALLOCATE( val_row(   nnz_row))
      
      ! Copy data for this row to temporary memory
      index_row = A%index( k1:k2)
      val_row   = A%val(   k1:k2)
      
      ! Sort the data
      DO k = 1, nnz_row
        jk = index_row( k)
        DO kk = k+1, nnz_row
          jkk = index_row( kk)
          IF (jkk < jk) THEN
            ! Switch columns
            index_row( kk) = index_row( kk) + index_row( k)
            index_row( k ) = index_row( kk) - index_row( k)
            index_row( kk) = index_row( kk) - index_row( k)
            val_row(   kk) = val_row(   kk) + val_row(   k)
            val_row(   k ) = val_row(   kk) - val_row(   k)
            val_row(   kk) = val_row(   kk) - val_row(   k)
          END IF
        END DO
      END DO ! DO k = 1, nnz_row
      
      ! Copy sorted data back
      A%index( k1:k2) = index_row
      A%val(   k1:k2) = val_row
      
      ! Clean up after yourself
      DEALLOCATE( index_row)
      DEALLOCATE( val_row  )
      
    END DO
    CALL sync
    
  END SUBROUTINE sort_columns_in_CSR_dist
  SUBROUTINE finalise_matrix_CSR_dist( A, i1, i2)
    ! Finalise a CSR-format sparse m-by-n matrix A from distributed to shared memory
    !
    ! NOTE: each process has data for rows i1-i2
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    INTEGER,                             INTENT(IN)    :: i1, i2
    
    ! Local variables:
    INTEGER                                            :: status(MPI_STATUS_SIZE)
    INTEGER                                            :: nnz_tot
    INTEGER,  DIMENSION(:    ), POINTER                :: ptr, index
    REAL(dp), DIMENSION(:    ), POINTER                :: val
    INTEGER                                            :: wptr, windex, wval
    INTEGER                                            :: p, k1, k2, k_loc, k_sum
    INTEGER                                            :: m, n
    
    ! Determine total number of non-zero elements in A
    CALL MPI_ALLREDUCE( A%nnz, nnz_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    m = A%m
    n = A%n
    
    ! Allocate temporary shared memory
    CALL allocate_shared_int_1D( A%m+1,   ptr,   wptr  )
    CALL allocate_shared_int_1D( nnz_tot, index, windex)
    CALL allocate_shared_dp_1D(  nnz_tot, val,   wval  )
    
    ! Determine range of indices for each process
    k1 = 0
    k2 = 0
    IF (par%master) THEN
      k1    = 1
      k2    = A%nnz
      k_sum = A%nnz
    END IF
    DO p = 1, par%n-1
      IF (par%master) THEN
        CALL MPI_SEND( k_sum, 1, MPI_INTEGER, p, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_RECV( k_loc, 1, MPI_INTEGER, p, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        k_sum = k_sum + k_loc
      ELSEIF (p == par%i) THEN
        CALL MPI_RECV( k_sum, 1, MPI_INTEGER, 0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
        CALL MPI_SEND( A%nnz, 1, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
        k1 = k_sum + 1
        k2 = k_sum + A%nnz
        A%ptr = A%ptr + k_sum
      END IF
      CALL sync
    END DO
    
    ! Copy local data to temporary shared memory
    ptr(   i1:i2) = A%ptr(   i1:i2  )
    index( k1:k2) = A%index( 1:A%nnz)
    val(   k1:k2) = A%val(   1:A%nnz)
    CALL sync
    
    ! Deallocate distributed shared memory for A
    CALL deallocate_matrix_CSR( A)
    
    ! Allocate shared memory for A
    CALL allocate_matrix_CSR_shared( A, m, n, nnz_tot)
    
    ! Copy data back from temporary memory
    A%ptr(   i1:i2) = ptr(   i1:i2)
    A%index( k1:k2) = index( k1:k2)
    A%val(   k1:k2) = val(   k1:k2)
    IF (par%master) A%m   = m
    IF (par%master) A%n   = n
    IF (par%master) A%nnz = nnz_tot
    IF (par%master) A%ptr( A%m+1) = A%nnz+1
    CALL sync
    
    ! Deallocate temporary memory
    CALL deallocate_shared( wptr  )
    CALL deallocate_shared( windex)
    CALL deallocate_shared( wval  )
  
  END SUBROUTINE finalise_matrix_CSR_dist
  SUBROUTINE deallocate_matrix_CSR( A)
    ! Deallocate the memory used by the CSR-format sparse m-by-n matrix A
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: A
    
    CALL deallocate_shared( A%wm)
    CALL deallocate_shared( A%wn)
    CALL deallocate_shared( A%wnnz_max)
    CALL deallocate_shared( A%wnnz)
    CALL deallocate_shared( A%wptr)
    CALL deallocate_shared( A%windex)
    CALL deallocate_shared( A%wval)
    
    NULLIFY( A%m      )
    NULLIFY( A%n      )
    NULLIFY( A%nnz_max)
    NULLIFY( A%nnz    )
    NULLIFY( A%ptr    )
    NULLIFY( A%index  )
    NULLIFY( A%val    )
  
  END SUBROUTINE deallocate_matrix_CSR
  SUBROUTINE check_CSR_for_double_entries( A)
    ! Check a CSR matrix representation for double entries
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: A
    
    ! Local variables:
    INTEGER                                            :: i,j,k,i1,i2,k2,j2
    
    ! Partition equations over the processors
    CALL partition_list( A%m, par%i, par%n, i1, i2)

    DO i = i1, i2  
    
      DO k = A%ptr( i), A%ptr( i+1)-1
        j = A%index( k)
        
        DO k2 = A%ptr( i), A%ptr( i+1)-1
          IF (k2 == k) CYCLE
          j2 = A%index( k2)
          IF (j2 == j) THEN
            WRITE(0,*) 'check_CSR_for_double_entries - ERROR: double entry detected!'
            CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
          END IF
        END DO
        
      END DO
      
    END DO ! DO 
    
  END SUBROUTINE check_CSR_for_double_entries
  
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

    ! External subroutines:      
    EXTERNAL DGTSV ! Lapack routine that solves tridiagonal systems (in double precision).

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
  
! == Debugging
  SUBROUTINE check_for_NaN_dp_1D(  d, d_name, routine_name)
    ! Check if NaN values occur in the 1-D dp data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:    ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name, routine_name
    
    ! Local variables:
    INTEGER                                                :: nx,i,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc, routine_name_loc
    
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
    IF (PRESENT( routine_name)) THEN
      routine_name_loc = TRIM(routine_name)
    ELSE
      routine_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF (d( i) /= d( i)) THEN
        WRITE(0,'(A,I4,A,A,A,A,A,A)') 'check_for_NaN_dp_1D - NaN detected at [', i, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d (i) > HUGE( d( i))) THEN
        WRITE(0,'(A,I4,A,A,A,A,A,A)') 'check_for_NaN_dp_1D - Inf detected at [', i, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d (i) < -HUGE( d( i))) THEN
        WRITE(0,'(A,I4,A,A,A,A,A,A)') 'check_for_NaN_dp_1D - -Inf detected at [', i, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_dp_1D
  SUBROUTINE check_for_NaN_dp_2D(  d, d_name, routine_name)
    ! Check if NaN values occur in the 2-D dp data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:,:  ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name, routine_name
    
    ! Local variables:
    INTEGER                                                :: nx,ny,i,j,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc, routine_name_loc
    
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
    IF (PRESENT( routine_name)) THEN
      routine_name_loc = TRIM(routine_name)
    ELSE
      routine_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF (d( j,i) /= d( j,i)) THEN
        WRITE(0,'(A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_dp_2D - NaN detected at [', i, ',', j, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( j,i) > HUGE( d( j,i))) THEN
        WRITE(0,'(A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_dp_2D - Inf detected at [', i, ',', j, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( j,i) < -HUGE( d( j,i))) THEN
        WRITE(0,'(A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_dp_2D - Inf detected at [', i, ',', j, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_dp_2D
  SUBROUTINE check_for_NaN_dp_3D(  d, d_name, routine_name)
    ! Check if NaN values occur in the 3-D dp data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    REAL(dp), DIMENSION(:,:,:),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name, routine_name
    
    ! Local variables:
    INTEGER                                                :: nx,ny,nz,i,j,k,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc, routine_name_loc
    
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
    IF (PRESENT( routine_name)) THEN
      routine_name_loc = TRIM(routine_name)
    ELSE
      routine_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    DO k = 1, nz
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF (d( k,j,i) /= d( k,j,i)) THEN
        WRITE(0,'(A,I4,A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_dp_3D - NaN detected at [', i, ',', j, ',', k, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( k,j,i) > HUGE( d( k,j,i))) THEN
        WRITE(0,'(A,I4,A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_dp_3D - Inf detected at [', i, ',', j, ',', k, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( k,j,i) < -HUGE( d( k,j,i))) THEN
        WRITE(0,'(A,I4,A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_dp_3D - -Inf detected at [', i, ',', j, ',', k, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_dp_3D
  SUBROUTINE check_for_NaN_int_1D( d, d_name, routine_name)
    ! Check if NaN values occur in the 1-D int data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,  DIMENSION(:    ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name, routine_name
    
    ! Local variables:
    INTEGER                                                :: nx,i,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc, routine_name_loc
    
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
    IF (PRESENT( routine_name)) THEN
      routine_name_loc = TRIM(routine_name)
    ELSE
      routine_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF (d( i) /= d( i)) THEN
        WRITE(0,'(A,I4,A,A,A,A,A,A)') 'check_for_NaN_int_1D - NaN detected at [', i, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( i) > HUGE( d( i))) THEN
        WRITE(0,'(A,I4,A,A,A,A,A,A)') 'check_for_NaN_int_1D - Inf detected at [', i, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( i) < -HUGE( d( i))) THEN
        WRITE(0,'(A,I4,A,A,A,A,A,A)') 'check_for_NaN_int_1D - -Inf detected at [', i, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_int_1D
  SUBROUTINE check_for_NaN_int_2D( d, d_name, routine_name)
    ! Check if NaN values occur in the 2-D int data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,  DIMENSION(:,:  ),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name, routine_name
    
    ! Local variables:
    INTEGER                                                :: nx,ny,i,j,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc, routine_name_loc
    
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
    IF (PRESENT( routine_name)) THEN
      routine_name_loc = TRIM(routine_name)
    ELSE
      routine_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF (d( j,i) /= d( j,i)) THEN
        WRITE(0,'(A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_int_2D - NaN detected at [', i, ',', j, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( j,i) > HUGE( d( j,i))) THEN
        WRITE(0,'(A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_int_2D - Inf detected at [', i, ',', j, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( j,i) < -HUGE( d( j,i))) THEN
        WRITE(0,'(A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_int_2D - -Inf detected at [', i, ',', j, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_int_2D
  SUBROUTINE check_for_NaN_int_3D( d, d_name, routine_name)
    ! Check if NaN values occur in the 3-D int data field d
    ! NOTE: parallelised!
    
    IMPLICIT NONE
    
    ! In/output variables:
    INTEGER,  DIMENSION(:,:,:),              INTENT(IN)    :: d
    CHARACTER(LEN=*),           OPTIONAL,    INTENT(IN)    :: d_name, routine_name
    
    ! Local variables:
    INTEGER                                                :: nx,ny,nz,i,j,k,i1,i2
    CHARACTER(LEN=256)                                     :: d_name_loc, routine_name_loc
    
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
    IF (PRESENT( routine_name)) THEN
      routine_name_loc = TRIM(routine_name)
    ELSE
      routine_name_loc = '?'
    END IF
    
    ! Inspect data field
    DO i = i1, i2
    DO j = 1, ny
    DO k = 1, nz
    
      ! Strangely enough, Fortran doesn't have an "isnan" function; instead,
      ! you use the property that a NaN is never equal to anything, including itself...
      
      IF (d( k,j,i) /= d( k,j,i)) THEN
        WRITE(0,'(A,I4,A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_int_3D - NaN detected at [', i, ',', j, ',', k, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( k,j,i) > HUGE( d( k,j,i))) THEN
        WRITE(0,'(A,I4,A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_int_3D - Inf detected at [', i, ',', j, ',', k, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      ELSEIF (d( k,j,i) < -HUGE( d( k,j,i))) THEN
        WRITE(0,'(A,I4,A,I4,A,I4,A,A,A,A,A)') 'check_for_NaN_int_3D - -Inf detected at [', i, ',', j, ',', k, &
          '] in variable "', TRIM(d_name_loc), '" from routine "', TRIM(routine_name_loc), '"'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
      
    END DO
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE check_for_NaN_int_3D

END MODULE utilities_module
