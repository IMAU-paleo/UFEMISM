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
  SUBROUTINE add_matrix_matrix_CSR( AA, BB, CC)
    ! Perform the matrix multiplication C = A+B, where all three
    ! matrices are provided in CSR format (parallelised)
      
    IMPLICIT NONE
    
    ! In- and output variables:
    TYPE(type_sparse_matrix_CSR),        INTENT(IN)    :: AA, BB
    TYPE(type_sparse_matrix_CSR),        INTENT(INOUT) :: CC
    
    ! Local variables:
    INTEGER                                            :: i1, i2, ic, jc
    INTEGER                                            :: ka1, ka2, ka, ja
    INTEGER                                            :: kb1, kb2, kb, jb
    LOGICAL                                            :: is_nonzero
    REAL(dp)                                           :: Cij
    
    ! Safety
    IF (AA%m /= BB%m .OR. AA%n /= BB%n) THEN
      IF (par%master) WRITE(0,*) 'append_overwrite_CSR - ERROR: A and B are not of the same size!'
      CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
    END IF
    
    ! Partition rows over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    ! Allocate distributed shared memory for C
    CALL allocate_matrix_CSR_dist( CC, AA%m, AA%n, AA%nnz + BB%nnz)
    
    ! Initialise
    CC%ptr = 1
    CC%nnz = 0
    
    DO ic = i1, i2
      
      DO jc = 1, CC%n
        
        is_nonzero = .FALSE.
        Cij        = 0._dp
        
        ! Check if this column has an entry in A; if so, add it
        ka1 = AA%ptr( ic)
        ka2 = AA%ptr( ic+1) - 1
        DO ka = ka1, ka2
          ja = AA%index( ka)
          IF (ja == jc) THEN
            is_nonzero = .TRUE.
            Cij = Cij + AA%val( ka)
            EXIT
          END IF
        END DO
        
        ! Check if this column has an entry in B; if so, add it
        kb1 = BB%ptr( ic)
        kb2 = BB%ptr( ic+1) - 1
        DO kb = kb1, kb2
          jb = BB%index( kb)
          IF (jb == jc) THEN
            is_nonzero = .TRUE.
            Cij = Cij + BB%val( kb)
            EXIT
          END IF
        END DO
        
        ! If the result is non-zero, list it in C
        IF (is_nonzero) THEN
          CC%nnz  = CC%nnz + 1
          CC%index( CC%nnz) = jc
          CC%val(   CC%nnz) = Cij
        END IF
        
      END DO ! DO j = 1, CC%n
      
      ! Finalise this row of C
      CC%ptr( ic+1) = CC%nnz+1
      
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
      IF (par%master) WRITE(0,*) 'append_overwrite_CSR - ERROR: A and B are not of the same size!'
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
  
  SUBROUTINE solve_matrix_equation_CSR( AA, b, x, choice_matrix_solver, SOR_nit, SOR_tol, SOR_omega, PETSc_rtol, PETSc_abstol)
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
      
      CALL solve_matrix_equation_CSR_SOR( AA, b, x, SOR_nit, SOR_tol, SOR_omega)
      
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
  SUBROUTINE solve_matrix_equation_CSR_SOR( AA, b, x, nit, tol, omega)
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
    REAL(dp)                                           :: lhs, res, cij, res_max, omega_dyn
    
    ! Partition equations over the processors
    CALL partition_list( AA%m, par%i, par%n, i1, i2)
    
    omega_dyn = omega
    
    res_max = tol * 2._dp
    it = 0
    SOR_iterate: DO WHILE (res_max > tol .AND. it < nit)
      it = it+1
      
      res_max = 0._dp

      DO i = i1, i2
      
        lhs = 0._dp
        cij = 0._dp
        DO k = AA%ptr( i), AA%ptr( i+1)-1
          j = AA%index( k)
          lhs = lhs + AA%val( k) * x( j)
          IF (j == i) cij = AA%val( k)
        END DO
        
        res = (lhs - b( i)) / cij
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
    
  END SUBROUTINE solve_matrix_equation_CSR_SOR
  
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
