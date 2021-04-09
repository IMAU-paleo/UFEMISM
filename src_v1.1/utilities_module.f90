MODULE utilities_module

  USE mpi
  USE parallel_module,                 ONLY: par, sync, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE configuration_module,            ONLY: dp, C
  USE parameters_module,               ONLY: pi, earth_radius
  USE data_types_module,               ONLY: type_mesh, type_grid

CONTAINS
  
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
  
! == Smoothing operations on the mesh
  SUBROUTINE smooth_Gaussian_2D( mesh, grid, d_mesh, r)
    ! Use pseudo-conservative remapping to map the data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.
    
    USE mesh_mapping_module, ONLY: map_mesh2grid_2D, map_grid2mesh_2D
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:    ),              INTENT(INOUT) :: d_mesh
    REAL(dp),                                INTENT(IN)    :: r
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:  ), POINTER                    :: d_grid
    INTEGER                                                :: wd_grid
    
    ! Allocate shared memory for the gridded data
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, d_grid, wd_grid)
    
    ! Map data from the mesh to the grid
    CALL map_mesh2grid_2D( mesh, grid, d_mesh, d_grid)
    
    ! Smooth data on the grid
    CALL smooth_Gaussian_2D_grid( grid, d_grid, r)
    
    ! Map smoothed data back to the mesh
    CALL map_grid2mesh_2D( mesh, grid, d_grid, d_mesh)
    
    ! Deallocate gridded data
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE smooth_Gaussian_2D
  SUBROUTINE smooth_Gaussian_3D( mesh, grid, d_mesh, r, nz)
    ! Use pseudo-conservative remapping to map the data from the mesh
    ! to the square grid. Apply the smoothing on the gridded data, then map
    ! it back to the mesh. The numerical diffusion arising from the two mapping
    ! operations is not a problem since we're smoothing the data anyway.
    
    USE mesh_mapping_module, ONLY: map_mesh2grid_3D, map_grid2mesh_3D
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                         INTENT(IN)    :: mesh
    TYPE(type_grid),                         INTENT(IN)    :: grid
    REAL(dp), DIMENSION(:,:  ),              INTENT(INOUT) :: d_mesh
    REAL(dp),                                INTENT(IN)    :: r
    INTEGER,                                 INTENT(IN)    :: nz
    
    ! Local variables:
    REAL(dp), DIMENSION(:,:,:), POINTER                    :: d_grid
    INTEGER                                                :: wd_grid
    
    ! Allocate shared memory for the gridded data
    CALL allocate_shared_dp_3D( grid%nx, grid%ny, nz, d_grid, wd_grid)
    
    ! Map data from the mesh to the grid
    CALL map_mesh2grid_3D( mesh, grid, d_mesh, d_grid)
    
    ! Smooth data on the grid
    CALL smooth_Gaussian_3D_grid( grid, d_grid, r, nz)
    
    ! Map smoothed data back to the mesh
    CALL map_grid2mesh_3D( mesh, grid, d_grid, d_mesh)
    
    ! Deallocate gridded data
    CALL deallocate_shared( wd_grid)
    
  END SUBROUTINE smooth_Gaussian_3D

! == Smoothing operations on a regular square grid
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
    DO i = -n, n
      f(i) = EXP( -0.5_dp * (REAL(i,dp) * grid%dx/r)**2)
    END DO
    f = f / SUM(f)
    
    ! Allocate temporary shared memory for the extended and smoothed data fields
    CALL allocate_shared_dp_2D( grid%nx + 2*n, grid%ny + 2*n, d_ext,        wd_ext       )
    CALL allocate_shared_dp_2D( grid%nx + 2*n, grid%ny + 2*n, d_ext_smooth, wd_ext_smooth)
    
    ! Copy data to the extended array and fill in the margins
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      d_ext( i+n,j+n) = d( i,j)
    END DO
    END DO
    CALL sync
    DO j = grid%j1, grid%j2
      d_ext(           1:          n, j+n) = d( 1      ,j)
      d_ext( grid%nx+n+1:grid%nx+2*n, j+n) = d( grid%nx,j)
    END DO
    CALL sync
    DO i = grid%i1, grid%i2
      d_ext( i+n,           1:          n) = d( i,1      )
      d_ext( i+n, grid%ny+n+1:grid%ny+2*n) = d( i,grid%ny)
    END DO
    IF (par%master) THEN
      d_ext( 1:n,                     1:n                    ) = d( 1      ,1      )
      d_ext( 1:n,                     grid%ny+n+1:grid%ny+2*n) = d( 1      ,grid%ny)
      d_ext( grid%nx+n+1:grid%nx+2*n, 1:n                    ) = d( grid%nx,1      )
      d_ext( grid%nx+n+1:grid%nx+2*n, grid%ny+n+1:grid%ny+2*n) = d( grid%nx,grid%ny)
    END IF
    CALL sync
    
    ! Convolute extended data with the smoothing filter
    d_ext_smooth( grid%i1+n:grid%i2+n,:) = 0._dp
    DO i = grid%i1, grid%i2
    DO j = 1,       grid%ny
      DO ii = -n, n
        d_ext_smooth( i+n,j+n) = d_ext_smooth( i+n,j+n) + d_ext( i+n+ii,j+n) * f(ii)
      END DO
    END DO
    END DO
    CALL sync
    
    d_ext( grid%i1+n:grid%i2+n,:) = d_ext_smooth( grid%i1+n:grid%i2+n,:)
    CALL sync
    DO i = grid%i1, grid%i2
      d_ext( i+n,           1:          n) = d( i,1      )
      d_ext( i+n, grid%ny+n+1:grid%ny+2*n) = d( i,grid%ny)
    END DO
    CALL sync
    
    d_ext_smooth( grid%i1+n:grid%i2+n,:) = 0._dp
    DO j = grid%j1, grid%j2
    DO i = 1,       grid%nx
      DO jj = -n, n
        d_ext_smooth( i+n,j+n) = d_ext_smooth( i+n,j+n) + d_ext( i+n,j+n+jj) * f(jj)
      END DO
    END DO
    END DO
    CALL sync
    
    ! Copy data back
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      d( i,j) = d_ext_smooth( i+n, j+n)
    END DO
    END DO
    
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
    CALL sync
    DO i = grid%i1, grid%i2
      d_ext(           1:          n, i+n) = d( 1      ,i)
      d_ext( grid%ny+n+1:grid%ny+2*n, i+n) = d( grid%ny,i)
    END DO
    CALL sync
    DO j = grid%j1, grid%j2
      d_ext( j+n,           1:          n) = d( j,1      )
      d_ext( j+n, grid%nx+n+1:grid%nx+2*n) = d( j,grid%nx)
    END DO
    IF (par%master) THEN
      d_ext( 1:n,                     1:n                    ) = d( 1      ,1      )
      d_ext( 1:n,                     grid%nx+n+1:grid%nx+2*n) = d( 1      ,grid%nx)
      d_ext( grid%ny+n+1:grid%ny+2*n, 1:n                    ) = d( grid%ny,1      )
      d_ext( grid%ny+n+1:grid%ny+2*n, grid%nx+n+1:grid%nx+2*n) = d( grid%ny,grid%nx)
    END IF
    CALL sync
    
    ! Convolute extended data with the smoothing filter
    d_ext_smooth( :,grid%i1+n:grid%i2+n) = 0._dp
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

END MODULE utilities_module
