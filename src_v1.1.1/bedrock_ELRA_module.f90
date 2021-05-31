MODULE bedrock_ELRA_module

  USE mpi
  USE configuration_module,            ONLY: dp, C
  USE parallel_module,                 ONLY: par, sync, &
                                             allocate_shared_int_0D, allocate_shared_dp_0D, &
                                             allocate_shared_int_1D, allocate_shared_dp_1D, &
                                             allocate_shared_int_2D, allocate_shared_dp_2D, &
                                             allocate_shared_int_3D, allocate_shared_dp_3D, &
                                             deallocate_shared
  USE data_types_module,               ONLY: type_model_region, type_mesh, type_grid, type_ice_model, type_PD_data_fields, type_remapping
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  USE parameters_module,               ONLY: grav, seawater_density, ice_density, pi
  USE general_ice_model_data_module,   ONLY: is_floating
  USE mesh_mapping_module,             ONLY: reallocate_field_dp, map_mesh2grid_2D, map_grid2mesh_2D

  IMPLICIT NONE
    
CONTAINS

  ! The ELRA GIA model
  SUBROUTINE run_ELRA_model( region)
    ! Use the ELRA model to update bedrock elevation. Once every (dt_bedrock_ELRA) years,
    ! update deformation rates. In all other time steps, just incrementally add deformation.
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Local variables:
    INTEGER                                            :: cerr,ierr
    INTEGER                                            :: vi
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        region%t0_ELRA = region%time
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_ELRA_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! If needed, update the bedrock deformation rate
    IF (region%do_ELRA) THEN
      CALL calculate_ELRA_bedrock_deformation_rate( region%mesh, region%grid_GIA, region%ice)
      region%t0_ELRA = region%time
    END IF
    
    ! Update bedrock with last calculated deformation rate
    DO vi = region%mesh%v1, region%mesh%v2
      region%ice%Hb(  vi) = region%ice%Hb( vi) + region%ice%dHb_dt( vi) * region%dt
      region%ice%dHb( vi) = region%ice%Hb( vi) - region%PD%Hb( vi)
    END DO
    CALL sync
    
  END SUBROUTINE run_ELRA_model
  SUBROUTINE calculate_ELRA_bedrock_deformation_rate( mesh, grid, ice)
    ! Use the ELRA model to update bedrock deformation rates.
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    
    ! Local variables:
    INTEGER                                            :: vi,i,j,n,k,l
    REAL(dp)                                           :: Lr
    
    ! Influence radius of the lithospheric rigidity
    Lr = (C%ELRA_lithosphere_flex_rigidity / (C%ELRA_mantle_density * grav))**0.25_dp
    
    ! Calculate the absolute and relative surface loads on the mesh
    DO vi = mesh%v1, mesh%v2
    
      ! Absolute surface load
      IF (is_floating( ice%Hi( vi), ice%Hb( vi), ice%SL( vi))) THEN
        ice%surface_load_mesh( vi) = (ice%SL( vi) - ice%Hb( vi)) * grid%dx**2 * seawater_density
      ELSEIF (ice%Hi( vi) > 0._dp) THEN
        ice%surface_load_mesh( vi) =  ice%Hi( vi)                * grid%dx**2 * ice_density
      ELSE
        ice%surface_load_mesh( vi) = 0._dp
      END IF
      
      ! Relative surface load
      ice%surface_load_rel_mesh( vi) = ice%surface_load_mesh( vi) - ice%surface_load_PD_mesh( vi)
      
    END DO
    CALL sync
    
    ! Map relative surface load to the GIA grid
    CALL map_mesh2grid_2D( mesh, grid, ice%surface_load_rel_mesh, ice%surface_load_rel_grid)
    
    ! Fill in the "extended" relative surface load (extrapolating the domain so we can do the 2D convolution)
    
    n = ice%flex_prof_rad
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%surface_load_rel_ext_grid( i+n,j+n) = ice%surface_load_rel_grid( i,j)
    END DO
    END DO
    CALL sync
    DO i = grid%i1, grid%i2
      ice%surface_load_rel_ext_grid( i,           1:          n) = ice%surface_load_rel_grid( i,1      )
      ice%surface_load_rel_ext_grid( i, grid%ny+n+1:grid%ny+2*n) = ice%surface_load_rel_grid( i,grid%ny)
    END DO
    CALL sync
    DO j = grid%j1, grid%j2
      ice%surface_load_rel_ext_grid(           1:          n, j) = ice%surface_load_rel_grid( 1      ,j)
      ice%surface_load_rel_ext_grid( grid%nx+n+1:grid%nx+2*n, j) = ice%surface_load_rel_grid( grid%nx,j)
    END DO
    CALL sync
    
    ! Calculate the equilibrium bedrock deformation for this surface load
    ! ( = convolute([surface load],[flexural profile]))
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ice%dHb_eq_grid( i,j) = 0._dp
      
      DO k = -n,n
      DO l = -n,n
        
        ice%dHb_eq_grid( i,j) = ice%dHb_eq_grid( i,j) + &
          (0.5_dp * grav * Lr**2 /(pi * C%ELRA_lithosphere_flex_rigidity) * ice%surface_load_rel_ext_grid( i+n+l, j+n+k) * ice%flex_prof_grid( k+n+1,l+n+1))
        
      END DO
      END DO
      
    END DO
    END DO
    CALL sync
    
    ! Map the actual bedrock deformation to the grid
    CALL map_mesh2grid_2D( mesh, grid, ice%dHb, ice%dHb_grid)
    
    ! Calculate the bedrock deformation rate from the difference between the current and the equilibrium deformation
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      ice%dHb_dt_grid( i,j) = (ice%dHb_eq_grid( i,j) - ice%dHb_grid( i,j)) / C%ELRA_bedrock_relaxation_time
    END DO
    END DO
    CALL sync
    
    ! Map the deformation rate back to the model mesh
    CALL map_grid2mesh_2D( mesh, grid, ice%dHb_dt_grid, ice%dHb_dt)
    
  END SUBROUTINE calculate_ELRA_bedrock_deformation_rate
  
  SUBROUTINE initialise_ELRA_model( mesh, grid, ice, PD)
    ! Allocate and initialise the ELRA GIA model
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_PD_data_fields),           INTENT(IN)    :: PD
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: i,j,n,k,l
    REAL(dp)                                           :: Lr, r
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_ELRA_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    IF (par%master) WRITE (0,*) '  Initialising ELRA GIA model...'
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV,          ice%surface_load_PD_mesh,  ice%wsurface_load_PD_mesh )
    CALL allocate_shared_dp_1D( mesh%nV,          ice%surface_load_mesh,     ice%wsurface_load_mesh    )
    CALL allocate_shared_dp_1D( mesh%nV,          ice%surface_load_rel_mesh, ice%wsurface_load_rel_mesh)
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, ice%surface_load_rel_grid, ice%wsurface_load_rel_grid)
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, ice%dHb_eq_grid,           ice%wdHb_eq_grid          )
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, ice%dHb_grid,              ice%wdHb_grid             )
    CALL allocate_shared_dp_2D( grid%nx, grid%ny, ice%dHb_dt_grid,           ice%wdHb_dt_grid          )
    
    ! Fill in the 2D flexural profile (= Kelvin function), with which
    ! a surface load is convoluted to find the surface deformation
    ! ============================================================
    
    ! Influence radius of the lithospheric rigidity
    Lr = (C%ELRA_lithosphere_flex_rigidity / (C%ELRA_mantle_density * grav))**0.25_dp
    
    ! Calculate radius (in number of grid cells) of the flexural profile
    CALL allocate_shared_int_0D( ice%flex_prof_rad, ice%wflex_prof_rad)
    IF (par%master) ice%flex_prof_rad = MIN( CEILING(grid%dx/2._dp), MAX(1, INT(6._dp * Lr / grid%dx) - 1))
    CALL sync
    
    n = 2 * ice%flex_prof_rad + 1
    CALL allocate_shared_dp_2D( n,         n,         ice%flex_prof_grid,            ice%wflex_prof_grid           )
    CALL allocate_shared_dp_2D( grid%nx+n, grid%ny+n, ice%surface_load_rel_ext_grid, ice%wsurface_load_rel_ext_grid)
    
    ! Calculate flexural profile
    IF (par%master) THEN
      DO i = -ice%flex_prof_rad, ice%flex_prof_rad
      DO j = -ice%flex_prof_rad, ice%flex_prof_rad
        l = i+ice%flex_prof_rad+1
        k = j+ice%flex_prof_rad+1
        r = grid%dx * SQRT( (REAL(i,dp))**2 + (REAL(j,dp))**2)
        ice%flex_prof_grid( l,k) = kelvin(r / Lr)
      END DO
      END DO
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Calculate the PD reference load
    ! ===============================
    
    CALL initialise_ELRA_PD_reference_load( mesh, grid, ice, PD)
        
  END SUBROUTINE initialise_ELRA_model
  SUBROUTINE initialise_ELRA_PD_reference_load( mesh, grid, ice, PD)
      
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    TYPE(type_grid),                     INTENT(IN)    :: grid
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_PD_data_fields),           INTENT(IN)    :: PD
    
    ! Local variables:
    INTEGER                                            :: vi
    
    ! Calculate PD reference load on the mesh
    DO vi = mesh%v1, mesh%v2
      IF (is_floating( PD%Hi( vi), PD%Hb( vi), 0._dp)) THEN
        ice%surface_load_PD_mesh( vi) = -PD%Hb( vi) * grid%dx**2 * seawater_density
      ELSEIF (PD%Hi( vi) > 0._dp) THEN
        ice%surface_load_PD_mesh( vi) =  PD%Hi( vi) * grid%dx**2 * ice_density
      END IF
    END DO
    CALL sync
    
  END SUBROUTINE initialise_ELRA_PD_reference_load
  SUBROUTINE remap_ELRA_model( mesh_old, mesh_new, map, ice, PD, grid)
    ! Remap or reallocate all the data fields
  
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_PD_data_fields),           INTENT(IN)    :: PD
    TYPE(type_grid),                     INTENT(IN)    :: grid
    
    ! Local variables:
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: int_dummy
    
    ! Exceptions for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream' .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        RETURN
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in remap_ELRA_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    int_dummy = mesh_old%nV
    int_dummy = map%trilin%vi( 1,1)
    
    CALL reallocate_field_dp( mesh_new%nV, ice%surface_load_PD_mesh,  ice%wsurface_load_PD_mesh )
    CALL reallocate_field_dp( mesh_new%nV, ice%surface_load_mesh,     ice%wsurface_load_mesh    )
    CALL reallocate_field_dp( mesh_new%nV, ice%surface_load_rel_mesh, ice%wsurface_load_rel_mesh)
    
    ! Recalculate the PD reference load on the GIA grid
    CALL initialise_ELRA_PD_reference_load( mesh_new, grid, ice, PD)
    
  END SUBROUTINE remap_ELRA_model

  ! The Kelvin function (just a wrapper for the Zhang & Jin "klvna" subroutine)
  FUNCTION kelvin(x) RESULT(kei)
    ! Evaluate the Kelvin function in a given point at a distance x from the load. Currently this is implemented 
    ! using the klvna routine in mklvna.f (see that routine for details).
    
    IMPLICIT NONE

    ! Input variables:
    ! x is the distance between a considered load and the point in which we want to know the deformation, x is in units Lr
    REAL(dp), INTENT(IN) :: x  
 
    ! Result variables:
    REAL(dp)             :: kei
 
    ! Local variables:
    REAL(dp)     :: xdp,berdp,beidp,gerdp,geidp,derdp,deidp,herdp,heidp

    xdp = x
    CALL klvna( xdp,berdp,beidp,gerdp,geidp,derdp,deidp,herdp,heidp)
    kei = geidp
    
  END FUNCTION kelvin
  SUBROUTINE klvna ( x, ber, bei, ger, gei, der, dei, her, hei )

!*****************************************************************************80
!
!! KLVNA: Kelvin functions ber(x), bei(x), ker(x), and kei(x), and derivatives.
!
!  Licensing:
!
!    This routine is copyrighted by Shanjie Zhang and Jianming Jin.  However, 
!    they give permission to incorporate this routine into a user program 
!    provided that the copyright is acknowledged.
!
!  Modified:
!
!    03 August 2012
!
!  Author:
!
!    Shanjie Zhang, Jianming Jin
!
!  Reference:
!
!    Shanjie Zhang, Jianming Jin,
!    Computation of Special Functions,
!    Wiley, 1996,
!    ISBN: 0-471-11963-6,
!    LC: QA351.C45.
!
!  Parameters:
!
!    Input, real (dp) X, the argument.
!
!    Output, real (dp) BER, BEI, GER, GEI, DER, DEI, HER, HEI, 
!    the values of ber x, bei x, ker x, kei x, ber'x, bei'x, ker'x, kei'x.
!

  USE parameters_module, ONLY: pi
  
  IMPLICIT NONE
  
  ! In/output variables:
  REAL(dp), INTENT(IN )   :: x
  REAL(dp), INTENT(OUT)   :: ber
  REAL(dp), INTENT(OUT)   :: bei
  REAL(dp), INTENT(OUT)   :: ger
  REAL(dp), INTENT(OUT)   :: gei
  REAL(dp), INTENT(OUT)   :: der
  REAL(dp), INTENT(OUT)   :: dei
  REAL(dp), INTENT(OUT)   :: her
  REAL(dp), INTENT(OUT)   :: hei

  ! Local variables:
  REAL(dp)   :: cn0
  REAL(dp)   :: cp0
  REAL(dp)   :: cs
  REAL(dp)   :: el
  REAL(dp)   :: eps
  REAL(dp)   :: fac
  REAL(dp)   :: gs

  INTEGER    :: k
  INTEGER    :: Km
  INTEGER    :: m

  REAL(dp)   :: pn0
  REAL(dp)   :: pn1
  REAL(dp)   :: pp0
  REAL(dp)   :: pp1
  REAL(dp)   :: qn0
  REAL(dp)   :: qn1
  REAL(dp)   :: qp0
  REAL(dp)   :: qp1
  REAL(dp)   :: r
  REAL(dp)   :: r0
  REAL(dp)   :: r1
  REAL(dp)   :: rc
  REAL(dp)   :: rs
  REAL(dp)   :: sn0
  REAL(dp)   :: sp0
  REAL(dp)   :: ss
  REAL(dp)   :: x2
  REAL(dp)   :: x4
  REAL(dp)   :: xc1
  REAL(dp)   :: xc2
  REAL(dp)   :: xd
  REAL(dp)   :: xe1
  REAL(dp)   :: xe2
  REAL(dp)   :: xt

  el = 0.5772156649015329_dp
  eps = 1.0D-15

  if ( x == 0.0_dp ) then
    ber = 1.0_dp
    bei = 0.0_dp
    ger = 1.0D+300
    gei = -0.25_dp * pi
    der = 0.0_dp
    dei = 0.0_dp
    her = -1.0D+300
    hei = 0.0_dp
    return
  end if

  x2 = 0.25_dp * x * x
  x4 = x2 * x2

  if ( abs ( x ) < 10.0_dp ) then

    ber = 1.0_dp
    r = 1.0_dp
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m - 1.0_dp ) ** 2 * x4
      ber = ber + r
      if ( abs ( r ) < abs ( ber ) * eps ) then
        exit
      end if
    end do

    bei = x2
    r = x2
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      bei = bei + r
      if ( abs ( r ) < abs ( bei ) * eps ) then
        exit
      end if
    end do

    ger = - ( log ( x / 2.0_dp ) + el ) * ber + 0.25_dp * pi * bei
    r = 1.0_dp
    gs = 0.0_dp
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m - 1.0_dp ) ** 2 * x4
      gs = gs + 1.0_dp / ( 2.0_dp * m - 1.0_dp ) + 1.0_dp / ( 2.0_dp * m )
      ger = ger + r * gs
      if ( abs ( r * gs ) < abs ( ger ) * eps ) then
        exit
      end if
    end do

    gei = x2 - ( log ( x / 2.0_dp ) + el ) * bei - 0.25_dp * pi * ber
    r = x2
    gs = 1.0_dp
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      gs = gs + 1.0_dp / ( 2.0_dp * m ) + 1.0_dp / ( 2.0_dp * m + 1.0_dp )
      gei = gei + r * gs
      if ( abs ( r * gs ) < abs ( gei ) * eps ) then
        exit
      end if
    end do

    der = -0.25_dp * x * x2
    r = der
    do m = 1, 60
      r = -0.25_dp * r / m / ( m + 1.0_dp ) &
        / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      der = der + r
      if ( abs ( r ) < abs ( der ) * eps ) then
        exit
      end if
    end do

    dei = 0.5_dp * x
    r = dei
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2.0_dp * m - 1.0_dp ) &
        / ( 2.0_dp * m + 1.0_dp ) * x4
      dei = dei + r
      if ( abs ( r ) < abs ( dei ) * eps ) then
        exit
      end if
    end do

    r = -0.25_dp * x * x2
    gs = 1.5_dp
    her = 1.5_dp * r - ber / x &
      - ( log ( x / 2.0_dp ) + el ) * der + 0.25_dp * pi * dei
    do m = 1, 60
      r = -0.25_dp * r / m / ( m + 1.0_dp ) &
        / ( 2.0_dp * m + 1.0_dp ) ** 2 * x4
      gs = gs + 1.0_dp / ( 2 * m + 1.0_dp ) + 1.0_dp &
        / ( 2 * m + 2.0_dp )
      her = her + r * gs
      if ( abs ( r * gs ) < abs ( her ) * eps ) then
        exit
      end if
    end do

    r = 0.5_dp * x
    gs = 1.0_dp
    hei = 0.5_dp * x - bei / x &
      - ( log ( x / 2.0_dp ) + el ) * dei - 0.25_dp * pi * der
    do m = 1, 60
      r = -0.25_dp * r / ( m * m ) / ( 2 * m - 1.0_dp ) &
        / ( 2 * m + 1.0_dp ) * x4
      gs = gs + 1.0_dp / ( 2.0_dp * m ) + 1.0_dp &
        / ( 2 * m + 1.0_dp )
      hei = hei + r * gs
      if ( abs ( r * gs ) < abs ( hei ) * eps ) then 
        return
      end if
    end do

  else

    pp0 = 1.0_dp
    pn0 = 1.0_dp
    qp0 = 0.0_dp
    qn0 = 0.0_dp
    r0 = 1.0_dp

    if ( abs ( x ) < 40.0_dp ) then
      km = 18
    else
      km = 10
    end if

    fac = 1.0_dp
    do k = 1, km
      fac = -fac
      xt = 0.25_dp * k * pi - int ( 0.125_dp * k ) * 2.0_dp * pi
      cs = cos ( xt )
      ss = sin ( xt )
      r0 = 0.125_dp * r0 * ( 2.0_dp * k - 1.0_dp ) ** 2 / k / x
      rc = r0 * cs
      rs = r0 * ss
      pp0 = pp0 + rc
      pn0 = pn0 + fac * rc
      qp0 = qp0 + rs
      qn0 = qn0 + fac * rs
    end do

    xd = x / sqrt (2.0_dp )
    xe1 = exp ( xd )
    xe2 = exp ( - xd )
    xc1 = 1.0_dp / sqrt ( 2.0_dp * pi * x )
    xc2 = sqrt ( 0.5_dp * pi / x )
    cp0 = cos ( xd + 0.125_dp * pi )
    cn0 = cos ( xd - 0.125_dp * pi )
    sp0 = sin ( xd + 0.125_dp * pi )
    sn0 = sin ( xd - 0.125_dp * pi )
    ger = xc2 * xe2 * (  pn0 * cp0 - qn0 * sp0 )
    gei = xc2 * xe2 * ( -pn0 * sp0 - qn0 * cp0 )
    ber = xc1 * xe1 * (  pp0 * cn0 + qp0 * sn0 ) - gei / pi
    bei = xc1 * xe1 * (  pp0 * sn0 - qp0 * cn0 ) + ger / pi
    pp1 = 1.0_dp
    pn1 = 1.0_dp
    qp1 = 0.0_dp
    qn1 = 0.0_dp
    r1 = 1.0_dp
    fac = 1.0_dp

    do k = 1, km
      fac = -fac
      xt = 0.25_dp * k * pi - int ( 0.125_dp * k ) * 2.0_dp * pi
      cs = cos ( xt )
      ss = sin ( xt )
      r1 = 0.125_dp * r1 &
        * ( 4.0_dp - ( 2.0_dp * k - 1.0_dp ) ** 2 ) / k / x
      rc = r1 * cs
      rs = r1 * ss
      pp1 = pp1 + fac * rc
      pn1 = pn1 + rc
      qp1 = qp1 + fac * rs
      qn1 = qn1 + rs
    end do

    her = xc2 * xe2 * ( - pn1 * cn0 + qn1 * sn0 )
    hei = xc2 * xe2 * (   pn1 * sn0 + qn1 * cn0 )
    der = xc1 * xe1 * (   pp1 * cp0 + qp1 * sp0 ) - hei / pi
    dei = xc1 * xe1 * (   pp1 * sp0 - qp1 * cp0 ) + her / pi

  end if

  END SUBROUTINE klvna

END MODULE bedrock_ELRA_module
