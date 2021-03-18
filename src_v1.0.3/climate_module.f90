MODULE climate_module

  USE mpi
  USE configuration_module,          ONLY: dp, C
  USE parallel_module,               ONLY: par, sync, &
                                           allocate_shared_int_0D, allocate_shared_dp_0D, &
                                           allocate_shared_int_1D, allocate_shared_dp_1D, &
                                           allocate_shared_int_2D, allocate_shared_dp_2D, &
                                           allocate_shared_int_3D, allocate_shared_dp_3D, &
                                           deallocate_shared
  USE data_types_module,             ONLY: type_mesh, type_ice_model, type_climate_model, type_init_data_fields, &
                                           type_climate_matrix, type_subclimate_global, type_subclimate_mesh, type_remapping
  USE netcdf_module,                 ONLY: debug, write_to_debug_file, inquire_PD_obs_data_file, read_PD_obs_data_file
  USE mesh_mapping_module,           ONLY: get_glob_to_mesh_map, map_cart_to_mesh, map_cart_to_mesh_3D, Remap_field_dp, Remap_field_dp_3D, &
                                           Remap_field_dp_monthly, Reallocate_field_dp, Reallocate_field_dp_3D, Reallocate_field_int

  IMPLICIT NONE
    
CONTAINS

  ! Run the matrix climate model on a region mesh
  SUBROUTINE run_climate_model( mesh, ice, climate, time)
    ! Run the regional climate model
    
    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh   
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    REAL(dp),                            INTENT(IN)    :: time
    
    ! Local variables:
    REAL(dp), PARAMETER                                :: lambda = -0.008_dp    ! Constant lapse rate of 8 K/km
    INTEGER                                            :: cerr, ierr
    INTEGER                                            :: vi, m
    REAL(dp)                                           :: dT_lapse
    
    ! Check if we need to apply any special benchmark experiment climate
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6') THEN
          
        CALL EISMINT_climate( mesh, ice, climate, time)
        RETURN
        
      ELSEIF (C%choice_benchmark_experiment == 'Halfar') THEN
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        RETURN
      ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        RETURN
      ELSE
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_climate_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! For now, only apply a basic lapse-rate based temperature correction
    ! to the PD observations.
    climate%applied%Precip(mesh%v1:mesh%v2,:) = climate%PD_obs%Precip(mesh%v1:mesh%v2,:)
    
    DO vi = mesh%v1, mesh%v2
    
      dT_lapse = (ice%Hs(vi) - climate%PD_obs%Hs(vi)) * lambda
    
      DO m = 1, 12
        climate%applied%T2m(vi,m) = climate%PD_obs%T2m(vi,m) + dT_lapse
      END DO
    END DO
    CALL sync
    
!   IF (par%master) WRITE(0,*) '   DENK DROM: run_climate_model - increasing temperature by 5 degrees'   
   climate%applied%T2m(mesh%v1:mesh%v2,:) = climate%applied%T2m(mesh%v1:mesh%v2,:) + 5._dp
   CALL sync
    
  END SUBROUTINE run_climate_model
  
  ! Temperature parameterisation for the EISMINT experiments
  SUBROUTINE EISMINT_climate( mesh, ice, climate, time)
    ! Simple lapse-rate temperature parameterisation
    
    USe parameters_module,           ONLY: pi
    
    IMPLICIT NONE

    TYPE(type_mesh),                     INTENT(IN)    :: mesh   
    TYPE(type_ice_model),                INTENT(IN)    :: ice 
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    REAL(dp),                            INTENT(IN)    :: time
    
    REAL(dp), PARAMETER                                :: lambda = -0.010_dp
    
    INTEGER                                            :: vi, m
    REAL(dp)                                           :: dT_lapse, d, dT
    
    ! Set precipitation to zero - SMB is parameterised anyway...
    climate%applied%Precip(mesh%v1:mesh%v2,:) = 0._dp
    
    ! Surface temperature for fixed or moving margin experiments
    IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_3') THEN
      ! Moving margin
          
      DO vi = mesh%v1, mesh%v2    
      
        dT_lapse = ice%Hs(vi) * lambda
          
        DO m = 1, 12
          climate%applied%T2m(vi,m) = 270._dp + dT_lapse
        END DO
      END DO
      
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_6') THEN
      ! Fixed margin
    
      DO vi = mesh%v1, mesh%v2    
        d = MAX( ABS(mesh%V(vi,1)/1000._dp), ABS(mesh%V(vi,2)/1000._dp))
    
        DO m = 1, 12
          climate%applied%T2m(vi,m) = 239._dp + (8.0E-08_dp * d**3)
        END DO
      END DO
      
    END IF
    CALL sync
    
    ! Glacial cycles
    IF     (C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_5') THEN
      IF (time > 0._dp) THEN
        dT = 10._dp * SIN(2 * pi * time / 20000._dp)
        DO vi = mesh%v1, mesh%v2 
          climate%applied%T2m(vi,:) = climate%applied%T2m(vi,:) + dT
        END DO
      END IF
    ELSEIF (C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_6') THEN
      IF (time > 0._dp) THEN
        dT = 10._dp * SIN(2 * pi * time / 40000._dp)
        DO vi = mesh%v1, mesh%v2 
          climate%applied%T2m(vi,:) = climate%applied%T2m(vi,:) + dT
        END DO
      END IF
    END IF
    CALL sync
    
  END SUBROUTINE EISMINT_climate
  
  ! Initialising the region-specific climate model, containing all the subclimates
  ! (PD observations, GCM snapshots and the applied climate) on the mesh
  SUBROUTINE initialise_climate_model( climate, matrix, mesh)
    ! Allocate shared memory for the regional climate models, containing the PD observed,
    ! GCM snapshots and applied climates as "subclimates"
    
    IMPLICIT NONE
    
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    TYPE(type_climate_matrix),           INTENT(IN)    :: matrix
    TYPE(type_mesh),                     INTENT(IN)    :: mesh  
        
    ! Allocate shared memory
    CALL initialise_subclimate( mesh, climate%PD_obs)
    CALL initialise_subclimate( mesh, climate%GCM_PI)
    CALL initialise_subclimate( mesh, climate%GCM_LGM)
    CALL initialise_subclimate( mesh, climate%applied)
    
    ! Map subclimates from global grid to mesh
    CALL map_subclimate_to_mesh( mesh,  matrix%PD_obs, climate%PD_obs)
  
    ! Initialise applied climate with present-day observations
    climate%applied%T2m(     mesh%v1:mesh%v2,:) = climate%PD_obs%T2m(     mesh%v1:mesh%v2,:)
    climate%applied%Precip(  mesh%v1:mesh%v2,:) = climate%PD_obs%Precip(  mesh%v1:mesh%v2,:)
    climate%applied%Hs(      mesh%v1:mesh%v2  ) = climate%PD_obs%Hs(      mesh%v1:mesh%v2  )
    climate%applied%Albedo(  mesh%v1:mesh%v2,:) = climate%PD_obs%Albedo(  mesh%v1:mesh%v2,:)
    climate%applied%Wind_WE( mesh%v1:mesh%v2,:) = climate%PD_obs%Wind_WE( mesh%v1:mesh%v2,:)
    climate%applied%Wind_SN( mesh%v1:mesh%v2,:) = climate%PD_obs%Wind_SN( mesh%v1:mesh%v2,:)
  
  END SUBROUTINE initialise_climate_model  
  SUBROUTINE initialise_subclimate( mesh, subclimate)
    ! Allocate shared memory for a "subclimate" (PD observed, GCM snapshot or applied climate) on the mesh
    
    IMPLICIT NONE
    
    TYPE(type_mesh),                     INTENT(IN)    :: mesh  
    TYPE(type_subclimate_mesh),          INTENT(INOUT) :: subclimate
    
    ! General forcing info (not relevant for PD observations but required for compatibility with GCM snapshots)
    CALL allocate_shared_dp_0D(              subclimate%CO2,        subclimate%wCO2       )
    CALL allocate_shared_dp_0D(              subclimate%orbit_time, subclimate%worbit_time)
    CALL allocate_shared_dp_0D(              subclimate%orbit_ecc,  subclimate%worbit_ecc )
    CALL allocate_shared_dp_0D(              subclimate%orbit_obl,  subclimate%worbit_obl )
    CALL allocate_shared_dp_0D(              subclimate%orbit_pre,  subclimate%worbit_pre )
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( mesh%nV, 12, subclimate%T2m,     subclimate%wT2m    )
    CALL allocate_shared_dp_2D( mesh%nV, 12, subclimate%Precip,  subclimate%wPrecip )
    CALL allocate_shared_dp_1D( mesh%nV,     subclimate%Hs,      subclimate%wHs     )
    CALL allocate_shared_dp_2D( mesh%nV, 12, subclimate%Albedo,  subclimate%wAlbedo )
    CALL allocate_shared_dp_2D( mesh%nV, 12, subclimate%Q_TOA,   subclimate%wQ_TOA  )
    CALL allocate_shared_dp_2D( mesh%nV, 12, subclimate%Wind_WE, subclimate%wWind_WE)
    CALL allocate_shared_dp_2D( mesh%nV, 12, subclimate%Wind_SN, subclimate%wWind_SN)
    
  END SUBROUTINE initialise_subclimate
  
  ! Map a global subclimate from the matrix (PD observed or GCM snapshot) to a region mesh
  SUBROUTINE map_subclimate_to_mesh( mesh,  gc, mc)
    ! Map data from a global "subclimate" (PD observed or GCM snapshot) to the mesh
    
    IMPLICIT NONE
    
    TYPE(type_mesh),                     INTENT(IN)    :: mesh 
    TYPE(type_subclimate_global),        INTENT(IN)    :: gc     ! Global climate
    TYPE(type_subclimate_mesh),          INTENT(INOUT) :: mc     ! Mesh   climate
    
    ! Local variables
    INTEGER,  DIMENSION(:,:  ), POINTER           :: map_mean_i  ! On cart grid: vertex indices
    REAL(dp), DIMENSION(:,:  ), POINTER           :: map_mean_w  ! On cart grid: weights
    INTEGER,  DIMENSION(:,:  ), POINTER           :: map_bilin_i ! On mesh     : cart indices (4x i and j))
    REAL(dp), DIMENSION(:,:  ), POINTER           :: map_bilin_w ! On mesh     : weights      (4x)
    INTEGER                                       :: wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w  
    
    IF (par%master) WRITE(0,*) '   Mapping ', TRIM(gc%name), ' data from global grid to mesh...'
    
    CALL get_glob_to_mesh_map(mesh, gc%lat, gc%lon, gc%nlat, gc%nlon,  map_mean_i,  map_mean_w,  map_bilin_i,  map_bilin_w, &
                                                                 wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w)
    CALL sync
    
    ! Map global climate data to the mesh
    CALL map_cart_to_mesh_3D( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, gc%nlat, 12, gc%i1, gc%i2,   gc%T2m,      mc%T2m     )
    CALL map_cart_to_mesh_3D( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, gc%nlat, 12, gc%i1, gc%i2,   gc%Precip,   mc%Precip  )
    CALL map_cart_to_mesh(    mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, gc%nlat,     gc%i1, gc%i2,   gc%Hs,       mc%Hs      )
    CALL map_cart_to_mesh_3D( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, gc%nlat, 12, gc%i1, gc%i2,   gc%Albedo,   mc%Albedo  )
    CALL map_cart_to_mesh_3D( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, gc%nlat, 12, gc%i1, gc%i2,   gc%Wind_WE,  mc%Wind_WE )
    CALL map_cart_to_mesh_3D( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, gc%nlat, 12, gc%i1, gc%i2,   gc%Wind_SN,  mc%Wind_SN )
    
    ! Deallocate shared memory for the mapping arrays
    CALL deallocate_shared( wmap_mean_i)
    CALL deallocate_shared( wmap_mean_w)
    CALL deallocate_shared( wmap_bilin_i)
    CALL deallocate_shared( wmap_bilin_w)
  
  END SUBROUTINE map_subclimate_to_mesh

  ! Initialising the climate matrix, containing all the global subclimates
  ! (PD observations and GCM snapshots)
  SUBROUTINE initialise_climate_matrix( matrix)
    ! Allocate shared memory for the global climate matrix
  
    IMPLICIT NONE
    
    TYPE(type_climate_matrix),      INTENT(INOUT) :: matrix
    
    IF (par%master) WRITE(0,*) ' Initialising the climate matrix...'
    
    ! The global ERA40 climate
    CALL initialise_PD_obs_data_fields(matrix%PD_obs, 'ERA40')
    
    ! The differenct GCM snapshots
!    CALL initialise_snapshot(nlat, nlon, matrix%PI)
!    CALL initialise_snapshot(nlat, nlon, matrix%LGM)
    
  END SUBROUTINE initialise_climate_matrix  
  SUBROUTINE initialise_PD_obs_data_fields( PD_obs, name)
    ! Allocate shared memory for the global PD observed climate data fields (stored in the climate matrix),
    ! read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_subclimate_global),   INTENT(INOUT) :: PD_obs
    CHARACTER(LEN=*),               INTENT(IN)    :: name
    
    ! Local variables:
    INTEGER                                       :: cerr, ierr
    
    PD_obs%name = name 
    
    ! Allocate memory for general forcing info (not relevant for PD but required for compatibility with GCM snapshots)
    CALL allocate_shared_dp_0D(                PD_obs%CO2,        PD_obs%wCO2       )
    CALL allocate_shared_dp_0D(                PD_obs%orbit_time, PD_obs%worbit_time)
    CALL allocate_shared_dp_0D(                PD_obs%orbit_ecc,  PD_obs%worbit_ecc )
    CALL allocate_shared_dp_0D(                PD_obs%orbit_obl,  PD_obs%worbit_obl )
    CALL allocate_shared_dp_0D(                PD_obs%orbit_pre,  PD_obs%worbit_pre )
        
    PD_obs%netcdf%filename   = C%filename_PD_obs_climate 
    
    ! Determine size of Cartesian grid, so that memory can be allocated accordingly.
    CALL allocate_shared_int_0D(       PD_obs%nlon, PD_obs%wnlon     )
    CALL allocate_shared_int_0D(       PD_obs%nlat, PD_obs%wnlat     )
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    IF (par%master) THEN
    
      IF (C%do_benchmark_experiment) THEN
        IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
            C%choice_benchmark_experiment == 'Halfar' .OR. &
            C%choice_benchmark_experiment == 'Bueler' .OR. &
            C%choice_benchmark_experiment == 'MISMIP_mod'.OR. &
            C%choice_benchmark_experiment == 'mesh_generation_test') THEN
      
          PD_obs%nlon = 10
          PD_obs%nlat = 10
          
        ELSE
          WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_obs_data_fields!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
      ELSE ! IF (C%do_benchmark_experiment) THEN
      
        ! Read data from input file
        CALL inquire_PD_obs_data_file(PD_obs)
        
      END IF ! IF (C%do_benchmark_experiment) THEN
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate memory  
    CALL allocate_shared_dp_1D( PD_obs%nlon,                  PD_obs%lon,     PD_obs%wlon    )
    CALL allocate_shared_dp_1D(              PD_obs%nlat,     PD_obs%lat,     PD_obs%wlat    )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%T2m,     PD_obs%wT2m    )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Precip,  PD_obs%wPrecip )
    CALL allocate_shared_dp_2D( PD_obs%nlon, PD_obs%nlat,     PD_obs%Hs,      PD_obs%wHs     )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Q_TOA,   PD_obs%wQ_TOA  )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Albedo,  PD_obs%wAlbedo )
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Wind_WE, PD_obs%wWind_WE)
    CALL allocate_shared_dp_3D( PD_obs%nlon, PD_obs%nlat, 12, PD_obs%Wind_SN, PD_obs%wWind_SN)
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    IF (par%master) THEN
    
      IF (C%do_benchmark_experiment) THEN
        IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
            C%choice_benchmark_experiment == 'Halfar' .OR. &
            C%choice_benchmark_experiment == 'Bueler' .OR. &
            C%choice_benchmark_experiment == 'MISMIP_mod'.OR. &
            C%choice_benchmark_experiment == 'mesh_generation_test') THEN
          ! No need to read any climate data
        ELSE
          WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_obs_data_fields!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
      ELSE ! IF (C%do_benchmark_experiment) THEN
      
        ! Read data from input file
        CALL read_PD_obs_data_file(PD_obs)
        
      END IF ! IF (C%do_benchmark_experiment) THEN
      
    END IF ! IF (par%master) THEN
    CALL sync
      
    ! Determine process domains
    PD_obs%i1 = MAX(1,           FLOOR(REAL(PD_obs%nlon *  par%i      / par%n)) + 1)
    PD_obs%i2 = MIN(PD_obs%nlon, FLOOR(REAL(PD_obs%nlon * (par%i + 1) / par%n)))
    
  END SUBROUTINE initialise_PD_obs_data_fields  
  SUBROUTINE initialise_snapshot( nlat, nlon, snapshot)
    ! Allocate shared memory and read data for the individual GCM snapshots (global),
    ! stored in the climate matrix.
    
    IMPLICIT NONE
    
    TYPE(type_subclimate_global),   INTENT(INOUT) :: snapshot
    INTEGER,                        INTENT(OUT)   :: nlat, nlon
    
    nlat = 1
    nlon = 1
    
    CALL allocate_shared_dp_0D(                snapshot%CO2,        snapshot%wCO2       )
    CALL allocate_shared_dp_0D(                snapshot%orbit_time, snapshot%worbit_time)
    CALL allocate_shared_dp_0D(                snapshot%orbit_ecc,  snapshot%worbit_ecc )
    CALL allocate_shared_dp_0D(                snapshot%orbit_obl,  snapshot%worbit_obl )
    CALL allocate_shared_dp_0D(                snapshot%orbit_pre,  snapshot%worbit_pre )
    
    CALL allocate_shared_dp_1D(nlon,           snapshot%lon,        snapshot%wlon       )
    CALL allocate_shared_dp_1D(      nlat,     snapshot%lat,        snapshot%wlat       )
    CALL allocate_shared_dp_3D(nlon, nlat, 12, snapshot%Q_TOA,      snapshot%wQ_TOA     )
    CALL allocate_shared_dp_3D(nlon, nlat, 12, snapshot%Albedo,     snapshot%wAlbedo    )
    CALL allocate_shared_dp_3D(nlon, nlat, 12, snapshot%T2m,        snapshot%wT2m       )
    CALL allocate_shared_dp_3D(nlon, nlat, 12, snapshot%Precip,     snapshot%wPrecip    )
    CALL allocate_shared_dp_2D(nlon, nlat,     snapshot%Hs,         snapshot%wHs        )
    CALL allocate_shared_dp_3D(nlon, nlat, 12, snapshot%Wind_WE,    snapshot%wWind_WE   )
    CALL allocate_shared_dp_3D(nlon, nlat, 12, snapshot%Wind_SN,    snapshot%wWind_SN   )  
    CALL allocate_shared_dp_3D(nlon, nlat, 12, snapshot%Albedo,     snapshot%wAlbedo    )
  
  END SUBROUTINE initialise_snapshot
  
  ! Remapping
  SUBROUTINE remap_climate_model( mesh_old, mesh_new, map, climate, matrix)
    ! Remap or reallocate all the data fields
  
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    TYPE(type_climate_model),            INTENT(INOUT) :: climate
    TYPE(type_climate_matrix),           INTENT(IN)    :: matrix
        
    ! All applied climate fields are remapped, since they are not updated in every timestep.
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, climate%applied%T2m,     climate%applied%wT2m,     'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, climate%applied%Precip,  climate%applied%wPrecip,  'trilin'        )
    CALL Remap_field_dp(         mesh_old, mesh_new, map, climate%applied%Hs,      climate%applied%wHs,      'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, climate%applied%Q_TOA,   climate%applied%wQ_TOA,   'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, climate%applied%Albedo,  climate%applied%wAlbedo,  'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, climate%applied%Wind_WE, climate%applied%wWind_WE, 'trilin'        )
    CALL Remap_field_dp_monthly( mesh_old, mesh_new, map, climate%applied%Wind_SN, climate%applied%wWind_SN, 'trilin'        )
    
    ! The subclimates are reallocated, then remapped from the global fields.
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%PD_obs%T2m,         climate%PD_obs%wT2m       , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%PD_obs%Precip,      climate%PD_obs%wPrecip    , 12)
    CALL Reallocate_field_dp(    mesh_new%nV,  climate%PD_obs%Hs,          climate%PD_obs%wHs            )
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%PD_obs%Q_TOA,       climate%PD_obs%wQ_TOA     , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%PD_obs%Albedo,      climate%PD_obs%wAlbedo    , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%PD_obs%Wind_WE,     climate%PD_obs%wWind_WE   , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%PD_obs%Wind_SN,     climate%PD_obs%wWind_SN   , 12)
    
    CALL map_subclimate_to_mesh( mesh_new,  matrix%PD_obs, climate%PD_obs)
    
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_PI%T2m,         climate%GCM_PI%wT2m       , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_PI%Precip,      climate%GCM_PI%wPrecip    , 12)
    CALL Reallocate_field_dp(    mesh_new%nV,  climate%GCM_PI%Hs,          climate%GCM_PI%wHs            )
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_PI%Q_TOA,       climate%GCM_PI%wQ_TOA     , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_PI%Albedo,      climate%GCM_PI%wAlbedo    , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_PI%Wind_WE,     climate%GCM_PI%wWind_WE   , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_PI%Wind_SN,     climate%GCM_PI%wWind_SN   , 12)
    
    !CALL map_subclimate_to_mesh( mesh_new,  matrix%GCM_PI, climate%GCM_PI)
    
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_LGM%T2m,        climate%GCM_LGM%wT2m      , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_LGM%Precip,     climate%GCM_LGM%wPrecip   , 12)
    CALL Reallocate_field_dp(    mesh_new%nV,  climate%GCM_LGM%Hs,         climate%GCM_LGM%wHs           )
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_LGM%Q_TOA,      climate%GCM_LGM%wQ_TOA    , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_LGM%Albedo,     climate%GCM_LGM%wAlbedo   , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_LGM%Wind_WE,    climate%GCM_LGM%wWind_WE  , 12)
    CALL Reallocate_field_dp_3D( mesh_new%nV,  climate%GCM_LGM%Wind_SN,    climate%GCM_LGM%wWind_SN  , 12)
    
    !CALL map_subclimate_to_mesh( mesh_new,  matrix%GCM_LGM, climate%GCM_LGM)
    
  END SUBROUTINE remap_climate_model

END MODULE climate_module
