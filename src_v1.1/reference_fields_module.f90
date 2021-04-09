MODULE reference_fields_module

  USE mpi
  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, &
                                         allocate_shared_int_0D, allocate_shared_dp_0D, &
                                         allocate_shared_int_1D, allocate_shared_dp_1D, &
                                         allocate_shared_int_2D, allocate_shared_dp_2D, &
                                         allocate_shared_int_3D, allocate_shared_dp_3D, &
                                         deallocate_shared
  USE data_types_module,           ONLY: type_mesh, type_init_data_fields, type_PD_data_fields
  USE parameters_module,           ONLY: seawater_density, ice_density
  USE netcdf_module,               ONLY: debug, write_to_debug_file, inquire_PD_data_file, inquire_init_data_file, &
                                         read_PD_data_file, read_init_data_file
  USE mesh_help_functions_module,  ONLY: partition_list
  USE mesh_mapping_module,         ONLY: create_remapping_arrays_mesh_grid, map_grid2mesh_2D, deallocate_remapping_arrays_mesh_grid

  IMPLICIT NONE
  
CONTAINS

  ! == Initialise init and PD reference data on their supplied grids
  SUBROUTINE initialise_PD_data_fields( PD, region_name)
    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
    INTEGER                                       :: i,j
    
    IF (C%do_benchmark_experiment) THEN
      ! Benchmark experiments: don't read data from an input file, but generate it internally
      CALL initialise_PD_data_fields_benchmarks( PD)
    ELSE
      ! Realistic experiments: read data from an input file
    
      IF      (region_name == 'NAM') THEN
        PD%netcdf%filename   = C%filename_PD_NAM
      ELSE IF (region_name == 'EAS') THEN
        PD%netcdf%filename   = C%filename_PD_EAS
      ELSE IF (region_name == 'GRL') THEN
        PD%netcdf%filename   = C%filename_PD_GRL
      ELSE IF (region_name == 'ANT') THEN
        PD%netcdf%filename   = C%filename_PD_ANT
      END IF
      
      IF (par%master) WRITE(0,*) '  Reading PD   reference data from file "', TRIM(PD%netcdf%filename), '"'
    
      ! Determine size of Cartesian grid, so that memory can be allocated accordingly.
      CALL allocate_shared_int_0D( PD%grid%nx,   PD%grid%wnx  )
      CALL allocate_shared_int_0D( PD%grid%ny,   PD%grid%wny  )
      
      IF (par%master) CALL inquire_PD_data_file(PD)
      CALL sync
    
      ! Allocate memory 
      CALL allocate_shared_dp_0D(                         PD%grid%dx,   PD%grid%wdx  )
      CALL allocate_shared_dp_0D(                         PD%grid%xmin, PD%grid%wxmin)
      CALL allocate_shared_dp_0D(                         PD%grid%xmax, PD%grid%wxmax)
      CALL allocate_shared_dp_0D(                         PD%grid%ymin, PD%grid%wymin)
      CALL allocate_shared_dp_0D(                         PD%grid%ymax, PD%grid%wymax)
      CALL allocate_shared_dp_1D( PD%grid%nx,             PD%grid%x,    PD%grid%wx   )
      CALL allocate_shared_dp_1D( PD%grid%ny,             PD%grid%y,    PD%grid%wy   )
      
      CALL allocate_shared_dp_2D( PD%grid%nx, PD%grid%ny, PD%Hi_grid,   PD%wHi_grid  )
      CALL allocate_shared_dp_2D( PD%grid%nx, PD%grid%ny, PD%Hb_grid,   PD%wHb_grid  )
      CALL allocate_shared_dp_2D( PD%grid%nx, PD%grid%ny, PD%Hs_grid,   PD%wHs_grid  )
      CALL allocate_shared_int_2D(PD%grid%nx, PD%grid%ny, PD%mask_grid, PD%wmask_grid)
      
      ! Read data from input file
      IF (par%master) CALL read_PD_data_file(PD)
      CALL sync
      
    END IF ! IF (C%do_benchmark_experiment) THEN
      
    ! Determine vertex domains
    CALL partition_list( PD%grid%nx, par%i, par%n, PD%grid%i1, PD%grid%i2)
    CALL partition_list( PD%grid%ny, par%i, par%n, PD%grid%j1, PD%grid%j2)
    
    ! Determine the masks
  
    DO i = PD%grid%i1, PD%grid%i2
    DO j = 1, PD%grid%ny
    
      PD%mask_grid( i,j) = 0
      
      IF (PD%Hi_grid(i,j) <= 0._dp) THEN
        IF (PD%Hb_grid(i,j) < 0._dp) THEN
          PD%mask_grid(i,j) = C%type_ocean
        ELSE
          PD%mask_grid(i,j) = C%type_land
        END IF
      ELSE
        ! Either sheet or shelf, depending on flotation criterion
        IF (PD%Hi_grid(i,j) > (-1._dp * PD%Hb_grid(i,j)) * seawater_density / ice_density) THEN
          PD%mask_grid(i,j) = C%type_sheet
        ELSE
          PD%mask_grid(i,j) = C%type_shelf
        END IF
      END IF
      
    END DO
    END DO
    CALL sync
  
    DO i = MAX(2,PD%grid%i1), MIN(PD%grid%nx-1,PD%grid%i2)
    DO j = 2, PD%grid%ny-1
    
      IF (PD%mask_grid(i,j) == C%type_sheet) THEN
      IF (PD%mask_grid(i-1,j) == C%type_shelf .OR. &
          PD%mask_grid(i+1,j) == C%type_shelf .OR. &
          PD%mask_grid(i,j-1) == C%type_shelf .OR. &
          PD%mask_grid(i,j+1) == C%type_shelf) THEN
        PD%mask_grid(i,j) = C%type_groundingline
      END IF
      END IF
  
      IF (PD%mask_grid(i,j) == C%type_sheet .OR. PD%mask_grid(i,j) == C%type_shelf) THEN
      IF ((PD%mask_grid(i-1,j) == C%type_ocean) .OR. &
          (PD%mask_grid(i+1,j) == C%type_ocean) .OR. &
          (PD%mask_grid(i,j-1) == C%type_ocean) .OR. &
          (PD%mask_grid(i,j+1) == C%type_ocean)) THEN
        PD%mask_grid(i,j) = C%type_calvingfront
      END IF
      END IF
      
    END DO 
    END DO ! DO i = MAX(2,PD%i1), MIN(PD%nx-1,PD%i2)
    CALL sync
    
  END SUBROUTINE initialise_PD_data_fields
  SUBROUTINE initialise_PD_data_fields_benchmarks( PD)
    ! Initialising PD reference data for the different benchmark experiments
    ! (needed because during mesh updates, bedrock is updated from the PD
    ! fields rather than from the init fields)
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    
    ! Local variables:
    INTEGER                                       :: i, j, cerr, ierr
    REAL(dp)                                      :: R
    
    REAL(dp), PARAMETER                           :: EISMINT_xmin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_xmax =  750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymax =  750000._dp
    
    REAL(dp), PARAMETER                           :: MISMIP_xmin = -1800000._dp
    REAL(dp), PARAMETER                           :: MISMIP_xmax =  1800000._dp
    REAL(dp), PARAMETER                           :: MISMIP_ymin = -1800000._dp
    REAL(dp), PARAMETER                           :: MISMIP_ymax =  1800000._dp
    
    ! Determine size of Cartesian grid, so that memory can be allocated accordingly.
    CALL allocate_shared_int_0D( PD%grid%nx, PD%grid%wnx)
    CALL allocate_shared_int_0D( PD%grid%ny, PD%grid%wny)
    CALL allocate_shared_dp_0D(  PD%grid%dx, PD%grid%wdx)
    
    IF (par%master) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler') THEN
        PD%grid%nx = 51
        PD%grid%ny = 51
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        PD%grid%nx = 900
        PD%grid%ny = 900
      ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        PD%grid%nx = 500
        PD%grid%ny = 500
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_data_fields!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate memory  
    CALL allocate_shared_dp_1D( PD%grid%nx,             PD%grid%x,    PD%grid%wx   )
    CALL allocate_shared_dp_1D( PD%grid%ny,             PD%grid%y,    PD%grid%wy   )
    CALL allocate_shared_dp_2D( PD%grid%nx, PD%grid%ny, PD%Hi_grid,   PD%wHi_grid  )
    CALL allocate_shared_dp_2D( PD%grid%nx, PD%grid%ny, PD%Hb_grid,   PD%wHb_grid  )
    CALL allocate_shared_dp_2D( PD%grid%nx, PD%grid%ny, PD%Hs_grid,   PD%wHs_grid  )
    CALL allocate_shared_int_2D(PD%grid%nx, PD%grid%ny, PD%mask_grid, PD%wmask_grid)
    
    IF (par%master) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
          C%choice_benchmark_experiment == 'Halfar' .OR. &
          C%choice_benchmark_experiment == 'Bueler') THEN
    
        ! Simple square grid with no data
        DO i = 1, PD%grid%nx
          PD%grid%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (PD%grid%nx-1)
        END DO
        DO j = 1, PD%grid%ny
          PD%grid%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (PD%grid%ny-1)
        END DO
        
        PD%grid%dx = PD%grid%x(2) - PD%grid%x(1)
        
        PD%Hi_grid = 0._dp
        PD%Hb_grid = 0._dp
        PD%Hs_grid = 0._dp
        
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        ! Bedrock slope as described by Pattyn et al., 2012
      
        DO i = 1, PD%grid%nx
        DO j = 1, PD%grid%ny
          PD%grid%x(i) = MISMIP_xmin + (MISMIP_xmax - MISMIP_xmin) * (i-1) / (PD%grid%nx-1)
          PD%grid%y(j) = MISMIP_ymin + (MISMIP_ymax - MISMIP_ymin) * (j-1) / (PD%grid%ny-1)
          R = NORM2([ PD%grid%x(i), PD%grid%y(j)])
          PD%Hb_grid( i,j) = 720._dp - 778.5_dp * (R / 750000._dp)
          PD%Hs_grid( i,j) = MAX( 0._dp, PD%Hb_grid( i,j))
        END DO
        END DO
        
        PD%grid%dx = PD%grid%x(2) - PD%grid%x(1)
        
        PD%Hi_grid = 0._dp
        
      ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        ! A small circular island with a little bay on the southeast side
      
        DO i = 1, PD%grid%nx
        DO j = 1, PD%grid%ny
          PD%grid%x(i) = -1000000._dp + 2000000._dp * (i-1) / (PD%grid%nx-1)
          PD%grid%y(j) = -1000000._dp + 2000000._dp * (j-1) / (PD%grid%ny-1)
          R = NORM2([ PD%grid%x( i), PD%grid%y( j)]) / 1000._dp
          PD%Hb_grid( i,j) = -450._dp - (400._dp * ATAN(( R - 600._dp) / 180._dp))
          R = NORM2([ (PD%grid%x(i) - 200000._dp), (PD%grid%y(j) + 100000._dp)]) / 1000._dp
          PD%Hb_grid( i,j) = PD%Hb_grid( i,j) + (350._dp * ATAN(( R - 240._dp) / 160._dp))
          PD%Hs_grid( i,j) = MAX( 0._dp, PD%Hb_grid( i,j))
        END DO
        END DO
        
        PD%grid%dx = PD%grid%x(2) - PD%grid%x(1)
        
        PD%Hi_grid = 0._dp
        
      ELSE
      
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_data_fields!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
          
      END IF ! IF (C%choice_benchmark_experiment == ... )
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE initialise_PD_data_fields_benchmarks
  
  SUBROUTINE initialise_init_data_fields( init, region_name)
    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
    INTEGER                                       :: i,j
    REAL(dp)                                      :: dx, dy, d2dx2, d2dy2, d2dxdy
    
    IF (C%do_benchmark_experiment) THEN
      ! Benchmark experiments: don't read data from an input file, but generate it internally
      CALL initialise_init_data_fields_benchmarks( init)
    ELSE
      ! Realistic experiments: read data from an input file
    
      IF      (region_name == 'NAM') THEN
        init%netcdf%filename   = C%filename_init_NAM
      ELSE IF (region_name == 'EAS') THEN
        init%netcdf%filename   = C%filename_init_EAS
      ELSE IF (region_name == 'GRL') THEN
        init%netcdf%filename   = C%filename_init_GRL
      ELSE IF (region_name == 'ANT') THEN
        init%netcdf%filename   = C%filename_init_ANT
      END IF
      
      IF (par%master) WRITE(0,*) '  Reading init reference data from file "', TRIM(init%netcdf%filename), '"'  
      
      ! Determine size of Cartesian grid, so that memory can be allocated accordingly.
      CALL allocate_shared_int_0D( init%grid%nx, init%grid%wnx)
      CALL allocate_shared_int_0D( init%grid%ny, init%grid%wny)
      IF (par%master) CALL inquire_init_data_file( init)
      CALL sync
    
      ! Allocate memory - init    
      CALL allocate_shared_dp_0D(                              init%grid%dx,                init%grid%wdx               )
      CALL allocate_shared_dp_0D(                              init%grid%xmin,              init%grid%wxmin             )
      CALL allocate_shared_dp_0D(                              init%grid%xmax,              init%grid%wxmax             )
      CALL allocate_shared_dp_0D(                              init%grid%ymin,              init%grid%wymin             )
      CALL allocate_shared_dp_0D(                              init%grid%ymax,              init%grid%wymax             )
      CALL allocate_shared_dp_1D(  init%grid%nx,               init%grid%x,                 init%grid%wx                )
      CALL allocate_shared_dp_1D(  init%grid%ny,               init%grid%y,                 init%grid%wy                )
      
      CALL allocate_shared_dp_2D(  init%grid%nx, init%grid%ny, init%Hi_grid,                init%wHi_grid               )
      CALL allocate_shared_dp_2D(  init%grid%nx, init%grid%ny, init%Hb_grid,                init%wHb_grid               )
      CALL allocate_shared_dp_2D(  init%grid%nx, init%grid%ny, init%Hs_grid,                init%wHs_grid               )
      
      CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_grid,              init%wmask_grid             )
      CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_ice_grid,          init%wmask_ice_grid         )
      CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_gl_grid,           init%wmask_gl_grid          )
      CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_cf_grid,           init%wmask_cf_grid          )
      CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_coast_grid,        init%wmask_coast_grid       )
      CALL allocate_shared_dp_2D(  init%grid%nx, init%grid%ny, init%surface_curvature_grid, init%wsurface_curvature_grid)
    
      ! Read data
      IF (par%master) CALL read_init_data_file(init) 
      
    END IF ! IF (par%master) THEN
    CALL sync
      
    ! Determine vertex domains
    CALL partition_list( init%grid%nx, par%i, par%n, init%grid%i1, init%grid%i2)
    CALL partition_list( init%grid%ny, par%i, par%n, init%grid%j1, init%grid%j2)
    
    ! Surface curvature
    init%surface_curvature_grid(init%grid%i1:init%grid%i2,:) = 0._dp
    dx = ABS(init%grid%x(2) - init%grid%x(1))
    dy = ABS(init%grid%y(2) - init%grid%y(1))
    DO i = MAX(2,init%grid%i1), MIN(init%grid%nx-1,init%grid%i2)
    DO j = 2, init%grid%ny-1
      d2dx2  = (init%Hs_grid(i+1,j  ) + init%Hs_grid(i-1,j  ) - 2._dp*init%Hs_grid(i,j)) / (dx**2)
      d2dy2  = (init%Hs_grid(i  ,j+1) + init%Hs_grid(i  ,j-1) - 2._dp*init%Hs_grid(i,j)) / (dy**2)
      d2dxdy = (init%Hs_grid(i+1,j+1) + init%Hs_grid(i-1,j-1) - init%Hs_grid(i+1,j-1) - init%Hs_grid(i-1,j+1)) / (4*dx*dy)
      init%surface_curvature_grid(i,j) = MAX(-1E-6, MIN(1E-6, SQRT(d2dx2**2 + d2dy2**2 + d2dxdy**2)))
    END DO
    END DO
    
    ! Determine the masks
    init%mask_grid(       init%grid%i1:init%grid%i2,:) = 0
    init%mask_ice_grid(   init%grid%i1:init%grid%i2,:) = 0
    init%mask_gl_grid(    init%grid%i1:init%grid%i2,:) = 0
    init%mask_cf_grid(    init%grid%i1:init%grid%i2,:) = 0
    init%mask_coast_grid( init%grid%i1:init%grid%i2,:) = 0
  
    DO i = init%grid%i1, init%grid%i2
    DO j = 1, init%grid%ny
      IF (init%Hi_grid(i,j) == 0._dp) THEN
        IF (init%Hb_grid(i,j) < 0._dp) THEN
          init%mask_grid(i,j) = C%type_ocean
        ELSE
          init%mask_grid(i,j) = C%type_land
        END IF
      ELSE
        ! Either sheet or shelf, depending on flotation criterion
        init%mask_ice_grid(i,j) = 1
        IF (init%Hi_grid(i,j) > (-1._dp * init%Hb_grid(i,j)) * seawater_density / ice_density) THEN
          init%mask_grid(i,j) = C%type_sheet
        ELSE
          init%mask_grid(i,j) = C%type_shelf
        END IF
      END IF
    END DO
    END DO
    CALL sync
  
    DO i = MAX(2,init%grid%i1), MIN(init%grid%nx-1,init%grid%i2)
    DO j = 2, init%grid%ny-1
      IF (init%mask_grid(i,j) == C%type_sheet) THEN
      IF (init%mask_grid(i-1,j) == C%type_shelf .OR. &
          init%mask_grid(i+1,j) == C%type_shelf .OR. &
          init%mask_grid(i,j-1) == C%type_shelf .OR. &
          init%mask_grid(i,j+1) == C%type_shelf) THEN
        init%mask_gl_grid(i,j) = 1
        init%mask_grid(i,j) = C%type_groundingline
      END IF
      END IF
  
      IF (init%mask_grid(i,j) == C%type_sheet .OR. init%mask_grid(i,j) == C%type_shelf) THEN
      IF ((init%mask_grid(i-1,j) == C%type_ocean) .OR. &
          (init%mask_grid(i+1,j) == C%type_ocean) .OR. &
          (init%mask_grid(i,j-1) == C%type_ocean) .OR. &
          (init%mask_grid(i,j+1) == C%type_ocean)) THEN
        init%mask_cf_grid(i,j) = 1
        init%mask_grid(i,j) = C%type_calvingfront
      END IF
      END IF
      
      IF (init%mask_grid(i,j) == C%type_land) THEN
      IF ((init%mask_grid(i-1,j) == C%type_ocean) .OR. &
          (init%mask_grid(i+1,j) == C%type_ocean) .OR. &
          (init%mask_grid(i,j-1) == C%type_ocean) .OR. &
          (init%mask_grid(i,j+1) == C%type_ocean)) THEN
        init%mask_coast_grid(i,j) = 1
      END IF
      END IF
    END DO
    END DO
    CALL sync
    
  END SUBROUTINE initialise_init_data_fields
  SUBROUTINE initialise_init_data_fields_benchmarks( init)
    ! Initialising init reference data for the different benchmark experiments
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    
    ! Local variables:
    INTEGER                                       :: i, j, cerr, ierr
    REAL(dp)                                      :: R
    
    REAL(dp), PARAMETER                           :: EISMINT_xmin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_xmax =  750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymax =  750000._dp
    
    REAL(dp), PARAMETER                           :: MISMIP_xmin = -1800000._dp
    REAL(dp), PARAMETER                           :: MISMIP_xmax =  1800000._dp
    REAL(dp), PARAMETER                           :: MISMIP_ymin = -1800000._dp
    REAL(dp), PARAMETER                           :: MISMIP_ymax =  1800000._dp
    
    ! Determine grid size
    ! ===================
    
    CALL allocate_shared_int_0D( init%grid%nx, init%grid%wnx)
    CALL allocate_shared_int_0D( init%grid%ny, init%grid%wny)
    CALL allocate_shared_dp_0D(  init%grid%dx, init%grid%wdx)
    
    IF (par%master) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6') THEN
        init%grid%nx = 51
        init%grid%ny = 51
      ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
              C%choice_benchmark_experiment == 'Bueler') THEN
        init%grid%nx = 2501
        init%grid%ny = 2501
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        init%grid%nx = 900
        init%grid%ny = 900
      ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        init%grid%nx = 500
        init%grid%ny = 500
      ELSE
        IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_init_data_fields_cart!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate memory
    ! ===============
    
    CALL allocate_shared_dp_1D(  init%grid%nx,               init%grid%x,                 init%grid%wx                )
    CALL allocate_shared_dp_1D(  init%grid%ny,               init%grid%y,                 init%grid%wy                )
    CALL allocate_shared_dp_2D(  init%grid%nx, init%grid%ny, init%Hi_grid,                init%wHi_grid               )
    CALL allocate_shared_dp_2D(  init%grid%nx, init%grid%ny, init%Hb_grid,                init%wHb_grid               )
    CALL allocate_shared_dp_2D(  init%grid%nx, init%grid%ny, init%Hs_grid,                init%wHs_grid               )
    
    CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_grid,              init%wmask_grid             )
    CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_ice_grid,          init%wmask_ice_grid         )
    CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_gl_grid,           init%wmask_gl_grid          )
    CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_cf_grid,           init%wmask_cf_grid          )
    CALL allocate_shared_int_2D( init%grid%nx, init%grid%ny, init%mask_coast_grid,        init%wmask_coast_grid       )
    CALL allocate_shared_dp_2D(  init%grid%nx, init%grid%ny, init%surface_curvature_grid, init%wsurface_curvature_grid)
    
    ! Fill in data
    ! ============
    
    IF (par%master) THEN
      IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
              C%choice_benchmark_experiment == 'EISMINT_6') THEN
              
        ! Simple square grid with no data
        DO i = 1, init%grid%nx
          init%grid%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (init%grid%nx-1)
        END DO
        DO j = 1, init%grid%ny
          init%grid%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (init%grid%ny-1)
        END DO
        
        init%grid%dx = init%grid%x(2) - init%grid%x(1)
        
        init%Hi_grid = 0._dp
        init%Hb_grid = 0._dp
        init%Hs_grid = 0._dp
        
      ELSEIF (C%choice_benchmark_experiment == 'Halfar') THEN
      
        ! Start with the Halfar solution for the given parameters
        DO i = 1, init%grid%nx
        DO j = 1, init%grid%ny
          init%grid%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (init%grid%nx-1)
          init%grid%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (init%grid%ny-1)
          init%Hi_grid(i,j) = Halfar_solution( C%halfar_solution_H0, C%halfar_solution_R0, init%grid%x(i), init%grid%y(j), 0._dp)
        END DO
        END DO
        init%Hs_grid = init%Hi_grid
        
        init%grid%dx = init%grid%x(2) - init%grid%x(1)
      
      ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
      
        ! Start with the Bueler solution for the given parameters
        DO i = 1, init%grid%nx
        DO j = 1, init%grid%ny
          init%grid%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (init%grid%nx-1)
          init%grid%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (init%grid%ny-1)
          init%Hi_grid(i,j) = Bueler_solution( C%halfar_solution_H0, C%halfar_solution_R0, C%bueler_solution_lambda, init%grid%x(i), init%grid%y(j), C%start_time_of_run)
        END DO
        END DO
        init%Hs_grid = init%Hi_grid
        
        init%grid%dx = init%grid%x(2) - init%grid%x(1)
        
      ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
        ! Bedrock slope as described by Pattyn et al., 2012
      
        DO i = 1, init%grid%nx
        DO j = 1, init%grid%ny
          init%grid%x(i) = MISMIP_xmin + (MISMIP_xmax - MISMIP_xmin) * (i-1) / (init%grid%nx-1)
          init%grid%y(j) = MISMIP_ymin + (MISMIP_ymax - MISMIP_ymin) * (j-1) / (init%grid%ny-1)
          R = NORM2([ init%grid%x(i), init%grid%y(j)])
          init%Hb_grid(i,j) = 720._dp - 778.5_dp * (R / 750000._dp)
          init%Hs_grid(i,j) = MAX( 0._dp, init%Hb_grid(i,j))
        END DO
        END DO
        
        init%grid%dx = init%grid%x(2) - init%grid%x(1)
        
        init%Hi_grid = 100._dp
        
      ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
        ! A small circular island with a little bay on the southeast side
      
        DO i = 1, init%grid%nx
        DO j = 1, init%grid%ny
          init%grid%x(i) = -1000000._dp + 2000000._dp * (i-1) / (init%grid%nx-1)
          init%grid%y(j) = -1000000._dp + 2000000._dp * (j-1) / (init%grid%ny-1)
          R = NORM2([ init%grid%x(i), init%grid%y(j)]) / 1000._dp
          init%Hb_grid(i,j) = -450._dp - (400._dp * ATAN(( R - 600._dp) / 180._dp))
          R = NORM2([ (init%grid%x(i) - 200000._dp), (init%grid%y(j) + 100000._dp)]) / 1000._dp
          init%Hb_grid(i,j) = init%Hb_grid(i,j) + (350._dp * ATAN(( R - 240._dp) / 160._dp))
          init%Hs_grid(i,j) = MAX( 0._dp, init%Hb_grid(i,j))
        END DO
        END DO
        
        init%grid%dx = init%grid%x(2) - init%grid%x(1)
      
      ELSE
      
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_init_data_fields_cart!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      
      END IF
    END IF ! IF (par%master) THEN
    CALL sync
    
  END SUBROUTINE initialise_init_data_fields_benchmarks
  SUBROUTINE deallocate_init_data( init)
    
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    
    CALL deallocate_shared( init%grid%wnx)
    CALL deallocate_shared( init%grid%wny)
    CALL deallocate_shared( init%grid%wx)
    CALL deallocate_shared( init%grid%wy)
    CALL deallocate_shared( init%wHi_grid)
    CALL deallocate_shared( init%wHb_grid)
    CALL deallocate_shared( init%wHs_grid)
    CALL deallocate_shared( init%wsurface_curvature_grid)
    CALL deallocate_shared( init%wmask_grid)
    CALL deallocate_shared( init%wmask_ice_grid)
    CALL deallocate_shared( init%wmask_gl_grid)
    CALL deallocate_shared( init%wmask_cf_grid)
    CALL deallocate_shared( init%wmask_coast_grid)
    CALL deallocate_shared( init%wHi)
    CALL deallocate_shared( init%wHb)
    CALL deallocate_shared( init%wHs)
    
  END SUBROUTINE deallocate_init_data
  
  ! == Map init and PD references data from their supplied grids to the model mesh
  SUBROUTINE map_PD_data_to_mesh(   mesh, PD)
    ! Calculate the mapping arrays, map data to the mesh. Allocate and deallocate memory for the mapping arrays within this subroutine.
    
    IMPLICIT NONE
  
    ! Input and output variables
    TYPE(type_mesh),                INTENT(INOUT) :: mesh
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    
    IF (par%master) WRITE(0,*) '  Mapping PD   reference data to the mesh...'
    
    ! Calculate mapping arrays
    CALL create_remapping_arrays_mesh_grid( mesh, PD%grid)
    
    ! Allocate memory for reference data on the mesh
    CALL allocate_shared_dp_1D( mesh%nV, PD%Hi, PD%wHi)
    CALL allocate_shared_dp_1D( mesh%nV, PD%Hb, PD%wHb)
    CALL allocate_shared_dp_1D( mesh%nV, PD%Hs, PD%wHs)
    
    ! Map PD data to the mesh
    CALL map_grid2mesh_2D( mesh, PD%grid, PD%Hi_grid, PD%Hi)
    CALL map_grid2mesh_2D( mesh, PD%grid, PD%Hb_grid, PD%Hb)
    PD%Hs( mesh%v1:mesh%v2) = MAX(0._dp, PD%Hb( mesh%v1:mesh%v2) + PD%Hi( mesh%v1:mesh%v2))
    
    ! Deallocate remapping arrays
    CALL deallocate_remapping_arrays_mesh_grid( PD%grid)
  
  END SUBROUTINE map_PD_data_to_mesh
  SUBROUTINE map_init_data_to_mesh( mesh, init)
    ! Calculate the mapping arrays, map data to the mesh. Allocate and deallocate memory for the mapping arrays within this subroutine.
    
    IMPLICIT NONE
  
    ! Input and output variables
    TYPE(type_mesh),                INTENT(INOUT) :: mesh
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    
    ! Local variables
    INTEGER                                       :: vi
    
    IF (par%master) WRITE(0,*) '  Mapping init reference data to the mesh...'
    
    ! Calculate mapping arrays
    CALL create_remapping_arrays_mesh_grid( mesh, init%grid)
    
    ! Allocate memory for reference data on the mesh
    CALL allocate_shared_dp_1D( mesh%nV, init%Hi, init%wHi)
    CALL allocate_shared_dp_1D( mesh%nV, init%Hb, init%wHb)
    CALL allocate_shared_dp_1D( mesh%nV, init%Hs, init%wHs)
    
    ! Map init data to the mesh
    CALL map_grid2mesh_2D( mesh, init%grid, init%Hi_grid, init%Hi)
    CALL map_grid2mesh_2D( mesh, init%grid, init%Hb_grid, init%Hb)
    init%Hs( mesh%v1:mesh%v2) = MAX(0._dp, init%Hb( mesh%v1:mesh%v2) + init%Hi( mesh%v1:mesh%v2))  
    
    ! For the Bueler benchmark experiment, recalculate initial ice thickness on the mesh,
    ! to make sure we have the correct starting conditions.
    IF (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'Bueler') THEN
      DO vi = mesh%v1, mesh%v2
        init%Hi(vi) = Bueler_solution( C%halfar_solution_H0, C%halfar_solution_R0, C%bueler_solution_lambda, mesh%V(vi,1), mesh%V(vi,2), C%start_time_of_run)
      END DO      
      init%Hs( mesh%v1:mesh%v2) = init%Hi( mesh%v1:mesh%v2)
      CALL sync
    END IF
    
    ! Deallocate remapping arrays
    CALL deallocate_remapping_arrays_mesh_grid( init%grid)
  
  END SUBROUTINE map_init_data_to_mesh
  
  ! == The Halfar and Bueler analytical solutions, used to initialise those benchmark experiments
  FUNCTION Halfar_solution( H0, R0, x, y, t) RESULT(H)
    ! Describes an ice-sheet at time t (in years) conforming to the Halfar similarity function 
    ! with dome thickness H0 and margin radius R0 at t0. Used to initialise the model
    ! for the Halfar solution test run
    
    USE parameters_module, ONLY: sec_per_year
    
    IMPLICIT NONE
    
    ! Input variables
    REAL(dp), INTENT(IN) :: H0 ! Ice dome thickness at t=0 [m]
    REAL(dp), INTENT(IN) :: R0 ! Ice margin radius  at t=0 [m]
    REAL(dp), INTENT(IN) :: x  ! x coordinate [m]
    REAL(dp), INTENT(IN) :: y  ! y coordinate [m]
    REAL(dp), INTENT(IN) :: t  ! Time from t0 [years]
    
    ! Result
    REAL(dp)             :: H  ! Ice thickness at [x,y] at t=0 [m]
    
    ! Local variables
    REAL(dp) :: A_flow, rho, g, Gamma, t0, r, f1, f2, f3, tp
    
    A_flow  = 1E-16_dp
    rho     = 910._dp
    g       = 9.81_dp
  
    Gamma = (2._dp / 5._dp) * (A_flow / sec_per_year) * (rho * g)**3._dp
    t0 = 1._dp / (18._dp * Gamma) * (7._dp/4._dp)**3._dp * (R0**4._dp)/(H0**7._dp)
  
    tp = (t * sec_per_year) + t0
  
    r = SQRT(x**2._dp + y**2._dp)
  
    f1 = (t0/tp)**(1._dp/9._dp)
    f2 = (t0/tp)**(1._dp/18._dp)
    f3 = (r/R0)
  
    H = H0 * f1 * MAX(0._dp, (1._dp - (f2*f3)**(4._dp/3._dp)))**(3._dp/7._dp)
  
  END FUNCTION Halfar_solution
  FUNCTION Bueler_solution( H0, R0, lambda, x, y, t) RESULT(H)
    ! Describes an ice-sheet at time t (in years) conforming to the Bueler solution
    ! with dome thickness H0 and margin radius R0 at t0, with a surface mass balance
    ! determined by lambda. Used to intialise the model for the Bueler solution test run
    
    USE parameters_module, ONLY: sec_per_year
    
    IMPLICIT NONE
    
    ! Input variables
    REAL(dp), INTENT(IN) :: H0      ! Ice dome thickness at t=0 [m]
    REAL(dp), INTENT(IN) :: R0      ! Ice margin radius  at t=0 [m]
    REAL(dp), INTENT(IN) :: lambda  ! Mass balance parameter
    REAL(dp), INTENT(IN) :: x       ! x coordinate [m]
    REAL(dp), INTENT(IN) :: y       ! y coordinate [m]
    REAL(dp), INTENT(IN) :: t       ! Time from t0 [years]
    
    ! Result
    REAL(dp)             :: H  ! Ice thickness at [x,y] at t=0 [m]
    
    ! Local variables
    REAL(dp) :: A_flow, rho, g, n, alpha, beta, Gamma, f1, f2, t0, tp, f3, f4
  
    A_flow  = 1E-16_dp
    rho     = 910._dp
    g       = 9.81_dp
    n       = 3._dp
    
    alpha = (2._dp - (n+1._dp)*lambda) / ((5._dp*n)+3._dp)
    beta  = (1._dp + ((2._dp*n)+1._dp)*lambda) / ((5._dp*n)+3._dp)
    Gamma = 2._dp/5._dp * (A_flow/sec_per_year) * (rho * g)**n
    
    f1 = ((2._dp*n)+1)/(n+1._dp)
    f2 = (R0**(n+1._dp))/(H0**((2._dp*n)+1._dp))
    t0 = (beta / Gamma) * (f1**n) * f2 
    
    !tp = (t * sec_per_year) + t0; % Acutal equation needs t in seconds from zero , but we want to supply t in years from t0
    tp = t * sec_per_year
    
    f1 = (tp / t0)**(-alpha)
    f2 = (tp / t0)**(-beta)
    f3 = SQRT( (x**2._dp) + (y**2._dp) )/R0
    f4 = MAX(0._dp, 1._dp - (f2*f3)**((n+1._dp)/n))
    H = H0 * f1 * f4**(n/((2._dp*n)+1._dp))
    
    !M = (lambda / tp) * H * sec_per_year
  
  END FUNCTION Bueler_solution

END MODULE reference_fields_module
