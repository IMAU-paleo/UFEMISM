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
  USE netcdf_module,               ONLY: debug, write_to_debug_file, inquire_PD_data_file, inquire_init_data_file_cart, &
                                         read_PD_data_file, read_init_data_file_cart!, inquire_init_data_file_mesh, read_init_data_file_mesh
  USE mesh_help_functions_module,  ONLY: partition_list
  USE mesh_mapping_module,         ONLY: get_cart_to_mesh_map, map_cart_to_mesh

  IMPLICIT NONE
  
CONTAINS

  SUBROUTINE initialise_PD_data_fields(        PD,   region_name)
    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
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

    INTEGER, PARAMETER                            :: type_land           = 0
    INTEGER, PARAMETER                            :: type_ocean          = 1
    INTEGER, PARAMETER                            :: type_lake           = 2
    INTEGER, PARAMETER                            :: type_sheet          = 3
    INTEGER, PARAMETER                            :: type_shelf          = 4
    INTEGER, PARAMETER                            :: type_coast          = 5
    INTEGER, PARAMETER                            :: type_margin         = 6
    INTEGER, PARAMETER                            :: type_groundingline  = 7
    INTEGER, PARAMETER                            :: type_calvingfront   = 8
    
    IF      (region_name == 'NAM') THEN
      PD%netcdf%filename   = C%filename_PD_NAM
    ELSE IF (region_name == 'EAS') THEN
      PD%netcdf%filename   = C%filename_PD_EAS
    ELSE IF (region_name == 'GRL') THEN
      PD%netcdf%filename   = C%filename_PD_GRL
    ELSE IF (region_name == 'ANT') THEN
      PD%netcdf%filename   = C%filename_PD_ANT
    END IF  
    
    ! Determine size of Cartesian grid, so that memory can be allocated accordingly.
    CALL allocate_shared_int_0D(            PD%nx,      PD%wnx     )
    CALL allocate_shared_int_0D(            PD%ny,      PD%wny     )
    
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
            C%choice_benchmark_experiment == 'Bueler') THEN
          PD%nx = 51
          PD%ny = 51
        ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
          PD%nx = 900
          PD%ny = 900
        ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
          PD%nx = 500
          PD%ny = 500
        ELSE
          IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_data_fields!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
        
      ELSE ! IF (C%do_benchmark_experiment) THEN
      
        ! Read data from input file
        CALL inquire_PD_data_file(PD)
        
      END IF ! IF (C%do_benchmark_experiment) THEN
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate memory - PD    
    CALL allocate_shared_dp_1D( PD%nx,            PD%x,            PD%wx           )
    CALL allocate_shared_dp_1D( PD%ny,            PD%y,            PD%wy           )
    CALL allocate_shared_dp_2D( PD%nx, PD%ny,     PD%Hi_cart,      PD%wHi_cart     )
    CALL allocate_shared_dp_2D( PD%nx, PD%ny,     PD%Hb_cart,      PD%wHb_cart     )
    CALL allocate_shared_dp_2D( PD%nx, PD%ny,     PD%Hs_cart,      PD%wHs_cart     )
    CALL allocate_shared_int_2D(PD%nx, PD%ny,     PD%mask_cart,    PD%wmask_cart   )
    
    ! Read PD data
    IF (par%master) THEN
    
      IF (C%do_benchmark_experiment) THEN
        IF (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
            C%choice_benchmark_experiment == 'EISMINT_6' .OR. &
            C%choice_benchmark_experiment == 'Halfar' .OR. &
            C%choice_benchmark_experiment == 'Bueler') THEN
      
          ! Simple square grid with no data
          DO i = 1, PD%nx
            PD%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (PD%nx-1)
          END DO
          DO j = 1, PD%ny
            PD%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (PD%ny-1)
          END DO
          
          PD%Hi_cart = 0._dp
          PD%Hb_cart = 0._dp
          PD%Hs_cart = 0._dp
          
        ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
          ! Bedrock slope as described by Pattyn et al., 2012
        
          DO i = 1, PD%nx
          DO j = 1, PD%ny
            PD%x(i) = MISMIP_xmin + (MISMIP_xmax - MISMIP_xmin) * (i-1) / (PD%nx-1)
            PD%y(j) = MISMIP_ymin + (MISMIP_ymax - MISMIP_ymin) * (j-1) / (PD%ny-1)
            R = NORM2([ PD%x(i), PD%y(j)])
            PD%Hb_cart(i,j) = 720._dp - 778.5_dp * (R / 750000._dp)
            PD%Hs_cart(i,j) = MAX( 0._dp, PD%Hb_cart(i,j))
          END DO
          END DO
          
          PD%Hi_cart = 0._dp
          
        ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
          ! A small circular island with a little bay on the southeast side
        
          DO i = 1, PD%nx
          DO j = 1, PD%ny
            PD%x(i) = -1000000._dp + 2000000._dp * (i-1) / (PD%nx-1)
            PD%y(j) = -1000000._dp + 2000000._dp * (j-1) / (PD%ny-1)
            R = NORM2([ PD%x(i), PD%y(j)]) / 1000._dp
            PD%Hb_cart(i,j) = -450._dp - (400._dp * ATAN(( R - 600._dp) / 180._dp))
            R = NORM2([ (PD%x(i) - 200000._dp), (PD%y(j) + 100000._dp)]) / 1000._dp
            PD%Hb_cart(i,j) = PD%Hb_cart(i,j) + (350._dp * ATAN(( R - 240._dp) / 160._dp))
            PD%Hs_cart(i,j) = MAX( 0._dp, PD%Hb_cart(i,j))
          END DO
          END DO
          
          PD%Hi_cart = 0._dp
          
        ELSE
        
          WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_PD_data_fields!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
            
        END IF
        
      ELSE ! IF (C%do_benchmark_experiment) THEN
      
        ! Read data from input file
        CALL read_PD_data_file(PD)
        
      END IF ! IF (C%do_benchmark_experiment) THEN
      
    END IF
    CALL sync
      
    ! Determine vertex domains
    CALL partition_list( PD%nx, par%i, par%n, PD%i1, PD%i2)
    
    ! Determine the masks
    PD%mask_cart( PD%i1:PD%i2,:) = 0
  
    DO i = PD%i1, PD%i2
    DO j = 1, PD%ny
      IF (PD%Hi_cart(i,j) <= 0._dp) THEN
        IF (PD%Hb_cart(i,j) < 0._dp) THEN
          PD%mask_cart(i,j) = type_ocean
        ELSE
          PD%mask_cart(i,j) = type_land
        END IF
      ELSE
        ! Either sheet or shelf, depending on flotation criterion
        IF (PD%Hi_cart(i,j) > (-1._dp * PD%Hb_cart(i,j)) * seawater_density / ice_density) THEN
          PD%mask_cart(i,j) = type_sheet
        ELSE
          PD%mask_cart(i,j) = type_shelf
        END IF
      END IF
    END DO
    END DO
    CALL sync
  
    DO i = MAX(2,PD%i1), MIN(PD%nx-1,PD%i2)
    DO j = 2, PD%ny-1
      IF (PD%mask_cart(i,j) == type_sheet) THEN
      IF ((PD%mask_cart(i-1,j) == type_shelf .OR. PD%mask_cart(i-1,j) == type_ocean) .OR. &
          (PD%mask_cart(i+1,j) == type_shelf .OR. PD%mask_cart(i+1,j) == type_ocean) .OR. &
          (PD%mask_cart(i,j-1) == type_shelf .OR. PD%mask_cart(i,j-1) == type_ocean) .OR. &
          (PD%mask_cart(i,j+1) == type_shelf .OR. PD%mask_cart(i,j+1) == type_ocean)) THEN
        PD%mask_cart(i,j) = type_groundingline
      END IF
      END IF
  
      IF (PD%mask_cart(i,j) == type_shelf) THEN
      IF ((PD%mask_cart(i-1,j) == type_ocean) .OR. &
          (PD%mask_cart(i+1,j) == type_ocean) .OR. &
          (PD%mask_cart(i,j-1) == type_ocean) .OR. &
          (PD%mask_cart(i,j+1) == type_ocean)) THEN
        PD%mask_cart(i,j) = type_calvingfront
      END IF
      END IF
    END DO 
    END DO ! DO i = MAX(2,PD%i1), MIN(PD%nx-1,PD%i2)
    CALL sync
    
  END SUBROUTINE initialise_PD_data_fields
  SUBROUTINE initialise_init_data_fields_cart( init, region_name)
    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
     
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
    INTEGER                                       :: i, j, cerr, ierr
    REAL(dp)                                      :: R
    REAL(dp)                                      :: dx, dy, d2dx2, d2dy2, d2dxdy
    
    REAL(dp), PARAMETER                           :: EISMINT_xmin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_xmax =  750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymin = -750000._dp
    REAL(dp), PARAMETER                           :: EISMINT_ymax =  750000._dp
    
    REAL(dp), PARAMETER                           :: MISMIP_xmin = -1800000._dp
    REAL(dp), PARAMETER                           :: MISMIP_xmax =  1800000._dp
    REAL(dp), PARAMETER                           :: MISMIP_ymin = -1800000._dp
    REAL(dp), PARAMETER                           :: MISMIP_ymax =  1800000._dp

    INTEGER, PARAMETER                            :: type_land           = 0
    INTEGER, PARAMETER                            :: type_ocean          = 1
    INTEGER, PARAMETER                            :: type_lake           = 2
    INTEGER, PARAMETER                            :: type_sheet          = 3
    INTEGER, PARAMETER                            :: type_shelf          = 4
    INTEGER, PARAMETER                            :: type_coast          = 5
    INTEGER, PARAMETER                            :: type_margin         = 6
    INTEGER, PARAMETER                            :: type_groundingline  = 7
    INTEGER, PARAMETER                            :: type_calvingfront   = 8
    
    IF      (region_name == 'NAM') THEN
      init%netcdf%filename   = C%filename_init_NAM
    ELSE IF (region_name == 'EAS') THEN
      init%netcdf%filename   = C%filename_init_EAS
    ELSE IF (region_name == 'GRL') THEN
      init%netcdf%filename   = C%filename_init_GRL
    ELSE IF (region_name == 'ANT') THEN
      init%netcdf%filename   = C%filename_init_ANT
    END IF  
    
    ! Determine size of Cartesian grid, so that memory can be allocated accordingly.
    CALL allocate_shared_int_0D( init%nx, init%wnx)
    CALL allocate_shared_int_0D( init%ny, init%wny)
    
    ! For the benchmark experiments, use dummy input data.
    ! For realistic experiments, read the provided input file.
    IF (par%master) THEN
    
      IF (C%do_benchmark_experiment) THEN
        IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_6') THEN
          init%nx = 51
          init%ny = 51
        ELSEIF (C%choice_benchmark_experiment == 'Halfar' .OR. &
                C%choice_benchmark_experiment == 'Bueler') THEN
          init%nx = 2501
          init%ny = 2501
        ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
          init%nx = 900
          init%ny = 900
        ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
          init%nx = 500
          init%ny = 500
        ELSE
          IF (par%master) WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_init_data_fields_cart!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        END IF
      ELSE ! IF (C%do_benchmark_experiment) THEN
      
        ! Read data from input file
        CALL inquire_init_data_file_cart( init)
        
      END IF ! IF (C%do_benchmark_experiment) THEN
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate memory - init    
    CALL allocate_shared_dp_1D(  init%nx,          init%x,                      init%wx                     )
    CALL allocate_shared_dp_1D(  init%ny,          init%y,                      init%wy                     )
    CALL allocate_shared_dp_2D(  init%nx, init%ny, init%Hi_cart,                init%wHi_cart               )
    CALL allocate_shared_dp_2D(  init%nx, init%ny, init%Hb_cart,                init%wHb_cart               )
    CALL allocate_shared_dp_2D(  init%nx, init%ny, init%Hs_cart,                init%wHs_cart               )
    CALL allocate_shared_dp_2D(  init%nx, init%ny, init%surface_curvature_cart, init%wsurface_curvature_cart)
    CALL allocate_shared_int_2D( init%nx, init%ny, init%mask_cart,              init%wmask_cart             )
    CALL allocate_shared_int_2D( init%nx, init%ny, init%mask_ice_cart,          init%wmask_ice_cart         )
    CALL allocate_shared_int_2D( init%nx, init%ny, init%mask_gl_cart,           init%wmask_gl_cart          )
    CALL allocate_shared_int_2D( init%nx, init%ny, init%mask_cf_cart,           init%wmask_cf_cart          )
    CALL allocate_shared_int_2D( init%nx, init%ny, init%mask_coast_cart,        init%wmask_coast_cart       )
    
    ! Read init data
    IF (par%master) THEN
    
      IF (C%do_benchmark_experiment) THEN
        IF     (C%choice_benchmark_experiment == 'EISMINT_1' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_2' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_3' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_4' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_5' .OR. &
                C%choice_benchmark_experiment == 'EISMINT_6') THEN
                
          ! Simple square grid with no data
          DO i = 1, init%nx
            init%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (init%nx-1)
          END DO
          DO j = 1, init%ny
            init%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (init%ny-1)
          END DO
          
          init%Hi_cart = 0._dp
          init%Hb_cart = 0._dp
          init%Hs_cart = 0._dp
          
        ELSEIF (C%choice_benchmark_experiment == 'Halfar') THEN
        
          ! Start with the Halfar solution for the given parameters
          DO i = 1, init%nx
          DO j = 1, init%ny
            init%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (init%nx-1)
            init%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (init%ny-1)
            init%Hi_cart(i,j) = Halfar_solution( C%halfar_solution_H0, C%halfar_solution_R0, init%x(i), init%y(j), 0._dp)
          END DO
          END DO
          init%Hs_cart = init%Hi_cart
        
        ELSEIF (C%choice_benchmark_experiment == 'Bueler') THEN
        
          ! Start with the Bueler solution for the given parameters
          DO i = 1, init%nx
          DO j = 1, init%ny
            init%x(i) = EISMINT_xmin + (EISMINT_xmax - EISMINT_xmin) * (i-1) / (init%nx-1)
            init%y(j) = EISMINT_ymin + (EISMINT_ymax - EISMINT_ymin) * (j-1) / (init%ny-1)
            init%Hi_cart(i,j) = Bueler_solution( C%halfar_solution_H0, C%halfar_solution_R0, C%bueler_solution_lambda, init%x(i), init%y(j), C%start_time_of_run)
          END DO
          END DO
          init%Hs_cart = init%Hi_cart
          
        ELSEIF (C%choice_benchmark_experiment == 'MISMIP_mod') THEN
          ! Bedrock slope as described by Pattyn et al., 2012
        
          DO i = 1, init%nx
          DO j = 1, init%ny
            init%x(i) = MISMIP_xmin + (MISMIP_xmax - MISMIP_xmin) * (i-1) / (init%nx-1)
            init%y(j) = MISMIP_ymin + (MISMIP_ymax - MISMIP_ymin) * (j-1) / (init%ny-1)
            R = NORM2([ init%x(i), init%y(j)])
            init%Hb_cart(i,j) = 720._dp - 778.5_dp * (R / 750000._dp)
            init%Hs_cart(i,j) = MAX( 0._dp, init%Hb_cart(i,j))
          END DO
          END DO
          
          init%Hi_cart = 100._dp
          
        ELSEIF (C%choice_benchmark_experiment == 'mesh_generation_test') THEN
          ! A small circular island with a little bay on the southeast side
        
          DO i = 1, init%nx
          DO j = 1, init%ny
            init%x(i) = -1000000._dp + 2000000._dp * (i-1) / (init%nx-1)
            init%y(j) = -1000000._dp + 2000000._dp * (j-1) / (init%ny-1)
            R = NORM2([ init%x(i), init%y(j)]) / 1000._dp
            init%Hb_cart(i,j) = -450._dp - (400._dp * ATAN(( R - 600._dp) / 180._dp))
            R = NORM2([ (init%x(i) - 200000._dp), (init%y(j) + 100000._dp)]) / 1000._dp
            init%Hb_cart(i,j) = init%Hb_cart(i,j) + (350._dp * ATAN(( R - 240._dp) / 160._dp))
            init%Hs_cart(i,j) = MAX( 0._dp, init%Hb_cart(i,j))
          END DO
          END DO
        
        ELSE
        
          WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_init_data_fields_cart!'
          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
        
        END IF
      ELSE ! IF (C%do_benchmark_experiment) THEN
      
        ! Read data from input file
        CALL read_init_data_file_cart(init) 
        
      END IF ! IF (C%do_benchmark_experiment) THEN
      
    END IF ! IF (par%master) THEN
    CALL sync
      
    ! Determine vertex domains
    CALL partition_list( init%nx, par%i, par%n, init%i1, init%i2)
    
    ! Surface curvature
    init%surface_curvature_cart(init%i1:init%i2,:) = 0._dp
    dx = ABS(init%x(2) - init%x(1))
    dy = ABS(init%y(2) - init%y(1))
    DO i = MAX(2,init%i1), MIN(init%nx-1,init%i2)
    DO j = 2, init%ny-1
      d2dx2  = (init%Hs_cart(i+1,j  ) + init%Hs_cart(i-1,j  ) - 2._dp*init%Hs_cart(i,j)) / (dx**2)
      d2dy2  = (init%Hs_cart(i  ,j+1) + init%Hs_cart(i  ,j-1) - 2._dp*init%Hs_cart(i,j)) / (dy**2)
      d2dxdy = (init%Hs_cart(i+1,j+1) + init%Hs_cart(i-1,j-1) - init%Hs_cart(i+1,j-1) - init%Hs_cart(i-1,j+1)) / (4*dx*dy)
      init%surface_curvature_cart(i,j) = MAX(-1E-6, MIN(1E-6, SQRT(d2dx2**2 + d2dy2**2 + d2dxdy**2)))
    END DO
    END DO
    
    ! Determine the masks
    init%mask_cart(       init%i1:init%i2,:) = 0
    init%mask_ice_cart(   init%i1:init%i2,:) = 0
    init%mask_gl_cart(    init%i1:init%i2,:) = 0
    init%mask_cf_cart(    init%i1:init%i2,:) = 0
    init%mask_coast_cart( init%i1:init%i2,:) = 0
  
    DO i = init%i1, init%i2
    DO j = 1, init%ny
      IF (init%Hi_cart(i,j) == 0._dp) THEN
        IF (init%Hb_cart(i,j) < 0._dp) THEN
          init%mask_cart(i,j) = type_ocean
        ELSE
          init%mask_cart(i,j) = type_land
        END IF
      ELSE
        ! Either sheet or shelf, depending on flotation criterion
        init%mask_ice_cart(i,j) = 1
        IF (init%Hi_cart(i,j) > (-1._dp * init%Hb_cart(i,j)) * seawater_density / ice_density) THEN
          init%mask_cart(i,j) = type_sheet
        ELSE
          init%mask_cart(i,j) = type_shelf
        END IF
      END IF
    END DO
    END DO
    CALL sync
  
    DO i = MAX(2,init%i1), MIN(init%nx-1,init%i2)
    DO j = 2, init%ny-1
      IF (init%mask_cart(i,j) == type_sheet) THEN
      IF ((init%mask_cart(i-1,j) == type_shelf .OR. init%mask_cart(i-1,j) == type_ocean) .OR. &
          (init%mask_cart(i+1,j) == type_shelf .OR. init%mask_cart(i+1,j) == type_ocean) .OR. &
          (init%mask_cart(i,j-1) == type_shelf .OR. init%mask_cart(i,j-1) == type_ocean) .OR. &
          (init%mask_cart(i,j+1) == type_shelf .OR. init%mask_cart(i,j+1) == type_ocean)) THEN
        init%mask_gl_cart(i,j) = 1
        init%mask_cart(i,j) = type_groundingline
      END IF
      END IF
  
      IF (init%mask_cart(i,j) == type_shelf) THEN
      IF ((init%mask_cart(i-1,j) == type_ocean) .OR. &
          (init%mask_cart(i+1,j) == type_ocean) .OR. &
          (init%mask_cart(i,j-1) == type_ocean) .OR. &
          (init%mask_cart(i,j+1) == type_ocean)) THEN
        init%mask_cf_cart(i,j) = 1
        init%mask_cart(i,j) = type_calvingfront
      END IF
      END IF
      
      IF (init%mask_cart(i,j) == type_land) THEN
      IF ((init%mask_cart(i-1,j) == type_ocean) .OR. &
          (init%mask_cart(i+1,j) == type_ocean) .OR. &
          (init%mask_cart(i,j-1) == type_ocean) .OR. &
          (init%mask_cart(i,j+1) == type_ocean)) THEN
        init%mask_coast_cart(i,j) = 1
      END IF
      END IF
    END DO
    END DO
    CALL sync      
    
  END SUBROUTINE initialise_init_data_fields_cart
!  SUBROUTINE initialise_init_data_fields_mesh( init, region_name, mesh)
!    ! Allocate memory for the reference data fields, read them from the specified NetCDF file (latter only done by master process).
!     
!    IMPLICIT NONE
!      
!    ! Input variables:
!    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
!    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
!    TYPE(type_mesh),                INTENT(INOUT) :: mesh
!    
!    INTEGER                                       :: nt
!    TYPE(type_netcdf_output)                      :: netcdf_placeholder
!    
!    IF      (region_name == 'NAM') THEN
!      init%netcdf%filename   = C%filename_init_NAM
!    ELSE IF (region_name == 'EAS') THEN
!      init%netcdf%filename   = C%filename_init_EAS
!    ELSE IF (region_name == 'GRL') THEN
!      init%netcdf%filename   = C%filename_init_GRL
!    ELSE IF (region_name == 'ANT') THEN
!      init%netcdf%filename   = C%filename_init_ANT
!    END IF
!    
!    ! Allocate memory for mesh metadata
!    CALL allocate_shared_dp_0D(  mesh%xmin,             mesh%wxmin           )
!    CALL allocate_shared_dp_0D(  mesh%xmax,             mesh%wxmax           )
!    CALL allocate_shared_dp_0D(  mesh%ymin,             mesh%wymin           )
!    CALL allocate_shared_dp_0D(  mesh%ymax,             mesh%wymax           )
!    CALL allocate_shared_int_0D( mesh%nV_mem,           mesh%wnV_mem         )
!    CALL allocate_shared_int_0D( mesh%nTri_mem,         mesh%wnTri_mem       )
!    CALL allocate_shared_int_0D( mesh%nC_mem,           mesh%wnC_mem         )
!    CALL allocate_shared_int_0D( mesh%nV,               mesh%wnV             )
!    CALL allocate_shared_int_0D( mesh%nTri,             mesh%wnTri           )
!    CALL allocate_shared_dp_0D(  mesh%alpha_min,        mesh%walpha_min      )
!    CALL allocate_shared_dp_0D(  mesh%dz_max_ice,       mesh%wdz_max_ice     )
!    CALL allocate_shared_dp_0D(  mesh%resolution_min,   mesh%wresolution_min )
!    CALL allocate_shared_dp_0D(  mesh%resolution_max,   mesh%wresolution_max )
!    
!    ! Check if all the required data fields are there
!    IF (par%master) THEN
!      CALL inquire_init_data_file_mesh(init, mesh, netcdf_placeholder, nt)      
!    END IF
!    CALL sync
!                
!    ! Determine vertex and triangle domains
!    mesh%v1 = MAX(1,         FLOOR(REAL(mesh%nV   *  par%i      / par%n)) + 1)
!    mesh%v2 = MIN(mesh%nV,   FLOOR(REAL(mesh%nV   * (par%i + 1) / par%n)))
!    mesh%t1 = MAX(1,         FLOOR(REAL(mesh%nTri *  par%i      / par%n)) + 1)
!    mesh%t2 = MIN(mesh%nTri, FLOOR(REAL(mesh%nTri * (par%i + 1) / par%n)))
!    CALL sync
!    
!    ! Allocate memory for mesh data
!    CALL allocate_shared_dp_2D(  mesh%nV,      2,           mesh%V,                mesh%wV               )
!    CALL allocate_shared_int_1D( mesh%nV,                   mesh%nC,               mesh%wnC              )
!    CALL allocate_shared_int_2D( mesh%nV,      mesh%nC_mem, mesh%C,                mesh%wC               )
!    CALL allocate_shared_int_1D( mesh%nV,                   mesh%niTri,            mesh%wniTri           )
!    CALL allocate_shared_int_2D( mesh%nV,      mesh%nC_mem, mesh%iTri,             mesh%wiTri            )
!    CALL allocate_shared_int_1D( mesh%nV,                   mesh%edge_index,       mesh%wedge_index      )
!    CALL allocate_shared_int_1D( mesh%nV,                   mesh%mesh_old_ti_in,   mesh%wmesh_old_ti_in  )
!    
!    CALL allocate_shared_int_2D( mesh%nTri,    3,           mesh%Tri,              mesh%wTri             )
!    CALL allocate_shared_int_2D( mesh%nTri,    3,           mesh%TriC,             mesh%wTriC            )
!    CALL allocate_shared_dp_2D(  mesh%nTri,    2,           mesh%Tricc,            mesh%wTricc           )
!    CALL allocate_shared_int_1D( mesh%nTri,                 mesh%Tri_edge_index,   mesh%wTri_edge_index  )
!    
!    CALL allocate_shared_int_2D( mesh%nTri,    2,           mesh%Triflip,          mesh%wTriflip         )
!    CALL allocate_shared_int_1D( mesh%nTri,                 mesh%RefMap,           mesh%wRefMap          )
!    CALL allocate_shared_int_1D( mesh%nTri,                 mesh%RefStack,         mesh%wRefStack        )
!    
!    ! Allocate memory for initial data
!    CALL allocate_shared_dp_1D(  mesh%nV,                   init%Hi,               init%wHi              )
!    CALL allocate_shared_dp_1D(  mesh%nV,                   init%Hb,               init%wHb              )
!    CALL allocate_shared_dp_1D(  mesh%nV,                   init%Hs,               init%wHs              )
!    CALL allocate_shared_int_1D( mesh%nV,                   init%mask,             init%wmask            )
!    CALL allocate_shared_dp_1D(  mesh%nV,                   init%U_SSA,            init%wU_SSA           )
!    CALL allocate_shared_dp_1D(  mesh%nV,                   init%V_SSA,            init%wV_SSA           )
!    CALL allocate_shared_dp_2D(  mesh%nV,      C%nZ,        init%Ti,               init%wTi              )    
!    CALL allocate_shared_dp_2D(  mesh%nV,      12,          init%Albedo,           init%wAlbedo          )
!    CALL allocate_shared_dp_2D(  mesh%nV,      12,          init%FirnDepth,        init%wFirnDepth       )
!    CALL allocate_shared_dp_1D(  mesh%nV,                   init%MeltPreviousYear, init%wMeltPreviousYear)
!    
!    ! Read mesh data and initial data
!    IF (par%master) THEN
!      CALL read_init_data_file_mesh(init, mesh, netcdf_placeholder, nt)
!    END IF
!    CALL sync
!    
!  END SUBROUTINE initialise_init_data_fields_mesh
  
  SUBROUTINE deallocate_init_data( init, filetype_init)
    
    IMPLICIT NONE
      
    ! Input variables:
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    CHARACTER(LEN=4),               INTENT(IN)    :: filetype_init
    
    IF (filetype_init == 'cart') THEN
    
      CALL deallocate_shared( init%wnx)
      CALL deallocate_shared( init%wny)
      CALL deallocate_shared( init%wx)
      CALL deallocate_shared( init%wy)
      CALL deallocate_shared( init%wHi_cart)
      CALL deallocate_shared( init%wHb_cart)
      CALL deallocate_shared( init%wHs_cart)
      CALL deallocate_shared( init%wsurface_curvature_cart)
      CALL deallocate_shared( init%wmask_cart)
      CALL deallocate_shared( init%wmask_ice_cart)
      CALL deallocate_shared( init%wmask_gl_cart)
      CALL deallocate_shared( init%wmask_cf_cart)
      CALL deallocate_shared( init%wmask_coast_cart)
      CALL deallocate_shared( init%wHi)
      CALL deallocate_shared( init%wHb)
      CALL deallocate_shared( init%wHs)
      CALL deallocate_shared( init%wmask)
      
    ELSE
    
      CALL deallocate_shared( init%wHi)
      CALL deallocate_shared( init%wHb)
      CALL deallocate_shared( init%wHs)
      CALL deallocate_shared( init%wmask)
      CALL deallocate_shared( init%wU_SSA)
      CALL deallocate_shared( init%wV_SSA)
      CALL deallocate_shared( init%wTi)
      CALL deallocate_shared( init%wAlbedo)
      CALL deallocate_shared( init%wFirnDepth)
      CALL deallocate_shared( init%wMeltPreviousYear)    
    
    END IF
    
  END SUBROUTINE deallocate_init_data
  
  SUBROUTINE map_PD_data_to_mesh(   mesh, PD)
    ! Calculate the mapping arrays, map data to the mesh. Allocate and deallocate memory for the mapping arrays within this subroutine.
    
    IMPLICIT NONE
  
    ! Input and output variables
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_PD_data_fields),      INTENT(INOUT) :: PD
    
    ! Local variables
    INTEGER,  DIMENSION(:,:  ), POINTER           :: map_mean_i  ! On cart grid: vertex indices
    REAL(dp), DIMENSION(:,:  ), POINTER           :: map_mean_w  ! On cart grid: weights
    INTEGER,  DIMENSION(:,:  ), POINTER           :: map_bilin_i ! On mesh     : cart indices (4x i and j))
    REAL(dp), DIMENSION(:,:  ), POINTER           :: map_bilin_w ! On mesh     : weights      (4x)
    INTEGER                                       :: wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w  
    
    ! Allocate shared memory for the mapped data
    CALL allocate_shared_dp_1D( mesh%nV,     PD%Hi,      PD%wHi     )
    CALL allocate_shared_dp_1D( mesh%nV,     PD%Hb,      PD%wHb     )
    CALL allocate_shared_dp_1D( mesh%nV,     PD%Hs,      PD%wHs     )
    CALL allocate_shared_int_1D(mesh%nV,     PD%mask,    PD%wmask   )
    
    ! Calculate mapping arrays
    CALL get_cart_to_mesh_map( mesh, PD%x, PD%y, PD%nx, PD%ny, map_mean_i,  map_mean_w,  map_bilin_i,  map_bilin_w, &
                                                         wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w)
    
    ! Map PD data to the mesh
    CALL map_cart_to_mesh( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, PD%ny, PD%i1, PD%i2,   PD%Hi_cart,    PD%Hi   )
    CALL map_cart_to_mesh( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, PD%ny, PD%i1, PD%i2,   PD%Hb_cart,    PD%Hb   )
    PD%Hs(mesh%v1:mesh%v2) = MAX(0._dp, PD%Hb(mesh%v1:mesh%v2) + PD%Hi(mesh%v1:mesh%v2))
    
    ! Deallocate shared memory for the mapping arrays
    CALL deallocate_shared( wmap_mean_i)
    CALL deallocate_shared( wmap_mean_w)
    CALL deallocate_shared( wmap_bilin_i)
    CALL deallocate_shared( wmap_bilin_w)    
  
  END SUBROUTINE map_PD_data_to_mesh
  SUBROUTINE map_init_data_to_mesh( mesh, init)
    ! Calculate the mapping arrays, map data to the mesh. Allocate and deallocate memory for the mapping arrays within this subroutine.
    
    IMPLICIT NONE
  
    ! Input and output variables
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_init_data_fields),    INTENT(INOUT) :: init
    
    ! Local variables
    INTEGER,  DIMENSION(:,:  ), POINTER           :: map_mean_i  ! On cart grid: vertex indices
    REAL(dp), DIMENSION(:,:  ), POINTER           :: map_mean_w  ! On cart grid: weights
    INTEGER,  DIMENSION(:,:  ), POINTER           :: map_bilin_i ! On mesh     : cart indices (4x i and j))
    REAL(dp), DIMENSION(:,:  ), POINTER           :: map_bilin_w ! On mesh     : weights      (4x)
    INTEGER                                       :: wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w
    
    INTEGER :: vi
    
    ! Allocate shared memory for the mapped data
    CALL allocate_shared_dp_1D( mesh%nV,     init%Hi,      init%wHi     )
    CALL allocate_shared_dp_1D( mesh%nV,     init%Hb,      init%wHb     )
    CALL allocate_shared_dp_1D( mesh%nV,     init%Hs,      init%wHs     )
    CALL allocate_shared_int_1D(mesh%nV,     init%mask,    init%wmask   )
    
    ! Calculate mapping arrays
    CALL get_cart_to_mesh_map( mesh, init%x, init%y, init%nx, init%ny,  map_mean_i,  map_mean_w,  map_bilin_i,  map_bilin_w, &
                                                                  wmap_mean_i, wmap_mean_w, wmap_bilin_i, wmap_bilin_w)
    CALL sync
    
    ! Map init data to the mesh
    CALL map_cart_to_mesh( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, init%ny, init%i1, init%i2,   init%Hi_cart,    init%Hi   )
    CALL map_cart_to_mesh( mesh, map_mean_i, map_mean_w, map_bilin_i, map_bilin_w, init%ny, init%i1, init%i2,   init%Hb_cart,    init%Hb   )
    init%Hs(mesh%v1:mesh%v2) = MAX(0._dp, init%Hb(mesh%v1:mesh%v2) + init%Hi(mesh%v1:mesh%v2))
    
    ! Deallocate shared memory for the mapping arrays
    CALL deallocate_shared( wmap_mean_i)
    CALL deallocate_shared( wmap_mean_w)
    CALL deallocate_shared( wmap_bilin_i)
    CALL deallocate_shared( wmap_bilin_w)
    
    ! For the Bueler benchmark experiment, recalculate initial ice thickness on the mesh,
    ! to make sure we have the correct starting conditions.
    IF (C%do_benchmark_experiment .AND. C%choice_benchmark_experiment == 'Bueler') THEN
    
      DO vi = mesh%v1, mesh%v2
        init%Hi(vi) = Bueler_solution( C%halfar_solution_H0, C%halfar_solution_R0, C%bueler_solution_lambda, mesh%V(vi,1), mesh%V(vi,2), C%start_time_of_run)
      END DO      
      init%Hs(mesh%v1:mesh%v2) = init%Hi(mesh%v1:mesh%v2)
      
    END IF   
  
  END SUBROUTINE map_init_data_to_mesh
  
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
