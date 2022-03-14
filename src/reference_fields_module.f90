MODULE reference_fields_module
  
  ! Contains the routines for setting up the three "reference geometries":
  ! - refgeo_PD:     present-day, used to calculate sea-level contribution, isotope change, and more
  ! - refgeo_init:   initial, used to initialise the simulation
  ! - refgeo_GIA_eq: GIA equilibrium, used for the GIA model

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE petscksp
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
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  
  ! Import specific functionality
  USE data_types_module,               ONLY: type_model_region, type_grid, type_reference_geometry, type_mesh
  USE netcdf_module,                   ONLY: inquire_reference_geometry_file, read_reference_geometry_file
  USE mesh_mapping_module,             ONLY: calc_remapping_operator_grid2mesh, map_grid2mesh_2D, deallocate_remapping_operators_grid2mesh
  USE utilities_module,                ONLY: is_floating, surface_elevation, remove_Lake_Vostok

  IMPLICIT NONE
  
CONTAINS

  ! Initialise all three reference geometries
  SUBROUTINE initialise_reference_geometries( refgeo_init, refgeo_PD, refgeo_GIAeq, region_name)
    ! Initialise all three reference geometries
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo_init
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo_PD
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo_GIAeq
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometries'
    CHARACTER(LEN=256)                            :: choice_refgeo_init, choice_refgeo_PD, choice_refgeo_GIAeq
    CHARACTER(LEN=256)                            :: filename_refgeo_init, filename_refgeo_PD, filename_refgeo_GIAeq
    REAL(dp)                                      :: time_to_restart_from
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Determine parameters for this region
    IF     (region_name == 'NAM') THEN
      choice_refgeo_init    = C%choice_refgeo_init_NAM
      choice_refgeo_PD      = C%choice_refgeo_PD_NAM
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_NAM
      filename_refgeo_init  = C%filename_refgeo_init_NAM
      filename_refgeo_PD    = C%filename_refgeo_PD_NAM
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_NAM
      time_to_restart_from  = C%time_to_restart_from_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_refgeo_init    = C%choice_refgeo_init_EAS
      choice_refgeo_PD      = C%choice_refgeo_PD_EAS
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_EAS
      filename_refgeo_init  = C%filename_refgeo_init_EAS
      filename_refgeo_PD    = C%filename_refgeo_PD_EAS
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_EAS
      time_to_restart_from  = C%time_to_restart_from_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_refgeo_init    = C%choice_refgeo_init_GRL
      choice_refgeo_PD      = C%choice_refgeo_PD_GRL
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_GRL
      filename_refgeo_init  = C%filename_refgeo_init_GRL
      filename_refgeo_PD    = C%filename_refgeo_PD_GRL
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_GRL
      time_to_restart_from  = C%time_to_restart_from_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_refgeo_init    = C%choice_refgeo_init_ANT
      choice_refgeo_PD      = C%choice_refgeo_PD_ANT
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_ANT
      filename_refgeo_init  = C%filename_refgeo_init_ANT
      filename_refgeo_PD    = C%filename_refgeo_PD_ANT
      filename_refgeo_GIAeq = C%filename_refgeo_GIAeq_ANT
      time_to_restart_from  = C%time_to_restart_from_ANT
    END IF
    
    ! Initial ice-sheet geometry
    ! ==========================
    
    IF     (choice_refgeo_init == 'idealised') THEN
      IF (par%master) WRITE(0,*) '  Initialising initial         reference geometry from idealised case "', TRIM(C%choice_refgeo_init_idealised), '"...'
      CALL initialise_reference_geometry_idealised_grid( refgeo_init, C%choice_refgeo_init_idealised, region_name, C%dx_refgeo_init_idealised)
    ELSEIF (choice_refgeo_init == 'realistic') THEN
      IF (par%master) WRITE(0,*) '  Initialising initial         reference geometry from file ', TRIM( filename_refgeo_init), '...'
      CALL initialise_reference_geometry_from_file( refgeo_init, filename_refgeo_init, region_name)
    ELSE
      CALL crash('unknown choice_refgeo_init "' // TRIM( choice_refgeo_init) // '"!')
    END IF
    
    ! Present-day ice-sheet geometry
    ! ==============================
    
    IF     (choice_refgeo_PD == 'idealised') THEN
      IF (par%master) WRITE(0,*) '  Initialising present-day     reference geometry from idealised case "', TRIM(C%choice_refgeo_PD_idealised), '"...'
      CALL initialise_reference_geometry_idealised_grid( refgeo_PD, C%choice_refgeo_PD_idealised, region_name, C%dx_refgeo_PD_idealised)
    ELSEIF (choice_refgeo_PD == 'realistic') THEN
      IF (par%master) WRITE(0,*) '  Initialising present-day     reference geometry from file ', TRIM( filename_refgeo_PD), '...'
      CALL initialise_reference_geometry_from_file( refgeo_PD, filename_refgeo_PD, region_name)
    ELSE
      CALL crash('unknown choice_refgeo_PD "' // TRIM( choice_refgeo_PD) // '"!')
    END IF
    
    ! GIA equilibrium ice-sheet geometry
    ! ==================================
    
    IF     (choice_refgeo_GIAeq == 'idealised') THEN
      IF (par%master) WRITE(0,*) '  Initialising GIA equilibrium reference geometry from idealised case "', TRIM(C%choice_refgeo_GIAeq_idealised), '"...'
      CALL initialise_reference_geometry_idealised_grid( refgeo_GIAeq, C%choice_refgeo_GIAeq_idealised, region_name, C%dx_refgeo_GIAeq_idealised)
    ELSEIF (choice_refgeo_GIAeq == 'realistic') THEN
      IF (par%master) WRITE(0,*) '  Initialising GIA equilibrium reference geometry from file ', TRIM( filename_refgeo_GIAeq), '...'
      CALL initialise_reference_geometry_from_file( refgeo_GIAeq, filename_refgeo_GIAeq, region_name)
    ELSE
      CALL crash('unknown choice_refgeo_GIAeq "' // TRIM( choice_refgeo_GIAeq) // '"!')
    END IF
    
    ! Fill in secondary data for the reference geometry (used to force mesh creation)
    CALL calc_reference_geometry_secondary_data( refgeo_init%grid , refgeo_init )
    CALL calc_reference_geometry_secondary_data( refgeo_PD%grid   , refgeo_PD   )
    CALL calc_reference_geometry_secondary_data( refgeo_GIAeq%grid, refgeo_GIAeq)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 78)
    
  END SUBROUTINE initialise_reference_geometries
  
  ! Initialise a reference geometry with data from a (timeless) NetCDF file (e.g. BedMachine)
  SUBROUTINE initialise_reference_geometry_from_file( refgeo, filename_refgeo, region_name)
    ! Initialise a reference geometry with data from a NetCDF file
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    CHARACTER(LEN=256),             INTENT(IN)    :: filename_refgeo
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_from_file'
    INTEGER                                       :: i,j,n
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp
    
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Inquire if all the required fields are present in the specified NetCDF file,
    ! and determine the dimensions of the memory to be allocated.
    CALL allocate_shared_int_0D( refgeo%grid%nx, refgeo%grid%wnx)
    CALL allocate_shared_int_0D( refgeo%grid%ny, refgeo%grid%wny)
    IF (par%master) THEN
      refgeo%netcdf%filename = filename_refgeo
      CALL inquire_reference_geometry_file( refgeo)
    END IF
    CALL sync

    ! Assign range to each processor
    CALL partition_list( refgeo%grid%nx, par%i, par%n, refgeo%grid%i1, refgeo%grid%i2)
    CALL partition_list( refgeo%grid%ny, par%i, par%n, refgeo%grid%j1, refgeo%grid%j2)

    ! Allocate memory for raw data
    CALL allocate_shared_dp_1D( refgeo%grid%nx,                 refgeo%grid%x,  refgeo%grid%wx )
    CALL allocate_shared_dp_1D(                 refgeo%grid%ny, refgeo%grid%y,  refgeo%grid%wy )

    CALL allocate_shared_dp_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%Hi_grid, refgeo%wHi_grid)
    CALL allocate_shared_dp_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%Hb_grid, refgeo%wHb_grid)
    CALL allocate_shared_dp_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%Hs_grid, refgeo%wHs_grid)

    ! Read data from input file
    IF (par%master) CALL read_reference_geometry_file( refgeo)
    CALL sync

    ! Fill in secondary grid parameters
    CALL allocate_shared_dp_0D( refgeo%grid%dx,   refgeo%grid%wdx  )
    CALL allocate_shared_dp_0D( refgeo%grid%xmin, refgeo%grid%wxmin)
    CALL allocate_shared_dp_0D( refgeo%grid%xmax, refgeo%grid%wxmax)
    CALL allocate_shared_dp_0D( refgeo%grid%ymin, refgeo%grid%wymin)
    CALL allocate_shared_dp_0D( refgeo%grid%ymax, refgeo%grid%wymax)
    IF (par%master) THEN
      refgeo%grid%dx = refgeo%grid%x( 2) - refgeo%grid%x( 1)
      refgeo%grid%xmin = refgeo%grid%x( 1             )
      refgeo%grid%xmax = refgeo%grid%x( refgeo%grid%nx)
      refgeo%grid%ymin = refgeo%grid%y( 1             )
      refgeo%grid%ymax = refgeo%grid%y( refgeo%grid%ny)
    END IF
    CALL sync

    ! Tolerance; points lying within this distance of each other are treated as identical
    CALL allocate_shared_dp_0D( refgeo%grid%tol_dist, refgeo%grid%wtol_dist)
    IF (par%master) refgeo%grid%tol_dist = ((refgeo%grid%xmax - refgeo%grid%xmin) + (refgeo%grid%ymax - refgeo%grid%ymin)) * tol / 2._dp
    
    ! Set up grid-to-vector translation tables
    CALL allocate_shared_int_0D(                   refgeo%grid%n           , refgeo%grid%wn           )
    IF (par%master) refgeo%grid%n  = refgeo%grid%nx * refgeo%grid%ny
    CALL sync
    CALL allocate_shared_int_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%grid%ij2n        , refgeo%grid%wij2n        )
    CALL allocate_shared_int_2D( refgeo%grid%n , 2      , refgeo%grid%n2ij        , refgeo%grid%wn2ij        )
    IF (par%master) THEN
      n = 0
      DO i = 1, refgeo%grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, refgeo%grid%ny
            n = n+1
            refgeo%grid%ij2n( i,j) = n
            refgeo%grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = refgeo%grid%ny, 1, -1
            n = n+1
            refgeo%grid%ij2n( i,j) = n
            refgeo%grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO
    END IF
    CALL sync
    
    ! Safety
    CALL check_for_NaN_dp_2D( refgeo%Hi_grid, 'refgeo%Hi_grid')
    CALL check_for_NaN_dp_2D( refgeo%Hb_grid, 'refgeo%Hb_grid')
    CALL check_for_NaN_dp_2D( refgeo%Hs_grid, 'refgeo%Hs_grid')
    
    ! Remove Lake Vostok from Antarctica (because it's annoying)
    IF (region_name == 'ANT'.AND. C%remove_Lake_Vostok) THEN
      CALL remove_Lake_Vostok( refgeo%grid%x, refgeo%grid%y, refgeo%Hi_grid, refgeo%Hb_grid, refgeo%Hs_grid)
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 16)
    
  END SUBROUTINE initialise_reference_geometry_from_file
  
  ! Initialise a reference geometry according to an idealised world on the initial grid
  SUBROUTINE initialise_reference_geometry_idealised_grid( refgeo, choice_refgeo_idealised, region_name, dx)
    ! Initialise a reference geometry according to an idealised world
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    CHARACTER(LEN=256),             INTENT(IN)    :: choice_refgeo_idealised
    CHARACTER(LEN=3),               INTENT(IN)    :: region_name
    REAL(dp),                       INTENT(IN)    :: dx
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Set up a grid
    CALL setup_idealised_geometry_grid( refgeo%grid, region_name, dx)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%Hi_grid, refgeo%wHi_grid)
    CALL allocate_shared_dp_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%Hb_grid, refgeo%wHb_grid)
    CALL allocate_shared_dp_2D( refgeo%grid%nx, refgeo%grid%ny, refgeo%Hs_grid, refgeo%wHs_grid)
    
    IF     (choice_refgeo_idealised == 'flatearth') THEN
      ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments
      CALL initialise_reference_geometry_idealised_grid_flatearth(     refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'Halfar') THEN
      ! The Halfar dome solution at t = 0
      CALL initialise_reference_geometry_idealised_grid_Halfar(        refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'Bueler') THEN
      ! The Bueler dome solution at t = 0
      CALL initialise_reference_geometry_idealised_grid_Bueler(        refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'SSA_icestream') THEN
      ! The SSA_icestream infinite slab on a flat slope
      CALL initialise_reference_geometry_idealised_grid_SSA_icestream( refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'MISMIP_mod') THEN
      ! The MISMIP_mod cone-shaped island
      CALL initialise_reference_geometry_idealised_grid_MISMIP_mod(    refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_A') THEN
      ! The ISMIP-HOM A bumpy slope
      CALL initialise_reference_geometry_idealised_grid_ISMIP_HOM_A(   refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_B') THEN
      ! The ISMIP-HOM B bumpy slope
      CALL initialise_reference_geometry_idealised_grid_ISMIP_HOM_B(   refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_C' .OR. &
            choice_refgeo_idealised == 'ISMIP_HOM_D') THEN
      ! The ISMIP-HOM C/D bumpy slope
      CALL initialise_reference_geometry_idealised_grid_ISMIP_HOM_CD(  refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_E') THEN
      ! The ISMIP-HOM E Glacier d'Arolla geometry
      CALL initialise_reference_geometry_idealised_grid_ISMIP_HOM_E(   refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_F') THEN
      ! The ISMIP-HOM A bumpy slope
      CALL initialise_reference_geometry_idealised_grid_ISMIP_HOM_F(   refgeo%grid, refgeo)
    ELSEIF (choice_refgeo_idealised == 'MISMIP+') THEN
      ! The MISMIP+ fjord geometry
      CALL initialise_reference_geometry_idealised_grid_MISMIPplus(    refgeo%grid, refgeo)
    ELSE
      CALL crash('unknown choice_refgeo_idealised "' // TRIM( choice_refgeo_idealised) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 16)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid
  SUBROUTINE initialise_reference_geometry_idealised_grid_flatearth(     grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_flatearth'
    INTEGER                                       :: i,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hi_grid( i,j) = 0._dp
      refgeo%Hb_grid( i,j) = 0._dp
      refgeo%Hs_grid( i,j) = 0._dp
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_flatearth
  SUBROUTINE initialise_reference_geometry_idealised_grid_Halfar(        grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The Halfar dome solution at t = 0
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_Halfar'
    INTEGER                                       :: i,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hi_grid( i,j) = Halfar_solution( grid%x( i), grid%y( j), C%start_time_of_run)
      refgeo%Hb_grid( i,j) = 0._dp
      refgeo%Hs_grid( i,j) = refgeo%Hi_grid( i,j)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_Halfar
  SUBROUTINE initialise_reference_geometry_idealised_grid_Bueler(        grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The Bueler dome solution at t = 0
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_Bueler'
    INTEGER                                       :: i,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hi_grid( i,j) = Bueler_solution( grid%x( i), grid%y( j), C%start_time_of_run)
      refgeo%Hb_grid( i,j) = 0._dp
      refgeo%Hs_grid( i,j) = refgeo%Hi_grid( i,j)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_Bueler
  SUBROUTINE initialise_reference_geometry_idealised_grid_SSA_icestream( grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The SSA_icestream infinite slab on a flat slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_SSA_icestream'
    INTEGER                                       :: i,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hi_grid( i,j) = 2000._dp
      refgeo%Hb_grid( i,j) = -0.001_dp * grid%x( i)
      refgeo%Hs_grid( i,j) = refgeo%Hb_grid( i,j) + refgeo%Hi_grid( i,j)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_SSA_icestream
  SUBROUTINE initialise_reference_geometry_idealised_grid_MISMIP_mod(    grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The MISMIP_mod cone-shaped island
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_MISMIP_mod'
    INTEGER                                       :: i,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Create a nice circular ice shelf
      IF (SQRT(grid%x( i)**2 + grid%y( j)**2) < grid%xmax * 0.95_dp) THEN
        refgeo%Hi_grid( i,j) = 100._dp
      ELSE
        refgeo%Hi_grid( i,j) = 0._dp
      END IF
      
      refgeo%Hb_grid( i,j) = 720._dp - 778.5_dp * SQRT( grid%x(i)**2 + grid%y(j)**2)/ 750000._dp
      refgeo%Hs_grid( i,j) = surface_elevation( refgeo%Hi_grid( i,j), refgeo%Hb_grid( i,j), 0._dp)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_MISMIP_mod
  SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_A(   grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM A bumpy slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_A'
    INTEGER                                       :: i,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hs_grid( i,j) = 2000._dp - grid%x( i) * TAN( 0.5_dp * pi / 180._dp)
      refgeo%Hb_grid( i,j) = refgeo%Hs_grid( i,j) - 1000._dp + 500._dp * SIN( grid%x( i) * 2._dp * pi / C%ISMIP_HOM_L) * SIN( grid%y( j) * 2._dp * pi / C%ISMIP_HOM_L)
      refgeo%Hi_grid( i,j) = refgeo%Hs_grid( i,j) - refgeo%Hb_grid( i,j)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_A
  SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_B(   grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM B bumpy slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_B'
    INTEGER                                       :: i,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hs_grid( i,j) = 2000._dp - grid%x( i) * TAN( 0.5_dp * pi / 180._dp)
      refgeo%Hb_grid( i,j) = refgeo%Hs_grid( i,j) - 1000._dp + 500._dp * SIN( grid%x(i) * 2._dp * pi / C%ISMIP_HOM_L)
      refgeo%Hi_grid( i,j) = refgeo%Hs_grid( i,j) - refgeo%Hb_grid( i,j)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_B
  SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_CD(  grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM C/D bumpy slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_CD'
    INTEGER                                       :: i,j
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hs_grid( i,j) = 2000._dp - grid%x( i) * TAN( 0.1_dp * pi / 180._dp)
      refgeo%Hb_grid( i,j) = refgeo%Hs_grid( i,j) - 1000._dp
      refgeo%Hi_grid( i,j) = refgeo%Hs_grid( i,j) - refgeo%Hb_grid( i,j)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_CD
  SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_E(   grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM E Glacier d'Arolla geometry
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_E'
    INTEGER                                       :: i,j
    REAL(dp)                                      :: x,Hs,Hb
    INTEGER                                       :: ios,slides
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! To prevent compiler warnings from unused variables
    i = grid%nx
      
    ! Read data from external file
    IF (par%master) THEN
      
      OPEN( UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION = 'READ')
      DO i = 1, 51
        READ( UNIT = 1337, FMT=*, IOSTAT=ios) x, Hb, Hs, slides
        DO j = 1, refgeo%grid%ny
          refgeo%Hb_grid( i,j) = Hb
          refgeo%Hi_grid( i,j) = Hs - Hb
          refgeo%Hs_grid( i,j) = Hs
        END DO
        IF (ios /= 0) THEN
          CALL crash(' length of text file "' // TRIM( C%ISMIP_HOM_E_Arolla_filename) // '" should be 51 lines!')
        END IF
      END DO
      CLOSE( UNIT  = 1337)
      
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_E
  SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_F(   grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM A bumpy slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_ISMIP_HOM_F'
    INTEGER                                       :: i,j
    
    REAL(dp), PARAMETER                           :: H0    = 1000._dp
    REAL(dp), PARAMETER                           :: a0    = 100._dp
    REAL(dp), PARAMETER                           :: sigma = 10000._dp
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      refgeo%Hs_grid( i,j) = 5000._dp - grid%x( i) * TAN( 3._dp * pi / 180._dp)
      refgeo%Hb_grid( i,j) = refgeo%Hs_grid( i,j) - H0 + a0 * EXP( -((grid%x( i) - 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                                       + a0 * EXP( -((grid%x( i) - 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                                       + a0 * EXP( -((grid%x( i) - 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                                       + a0 * EXP( -((grid%x( i) - 0._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                                       + a0 * EXP( -((grid%x( i) - 0._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                                       + a0 * EXP( -((grid%x( i) - 0._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                                       + a0 * EXP( -((grid%x( i) + 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                                       + a0 * EXP( -((grid%x( i) + 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                                       + a0 * EXP( -((grid%x( i) + 1._dp * C%ISMIP_HOM_L)**2 + (grid%y( j) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2)
      refgeo%Hi_grid( i,j) = refgeo%Hs_grid( i,j) - refgeo%Hb_grid( i,j)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_ISMIP_HOM_F
  SUBROUTINE initialise_reference_geometry_idealised_grid_MISMIPplus(    grid, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The MISMIpplus fjord geometry
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_grid_MISMIPplus'
    INTEGER                                       :: i,j
    
    REAL(dp)                                      :: x,y,xtilde,Bx,By
    REAL(dp), PARAMETER                           :: B0     = -150._dp
    REAL(dp), PARAMETER                           :: B2     = -728.8_dp
    REAL(dp), PARAMETER                           :: B4     = 343.91_dp
    REAL(dp), PARAMETER                           :: B6     = -50.57_dp
    REAL(dp), PARAMETER                           :: xbar   = 300000._dp
    REAL(dp), PARAMETER                           :: fc     = 4000._dp
    REAL(dp), PARAMETER                           :: dc     = 500._dp
    REAL(dp), PARAMETER                           :: wc     = 24000._dp
    REAL(dp), PARAMETER                           :: zbdeep = -720._dp
    
    ! Add routine to path
    CALL init_routine( routine_name)
  
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      x = grid%x( i) + 400000._dp
      y = -40000._dp +  80000._dp * REAL(j-1,dp) / REAL(grid%ny-1,dp)
      xtilde = x / xbar
      Bx = B0 + (B2 * xtilde**2._dp) + (B4 * xtilde**4._dp) + (B6 * xtilde**6._dp)
      By = (dc / (1 + EXP(-2._dp*(y - wc)/fc))) + &
           (dc / (1 + EXP( 2._dp*(y + wc)/fc)))
      refgeo%Hi_grid( i,j) = 100._dp
      refgeo%Hb_grid( i,j) = MAX( Bx + By, zbdeep)
      refgeo%Hs_grid( i,j) = surface_elevation( refgeo%Hi_grid( i,j), refgeo%Hb_grid( i,j), 0._dp)
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_grid_MISMIPplus
  
  ! Set up a square grid for the idealised reference geometry (since it is now not provided externally)
  SUBROUTINE setup_idealised_geometry_grid( grid, region_name, dx)
    ! Set up a square grid for the idealised reference geometry (since it is now not provided externally)
  
    IMPLICIT NONE  
    
    ! In/output variables:
    TYPE(type_grid),                 INTENT(INOUT)     :: grid
    CHARACTER(LEN=3),                INTENT(IN)        :: region_name
    REAL(dp),                        INTENT(IN)        :: dx
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                     :: routine_name = 'setup_idealised_geometry_grid'
    INTEGER                                            :: nsx, nsy, i, j, n
    REAL(dp)                                           :: xmin, xmax, ymin, ymax, xmid, ymid
    REAL(dp), PARAMETER                                :: tol = 1E-9_dp
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Assign dummy values to suppress compiler warnings
    xmin = 0._dp
    xmax = 0._dp
    ymin = 0._dp
    ymax = 0._dp
    xmid = 0._dp
    ymid = 0._dp
    nsx  = 0
    nsy  = 0
    
    ! Allocate shared memory
    CALL allocate_shared_int_0D(                   grid%nx          , grid%wnx          )
    CALL allocate_shared_int_0D(                   grid%ny          , grid%wny          )
    CALL allocate_shared_dp_0D(                    grid%dx          , grid%wdx          )
    CALL allocate_shared_dp_0D(                    grid%xmin        , grid%wxmin        )
    CALL allocate_shared_dp_0D(                    grid%xmax        , grid%wxmax        )
    CALL allocate_shared_dp_0D(                    grid%ymin        , grid%wymin        )
    CALL allocate_shared_dp_0D(                    grid%ymax        , grid%wymax        )
    
    ! Let the Master do the work
    IF (par%master) THEN
    
      grid%dx = dx
  
      ! Resolution, domain size, and projection parameters for this region are determined from the config
      IF     (region_name == 'NAM') THEN
        xmin = C%xmin_NAM
        xmax = C%xmax_NAM
        ymin = C%ymin_NAM
        ymax = C%ymax_NAM
      ELSEIF (region_name == 'EAS') THEN
        xmin = C%xmin_EAS
        xmax = C%xmax_EAS
        ymin = C%ymin_EAS
        ymax = C%ymax_EAS
      ELSEIF (region_name == 'GRL') THEN
        xmin = C%xmin_GRL
        xmax = C%xmax_GRL
        ymin = C%ymin_GRL
        ymax = C%ymax_GRL
      ELSEIF (region_name == 'ANT') THEN
        xmin = C%xmin_ANT
        xmax = C%xmax_ANT
        ymin = C%ymin_ANT
        ymax = C%ymax_ANT
      END IF
      
      ! Determine the number of grid cells we can fit in this domain
      xmid = (xmax + xmin) / 2._dp
      ymid = (ymax + ymin) / 2._dp
      nsx  = FLOOR( (xmax - xmid) / grid%dx)
      nsy  = FLOOR( (ymax - ymid) / grid%dx)
      
      ! Small exceptions for very weird benchmark experiments
      IF (C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'SSA_icestream') nsx = 3
      IF (C%choice_refgeo_init_ANT == 'idealised' .AND. C%choice_refgeo_init_idealised == 'ISMIP_HOM_E')   nsx = 25
      
      grid%nx = 1 + 2*nsx
      grid%ny = 1 + 2*nsy
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Assign range to each processor
    CALL partition_list( grid%nx, par%i, par%n, grid%i1, grid%i2)
    CALL partition_list( grid%ny, par%i, par%n, grid%j1, grid%j2)
    
    ! Allocate shared memory for x and y
    CALL allocate_shared_dp_1D( grid%nx, grid%x, grid%wx)
    CALL allocate_shared_dp_1D( grid%ny, grid%y, grid%wy)
    
    ! Fill in x and y
    IF (par%master) THEN
      DO i = 1, grid%nx
        grid%x( i) = -nsx*grid%dx + (i-1)*grid%dx
      END DO
      DO j = 1, grid%ny
        grid%y( j) = -nsy*grid%dx + (j-1)*grid%dx
      END DO
      
      grid%xmin = MINVAL(grid%x)
      grid%xmax = MAXVAL(grid%x)
      grid%ymin = MINVAL(grid%y)
      grid%ymax = MAXVAL(grid%y)
    END IF ! IF (par%master) THEN
    CALL sync
      
    ! Tolerance; points lying within this distance of each other are treated as identical
    CALL allocate_shared_dp_0D( grid%tol_dist, grid%wtol_dist)
    IF (par%master) grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp
    
    ! Set up grid-to-vector translation tables
    CALL allocate_shared_int_0D(                   grid%n           , grid%wn           )
    IF (par%master) grid%n  = grid%nx * grid%ny
    CALL sync
    CALL allocate_shared_int_2D( grid%nx, grid%ny, grid%ij2n        , grid%wij2n        )
    CALL allocate_shared_int_2D( grid%n , 2      , grid%n2ij        , grid%wn2ij        )
    IF (par%master) THEN
      n = 0
      DO i = 1, grid%nx
        IF (MOD(i,2) == 1) THEN
          DO j = 1, grid%ny
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        ELSE
          DO j = grid%ny, 1, -1
            n = n+1
            grid%ij2n( i,j) = n
            grid%n2ij( n,:) = [i,j]
          END DO
        END IF
      END DO
    END IF
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 13)
    
  END SUBROUTINE setup_idealised_geometry_grid
  
  ! Fill in secondary data for the reference geometry (used to force mesh creation)
  SUBROUTINE calc_reference_geometry_secondary_data( grid, refgeo)
    ! Fill in secondary data for the reference geometry (used to force mesh creation)
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_grid),                INTENT(IN)    :: grid
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'calc_reference_geometry_secondary_data'
    INTEGER                                       :: i,j,ii,jj
    REAL(dp)                                      :: d2Hs_dx2, d2Hs_dxdy, d2Hs_dy2
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_2D(  grid%nx, grid%ny, refgeo%surf_curv  , refgeo%wsurf_curv  )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, refgeo%mask_land  , refgeo%wmask_land  )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, refgeo%mask_ocean , refgeo%wmask_ocean )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, refgeo%mask_ice   , refgeo%wmask_ice   )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, refgeo%mask_sheet , refgeo%wmask_sheet )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, refgeo%mask_shelf , refgeo%wmask_shelf )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, refgeo%mask_margin, refgeo%wmask_margin)
    CALL allocate_shared_int_2D( grid%nx, grid%ny, refgeo%mask_gl    , refgeo%wmask_gl    )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, refgeo%mask_cf    , refgeo%wmask_cf    )
    CALL allocate_shared_int_2D( grid%nx, grid%ny, refgeo%mask_coast , refgeo%wmask_coast )
    
    ! Calculate surface curvature
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      IF (i == 1 .OR. i == grid%nx .OR. j == 1 .OR. j == grid%ny) THEN
        d2Hs_dx2  = 0._dp
        d2Hs_dxdy = 0._dp
        d2Hs_dy2  = 0._dp
      ELSE
        d2Hs_dx2  = (refgeo%Hs_grid( i+1,j) + refgeo%Hs_grid( i-1,j) - 2._dp * refgeo%Hs_grid( i,j)) / (grid%dx**2)
        d2Hs_dxdy = (refgeo%Hs_grid( i+1,j+1) + refgeo%Hs_grid( i-1,j-1) - refgeo%Hs_grid( i+1,j-1) - refgeo%Hs_grid( i-1,j+1)) / (4._dp * grid%dx * grid%dx)
        d2Hs_dy2  = (refgeo%Hs_grid( i,j+1) + refgeo%Hs_grid( i,j-1) - 2._dp * refgeo%Hs_grid( i,j)) / (grid%dx**2)
      END IF
      
      refgeo%surf_curv( i,j) = MAX( -1E-6_dp, MIN( 1E-6_dp, SQRT( d2Hs_dx2**2 + d2Hs_dxdy**2 + d2Hs_dy2**2)))
      
    END DO
    END DO
    CALL sync
    
    ! Fill in masks
    
    ! Land/ocean
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      IF (is_floating( refgeo%Hi_grid( i,j), refgeo%Hb_grid( i,j), 0._dp)) THEN
        refgeo%mask_land(  i,j) = 1
      ELSE
        refgeo%mask_ocean( i,j) = 1
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Ice/sheet/shelf
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      IF (refgeo%Hi_grid( i,j) > 0._dp) THEN
      
        refgeo%mask_ice(  i,j) = 1
        
        IF (refgeo%mask_land( i,j) == 1) THEN
          refgeo%mask_sheet( i,j) = 1
        ELSE
          refgeo%mask_shelf( i,j) = 1
        END IF
        
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Transitions (margin, grounding line, calving front, coastline)
    DO i = grid%i1, grid%i2
    DO j = 1, grid%ny
      
      ! Ice next to non-ice equals ice margin
      IF (refgeo%mask_ice( i,j) == 1) THEN
        DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          IF (refgeo%mask_ice( ii,jj) == 0) THEN
            refgeo%mask_margin( i,j) = 1
          END IF
        END DO
        END DO
      END IF
      
      ! Sheet next to shelf equals grounding line
      IF (refgeo%mask_sheet( i,j) == 1) THEN
        DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          IF (refgeo%mask_shelf( ii,jj) == 1) THEN
            refgeo%mask_gl( i,j) = 1
          END IF
        END DO
        END DO
      END IF
      
      ! Ice next to open ocean equals calving front
      IF (refgeo%mask_ice( i,j) == 1) THEN
        DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          IF (refgeo%mask_ocean( ii,jj) == 1 .AND. refgeo%mask_ice( ii,jj) == 0) THEN
            refgeo%mask_cf( i,j) = 1
          END IF
        END DO
        END DO
      END IF
      
      ! Dry land next to open ocean equals coastline
      IF (refgeo%mask_land( i,j) == 1) THEN
        DO ii = MAX( 1, i-1), MIN( grid%nx, i+1)
        DO jj = MAX( 1, j-1), MIN( grid%ny, j+1)
          IF (refgeo%mask_ocean( ii,jj) == 1 .AND. refgeo%mask_ice( ii,jj) == 0) THEN
            refgeo%mask_coast( i,j) = 1
          END IF
        END DO
        END DO
      END IF
      
    END DO
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 10)
    
  END SUBROUTINE calc_reference_geometry_secondary_data
  
  ! Initialise a reference geometry according to an idealised world on the model mesh
  SUBROUTINE initialise_reference_geometry_idealised_mesh( mesh, refgeo, choice_refgeo_idealised)
    ! Initialise a reference geometry according to an idealised world
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    CHARACTER(LEN=256),             INTENT(IN)    :: choice_refgeo_idealised
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Allocate shared memory
    CALL allocate_shared_dp_1D( mesh%nV, refgeo%Hi, refgeo%wHi)
    CALL allocate_shared_dp_1D( mesh%nV, refgeo%Hb, refgeo%wHb)
    CALL allocate_shared_dp_1D( mesh%nV, refgeo%Hs, refgeo%wHs)
    
    IF     (choice_refgeo_idealised == 'flatearth') THEN
      ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments
      CALL initialise_reference_geometry_idealised_mesh_flatearth(     mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'Halfar') THEN
      ! The Halfar dome solution at t = 0
      CALL initialise_reference_geometry_idealised_mesh_Halfar(        mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'Bueler') THEN
      ! The Bueler dome solution at t = 0
      CALL initialise_reference_geometry_idealised_mesh_Bueler(        mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'SSA_icestream') THEN
      ! The SSA_icestream infinite slab on a flat slope
      CALL initialise_reference_geometry_idealised_mesh_SSA_icestream( mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'MISMIP_mod') THEN
      ! The MISMIP_mod cone-shaped island
      CALL initialise_reference_geometry_idealised_mesh_MISMIP_mod(    mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_A') THEN
      ! The ISMIP-HOM A bumpy slope
      CALL initialise_reference_geometry_idealised_mesh_ISMIP_HOM_A(   mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_B') THEN
      ! The ISMIP-HOM B bumpy slope
      CALL initialise_reference_geometry_idealised_mesh_ISMIP_HOM_B(   mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_C' .OR. &
            choice_refgeo_idealised == 'ISMIP_HOM_D') THEN
      ! The ISMIP-HOM C/D bumpy slope
      CALL initialise_reference_geometry_idealised_mesh_ISMIP_HOM_CD(  mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_E') THEN
      ! The ISMIP-HOM E Glacier d'Arolla geometry
      CALL initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E(   mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'ISMIP_HOM_F') THEN
      ! The ISMIP-HOM A bumpy slope
      CALL initialise_reference_geometry_idealised_mesh_ISMIP_HOM_F(   mesh, refgeo)
    ELSEIF (choice_refgeo_idealised == 'MISMIP+') THEN
      ! The MISMIP+ fjord geometry
      CALL initialise_reference_geometry_idealised_mesh_MISMIPplus(    mesh, refgeo)
    ELSE
      CALL crash('unknown choice_refgeo_idealised "' // TRIM( choice_refgeo_idealised) // '"!')
    END IF
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 3)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh
  SUBROUTINE initialise_reference_geometry_idealised_mesh_flatearth(     mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! Simply a flat, empty earth. Used for example in the EISMINT-1 benchmark experiments
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_flatearth'
    INTEGER                                       :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hi( vi) = 0._dp
      refgeo%Hb( vi) = 0._dp
      refgeo%Hs( vi) = 0._dp
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_flatearth
  SUBROUTINE initialise_reference_geometry_idealised_mesh_Halfar(        mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The Halfar dome solution at t = 0
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_Halfar'
    INTEGER                                       :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hi( vi) = Halfar_solution( mesh%V( vi,1), mesh%V( vi,2), C%start_time_of_run)
      refgeo%Hb( vi) = 0._dp
      refgeo%Hs( vi) = refgeo%Hi( vi)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_Halfar
  SUBROUTINE initialise_reference_geometry_idealised_mesh_Bueler(        mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The Bueler dome solution at t = 0
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_Bueler'
    INTEGER                                       :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hi( vi) = Bueler_solution( mesh%V( vi,1), mesh%V( vi,2), C%start_time_of_run)
      refgeo%Hb( vi) = 0._dp
      refgeo%Hs( vi) = refgeo%Hi( vi)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_Bueler
  SUBROUTINE initialise_reference_geometry_idealised_mesh_SSA_icestream( mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The SSA_icestream infinite slab on a flat slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_SSA_icestream'
    INTEGER                                       :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hi( vi) = 2000._dp
      refgeo%Hb( vi) = -0.001_dp * mesh%V( vi,1)
      refgeo%Hs( vi) = refgeo%Hb( vi) + refgeo%Hi( vi)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_SSA_icestream
  SUBROUTINE initialise_reference_geometry_idealised_mesh_MISMIP_mod(    mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The MISMIP_mod cone-shaped island
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_MISMIP_mod'
    INTEGER                                       :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO vi = mesh%vi1, mesh%vi2
      
      ! Create a nice circular ice shelf
      IF (SQRT(mesh%V( vi,1)**2 + mesh%V( vi,2)**2) < mesh%xmax * 0.95_dp) THEN
        refgeo%Hi( vi) = 100._dp
      ELSE
        refgeo%Hi( vi) = 0._dp
      END IF
      
      refgeo%Hb( vi) = 720._dp - 778.5_dp * SQRT( mesh%V( vi,1)**2 + mesh%V( vi,2)**2)/ 750000._dp
      refgeo%Hs( vi) = surface_elevation( refgeo%Hi( vi), refgeo%Hb( vi), 0._dp)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_MISMIP_mod
  SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_A(   mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM A bumpy slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_A'
    INTEGER                                       :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hs( vi) = 2000._dp - mesh%V( vi,1) * TAN( 0.5_dp * pi / 180._dp)
      refgeo%Hb( vi) = refgeo%Hs( vi) - 1000._dp + 500._dp * SIN( mesh%V( vi,1) * 2._dp * pi / C%ISMIP_HOM_L) * SIN( mesh%V( vi,2) * 2._dp * pi / C%ISMIP_HOM_L)
      refgeo%Hi( vi) = refgeo%Hs( vi) - refgeo%Hb( vi)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_A
  SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_B(   mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM B bumpy slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_B'
    INTEGER                                       :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hs( vi) = 2000._dp - mesh%V( vi,1) * TAN( 0.5_dp * pi / 180._dp)
      refgeo%Hb( vi) = refgeo%Hs( vi) - 1000._dp + 500._dp * SIN( mesh%V( vi,1) * 2._dp * pi / C%ISMIP_HOM_L)
      refgeo%Hi( vi) = refgeo%Hs( vi) - refgeo%Hb( vi)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_B
  SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_CD(  mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM C/D bumpy slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_CD'
    INTEGER                                       :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
      
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hs( vi) = 2000._dp - mesh%V( vi,1) * TAN( 0.1_dp * pi / 180._dp)
      refgeo%Hb( vi) = refgeo%Hs( vi) - 1000._dp
      refgeo%Hi( vi) = refgeo%Hs( vi) - refgeo%Hb( vi)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_CD
  SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E(   mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM E Glacier d'Arolla geometry
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E'
    
    ! DENK DROM
    REAL(dp) :: dp_dummy
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! DENK DROM
    dp_dummy = mesh%V( 1,1)
    dp_dummy = refgeo%Hi( 1)
    CALL crash('FIXME!')
    
!    ! Local variables
!    INTEGER                                       :: vi
!    REAL(dp)                                      :: x,Hs,Hb
!    INTEGER                                       :: ios,slides
!    
!    ! To prevent compiler warnings from unused variables
!    i = grid%nx
!      
!    ! Read data from external file
!    IF (par%master) THEN
!      
!      OPEN( UNIT = 1337, FILE=C%ISMIP_HOM_E_Arolla_filename, ACTION = 'READ')
!      DO i = 1, 51
!        READ( UNIT = 1337, FMT=*, IOSTAT=ios) x, Hb, Hs, slides
!        DO j = 1, refgeo%grid%ny
!          refgeo%Hb( vi) = Hb
!          refgeo%Hi( vi) = Hs - Hb
!          refgeo%Hs( vi) = Hs
!        END DO
!        IF (ios /= 0) THEN
!          WRITE(0,*) ' initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E - ERROR: length of text file "', TRIM(C%ISMIP_HOM_E_Arolla_filename), '" should be 51 lines!'
!          CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
!        END IF
!      END DO
!      CLOSE( UNIT  = 1337)
!      
!    END IF ! IF (par%master) THEN
!    CALL sync
!    
!    ! Finalise routine path
!    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_E
  SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_F(   mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The ISMIP-HOM A bumpy slope
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_ISMIP_HOM_F'
    INTEGER                                       :: vi
    
    REAL(dp), PARAMETER                           :: H0    = 1000._dp
    REAL(dp), PARAMETER                           :: a0    = 100._dp
    REAL(dp), PARAMETER                           :: sigma = 10000._dp
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hs( vi) = 5000._dp - mesh%V( vi,1) * TAN( 3._dp * pi / 180._dp)
      refgeo%Hb( vi) = refgeo%Hs( vi) - H0 + a0 * EXP( -((mesh%V( vi,1) - 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                           + a0 * EXP( -((mesh%V( vi,1) - 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                           + a0 * EXP( -((mesh%V( vi,1) - 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                           + a0 * EXP( -((mesh%V( vi,1) - 0._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                           + a0 * EXP( -((mesh%V( vi,1) - 0._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                           + a0 * EXP( -((mesh%V( vi,1) - 0._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                           + a0 * EXP( -((mesh%V( vi,1) + 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 1._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                           + a0 * EXP( -((mesh%V( vi,1) + 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) - 0._dp * C%ISMIP_HOM_L)**2) / sigma**2) &
                                           + a0 * EXP( -((mesh%V( vi,1) + 1._dp * C%ISMIP_HOM_L)**2 + (mesh%V( vi,2) + 1._dp * C%ISMIP_HOM_L)**2) / sigma**2)
      refgeo%Hi( vi) = refgeo%Hs( vi) - refgeo%Hb( vi)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_ISMIP_HOM_F
  SUBROUTINE initialise_reference_geometry_idealised_mesh_MISMIPplus(    mesh, refgeo)
    ! Initialise reference geometry according to an idealised world
    !
    ! The MISMIpplus fjord geometry
     
    IMPLICIT NONE
      
    ! In/output variables:
    TYPE(type_mesh),                INTENT(IN)    :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'initialise_reference_geometry_idealised_mesh_MISMIPplus'
    INTEGER                                       :: vi
    
    REAL(dp)                                      :: x,y,xtilde,Bx,By
    REAL(dp), PARAMETER                           :: B0     = -150._dp
    REAL(dp), PARAMETER                           :: B2     = -728.8_dp
    REAL(dp), PARAMETER                           :: B4     = 343.91_dp
    REAL(dp), PARAMETER                           :: B6     = -50.57_dp
    REAL(dp), PARAMETER                           :: xbar   = 300000._dp
    REAL(dp), PARAMETER                           :: fc     = 4000._dp
    REAL(dp), PARAMETER                           :: dc     = 500._dp
    REAL(dp), PARAMETER                           :: wc     = 24000._dp
    REAL(dp), PARAMETER                           :: zbdeep = -720._dp
    
    ! Add routine to path
    CALL init_routine( routine_name)
  
    DO vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1) + 400000._dp
      y = mesh%V( vi,2)
      xtilde = x / xbar
      Bx = B0 + (B2 * xtilde**2._dp) + (B4 * xtilde**4._dp) + (B6 * xtilde**6._dp)
      By = (dc / (1 + EXP(-2._dp*(y - wc)/fc))) + &
           (dc / (1 + EXP( 2._dp*(y + wc)/fc)))
      refgeo%Hi( vi) = 100._dp
      refgeo%Hb( vi) = MAX( Bx + By, zbdeep)
      refgeo%Hs( vi) = surface_elevation( refgeo%Hi( vi), refgeo%Hb( vi), 0._dp)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE initialise_reference_geometry_idealised_mesh_MISMIPplus

  ! Map init and PD references data from their supplied grids to the model mesh
  SUBROUTINE map_reference_geometries_to_mesh( region, mesh)
    ! Map the initial, present-day, and GIAeq reference geometries from their original
    ! square grids to the model mesh.
    !
    ! Since calculating remapping operators can take some time, and since the three
    ! square grids are often identical, some time can be saved by calculating the
    ! operators only once.
    
    IMPLICIT NONE
  
    ! Input and output variables
    TYPE(type_model_region),        INTENT(INOUT) :: region
    TYPE(type_mesh),                INTENT(INOUT) :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'map_reference_geometries_to_mesh'
    CHARACTER(LEN=256)                            :: choice_refgeo_init, choice_refgeo_PD, choice_refgeo_GIAeq
    LOGICAL                                       :: did_remap_init, did_remap_PD, did_remap_GIAeq, do_reuse_init_map, do_reuse_PD_map
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Determine parameters for this region
    IF     (region%name == 'NAM') THEN
      choice_refgeo_init    = C%choice_refgeo_init_NAM
      choice_refgeo_PD      = C%choice_refgeo_PD_NAM
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_NAM
    ELSEIF (region%name == 'EAS') THEN
      choice_refgeo_init    = C%choice_refgeo_init_EAS
      choice_refgeo_PD      = C%choice_refgeo_PD_EAS
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_EAS
    ELSEIF (region%name == 'GRL') THEN
      choice_refgeo_init    = C%choice_refgeo_init_GRL
      choice_refgeo_PD      = C%choice_refgeo_PD_GRL
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_GRL
    ELSEIF (region%name == 'ANT') THEN
      choice_refgeo_init    = C%choice_refgeo_init_ANT
      choice_refgeo_PD      = C%choice_refgeo_PD_ANT
      choice_refgeo_GIAeq   = C%choice_refgeo_GIAeq_ANT
    END IF
    
    ! Initial ice-sheet geometry
    ! ==========================
    
    IF (par%master) WRITE(0,*) '  Mapping initial reference geometry to the mesh...'
    
    did_remap_init = .FALSE.
    
    IF     (choice_refgeo_init == 'idealised') THEN
      ! For idealised geometries, calculate the exact solution directly on the mesh instead of remapping it from the grid
      
      CALL initialise_reference_geometry_idealised_mesh( mesh, region%refgeo_init, C%choice_refgeo_init_idealised)
      
    ELSEIF (choice_refgeo_init == 'realistic') THEN
      ! For realistic geometries, remap the data from the grid
      
      CALL calc_remapping_operator_grid2mesh( region%refgeo_init%grid, mesh)
      CALL map_reference_geometry_to_mesh( mesh, region%refgeo_init)
      did_remap_init = .TRUE.
      
    ELSE
      CALL crash('unknown choice_refgeo_init "' // TRIM( choice_refgeo_init) // '"!')
    END IF
    
    ! Present-day ice-sheet geometry
    ! ==============================
    
    IF (par%master) WRITE(0,*) '  Mapping present-day reference geometry to the mesh...'
    
    did_remap_PD = .FALSE.
    
    IF     (choice_refgeo_PD == 'idealised') THEN
      ! For idealised geometries, calculate the exact solution directly on the mesh instead of remapping it from the grid
      
      CALL initialise_reference_geometry_idealised_mesh( mesh, region%refgeo_PD, C%choice_refgeo_PD_idealised)
      
    ELSEIF (choice_refgeo_PD == 'realistic') THEN
      ! For realistic geometries, remap the data from the grid
      
      ! Check if we can re-use the remapping arrays from the initial geometry
      do_reuse_init_map = .FALSE.
      IF (did_remap_init) THEN
        IF (region%refgeo_PD%grid%xmin == region%refgeo_init%grid%xmin .AND. &
            region%refgeo_PD%grid%xmax == region%refgeo_init%grid%xmax .AND. &
            region%refgeo_PD%grid%ymin == region%refgeo_init%grid%ymin .AND. &
            region%refgeo_PD%grid%ymax == region%refgeo_init%grid%ymax .AND. &
            region%refgeo_PD%grid%dx   == region%refgeo_init%grid%dx   .AND. &
            region%refgeo_PD%grid%nx   == region%refgeo_init%grid%nx   .AND. &
            region%refgeo_PD%grid%ny   == region%refgeo_init%grid%ny) THEN
          do_reuse_init_map = .TRUE.
        END IF
      END IF
      IF (do_reuse_init_map) THEN
        CALL MatDuplicate( region%refgeo_init%grid%M_map_grid2mesh, MAT_COPY_VALUES, region%refgeo_PD%grid%M_map_grid2mesh, perr)
      ELSE
        CALL calc_remapping_operator_grid2mesh( region%refgeo_PD%grid, mesh)
      END IF
      
      CALL map_reference_geometry_to_mesh( mesh, region%refgeo_PD)
      did_remap_PD = .TRUE.
      
    ELSE
      CALL crash('unknown choice_refgeo_PD "' // TRIM( choice_refgeo_PD) // '"!')
    END IF
    
    ! GIA equilibrium ice-sheet geometry
    ! ==================================
    
    IF (par%master) WRITE(0,*) '  Mapping GIA equilibrium reference geometry to the mesh...'
    
    did_remap_GIAeq = .FALSE.
    
    IF     (choice_refgeo_GIAeq == 'idealised') THEN
      ! For idealised geometries, calculate the exact solution directly on the mesh instead of remapping it from the grid
      
      CALL initialise_reference_geometry_idealised_mesh( mesh, region%refgeo_GIAeq, C%choice_refgeo_GIAeq_idealised)
      
    ELSEIF (choice_refgeo_GIAeq == 'realistic') THEN
      ! For realistic geometries, remap the data from the grid
      
      ! Check if we can re-use the remapping arrays from the initial/PD geometry
      do_reuse_init_map = .FALSE.
      IF (did_remap_init) THEN
        IF (region%refgeo_GIAeq%grid%xmin == region%refgeo_init%grid%xmin .AND. &
            region%refgeo_GIAeq%grid%xmax == region%refgeo_init%grid%xmax .AND. &
            region%refgeo_GIAeq%grid%ymin == region%refgeo_init%grid%ymin .AND. &
            region%refgeo_GIAeq%grid%ymax == region%refgeo_init%grid%ymax .AND. &
            region%refgeo_GIAeq%grid%dx   == region%refgeo_init%grid%dx   .AND. &
            region%refgeo_GIAeq%grid%nx   == region%refgeo_init%grid%nx   .AND. &
            region%refgeo_GIAeq%grid%ny   == region%refgeo_init%grid%ny) THEN
          do_reuse_init_map = .TRUE.
        END IF
      END IF
      do_reuse_PD_map = .FALSE.
      IF (did_remap_init) THEN
        IF (region%refgeo_GIAeq%grid%xmin == region%refgeo_PD%grid%xmin .AND. &
            region%refgeo_GIAeq%grid%xmax == region%refgeo_PD%grid%xmax .AND. &
            region%refgeo_GIAeq%grid%ymin == region%refgeo_PD%grid%ymin .AND. &
            region%refgeo_GIAeq%grid%ymax == region%refgeo_PD%grid%ymax .AND. &
            region%refgeo_GIAeq%grid%dx   == region%refgeo_PD%grid%dx   .AND. &
            region%refgeo_GIAeq%grid%nx   == region%refgeo_PD%grid%nx   .AND. &
            region%refgeo_GIAeq%grid%ny   == region%refgeo_PD%grid%ny) THEN
          do_reuse_PD_map = .TRUE.
        END IF
      END IF
      IF     (do_reuse_init_map) THEN
        CALL MatDuplicate( region%refgeo_init%grid%M_map_grid2mesh, MAT_COPY_VALUES, region%refgeo_GIAeq%grid%M_map_grid2mesh, perr)
      ELSEIF (do_reuse_PD_map) THEN
        CALL MatDuplicate( region%refgeo_PD%grid%M_map_grid2mesh  , MAT_COPY_VALUES, region%refgeo_GIAeq%grid%M_map_grid2mesh, perr)
      ELSE
        CALL calc_remapping_operator_grid2mesh( region%refgeo_GIAeq%grid, mesh)
      END IF
      
      CALL map_reference_geometry_to_mesh( mesh, region%refgeo_GIAeq)
      did_remap_GIAeq = .TRUE.
      
    ELSE
      CALL crash('unknown choice_refgeo_GIAeq "' // TRIM( choice_refgeo_GIAeq) // '"!')
    END IF
    
    ! Clean up after yourself
    IF (did_remap_init ) CALL deallocate_remapping_operators_grid2mesh( region%refgeo_init%grid )
    IF (did_remap_PD   ) CALL deallocate_remapping_operators_grid2mesh( region%refgeo_PD%grid   )
    IF (did_remap_GIAeq) CALL deallocate_remapping_operators_grid2mesh( region%refgeo_GIAeq%grid)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 9)
    
  END SUBROUTINE map_reference_geometries_to_mesh
  SUBROUTINE map_reference_geometry_to_mesh( mesh, refgeo)
    ! Map data for a single reference geometry from its original square grid to the model mesh.
    
    IMPLICIT NONE
  
    ! Input and output variables
    TYPE(type_mesh),                INTENT(INOUT) :: mesh
    TYPE(type_reference_geometry),  INTENT(INOUT) :: refgeo
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'map_reference_geometry_to_mesh'
    INTEGER                                       :: vi
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    ! Map PD data to the mesh
    CALL map_grid2mesh_2D( refgeo%grid, mesh, refgeo%Hi_grid, refgeo%Hi)
    CALL map_grid2mesh_2D( refgeo%grid, mesh, refgeo%Hb_grid, refgeo%Hb)
    
    DO vi = mesh%vi1, mesh%vi2
      refgeo%Hs( vi) = surface_elevation( refgeo%Hi( vi), refgeo%Hb( vi), 0._dp)
    END DO
    CALL sync
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
  
  END SUBROUTINE map_reference_geometry_to_mesh

  ! Analytical solutions used to initialise some benchmark experiments
  FUNCTION Halfar_solution( x, y, t) RESULT(H)
    ! Describes an ice-sheet at time t (in years) conforming to the Halfar similarity function 
    ! with dome thickness H0 and margin radius R0 at t0. Used to initialise the model
    ! for the Halfar solution test run
    
    IMPLICIT NONE
    
    ! Input variables
    REAL(dp), INTENT(IN) :: x  ! x coordinate [m]
    REAL(dp), INTENT(IN) :: y  ! y coordinate [m]
    REAL(dp), INTENT(IN) :: t  ! Time from t0 [years]
    
    ! Result
    REAL(dp)             :: H  ! Ice thickness at [x,y] at t=0 [m]
    
    ! Local variables
    REAL(dp) :: A_flow, rho, g, Gamma, t0, r, f1, f2, f3, tp
    
    REAL(dp), PARAMETER :: H0 = 5000._dp   ! Ice dome thickness at t=0 [m]
    REAL(dp), PARAMETER :: R0 = 300000._dp ! Ice margin radius  at t=0 [m]
    
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
  FUNCTION Bueler_solution( x, y, t) RESULT(H)
    ! Describes an ice-sheet at time t (in years) conforming to the Bueler solution
    ! with dome thickness H0 and margin radius R0 at t0, with a surface mass balance
    ! determined by lambda. Used to intialise the model for the Bueler solution test run
    
    IMPLICIT NONE
    
    ! Input variables
    REAL(dp), INTENT(IN) :: x       ! x coordinate [m]
    REAL(dp), INTENT(IN) :: y       ! y coordinate [m]
    REAL(dp), INTENT(IN) :: t       ! Time from t0 [years]
    
    ! Result
    REAL(dp)             :: H  ! Ice thickness at [x,y] at t=0 [m]
    
    ! Local variables
    REAL(dp) :: A_flow, rho, g, n, alpha, beta, Gamma, f1, f2, t0, tp, f3, f4
    
    REAL(dp), PARAMETER :: H0     = 3000._dp    ! Ice dome thickness at t=0 [m]
    REAL(dp), PARAMETER :: R0     = 500000._dp  ! Ice margin radius  at t=0 [m]
    REAL(dp), PARAMETER :: lambda = 5.0_dp      ! Mass balance parameter
  
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
    
    !tp = (t * sec_per_year) + t0; % Actual equation needs t in seconds from zero , but we want to supply t in years from t0
    tp = t * sec_per_year
    
    f1 = (tp / t0)**(-alpha)
    f2 = (tp / t0)**(-beta)
    f3 = SQRT( (x**2._dp) + (y**2._dp) )/R0
    f4 = MAX(0._dp, 1._dp - (f2*f3)**((n+1._dp)/n))
    H = H0 * f1 * f4**(n/((2._dp*n)+1._dp))
    
    !M = (lambda / tp) * H * sec_per_year
  
  END FUNCTION Bueler_solution

END MODULE reference_fields_module
