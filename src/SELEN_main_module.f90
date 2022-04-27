MODULE SELEN_main_module

  ! The "coupling" module that contains all the stuff for coupling SELEN to IMAU-ICE.
  ! The two main routines that are called from the IMAU-ICE program are "initialise_SELEN"
  ! and "run_SELEN". initialise_SELEN does all the pre-work of allocating memory, reading
  ! the external LJ and MJ values, initialising the TABOO GIA model (for the Love numbers),
  ! setting up the mapping between the four regional ice-model grids and the SELEN global
  ! grid, initialising (either by reading or by calculating) the spherical harmonics
  ! transfer functions, initialising the ice loading history, and creating a SELEN output
  ! NetCDF file. run_SELEN updates the 80 kyr ice loading history (regional and global),
  ! runs SELEN to solve the SLE, maps the results back to the ice model grid, and writes
  ! the (global) results to the output NetCDF file.

  USE mpi
  USE configuration_module,              ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module,                 ONLY: pi
  USE parallel_module,                   ONLY: par, sync, ierr, partition_list, &
                                               allocate_shared_int_0D, allocate_shared_dp_0D, &
                                               allocate_shared_int_1D, allocate_shared_dp_1D, &
                                               allocate_shared_int_2D, allocate_shared_dp_2D, &
                                               allocate_shared_int_3D, allocate_shared_dp_3D, &
                                               allocate_shared_complex_1D, deallocate_shared
  USE data_types_module,                 ONLY: type_model_region, type_SELEN_global
  USE netcdf_module,                     ONLY: inquire_SELEN_global_topo_file, read_SELEN_global_topo_file, &
                                               create_SELEN_output_file
  USE mesh_mapping_module,               ONLY: calc_remapping_operator_mesh2grid, map_mesh2grid_2D
  USE general_ice_model_data_module,     ONLY: is_floating
  USE SELEN_mapping_module,              ONLY: create_GIA_grid_to_SELEN_maps, map_GIA_grid_to_SELEN
  USE SELEN_harmonics_module,            ONLY: DOM, PLEG, HARMO, PLMBAR_MOD
  USE SELEN_taboo_hp_module,             ONLY: taboo_hp_model

  IMPLICIT NONE

CONTAINS

! == Initialisation ==
! ====================

! == Initialise all the SELEN stuff.
  SUBROUTINE initialise_SELEN( SELEN, NAM, EAS, GRL, ANT, version_number)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN
    TYPE(type_model_region),             INTENT(INOUT) :: NAM, EAS, GRL, ANT
    CHARACTER(LEN=256),                  INTENT(IN)    :: version_number

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SELEN'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) ''
    IF (par%master) WRITE(0,*) ' Initialising SELEN...'

    ! Initialise the global irregular mesh and reference topography, ice loading & ocean function
    IF (par%master) WRITE(0,*) '  Initialising SELEN global mesh and topography...'
    CALL initialise_SELEN_global_data( SELEN)

    ! Initialise the reference load and topography for the four ice-model regions
    IF (par%master) WRITE(0,*) '  Initialising SELEN regional reference fields...'
    IF (C%do_NAM) CALL initialise_reference_fields_regional( NAM)
    IF (C%do_EAS) CALL initialise_reference_fields_regional( EAS)
    IF (C%do_GRL) CALL initialise_reference_fields_regional( GRL)
    IF (C%do_ANT) CALL initialise_reference_fields_regional( ANT)

    ! Link ice-model regional GIA grid points to SELEN global grid points.
    ! Mind the sequence! GRL first, then EAS, then NAM.
    IF (C%do_GRL) CALL create_GIA_grid_to_SELEN_maps( SELEN, GRL, 3)
    IF (C%do_EAS) CALL create_GIA_grid_to_SELEN_maps( SELEN, EAS, 2)
    IF (C%do_NAM) CALL create_GIA_grid_to_SELEN_maps( SELEN, NAM, 1)
    IF (C%do_ANT) CALL create_GIA_grid_to_SELEN_maps( SELEN, ANT, 4)

    ! Create regional pixel lists
    IF (C%do_NAM) CALL create_icemodel_regional_pixel_list( SELEN, NAM, 1)
    IF (C%do_EAS) CALL create_icemodel_regional_pixel_list( SELEN, EAS, 2)
    IF (C%do_GRL) CALL create_icemodel_regional_pixel_list( SELEN, GRL, 3)
    IF (C%do_ANT) CALL create_icemodel_regional_pixel_list( SELEN, ANT, 4)

    ! Create list of SELEN global pixels that belong to ice model regions
    CALL create_icemodel_pixel_list( SELEN)

    ! Map present-day topography and ice loading from the regional GIA grids to the SELEN global mesh
    CALL adapt_SELEN_global_data_to_ice_model( SELEN, NAM, EAS, GRL, ANT)

    ! Read Love numbers from the provided external file
    IF (par%master) WRITE(0,*) '  Reading Love numbers...'
    CALL read_Love_numbers( SELEN)

    ! Initialise the spherical harmonics files (either read from files, or create from scratch)
    IF (par%master) WRITE(0,*) '  Initialising spherical harmonics...'
    IF (C%do_NAM) CALL initialise_spherical_harmonics( SELEN, NAM, version_number)
    IF (C%do_EAS) CALL initialise_spherical_harmonics( SELEN, EAS, version_number)
    IF (C%do_GRL) CALL initialise_spherical_harmonics( SELEN, GRL, version_number)
    IF (C%do_ANT) CALL initialise_spherical_harmonics( SELEN, ANT, version_number)

    ! Initialise ALF and LONG_TABLE
    IF (par%master) WRITE(0,*) '  Initialising ALF and LONG_TABLE...'
    CALL initialise_ALF_and_LONG_TABLE( SELEN)

    ! Initialise ice loading history
    IF (par%master) WRITE(0,*) '  Initialising ice loading history...'
    CALL initialise_ice_loading_history( SELEN, NAM, EAS, GRL, ANT)

    ! Calculating Green's functions (from load Love numbers)
    IF (par%master) WRITE(0,*) '  Calling TABOO and calculate Green functions...'
    IF (par%master) CALL taboo_hp_model

    ! Allocate memory for SELEN end results
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%MEM_S,       SELEN%wMEM_S      )
    CALL allocate_shared_complex_1D( C%SELEN_jmax * (C%SELEN_irreg_time_n+1), SELEN%MEM_U,       SELEN%wMEM_U      )
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV,                           SELEN%Hi_glob,     SELEN%wHi_glob    )
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV,                           SELEN%Hi_rel_glob, SELEN%wHi_rel_glob)
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV,                           SELEN%U_glob,      SELEN%wU_glob     )
    CALL allocate_shared_dp_1D(      SELEN%mesh%nV,                           SELEN%N_glob,      SELEN%wN_glob     )
    CALL allocate_shared_int_1D(     SELEN%mesh%nV,                           SELEN%of_glob,     SELEN%wof_glob    )

    ! ! Create NetCDF output file
    CALL create_SELEN_output_file( SELEN)

    IF (par%master) WRITE(0,*) ' Finished initialising SELEN.'
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=19)

  END SUBROUTINE initialise_SELEN
  SUBROUTINE read_Love_numbers( SELEN)
    ! Read the Love numbers from the provided external file

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_Love_numbers'
    INTEGER                                            :: j
    LOGICAL                                            :: exists

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Read the MJ and LJ values from a binary file
    ! ============================================

    INQUIRE( FILE = TRIM(C%selen_dir)//'/'//TRIM(C%SELEN_LMJ_VALUES_filename), EXIST=exists)
    IF (.NOT. exists) THEN
      CALL crash('LMJ file "' // TRIM( C%SELEN_LMJ_VALUES_filename) // '" not found!')
    END IF

    CALL allocate_shared_int_1D( C%SELEN_jmax, SELEN%MJ_VAL, SELEN%wMJ_VAL)
    CALL allocate_shared_int_1D( C%SELEN_jmax, SELEN%LJ_VAL, SELEN%wLJ_VAL)

    IF (par%master) THEN
      OPEN(   UNIT = 6004, FILE = TRIM(C%selen_dir)//'/'//TRIM(C%SELEN_LMJ_VALUES_filename), STATUS = 'OLD',FORM = 'UNFORMATTED')
      DO j = 1, C%SELEN_jmax
        READ( UNIT = 6004) SELEN%LJ_VAL(j), SELEN%MJ_VAL(j)
      END DO
      CLOSE(  UNIT = 6004)
    END IF
    CALL sync

    CALL partition_list( C%SELEN_jmax, par%i, par%n, C%SELEN_j1, C%SELEN_j2)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=2)

  END SUBROUTINE read_Love_numbers

  ! == Initialise ice loading history
  SUBROUTINE initialise_ice_loading_history( SELEN, NAM, EAS, GRL, ANT)
    ! Allocate memory for the global and regional ice loading histories

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN
    TYPE(type_model_region),             INTENT(INOUT) :: NAM, EAS, GRL, ANT

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_loading_history'
    INTEGER                                            :: ki

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory:
    CALL allocate_shared_dp_1D( SELEN%mesh%nV * (C%SELEN_irreg_time_n+1), SELEN%dice_loading_history_irreg_glob, SELEN%wdice_loading_history_irreg_glob)
    SELEN%ice_loading_history_irreg_glob( 1:SELEN%mesh%nV,0:C%SELEN_irreg_time_n) => SELEN%dice_loading_history_irreg_glob

    ! Assume the world is in a steady state
    DO ki = 0, C%SELEN_irreg_time_n
      SELEN%ice_loading_history_irreg_glob( C%SELEN_i1:C%SELEN_i2,ki) = SELEN%load_ref( C%SELEN_i1:C%SELEN_i2)
    END DO
    CALL sync

    ! Initialise for the ice model regions
    IF (C%do_NAM) CALL initialise_ice_loading_history_region( NAM)
    IF (C%do_EAS) CALL initialise_ice_loading_history_region( EAS)
    IF (C%do_GRL) CALL initialise_ice_loading_history_region( GRL)
    IF (C%do_ANT) CALL initialise_ice_loading_history_region( ANT)

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=3)

  END SUBROUTINE initialise_ice_loading_history
  SUBROUTINE initialise_ice_loading_history_region( region)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ice_loading_history_region'
    INTEGER                                            :: k, ki

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory:
    CALL allocate_shared_dp_3D( region%grid_GIA%nx, region%grid_GIA%ny, C%SELEN_reg_time_n, region%SELEN%ice_loading_history_reg_sq,   region%SELEN%wice_loading_history_reg_sq  )
    CALL allocate_shared_dp_3D( region%grid_GIA%nx, region%grid_GIA%ny, C%SELEN_irreg_time_n, region%SELEN%ice_loading_history_irreg_sq, region%SELEN%wice_loading_history_irreg_sq)

    ! Assume the world is in a steady state
    DO k = 1, C%SELEN_reg_time_n
      region%SELEN%ice_loading_history_reg_sq(   region%grid_GIA%i1:region%grid_GIA%i2,:,k ) = region%SELEN%load_ref( region%grid_GIA%i1:region%grid_GIA%i2,:)
    END DO
    DO ki = 1, C%SELEN_irreg_time_n
      region%SELEN%ice_loading_history_irreg_sq( region%grid_GIA%i1:region%grid_GIA%i2,:,ki) = region%SELEN%load_ref( region%grid_GIA%i1:region%grid_GIA%i2,:)
    END DO
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=2)

  END SUBROUTINE initialise_ice_loading_history_region

  ! == Initialise ALF and LONG_TABLE
  SUBROUTINE initialise_ALF_and_LONG_TABLE( SELEN)
    ! ALF contains some data regarding spherical harmonics of the anchor points
    ! LONG_TABLE contains some data regarding the entire global SELEN grid
    ! However, since it's been indexed starting at zero (WHY?!), storing it as
    ! shared memory is a bit tricky. This solved by creating a "wrapper function"
    ! to access the shared memory data, where the function simply adds 1 to the queried index.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_ALF_and_LONG_TABLE'
    INTEGER                                            :: i, i1, i2, li

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_2D(       C%SELEN_jmax,     SELEN%mesh%nanc, SELEN%ALF,         SELEN%wALF       )
    CALL allocate_shared_complex_1D( (C%SELEN_n_harmonics+1) * SELEN%mesh%nV,    SELEN%dLONG_TABLE, SELEN%wdLONG_TABLE)

    SELEN%LONG_TABLE(0:C%SELEN_n_harmonics,1:SELEN%mesh%nV) => SELEN%dLONG_TABLE

    ! Calculate ALF
    CALL partition_list( SELEN%mesh%nanc, par%i, par%n, i1, i2)
    DO i = i1, i2
      CALL PLMBAR_MOD( SELEN, C%SELEN_n_harmonics, SELEN%mesh%ancplist_lat( i), SELEN%ALF( :,i))
    END DO ! DO i = i1, i2
    CALL sync

    ! Calculate LONG_TABLE
    DO i = C%SELEN_i1, C%SELEN_i2
      DO li = 0, C%SELEN_n_harmonics
        SELEN%LONG_TABLE( li,i) = cmplx(COS( dble( li) * SELEN%mesh%lon( i) * (pi / 180._dp)), &
                                        SIN( dble( li) * SELEN%mesh%lon( i) * (pi / 180._dp))  )
      END DO
    END DO ! DO i = C%SELEN_i1, C%SELEN_i2
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=2)

  END SUBROUTINE initialise_ALF_and_LONG_TABLE

  ! == Initialise regional and global reference fields
  SUBROUTINE initialise_SELEN_global_data( SELEN)
    ! Read the SELEN global irregular mesh and topography on that mesh from the provided external file.
    ! Map PD ice loading from the ice model mesh(es) to the SELEN global mesh.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SELEN_global_data'
    INTEGER                                            :: i1, i2, iap

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for mesh size variables
    CALL allocate_shared_int_0D( SELEN%mesh%nV,     SELEN%mesh%wnV    )
    CALL allocate_shared_int_0D( SELEN%mesh%nTri,   SELEN%mesh%wnTri  )
    CALL allocate_shared_int_0D( SELEN%mesh%nC_mem, SELEN%mesh%wnC_mem)

    ! Inquire into the provided NetCDF file, obtain mesh size
    IF (par%master) SELEN%netcdf_topo%filename = TRIM(C%SELEN_dir)//'/'//TRIM(C%SELEN_global_topo_filename)
    IF (par%master) WRITE(0,*) '  Reading SELEN global mesh and topography data from file "', TRIM(SELEN%netcdf_topo%filename), '"...'
    IF (par%master) CALL inquire_SELEN_global_topo_file( SELEN)
    CALL sync
    C%SELEN_alfa = SQRT(4._dp / DBLE( SELEN%mesh%nV))
    CALL partition_list( SELEN%mesh%nV, par%i, par%n, C%SELEN_i1, C%SELEN_i2)

    ! Allocate memory for the mesh and the global topography data
    CALL allocate_shared_dp_2D(  SELEN%mesh%nV,   3,                 SELEN%mesh%V,          SELEN%mesh%wV         )
    CALL allocate_shared_int_2D( SELEN%mesh%nTri, 3,                 SELEN%mesh%Tri,        SELEN%mesh%wTri       )
    CALL allocate_shared_int_1D( SELEN%mesh%nV,                      SELEN%mesh%nC,         SELEN%mesh%wnC        )
    CALL allocate_shared_int_2D( SELEN%mesh%nV,   SELEN%mesh%nC_mem, SELEN%mesh%C,          SELEN%mesh%wC         )
    CALL allocate_shared_int_1D( SELEN%mesh%nV,                      SELEN%mesh%niTri,      SELEN%mesh%wniTri     )
    CALL allocate_shared_int_2D( SELEN%mesh%nV,   SELEN%mesh%nC_mem, SELEN%mesh%iTri,       SELEN%mesh%wiTri      )
    CALL allocate_shared_dp_1D(  SELEN%mesh%nV,                      SELEN%mesh%lat,        SELEN%mesh%wlat       )
    CALL allocate_shared_dp_1D(  SELEN%mesh%nV,                      SELEN%mesh%lon,        SELEN%mesh%wlon       )
    CALL allocate_shared_int_1D( SELEN%mesh%nV,                      SELEN%mesh%ianc,       SELEN%mesh%wianc      )

    CALL allocate_shared_dp_1D(  SELEN%mesh%nV,                      SELEN%topo_ref,        SELEN%wtopo_ref       )
    CALL allocate_shared_dp_1D(  SELEN%mesh%nV,                      SELEN%load_ref,        SELEN%wload_ref       )
    CALL allocate_shared_int_1D( SELEN%mesh%nV,                      SELEN%of_ref,          SELEN%wof_ref         )
    CALL allocate_shared_int_2D( SELEN%mesh%nV,   2,                 SELEN%icemodel_region, SELEN%wicemodel_region)

    ! Read the mesh and global topography data from the provided NetCDF file
    IF (par%master) CALL read_SELEN_global_topo_file( SELEN)
    IF (par%master) WRITE(0,*) '  Initialised a global mesh with ', SELEN%mesh%nV, ' grid cells.'
    CALL sync

    ! Create the anchor point lists
    ! =============================

    ! Allocate memory
    CALL allocate_shared_int_0D( SELEN%mesh%nanc, SELEN%mesh%wnanc)
    IF (par%master) SELEN%mesh%nanc = SELEN%mesh%ianc( SELEN%mesh%nV)
    CALL sync
    CALL allocate_shared_int_2D( SELEN%mesh%nanc, 2, SELEN%mesh%ancplist_isl, SELEN%mesh%wancplist_isl)
    CALL allocate_shared_dp_1D(  SELEN%mesh%nanc,    SELEN%mesh%ancplist_lat, SELEN%mesh%wancplist_lat)

    IF (par%master) THEN

      i1 = 1
      i2 = 1

      DO iap = 1, SELEN%mesh%nanc

        DO WHILE (SELEN%mesh%ianc( i1) < iap)
          i1 = i1+1
        END DO
        i2 = i1
        DO WHILE (SELEN%mesh%ianc( i2) == iap)
          i2 = i2+1
          IF (i2 == SELEN%mesh%nV+1) EXIT
        END DO
        i2 = i2-1

        SELEN%mesh%ancplist_isl( iap,:) = [i1,i2]
        SELEN%mesh%ancplist_lat( iap  ) = SELEN%mesh%lat( i1)

      END DO ! DO iap = 1, SELEN%mesh%nanc

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=19)

  END SUBROUTINE initialise_SELEN_global_data
  SUBROUTINE initialise_reference_fields_regional( region)
    ! Initialise the reference fields needed for SELEN on UFEMISM's regional square grids

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_reference_fields_regional'
    INTEGER                                            :: i,j
    REAL(dp), DIMENSION(:,:  ), POINTER                ::  Hi_PD_grid_GIA,  Hb_PD_grid_GIA
    INTEGER                                            :: wHi_PD_grid_GIA, wHb_PD_grid_GIA

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) WRITE(0,*) '   Initialising SELEN reference fields for region ', region%name, '...'

    CALL allocate_shared_dp_2D(  region%grid_GIA%nx, region%grid_GIA%ny, region%SELEN%load_ref, region%SELEN%wload_ref)
    CALL allocate_shared_dp_2D(  region%grid_GIA%nx, region%grid_GIA%ny, region%SELEN%topo_ref, region%SELEN%wtopo_ref)

    ! Map PD data from the ice model grid to the GIA grid

    CALL allocate_shared_dp_2D( region%grid_GIA%nx, region%grid_GIA%ny, Hi_PD_grid_GIA, wHi_PD_grid_GIA)
    CALL allocate_shared_dp_2D( region%grid_GIA%nx, region%grid_GIA%ny, Hb_PD_grid_GIA, wHb_PD_grid_GIA)

    ! (remapping operators between mesh and grid, and viceversa, already computed
    ! during initialisation of grid, and recomputed after mesh update)
    CALL map_mesh2grid_2D( region%mesh, region%grid_GIA, region%refgeo_PD%Hi, Hi_PD_grid_GIA)
    CALL map_mesh2grid_2D( region%mesh, region%grid_GIA, region%refgeo_PD%Hb, Hb_PD_grid_GIA)

    ! Compute the PD load and topo on the regional GIA grid and save them for SELEN
    DO j = 1, region%grid_GIA%ny
    DO i = region%grid_GIA%i1, region%grid_GIA%i2

      ! Only grounded ice counts as loading
      IF (.NOT. is_floating( Hi_PD_grid_GIA( i,j), Hb_PD_grid_GIA( i,j), 0._dp)) THEN
        region%SELEN%load_ref( i,j) = Hi_PD_grid_GIA( i,j)
      ELSE
        region%SELEN%load_ref( i,j) = 0._dp
      END IF

      ! Topography
      region%SELEN%topo_ref( i,j) = Hb_PD_grid_GIA( i,j)

    END DO
    END DO
    CALL sync

    ! Clean up after yourself
    CALL deallocate_shared( wHi_PD_grid_GIA)
    CALL deallocate_shared( wHb_PD_grid_GIA)

    ! Initialise SELEN results on the regional GIA grid
    CALL allocate_shared_dp_2D( region%grid_GIA%nx, region%grid_GIA%ny, region%SELEN%dHb_t_grid_GIA,       region%SELEN%wdHb_t_grid_GIA      )
    CALL allocate_shared_dp_2D( region%grid_GIA%nx, region%grid_GIA%ny, region%SELEN%dHb_tplusdt_grid_GIA, region%SELEN%wdHb_tplusdt_grid_GIA)
    CALL allocate_shared_dp_2D( region%grid_GIA%nx, region%grid_GIA%ny, region%SELEN%SL_t_grid_GIA,        region%SELEN%wSL_t_grid_GIA       )
    CALL allocate_shared_dp_2D( region%grid_GIA%nx, region%grid_GIA%ny, region%SELEN%SL_tplusdt_grid_GIA,  region%SELEN%wSL_tplusdt_grid_GIA )

    ! Initialise SELEN results on the regional mesh
    CALL allocate_shared_dp_1D( region%mesh%nV,                         region%SELEN%dHb_t,                region%SELEN%wdHb_t               )
    CALL allocate_shared_dp_1D( region%mesh%nV,                         region%SELEN%dHb_tplusdt,          region%SELEN%wdHb_tplusdt         )
    CALL allocate_shared_dp_1D( region%mesh%nV,                         region%SELEN%SL_t,                 region%SELEN%wSL_t                )
    CALL allocate_shared_dp_1D( region%mesh%nV,                         region%SELEN%SL_tplusdt,           region%SELEN%wSL_tplusdt          )

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=10)

  END SUBROUTINE initialise_reference_fields_regional
  SUBROUTINE create_icemodel_regional_pixel_list( SELEN, region, region_label)
    ! Create a list of all SELEN global grid pixels that lie inside this ice model region

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(IN)    :: SELEN
    TYPE(type_model_region),             INTENT(INOUT) :: region
    INTEGER,                             INTENT(IN)    :: region_label

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_icemodel_regional_pixel_list'
    INTEGER                                            :: i, ir

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Create lists of pixels that lie inside of the icemodel regions
    ! Determine how many there are
    CALL allocate_shared_int_0D( region%SELEN%nr, region%SELEN%wnslr)

    IF (par%master) THEN
      region%SELEN%nr = 0
      DO i = 1, SELEN%mesh%nV
        IF (SELEN%icemodel_region( i,1) == region_label) THEN
          region%SELEN%nr = region%SELEN%nr + 1
        END IF
      END DO ! DO i = 1, SELEN%mesh%nV
    END IF ! IF (par%master) THEN
    CALL sync

    ! Determine parallelisation process domains
    CALL partition_list( region%SELEN%nr, par%i, par%n, region%SELEN%ir1, region%SELEN%ir2)

    ! Allocate shared memory
    CALL allocate_shared_int_1D( region%SELEN%nr, region%SELEN%map_isl_region2glob, region%SELEN%wmap_isl_region2glob)

    ! Let the master fill in the list
    IF (par%master) THEN
      ir = 0
      DO i = 1, SELEN%mesh%nV
        IF (SELEN%icemodel_region( i,1) == region_label) THEN
          ir = ir + 1
          region%SELEN%map_isl_region2glob( ir) = i
          SELEN%icemodel_region( i,2) = ir
        END IF
      END DO ! DO i = 1, SELEN%mesh%nV
    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=2)

  END SUBROUTINE create_icemodel_regional_pixel_list
  SUBROUTINE create_icemodel_pixel_list( SELEN)
    ! For each SELEN global grid pixel, list to which ice model region it belongs (if any)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_icemodel_pixel_list'
    INTEGER                                            :: i, ir

    ! Add routine to path
    CALL init_routine( routine_name)

    CALL allocate_shared_int_0D( SELEN%nel_icemodel, SELEN%wnel_icemodel)

    ! Let the master determine how many there are
    IF (par%master) THEN
      SELEN%nel_icemodel = 0
      DO i = 1, SELEN%mesh%nV
        IF (SELEN%icemodel_region( i,1) > 0) SELEN%nel_icemodel = SELEN%nel_icemodel + 1
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    CALL allocate_shared_int_2D( SELEN%nel_icemodel, 2, SELEN%isl_icemodel, SELEN%wisl_icemodel)

    ! Let the master fill in the data
    IF (par%master) THEN
      ir = 0
      DO i = 1, SELEN%mesh%nV
        IF (SELEN%icemodel_region( i,1) > 0) THEN
          ir = ir + 1
          SELEN%isl_icemodel( ir,:) = [i, SELEN%icemodel_region( i,1)]
        END IF
      END DO
    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected=2)

  END SUBROUTINE create_icemodel_pixel_list
  SUBROUTINE adapt_SELEN_global_data_to_ice_model( SELEN, NAM, EAS, GRL, ANT)
    ! Adapt SELEN's global reference fields with the initial data from the ice model regions

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(INOUT) :: SELEN
    TYPE(type_model_region),             INTENT(INOUT) :: NAM, EAS, GRL, ANT

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'adapt_SELEN_global_data_to_ice_model'
    REAL(dp), DIMENSION(:    ), POINTER                :: Hi_NAM, Hi_EAS, Hi_GRL, Hi_ANT ! ice model ice thickness  (mapped from regional square grid to SELEN mesh)
    REAL(dp), DIMENSION(:    ), POINTER                :: Hb_NAM, Hb_EAS, Hb_GRL, Hb_ANT ! ice model bed topography (mapped from regional square grid to SELEN mesh)
    INTEGER                                            :: wHi_NAM, wHi_EAS, wHi_GRL, wHi_ANT, wHb_NAM, wHb_EAS, wHb_GRL, wHb_ANT
    INTEGER                                            :: i

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate shared memory
    CALL allocate_shared_dp_1D( SELEN%mesh%nV, Hi_NAM, wHi_NAM)
    CALL allocate_shared_dp_1D( SELEN%mesh%nV, Hi_EAS, wHi_EAS)
    CALL allocate_shared_dp_1D( SELEN%mesh%nV, Hi_GRL, wHi_GRL)
    CALL allocate_shared_dp_1D( SELEN%mesh%nV, Hi_ANT, wHi_ANT)
    CALL allocate_shared_dp_1D( SELEN%mesh%nV, Hb_NAM, wHb_NAM)
    CALL allocate_shared_dp_1D( SELEN%mesh%nV, Hb_EAS, wHb_EAS)
    CALL allocate_shared_dp_1D( SELEN%mesh%nV, Hb_GRL, wHb_GRL)
    CALL allocate_shared_dp_1D( SELEN%mesh%nV, Hb_ANT, wHb_ANT)

  ! Map ice-model reference load from the regional GIA grids to the SELEN global grid
  ! =================================================================================

    IF (C%do_NAM) CALL map_GIA_grid_to_SELEN( SELEN, NAM, NAM%SELEN%load_ref, Hi_NAM, .TRUE.)
    IF (C%do_EAS) CALL map_GIA_grid_to_SELEN( SELEN, EAS, EAS%SELEN%load_ref, Hi_EAS, .TRUE.)
    IF (C%do_GRL) CALL map_GIA_grid_to_SELEN( SELEN, GRL, GRL%SELEN%load_ref, Hi_GRL, .TRUE.)
    IF (C%do_ANT) CALL map_GIA_grid_to_SELEN( SELEN, ANT, ANT%SELEN%load_ref, Hi_ANT, .TRUE.)

    SELEN%load_ref( C%SELEN_i1:C%SELEN_i2) = Hi_NAM( C%SELEN_i1:C%SELEN_i2) + &
                                             Hi_EAS( C%SELEN_i1:C%SELEN_i2) + &
                                             Hi_GRL( C%SELEN_i1:C%SELEN_i2) + &
                                             Hi_ANT( C%SELEN_i1:C%SELEN_i2)
    CALL sync

  ! Overwrite SELEN bed topography with data from the ice model
  ! (slightly trickier, as data outside of the ice model regions must be kept)
  ! ==========================================================================

    IF (C%do_NAM) CALL map_GIA_grid_to_SELEN( SELEN, NAM, NAM%SELEN%topo_ref, Hb_NAM, .FALSE.)
    IF (C%do_EAS) CALL map_GIA_grid_to_SELEN( SELEN, EAS, EAS%SELEN%topo_ref, Hb_EAS, .FALSE.)
    IF (C%do_GRL) CALL map_GIA_grid_to_SELEN( SELEN, GRL, GRL%SELEN%topo_ref, Hb_GRL, .FALSE.)
    IF (C%do_ANT) CALL map_GIA_grid_to_SELEN( SELEN, ANT, ANT%SELEN%topo_ref, Hb_ANT, .FALSE.)

    DO i = C%SELEN_i1, C%SELEN_i2
      IF     (SELEN%icemodel_region( i,1) == 1) THEN
        SELEN%topo_ref( i) = Hb_NAM( i)
      ELSEIF (SELEN%icemodel_region( i,1) == 2) THEN
        SELEN%topo_ref( i) = Hb_EAS( i)
      ELSEIF (SELEN%icemodel_region( i,1) == 3) THEN
        SELEN%topo_ref( i) = Hb_GRL( i)
      ELSEIF (SELEN%icemodel_region( i,1) == 4) THEN
        SELEN%topo_ref( i) = Hb_ANT( i)
      END IF
    END DO

  ! Recalculate reference ocean function on the SELEN global grid
  ! =============================================================

    DO i = C%SELEN_i1, C%SELEN_i2
      IF (is_floating( SELEN%load_ref( i), SELEN%topo_ref( i), 0._dp)) THEN
        SELEN%of_ref( i) = 1
      ELSE
        SELEN%of_ref( i) = 0
      END IF
    END DO

    ! Clean up after yourself
    CALL deallocate_shared( wHi_NAM)
    CALL deallocate_shared( wHi_EAS)
    CALL deallocate_shared( wHi_GRL)
    CALL deallocate_shared( wHi_ANT)
    CALL deallocate_shared( wHb_NAM)
    CALL deallocate_shared( wHb_EAS)
    CALL deallocate_shared( wHb_GRL)
    CALL deallocate_shared( wHb_ANT)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE adapt_SELEN_global_data_to_ice_model

  ! == Initialise spherical harmonics
  SUBROUTINE initialise_spherical_harmonics( SELEN, region, version_number)
    ! Check if spherical harmonics files for the current ice model setting exist. If so,
    ! read those. If not, create them from scratch (which is heavy, so hope they exist!)
    ! Note that the spherical harmonics are not stored in memory, but in the external
    ! files, with the paths to those files stored in region%SELEN%sh_filename

    ! When creating a set of SH files, a header file is created that describes
    ! the ice model settings for which those SH files were created. We check all
    ! existing header files, if any of them match the current settings, we read
    ! the SH files listed there. If none of them match, we create a set of SH
    ! files (and the accompanying header file) from scratch

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(IN)    :: SELEN
    TYPE(type_model_region),             INTENT(INOUT) :: region
    CHARACTER(LEN=256),                  INTENT(IN)    :: version_number

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_spherical_harmonics'
    LOGICAL                                            :: foundmatch
    INTEGER                                            :: newfolder_i

    ! Add routine to path
    CALL init_routine( routine_name)

    ! First, check if any existing header matches the current ice model set-up.
    CALL check_for_matching_sh_header( region, version_number, foundmatch, newfolder_i)

    IF (foundmatch) THEN
      IF (par%master) WRITE(0,*) '   Found valid  spherical harmonics in folder "', TRIM(region%SELEN%sh_foldername), '"'
    ELSE
      ! No header fitting the current ice model set-up was found. Create a new one describing
      ! the current set-up, and generate spherical harmonics files from scratch.
      CALL create_new_sh_folder( region, version_number, newfolder_i)
      IF (par%master) WRITE(0,*) '   Creating new spherical harmonics in folder "', TRIM(region%SELEN%sh_foldername), '"'
      CALL create_new_sh_files( SELEN, region)
    END IF ! IF (.NOT. foundmatch) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_spherical_harmonics
  SUBROUTINE create_new_sh_folder( region, version_number, newfolder_i)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    CHARACTER(LEN=256),                  INTENT(IN)    :: version_number
    INTEGER,                             INTENT(IN)    :: newfolder_i

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_new_sh_folder'
    INTEGER, DIMENSION(8)                              :: datevec
    CHARACTER(LEN=256)                                 :: header_filename

    ! Add routine to path
    CALL init_routine( routine_name)

  ! Create a new folder where the new spherical
  ! harmonics binary files will be stored.
  ! ==========================================

    IF (par%master) THEN

      IF     (newfolder_i < 10)   THEN
        WRITE( region%SELEN%sh_foldername,'(A,A,A,A,I1)') TRIM(C%selen_dir), '/spherical_harmonics_', region%name, '_00', newfolder_i
      ELSEIF (newfolder_i < 100)  THEN
        WRITE( region%SELEN%sh_foldername,'(A,A,A,A,I2)') TRIM(C%selen_dir), '/spherical_harmonics_', region%name, '_0',  newfolder_i
      ELSEIF (newfolder_i < 1000) THEN
        WRITE( region%SELEN%sh_foldername,'(A,A,A,A,I3)') TRIM(C%selen_dir), '/spherical_harmonics_', region%name, '_',   newfolder_i
      END IF

      CALL system('mkdir ' // TRIM(region%SELEN%sh_foldername))

    END IF ! IF (par%master) THEN
    CALL sync

    CALL MPI_BCAST( region%SELEN%sh_foldername, 256, MPI_CHAR,    0, MPI_COMM_WORLD, ierr)

  ! Create a header file describing the current UFEMISM set-up.
  ! ===========================================================

    IF (par%master) THEN
      header_filename = TRIM( region%SELEN%sh_foldername)//'/'//'header.txt'
      OPEN( UNIT = 1337, FILE = header_filename, STATUS = 'NEW')

      CALL date_and_time( VALUES = datevec)

      WRITE(UNIT = 1337, FMT = '(A)') '&HEADER'
      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A,I4,A,I2,A,I2)') '! Icemodel-SELEN header file, created on ', datevec(1), '-', datevec(2), '-', datevec(3)
      WRITE(UNIT = 1337, FMT = '(A)') '!'
      WRITE(UNIT = 1337, FMT = '(A)') '! This header describes the icemodel set-up that was used to created a set of'
      WRITE(UNIT = 1337, FMT = '(A)') '! spherical harmonics files. Since creating these is computationally intensive,'
      WRITE(UNIT = 1337, FMT = '(A)') '! reading them from files is preferred. These header files make ice model'
      WRITE(UNIT = 1337, FMT = '(A)') '! a bit more flexible when using different input files for the four model regions.'
      WRITE(UNIT = 1337, FMT = '(A)') ''

      WRITE(UNIT = 1337, FMT = '(A,A,A)')   'version_number       = "', TRIM(version_number), '"'
      WRITE(UNIT = 1337, FMT = '(A)') ''

      WRITE(UNIT = 1337, FMT = '(A,I5)')    'nx                    = ', region%grid_GIA%nx
      WRITE(UNIT = 1337, FMT = '(A,I5)')    'ny                    = ', region%grid_GIA%ny
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'dx                    = ', region%grid_GIA%dx
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'lambda_M              = ', region%grid_GIA%lambda_M
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'phi_M                 = ', region%grid_GIA%phi_M
      WRITE(UNIT = 1337, FMT = '(A,F14.4)') 'alpha_stereo          = ', region%grid_GIA%alpha_stereo
      WRITE(UNIT = 1337, FMT = '(A,I5)')    'SELEN_n_harmonics     = ', C%SELEN_n_harmonics

      WRITE(UNIT = 1337, FMT = '(A)') ''
      WRITE(UNIT = 1337, FMT = '(A)') '/'

      CLOSE(UNIT = 1337)

    END IF ! IF (par%master) THEN
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_new_sh_folder
  SUBROUTINE read_sh_header( header_filename, version_number, nx, ny, dx, lambda_M, phi_M, alpha_stereo, SELEN_n_harmonics)

    IMPLICIT NONE

    ! In/output variables:
    CHARACTER(LEN=256),                  INTENT(IN)    :: header_filename
    CHARACTER(LEN=256),                  INTENT(OUT)   :: version_number
    INTEGER,                             INTENT(OUT)   :: nx
    INTEGER,                             INTENT(OUT)   :: ny
    REAL(dp),                            INTENT(OUT)   :: dx
    REAL(dp),                            INTENT(OUT)   :: lambda_M
    REAL(dp),                            INTENT(OUT)   :: phi_M
    REAL(dp),                            INTENT(OUT)   :: alpha_stereo
    INTEGER,                             INTENT(OUT)   :: SELEN_n_harmonics

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'read_sh_header'
    INTEGER                                            :: ios

    ! The NAMELIST that's used to read the external header file.
    NAMELIST /HEADER/version_number,                      &
                     nx,                                  &
                     ny,                                  &
                     dx,                                  &
                     lambda_M,                            &
                     phi_M,                               &
                     alpha_stereo,                        &
                     SELEN_n_harmonics

    ! Add routine to path
    CALL init_routine( routine_name)

    OPEN( UNIT = 29, FILE = TRIM(header_filename), STATUS='OLD', ACTION='READ', iostat=ios)
    IF (ios /= 0) THEN
      CALL crash('could not open header file "' // TRIM( header_filename) // '"')
    END IF

    ! In the following statement the entire configuration file is read, using the namelist (NML=CONFIG)
    READ(  UNIT = 29, NML = HEADER, IOSTAT = ios)
    CLOSE( UNIT = 29)

    IF (ios /= 0) THEN
      CALL crash('could not read header file "' // TRIM( header_filename) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE read_sh_header
  SUBROUTINE check_for_matching_sh_header( region, version_number, foundmatch, newfolder_i)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_model_region),             INTENT(INOUT) :: region
    CHARACTER(LEN=256),                  INTENT(IN)    :: version_number
    LOGICAL,                             INTENT(OUT)   :: foundmatch
    INTEGER,                             INTENT(OUT)   :: newfolder_i

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'check_for_matching_sh_header'
    CHARACTER(LEN=256)                                 :: header_filename
    LOGICAL                                            :: header_exists

    CHARACTER(LEN=256)                                 :: version_number_read
    INTEGER                                            :: nx_read
    INTEGER                                            :: ny_read
    REAL(dp)                                           :: dx_read
    REAL(dp)                                           :: lambda_M_read
    REAL(dp)                                           :: phi_M_read
    REAL(dp)                                           :: alpha_stereo_read
    INTEGER                                            :: SELEN_n_harmonics_read

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) THEN

      foundmatch  = .FALSE.

      newfolder_i = 1
      DO WHILE (newfolder_i < 1000)

        ! Generate a foldername to inspect
        IF     (newfolder_i < 10)   THEN
          WRITE( region%SELEN%sh_foldername,'(A,A,A,A,I1)') TRIM(C%selen_dir), '/spherical_harmonics_', region%name, '_00', newfolder_i
        ELSEIF (newfolder_i < 100)  THEN
          WRITE( region%SELEN%sh_foldername,'(A,A,A,A,I2)') TRIM(C%selen_dir), '/spherical_harmonics_', region%name, '_0',  newfolder_i
        ELSEIF (newfolder_i < 1000) THEN
          WRITE( region%SELEN%sh_foldername,'(A,A,A,A,I3)') TRIM(C%selen_dir), '/spherical_harmonics_', region%name, '_',   newfolder_i
        ELSE
          CALL crash('tried a thousand folders!')
        END IF

        ! Check if a header in this folder exists. If not, then we've inspected all existing headers
        ! without finding the good one, so we must generate the SH files from scratch.
        header_filename = TRIM( region%SELEN%sh_foldername) // '/header.txt'
        INQUIRE( FILE = header_filename, EXIST = header_exists)

        IF (.NOT. header_exists) THEN
          ! No more headers exist to be inspected. Return foundmatch = .FALSE.,
          ! so that a new header file will be created at header_i, and a new
          ! spherical harmonics file will be calculated.

          EXIT

        ELSE
          ! If the header exists, read it and see if it fits the current UFEMISM set-up.

          CALL read_sh_header( header_filename, version_number_read, nx_read, ny_read, dx_read, lambda_M_read, phi_M_read, alpha_stereo_read, SELEN_n_harmonics_read)

          IF ( TRIM(version_number_read) == TRIM(version_number)         .AND. &
               nx_read                   == region%grid_GIA%nx           .AND. &
               ny_read                   == region%grid_GIA%ny           .AND. &
               dx_read                   == region%grid_GIA%dx           .AND. &
               lambda_M_read             == region%grid_GIA%lambda_M     .AND. &
               phi_M_read                == region%grid_GIA%phi_M        .AND. &
               alpha_stereo_read         == region%grid_GIA%alpha_stereo .AND. &
               SELEN_n_harmonics_read    == C%SELEN_n_harmonics) THEN
            ! This header matches the current model set-up!

            foundmatch = .TRUE.
            EXIT

          ELSE
            ! This header doesn't match the current model set-up. Try the next one.
            newfolder_i = newfolder_i + 1
          END IF

        END IF

      END DO ! DO WHILE (header_i < 1000)

    END IF ! IF (par%master) THEN
    CALL sync

    CALL MPI_BCAST( region%SELEN%sh_foldername, 256, MPI_CHAR,    0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( foundmatch,                 1,   MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( header_filename,            256, MPI_CHAR,    0, MPI_COMM_WORLD, ierr)
    CALL MPI_BCAST( newfolder_i,                1,   MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE check_for_matching_sh_header
  SUBROUTINE create_new_sh_files( SELEN, region)
    ! Calculate the spherical harmonics for this UFEMISM region's square grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_SELEN_global),             INTENT(IN)    :: SELEN
    TYPE(type_model_region),             INTENT(INOUT) :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'create_new_sh_files'
    INTEGER                                            :: i,j,li
    CHARACTER(LEN=256)                                 :: sh_filename
    COMPLEX*16, DIMENSION( C%SELEN_jmax)               :: Ypx
    INTEGER                                            :: nblocks, bi, b1, b2, ir, ir1, ir2
    REAL(dp),   DIMENSION(0:C%SELEN_n_harmonics+1)     :: LEG       ! Legendre polynome
    REAL(dp),   DIMENSION(0:C%SELEN_n_harmonics  )     :: T         ! for the shape factors

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Degree-dependent factors
    CALL PLEG( C%SELEN_n_harmonics+1, 90.- (180._dp / pi) * C%SELEN_alfa * 1., LEG)
    T(0) = (1._dp - COS(C%SELEN_alfa)) / 2._dp
    DO li = 1, C%SELEN_n_harmonics
      T( li) = (LEG( li-1) - LEG( li+1)) / (4._dp*li + 2._dp)
    END DO

  ! Calculate spherical harmonics for the GIA grid. Create a separate file
  ! for each column to keep file sizes and memory usage down.
  ! ===========================================================================

    ! Go over all the columns assigned to this process
    DO i = region%grid_GIA%i1, region%grid_GIA%i2

      ! Determine file name
      sh_filename = ''
      IF     (i < 10)    THEN
        WRITE( sh_filename,'(A,A,I1,A)') TRIM(region%SELEN%sh_foldername), '/sh_square_000', i, '.bin'
      ELSEIF (i < 100)   THEN
        WRITE( sh_filename,'(A,A,I2,A)') TRIM(region%SELEN%sh_foldername), '/sh_square_00',  i, '.bin'
      ELSEIF (i < 1000)  THEN
        WRITE( sh_filename,'(A,A,I3,A)') TRIM(region%SELEN%sh_foldername), '/sh_square_0',   i, '.bin'
      ELSEIF (i < 10000) THEN
        WRITE( sh_filename,'(A,A,I4,A)') TRIM(region%SELEN%sh_foldername), '/sh_square_',    i, '.bin'
      END IF

      ! Create a new binary file
      OPEN( UNIT = 1337+par%i, FILE = sh_filename, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')

      ! Calculate and write spherical harmonics for all rows in this column
      DO j = 1, region%grid_GIA%ny

        ! Calculate
        CALL harmo( SELEN, C%SELEN_n_harmonics, region%grid_GIA%lon( i,j), region%grid_GIA%lat( i,j), Ypx)

        ! Adjust
        DO li = 1, C%SELEN_jmax
          Ypx( li) = Ypx( li) * (2 - DOM( SELEN, li))
        END DO ! DO li = 1, C%SELEN_jmax

        ! Write
        WRITE( UNIT = 1337+par%i) Ypx

      END DO ! DO j = 1, region%grid_GIA%ny

      ! Close the binary file
      CLOSE( UNIT = 1337+par%i)

    END DO ! DO i = region%grid_GIA%i1, region%grid_GIA%i2
    CALL sync

  ! Calculate spherical harmonics for the SELEN global grid pixels inside this ice model region.
  ! Create a separate file for each 100 pixels to keep file sizes and memory usage down.
  ! ============================================================================================

    ! Determine how many blocks of 100 pixels there are, divide those over the processes
    nblocks = CEILING( REAL(region%SELEN%nr) / 100._dp)
    CALL partition_list( nblocks, par%i, par%n, b1, b2)

    ! Go over all the blocks assigned to this process
    DO bi = b1, b2

      ! Determine file name
      sh_filename = ''
      IF     (bi < 10)    THEN
        WRITE( sh_filename,'(A,A,I1,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_000', bi, '.bin'
      ELSEIF (bi < 100)   THEN
        WRITE( sh_filename,'(A,A,I2,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_00',  bi, '.bin'
      ELSEIF (bi < 1000)  THEN
        WRITE( sh_filename,'(A,A,I3,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_0',   bi, '.bin'
      ELSEIF (bi < 10000) THEN
        WRITE( sh_filename,'(A,A,I4,A)') TRIM(region%SELEN%sh_foldername), '/sh_glob_',    bi, '.bin'
      END IF

      ! Create a new binary file
      OPEN( UNIT = 1337+par%i, FILE = sh_filename, STATUS = 'UNKNOWN', FORM = 'UNFORMATTED')

      ! Calculate and write spherical harmonics for all pixels in this block
      ir1 = (bi-1)*100+1
      ir2 = MIN( bi*100, region%SELEN%nr)
      DO ir = ir1, ir2

        ! Look up this pixel's global index
        i = region%SELEN%map_isl_region2glob( ir)

        CALL harmo( SELEN, C%SELEN_n_harmonics, SELEN%mesh%lon( i), SELEN%mesh%lat( i), Ypx)
        DO li = 1, C%SELEN_jmax
          Ypx( li) = T( SELEN%LJ_VAL( li)) * CONJG( Ypx( li))
        END DO

        ! Write
        WRITE( UNIT = 1337+par%i) Ypx

      END DO ! DO isli = isli1, isli2

      ! Close the binary file
      CLOSE( UNIT = 1337+par%i)

    END DO ! DO bi = b1, b2
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_new_sh_files

END MODULE SELEN_main_module