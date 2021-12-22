MODULE isotopes_module

  ! Contains all the routines for calculating the isotope content of the ice sheet.

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
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D
  USE netcdf_module,                   ONLY: debug, write_to_debug_file
  
  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh, type_model_region, type_remapping
  USE mesh_mapping_module,             ONLY: remap_field_dp

  IMPLICIT NONE
    
CONTAINS

  ! Run the isotopes model
  SUBROUTINE run_isotopes_model( region)
  
    ! Based on the ANICE routines by Bas de Boer (November 2010).
    !
    ! using the implicit ice fluxes (vertical averaged) and applied mass balance fields from the 
    ! ice thickness subroutine to calculate the advected isotope flux from all 4 directions.
    !
    ! Because of the dynamic shelf, we can also calculate IsoIce over the shelf area and threat 
    ! all ice covered gridpoints the same since we are using vertically averaged velocities
    !
    ! IsoIce_adv should be zero for no ice (or at least have a small value, just like Hi)
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region
    
    ! Local variables
    INTEGER                                            :: vi, ci, vj, aci
    REAL(dp)                                           :: Ts, Ts_ref, Hs, Hs_ref
    REAL(dp)                                           :: IsoMin,IsoMax    ! minimum and maximum value of IsoIce
    
    REAL(dp)                                           :: Upar, dVIso
    REAL(dp)                                           ::  dIso_dt,  VIso_old,  VIso_new
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in run_isotopes_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Calculate the isotope content of annual mean precipitation
    ! ==========================================================
    
    IsoMax = -1E9_dp
    IsoMin =  1E9_dp
    
    DO vi = region%mesh%vi1, region%mesh%vi2
      
      Ts     = SUM( region%climate%applied%T2m( vi,:)) / 12._dp
      Ts_ref = SUM( region%climate%PD_obs%T2m(  vi,:)) / 12._dp
      Hs     = region%ice%Hs_a( vi)
      Hs_ref = region%climate%PD_obs%Hs( vi)
      
      region%ice%IsoSurf( vi) = region%ice%IsoRef( vi)                 &
                              + 0.35_dp              * (Ts - Ts_ref    &
                              - C%constant_lapserate * (Hs - Hs_ref))  &
                              - 0.0062_dp            * (Hs - Hs_ref)   ! from Clarke et al., 2005  
        
      IF (region%ice%mask_ice_a( vi) == 1) THEN
        IsoMax = MAX( IsoMax, region%ice%IsoSurf( vi))
        IsoMin = MIN( IsoMin, region%ice%IsoSurf( vi))
      END IF           

    END DO
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, IsoMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, IsoMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! Calculate the mass gain/loss of d18O
    DO vi = region%mesh%vi1, region%mesh%vi2
    
      region%ice%MB_iso( vi) = 0._dp
    
      IF ( region%SMB%SMB_year( vi) > 0._dp ) THEN
        ! Surface accumulation has the isotope content of precipitation
        region%ice%MB_iso( vi) = region%ice%MB_iso( vi) + region%SMB%SMB_year( vi) * region%ice%IsoSurf( vi) * region%mesh%A( vi)    ! (applied MB) mass gain, so d18O from precipitation
      ELSE
        ! Surface melt has the isotope content of the ice itself
        region%ice%MB_iso( vi) = region%ice%MB_iso( vi) + region%SMB%SMB_year( vi) * region%ice%IsoIce(  vi) * region%mesh%A( vi)    ! (applied MB) mass loss, so d18O from ice
      END IF
      
      ! Both basal melt and basal freezing have the isotope content of the ice itself (the latter
      ! is not really true, but it's the best we can do for now)
      region%ice%MB_iso(   vi) = region%ice%MB_iso( vi) + region%BMB%BMB(      vi) * region%ice%IsoIce(  vi) * region%mesh%A( vi)
      
    END DO
    CALL sync

    ! Calculate the new d18O_ice from the ice fluxes and applied mass balance
    ! =======================================================================
    
    IF (par%master) WRITE(0,*) 'run_isotopes_model - FIXME!'
    CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)

!    DO vi = region%mesh%vi1, region%mesh%vi2
!    
!      region%ice%IsoIce_new( vi) = region%ice%IsoIce( vi)
!    
!      ! Don't treat boundary vertices
!      IF (region%mesh%edge_index( vi) > 0) CYCLE
!    
!      ! Only treat ice-covered vertices
!      IF (region%ice%mask_ice_a( vi) == 1) THEN
!      
!        dIso_dt = 0._dp
!      
!        ! Calculate advective isotope fluxes
!        DO ci = 1, region%mesh%nC( vi)
!        
!          vj  = region%mesh%C(    vi,ci)
!          aci = region%mesh%iAci( vi,ci)
!          
!          IF (region%mesh%Aci( aci,1) == vi) THEN
!            ! Velocity defined from vi to vj
!            Upar =  (region%ice%Up_SIA_Ac( aci) + region%ice%Up_SSA_Ac( aci))
!          ELSE
!            ! Velocity defined from vj to vi
!            Upar = -(region%ice%Up_SIA_Ac( aci) + region%ice%Up_SSA_Ac( aci))
!          END IF
!          
!          IF (Upar > 0._dp) THEN
!            ! Ice flows from vi to vj
!            dVIso = Upar * region%ice%IsoIce( vi) * region%ice%Hi( vi) * region%mesh%Cw( vi,ci)
!          ELSE
!            ! Ice flows from vj to vi
!            dVIso = Upar * region%ice%IsoIce( vj) * region%ice%Hi( vj) * region%mesh%Cw( vi,ci)
!          END IF
!          
!          dIso_dt = dIso_dt + dViso
!          
!        END DO ! DO ci = 1, region%mesh%nC( vi)
!        
!        ! Add surface + basal mass balance
!        dIso_dt = dIso_dt + region%ice%MB_iso( vi)
!        
!        ! Update vertically averaged ice isotope content
!        VIso_old = region%ice%IsoIce( vi) * region%ice%Hi_prev( vi) * region%mesh%A( vi)
!        VIso_new = VIso_old + dIso_dt * region%dt
!        region%ice%IsoIce_new( vi) = MIN( IsoMax, MAX( IsoMin, VIso_new / (region%mesh%A( vi) * region%ice%Hi( vi)) ))
!
!      END IF ! IF (ice%mask_ice( vi) == 1) THEN
!      
!    END DO
!    CALL sync
    
    ! Update field
    region%ice%IsoIce( region%mesh%vi1:region%mesh%vi2) = region%ice%IsoIce_new( region%mesh%vi1:region%mesh%vi2)
    
    ! Calculate mean isotope content of the whole ice sheet
    CALL calculate_isotope_content( region%mesh, region%ice%Hi_a, region%ice%IsoIce, region%mean_isotope_content, region%d18O_contribution)
    
  END SUBROUTINE run_isotopes_model
  SUBROUTINE calculate_isotope_content( mesh, Hi, IsoIce, mean_isotope_content, d18O_contribution)
    ! Calculate mean isotope content of the whole ice sheet
    
    USE parameters_module, ONLY: ice_density, seawater_density, ocean_area, mean_ocean_depth
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: Hi
    REAL(dp), DIMENSION(:    ),          INTENT(IN)    :: IsoIce
    REAL(dp),                            INTENT(OUT)   :: mean_isotope_content
    REAL(dp),                            INTENT(OUT)   :: d18O_contribution
    
    ! Local variables
    INTEGER                                            :: vi
    REAL(dp)                                           :: Hi_msle
    REAL(dp)                                           :: total_isotope_content
    REAL(dp)                                           :: total_ice_volume_msle

    ! Calculate total isotope content
    ! ===============================
    
    total_isotope_content = 0._dp
    total_ice_volume_msle = 0._dp
    
    DO vi = mesh%vi1, mesh%vi2
          
      IF (Hi( vi) > 0._dp) THEN
        Hi_msle = Hi( vi) * mesh%A( vi) * ice_density / (seawater_density * ocean_area)
        total_isotope_content = total_isotope_content + Hi_msle * IsoIce( vi)
        total_ice_volume_msle = total_ice_volume_msle + Hi_msle
      END IF   

    END DO
    CALL sync
    
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_isotope_content, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    CALL MPI_ALLREDUCE( MPI_IN_PLACE, total_ice_volume_msle, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    
    ! Weighted average of isotope content with Hi
    IF (par%master) THEN
      IF (total_ice_volume_msle > 0._dp) THEN
        mean_isotope_content  = total_isotope_content / total_ice_volume_msle
      ELSE
        mean_isotope_content  = 0._dp
      END IF
    END IF
    CALL sync
    
    ! Contribution to benthic d18O
    d18O_contribution = -1._dp * mean_isotope_content * total_ice_volume_msle / mean_ocean_depth
    
  END SUBROUTINE calculate_isotope_content
  
  ! Initialise the isotopes model (allocating shared memory)
  SUBROUTINE initialise_isotopes_model( region)
    ! Allocate memory for the data fields of the isotopes model and initialise the reference fields.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region
    
    ! Local variables
    CHARACTER(LEN=64), PARAMETER                       :: routine_name = 'initialise_isotopes_model'
    INTEGER                                            :: n1, n2
    INTEGER                                            :: vi
    REAL(dp)                                           :: Ts, Ts_ref, Hs, Hs_ref
    
    n1 = par%mem%n
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in initialise_isotopes_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    IF (par%master) WRITE (0,*) '  Initialising isotopes model...'
    
    ! Allocate memory
    CALL allocate_shared_dp_1D( region%mesh%nV, region%ice%IsoRef    , region%ice%wIsoRef    )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%ice%IsoSurf   , region%ice%wIsoSurf   )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%ice%MB_iso    , region%ice%wMB_iso    )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%ice%IsoIce    , region%ice%wIsoIce    )
    CALL allocate_shared_dp_1D( region%mesh%nV, region%ice%IsoIce_new, region%ice%wIsoIce_new)
    
    ! Calculate present-day isotope content of precipitation    
    CALL calculate_reference_isotopes( region)
    
    ! Initialise ice sheet isotope content with the isotope content of present-day annual mean precipitation,
    ! so that we can calculate the total present-day isotope content of the ice-sheet
    DO vi = region%mesh%vi1, region%mesh%vi2
    
      IF (region%refgeo_PD%Hi( vi) > 0._dp) THEN
      
        Ts     = SUM( region%climate%PD_obs%T2m( vi,:)) / 12._dp
        Ts_ref = SUM( region%climate%PD_obs%T2m( vi,:)) / 12._dp
        Hs     = region%refgeo_PD%Hs( vi)
        Hs_ref = region%climate%PD_obs%Hs( vi)
        
        region%ice%IsoIce( vi) = region%ice%IsoRef( vi)                 &
                               + 0.35_dp              * (Ts - Ts_ref    &
                               - C%constant_lapserate * (Hs - Hs_ref))  &
                               - 0.0062_dp            * (Hs - Hs_ref)   ! from Clarke et al., 2005  
        
      ELSE
        region%ice%IsoIce( vi) = 0._dp ! = No ice
      END IF           

    END DO
    CALL sync
    
    ! Calculate mean isotope content of the whole ice sheet at present-day
    CALL calculate_isotope_content( region%mesh, region%refgeo_PD%Hi, region%ice%IsoIce, region%mean_isotope_content_PD, region%d18O_contribution_PD)
    
    ! Initialise ice sheet isotope content with the isotope content of annual mean precipitation at the start of the simulation
    ! (need not be the same as present-day conditions, that's why we need to repeat the calculation)
    DO vi = region%mesh%vi1, region%mesh%vi2
    
      IF (region%ice%mask_ice_a( vi) == 1) THEN
      
        Ts     = SUM( region%climate%applied%T2m( vi,:)) / 12._dp
        Ts_ref = SUM( region%climate%PD_obs%T2m(  vi,:)) / 12._dp
        Hs     = region%ice%Hs_a( vi)
        Hs_ref = region%climate%PD_obs%Hs( vi)
        
        region%ice%IsoIce( vi) = region%ice%IsoRef( vi)                 &
                               + 0.35_dp              * (Ts - Ts_ref    &
                               - C%constant_lapserate * (Hs - Hs_ref))  &
                               - 0.0062_dp            * (Hs - Hs_ref)   ! from Clarke et al., 2005  
        
      ELSE
        region%ice%IsoIce( vi) = 0._dp ! = No ice
      END IF           

    END DO
    CALL sync
    
    ! Calculate mean isotope content of the whole ice sheet at the start of the simulation
    CALL calculate_isotope_content( region%mesh, region%ice%Hi_a, region%ice%IsoIce, region%mean_isotope_content, region%d18O_contribution)
    
    n2 = par%mem%n
    CALL write_to_memory_log( routine_name, n1, n2)
    
  END SUBROUTINE initialise_isotopes_model
  SUBROUTINE calculate_reference_isotopes( region)
    ! Allocate memory for the data fields of the isotopes model.
    
    IMPLICIT NONE
    
    ! In/output variables
    TYPE(type_model_region),             INTENT(INOUT) :: region
    
    ! Local variables
    INTEGER                                            :: vi
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in calculate_reference_isotopes!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Calculate reference field of d18O of precipitation
    ! (Zwally, H. J. and Giovinetto, M. B.: Areal distribution of the oxygen-isotope ratio in Greenland, Annals of Glaciology 25, 208-213, 1997)
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%ice%IsoRef( vi) = 0.691_dp * SUM(region%climate%PD_obs%T2m( vi,:) / 12._dp) - 202.172_dp
    END DO
    CALL sync
    
  END SUBROUTINE calculate_reference_isotopes
  SUBROUTINE remap_isotopes_model( mesh_old, mesh_new, map, region)
    ! Remap or reallocate all the data fields
  
    ! In/output variables:
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                     INTENT(IN)    :: mesh_new
    TYPE(type_remapping),                INTENT(IN)    :: map
    TYPE(type_model_region),             INTENT(INOUT) :: region
    
    ! Not needed for benchmark experiments
    IF (C%do_benchmark_experiment) THEN
      IF (C%choice_benchmark_experiment == 'EISMINT_1'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_2'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_3'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_4'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_5'  .OR. &
          C%choice_benchmark_experiment == 'EISMINT_6'  .OR. &
          C%choice_benchmark_experiment == 'Halfar'     .OR. &
          C%choice_benchmark_experiment == 'Bueler'     .OR. &
          C%choice_benchmark_experiment == 'MISMIP_mod' .OR. &
          C%choice_benchmark_experiment == 'mesh_generation_test' .OR. &
          C%choice_benchmark_experiment == 'SSA_icestream') THEN
        RETURN
      ELSE 
        WRITE(0,*) '  ERROR: benchmark experiment "', TRIM(C%choice_benchmark_experiment), '" not implemented in remap_isotopes_model!'
        CALL MPI_ABORT( MPI_COMM_WORLD, cerr, ierr)
      END IF
    END IF ! IF (C%do_benchmark_experiment) THEN
    
    ! Reallocate memory
    CALL reallocate_shared_dp_1D( mesh_new%nV, region%ice%IsoRef,     region%ice%wIsoRef    )
    CALL reallocate_shared_dp_1D( mesh_new%nV, region%ice%IsoSurf,    region%ice%wIsoSurf   )
    CALL reallocate_shared_dp_1D( mesh_new%nV, region%ice%MB_iso,     region%ice%wMB_iso    )
    CALL reallocate_shared_dp_1D( mesh_new%nV, region%ice%IsoIce_new, region%ice%wIsoIce_new)
    
    ! Remap the actual isotope content of the ice sheet
    CALL remap_field_dp( mesh_old, mesh_new, map, region%ice%IsoIce, region%ice%wIsoIce, 'cons_1st_order')
    
  END SUBROUTINE remap_isotopes_model

END MODULE isotopes_module