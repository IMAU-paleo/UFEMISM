  SUBROUTINE initialise_lowres_SSA_multigrid_ice_model(ice, mesh)
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_ice_model),                INTENT(INOUT) :: ice
    TYPE(type_mesh),                     INTENT(IN)    :: mesh  
      
    ! Allocate memory
    ! ===============
    
!    ! Basic data - ice thickness, bedrock height, surface height, mask, 3D ice velocities and temperature
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Hi,                  ice%wHi                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Hb,                  ice%wHb                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Hs,                  ice%wHs                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%sealevel,            ice%wsealevel           )
!    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Us,                  ice%wUs                 )
!    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Vs,                  ice%wVs                 )
!    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%Ti,                  ice%wTi                 )  
    
    ! Different masks
!    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask,                ice%wmask               )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_ice,            ice%wmask_ice           )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_sheet,          ice%wmask_sheet         )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_shelf,          ice%wmask_shelf         )
    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_ocean,          ice%wmask_ocean         )
!    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_margin,         ice%wmask_margin        )
!    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_groundingline,  ice%wmask_groundingline )
!    CALL allocate_shared_memory_int_1D( mesh%nV,       ice%mask_calvingfront,   ice%wmask_calvingfront  )
    
!    ! Ice physical properties
!    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%Ti_pmp,              ice%wTi_pmp             )
!    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%A_flow,              ice%wA_flow             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%A_flow_mean,         ice%wA_flow_mean        )
!    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%Cpi,                 ice%wCpi                ) 
!    CALL allocate_shared_memory_dp_2D(  mesh%nV, C%nZ, ice%Ki,                  ice%wKi                 ) 
!    
!    ! Spatial derivatives and curvatures
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHi_dx,              ice%wdHi_dx             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHi_dy,              ice%wdHi_dy             )
!    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHi_dt,              ice%wdHi_dt             )
!    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHb_dx,              ice%wdHb_dx             )
!    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHb_dy,              ice%wdHb_dy             )
!    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHb_dt,              ice%wdHb_dt             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dx,              ice%wdHs_dx             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dy,              ice%wdHs_dy             )
!    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dt,              ice%wdHs_dt             )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dx_shelf,        ice%wdHs_dx_shelf       )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dHs_dy_shelf,        ice%wdHs_dy_shelf       )
    
!    ! Ice dynamics - SSA
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%U_SSA,               ice%wU_SSA              )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%V_SSA,               ice%wV_SSA              )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%tau_yield,           ice%wtau_yield          )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dUdx,                ice%wdUdx               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dUdy,                ice%wdUdy               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dVdx,                ice%wdVdx               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%dVdy,                ice%wdVdy               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%nu,                  ice%wnu                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%beta_base,           ice%wbeta_base          )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%rhs_x,               ice%wrhs_x              )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%rhs_y,               ice%wrhs_y              )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%eu_i,                ice%weu_i               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%ev_i,                ice%wev_i               )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Rx,                  ice%wRx                 )
    CALL allocate_shared_memory_dp_1D(  mesh%nV,       ice%Ry,                  ice%wRy                 )
    
    
  END SUBROUTINE initialise_lowres_SSA_multigrid_ice_model
    
  ! SSA solver - multigrid
  SUBROUTINE multigrid_SSA_solver(region)
    ! Use a sequential multigrid SOR scheme to solver the SSA really quickly.
  
    ! Input variables
    TYPE(type_model_region),         INTENT(INOUT)     :: region
    
    ! Local data    
    INTEGER                                            :: vi, ci, vj
    REAL(dp)                                           :: Ux, Uy, U
    LOGICAL, PARAMETER                                 :: debug = .FALSE.
    REAL(dp)                                           :: t1
    
  ! Upscale the input data from the final mesh to the relevant nested meshes
  ! ========================================================================
    
    ! First from the model mesh to the finest nested mesh
    IF     (region%n_nested_meshes == 7) THEN
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from model mesh to nested mesh 7'
      IF (par%master) t1 = MPI_WTIME()
      CALL map_input_data_between_nested_meshes( region%mesh, region%mesh_low_7, region%ice, region%ice_low_7, &
           region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    ELSEIF (region%n_nested_meshes == 6) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from model mesh to nested mesh 6'
      CALL map_input_data_between_nested_meshes( region%mesh, region%mesh_low_6, region%ice, region%ice_low_6, &
           region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    ELSEIF (region%n_nested_meshes == 5) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from model mesh to nested mesh 5'
      CALL map_input_data_between_nested_meshes( region%mesh, region%mesh_low_5, region%ice, region%ice_low_5, &
           region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    ELSEIF (region%n_nested_meshes == 4) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from model mesh to nested mesh 4'
      CALL map_input_data_between_nested_meshes( region%mesh, region%mesh_low_4, region%ice, region%ice_low_4, &
           region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    ELSEIF (region%n_nested_meshes == 3) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from model mesh to nested mesh 3'
      CALL map_input_data_between_nested_meshes( region%mesh, region%mesh_low_3, region%ice, region%ice_low_3, &
           region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    ELSEIF (region%n_nested_meshes == 2) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from model mesh to nested mesh 2'
      CALL map_input_data_between_nested_meshes( region%mesh, region%mesh_low_2, region%ice, region%ice_low_2, &
           region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    ELSEIF (region%n_nested_meshes == 1) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from model mesh to nested mesh 1'
      CALL map_input_data_between_nested_meshes( region%mesh, region%mesh_low_1, region%ice, region%ice_low_1, &
           region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    END IF
    
    ! Then between all the nested meshes
    IF     (region%n_nested_meshes >= 7) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from nested mesh 7 to nested mesh 6'
      CALL map_input_data_between_nested_meshes( region%mesh_low_7, region%mesh_low_6, region%ice_low_7, region%ice_low_6, &
           region%map_76_nV, region%map_76_vi, region%map_76_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    END IF
    IF     (region%n_nested_meshes >= 6) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from nested mesh 6 to nested mesh 5'
      CALL map_input_data_between_nested_meshes( region%mesh_low_6, region%mesh_low_5, region%ice_low_6, region%ice_low_5, &
           region%map_65_nV, region%map_65_vi, region%map_65_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    END IF
    IF     (region%n_nested_meshes >= 5) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from nested mesh 5 to nested mesh 4'
      CALL map_input_data_between_nested_meshes( region%mesh_low_5, region%mesh_low_4, region%ice_low_5, region%ice_low_4, &
           region%map_54_nV, region%map_54_vi, region%map_54_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    END IF
    IF     (region%n_nested_meshes >= 4) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from nested mesh 4 to nested mesh 3'
      CALL map_input_data_between_nested_meshes( region%mesh_low_4, region%mesh_low_3, region%ice_low_4, region%ice_low_3, &
           region%map_43_nV, region%map_43_vi, region%map_43_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    END IF
    IF     (region%n_nested_meshes >= 3) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from nested mesh 3 to nested mesh 2'
      CALL map_input_data_between_nested_meshes( region%mesh_low_3, region%mesh_low_2, region%ice_low_3, region%ice_low_2, &
           region%map_32_nV, region%map_32_vi, region%map_32_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    END IF
    IF     (region%n_nested_meshes >= 2) THEN
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping input data from nested mesh 2 to nested mesh 1'
      CALL map_input_data_between_nested_meshes( region%mesh_low_2, region%mesh_low_1, region%ice_low_2, region%ice_low_1, &
           region%map_21_nV, region%map_21_vi, region%map_21_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    END IF
    
! Solve the SSA on the nested meshes
! ==================================

    ! Nested mesh #1
    ! ==============
    
    IF (region%n_nested_meshes >= 1) THEN
    
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: solving SSA on nested mesh 1'
      CALL calculate_ice_velocities_SSA( region%ice_low_1, region%mesh_low_1)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      
      ! Either map to nested mesh #2 or to the model mesh
      IF (region%n_nested_meshes >= 2) THEN
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 1 to nested mesh 2'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_1, region%mesh_low_2, region%ice_low_1, region%ice_low_2, &
             region%map_12_nV, region%map_12_vi, region%map_12_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      ELSE
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 1 to model mesh'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_1, region%mesh, region%ice_low_1, region%ice, &
             region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      END IF
    
    END IF

    ! Nested mesh #2
    ! ==============
    
    IF (region%n_nested_meshes >= 2) THEN
      
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: solving SSA on nested mesh 2'
      CALL calculate_ice_velocities_SSA( region%ice_low_2, region%mesh_low_2)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      
      ! Either map to nested mesh #3 or to the model mesh
      IF (region%n_nested_meshes >= 3) THEN
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 2 to nested mesh 3'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_2, region%mesh_low_3, region%ice_low_2, region%ice_low_3, &
             region%map_23_nV, region%map_23_vi, region%map_23_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      ELSE
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 2 to model mesh'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_2, region%mesh, region%ice_low_2, region%ice, &
             region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      END IF
      
    END IF
  
    ! Nested mesh #3
    ! ==============
    
    IF (region%n_nested_meshes >= 3) THEN
    
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: solving SSA on nested mesh 3'
      CALL calculate_ice_velocities_SSA( region%ice_low_3, region%mesh_low_3)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      
      ! Either map to nested mesh #4 or to the model mesh
      IF (region%n_nested_meshes >= 4) THEN
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 3 to nested mesh 4'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_3, region%mesh_low_4, region%ice_low_3, region%ice_low_4, &
             region%map_34_nV, region%map_34_vi, region%map_34_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      ELSE
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 3 to model mesh'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_3, region%mesh, region%ice_low_3, region%ice, &
             region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      END IF
      
    END IF
  
    ! Nested mesh #4
    ! ==============
    
    IF (region%n_nested_meshes >= 4) THEN
    
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: solving SSA on nested mesh 4'
      CALL calculate_ice_velocities_SSA( region%ice_low_4, region%mesh_low_4)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      
      ! Either map to nested mesh #5 or to the model mesh
      IF (region%n_nested_meshes >= 5) THEN
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 4 to nested mesh 5'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_4, region%mesh_low_5, region%ice_low_4, region%ice_low_5, &
             region%map_45_nV, region%map_45_vi, region%map_45_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      ELSE
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 4 to model mesh'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_4, region%mesh, region%ice_low_4, region%ice, &
             region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      END IF
      
    END IF
  
    ! Nested mesh #5
    ! ==============
    
    IF (region%n_nested_meshes >= 5) THEN
    
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: solving SSA on nested mesh 5'
      CALL calculate_ice_velocities_SSA( region%ice_low_5, region%mesh_low_5)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      
      ! Either map to nested mesh #6 or to the model mesh
      IF (region%n_nested_meshes >= 6) THEN
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 5 to nested mesh 6'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_5, region%mesh_low_6, region%ice_low_5, region%ice_low_6, &
             region%map_56_nV, region%map_56_vi, region%map_56_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      ELSE
        IF (par%master) t1 = MPI_WTIME()
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 5 to model mesh'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_5, region%mesh, region%ice_low_5, region%ice, &
             region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      END IF
    
    END IF
  
    ! Nested mesh #6
    ! ==============
    
    IF (region%n_nested_meshes >= 6) THEN
    
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: solving SSA on nested mesh 6'
      CALL calculate_ice_velocities_SSA( region%ice_low_6, region%mesh_low_6)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      
      ! Either map to nested mesh #7 or to the model mesh
      IF (region%n_nested_meshes >= 7) THEN
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 6 to nested mesh 7'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_6, region%mesh_low_7, region%ice_low_6, region%ice_low_7, &
             region%map_67_nV, region%map_67_vi, region%map_67_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      ELSE
        IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 6 to model mesh'
        CALL map_SSA_solution_between_nested_meshes( region%mesh_low_6, region%mesh, region%ice_low_6, region%ice, &
             region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w)
        IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      END IF
      
    END IF
  
    ! Nested mesh #7
    ! ==============
    
    IF (region%n_nested_meshes >= 7) THEN
    
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: solving SSA on nested mesh 7'
      CALL calculate_ice_velocities_SSA( region%ice_low_7, region%mesh_low_7)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
      
      ! Map to the model mesh
      IF (par%master) t1 = MPI_WTIME()
      IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: mapping solution from nested mesh 7 to model mesh'
      CALL map_SSA_solution_between_nested_meshes( region%mesh_low_7, region%mesh, region%ice_low_7, region%ice, &
             region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w)
      IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
           
    END IF
    
  ! Finally, solve the SSA on the model mesh
  ! ===============================
    
    IF (par%master) t1 = MPI_WTIME()
    IF (par%master .AND. debug) WRITE(0,*) ' Multigrid SSA solver: solving SSA on model mesh'
    CALL calculate_ice_velocities_SSA( region%ice, region%mesh)
    IF (par%master .AND. debug) WRITE(0,'(A,F7.4,A)') '  Finished in ', MPI_WTIME() - t1, ' seconds.'
    
    ! Map data to Ac mesh for ice thickness calculation
    CALL MapAaToAc(region%mesh, region%ice%U_SSA, region%ice%Ux_SSA_ac)
    CALL MapAaToAc(region%mesh, region%ice%V_SSA, region%ice%Uy_SSA_ac)
    
    ! Get parallel and orthogonal components of SSA velocities for ice thickness calculation
    DO ci = region%mesh%c1, region%mesh%c2
      vi = region%mesh%Aci(ci,1)
      vj = region%mesh%Aci(ci,2)
      Ux = region%mesh%V(vj,1)-region%mesh%V(vi,1)
      Uy = region%mesh%V(vj,2)-region%mesh%V(vi,2)
      U  = SQRT(Ux**2+Uy**2)
      
      region%ice%Up_SSA_ac(ci) = region%ice%Ux_SSA_ac(ci) * Ux/U + region%ice%Uy_SSA_ac(ci) * Uy/U
      region%ice%Uo_SSA_ac(ci) = region%ice%Uy_SSA_ac(ci) * Ux/U - region%ice%Ux_SSA_ac(ci) * Uy/U
    END DO
    
    
  END SUBROUTINE multigrid_SSA_solver
  SUBROUTINE map_input_data_between_nested_meshes(mesh1, mesh2, ice1, ice2, map_nV, map_vi, map_w)
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh1
    TYPE(type_mesh),                     INTENT(IN)    :: mesh2
    TYPE(type_ice_model),                INTENT(INOUT) :: ice1
    TYPE(type_ice_model),                INTENT(INOUT) :: ice2
    INTEGER,  DIMENSION(:  ),            INTENT(IN)    :: map_nV
    INTEGER,  DIMENSION(:,:),            INTENT(IN)    :: map_vi
    REAL(dp), DIMENSION(:,:),            INTENT(IN)    :: map_w
  
    ! Local variables
    INTEGER                                       :: vi
    
    CALL MapMesh2Mesh( mesh2, map_nV, map_vi, map_w, ice1%Hi,           ice2%Hi          )
    CALL MapMesh2Mesh( mesh2, map_nV, map_vi, map_w, ice1%Hb,           ice2%Hb          )
    CALL MapMesh2Mesh( mesh2, map_nV, map_vi, map_w, ice1%Hs,           ice2%Hs          )
    CALL MapMesh2Mesh( mesh2, map_nV, map_vi, map_w, ice1%sealevel,     ice2%sealevel    )
    CALL MapMesh2Mesh( mesh2, map_nV, map_vi, map_w, ice1%A_flow_mean,  ice2%A_flow_mean ) 
  
    ! Update masks
    ! ============
    
    ice2%mask_ice(          mesh2%v1:mesh2%v2) = 0
    ice2%mask_sheet(        mesh2%v1:mesh2%v2) = 0
    ice2%mask_shelf(        mesh2%v1:mesh2%v2) = 0
    ice2%mask_ocean(        mesh2%v1:mesh2%v2) = 0
    CALL sync
  
    DO vi = mesh2%v1, mesh2%v2
    
      ! Determine ice
      IF (ice2%Hi(vi) > 0._dp) THEN
        ice2%mask_ice(vi) = 1
      END IF
      
      ! Determine sheet
      IF (ice2%mask_ice(vi) == 1 .AND. (ice2%Hi(vi) > (ice2%sealevel(vi) - ice2%Hb(vi)) * seawater_density/ice_density)) THEN
        ice2%mask_sheet(vi) = 1
      END IF
    
      ! Determine shelf
      IF (ice2%mask_ice(vi) == 1 .AND. ice2%mask_sheet(vi) == 0) THEN
        ice2%mask_shelf(vi) = 1
      END IF
      
      ! Determine open ocean
      IF ((ice2%Hb(vi) < ice2%sealevel(vi)) .AND. (ice2%Hi(vi) == 0._dp)) THEN
        ice2%mask_ocean(vi) = 1
      END IF
      
    END DO
    CALL sync  
    
    ! Calculate surface slopes
    CALL GetMeshDerivatives( mesh2, ice2%Hi, ice2%dHi_dx,    ice2%dHi_dy)
    CALL GetMeshDerivatives( mesh2, ice2%Hs, ice2%dHs_dx,    ice2%dHs_dy)
    
    ! Use a different surface slope for shelves. Don't really know why, but when you use the "regular" surface slope,
    ! the SSA gives ice velocities upwards of 1000 km/y. ANICE does the same thing, the comments there says they (who?)
    ! suspect it's an initialisation thing. Well, that's definitely not it since we initialise with observations rather
    ! than with output from another ice model. Anyway, it works, so leave it in...
    DO vi = mesh2%v1, mesh2%v2
      IF (ice2%mask_sheet(vi)==1) THEN
        ice2%dHs_dx_shelf(vi) = ice2%dHs_dx(vi)
        ice2%dHs_dy_shelf(vi) = ice2%dHs_dy(vi)
      ELSE
        ice2%dHs_dx_shelf(vi) = (1._dp - ice_density / seawater_density) * ice2%dHi_dx(vi)
        ice2%dHs_dy_shelf(vi) = (1._dp - ice_density / seawater_density) * ice2%dHi_dy(vi)
      END IF
    END DO
    
    
  END SUBROUTINE map_input_data_between_nested_meshes
  SUBROUTINE map_SSA_solution_between_nested_meshes(mesh1, mesh2, ice1, ice2, map_nV, map_vi, map_w)
      
    IMPLICIT NONE
    
    ! In- and output variables
    TYPE(type_mesh),                     INTENT(IN)    :: mesh1
    TYPE(type_mesh),                     INTENT(IN)    :: mesh2
    TYPE(type_ice_model),                INTENT(INOUT) :: ice1
    TYPE(type_ice_model),                INTENT(INOUT) :: ice2
    INTEGER,  DIMENSION(:  ),            INTENT(IN)    :: map_nV
    INTEGER,  DIMENSION(:,:),            INTENT(IN)    :: map_vi
    REAL(dp), DIMENSION(:,:),            INTENT(IN)    :: map_w
    
    ! Local variables
    REAL(dp), DIMENSION(:  ), POINTER                  ::  dU_SSA_1,  dV_SSA_1
    REAL(dp), DIMENSION(:  ), POINTER                  ::  dU_SSA_2,  dV_SSA_2
    INTEGER                                            :: wdU_SSA_1, wdV_SSA_1
    INTEGER                                            :: wdU_SSA_2, wdV_SSA_2
    
    CALL allocate_shared_memory_dp_1D( mesh1%nV, dU_SSA_1, wdU_SSA_1)
    CALL allocate_shared_memory_dp_1D( mesh1%nV, dV_SSA_1, wdV_SSA_1)
    CALL allocate_shared_memory_dp_1D( mesh2%nV, dU_SSA_2, wdU_SSA_2)
    CALL allocate_shared_memory_dp_1D( mesh2%nV, dV_SSA_2, wdV_SSA_2)
    
    dU_SSA_1(mesh1%v1:mesh1%v2) = ice1%U_SSA(mesh1%v1:mesh1%v2) - ice1%U_SSA_old(mesh1%v1:mesh1%v2)
    dV_SSA_1(mesh1%v1:mesh1%v2) = ice1%V_SSA(mesh1%v1:mesh1%v2) - ice1%V_SSA_old(mesh1%v1:mesh1%v2)
    
    CALL MapMesh2Mesh( mesh2, map_nV, map_vi, map_w, dU_SSA_1, dU_SSA_2)
    CALL MapMesh2Mesh( mesh2, map_nV, map_vi, map_w, dV_SSA_1, dV_SSA_2)
    
    ice2%U_SSA(mesh2%v1:mesh2%v2) = ice2%U_SSA(mesh2%v1:mesh2%v2) + dU_SSA_2(mesh2%v1:mesh2%v2)
    ice2%V_SSA(mesh2%v1:mesh2%v2) = ice2%V_SSA(mesh2%v1:mesh2%v2) + dV_SSA_2(mesh2%v1:mesh2%v2)
    
    CALL deallocate_shared_memory( wdU_SSA_1)
    CALL deallocate_shared_memory( wdV_SSA_1)
    CALL deallocate_shared_memory( wdU_SSA_2)
    CALL deallocate_shared_memory( wdV_SSA_2)
    
    NULLIFY( dU_SSA_1)
    NULLIFY( dV_SSA_1)
    NULLIFY( dU_SSA_2)
    NULLIFY( dV_SSA_2)
    
  END SUBROUTINE map_SSA_solution_between_nested_meshes
  
  
  
  SUBROUTINE CreateMeshFromCartData_Nested(region)
    ! Create the first mesh, using the data from the initial file to force the resolution.
  
    ! Input variables
    TYPE(type_model_region),    INTENT(INOUT)     :: region
    
    ! Local variables
    TYPE(type_mesh)                               :: submesh
    TYPE(type_mesh)                               :: submesh_merged
    INTEGER                                       :: orientation
    REAL(dp)                                      :: xmin, xmax, ymin, ymax
    TYPE(type_mesh_config)                        :: C_mesh
    REAL(dp), DIMENSION(6)                        :: res
    REAL(dp)                                      :: smallest_res
    INTEGER                                       :: nVmap
    LOGICAL, PARAMETER                            :: debug_nested_meshes = .TRUE.
    
    ! Paralellised mesh merging only works for (integer power of 2) processes
    IF (.NOT. (par%n==2 .OR. par%n==4 .OR. par%n==8 .OR. par%n==16)) THEN
      IF (par%master) WRITE(0,'(A,I2,A)') ' ERROR: parallelised mesh creation not implemented for ', par%n, ' processors!'
      IF (par%master) WRITE(0,'(A,I2,A)') '        Valid choices; 2, 4, 8, 16'
      STOP
    END IF
    
! Determine how many nested low-resolution meshes are required
! ============================================================

    res = [ C%res_max, C%res_max_margin, C%res_max_gl, C%res_max_cf, C%res_max_mountain, C%velocity_scaling_res(C%velocity_scaling_nvals) ]
    smallest_res = MINVAL(res)
    
    IF     (smallest_res > C%nested_mesh_resolutions(1)) THEN
      WRITE(0,*) ' Mesh creation ERROR - nested meshes need a smallest resolution of at most 200 km!'
    ELSEIF (smallest_res > C%nested_mesh_resolutions(2)) THEN
      region%n_nested_meshes = 1
    ELSEIF (smallest_res > C%nested_mesh_resolutions(3)) THEN
      region%n_nested_meshes = 2
    ELSEIF (smallest_res > C%nested_mesh_resolutions(4)) THEN
      region%n_nested_meshes = 3
    ELSEIF (smallest_res > C%nested_mesh_resolutions(5)) THEN
      region%n_nested_meshes = 4
    ELSEIF (smallest_res > C%nested_mesh_resolutions(6)) THEN
      region%n_nested_meshes = 5
    ELSEIF (smallest_res > C%nested_mesh_resolutions(7)) THEN
      region%n_nested_meshes = 6
    ELSE
      region%n_nested_meshes = 7
    END IF
    
! Initialise the dummy submeshes
! ==============================

    IF (region%name == 'GRL') THEN
      orientation = 1 ! 0 - eastwest, 1 = northsouth
    ELSE
      orientation = 0
    END IF
        
    ! Determine the domain of the submesh. No previous mesh exists,
    ! so no domain partition optimisation is possible.
    IF (orientation == 0) THEN
      xmin = MINVAL(region%init%x) + REAL(par%i    ) * (MAXVAL(region%init%x) - MINVAL(region%init%x)) / REAL(par%n)
      xmax = MINVAL(region%init%x) + REAL(par%i + 1) * (MAXVAL(region%init%x) - MINVAL(region%init%x)) / REAL(par%n)
      ymin = MINVAL(region%init%y)
      ymax = MAXVAL(region%init%y)
    ELSEIF (orientation == 1) THEN
      xmin = MINVAL(region%init%x)
      xmax = MAXVAL(region%init%x)
      ymin = MINVAL(region%init%y) + REAL(par%i    ) * (MAXVAL(region%init%y) - MINVAL(region%init%y)) / REAL(par%n)
      ymax = MINVAL(region%init%y) + REAL(par%i + 1) * (MAXVAL(region%init%y) - MINVAL(region%init%y)) / REAL(par%n)
    END IF
    
    ! Allocate memory, initialise a dummy mesh, refine it, and crop the memory   
    CALL AllocateSubmesh(     submesh, 1000)
    CALL InitialiseDummyMesh( submesh, xmin, xmax, ymin, ymax)
    
! Sequentially refine the nested low-resolution meshes
! ====================================================

    IF (region%n_nested_meshes >= 1) THEN
    
      C_mesh%dz_max_ice       = MAX(C%dz_max_ice,       30000._dp)
      C_mesh%res_max          = MAX(c%res_max,          1600._dp)
      C_mesh%res_max_margin   = MAX(C%res_max_margin,   C%nested_mesh_resolutions(1))
      C_mesh%res_max_gl       = MAX(C%res_max_gl,       C%nested_mesh_resolutions(1))
      C_mesh%res_max_cf       = MAX(C%res_max_cf,       C%nested_mesh_resolutions(1))
      C_mesh%res_max_mountain = MAX(C%res_max_mountain, C%nested_mesh_resolutions(1))
      C_mesh%res_max_coast    = MAX(C%res_max_coast,    C%nested_mesh_resolutions(1))
      
      CALL ExtendSubmesh(      submesh, submesh%nV+100)
      CALL RefineMesh(         submesh, region%init, C_mesh)      
      CALL MergeSubmeshesMeta( submesh, orientation, submesh_merged)
      CALL CreateFinalMeshFromMergedSubmesh( submesh_merged, region%name, region%mesh_low_1)

      IF (par%master .AND. debug_nested_meshes) THEN
        WRITE(0,'(A)')                '   Finished creating nested mesh #1.'
        WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh_low_1%nV
        WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh_low_1%nTri
        WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh_low_1%resolution_min/1000._dp, ' - ', region%mesh_low_1%resolution_max/1000._dp, ' km'
      END IF
    
    END IF
      
    IF (region%n_nested_meshes >= 2) THEN
    
      C_mesh%dz_max_ice       = MAX(C%dz_max_ice,       14000._dp)
      C_mesh%res_max          = MAX(c%res_max,          1400._dp)
      C_mesh%res_max_margin   = MAX(C%res_max_margin,   C%nested_mesh_resolutions(2))
      C_mesh%res_max_gl       = MAX(C%res_max_gl,       C%nested_mesh_resolutions(2))
      C_mesh%res_max_cf       = MAX(C%res_max_cf,       C%nested_mesh_resolutions(2))
      C_mesh%res_max_mountain = MAX(C%res_max_mountain, C%nested_mesh_resolutions(2))
      C_mesh%res_max_coast    = MAX(C%res_max_coast,    C%nested_mesh_resolutions(2))
      
      CALL ExtendSubmesh(      submesh, submesh%nV+100)
      CALL RefineMesh(         submesh, region%init, C_mesh)
      CALL MergeSubmeshesMeta( submesh, orientation, submesh_merged)
      CALL CreateFinalMeshFromMergedSubmesh( submesh_merged, region%name, region%mesh_low_2)
      
      ! Calculate the mapping arrays
      CALL GetMesh2MeshMap(region%mesh_low_1, region%mesh_low_2, region%map_12_nV, region%map_12_vi, region%map_12_w, region%wmap_12_nV, region%wmap_12_vi, region%wmap_12_w)
      CALL GetMesh2MeshMap(region%mesh_low_2, region%mesh_low_1, region%map_21_nV, region%map_21_vi, region%map_21_w, region%wmap_21_nV, region%wmap_21_vi, region%wmap_21_w)
      CALL sync

      IF (par%master .AND. debug_nested_meshes) THEN
        WRITE(0,'(A)')                '   Finished creating nested mesh #2.'
        WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh_low_2%nV
        WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh_low_2%nTri
        WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh_low_2%resolution_min/1000._dp, ' - ', region%mesh_low_2%resolution_max/1000._dp, ' km'
      END IF
    
    END IF
      
    IF (region%n_nested_meshes >= 3) THEN
    
      C_mesh%dz_max_ice       = MAX(C%dz_max_ice,       12000._dp)
      C_mesh%res_max          = MAX(c%res_max,          1200._dp)
      C_mesh%res_max_margin   = MAX(C%res_max_margin,   C%nested_mesh_resolutions(3))
      C_mesh%res_max_gl       = MAX(C%res_max_gl,       C%nested_mesh_resolutions(3))
      C_mesh%res_max_cf       = MAX(C%res_max_cf,       C%nested_mesh_resolutions(3))
      C_mesh%res_max_mountain = MAX(C%res_max_mountain, C%nested_mesh_resolutions(3))
      C_mesh%res_max_coast    = MAX(C%res_max_coast,    C%nested_mesh_resolutions(3))
      
      CALL ExtendSubmesh(      submesh, submesh%nV+100)
      CALL RefineMesh(         submesh, region%init, C_mesh)
      CALL MergeSubmeshesMeta( submesh, orientation, submesh_merged)
      CALL CreateFinalMeshFromMergedSubmesh( submesh_merged, region%name, region%mesh_low_3)
      
      ! Calculate the mapping arrays
      CALL GetMesh2Meshmap(region%mesh_low_2, region%mesh_low_3, region%map_23_nV, region%map_23_vi, region%map_23_w, region%wmap_23_nV, region%wmap_23_vi, region%wmap_23_w)
      CALL GetMesh2Meshmap(region%mesh_low_3, region%mesh_low_2, region%map_32_nV, region%map_32_vi, region%map_32_w, region%wmap_32_nV, region%wmap_32_vi, region%wmap_32_w)
      CALL sync

      IF (par%master .AND. debug_nested_meshes) THEN
        WRITE(0,'(A)')                '   Finished creating nested mesh #3.'
        WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh_low_3%nV
        WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh_low_3%nTri
        WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh_low_3%resolution_min/1000._dp, ' - ', region%mesh_low_3%resolution_max/1000._dp, ' km'
      END IF
    
    END IF
      
    IF (region%n_nested_meshes >= 4) THEN
    
      C_mesh%dz_max_ice       = MAX(C%dz_max_ice,       10000._dp)
      C_mesh%res_max          = MAX(c%res_max,          1000._dp)
      C_mesh%res_max_margin   = MAX(C%res_max_margin,   C%nested_mesh_resolutions(4))
      C_mesh%res_max_gl       = MAX(C%res_max_gl,       C%nested_mesh_resolutions(4))
      C_mesh%res_max_cf       = MAX(C%res_max_cf,       C%nested_mesh_resolutions(4))
      C_mesh%res_max_mountain = MAX(C%res_max_mountain, C%nested_mesh_resolutions(4))
      C_mesh%res_max_coast    = MAX(C%res_max_coast,    C%nested_mesh_resolutions(4))
      
      CALL ExtendSubmesh(      submesh, submesh%nV+100)
      CALL RefineMesh(         submesh, region%init, C_mesh)
      CALL MergeSubmeshesMeta( submesh, orientation, submesh_merged)
      CALL CreateFinalMeshFromMergedSubmesh( submesh_merged, region%name, region%mesh_low_4)
      
      ! Calculate the mapping arrays
      CALL GetMesh2MeshMap(region%mesh_low_3, region%mesh_low_4, region%map_34_nV, region%map_34_vi, region%map_34_w, region%wmap_34_nV, region%wmap_34_vi, region%wmap_34_w)
      CALL GetMesh2MeshMap(region%mesh_low_4, region%mesh_low_3, region%map_43_nV, region%map_43_vi, region%map_43_w, region%wmap_43_nV, region%wmap_43_vi, region%wmap_43_w)
      CALL sync

      IF (par%master .AND. debug_nested_meshes) THEN
        WRITE(0,'(A)')                '   Finished creating nested mesh #4.'
        WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh_low_4%nV
        WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh_low_4%nTri
        WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh_low_4%resolution_min/1000._dp, ' - ', region%mesh_low_4%resolution_max/1000._dp, ' km'
      END IF
    
    END IF
      
    IF (region%n_nested_meshes >= 5) THEN
    
      C_mesh%dz_max_ice       = MAX(C%dz_max_ice,       8000._dp)
      C_mesh%res_max          = MAX(c%res_max,          800._dp)
      C_mesh%res_max_margin   = MAX(C%res_max_margin,   C%nested_mesh_resolutions(5))
      C_mesh%res_max_gl       = MAX(C%res_max_gl,       C%nested_mesh_resolutions(5))
      C_mesh%res_max_cf       = MAX(C%res_max_cf,       C%nested_mesh_resolutions(5))
      C_mesh%res_max_mountain = MAX(C%res_max_mountain, C%nested_mesh_resolutions(5))
      C_mesh%res_max_coast    = MAX(C%res_max_coast,    C%nested_mesh_resolutions(5))
      
      CALL ExtendSubmesh(      submesh, submesh%nV+100)
      CALL RefineMesh(         submesh, region%init, C_mesh)
      CALL MergeSubmeshesMeta( submesh, orientation, submesh_merged)
      CALL CreateFinalMeshFromMergedSubmesh( submesh_merged, region%name, region%mesh_low_5)
      
      ! Calculate the mapping arrays
      CALL GetMesh2MeshMap(region%mesh_low_4, region%mesh_low_5, region%map_45_nV, region%map_45_vi, region%map_45_w, region%wmap_45_nV, region%wmap_45_vi, region%wmap_45_w)
      CALL GetMesh2MeshMap(region%mesh_low_5, region%mesh_low_4, region%map_54_nV, region%map_54_vi, region%map_54_w, region%wmap_54_nV, region%wmap_54_vi, region%wmap_54_w)
      CALL sync

      IF (par%master .AND. debug_nested_meshes) THEN
        WRITE(0,'(A)')                '   Finished creating nested mesh #5.'
        WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh_low_5%nV
        WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh_low_5%nTri
        WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh_low_5%resolution_min/1000._dp, ' - ', region%mesh_low_5%resolution_max/1000._dp, ' km'
      END IF
    
    END IF
      
    IF (region%n_nested_meshes >= 6) THEN
    
      C_mesh%dz_max_ice       = MAX(C%dz_max_ice,       6000._dp)
      C_mesh%res_max          = MAX(c%res_max,          600._dp)
      C_mesh%res_max_margin   = MAX(C%res_max_margin,   C%nested_mesh_resolutions(6))
      C_mesh%res_max_gl       = MAX(C%res_max_gl,       C%nested_mesh_resolutions(6))
      C_mesh%res_max_cf       = MAX(C%res_max_cf,       C%nested_mesh_resolutions(6))
      C_mesh%res_max_mountain = MAX(C%res_max_mountain, C%nested_mesh_resolutions(6))
      C_mesh%res_max_coast    = MAX(C%res_max_coast,    C%nested_mesh_resolutions(6))
      
      CALL ExtendSubmesh(      submesh, submesh%nV+100)
      CALL RefineMesh(         submesh, region%init, C_mesh)
      CALL MergeSubmeshesMeta( submesh, orientation, submesh_merged)
      CALL CreateFinalMeshFromMergedSubmesh( submesh_merged, region%name, region%mesh_low_6)
      
      ! Calculate the mapping arrays
      CALL GetMesh2MeshMap(region%mesh_low_5, region%mesh_low_6, region%map_56_nV, region%map_56_vi, region%map_56_w, region%wmap_56_nV, region%wmap_56_vi, region%wmap_56_w)
      CALL GetMesh2MeshMap(region%mesh_low_6, region%mesh_low_5, region%map_65_nV, region%map_65_vi, region%map_65_w, region%wmap_65_nV, region%wmap_65_vi, region%wmap_65_w)
      CALL sync

      IF (par%master .AND. debug_nested_meshes) THEN
        WRITE(0,'(A)')                '   Finished creating nested mesh #6.'
        WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh_low_6%nV
        WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh_low_6%nTri
        WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh_low_6%resolution_min/1000._dp, ' - ', region%mesh_low_6%resolution_max/1000._dp, ' km'
      END IF
    
    END IF
      
    IF (region%n_nested_meshes >= 7) THEN
    
      C_mesh%dz_max_ice       = MAX(C%dz_max_ice,       3000._dp)
      C_mesh%res_max          = MAX(c%res_max,          300._dp)
      C_mesh%res_max_margin   = MAX(C%res_max_margin,   C%nested_mesh_resolutions(7))
      C_mesh%res_max_gl       = MAX(C%res_max_gl,       C%nested_mesh_resolutions(7))
      C_mesh%res_max_cf       = MAX(C%res_max_cf,       C%nested_mesh_resolutions(7))
      C_mesh%res_max_mountain = MAX(C%res_max_mountain, C%nested_mesh_resolutions(7))
      C_mesh%res_max_coast    = MAX(C%res_max_coast,    C%nested_mesh_resolutions(7))
      
      CALL ExtendSubmesh(      submesh, submesh%nV+100)
      CALL RefineMesh(         submesh, region%init, C_mesh)
      CALL MergeSubmeshesMeta( submesh, orientation, submesh_merged)
      CALL CreateFinalMeshFromMergedSubmesh( submesh_merged, region%name, region%mesh_low_7)
      
      ! Calculate the mapping arrays
      CALL GetMesh2MeshMap(region%mesh_low_6, region%mesh_low_7, region%map_67_nV, region%map_67_vi, region%map_67_w, region%wmap_67_nV, region%wmap_67_vi, region%wmap_67_w)
      CALL GetMesh2MeshMap(region%mesh_low_7, region%mesh_low_6, region%map_76_nV, region%map_76_vi, region%map_76_w, region%wmap_76_nV, region%wmap_76_vi, region%wmap_76_w)
      CALL sync

      IF (par%master .AND. debug_nested_meshes) THEN
        WRITE(0,'(A)')                '   Finished creating nested mesh #7.'
        WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh_low_7%nV
        WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh_low_7%nTri
        WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh_low_7%resolution_min/1000._dp, ' - ', region%mesh_low_7%resolution_max/1000._dp, ' km'
      END IF
    
    END IF
    
! Final model mesh
! ================
    
    C_mesh%dz_max_ice       = C%dz_max_ice
    C_mesh%res_max          = C%res_max
    C_mesh%res_max_margin   = C%res_max_margin
    C_mesh%res_max_gl       = C%res_max_gl
    C_mesh%res_max_cf       = C%res_max_cf
    C_mesh%res_max_mountain = C%res_max_mountain
    C_mesh%res_max_coast    = C%res_max_coast
     
    CALL ExtendSubmesh(      submesh, submesh%nV + 100)
    CALL RefineMesh(         submesh, region%init, C_mesh)   
    CALL MergeSubmeshesMeta( submesh, orientation, submesh_merged)
    CALL CreateFinalMeshFromMergedSubmesh( submesh_merged, region%name, region%mesh)    
    CALL DeallocateSubmesh( submesh)
    
    ! Determine how much shared memory to allocate for the final set of mapping arrays    
    IF (region%n_nested_meshes == 1) nVmap = region%mesh_low_1%nV
    IF (region%n_nested_meshes == 2) nVmap = region%mesh_low_2%nV
    IF (region%n_nested_meshes == 3) nVmap = region%mesh_low_3%nV
    IF (region%n_nested_meshes == 4) nVmap = region%mesh_low_4%nV
    IF (region%n_nested_meshes == 5) nVmap = region%mesh_low_5%nV
    IF (region%n_nested_meshes == 6) nVmap = region%mesh_low_6%nV
    IF (region%n_nested_meshes == 7) nVmap = region%mesh_low_7%nV
    
    ! Calculate the final set of mapping arrays
    IF     (region%n_nested_meshes == 1 .AND. par%master) THEN
      CALL GetMesh2MeshMap(region%mesh_low_1, region%mesh, &
        region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w,  region%wmap_nest_final_nV, region%wmap_nest_final_vi, region%wmap_nest_final_w)
      CALL GetMesh2MeshMap(region%mesh, region%mesh_low_1, &
        region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w,  region%wmap_final_nest_nV, region%wmap_final_nest_vi, region%wmap_final_nest_w)
    ELSEIF (region%n_nested_meshes == 2 .AND. par%master) THEN
      CALL GetMesh2MeshMap(region%mesh_low_2, region%mesh, &
        region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w,  region%wmap_nest_final_nV, region%wmap_nest_final_vi, region%wmap_nest_final_w)
      CALL GetMesh2MeshMap(region%mesh, region%mesh_low_2, &
        region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w,  region%wmap_final_nest_nV, region%wmap_final_nest_vi, region%wmap_final_nest_w)
    ELSEIF (region%n_nested_meshes == 3 .AND. par%master) THEN
      CALL GetMesh2MeshMap(region%mesh_low_3, region%mesh, &
        region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w,  region%wmap_nest_final_nV, region%wmap_nest_final_vi, region%wmap_nest_final_w)
      CALL GetMesh2MeshMap(region%mesh, region%mesh_low_3, &
        region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w,  region%wmap_final_nest_nV, region%wmap_final_nest_vi, region%wmap_final_nest_w)
    ELSEIF (region%n_nested_meshes == 4 .AND. par%master) THEN
      CALL GetMesh2MeshMap(region%mesh_low_4, region%mesh, &
        region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w,  region%wmap_nest_final_nV, region%wmap_nest_final_vi, region%wmap_nest_final_w)
      CALL GetMesh2MeshMap(region%mesh, region%mesh_low_4, &
        region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w,  region%wmap_final_nest_nV, region%wmap_final_nest_vi, region%wmap_final_nest_w)
    ELSEIF (region%n_nested_meshes == 5 .AND. par%master) THEN
      CALL GetMesh2MeshMap(region%mesh_low_5, region%mesh, &
        region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w,  region%wmap_nest_final_nV, region%wmap_nest_final_vi, region%wmap_nest_final_w)
      CALL GetMesh2MeshMap(region%mesh, region%mesh_low_5, &
        region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w,  region%wmap_final_nest_nV, region%wmap_final_nest_vi, region%wmap_final_nest_w)
    ELSEIF (region%n_nested_meshes == 6 .AND. par%master) THEN
      CALL GetMesh2MeshMap(region%mesh_low_6, region%mesh, &
        region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w,  region%wmap_nest_final_nV, region%wmap_nest_final_vi, region%wmap_nest_final_w)
      CALL GetMesh2MeshMap(region%mesh, region%mesh_low_6, &
        region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w,  region%wmap_final_nest_nV, region%wmap_final_nest_vi, region%wmap_final_nest_w)
    ELSEIF (region%n_nested_meshes == 7 .AND. par%master) THEN
      CALL GetMesh2MeshMap(region%mesh_low_7, region%mesh, &
        region%map_nest_final_nV, region%map_nest_final_vi, region%map_nest_final_w,  region%wmap_nest_final_nV, region%wmap_nest_final_vi, region%wmap_nest_final_w)
      CALL GetMesh2MeshMap(region%mesh, region%mesh_low_7, &
        region%map_final_nest_nV, region%map_final_nest_vi, region%map_final_nest_w,  region%wmap_final_nest_nV, region%wmap_final_nest_vi, region%wmap_final_nest_w)
    END IF
    CALL sync

    IF (par%master) THEN
      WRITE(0,'(A)')                '   Finished creating final mesh.'
      WRITE(0,'(A,I6)')             '    Vertices  : ', region%mesh%nV
      WRITE(0,'(A,I6)')             '    Triangles : ', region%mesh%nTri
      WRITE(0,'(A,F7.1,A,F7.1,A)')  '    Resolution: ', region%mesh%resolution_min/1000._dp, ' - ', region%mesh%resolution_max/1000._dp, ' km'
    END IF
    
  END SUBROUTINE CreateMeshFromCartData_Nested
