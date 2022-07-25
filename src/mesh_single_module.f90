module mesh_single_module
contains
  ! == Mesh creation routines ==
  SUBROUTINE create_single_mesh_from_cart_data( region)
    ! Create the first mesh, using the data from the initial file to force the resolution.
    use mesh_creation_module

    IMPLICIT NONE

    ! Input variables
    TYPE(type_model_region),    INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_single_mesh_from_cart_data'
    INTEGER                                       :: orientation
    TYPE(type_mesh)                               :: submesh
    REAL(dp)                                      :: xmin, xmax, ymin, ymax
    REAL(dp)                                      :: res_min_inc
    CHARACTER(LEN=2)                              :: str_processid

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%master) then
      write(*,"(A)") '  Creating the first mesh...'
    end if

    ! Orientation of domain partitioning: east-west for GRL, north-south everywhere else
    IF (region%name == 'GRL') THEN
      orientation = 1
    ELSE
      orientation = 0
    END IF

    ! Determine the domain of this process' submesh.
    ymin = MINVAL(region%refgeo_init%grid%y)
    ymax = MAXVAL(region%refgeo_init%grid%y)
    xmin = MINVAL(region%refgeo_init%grid%x)
    xmax = MAXVAL(region%refgeo_init%grid%x)

    ! Allocate memory and initialise a dummy mesh
    CALL allocate_submesh_primary( submesh, region%name, 10, 20, C%nconmax)
    CALL initialise_dummy_mesh(    submesh, xmin, xmax, ymin, ymax)
    CALL perturb_dummy_mesh(       submesh, 0)


    if (par%master) then
      res_min_inc = C%res_max * 2._dp

      DO WHILE (res_min_inc > C%res_min)

        ! Increase resolution
        res_min_inc = res_min_inc / 2._dp

        ! Determine resolutions
        submesh%res_min          = MAX( C%res_min,          res_min_inc)
        submesh%res_max_margin   = MAX( C%res_max_margin,   res_min_inc)
        submesh%res_max_gl       = MAX( C%res_max_gl,       res_min_inc)
        submesh%res_max_cf       = MAX( C%res_max_cf,       res_min_inc)
        submesh%res_max_mountain = MAX( C%res_max_mountain, res_min_inc)
        submesh%res_max_coast    = MAX( C%res_max_coast,    res_min_inc)

        if (debug_mesh_creation) then
          write(*,"(A,I3,A,F4.1)") '  Process ', par%i, ' refining submesh to ', submesh%res_max_gl, ' km...'
        end if

        ! Refine the process submesh
        CALL refine_mesh( submesh, region%refgeo_init)

        ! Split any new triangles (added during alignment) that are too sharp
        ! CALL refine_submesh_geo_only( submesh)

        ! Smooth the submesh using Lloyd' algorithm
        CALL Lloyds_algorithm_single_iteration_submesh( submesh)

        ! After the last refinement step, apply Lloyds algorithm two more times, because we can.
        IF (res_min_inc <= C%res_min) THEN
        CALL Lloyds_algorithm_single_iteration_submesh( submesh)
        CALL Lloyds_algorithm_single_iteration_submesh( submesh)
        END IF

        ! Write submesh to text file for debugging
        write(str_processid,'(I2)') par%i; str_processid = adjustl(str_processid)
        if (debug_mesh_creation) then
          CALL write_mesh_to_text_file( submesh, 'submesh_proc_' // TRIM(str_processid) // '.txt')
        end if

        ! Check if everything went correctly
        CALL check_mesh( submesh)

      END DO

    end if

    IF (debug_mesh_creation .AND. par%master) WRITE(0,*) '  Creating final mesh...'
    CALL create_final_mesh_from_merged_submesh( submesh, region%mesh)

    IF (par%master) THEN
       write(*,"(A,I6)")             '   Vertices  : ', region%mesh%nV
       write(*,"(A,I6)")             '   Triangles : ', region%mesh%nTri
       write(*,"(A,F7.1,A,F7.1,A)")  '   Resolution: ', region%mesh%resolution_min/1000._dp, ' - ', region%mesh%resolution_max/1000._dp, ' km'
       write(*,"(A)")                '  Finished creating final mesh.'
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_single_mesh_from_cart_data

  SUBROUTINE create_new_mesh_single( region)
    use mesh_update_module ! use everything

    IMPLICIT NONE

    ! Input variables
    TYPE(type_model_region),    INTENT(INOUT)     :: region

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                 :: routine_name = 'create_new_mesh_single'
    INTEGER                                       :: vi
    REAL(dp), DIMENSION(:    ), POINTER           ::  d2dx2,  d2dxdy,  d2dy2
    INTEGER                                       :: wd2dx2, wd2dxdy, wd2dy2
    INTEGER                                       :: orientation, it
    TYPE(type_mesh)                               :: submesh
    REAL(dp)                                      :: xmin, xmax, ymin, ymax
    REAL(dp)                                      :: res_min_inc
    CHARACTER(LEN=2)                              :: str_processid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Screen meesage
    if (par%master) then
      if (C%do_time_display) then
        if (mod(region%time-region%dt,C%dt_output) /= 0._dp) then
          write(*,"(A)",advance="yes") repeat(c_backspace,17) // &
                                       ' - mesh time!    '
        else
          ! Output took care of advancing a newline.
        end if
      end if
      write(*,"(A)",advance="yes") '  Creating a new mesh for region ' &
                                   // TRIM(region%mesh%region_name) // '...'
    end if

    ! Orientation of domain partitioning: east-west for GRL, north-south everywhere else
    IF (region%name == 'GRL') THEN
      orientation = 1
    ELSE
      orientation = 0
    END IF

    ! Make the ice margin 1 row of elements wider. This slightly increases the vertex count of
    ! the new mesh, but means the new mesh stays Fit much longer.
    CALL widen_high_res_zones( region%ice, region%mesh, region%time)

    ! Calculate surface curvature
    allocate( d2dx2 (region%mesh%vi1:region%mesh%vi2))
    allocate( d2dxdy(region%mesh%vi1:region%mesh%vi2))
    allocate( d2dy2 (region%mesh%vi1:region%mesh%vi2))
    CALL d2dx2_a_to_a_2D(  region%mesh, region%ice%Hs_a, d2dx2 )
    CALL d2dxdy_a_to_a_2D( region%mesh, region%ice%Hs_a, d2dxdy)
    CALL d2dy2_a_to_a_2D(  region%mesh, region%ice%Hs_a, d2dy2 )
    DO vi = region%mesh%vi1, region%mesh%vi2
      region%ice%surf_curv( vi) = MAX(-1E-6, MIN(1E-6, SQRT(d2dx2( vi)**2 + d2dy2( vi)**2 + d2dxdy( vi)**2)))
    END DO
    call allgather_array(region%ice%surf_curv)
    deallocate( d2dx2)
    deallocate( d2dxdy)
    deallocate( d2dy2)

    ! Determine the domain of this process' submesh (based on distribution of vertices in
    ! the previous mesh, which works very well for workload balancing)
    ymin = region%mesh%ymin
    ymax = region%mesh%ymax
    xmin = region%mesh%xmin
    xmax = region%mesh%xmax

    ! Allocate memory, initialise a dummy mesh, refine it, and crop the memory
    CALL allocate_submesh_primary( submesh, region%name, 10, 20, C%nconmax)
    CALL initialise_dummy_mesh(    submesh, xmin, xmax, ymin, ymax)
    CALL perturb_dummy_mesh(       submesh, 1 - region%mesh%perturb_dir)

    if (par%master) then
      res_min_inc = C%res_max * 2._dp

      it = 0
      DO WHILE (res_min_inc > C%res_min)
        it = it + 1

        ! Increase resolution
        res_min_inc = res_min_inc / 2._dp

        ! Determine resolutions
        submesh%res_min          = MAX( C%res_min,          res_min_inc)
        submesh%res_max_margin   = MAX( C%res_max_margin,   res_min_inc)
        submesh%res_max_gl       = MAX( C%res_max_gl,       res_min_inc)
        submesh%res_max_cf       = MAX( C%res_max_cf,       res_min_inc)
        submesh%res_max_mountain = MAX( C%res_max_mountain, res_min_inc)
        submesh%res_max_coast    = MAX( C%res_max_coast,    res_min_inc)

        ! Refine the process submesh
        CALL refine_mesh( submesh, region%mesh, region%ice, region%refgeo_PD)

        ! Split any new triangles (added during alignment) that are too sharp
        !CALL refine_submesh_geo_only( submesh)

        ! Smooth the submesh using Lloyd' algorithm
        CALL Lloyds_algorithm_single_iteration_submesh( submesh)

        ! Write submesh to text file for debugging
        WRITE(str_processid,'(I2)') par%i;   str_processid = ADJUSTL(str_processid)
        IF (debug_mesh_creation) CALL write_mesh_to_text_file( submesh, 'submesh_proc_' // TRIM(str_processid) // '.txt')

        ! Check if everything went correctly
        CALL check_mesh( submesh)

      END DO
    end if

    IF (debug_mesh_creation .AND. par%master) WRITE(0,*) '  Creating final mesh...'
    CALL create_final_mesh_from_merged_submesh( submesh, region%mesh_new)

    IF (par%master) THEN
      write(*,"(A,I6)")            '   Vertices  : ', region%mesh_new%nV
      write(*,"(A,I6)")            '   Triangles : ', region%mesh_new%nTri
      write(*,"(A,F7.1,A,F7.1,A)") '   Resolution: ', region%mesh_new%resolution_min/1000._dp, ' - ', region%mesh_new%resolution_max/1000._dp, ' km'
      write(*,"(A)")               '  Finished creating final mesh.'
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_new_mesh_single
end module
