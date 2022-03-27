MODULE mesh_memory_module
  ! Routines for allocating, deallocating, extending and cropping the memory for the mesh data.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: tMat, MatDestroy, perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list, &
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

  ! Import specific functionality
  USE parallel_module,                 ONLY: adapt_shared_int_1D,    adapt_shared_dp_1D, &
                                             adapt_shared_int_2D,    adapt_shared_dp_2D, &
                                             adapt_shared_int_3D,    adapt_shared_dp_3D, &
                                             adapt_shared_bool_1D, &
                                             allocate_shared_dist_int_0D, allocate_shared_dist_dp_0D, &
                                             allocate_shared_dist_int_1D, allocate_shared_dist_dp_1D, &
                                             allocate_shared_dist_int_2D, allocate_shared_dist_dp_2D, &
                                             allocate_shared_dist_int_3D, allocate_shared_dist_dp_3D, &
                                             allocate_shared_dist_bool_1D, &
                                             adapt_shared_dist_int_1D,    adapt_shared_dist_dp_1D, &
                                             adapt_shared_dist_int_2D,    adapt_shared_dist_dp_2D, &
                                             adapt_shared_dist_int_3D,    adapt_shared_dist_dp_3D, &
                                             adapt_shared_dist_bool_1D, &
                                             share_memory_access_int_0D, share_memory_access_dp_0D, &
                                             share_memory_access_int_1D, share_memory_access_dp_1D, &
                                             share_memory_access_int_2D, share_memory_access_dp_2D, &
                                             share_memory_access_int_3D, share_memory_access_dp_3D
  USE data_types_module,               ONLY: type_mesh, type_mesh_new
  USE sparse_matrix_module,            ONLY: deallocate_matrix_CSR
  use reallocate_mod,                  only: reallocate

  IMPLICIT NONE

  interface extend_submesh_primary
    procedure :: extend_submesh_primary
    procedure :: extend_submesh_primary_new
  end interface
  interface crop_submesh_primary
    procedure :: crop_submesh_primary
    procedure :: crop_submesh_primary_new
  end interface crop_submesh_primary
  interface deallocate_submesh_primary
    procedure :: deallocate_submesh_primary
    procedure :: deallocate_submesh_primary_new
  end interface
  interface share_submesh_access
    procedure :: share_submesh_access
    procedure :: share_submesh_access_new
  end interface
  interface allocate_submesh_primary
    procedure :: allocate_submesh_primary
    procedure :: allocate_submesh_primary_new
  end interface
  interface move_data_from_submesh_to_mesh
    procedure :: move_data_from_submesh_to_mesh
    procedure :: move_data_from_submesh_to_mesh_new
  end interface
CONTAINS
  
  SUBROUTINE allocate_mesh_primary(       mesh, region_name, nV_mem, nTri_mem, nC_mem)
    ! Allocate memory for mesh data
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    CHARACTER(LEN=3),                INTENT(IN)        :: region_name
    INTEGER,                         INTENT(IN)        :: nV_mem, nTri_mem, nC_mem
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_mesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    mesh%region_name = region_name
    
    CALL allocate_shared_dp_0D(  mesh%lambda_M,         mesh%wlambda_M        )
    CALL allocate_shared_dp_0D(  mesh%phi_M,            mesh%wphi_M           )
    CALL allocate_shared_dp_0D(  mesh%alpha_stereo,     mesh%walpha_stereo    )
    
    IF (mesh%region_name     == 'NAM') THEN
      mesh%lambda_M     = C%lambda_M_NAM
      mesh%phi_M        = C%phi_M_NAM
      mesh%alpha_stereo = C%alpha_stereo_NAM
    ELSEIF (mesh%region_name == 'EAS') THEN
      mesh%lambda_M     = C%lambda_M_EAS
      mesh%phi_M        = C%phi_M_EAS
      mesh%alpha_stereo = C%alpha_stereo_EAS
    ELSEIF (mesh%region_name == 'GRL') THEN
      mesh%lambda_M     = C%lambda_M_GRL
      mesh%phi_M        = C%phi_M_GRL
      mesh%alpha_stereo = C%alpha_stereo_GRL
    ELSEIF (mesh%region_name == 'ANT') THEN
      mesh%lambda_M     = C%lambda_M_ANT
      mesh%phi_M        = C%phi_M_ANT
      mesh%alpha_stereo = C%alpha_stereo_ANT
    END IF
    
    CALL allocate_shared_dp_0D(  mesh%xmin,             mesh%wxmin            )
    CALL allocate_shared_dp_0D(  mesh%xmax,             mesh%wxmax            )
    CALL allocate_shared_dp_0D(  mesh%ymin,             mesh%wymin            )
    CALL allocate_shared_dp_0D(  mesh%ymax,             mesh%wymax            )
    CALL allocate_shared_dp_0D(  mesh%tol_dist,         mesh%wtol_dist        )
    CALL allocate_shared_int_0D( mesh%nC_mem,           mesh%wnC_mem          )
    CALL allocate_shared_int_0D( mesh%nV_mem,           mesh%wnV_mem          )
    CALL allocate_shared_int_0D( mesh%nTri_mem,         mesh%wnTri_mem        )
    CALL allocate_shared_int_0D( mesh%nV,               mesh%wnV              )
    CALL allocate_shared_int_0D( mesh%nTri,             mesh%wnTri            )
    CALL allocate_shared_int_0D( mesh%perturb_dir,      mesh%wperturb_dir     )
    CALL allocate_shared_dp_0D(  mesh%alpha_min,        mesh%walpha_min       )
    CALL allocate_shared_dp_0D(  mesh%dz_max_ice,       mesh%wdz_max_ice      )
    CALL allocate_shared_dp_0D(  mesh%res_max,          mesh%wres_max         )
    CALL allocate_shared_dp_0D(  mesh%res_max_margin,   mesh%wres_max_margin  )
    CALL allocate_shared_dp_0D(  mesh%res_max_gl,       mesh%wres_max_gl      )
    CALL allocate_shared_dp_0D(  mesh%res_max_cf,       mesh%wres_max_cf      )
    CALL allocate_shared_dp_0D(  mesh%res_max_mountain, mesh%wres_max_mountain)
    CALL allocate_shared_dp_0D(  mesh%res_max_coast,    mesh%wres_max_coast   )
    CALL allocate_shared_dp_0D(  mesh%res_min,          mesh%wres_min         )
    CALL allocate_shared_dp_0D(  mesh%resolution_min,   mesh%wresolution_min  )
    CALL allocate_shared_dp_0D(  mesh%resolution_max,   mesh%wresolution_max  )
    
    IF (par%master) THEN
      mesh%nV_mem   = nV_mem
      mesh%nTri_mem = nTri_mem
      mesh%nC_mem   = nC_mem
    END IF
    CALL sync
    
    CALL allocate_shared_dp_2D(  nV_mem,   2,      mesh%V,               mesh%wV              )
    CALL allocate_shared_int_1D( nV_mem,           mesh%nC,              mesh%wnC             )
    CALL allocate_shared_int_2D( nV_mem,   nC_mem, mesh%C,               mesh%wC              )
    CALL allocate_shared_int_1D( nV_mem,           mesh%niTri,           mesh%wniTri          )
    CALL allocate_shared_int_2D( nV_mem,   nC_mem, mesh%iTri,            mesh%wiTri           )
    CALL allocate_shared_int_1D( nV_mem,           mesh%edge_index,      mesh%wedge_index     )
    CALL allocate_shared_int_1D( nV_mem,           mesh%mesh_old_ti_in,  mesh%wmesh_old_ti_in )
    
    CALL allocate_shared_int_2D( nTri_mem, 3,      mesh%Tri,             mesh%wTri            )
    CALL allocate_shared_dp_2D(  nTri_mem, 2,      mesh%Tricc,           mesh%wTricc          )
    CALL allocate_shared_int_2D( nTri_mem, 3,      mesh%TriC,            mesh%wTriC           )
    CALL allocate_shared_int_1D( nTri_mem,         mesh%Tri_edge_index,  mesh%wTri_edge_index )
    
    CALL allocate_shared_int_2D( nTri_mem, 2,      mesh%Triflip,         mesh%wTriflip        )
    CALL allocate_shared_int_1D( nTri_mem,         mesh%RefMap,          mesh%wRefMap         )
    CALL allocate_shared_int_1D( nTri_mem,         mesh%RefStack,        mesh%wRefStack       )
    CALL allocate_shared_int_0D(                   mesh%RefStackN,       mesh%wRefStackN      )
    
    ! Distributed shared memory for FloodFill maps and stacks
    CALL allocate_shared_dist_int_1D( nV_mem,      mesh%VMap,            mesh%wVMap           )
    CALL allocate_shared_dist_int_1D( nV_mem,      mesh%VStack1,         mesh%wVStack1        )
    CALL allocate_shared_dist_int_1D( nV_mem,      mesh%VStack2,         mesh%wVStack2        )
    CALL allocate_shared_dist_int_1D( nTri_mem,    mesh%TriMap,          mesh%wTriMap         )
    CALL allocate_shared_dist_int_1D( nTri_mem,    mesh%TriStack1,       mesh%wTriStack1      )
    CALL allocate_shared_dist_int_1D( nTri_mem,    mesh%TriStack2,       mesh%wTriStack2      )
    
    ! POI stuff
    CALL allocate_shared_int_0D(  mesh%nPOI, mesh%wnPOI)
    
    IF (mesh%region_name     == 'NAM') THEN
      mesh%nPOI         = C%nPOI_NAM
    ELSEIF (mesh%region_name == 'EAS') THEN
      mesh%nPOI         = C%nPOI_EAS
    ELSEIF (mesh%region_name == 'GRL') THEN
      mesh%nPOI         = C%nPOI_GRL
    ELSEIF (mesh%region_name == 'ANT') THEN
      mesh%lambda_M     = C%lambda_M_ANT
      mesh%nPOI         = C%nPOI_ANT
    END IF
    
    CALL allocate_shared_dp_2D(  mesh%nPOI, 2, mesh%POI_coordinates,           mesh%wPOI_coordinates          )
    CALL allocate_shared_dp_2D(  mesh%nPOI, 2, mesh%POI_XY_coordinates,        mesh%wPOI_XY_coordinates       )
    CALL allocate_shared_dp_1D(  mesh%nPOI,    mesh%POI_resolutions,           mesh%wPOI_resolutions          )
    CALL allocate_shared_int_2D( mesh%nPOI, 3, mesh%POI_vi,                    mesh%wPOI_vi                   )
    CALL allocate_shared_dp_2D(  mesh%nPOI, 3, mesh%POI_w,                     mesh%wPOI_w                    )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 52)

  END SUBROUTINE allocate_mesh_primary
  SUBROUTINE extend_submesh_primary_new(      mesh, nV_mem_new, nTri_mem_new)
    ! For when we didn't allocate enough. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (bigger) one, and copy the data back.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh_new),                 INTENT(INOUT)     :: mesh
    INTEGER,                         INTENT(IN)        :: nV_mem_new, nTri_mem_new
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extend_submesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
 
    mesh%nV_mem   = nV_mem_new
    mesh%nTri_mem = nTri_mem_new
       
    call reallocate(mesh%V             , nv_mem_new,           2)
    call reallocate(mesh%V             , nV_mem_new,           2)
    call reallocate(mesh%nC            , nV_mem_new             )
    call reallocate(mesh%C             , nV_mem_new, mesh%nC_mem)
    call reallocate(mesh%niTri         , nV_mem_new             )
    call reallocate(mesh%iTri          , nV_mem_new, mesh%nC_mem)
    call reallocate(mesh%edge_index    , nV_mem_new             )
    call reallocate(mesh%mesh_old_ti_in, nV_mem_new             )
    
    mesh%mesh_old_ti_in(mesh%nV+1:nV_mem_new) = 1
    
    call reallocate( mesh%Tri,            nTri_mem_new, 3)
    call reallocate( mesh%Tricc,          nTri_mem_new, 2)
    call reallocate( mesh%TriC,           nTri_mem_new, 3)
    call reallocate( mesh%Tri_edge_index, nTri_mem_new   )

    call reallocate( mesh%Triflip,        nTri_mem_new, 2)
    call reallocate( mesh%RefMap,         nTri_mem_new   )
    call reallocate( mesh%RefStack,       nTri_mem_new   )
    
    ! Distributed shared memory for FloodFill maps and stacks
    call reallocate( mesh%VMap,      nV_mem_new   )
    call reallocate( mesh%VStack1,   nV_mem_new   )
    call reallocate( mesh%VStack2,   nV_mem_new   )
    call reallocate( mesh%TriMap,    nTri_mem_new )
    call reallocate( mesh%TriStack1, nTri_mem_new )
    call reallocate( mesh%TriStack2, nTri_mem_new )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE extend_submesh_primary_new
  SUBROUTINE extend_mesh_primary(         mesh, nV_mem_new, nTri_mem_new)
    ! For when we didn't allocate enough. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (bigger) one, and copy the data back.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    INTEGER,                         INTENT(IN)        :: nV_mem_new, nTri_mem_new
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extend_mesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    IF (par%master) mesh%nV_mem   = nV_mem_new
    IF (par%master) mesh%nTri_mem = nTri_mem_new
       
    CALL adapt_shared_dp_2D(       mesh%nV,   nV_mem_new,    2,           mesh%V,              mesh%wV             )
    CALL adapt_shared_dp_2D(       mesh%nV,   nV_mem_new,    2,           mesh%V,              mesh%wV             )
    CALL adapt_shared_dp_2D(       mesh%nV,   nV_mem_new,    2,           mesh%V,              mesh%wV             )
    CALL adapt_shared_int_1D(      mesh%nV,   nV_mem_new,                 mesh%nC,             mesh%wnC            )
    CALL adapt_shared_int_2D(      mesh%nV,   nV_mem_new,    mesh%nC_mem, mesh%C,              mesh%wC             )
    CALL adapt_shared_int_1D(      mesh%nV,   nV_mem_new,                 mesh%niTri,          mesh%wniTri         )
    CALL adapt_shared_int_2D(      mesh%nV,   nV_mem_new,    mesh%nC_mem, mesh%iTri,           mesh%wiTri          )
    CALL adapt_shared_int_1D(      mesh%nV,   nV_mem_new,                 mesh%edge_index,     mesh%wedge_index    )
    CALL adapt_shared_int_1D(      mesh%nV,   nV_mem_new,                 mesh%mesh_old_ti_in, mesh%wmesh_old_ti_in)
    
    IF (par%master) mesh%mesh_old_ti_in(mesh%nV+1:nV_mem_new) = 1
    
    CALL adapt_shared_int_2D(      mesh%nTri, nTri_mem_new, 3,            mesh%Tri,            mesh%wTri           )
    CALL adapt_shared_dp_2D(       mesh%nTri, nTri_mem_new, 2,            mesh%Tricc,          mesh%wTricc         )
    CALL adapt_shared_int_2D(      mesh%nTri, nTri_mem_new, 3,            mesh%TriC,           mesh%wTriC          )
    CALL adapt_shared_int_1D(      mesh%nTri, nTri_mem_new,               mesh%Tri_edge_index, mesh%wTri_edge_index)
    
    CALL adapt_shared_int_2D(      mesh%nTri, nTri_mem_new, 2,            mesh%Triflip,        mesh%wTriflip       )
    CALL adapt_shared_int_1D(      mesh%nTri, nTri_mem_new,               mesh%RefMap,         mesh%wRefMap        )
    CALL adapt_shared_int_1D(      mesh%nTri, nTri_mem_new,               mesh%RefStack,       mesh%wRefStack      )
    
    ! Distributed shared memory for FloodFill maps and stacks
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                 mesh%VMap,           mesh%wVMap          )
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                 mesh%VStack1,        mesh%wVStack1       )
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                 mesh%VStack2,        mesh%wVStack2       )
    CALL adapt_shared_dist_int_1D( mesh%nTri, nTri_mem_new,               mesh%TriMap,         mesh%wTriMap        )
    CALL adapt_shared_dist_int_1D( mesh%nTri, nTri_mem_new,               mesh%TriStack1,      mesh%wTriStack1     )
    CALL adapt_shared_dist_int_1D( mesh%nTri, nTri_mem_new,               mesh%TriStack2,      mesh%wTriStack2     )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE extend_mesh_primary
  SUBROUTINE crop_mesh_primary(           mesh)
    ! For when we allocated too much. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (smaller) one, and copy the data back.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'crop_mesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
   
    IF (par%master) mesh%nV_mem   = mesh%nV
    IF (par%master) mesh%nTri_mem = mesh%nTri
       
    CALL adapt_shared_dp_2D(       mesh%nV,   mesh%nV,   2,             mesh%V,              mesh%wV             )
    CALL adapt_shared_int_1D(      mesh%nV,   mesh%nV,                  mesh%nC,             mesh%wnC            )
    CALL adapt_shared_int_2D(      mesh%nV,   mesh%nV,   mesh%nC_mem,   mesh%C,              mesh%wC             )
    CALL adapt_shared_int_1D(      mesh%nV,   mesh%nV,                  mesh%niTri,          mesh%wniTri         )
    CALL adapt_shared_int_2D(      mesh%nV,   mesh%nV,   mesh%nC_mem,   mesh%iTri,           mesh%wiTri          )
    CALL adapt_shared_int_1D(      mesh%nV,   mesh%nV,                  mesh%edge_index,     mesh%wedge_index    )
    CALL adapt_shared_int_1D(      mesh%nV,   mesh%nV,                  mesh%mesh_old_ti_in, mesh%wmesh_old_ti_in)
    
    CALL adapt_shared_int_2D(      mesh%nTri, mesh%nTri, 3,             mesh%Tri,            mesh%wTri           )
    CALL adapt_shared_dp_2D(       mesh%nTri, mesh%nTri, 2,             mesh%Tricc,          mesh%wTricc         )
    CALL adapt_shared_int_2D(      mesh%nTri, mesh%nTri, 3,             mesh%TriC,           mesh%wTriC          )
    CALL adapt_shared_int_1D(      mesh%nTri, mesh%nTri,                mesh%Tri_edge_index, mesh%wTri_edge_index)
    
    CALL adapt_shared_int_2D(      mesh%nTri, mesh%nTri, 2,             mesh%Triflip,        mesh%wTriflip       )
    CALL adapt_shared_int_1D(      mesh%nTri, mesh%nTri,                mesh%RefMap,         mesh%wRefMap        )
    CALL adapt_shared_int_1D(      mesh%nTri, mesh%nTri,                mesh%RefStack,       mesh%wRefStack      )
    
    ! Distributed shared memory for FloodFill maps and stacks
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%VMap,           mesh%wVMap          )
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%VStack1,        mesh%wVStack1       )
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%VStack2,        mesh%wVStack2       )
    CALL adapt_shared_dist_int_1D( mesh%nTri, mesh%nTri,                mesh%TriMap,         mesh%wTriMap        )
    CALL adapt_shared_dist_int_1D( mesh%nTri, mesh%nTri,                mesh%TriStack1,      mesh%wTriStack1     )
    CALL adapt_shared_dist_int_1D( mesh%nTri, mesh%nTri,                mesh%TriStack2,      mesh%wTriStack2     )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE crop_mesh_primary
  SUBROUTINE allocate_mesh_secondary(     mesh)
    ! Allocate memory for mesh data
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_mesh_secondary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    CALL allocate_shared_dp_1D(  mesh%nV  ,              mesh%A,               mesh%wA              )
    CALL allocate_shared_dp_2D(  mesh%nV  , 2,           mesh%VorGC,           mesh%wVorGC          )
    CALL allocate_shared_dp_1D(  mesh%nV  ,              mesh%R,               mesh%wR              )
    CALL allocate_shared_dp_2D(  mesh%nV  , mesh%nC_mem, mesh%Cw,              mesh%wCw             )
        
    CALL allocate_shared_dp_1D(  mesh%nTri,              mesh%TriA,            mesh%wTriA           )
    CALL allocate_shared_dp_2D(  mesh%nTri, 2,           mesh%TriGC,           mesh%wTriGC          )
    
    CALL allocate_shared_dp_1D(  mesh%nV  ,              mesh%lat,             mesh%wlat            )
    CALL allocate_shared_dp_1D(  mesh%nV  ,              mesh%lon,             mesh%wlon            )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 9)
    
  END SUBROUTINE allocate_mesh_secondary
  SUBROUTINE deallocate_mesh_all(         mesh)
    ! Deallocate memory for mesh data
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_mesh_all'
    
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Basic meta properties
    ! =====================
    
    CALL deallocate_shared( mesh%wlambda_M)
    CALL deallocate_shared( mesh%wphi_M)
    CALL deallocate_shared( mesh%walpha_stereo)
    CALL deallocate_shared( mesh%wxmin)
    CALL deallocate_shared( mesh%wxmax)
    CALL deallocate_shared( mesh%wymin)
    CALL deallocate_shared( mesh%wymax)
    CALL deallocate_shared( mesh%wtol_dist)
    CALL deallocate_shared( mesh%wnV_mem)
    CALL deallocate_shared( mesh%wnTri_mem)
    CALL deallocate_shared( mesh%wnC_mem)
    CALL deallocate_shared( mesh%wnV)
    CALL deallocate_shared( mesh%wnTri)
    CALL deallocate_shared( mesh%wperturb_dir)
    CALL deallocate_shared( mesh%walpha_min)
    CALL deallocate_shared( mesh%wdz_max_ice)
    CALL deallocate_shared( mesh%wres_max)
    CALL deallocate_shared( mesh%wres_max_margin)
    CALL deallocate_shared( mesh%wres_max_gl)
    CALL deallocate_shared( mesh%wres_max_cf)
    CALL deallocate_shared( mesh%wres_max_mountain)
    CALL deallocate_shared( mesh%wres_max_coast)
    CALL deallocate_shared( mesh%wres_min)
    CALL deallocate_shared( mesh%wresolution_min)
    CALL deallocate_shared( mesh%wresolution_max)

    ! Primary mesh data (needed for mesh creation & refinement)
    ! =========================================================
    
    CALL deallocate_shared( mesh%wV)
    CALL deallocate_shared( mesh%wnC)
    CALL deallocate_shared( mesh%wC)
    CALL deallocate_shared( mesh%wniTri)
    CALL deallocate_shared( mesh%wiTri)
    CALL deallocate_shared( mesh%wedge_index)
    CALL deallocate_shared( mesh%wmesh_old_ti_in)
    
    CALL deallocate_shared( mesh%wTri)
    CALL deallocate_shared( mesh%wTricc)
    CALL deallocate_shared( mesh%wTriC)
    CALL deallocate_shared( mesh%wTri_edge_index)
    
    CALL deallocate_shared( mesh%wTriflip)
    CALL deallocate_shared( mesh%wRefMap)
    CALL deallocate_shared( mesh%wRefStack)
    CALL deallocate_shared( mesh%wRefStackN)
    
    CALL deallocate_shared( mesh%wVMap)
    CALL deallocate_shared( mesh%wVStack1)
    CALL deallocate_shared( mesh%wVStack2)
    CALL deallocate_shared( mesh%wTriMap)
    CALL deallocate_shared( mesh%wTriStack1)
    CALL deallocate_shared( mesh%wTriStack2)
    
    CALL deallocate_shared( mesh%wnPOI               )
    CALL deallocate_shared( mesh%wPOI_coordinates    )
    CALL deallocate_shared( mesh%wPOI_XY_coordinates )
    CALL deallocate_shared( mesh%wPOI_resolutions    )
    CALL deallocate_shared( mesh%wPOI_vi             )
    CALL deallocate_shared( mesh%wPOI_w              )
    
    ! Secondary mesh data
    ! ===================

    CALL deallocate_shared( mesh%wA)
    CALL deallocate_shared( mesh%wVorGC)
    CALL deallocate_shared( mesh%wR)
    CALL deallocate_shared( mesh%wCw)
    CALL deallocate_shared( mesh%wTriGC)
    CALL deallocate_shared( mesh%wTriA)
    
    CALL deallocate_shared( mesh%wlat)
    CALL deallocate_shared( mesh%wlon)
    
    CALL deallocate_shared( mesh%wnAc)
    CALL deallocate_shared( mesh%wVAc)
    CALL deallocate_shared( mesh%wAci)
    CALL deallocate_shared( mesh%wiAci)
    CALL deallocate_shared( mesh%wedge_index_Ac)

    CALL MatDestroy( mesh%M_map_a_b             , perr)
    CALL MatDestroy( mesh%M_map_a_c             , perr)
    CALL MatDestroy( mesh%M_map_b_a             , perr)
    CALL MatDestroy( mesh%M_map_b_c             , perr)
    CALL MatDestroy( mesh%M_map_c_a             , perr)
    CALL MatDestroy( mesh%M_map_c_b             , perr)
    
    CALL MatDestroy( mesh%M_ddx_a_a             , perr)
    CALL MatDestroy( mesh%M_ddx_a_b             , perr)
    CALL MatDestroy( mesh%M_ddx_a_c             , perr)
    CALL MatDestroy( mesh%M_ddx_b_a             , perr)
    CALL MatDestroy( mesh%M_ddx_b_b             , perr)
    CALL MatDestroy( mesh%M_ddx_b_c             , perr)
    CALL MatDestroy( mesh%M_ddx_c_a             , perr)
    CALL MatDestroy( mesh%M_ddx_c_b             , perr)
    CALL MatDestroy( mesh%M_ddx_c_c             , perr)
    
    CALL MatDestroy( mesh%M_ddy_a_a             , perr)
    CALL MatDestroy( mesh%M_ddy_a_b             , perr)
    CALL MatDestroy( mesh%M_ddy_a_c             , perr)
    CALL MatDestroy( mesh%M_ddy_b_a             , perr)
    CALL MatDestroy( mesh%M_ddy_b_b             , perr)
    CALL MatDestroy( mesh%M_ddy_b_c             , perr)
    CALL MatDestroy( mesh%M_ddy_c_a             , perr)
    CALL MatDestroy( mesh%M_ddy_c_b             , perr)
    CALL MatDestroy( mesh%M_ddy_c_c             , perr)
    
    CALL MatDestroy( mesh%M2_ddx_b_b            , perr)
    CALL MatDestroy( mesh%M2_ddy_b_b            , perr)
    CALL MatDestroy( mesh%M2_d2dx2_b_b          , perr)
    CALL MatDestroy( mesh%M2_d2dxdy_b_b         , perr)
    CALL MatDestroy( mesh%M2_d2dy2_b_b          , perr)
    
    CALL deallocate_matrix_CSR( mesh%M2_ddx_b_b_CSR    )
    CALL deallocate_matrix_CSR( mesh%M2_ddy_b_b_CSR    )
    CALL deallocate_matrix_CSR( mesh%M2_d2dx2_b_b_CSR  )
    CALL deallocate_matrix_CSR( mesh%M2_d2dxdy_b_b_CSR )
    CALL deallocate_matrix_CSR( mesh%M2_d2dy2_b_b_CSR  )
    
    CALL MatDestroy( mesh%M_Neumann_BC_b        , perr)
    CALL deallocate_matrix_CSR( mesh%M_Neumann_BC_b_CSR)
    
    CALL deallocate_shared( mesh%wnV_transect)
    CALL deallocate_shared( mesh%wvi_transect)
    CALL deallocate_shared( mesh%ww_transect )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
 
   END SUBROUTINE deallocate_mesh_all
  
  SUBROUTINE allocate_submesh_primary(    mesh, region_name, nV_mem, nTri_mem, nC_mem)
    ! Allocate memory for mesh data
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    CHARACTER(LEN=3),                INTENT(IN)        :: region_name
    INTEGER,                         INTENT(IN)        :: nV_mem, nTri_mem, nC_mem
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_submesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    mesh%region_name = region_name
    
    CALL allocate_shared_dist_dp_0D(  mesh%lambda_M,         mesh%wlambda_M        )
    CALL allocate_shared_dist_dp_0D(  mesh%phi_M,            mesh%wphi_M           )
    CALL allocate_shared_dist_dp_0D(  mesh%alpha_stereo,     mesh%walpha_stereo    )
    
    IF (mesh%region_name     == 'NAM') THEN
      mesh%lambda_M     = C%lambda_M_NAM
      mesh%phi_M        = C%phi_M_NAM
      mesh%alpha_stereo = C%alpha_stereo_NAM
    ELSEIF (mesh%region_name == 'EAS') THEN
      mesh%lambda_M     = C%lambda_M_EAS
      mesh%phi_M        = C%phi_M_EAS
      mesh%alpha_stereo = C%alpha_stereo_EAS
    ELSEIF (mesh%region_name == 'GRL') THEN
      mesh%lambda_M     = C%lambda_M_GRL
      mesh%phi_M        = C%phi_M_GRL
      mesh%alpha_stereo = C%alpha_stereo_GRL
    ELSEIF (mesh%region_name == 'ANT') THEN
      mesh%lambda_M     = C%lambda_M_ANT
      mesh%phi_M        = C%phi_M_ANT
      mesh%alpha_stereo = C%alpha_stereo_ANT
    END IF
    
    CALL allocate_shared_dist_dp_0D(  mesh%xmin,             mesh%wxmin            )
    CALL allocate_shared_dist_dp_0D(  mesh%xmax,             mesh%wxmax            )
    CALL allocate_shared_dist_dp_0D(  mesh%ymin,             mesh%wymin            )
    CALL allocate_shared_dist_dp_0D(  mesh%ymax,             mesh%wymax            )
    CALL allocate_shared_dist_dp_0D(  mesh%tol_dist,         mesh%wtol_dist        )
    CALL allocate_shared_dist_int_0D( mesh%nC_mem,           mesh%wnC_mem          )
    CALL allocate_shared_dist_int_0D( mesh%nV_mem,           mesh%wnV_mem          )
    CALL allocate_shared_dist_int_0D( mesh%nTri_mem,         mesh%wnTri_mem        )
    CALL allocate_shared_dist_int_0D( mesh%nV,               mesh%wnV              )
    CALL allocate_shared_dist_int_0D( mesh%nTri,             mesh%wnTri            )
    CALL allocate_shared_dist_int_0D( mesh%perturb_dir,      mesh%wperturb_dir     )
    CALL allocate_shared_dist_dp_0D(  mesh%alpha_min,        mesh%walpha_min       )
    CALL allocate_shared_dist_dp_0D(  mesh%dz_max_ice,       mesh%wdz_max_ice      )
    CALL allocate_shared_dist_dp_0D(  mesh%res_max,          mesh%wres_max         )
    CALL allocate_shared_dist_dp_0D(  mesh%res_max_margin,   mesh%wres_max_margin  )
    CALL allocate_shared_dist_dp_0D(  mesh%res_max_gl,       mesh%wres_max_gl      )
    CALL allocate_shared_dist_dp_0D(  mesh%res_max_cf,       mesh%wres_max_cf      )
    CALL allocate_shared_dist_dp_0D(  mesh%res_max_mountain, mesh%wres_max_mountain)
    CALL allocate_shared_dist_dp_0D(  mesh%res_max_coast,    mesh%wres_max_coast   )
    CALL allocate_shared_dist_dp_0D(  mesh%res_min,          mesh%wres_min         )
    CALL allocate_shared_dist_dp_0D(  mesh%resolution_min,   mesh%wresolution_min  )
    CALL allocate_shared_dist_dp_0D(  mesh%resolution_max,   mesh%wresolution_max  )
    
    mesh%nV_mem   = nV_mem
    mesh%nTri_mem = nTri_mem
    mesh%nC_mem   = nC_mem
    
    CALL allocate_shared_dist_dp_2D(  nV_mem,   2,      mesh%V,               mesh%wV              )
    CALL allocate_shared_dist_int_1D( nV_mem,           mesh%nC,              mesh%wnC             )
    CALL allocate_shared_dist_int_2D( nV_mem,   nC_mem, mesh%C,               mesh%wC              )
    CALL allocate_shared_dist_int_1D( nV_mem,           mesh%niTri,           mesh%wniTri          )
    CALL allocate_shared_dist_int_2D( nV_mem,   nC_mem, mesh%iTri,            mesh%wiTri           )
    CALL allocate_shared_dist_int_1D( nV_mem,           mesh%edge_index,      mesh%wedge_index     )
    CALL allocate_shared_dist_int_1D( nV_mem,           mesh%mesh_old_ti_in,  mesh%wmesh_old_ti_in )
    
    CALL allocate_shared_dist_int_2D( nTri_mem, 3,      mesh%Tri,             mesh%wTri            )
    CALL allocate_shared_dist_dp_2D(  nTri_mem, 2,      mesh%Tricc,           mesh%wTricc          )
    CALL allocate_shared_dist_int_2D( nTri_mem, 3,      mesh%TriC,            mesh%wTriC           )
    CALL allocate_shared_dist_int_1D( nTri_mem,         mesh%Tri_edge_index,  mesh%wTri_edge_index )
    
    CALL allocate_shared_dist_int_2D( nTri_mem, 2,      mesh%Triflip,         mesh%wTriflip        )
    CALL allocate_shared_dist_int_1D( nTri_mem,         mesh%RefMap,          mesh%wRefMap         )
    CALL allocate_shared_dist_int_1D( nTri_mem,         mesh%RefStack,        mesh%wRefStack       )
    CALL allocate_shared_dist_int_0D(                   mesh%RefStackN,       mesh%wRefStackN      )
    
    ! Distributed shared memory for FloodFill maps and stacks
    CALL allocate_shared_dist_int_1D( nV_mem,           mesh%VMap,            mesh%wVMap           )
    CALL allocate_shared_dist_int_1D( nV_mem,           mesh%VStack1,         mesh%wVStack1        )
    CALL allocate_shared_dist_int_1D( nV_mem,           mesh%VStack2,         mesh%wVStack2        )
    CALL allocate_shared_dist_int_1D( nTri_mem,         mesh%TriMap,          mesh%wTriMap         )
    CALL allocate_shared_dist_int_1D( nTri_mem,         mesh%TriStack1,       mesh%wTriStack1      )
    CALL allocate_shared_dist_int_1D( nTri_mem,         mesh%TriStack2,       mesh%wTriStack2      )
    
    ! POI stuff
    CALL allocate_shared_int_0D(  mesh%nPOI, mesh%wnPOI)
    
    IF (mesh%region_name == 'NAM') THEN
      mesh%nPOI = C%nPOI_NAM
    ELSEIF (mesh%region_name == 'EAS') THEN
      mesh%nPOI = C%nPOI_EAS
    ELSEIF (mesh%region_name == 'GRL') THEN
      mesh%nPOI = C%nPOI_GRL
    ELSEIF (mesh%region_name == 'ANT') THEN
      mesh%nPOI = C%nPOI_ANT
    END IF
    
    CALL allocate_shared_dp_2D(  mesh%nPOI, 2, mesh%POI_coordinates,           mesh%wPOI_coordinates          )
    CALL allocate_shared_dp_2D(  mesh%nPOI, 2, mesh%POI_XY_coordinates,        mesh%wPOI_XY_coordinates       )
    CALL allocate_shared_dp_1D(  mesh%nPOI,    mesh%POI_resolutions,           mesh%wPOI_resolutions          )
    CALL allocate_shared_int_2D( mesh%nPOI, 3, mesh%POI_vi,                    mesh%wPOI_vi                   )
    CALL allocate_shared_dp_2D(  mesh%nPOI, 3, mesh%POI_w,                     mesh%wPOI_w                    )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 52)
    
  END SUBROUTINE allocate_submesh_primary
  SUBROUTINE allocate_submesh_primary_new(    mesh, region_name, nV_mem, nTri_mem, nC_mem)
    ! Allocate memory for mesh data
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh_new),                 INTENT(INOUT)     :: mesh
    CHARACTER(LEN=3),                INTENT(IN)        :: region_name
    INTEGER,                         INTENT(IN)        :: nV_mem, nTri_mem, nC_mem
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_submesh_primary'
    
    ! Add routine to path
    mesh%region_name = region_name
    
    IF (mesh%region_name     == 'NAM') THEN
      mesh%lambda_M     = C%lambda_M_NAM
      mesh%phi_M        = C%phi_M_NAM
      mesh%alpha_stereo = C%alpha_stereo_NAM
    ELSEIF (mesh%region_name == 'EAS') THEN
      mesh%lambda_M     = C%lambda_M_EAS
      mesh%phi_M        = C%phi_M_EAS
      mesh%alpha_stereo = C%alpha_stereo_EAS
    ELSEIF (mesh%region_name == 'GRL') THEN
      mesh%lambda_M     = C%lambda_M_GRL
      mesh%phi_M        = C%phi_M_GRL
      mesh%alpha_stereo = C%alpha_stereo_GRL
    ELSEIF (mesh%region_name == 'ANT') THEN
      mesh%lambda_M     = C%lambda_M_ANT
      mesh%phi_M        = C%phi_M_ANT
      mesh%alpha_stereo = C%alpha_stereo_ANT
    END IF
    
    mesh%nV_mem   = nV_mem
    mesh%nTri_mem = nTri_mem
    mesh%nC_mem   = nC_mem
    allocate( mesh%V              (nV_mem,   2     ))
    allocate( mesh%nC             (nV_mem          ))
    allocate( mesh%C              (nV_mem,   nC_mem))
    allocate( mesh%niTri          (nV_mem          ))
    allocate( mesh%iTri           (nV_mem,   nC_mem))
    allocate( mesh%edge_index     (nV_mem          ))
    allocate( mesh%mesh_old_ti_in (nV_mem          ))

    allocate( mesh%Tri            (nTri_mem, 3     ))
    allocate( mesh%Tricc          (nTri_mem, 2     ))
    allocate( mesh%TriC           (nTri_mem, 3     )) 
    allocate( mesh%Tri_edge_index (nTri_mem        )) 
                                                     
    allocate( mesh%Triflip        (nTri_mem, 2     ))
    allocate( mesh%RefMap         (nTri_mem        ))
    allocate( mesh%RefStack       (nTri_mem        ))
    
    ! Distributed shared memory for FloodFill maps and stacks
    allocate( mesh%VMap           (nV_mem          ))
    allocate( mesh%VStack1        (nV_mem          ))
    allocate( mesh%VStack2        (nV_mem          ))
    allocate( mesh%TriMap         (nTri_mem        ))
    allocate( mesh%TriStack1      (nTri_mem        ))
    allocate( mesh%TriStack2      (nTri_mem        ))
    
    ! POI stuff
    IF (mesh%region_name == 'NAM') THEN
      mesh%nPOI = C%nPOI_NAM
    ELSEIF (mesh%region_name == 'EAS') THEN
      mesh%nPOI = C%nPOI_EAS
    ELSEIF (mesh%region_name == 'GRL') THEN
      mesh%nPOI = C%nPOI_GRL
    ELSEIF (mesh%region_name == 'ANT') THEN
      mesh%nPOI = C%nPOI_ANT
    END IF
    
    allocate( mesh%POI_coordinates    ( mesh%nPOI, 2 )) 
    allocate( mesh%POI_XY_coordinates ( mesh%nPOI, 2 )) 
    allocate( mesh%POI_resolutions    ( mesh%nPOI    ))
    allocate( mesh%POI_vi             ( mesh%nPOI, 3 )) 
    allocate( mesh%POI_w              ( mesh%nPOI, 3 )) 

  END SUBROUTINE allocate_submesh_primary_new
  SUBROUTINE extend_submesh_primary(      mesh, nV_mem_new, nTri_mem_new)
    ! For when we didn't allocate enough. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (bigger) one, and copy the data back.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    INTEGER,                         INTENT(IN)        :: nV_mem_new, nTri_mem_new
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extend_submesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
 
    mesh%nV_mem   = nV_mem_new
    mesh%nTri_mem = nTri_mem_new
       
    CALL adapt_shared_dist_dp_2D(  mesh%nV,   nV_mem_new,   2,             mesh%V,              mesh%wV             )
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                  mesh%nC,             mesh%wnC            )
    CALL adapt_shared_dist_int_2D( mesh%nV,   nV_mem_new,   mesh%nC_mem,   mesh%C,              mesh%wC             )
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                  mesh%niTri,          mesh%wniTri         )
    CALL adapt_shared_dist_int_2D( mesh%nV,   nV_mem_new,   mesh%nC_mem,   mesh%iTri,           mesh%wiTri          )
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                  mesh%edge_index,     mesh%wedge_index    )
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                  mesh%mesh_old_ti_in, mesh%wmesh_old_ti_in)
    
    mesh%mesh_old_ti_in(mesh%nV+1:nV_mem_new) = 1
    
    CALL adapt_shared_dist_int_2D( mesh%nTri, nTri_mem_new, 3,             mesh%Tri,            mesh%wTri           )
    CALL adapt_shared_dist_dp_2D(  mesh%nTri, nTri_mem_new, 2,             mesh%Tricc,          mesh%wTricc         )
    CALL adapt_shared_dist_int_2D( mesh%nTri, nTri_mem_new, 3,             mesh%TriC,           mesh%wTriC          )
    CALL adapt_shared_dist_int_1D( mesh%nTri, nTri_mem_new,                mesh%Tri_edge_index, mesh%wTri_edge_index)
    
    CALL adapt_shared_dist_int_2D( mesh%nTri, nTri_mem_new, 2,             mesh%Triflip,        mesh%wTriflip       )
    CALL adapt_shared_dist_int_1D( mesh%nTri, nTri_mem_new,                mesh%RefMap,         mesh%wRefMap        )
    CALL adapt_shared_dist_int_1D( mesh%nTri, nTri_mem_new,                mesh%RefStack,       mesh%wRefStack      )
    
    ! Distributed shared memory for FloodFill maps and stacks
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                  mesh%VMap,           mesh%wVMap          )
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                  mesh%VStack1,        mesh%wVStack1       )
    CALL adapt_shared_dist_int_1D( mesh%nV,   nV_mem_new,                  mesh%VStack2,        mesh%wVStack2       )
    CALL adapt_shared_dist_int_1D( mesh%nTri, nTri_mem_new,                mesh%TriMap,         mesh%wTriMap        )
    CALL adapt_shared_dist_int_1D( mesh%nTri, nTri_mem_new,                mesh%TriStack1,      mesh%wTriStack1     )
    CALL adapt_shared_dist_int_1D( mesh%nTri, nTri_mem_new,                mesh%TriStack2,      mesh%wTriStack2     )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE extend_submesh_primary
  SUBROUTINE crop_submesh_primary_new(        mesh)
    ! For when we allocated too much. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (smaller) one, and copy the data back.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh_new),                 INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'crop_submesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
   
    mesh%nV_mem = mesh%nV
    
    call reallocate( mesh%V             , mesh%nV,   2          )
    call reallocate( mesh%nC            , mesh%nV               )
    call reallocate( mesh%C             , mesh%nV,   mesh%nC_mem)
    call reallocate( mesh%niTri         , mesh%nV               )
    call reallocate( mesh%iTri          , mesh%nV,   mesh%nC_mem)
    call reallocate( mesh%edge_index    , mesh%nV               )
    call reallocate( mesh%mesh_old_ti_in, mesh%nV               )

    call reallocate( mesh%Tri           , mesh%nTri, 3          )
    call reallocate( mesh%Tricc         , mesh%nTri, 2          )
    call reallocate( mesh%TriC          , mesh%nTri, 3          )
    call reallocate( mesh%Tri_edge_index, mesh%nTri             )

    call reallocate( mesh%Triflip       , mesh%nTri, 2          )
    call reallocate( mesh%RefMap        , mesh%nTri             )
    call reallocate( mesh%RefStack      , mesh%nTri             )

    ! Distributed shared memory for,FloodFill maps and stacks
    call reallocate( mesh%VMap          , mesh%nV  )
    call reallocate( mesh%VStack1       , mesh%nV  )
    call reallocate( mesh%VStack2       , mesh%nV  )
    call reallocate( mesh%TriMap        , mesh%nTri)
    call reallocate( mesh%TriStack1     , mesh%nTri)
    call reallocate( mesh%TriStack2     , mesh%nTri)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE crop_submesh_primary_new
  SUBROUTINE crop_submesh_primary(        mesh)
    ! For when we allocated too much. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (smaller) one, and copy the data back.
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'crop_submesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
   
    mesh%nV_mem = mesh%nV
       
    CALL adapt_shared_dist_dp_2D(  mesh%nV,   mesh%nV,   2,             mesh%V,              mesh%wV             )
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%nC,             mesh%wnC            )
    CALL adapt_shared_dist_int_2D( mesh%nV,   mesh%nV,   mesh%nC_mem,   mesh%C,              mesh%wC             )
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%niTri,          mesh%wniTri         )
    CALL adapt_shared_dist_int_2D( mesh%nV,   mesh%nV,   mesh%nC_mem,   mesh%iTri,           mesh%wiTri          )
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%edge_index,     mesh%wedge_index    )
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%mesh_old_ti_in, mesh%wmesh_old_ti_in)
    
    CALL adapt_shared_dist_int_2D( mesh%nTri, mesh%nTri, 3,             mesh%Tri,            mesh%wTri           )
    CALL adapt_shared_dist_dp_2D(  mesh%nTri, mesh%nTri, 2,             mesh%Tricc,          mesh%wTricc         )
    CALL adapt_shared_dist_int_2D( mesh%nTri, mesh%nTri, 3,             mesh%TriC,           mesh%wTriC          )
    CALL adapt_shared_dist_int_1D( mesh%nTri, mesh%nTri,                mesh%Tri_edge_index, mesh%wTri_edge_index)
    
    CALL adapt_shared_dist_int_2D( mesh%nTri, mesh%nTri, 2,             mesh%Triflip,        mesh%wTriflip       )
    CALL adapt_shared_dist_int_1D( mesh%nTri, mesh%nTri,                mesh%RefMap,         mesh%wRefMap        )
    CALL adapt_shared_dist_int_1D( mesh%nTri, mesh%nTri,                mesh%RefStack,       mesh%wRefStack      )
    
    ! Distributed shared memory for FloodFill maps and stacks
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%VMap,           mesh%wVMap          )
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%VStack1,        mesh%wVStack1       )
    CALL adapt_shared_dist_int_1D( mesh%nV,   mesh%nV,                  mesh%VStack2,        mesh%wVStack2       )
    CALL adapt_shared_dist_int_1D( mesh%nTri, mesh%nTri,                mesh%TriMap,         mesh%wTriMap        )
    CALL adapt_shared_dist_int_1D( mesh%nTri, mesh%nTri,                mesh%TriStack1,      mesh%wTriStack1     )
    CALL adapt_shared_dist_int_1D( mesh%nTri, mesh%nTri,                mesh%TriStack2,      mesh%wTriStack2     )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE crop_submesh_primary
  SUBROUTINE deallocate_submesh_primary_new(  mesh)
    ! Deallocate memory for mesh data
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh_new),                 INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_submesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    deallocate(mesh%V)
    deallocate(mesh%nC)
    deallocate(mesh%C)
    deallocate(mesh%niTri)
    deallocate(mesh%iTri)
    deallocate(mesh%edge_index)
    deallocate(mesh%mesh_old_ti_in)
    
    deallocate(mesh%Tri)
    deallocate(mesh%Tricc)
    deallocate(mesh%TriC)
    deallocate(mesh%Tri_edge_index)

    deallocate(mesh%Triflip)
    deallocate(mesh%RefMap)
    deallocate(mesh%RefStack)
    deallocate(mesh%RefStackN)

    deallocate(mesh%VMap)
    deallocate(mesh%VStack1)
    deallocate(mesh%VStack2)
    deallocate(mesh%TriMap)
    deallocate(mesh%TriStack1)
    deallocate(mesh%TriStack2)

    deallocate(mesh%nPOI               )
    deallocate(mesh%POI_coordinates    )
    deallocate(mesh%POI_XY_coordinates )
    deallocate(mesh%POI_resolutions    )
    deallocate(mesh%POI_vi             )
    deallocate(mesh%POI_w              )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
 
   END SUBROUTINE deallocate_submesh_primary_new
  SUBROUTINE deallocate_submesh_primary(  mesh)
    ! Deallocate memory for mesh data
    
    IMPLICIT NONE
    
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'deallocate_submesh_primary'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    CALL deallocate_shared( mesh%wlambda_M)
    CALL deallocate_shared( mesh%wphi_M)
    CALL deallocate_shared( mesh%walpha_stereo)
    CALL deallocate_shared( mesh%wxmin)
    CALL deallocate_shared( mesh%wxmax)
    CALL deallocate_shared( mesh%wymin)
    CALL deallocate_shared( mesh%wymax)
    CALL deallocate_shared( mesh%wtol_dist)
    CALL deallocate_shared( mesh%wnC_mem)
    CALL deallocate_shared( mesh%wnV_mem)
    CALL deallocate_shared( mesh%wnTri_mem)
    CALL deallocate_shared( mesh%wnV)
    CALL deallocate_shared( mesh%wnTri)
    CALL deallocate_shared( mesh%wperturb_dir)
    CALL deallocate_shared( mesh%walpha_min)
    CALL deallocate_shared( mesh%wdz_max_ice)
    CALL deallocate_shared( mesh%wres_max)
    CALL deallocate_shared( mesh%wres_max_margin)
    CALL deallocate_shared( mesh%wres_max_gl)
    CALL deallocate_shared( mesh%wres_max_cf)
    CALL deallocate_shared( mesh%wres_max_mountain)
    CALL deallocate_shared( mesh%wres_max_coast)
    CALL deallocate_shared( mesh%wres_min)
    CALL deallocate_shared( mesh%wresolution_min)
    CALL deallocate_shared( mesh%wresolution_max)
    
    CALL deallocate_shared( mesh%wV)
    CALL deallocate_shared( mesh%wnC)
    CALL deallocate_shared( mesh%wC)
    CALL deallocate_shared( mesh%wniTri)
    CALL deallocate_shared( mesh%wiTri)
    CALL deallocate_shared( mesh%wedge_index)
    CALL deallocate_shared( mesh%wmesh_old_ti_in)
    
    CALL deallocate_shared( mesh%wTri)
    CALL deallocate_shared( mesh%wTricc)
    CALL deallocate_shared( mesh%wTriC)
    CALL deallocate_shared( mesh%wTri_edge_index)
    
    CALL deallocate_shared( mesh%wTriflip)
    CALL deallocate_shared( mesh%wRefMap)
    CALL deallocate_shared( mesh%wRefStack)
    CALL deallocate_shared( mesh%wRefStackN)
    
    CALL deallocate_shared( mesh%wVMap)
    CALL deallocate_shared( mesh%wVStack1)
    CALL deallocate_shared( mesh%wVStack2)
    CALL deallocate_shared( mesh%wTriMap)
    CALL deallocate_shared( mesh%wTriStack1)
    CALL deallocate_shared( mesh%wTriStack2)
    
    CALL deallocate_shared( mesh%wnPOI               )
    CALL deallocate_shared( mesh%wPOI_coordinates    )
    CALL deallocate_shared( mesh%wPOI_XY_coordinates )
    CALL deallocate_shared( mesh%wPOI_resolutions    )
    CALL deallocate_shared( mesh%wPOI_vi             )
    CALL deallocate_shared( mesh%wPOI_w              )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
 
   END SUBROUTINE deallocate_submesh_primary
  SUBROUTINE share_submesh_access( p_left, p_right, submesh, submesh_right)
    ! Give process p_left access to the submesh memory of p_right
    ! Used in submesh merging
    
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_LOC
    
    IMPLICIT NONE
 
    ! In/output variables:
    INTEGER,                         INTENT(IN)        :: p_left
    INTEGER,                         INTENT(IN)        :: p_right
    TYPE(type_mesh),                 INTENT(IN)        :: submesh
    TYPE(type_mesh),                 INTENT(INOUT)     :: submesh_right
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'share_submesh_access'
    INTEGER                                            :: nV, nTri, nconmax
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    
    CALL share_memory_access_dp_0D(  p_left, p_right, submesh_right%xmin,             submesh%wxmin,             submesh_right%wxmin            )
    CALL share_memory_access_dp_0D(  p_left, p_right, submesh_right%xmax,             submesh%wxmax,             submesh_right%wxmax            )
    CALL share_memory_access_dp_0D(  p_left, p_right, submesh_right%ymin,             submesh%wymin,             submesh_right%wymin            )
    CALL share_memory_access_dp_0D(  p_left, p_right, submesh_right%ymax,             submesh%wymax,             submesh_right%wymax            )
    CALL share_memory_access_int_0D( p_left, p_right, submesh_right%nC_mem,           submesh%wnC_mem,           submesh_right%wnC_mem          )
    CALL share_memory_access_int_0D( p_left, p_right, submesh_right%nV_mem,           submesh%wnV_mem,           submesh_right%wnV_mem          )
    CALL share_memory_access_int_0D( p_left, p_right, submesh_right%nTri_mem,         submesh%wnTri_mem,         submesh_right%wnTri_mem        )
    CALL share_memory_access_int_0D( p_left, p_right, submesh_right%nV,               submesh%wnV,               submesh_right%wnV              )
    CALL share_memory_access_int_0D( p_left, p_right, submesh_right%nTri,             submesh%wnTri,             submesh_right%wnTri            )
    
    IF (par%i == p_left) THEN
      nV      = submesh_right%nV_mem
      nTri    = submesh_right%nTri_mem
      nconmax = submesh_right%nC_mem
    ELSEIF (par%i == p_right) THEN
      nV      = submesh%nV_mem
      nTri    = submesh%nTri_mem
      nconmax = submesh%nC_mem
    END IF
    
    CALL share_memory_access_dp_2D(  p_left, p_right, submesh_right%V,                submesh%wV,                submesh_right%wV,                nV,   2      )
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%nC,               submesh%wnC,               submesh_right%wnC,               nV           )
    CALL share_memory_access_int_2D( p_left, p_right, submesh_right%C,                submesh%wC,                submesh_right%wC,                nV,   nconmax)
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%niTri,            submesh%wniTri,            submesh_right%wniTri,            nV           )
    CALL share_memory_access_int_2D( p_left, p_right, submesh_right%iTri,             submesh%witri,             submesh_right%wiTri,             nV,   nconmax)
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%edge_index,       submesh%wedge_index,       submesh_right%wedge_index,       nV           )
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%mesh_old_ti_in,   submesh%wmesh_old_ti_in,   submesh_right%wmesh_old_ti_in,   nV           )
    
    CALL share_memory_access_int_2D( p_left, p_right, submesh_right%Tri,              submesh%wTri,              submesh_right%wTri,              nTri, 3      )
    CALL share_memory_access_dp_2D(  p_left, p_right, submesh_right%Tricc,            submesh%wTricc,            submesh_right%wTricc,            nTri, 2      )
    CALL share_memory_access_int_2D( p_left, p_right, submesh_right%TriC,             submesh%wTriC,             submesh_right%wTriC,             nTri, 3      )
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%Tri_edge_index,   submesh%wTri_edge_index,   submesh_right%wTri_edge_index,   nTri         )
        
    CALL share_memory_access_int_2D( p_left, p_right, submesh_right%Triflip,          submesh%wTriflip,          submesh_right%wTriflip,          nTri, 2      )
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%RefMap,           submesh%wRefMap,           submesh_right%wRefMap,           nTri         )
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%RefStack,         submesh%wRefStack,         submesh_right%wRefStack,         nTri         )
    CALL share_memory_access_int_0D( p_left, p_right, submesh_right%RefStackN,        submesh%wRefStackN,        submesh_right%wRefStackN                      )
    
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%VMap,             submesh%wVMap,             submesh_right%wVMap,             nV           )
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%VStack1,          submesh%wVStack1,          submesh_right%wVStack1,          nV           )
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%VStack2,          submesh%wVStack2,          submesh_right%wVStack2,          nV           )
    
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%TriMap,           submesh%wTriMap,           submesh_right%wTriMap,           nTri         )
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%TriStack1,        submesh%wTriStack1,        submesh_right%wTriStack1,        nTri         )
    CALL share_memory_access_int_1D( p_left, p_right, submesh_right%TriStack2,        submesh%wTriStack2,        submesh_right%wTriStack2,        nTri         )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE share_submesh_access
  SUBROUTINE share_submesh_access_new( p_left, p_right, submesh, submesh_right)
    ! Give process p_left access to the submesh memory of p_right
    ! Used in submesh merging
    
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_F_POINTER, C_LOC
    
    IMPLICIT NONE
 
    ! In/output variables:
    INTEGER,                         INTENT(IN)        :: p_left
    INTEGER,                         INTENT(IN)        :: p_right
    INTEGER                                       :: status(MPI_STATUS_SIZE)
    TYPE(type_mesh_new),                 INTENT(IN)        :: submesh
    TYPE(type_mesh_new),                 INTENT(INOUT)     :: submesh_right
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'share_submesh_access'
    INTEGER                                            :: nV, nTri, nconmax
    
    ! Add routine to path
    CALL init_routine( routine_name)

    ! Left == receiving
    if (par%i == p_left) then
      ! these could be combined ...
      call mpi_recv( submesh_right%xmin,    1, MPI_REAL8  , p_right, 1, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%xmax,    1, MPI_REAL8  , p_right, 2, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%ymin,    1, MPI_REAL8  , p_right, 3, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%ymax,    1, MPI_REAL8  , p_right, 4, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%nC_mem,  1, MPI_INTEGER, p_right, 5, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%nV_mem,  1, MPI_INTEGER, p_right, 6, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%nTri_mem,1, MPI_INTEGER, p_right, 7, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%nV,      1, MPI_INTEGER, p_right, 8, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%nTri,    1, MPI_INTEGER, p_right, 9, MPI_COMM_WORLD, status, ierr )
    else
    ! Sending 
      ! these could be combined ...
      call mpi_send( submesh      %xmin,    1, MPI_REAL8  , p_left , 1, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %xmax,    1, MPI_REAL8  , p_left , 2, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %ymin,    1, MPI_REAL8  , p_left , 3, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %ymax,    1, MPI_REAL8  , p_left , 4, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %nC_mem,  1, MPI_INTEGER, p_left , 5, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %nV_mem,  1, MPI_INTEGER, p_left , 6, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %nTri_mem,1, MPI_INTEGER, p_left , 7, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %nV,      1, MPI_INTEGER, p_left , 8, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %nTri,    1, MPI_INTEGER, p_left , 9, MPI_COMM_WORLD,         ierr )
    end if

    IF (par%i == p_left) THEN
    ! Sending
      nV      = submesh_right%nV_mem
      nTri    = submesh_right%nTri_mem
      nconmax = submesh_right%nC_mem
    ELSE
    ! Receiving
      nV      = submesh%nV_mem
      nTri    = submesh%nTri_mem
      nconmax = submesh%nC_mem
    END IF

    if (par%i == p_left) then
      !MAKE ROOM
      allocate( submesh_right%V             (   nV,   2      ))
      allocate( submesh_right%nC            (   nV           ))
      allocate( submesh_right%C             (   nV,   nconmax))
      allocate( submesh_right%niTri         (   nV           ))
      allocate( submesh_right%iTri          (   nV,   nconmax))
      allocate( submesh_right%edge_index    (   nV           ))
      allocate( submesh_right%mesh_old_ti_in(   nV           ))

      allocate( submesh_right%Tri           (   nTri, 3      ))
      allocate( submesh_right%Tricc         (   nTri, 2      ))
      allocate( submesh_right%TriC          (   nTri, 3      ))
      allocate( submesh_right%Tri_edge_index(   nTri         ))

      allocate( submesh_right%Triflip       (   nTri, 2      ))
      allocate( submesh_right%RefMap        (   nTri         ))
      allocate( submesh_right%RefStack      (   nTri         ))

      allocate( submesh_right%VMap          (   nV           ))
      allocate( submesh_right%VStack1       (   nV           ))
      allocate( submesh_right%VStack2       (   nV           ))

      allocate( submesh_right%TriMap        (   nTri         ))
      allocate( submesh_right%TriStack1     (   nTri         ))
      allocate( submesh_right%TriStack2     (   nTri         ))
      
      ! ok now we can receive
      call mpi_recv( submesh_right%V,              nV*2      , mpi_real8  , p_right, 1, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%nC,             nV        , mpi_integer, p_right, 2, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%C,              nV*nconmax, mpi_integer, p_right, 3, MPI_COMM_WORLD, status, ierr ) 
      call mpi_recv( submesh_right%niTri,          nV        , mpi_integer, p_right, 4, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%iTri,           nV*nconmax, mpi_integer, p_right, 5, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%edge_index,     nV        , mpi_integer, p_right, 6, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%mesh_old_ti_in, nV        , mpi_integer, p_right, 7, MPI_COMM_WORLD, status, ierr )

      call mpi_recv( submesh_right%Tri,            nTri*3    , mpi_integer, p_right, 9, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%Tricc,          nTri*2    , mpi_real8  , p_right,10, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%TriC,           nTri*3    , mpi_integer, p_right,11, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%Tri_edge_index, nTri      , mpi_integer, p_right,12, MPI_COMM_WORLD, status, ierr )

      call mpi_recv( submesh_right%Triflip,        nTri*2    , mpi_integer, p_right,13, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%RefMap,         nTri      , mpi_integer, p_right,14, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%RefStack,       nTri      , mpi_integer, p_right,15, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%RefStackN,      1         , mpi_integer, p_right,16, MPI_COMM_WORLD, status, ierr )

      call mpi_recv( submesh_right%VMap,           nV        , mpi_integer, p_right,17, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%VStack1,        nV        , mpi_integer, p_right,18, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%VStack2,        nV        , mpi_integer, p_right,19, MPI_COMM_WORLD, status, ierr )

      call mpi_recv( submesh_right%TriMap,         nTri      , mpi_integer, p_right,20, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%TriStack1,      nTri      , mpi_integer, p_right,21, MPI_COMM_WORLD, status, ierr )
      call mpi_recv( submesh_right%TriStack2,      nTri      , mpi_integer, p_right,22, MPI_COMM_WORLD, status, ierr )
    else 
   ! the right must send
      call mpi_send( submesh      %V,              nV*2      , mpi_real8  , p_left , 1, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %nC,             nV        , mpi_integer, p_left , 2, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %C,              nV*nconmax, mpi_integer, p_left , 3, MPI_COMM_WORLD,         ierr ) 
      call mpi_send( submesh      %niTri,          nV        , mpi_integer, p_left , 4, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %iTri,           nV*nconmax, mpi_integer, p_left , 5, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %edge_index,     nV        , mpi_integer, p_left , 6, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %mesh_old_ti_in, nV        , mpi_integer, p_left , 7, MPI_COMM_WORLD,         ierr )

      call mpi_send( submesh      %Tri,            nTri*3    , mpi_integer, p_left , 9, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %Tricc,          nTri*2    , mpi_real8  , p_left ,10, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %TriC,           nTri*3    , mpi_integer, p_left ,11, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %Tri_edge_index, nTri      , mpi_integer, p_left ,12, MPI_COMM_WORLD,         ierr )

      call mpi_send( submesh      %Triflip,        nTri*2    , mpi_integer, p_left ,13, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %RefMap,         nTri      , mpi_integer, p_left ,14, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %RefStack,       nTri      , mpi_integer, p_left ,15, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %RefStackN,      1         , mpi_integer, p_left ,16, MPI_COMM_WORLD,         ierr )

      call mpi_send( submesh      %VMap,           nV        , mpi_integer, p_left ,17, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %VStack1,        nV        , mpi_integer, p_left ,18, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %VStack2,        nV        , mpi_integer, p_left ,19, MPI_COMM_WORLD,         ierr )

      call mpi_send( submesh      %TriMap,         nTri      , mpi_integer, p_left ,20, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %TriStack1,      nTri      , mpi_integer, p_left ,21, MPI_COMM_WORLD,         ierr )
      call mpi_send( submesh      %TriStack2,      nTri      , mpi_integer, p_left ,22, MPI_COMM_WORLD,         ierr )
    endif
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE share_submesh_access_new
  
  SUBROUTINE move_data_from_submesh_to_mesh( mesh, submesh)
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    TYPE(type_mesh),                 INTENT(IN)        :: submesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'move_data_from_submesh_to_mesh'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    mesh%nV          = submesh%nV
    mesh%nTri        = submesh%nTri
    mesh%nV_mem      = submesh%nV_mem
    mesh%nTri_mem    = submesh%nTri_mem
    
    mesh%xmin        = submesh%xmin
    mesh%xmax        = submesh%xmax
    mesh%ymin        = submesh%ymin
    mesh%ymax        = submesh%ymax
    mesh%tol_dist    = submesh%tol_dist
    
    mesh%perturb_dir = submesh%perturb_dir
        
    mesh%V(              1:submesh%nV  ,:) = submesh%V(              1:submesh%nV  ,:)
    mesh%nC(             1:submesh%nV    ) = submesh%nC(             1:submesh%nV    )
    mesh%C(              1:submesh%nV  ,:) = submesh%C(              1:submesh%nV  ,:)
    mesh%niTri(          1:submesh%nV    ) = submesh%niTri(          1:submesh%nV    )
    mesh%iTri(           1:submesh%nV  ,:) = submesh%iTri(           1:submesh%nV  ,:)
    mesh%edge_index(     1:submesh%nV    ) = submesh%edge_index(     1:submesh%nV    )
    mesh%mesh_old_ti_in( 1:submesh%nV    ) = submesh%mesh_old_ti_in( 1:submesh%nV    )
    
    mesh%Tri(            1:submesh%nTri,:) = submesh%Tri(            1:submesh%nTri,:)
    mesh%Tricc(          1:submesh%nTri,:) = submesh%Tricc(          1:submesh%nTri,:)
    mesh%TriC(           1:submesh%nTri,:) = submesh%TriC(           1:submesh%nTri,:)
    mesh%Tri_edge_index( 1:submesh%nTri  ) = submesh%Tri_edge_index( 1:submesh%nTri  )
    
    mesh%Triflip(        1:submesh%nTri,:) = submesh%Triflip(        1:submesh%nTri,:)
    mesh%RefMap(         1:submesh%nTri  ) = submesh%RefMap(         1:submesh%nTri  )
    mesh%RefStack(       1:submesh%nTri  ) = submesh%RefStack(       1:submesh%nTri  )
    mesh%RefStackN                         = submesh%RefStackN
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE move_data_from_submesh_to_mesh
  SUBROUTINE move_data_from_submesh_to_mesh_new( mesh, submesh)
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    TYPE(type_mesh_new),                 INTENT(IN)        :: submesh
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'move_data_from_submesh_to_mesh'
    
    ! Add routine to path
    CALL init_routine( routine_name)
    
    mesh%nV          = submesh%nV
    mesh%nTri        = submesh%nTri
    mesh%nV_mem      = submesh%nV_mem
    mesh%nTri_mem    = submesh%nTri_mem
    
    mesh%xmin        = submesh%xmin
    mesh%xmax        = submesh%xmax
    mesh%ymin        = submesh%ymin
    mesh%ymax        = submesh%ymax
    mesh%tol_dist    = submesh%tol_dist
    
    mesh%perturb_dir = submesh%perturb_dir
        
    mesh%V(              1:submesh%nV  ,:) = submesh%V(              1:submesh%nV  ,:)
    mesh%nC(             1:submesh%nV    ) = submesh%nC(             1:submesh%nV    )
    mesh%C(              1:submesh%nV  ,:) = submesh%C(              1:submesh%nV  ,:)
    mesh%niTri(          1:submesh%nV    ) = submesh%niTri(          1:submesh%nV    )
    mesh%iTri(           1:submesh%nV  ,:) = submesh%iTri(           1:submesh%nV  ,:)
    mesh%edge_index(     1:submesh%nV    ) = submesh%edge_index(     1:submesh%nV    )
    mesh%mesh_old_ti_in( 1:submesh%nV    ) = submesh%mesh_old_ti_in( 1:submesh%nV    )
    
    mesh%Tri(            1:submesh%nTri,:) = submesh%Tri(            1:submesh%nTri,:)
    mesh%Tricc(          1:submesh%nTri,:) = submesh%Tricc(          1:submesh%nTri,:)
    mesh%TriC(           1:submesh%nTri,:) = submesh%TriC(           1:submesh%nTri,:)
    mesh%Tri_edge_index( 1:submesh%nTri  ) = submesh%Tri_edge_index( 1:submesh%nTri  )
    
    mesh%Triflip(        1:submesh%nTri,:) = submesh%Triflip(        1:submesh%nTri,:)
    mesh%RefMap(         1:submesh%nTri  ) = submesh%RefMap(         1:submesh%nTri  )
    mesh%RefStack(       1:submesh%nTri  ) = submesh%RefStack(       1:submesh%nTri  )
    mesh%RefStackN                         = submesh%RefStackN
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE move_data_from_submesh_to_mesh_new

END MODULE mesh_memory_module
