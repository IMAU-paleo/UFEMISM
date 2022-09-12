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
  USE utilities_module,                ONLY: check_for_NaN_dp_1D,  check_for_NaN_dp_2D,  check_for_NaN_dp_3D, &
                                             check_for_NaN_int_1D, check_for_NaN_int_2D, check_for_NaN_int_3D

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
  USE data_types_module,               ONLY: type_mesh
  USE sparse_matrix_module,            ONLY: deallocate_matrix_CSR

  IMPLICIT NONE

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
    CALL allocate_shared_dp_0D(  mesh%beta_stereo,      mesh%wbeta_stereo     )

    IF (mesh%region_name     == 'NAM') THEN
      mesh%lambda_M     = C%lambda_M_NAM
      mesh%phi_M        = C%phi_M_NAM
      mesh%beta_stereo  = C%beta_stereo_NAM
    ELSEIF (mesh%region_name == 'EAS') THEN
      mesh%lambda_M     = C%lambda_M_EAS
      mesh%phi_M        = C%phi_M_EAS
      mesh%beta_stereo  = C%beta_stereo_EAS
    ELSEIF (mesh%region_name == 'GRL') THEN
      mesh%lambda_M     = C%lambda_M_GRL
      mesh%phi_M        = C%phi_M_GRL
      mesh%beta_stereo  = C%beta_stereo_GRL
    ELSEIF (mesh%region_name == 'ANT') THEN
      mesh%lambda_M     = C%lambda_M_ANT
      mesh%phi_M        = C%phi_M_ANT
      mesh%beta_stereo  = C%beta_stereo_ANT
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
    CALL deallocate_shared( mesh%wbeta_stereo)
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
    CALL allocate_shared_dist_dp_0D(  mesh%beta_stereo,      mesh%wbeta_stereo     )

    IF (mesh%region_name     == 'NAM') THEN
      mesh%lambda_M     = C%lambda_M_NAM
      mesh%phi_M        = C%phi_M_NAM
      mesh%beta_stereo  = C%beta_stereo_NAM
    ELSEIF (mesh%region_name == 'EAS') THEN
      mesh%lambda_M     = C%lambda_M_EAS
      mesh%phi_M        = C%phi_M_EAS
      mesh%beta_stereo  = C%beta_stereo_EAS
    ELSEIF (mesh%region_name == 'GRL') THEN
      mesh%lambda_M     = C%lambda_M_GRL
      mesh%phi_M        = C%phi_M_GRL
      mesh%beta_stereo  = C%beta_stereo_GRL
    ELSEIF (mesh%region_name == 'ANT') THEN
      mesh%lambda_M     = C%lambda_M_ANT
      mesh%phi_M        = C%phi_M_ANT
      mesh%beta_stereo  = C%beta_stereo_ANT
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
    CALL deallocate_shared( mesh%wbeta_stereo)
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

END MODULE mesh_memory_module
