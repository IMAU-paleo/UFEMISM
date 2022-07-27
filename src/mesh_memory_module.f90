MODULE mesh_memory_module
  ! Routines for allocating, deallocating, extending and cropping the memory for the mesh data.

  ! Import basic functionality
#include <petsc/finclude/petscksp.h>
  USE mpi
  USE configuration_module,            ONLY: dp, C, routine_path, init_routine, finalise_routine, crash, warning
  USE parameters_module
  USE petsc_module,                    ONLY: tMat, MatDestroy, perr
  USE parallel_module,                 ONLY: par, sync, ierr, cerr, partition_list

  ! Import specific functionality
  USE data_types_module,               ONLY: type_mesh
  USE sparse_matrix_module,            ONLY: deallocate_matrix_CSR
  use reallocate_mod,                  only: reallocate

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

    ! Something goes wrong and allocated stuff is passed here, so deallocate it
    ! TODO: can be fixed by making the mesh intent out and fixing the bugs leading to this behaviour
    if ( allocated(mesh%V             )) deallocate( mesh%V             )
    if ( allocated(mesh%nC            )) deallocate( mesh%nC            )
    if ( allocated(mesh%C             )) deallocate( mesh%C             )
    if ( allocated(mesh%niTri         )) deallocate( mesh%niTri         )
    if ( allocated(mesh%iTri          )) deallocate( mesh%iTri          )
    if ( allocated(mesh%edge_index    )) deallocate( mesh%edge_index    )
    if ( allocated(mesh%mesh_old_ti_in)) deallocate( mesh%mesh_old_ti_in)
    if ( allocated(mesh%Tri           )) deallocate( mesh%Tri           )
    if ( allocated(mesh%Tricc         )) deallocate( mesh%Tricc         )
    if ( allocated(mesh%TriC          )) deallocate( mesh%TriC          )
    if ( allocated(mesh%Tri_edge_index)) deallocate( mesh%Tri_edge_index)
    if ( allocated(mesh%Triflip       )) deallocate( mesh%Triflip       )
    if ( allocated(mesh%RefMap        )) deallocate( mesh%RefMap        )
    if ( allocated(mesh%RefStack      )) deallocate( mesh%RefStack      )
    if ( allocated(mesh%VMap          )) deallocate( mesh%VMap          )
    if ( allocated(mesh%VStack1       )) deallocate( mesh%VStack1       )
    if ( allocated(mesh%VStack2       )) deallocate( mesh%VStack2       )
    if ( allocated(mesh%TriMap        )) deallocate( mesh%TriMap        )
    if ( allocated(mesh%TriStack1     )) deallocate( mesh%TriStack1     )
    if ( allocated(mesh%TriStack2     )) deallocate( mesh%TriStack2     )
    if ( allocated(mesh%POI_coordinates   )) deallocate(mesh%POI_coordinates   )
    if ( allocated(mesh%POI_XY_coordinates)) deallocate(mesh%POI_XY_coordinates)
    if ( allocated(mesh%POI_resolutions   )) deallocate(mesh%POI_resolutions   )
    if ( allocated(mesh%POI_vi            )) deallocate(mesh%POI_vi            )
    if ( allocated(mesh%POI_w             )) deallocate(mesh%POI_w             )
    
    allocate( mesh%V              (nV_mem,   2     ), source=0._dp)
    allocate( mesh%nC             (nV_mem          ), source=0)
    allocate( mesh%C              (nV_mem,   nC_mem), source=0)
    allocate( mesh%niTri          (nV_mem          ), source=0)
    allocate( mesh%iTri           (nV_mem,   nC_mem), source=0)
    allocate( mesh%edge_index     (nV_mem          ), source=0)
    allocate( mesh%mesh_old_ti_in (nV_mem          ), source=0)

    allocate( mesh%Tri            (nTri_mem, 3     ), source=0)
    allocate( mesh%Tricc          (nTri_mem, 2     ), source=0._dp)
    allocate( mesh%TriC           (nTri_mem, 3     ), source=0) 
    allocate( mesh%Tri_edge_index (nTri_mem        ), source=0) 
                                                     
    allocate( mesh%Triflip        (nTri_mem, 2     ), source=0)
    allocate( mesh%RefMap         (nTri_mem        ), source=0)
    allocate( mesh%RefStack       (nTri_mem        ), source=0)
    
    ! Distributed shared memory for FloodFill maps and stacks
    allocate( mesh%VMap           (nV_mem          ), source=0)
    allocate( mesh%VStack1        (nV_mem          ), source=0)
    allocate( mesh%VStack2        (nV_mem          ), source=0)
    allocate( mesh%TriMap         (nTri_mem        ), source=0)
    allocate( mesh%TriStack1      (nTri_mem        ), source=0)
    allocate( mesh%TriStack2      (nTri_mem        ), source=0)
    
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
   
    mesh%nV_mem   = mesh%nV
    mesh%nTri_mem = mesh%nTri
       
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
    
    if (allocated(mesh% A    )) deallocate(mesh% A    ) 
    if (allocated(mesh% VorGC)) deallocate(mesh% VorGC) 
    if (allocated(mesh% R    )) deallocate(mesh% R    ) 
    if (allocated(mesh% Cw   )) deallocate(mesh% Cw   ) 
    if (allocated(mesh% TriA )) deallocate(mesh% TriA ) 
    if (allocated(mesh% TriGC)) deallocate(mesh% TriGC) 
    if (allocated(mesh% lat  )) deallocate(mesh% lat  ) 
    if (allocated(mesh% lon  )) deallocate(mesh% lon  ) 

    allocate(mesh% A    ( 1       :mesh%nV    )) ! Globally used kinda
    allocate(mesh% VorGC( mesh%vi1:mesh%vi2,2 ))
    allocate(mesh% R    ( 1       :mesh%nV    )) ! Globally used kinda
    allocate(mesh% Cw   ( mesh%vi1:mesh%vi2,mesh%nC_mem ))
        
    allocate(mesh% TriA ( mesh%ti1:mesh%ti2   ))
    allocate(mesh% TriGC( 1:mesh%nTri,2       )) ! Globally needed

    allocate(mesh% lat  ( mesh%vi1:mesh%vi2   ))
    allocate(mesh% lon  ( mesh%vi1:mesh%vi2   ))
    
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
    
!   CALL deallocate_shared( mesh%wlambda_M)
!   CALL deallocate_shared( mesh%wphi_M)
!   CALL deallocate_shared( mesh%walpha_stereo)
!   CALL deallocate_shared( mesh%wxmin)
!   CALL deallocate_shared( mesh%wxmax)
!   CALL deallocate_shared( mesh%wymin)
!   CALL deallocate_shared( mesh%wymax)
!   CALL deallocate_shared( mesh%wtol_dist)
!   CALL deallocate_shared( mesh%wnV_mem)
!   CALL deallocate_shared( mesh%wnTri_mem)
!   CALL deallocate_shared( mesh%wnC_mem)
!   CALL deallocate_shared( mesh%wnV)
!   CALL deallocate_shared( mesh%wnTri)
!   CALL deallocate_shared( mesh%wperturb_dir)
!   CALL deallocate_shared( mesh%walpha_min)
!   CALL deallocate_shared( mesh%wdz_max_ice)
!   CALL deallocate_shared( mesh%wres_max)
!   CALL deallocate_shared( mesh%wres_max_margin)
!   CALL deallocate_shared( mesh%wres_max_gl)
!   CALL deallocate_shared( mesh%wres_max_cf)
!   CALL deallocate_shared( mesh%wres_max_mountain)
!   CALL deallocate_shared( mesh%wres_max_coast)
!   CALL deallocate_shared( mesh%wres_min)
!   CALL deallocate_shared( mesh%wresolution_min)
!   CALL deallocate_shared( mesh%wresolution_max)

!   ! Primary mesh data (needed for mesh creation & refinement)
!   ! =========================================================
!   
!   CALL deallocate_shared( mesh%wV)
!   CALL deallocate_shared( mesh%wnC)
!   CALL deallocate_shared( mesh%wC)
!   CALL deallocate_shared( mesh%wniTri)
!   CALL deallocate_shared( mesh%wiTri)
!   CALL deallocate_shared( mesh%wedge_index)
!   CALL deallocate_shared( mesh%wmesh_old_ti_in)
!   
!   CALL deallocate_shared( mesh%wTri)
!   CALL deallocate_shared( mesh%wTricc)
!   CALL deallocate_shared( mesh%wTriC)
!   CALL deallocate_shared( mesh%wTri_edge_index)
!   
!   CALL deallocate_shared( mesh%wTriflip)
!   CALL deallocate_shared( mesh%wRefMap)
!   CALL deallocate_shared( mesh%wRefStack)
!   CALL deallocate_shared( mesh%wRefStackN)
!   
!   CALL deallocate_shared( mesh%wVMap)
!   CALL deallocate_shared( mesh%wVStack1)
!   CALL deallocate_shared( mesh%wVStack2)
!   CALL deallocate_shared( mesh%wTriMap)
!   CALL deallocate_shared( mesh%wTriStack1)
!   CALL deallocate_shared( mesh%wTriStack2)
!   
!   CALL deallocate_shared( mesh%wnPOI               )
!   CALL deallocate_shared( mesh%wPOI_coordinates    )
!   CALL deallocate_shared( mesh%wPOI_XY_coordinates )
!   CALL deallocate_shared( mesh%wPOI_resolutions    )
!   CALL deallocate_shared( mesh%wPOI_vi             )
!   CALL deallocate_shared( mesh%wPOI_w              )
!   
!   ! Secondary mesh data
!   ! ===================

!   CALL deallocate_shared( mesh%wA)
!   CALL deallocate_shared( mesh%wVorGC)
!   CALL deallocate_shared( mesh%wR)
!   CALL deallocate_shared( mesh%wCw)
!   CALL deallocate_shared( mesh%wTriGC)
!   CALL deallocate_shared( mesh%wTriA)
!   
!   CALL deallocate_shared( mesh%wlat)
!   CALL deallocate_shared( mesh%wlon)
!   
!   CALL deallocate_shared( mesh%wnAc)
!   CALL deallocate_shared( mesh%wVAc)
!   CALL deallocate_shared( mesh%wAci)
!   CALL deallocate_shared( mesh%wiAci)
!   CALL deallocate_shared( mesh%wedge_index_Ac)

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
    
!    CALL deallocate_shared( mesh%wnV_transect)
!    CALL deallocate_shared( mesh%wvi_transect)
!    CALL deallocate_shared( mesh%ww_transect )
    
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
    allocate( mesh%V              (nV_mem,   2     ), source=0._dp)
    allocate( mesh%nC             (nV_mem          ), source=0)
    allocate( mesh%C              (nV_mem,   nC_mem), source=0)
    allocate( mesh%niTri          (nV_mem          ), source=0)
    allocate( mesh%iTri           (nV_mem,   nC_mem), source=0)
    allocate( mesh%edge_index     (nV_mem          ), source=0)
    allocate( mesh%mesh_old_ti_in (nV_mem          ), source=0)

    allocate( mesh%Tri            (nTri_mem, 3     ), source=0)
    allocate( mesh%Tricc          (nTri_mem, 2     ), source=0._dp)
    allocate( mesh%TriC           (nTri_mem, 3     ), source=0) 
    allocate( mesh%Tri_edge_index (nTri_mem        ), source=0) 
                                                     
    allocate( mesh%Triflip        (nTri_mem, 2     ), source=0)
    allocate( mesh%RefMap         (nTri_mem        ), source=0)
    allocate( mesh%RefStack       (nTri_mem        ), source=0)
    
    ! Distributed shared memory for FloodFill maps and stacks
    allocate( mesh%VMap           (nV_mem          ), source=0)
    allocate( mesh%VStack1        (nV_mem          ), source=0)
    allocate( mesh%VStack2        (nV_mem          ), source=0)
    allocate( mesh%TriMap         (nTri_mem        ), source=0)
    allocate( mesh%TriStack1      (nTri_mem        ), source=0)
    allocate( mesh%TriStack2      (nTri_mem        ), source=0)
    
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
    mesh%nTri_mem = mesh%nTri
    
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
    
! Because we switch to allocatable, it will be automatically deallocated when going out of scope, bonus!
!   deallocate(mesh%V)
!   deallocate(mesh%nC)
!   deallocate(mesh%C)
!   deallocate(mesh%niTri)
!   deallocate(mesh%iTri)
!   deallocate(mesh%edge_index)
!   deallocate(mesh%mesh_old_ti_in)
!   
!   deallocate(mesh%Tri)
!   deallocate(mesh%Tricc)
!   deallocate(mesh%TriC)
!   deallocate(mesh%Tri_edge_index)

!   deallocate(mesh%Triflip)
!   deallocate(mesh%RefMap)
!   deallocate(mesh%RefStack)

!   deallocate(mesh%VMap)
!   deallocate(mesh%VStack1)
!   deallocate(mesh%VStack2)
!   deallocate(mesh%TriMap)
!   deallocate(mesh%TriStack1)
!   deallocate(mesh%TriStack2)

!   deallocate(mesh%nPOI               )
!   deallocate(mesh%POI_coordinates    )
!   deallocate(mesh%POI_XY_coordinates )
!   deallocate(mesh%POI_resolutions    )
!   deallocate(mesh%POI_vi             )
!   deallocate(mesh%POI_w              )
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
 
  END SUBROUTINE deallocate_submesh_primary
  SUBROUTINE share_submesh_access( p_left, p_right, submesh, submesh_right)
    ! Give process p_left access to the submesh memory of p_right
    ! Used in submesh merging
    
    IMPLICIT NONE
 
    ! In/output variables:
    INTEGER,                         INTENT(IN)        :: p_left
    INTEGER,                         INTENT(IN)        :: p_right
    INTEGER                                       :: status(MPI_STATUS_SIZE)
    TYPE(type_mesh),                 INTENT(IN)        :: submesh
    TYPE(type_mesh),                 INTENT(INOUT)     :: submesh_right
    
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
      allocate( submesh_right%V             (   nV,   2      ), source=0._dp)
      allocate( submesh_right%nC            (   nV           ), source=0)
      allocate( submesh_right%C             (   nV,   nconmax), source=0)
      allocate( submesh_right%niTri         (   nV           ), source=0)
      allocate( submesh_right%iTri          (   nV,   nconmax), source=0)
      allocate( submesh_right%edge_index    (   nV           ), source=0)
      allocate( submesh_right%mesh_old_ti_in(   nV           ), source=0)

      allocate( submesh_right%Tri           (   nTri, 3      ), source=0)
      allocate( submesh_right%Tricc         (   nTri, 2      ), source=0._dp)
      allocate( submesh_right%TriC          (   nTri, 3      ), source=0)
      allocate( submesh_right%Tri_edge_index(   nTri         ), source=0)

      allocate( submesh_right%Triflip       (   nTri, 2      ), source=0)
      allocate( submesh_right%RefMap        (   nTri         ), source=0)
      allocate( submesh_right%RefStack      (   nTri         ), source=0)

      allocate( submesh_right%VMap          (   nV           ), source=0)
      allocate( submesh_right%VStack1       (   nV           ), source=0)
      allocate( submesh_right%VStack2       (   nV           ), source=0)

      allocate( submesh_right%TriMap        (   nTri         ), source=0)
      allocate( submesh_right%TriStack1     (   nTri         ), source=0)
      allocate( submesh_right%TriStack2     (   nTri         ), source=0)
      
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
    
  END SUBROUTINE share_submesh_access
  
  SUBROUTINE move_data_from_submesh_to_mesh( mesh, submesh)
    use mpi_f08
    
    IMPLICIT NONE
  
    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    TYPE(type_mesh),                 INTENT(INOUT)     :: submesh ! We pointer move from here, so must be inout
    
    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'move_data_from_submesh_to_mesh'
    type(mpi_request)                                  :: reqs(15)
    integer                                            :: nconmax

    nconmax = C%nconmax
    
    ! Add routine to path
    CALL init_routine( routine_name)
    if (par% master) then  

      call crop_submesh_primary(submesh) !make sure nv_mem == nv and ntri == ntri_mem
      if ((submesh%nv /= submesh%nv_mem) .or. (submesh%ntri /= submesh%ntri_mem)) then
        write(0,*) "Error in cropping, memory allocated is still larger then expected"
        error stop
      end if

      mesh%nV          = submesh%nV
      mesh%nTri        = submesh%nTri
      mesh%nV_mem      = submesh%nV_mem
      mesh%nTri_mem    = submesh%nTri_mem
      mesh%nC_mem      = submesh%nC_mem
      
      mesh%xmin        = submesh%xmin
      mesh%xmax        = submesh%xmax
      mesh%ymin        = submesh%ymin
      mesh%ymax        = submesh%ymax
      mesh%tol_dist    = submesh%tol_dist
      
      mesh%perturb_dir = submesh%perturb_dir
    endif

    call mpi_ibcast( mesh%nV         , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 1), ierr )
    call mpi_ibcast( mesh%nTri       , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 2), ierr )
    call mpi_ibcast( mesh%nV_mem     , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 3), ierr )
    call mpi_ibcast( mesh%nTri_mem   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 4), ierr )
    call mpi_ibcast( mesh%nC_mem     , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 5), ierr )

    call mpi_ibcast( mesh%xmin       , 1, MPI_REAL8  , 0, MPI_COMM_WORLD, reqs( 6), ierr )
    call mpi_ibcast( mesh%xmax       , 1, MPI_REAL8  , 0, MPI_COMM_WORLD, reqs( 7), ierr )
    call mpi_ibcast( mesh%ymin       , 1, MPI_REAL8  , 0, MPI_COMM_WORLD, reqs( 8), ierr )
    call mpi_ibcast( mesh%ymax       , 1, MPI_REAL8  , 0, MPI_COMM_WORLD, reqs( 9), ierr )
    call mpi_ibcast( mesh%tol_dist   , 1, MPI_REAL8  , 0, MPI_COMM_WORLD, reqs(10), ierr )

    call mpi_ibcast( mesh%perturb_dir, 1, MPI_REAL8  , 0, MPI_COMM_WORLD, reqs(11), ierr )

    call mpi_waitall( 11, reqs, MPI_STATUSES_IGNORE, ierr)

    CALL allocate_mesh_primary( mesh, submesh%region_name, mesh%nV, mesh%nTri, mesh%nC_mem)

    ! Fill the good data 
    if (par% master) then  
      call move_alloc(submesh%V, mesh%V)
      call move_alloc(submesh%nC, mesh%nC)
      call move_alloc(submesh%C, mesh%C)
      call move_alloc(submesh%niTri, mesh%niTri)
      call move_alloc(submesh%iTri, mesh%iTri)
      call move_alloc(submesh%edge_index, mesh%edge_index)
      call move_alloc(submesh%mesh_old_ti_in, mesh%mesh_old_ti_in)

      call move_alloc(submesh%Tri, mesh%Tri)
      call move_alloc(submesh%tricC, mesh%Tricc)
      call move_alloc(submesh%triC, mesh%triC)
      call move_alloc(submesh%Tri_edge_index, mesh%tri_edge_index)

      call move_alloc(submesh%Triflip, mesh%Triflip)
      call move_alloc(submesh%Refmap, mesh%Refmap)
      call move_alloc(submesh%RefStack, mesh%RefStack)

      mesh% RefStackN = submesh% RefStackN
    end if

    call mpi_ibcast( mesh%V          , mesh%nV*2      , MPI_REAL8  , 0, MPI_COMM_WORLD, reqs( 1), ierr )
    call mpi_ibcast( mesh%nC         , mesh%nV        , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 2), ierr )
    call mpi_ibcast( mesh%C          , mesh%nV*nconmax, MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 3), ierr )
    call mpi_ibcast( mesh%niTri      , mesh%nV        , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 4), ierr )
    call mpi_ibcast( mesh%iTri       , mesh%nV*nconmax, MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 5), ierr )
    call mpi_ibcast( mesh%edge_index , mesh%nV        , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 6), ierr )
    call mpi_ibcast( mesh%mesh_old_ti_in, mesh%nV     , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 7), ierr )

    call mpi_ibcast( mesh%Tri        , mesh%nTri*3    , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs( 8), ierr )
    call mpi_ibcast( mesh%Tricc      , mesh%nTri*2    , MPI_REAL8  , 0, MPI_COMM_WORLD, reqs( 9), ierr )
    call mpi_ibcast( mesh%Tric       , mesh%nTri*3    , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs(10), ierr )
    call mpi_ibcast( mesh%Tri_edge_index, mesh%nTri   , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs(11), ierr )

    call mpi_ibcast( mesh%Triflip    , mesh%nTri*2    , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs(12), ierr )
    call mpi_ibcast( mesh%Refmap     , mesh%nTri      , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs(13), ierr )
    call mpi_ibcast( mesh%RefStack   , mesh%nTri      , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs(14), ierr )
    call mpi_ibcast( mesh%RefStackN  , 1              , MPI_INTEGER, 0, MPI_COMM_WORLD, reqs(15), ierr )

    call mpi_waitall( 15, reqs, MPI_STATUSES_IGNORE, ierr)
    
    ! Finalise routine path
    CALL finalise_routine( routine_name)
    
  END SUBROUTINE move_data_from_submesh_to_mesh

END MODULE mesh_memory_module
