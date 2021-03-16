MODULE mesh_rectangular_module

  USE configuration_module,        ONLY: dp, C
  USE parallel_module,             ONLY: par, sync, allocate_shared_memory_int_0D, allocate_shared_memory_dp_0D, &
                                                    allocate_shared_memory_int_1D, allocate_shared_memory_dp_1D, &
                                                    allocate_shared_memory_int_2D, allocate_shared_memory_dp_2D, &
                                                    allocate_shared_memory_int_3D, allocate_shared_memory_dp_3D, &
                                                    allocate_shared_memory_bool_1D, deallocate_shared_memory, &
                                                    adapt_shared_memory_int_1D,    adapt_shared_memory_dp_1D, &
                                                    adapt_shared_memory_int_2D,    adapt_shared_memory_dp_2D, &
                                                    adapt_shared_memory_int_3D,    adapt_shared_memory_dp_3D, &
                                                    adapt_shared_memory_bool_1D                   
  USE data_types_module,           ONLY: type_mesh
  USE mesh_help_functions_module,  ONLY: FindConnectionWidths, FindTriangleAreas, FindVoronoiCellAreas, GetLatLonCoordinates, &
                                         DetermineMeshResolution, oblique_sg_projection, WriteMeshToScreen, CheckMesh, &
                                         SwitchVertices, FindPOIXYCoordinates, FindPOIVerticesAndWeights, FindVoronoiCellGeometricCentres
  USE mesh_memory_module,          ONLY: AllocateMesh, CropMeshMemory, AllocateMesh_extra
  USE mesh_derivatives_module,     ONLY: GetNeighbourFunctions, GetMeshCurvatures, GetMeshCurvaturesVertex
  USE mesh_ArakawaC_module,        ONLY: MakeAcMesh

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
  CONTAINS
  
  SUBROUTINE CreateSemiRectangularMesh(mesh, region_name, xmin, xmax, ymin, ymax)
    ! Create a structured "semi-rectangular" mesh
  
    ! Input variables (initial data)
    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    CHARACTER(LEN=3),           INTENT(IN)        :: region_name
    REAL(dp),                   INTENT(IN)        :: xmin
    REAL(dp),                   INTENT(IN)        :: xmax
    REAL(dp),                   INTENT(IN)        :: ymin
    REAL(dp),                   INTENT(IN)        :: ymax
    
    ! Local variables
    INTEGER                                       :: nx, ny
    REAL(dp)                                      :: dx, dy
    INTEGER                                       :: nV, i, j, nTri_per_col, til, tir, vi, v1, v2, v3, v4
    INTEGER                                       :: vi_SE, vi_NE, vi_NW
    REAL(dp), PARAMETER                           :: tol = 1E-9_dp
    
    nx = NINT((xmax - xmin) / (C%res_rectangular * 1000._dp))+1
    ny = nx
    
    nV = nx*ny
    CALL AllocateMesh(mesh, region_name, nV+100)
        
    mesh%xmin = xmin
    mesh%xmax = xmax
    mesh%ymin = ymin
    mesh%ymax = ymax
    
    ! Horizontal distance of tolerance. Used for several small routines - points that lie within
    ! this distance of each other (vertices, triangle circumcenters, etc.) are treated as identical.
    mesh%tol_dist   = ((mesh%xmax - mesh%xmin) + (mesh%ymax-mesh%ymin)) * tol / 2._dp

    dx = (mesh%xmax - mesh%xmin) / REAL(nx-1)
    dy = (mesh%ymax - mesh%ymin) / REAL(ny-1)

    mesh%nV      = nx * ny
    mesh%nTri    = (nx-1) * (ny-1) * 2
        
    nTri_per_col = 2 * (ny-1)
    
    IF (par%master) THEN
    
    DO i = 1, nx
      DO j = 1, ny

        ! Vertex index
        vi = (i-1)*ny + j

        ! Vertex coordinates
        mesh%V(vi,:) = [mesh%xmin + REAL(i-1)*dx, mesh%ymin + REAL(j-1)*dy]

        ! Triangle indices
        til = (i-1)*nTri_per_col + (j-1)*2 + 1
        tir = (i-1)*nTri_per_col + (j-1)*2 + 2

        v1 = vi
        v2 = vi + ny
        v3 = v2 + 1
        v4 = v1 + 1

        IF (i<nx .AND. j<ny) THEN
          mesh%Tri(til,:) = [v1, v2, v4]
          mesh%Tri(tir,:) = [v2, v3, v4]
        END IF

        ! Determine other properties - vertex connectivity, inverse triangles,
        ! edge index, etc.
        IF (i==1) THEN
          IF (j==1) THEN
            ! SW corner
            mesh%nC(              vi    ) = 2
            mesh%C(               vi,1:2) = [ny+1, 2]
            mesh%niTri(           vi    ) = 1
            mesh%iTri(            vi,1  ) = 1
            mesh%edge_index(      vi    ) = 6

            mesh%Tricc(           til,:)  = mesh%V(vi,:) + [dx/2, dy/2]
            mesh%TriC(            til,:)  = [2, 0, 0]
            mesh%Tri_edge_index(  til  )  = 5

            mesh%Tricc(           tir,:)  = mesh%V(vi,:) + [dx/2, dy/2]
            mesh%TriC(            tir,:)  = [3, 1, nTri_per_col+1]
            mesh%Tri_edge_index(  tir  )  = 0
          ELSEIF (j==ny) THEN
            ! NW corner
            vi_NW = vi
            mesh%nC(              vi    ) = 3
            mesh%C(               vi,1:3) = [ny-1, 2*ny-1, 2*ny]
            mesh%niTri(           vi    ) = 2
            mesh%iTri(            vi,1:2) = [nTri_per_col-1, nTri_per_col]
            mesh%edge_index(      vi    ) = 8
          ELSE
            ! W edge
            mesh%nC(              vi    ) = 4
            mesh%C(               vi,1:4) = [vi-1, vi+ny-1, vi+ny, vi+1]
            mesh%niTri(           vi    ) = 3
            mesh%iTri(            vi,1:3) = [vi*2-3, vi*2-2, vi*2-1]
            mesh%edge_index(      vi    ) = 7

            mesh%Tricc(           til,:)  = mesh%V(vi,:) + [dx/2, dy/2]
            mesh%TriC(            til,:)  = [til+1, 0, til-1]
            mesh%Tri_edge_index(  til  )  = 7

            mesh%Tricc(           tir,:)  = mesh%V(vi,:) + [dx/2, dy/2]
            mesh%TriC(            tir,:)  = [tir+1, tir-1, tir+nTri_per_col-1]
            mesh%Tri_edge_index(  tir  )  = 0

            IF (j==ny-1) THEN
              mesh%TriC(          tir,:)  = [0, tir-1, tir+nTri_per_col-1]
              mesh%Tri_edge_index(tir  )  = 1
            END IF
          END IF
        ELSEIF (i==nx) THEN
          IF (j==1) THEN
            ! SE corner
            vi_SE = vi
            mesh%nC(              vi    ) = 3
            mesh%C(               vi,1:3) = [vi+1, vi-ny+1, vi-ny]
            mesh%niTri(           vi    ) = 2
            mesh%iTri(            vi,1:2) = [nTri_per_col*(nx-2)+2, nTri_per_col*(nx-2)+1]
            mesh%edge_index(      vi    ) = 4
          ELSEIF (j==ny) THEN
            ! NE corner
            vi_NE = vi
            mesh%nC(              vi    ) = 2
            mesh%C(               vi,1:2) = [vi-ny, vi-1]
            mesh%niTri(           vi    ) = 1
            mesh%iTri(            vi,1  ) = mesh%nTri
            mesh%edge_index(      vi    ) = 2
          ELSE
            ! E edge
            mesh%nC(              vi    ) = 4
            mesh%C(               vi,1:4) = [vi+1, vi-ny+1, vi-ny, vi-1]
            mesh%niTri(           vi    ) = 3
            mesh%iTri(            vi,1:3) = [nTri_per_col*(nx-2)+(j*2), nTri_per_col*(nx-2)+(j*2)-1, nTri_per_col*(nx-2)+(j*2)-2]
            mesh%edge_index(      vi    ) = 3
          END IF
        ELSE
          if (j==1) THEN
            ! S edge
            mesh%nC(              vi    ) = 4
            mesh%C(               vi,1:4) = [vi+ny, vi+1, vi-ny+1, vi-ny]
            mesh%niTri(           vi    ) = 3
            mesh%iTri(            vi,1:3) = [nTri_per_col*(i-1)+1, nTri_per_col*(i-2)+2, nTri_per_col*(i-2)+1]
            mesh%edge_index(      vi    ) = 5

            mesh%Tricc(           til,:)  = mesh%V(vi,:) + [dx/2, dy/2]
            mesh%TriC(            til,:)  = [til+1, til-nTri_per_col+1, 0]
            mesh%Tri_edge_index(  til  )  = 5

            mesh%Tricc(           tir,:)  = mesh%V(vi,:) + [dx/2, dy/2]
            mesh%TriC(            tir,:)  = [tir+1, tir-1, tir+nTri_per_col-1]
            mesh%Tri_edge_index(  tir  )  = 0

            IF (i==nx-1) THEN
              mesh%TriC(          tir,:)  = [tir+1, tir-1, 0]
              mesh%Tri_edge_index(tir  )  = 3
            END IF
          ELSEIF (j==ny) THEN
            ! N edge
            mesh%nC(              vi    ) = 4
            mesh%C(               vi,1:4) = [vi-ny, vi-1, vi+ny-1, vi+ny]
            mesh%niTri(           vi    ) = 3
            mesh%iTri(            vi,1:3) = [nTri_per_col*(i-1), nTri_per_col*i-1, nTri_per_col*i]
            mesh%edge_index(      vi    ) = 1
          ELSE
            ! non-edge vertex
            mesh%nC(              vi    ) = 6
            mesh%C(               vi,1:6) = [vi+ny, vi+1, vi-ny+1, vi-ny, vi-1, vi+ny-1]
            mesh%niTri(           vi    ) = 6
            mesh%iTri(            vi,1:6) = [ nTri_per_col*(i-1)+j*2-1, &
                                              nTri_per_col*(i-2)+j*2, &
                                              nTri_per_col*(i-2)+j*2-1, &
                                              nTri_per_col*(i-2)+j*2-2, &
                                              nTri_per_col*(i-1)+j*2-3, &
                                              nTri_per_col*(i-1)+j*2-2]
            mesh%edge_index(      vi    ) = 0

            mesh%Tricc(           til,:)  = mesh%V(vi,:) + [dx/2, dy/2]
            mesh%TriC(            til,:)  = [til+1, til-nTri_per_col+1, til-1]
            mesh%Tri_edge_index(  til  )  = 0

            mesh%Tricc(           tir,:)  = mesh%V(vi,:) + [dx/2, dy/2]
            mesh%TriC(            tir,:)  = [tir+1, tir-1, tir+nTri_per_col-1]
            mesh%Tri_edge_index(  tir  )  = 0

            IF (i<nx-1 .AND. j==ny-1) THEN
              mesh%TriC(          tir,:)  = [0, tir-1, tir+nTri_per_col-1]
              mesh%Tri_edge_index(tir  )  = 1
            ELSEIF (i==nx-1 .AND. j<ny-1) THEN
              mesh%TriC(          tir,:)  = [tir+1, tir-1, 0]
              mesh%Tri_edge_index(tir  )  = 3
            ELSEIF (i==nx-1 .AND. j==ny-1) THEN
              mesh%TriC(          tir,:)  = [0, tir-1, 0]
              mesh%Tri_edge_index(tir  )  = 1
            END IF
          END IF
        END IF 

      END DO
    END DO
    
    ! Make sure vertices 2, 3 and 4 are at the SE, NE and NW corners (some routines expect this)
    CALL SwitchVertices( mesh, 2, vi_SE)
    CALL SwitchVertices( mesh, 3, vi_NE)
    CALL SwitchVertices( mesh, 4, vi_NW)
    
    END IF ! IF (par%master) THEN

    ! Crop surplus allocated memory
    CALL CropMeshMemory(mesh)
    
    ! Determine vertex and triangle domains
    mesh%v1 = MAX(1,         FLOOR(REAL(mesh%nV   *  par%i      / par%n)) + 1)
    mesh%v2 = MIN(mesh%nV,   FLOOR(REAL(mesh%nV   * (par%i + 1) / par%n)))
    mesh%t1 = MAX(1,         FLOOR(REAL(mesh%nTri *  par%i      / par%n)) + 1)
    mesh%t2 = MIN(mesh%nTri, FLOOR(REAL(mesh%nTri * (par%i + 1) / par%n)))
    CALL sync
    mesh%dz_max_ice     = C%dz_max_ice
    mesh%alpha_min      = C%mesh_alpha_min
    
    CALL CheckMesh( mesh)
    
    ! Calculate extra mesh data
    CALL AllocateMesh_extra(      mesh)       
    CALL GetLatLonCoordinates(    mesh)
    CALL FindVoronoiCellAreas(    mesh)
    CALL FindTriangleAreas(       mesh)
    CALL FindConnectionWidths(    mesh)
    CALL MakeAcMesh(              mesh)
    CALL GetNeighbourFunctions(   mesh)
    CALL DetermineMeshResolution( mesh)
    CALL FindPOIXYCoordinates(            mesh)
    CALL FindPOIVerticesAndWeights(       mesh)
    CALL FindVoronoiCellGeometricCentres( mesh)

    IF (par%master) THEN
      WRITE(0,'(A)')                '   Finished creating mesh.'
      WRITE(0,'(A,I6)')             '    Vertices  : ', mesh%nV
      WRITE(0,'(A,I6)')             '    Triangles : ', mesh%nTri
      WRITE(0,'(A,F5.1,A,F5.1,A)')  '    Resolution: ', mesh%resolution_min/1000._dp, ' - ', mesh%resolution_max/1000._dp, ' km'
    END IF
    
  END SUBROUTINE CreateSemiRectangularMesh

END MODULE mesh_rectangular_module
