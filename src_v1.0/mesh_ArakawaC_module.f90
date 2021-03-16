MODULE mesh_ArakawaC_module
  ! Routines used in creating the Arakawa C staggered mesh.

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
  USE mesh_help_functions_module,  ONLY: IsBoundarySegment

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'

  CONTAINS
    
  ! == Subroutines for the Arakawa C mesh
  SUBROUTINE MakeAcMesh(mesh)
    ! Determine the Arakawa C mesh vertices and their neighbour functions.
    ! For a detailed explanation, read the documentation (not kidding, its in there).

    TYPE(type_mesh),            INTENT(INOUT)     :: mesh
    
    INTEGER                                       :: vi, ci, vj, cj, iti, ti, n1, n2, n3, vl, vr
    INTEGER,  DIMENSION(:,:), ALLOCATABLE         :: Aci_temp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: VAc_temp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Nx_ac_temp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: Ny_ac_temp
    REAL(dp), DIMENSION(:  ), ALLOCATABLE         :: Np_ac_temp
    REAL(dp), DIMENSION(:,:), ALLOCATABLE         :: No_ac_temp
    REAL(dp), DIMENSION(4)                        :: Nxl, Nxr, Nyl, Nyr
    REAL(dp)                                      :: Nzl, Nzr
    REAL(dp), DIMENSION(4)                        :: Nx, Ny
    REAL(dp)                                      :: Ux, Uy, U, Np
    REAL(dp), DIMENSION(4)                        :: No
                
    ! Allocate memory for the mapping arrays. Aci is stored in temporary memory on the master process first,
    ! once we know how much is needed, we will allocate shared memory and transfer the data there.
    CALL allocate_shared_memory_int_2D( mesh%nV,     mesh%nconmax, mesh%iAci,            mesh%wiAci          )
    CALL allocate_shared_memory_int_0D(                            mesh%nAC,             mesh%wnAC           )
    
    ! Go through all vertex connections, determine the location of the Ac vertex (on the connection midpoint),
    ! the indices of the left and right Aa vertices, and the neighbour functions on the Ac vertex
    ! Only done by master
    ! ====================================
    
    vr = 0
    vl = 0
    
    IF (.NOT. par%master) THEN
    
      ! Allocate a little memory to prevent compiler warnings
      ALLOCATE(VAc_temp(   1,1))
      ALLOCATE(Aci_temp(   1,1))
      ALLOCATE(Nx_ac_temp( 1,1))
      ALLOCATE(Ny_ac_temp( 1,1))
      ALLOCATE(Np_ac_temp( 1  ))
      ALLOCATE(No_ac_temp( 1,1))
    
    ELSE
    
      ! Allocate temporary memory for VAc, Aci and neighbour functions
      ALLOCATE(VAc_temp(   mesh%nTri * 3,2))
      ALLOCATE(Aci_temp(   mesh%nTri * 3,4))
      ALLOCATE(Nx_ac_temp( mesh%nTri * 3,4))
      ALLOCATE(Ny_ac_temp( mesh%nTri * 3,4))
      ALLOCATE(Np_ac_temp( mesh%nTri * 3  ))
      ALLOCATE(No_ac_temp( mesh%nTri * 3,4))
      
      mesh%nAc   = 0
      mesh%iAci  = 0
      Aci_temp   = 0
      VAc_temp   = 0._dp
      Nx_ac_temp = 0._dp
      Ny_ac_temp = 0._dp
      Np_ac_temp = 0._dp
      No_ac_temp = 0._dp
      
      ! Go over all Aa vertices
      DO vi = 1, mesh%nV
      
        ! Go over all the vertex connections
        DO ci = 1, mesh%nC(vi)
          vj = mesh%C(vi,ci)
          
          ! Skip connections that were already considered in the opposite direction
          IF (mesh%iAci(vi,ci)>0) CYCLE
          
          ! List Ac vertex number in reference array
          mesh%nAc = mesh%nAc + 1
          mesh%iAci(vi,ci) = mesh%nAc
          VAc_temp(mesh%nAc,:) = (mesh%V(vi,:) + mesh%V(vj,:)) / 2._dp
            
          ! Find reverse connection
          DO cj = 1, mesh%nC(vj)
            IF (mesh%C(vj,cj)==vi) THEN
              mesh%iAci(vj,cj) = mesh%nAC
              EXIT
            END IF
          END DO
          
          IF (.NOT. IsBoundarySegment( mesh, vi, vj)) THEN
              
            ! Find indices of vl and vr
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO n1 = 1, 3
                n2 = n1+1
                IF (n2==4) n2=1
                n3 = n2+1
                IF (n3==4) n3=1
                
                IF (mesh%Tri(ti,n1)==vi .AND. mesh%Tri(ti,n2)==vj) THEN
                  vl = mesh%Tri(ti,n3)
                ELSE IF (mesh%Tri(ti,n1)==vj .AND. mesh%Tri(ti,n2)==vi) THEN
                  vr = mesh%Tri(ti,n3)
                END IF
              END DO
            END DO
            
            ! List relevant Aa vertex indices vor Ac vertex in reference array
            Aci_temp(mesh%nAc,:) = [vi, vj, vl, vr]
            
            ! Calculate neighbour functions - first on adjacent triangles, then on Ac vertex
            ! x an y just as average of the left and right triangles
            Nxl = [mesh%V(vl,2)-mesh%V(vj,2), mesh%V(vi,2)-mesh%V(vl,2), mesh%V(vj,2)-mesh%V(vi,2), 0._dp]
            Nyl = [mesh%V(vj,1)-mesh%V(vl,1), mesh%V(vl,1)-mesh%V(vi,1), mesh%V(vi,1)-mesh%V(vj,1), 0._dp]
            Nxr = [mesh%V(vj,2)-mesh%V(vr,2), mesh%V(vr,2)-mesh%V(vi,2), 0._dp, mesh%V(vi,2)-mesh%V(vj,2)]
            Nyr = [mesh%V(vr,1)-mesh%V(vj,1), mesh%V(vi,1)-mesh%V(vr,1), 0._dp, mesh%V(vj,1)-mesh%V(vi,1)]
            Nzl = ( (mesh%V(vj,1)-mesh%V(vi,1)) * (mesh%V(vl,2)-mesh%V(vi,2)) ) - &
                  ( (mesh%V(vj,2)-mesh%V(vi,2)) * (mesh%V(vl,1)-mesh%V(vi,1)) )
            Nzr = ( (mesh%V(vr,1)-mesh%V(vi,1)) * (mesh%V(vj,2)-mesh%V(vi,2)) ) - &
                  ( (mesh%V(vr,2)-mesh%V(vi,2)) * (mesh%V(vj,1)-mesh%V(vi,1)) )
                  
            Nx = -((Nxl / Nzl) + (Nxr / Nzr))/2._dp
            Ny = -((Nyl / Nzl) + (Nyr / Nzr))/2._dp
            
          ELSE
            ! Surface slope on an edge segment Ac vertex equals that of the (single) triangle
              
            ! Find index of vl
            DO iti = 1, mesh%niTri(vi)
              ti = mesh%iTri(vi,iti)
              DO n1 = 1, 3
                n2 = n1+1
                IF (n2==4) n2=1
                n3 = n2+1
                IF (n3==4) n3=1
                
                IF ((mesh%Tri(ti,n1)==vi .AND. mesh%Tri(ti,n2)==vj) .OR. (mesh%Tri(ti,n1)==vj .AND. mesh%Tri(ti,n2)==vi)) THEN
                  vl = mesh%Tri(ti,n3)
                END IF
              END DO
            END DO
            
            ! List relevant Aa vertex indices vor Ac vertex in reference array
            ! NOTE: fourth entry is set to 1 because it cant be zero, but the fourth
            !       neighbour function is set to zero, so it doesn't matter...
            Aci_temp(mesh%nAc,:) = [vi, vj, vl, 1]
            
            ! Calculate neighbour functions on single adjacent triangle (equal to Ac vertex values)
            ! x an y just as average of the left and right triangles
            Nxl = [mesh%V(vl,2)-mesh%V(vj,2), mesh%V(vi,2)-mesh%V(vl,2), mesh%V(vj,2)-mesh%V(vi,2), 0._dp]
            Nyl = [mesh%V(vj,1)-mesh%V(vl,1), mesh%V(vl,1)-mesh%V(vi,1), mesh%V(vi,1)-mesh%V(vj,1), 0._dp]
            Nzl = ( (mesh%V(vj,1)-mesh%V(vi,1)) * (mesh%V(vl,2)-mesh%V(vi,2)) ) - &
                  ( (mesh%V(vj,2)-mesh%V(vi,2)) * (mesh%V(vl,1)-mesh%V(vi,1)) )
                  
            Nx = -Nxl / Nzl
            Ny = -Nyl / Nzl           
          
          END IF
            
          ! Then p and o by rotating x and y
          Ux = mesh%V(vj,1)-mesh%V(vi,1)
          Uy = mesh%V(vj,2)-mesh%V(vi,2)
          U  = SQRT(Ux**2+Uy**2)
          
          Np = 1._dp / U          ! Slope along line from j to i is just f(j)-f(i)/dist(ji), only need to save that single number
                                  ! Checked that this gives the same result as just rotating Nx, Ny
          No = (Ny*Ux - Nx*Uy)/U  ! This one still depends on all four vertices, no symmetries...
          
          ! Save everything in the temporary memory
          Nx_ac_temp(mesh%nAc,:) = Nx
          Ny_ac_temp(mesh%nAc,:) = Ny
          Np_ac_temp(mesh%nAc  ) = Np
          No_ac_temp(mesh%nAc,:) = No 
            
        END DO ! DO ci = 1, mesh%nC(vi)
      END DO ! DO vi = 1, mesh%nV
    
    END IF ! IF (par%master) THEN
    CALL sync
    
    ! Allocate shared memory, move data there
    ! =======================================
    
    CALL allocate_shared_memory_dp_2D(  mesh%nAc,       2,            mesh%VAc,             mesh%wVAc           )
    CALL allocate_shared_memory_int_2D( mesh%nAc,       4,            mesh%Aci,             mesh%wAci           )
    CALL allocate_shared_memory_dp_2D(  mesh%nAc,       4,            mesh%Nx_ac,           mesh%wNx_ac         )
    CALL allocate_shared_memory_dp_2D(  mesh%nAc,       4,            mesh%Ny_ac,           mesh%wNy_ac         )
    CALL allocate_shared_memory_dp_1D(  mesh%nAc,                     mesh%Np_ac,           mesh%wNp_ac         )
    CALL allocate_shared_memory_dp_2D(  mesh%nAc,       4,            mesh%No_ac,           mesh%wNo_ac         )
    
    IF (par%master) THEN
      mesh%VAc   = VAc_temp(  1:mesh%nAc,:)
      mesh%Aci   = Aci_temp(  1:mesh%nAc,:)
      mesh%Nx_ac = Nx_ac_temp(1:mesh%nAc,:)
      mesh%Ny_ac = Ny_ac_temp(1:mesh%nAc,:)
      mesh%Np_ac = Np_ac_temp(1:mesh%nAc  )
      mesh%No_ac = No_ac_temp(1:mesh%nAc,:)
    END IF
    
    ! Deallocate temporary memory
    ! ===========================
    
    DEALLOCATE(VAc_temp)
    DEALLOCATE(Aci_temp)
    DEALLOCATE(Nx_ac_temp)
    DEALLOCATE(Ny_ac_temp)
    DEALLOCATE(Np_ac_temp)
    DEALLOCATE(No_ac_temp)
    CALL sync
    
    ! Determine Ac vertex domains
    mesh%c1 = MAX(1,        FLOOR(REAL(mesh%nAc *  par%i      / par%n)) + 1)
    mesh%c2 = MIN(mesh%nAc, FLOOR(REAL(mesh%nAc * (par%i + 1) / par%n)))
    CALL sync
    
  END SUBROUTINE MakeAcMesh
  
  SUBROUTINE GetAcMeshDerivatives(mesh,meshdata,ddx_ac,ddy_ac,ddp_ac,ddo_ac)

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: meshdata
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: ddx_ac, ddy_ac, ddp_ac, ddo_ac

    INTEGER                                       :: ci, n

    ddx_ac(mesh%c1:mesh%c2)  = 0._dp
    ddy_ac(mesh%c1:mesh%c2)  = 0._dp
    ddp_ac(mesh%c1:mesh%c2)  = 0._dp
    ddo_ac(mesh%c1:mesh%c2)  = 0._dp

    DO ci = mesh%c1, mesh%c2
      DO n = 1, 4
        ddx_ac(ci)  = ddx_ac(ci)  + mesh%Nx_ac(ci,n) * meshdata(mesh%Aci(ci,n))
        ddy_ac(ci)  = ddy_ac(ci)  + mesh%Ny_ac(ci,n) * meshdata(mesh%Aci(ci,n))
        ddo_ac(ci)  = ddo_ac(ci)  + mesh%No_ac(ci,n) * meshdata(mesh%Aci(ci,n))
      END DO
      ddp_ac(ci)  = mesh%Np_ac(ci) * (meshdata(mesh%Aci(ci,2)) - meshdata(mesh%Aci(ci,1)))
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
    
  END SUBROUTINE GetAcMeshDerivatives
  
  SUBROUTINE MapAaToAc(mesh,meshdata_aa,meshdata_ac)
    ! Map data from the regular Aa mesh to the staggered Ac mesh

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: meshdata_aa
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: meshdata_ac

    INTEGER                                       :: ci, vi, vj
    
    meshdata_ac(mesh%c1:mesh%c2) = 0._dp
    
    DO ci = mesh%c1, mesh%c2
    
      vi = mesh%Aci(ci,1)
      vj = mesh%Aci(ci,2)
      
      meshdata_ac(ci) = (meshdata_aa(vi) + meshdata_aa(vj)) / 2._dp
    
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
    
  END SUBROUTINE MapAaToAc
  
  SUBROUTINE MapAaToAc_3D(mesh,meshdata_aa,nz,meshdata_ac)
    ! Map data from the regular Aa mesh to the staggered Ac mesh

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: meshdata_aa
    INTEGER,                    INTENT(IN)        :: nz
    REAL(dp), DIMENSION(:,:),   INTENT(OUT)       :: meshdata_ac

    INTEGER                                       :: ci, vi, vj, k
    
    meshdata_ac(mesh%c1:mesh%c2,:) = 0._dp
    
    DO ci = mesh%c1, mesh%c2
    
      vi = mesh%Aci(ci,1)
      vj = mesh%Aci(ci,2)
      
      DO k = 1, nz
        meshdata_ac(ci,k) = (meshdata_aa(vi,k) + meshdata_aa(vj,k)) / 2._dp
      END DO
    
    END DO ! DO ci = mesh%c1, mesh%c2
    CALL sync
    
  END SUBROUTINE MapAaToAc_3D
  
  SUBROUTINE MapAcToAa(mesh,meshdata_ac,meshdata_aa)
    ! Map data from the regular Aa mesh to the staggered Ac mesh

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:),     INTENT(IN)        :: meshdata_ac
    REAL(dp), DIMENSION(:),     INTENT(OUT)       :: meshdata_aa

    INTEGER                                       :: vi, ci, aci
    
    meshdata_aa(mesh%v1:mesh%v2) = 0._dp
    
    DO vi = mesh%v1, mesh%v2
      DO ci = 1, mesh%nC(vi)
        aci = mesh%iAci(vi,ci)
        meshdata_aa(vi) = meshdata_aa(vi) + meshdata_ac(aci) / mesh%nC(vi)
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
  END SUBROUTINE MapAcToAa
  
  SUBROUTINE MapAcToAa_3D(mesh,meshdata_ac,nz,meshdata_aa)
    ! Map data from the regular Aa mesh to the staggered Ac mesh

    TYPE(type_mesh),            INTENT(IN)        :: mesh
    REAL(dp), DIMENSION(:,:),   INTENT(IN)        :: meshdata_ac
    INTEGER,                    INTENT(IN)        :: nz
    REAL(dp), DIMENSION(:,:),   INTENT(OUT)       :: meshdata_aa

    INTEGER                                       :: vi, ci, aci, k
    
    meshdata_aa(mesh%v1:mesh%v2,:) = 0._dp
    
    DO vi = mesh%v1, mesh%v2
      DO ci = 1, mesh%nC(vi)
        aci = mesh%iAci(vi,ci)
        DO k = 1, nz
          meshdata_aa(vi,k) = meshdata_aa(vi,k) + meshdata_ac(aci,k) / mesh%nC(vi)
        END DO
      END DO
    END DO ! DO vi = mesh%v1, mesh%v2
    CALL sync
    
  END SUBROUTINE MapAcToAa_3D


END MODULE mesh_ArakawaC_module