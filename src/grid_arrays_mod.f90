MODULE grid_arrays_mod

USE constants_mod

IMPLICIT NONE


TYPE particle
    REAL(4) :: x(3) = 0.0  ! position
    REAL(4) :: v(3) = 0.0  ! velocity
END TYPE particle

REAL(4) :: dx, Lx, Ly, Lz
INTEGER :: Ne, Np, tstep
REAL(4), ALLOCATABLE :: B(:,:,:,:)
TYPE(particle), ALLOCATABLE :: electron(:), proton(:)



CONTAINS



SUBROUTINE create_grid_arrays()

    REAL(4) :: mem
    
    ! initialize particle numbers
    Ne = 0
    Np = 0

    ! initilize grid spacing and size
    dx = 1.0 / MAX(nx,ny,nz)
    Lx = nx * dx
    Ly = ny * dx
    Lz = nz * dx

    ! compute storage requirements 
    mem = 4.0 * ((nx+2.0*nb) * (ny+2.0*nb) * (nz+2.0*nb) + 6.0 * (Ne + Np)) * 1.e-6

    PRINT*,'Total memory required for grid array allocation (Mb) = ', mem

    ALLOCATE(B(1-nb:nx+nb,1-nb:ny+nb,1-nb:nz+nb,3))
    ALLOCATE(electron(1:Ne_max), proton(1:Np_max))    


    PRINT*,''
    PRINT*,'dx, Lx, Ly, Lz = ', dx, Lx, Ly, Lz
    

END SUBROUTINE create_grid_arrays



SUBROUTINE destroy_grid_arrays()

    DEALLOCATE(B)
    DEALLOCATE(electron, proton)
    
END SUBROUTINE destroy_grid_arrays



END MODULE grid_arrays_mod