MODULE generate_bfield_mod

USE constants_mod
USE grid_arrays_mod

IMPLICIT NONE

REAL(4), PARAMETER :: gamma = 11.0 / 3.0  ! kolmogorov spectral index [k^2 * k^(-5/3) = k^(11/3)]
INTEGER, PARAMETER :: Nmodes = 100


CONTAINS


! generates B field on the spatial grid
SUBROUTINE generate_bfield_grid_ver1()

    INTEGER :: n, ix, iy, iz
    REAL(4) :: rn(5), theta, phi, alf, beta, zeta(3), k, dk, khat(3), kdotx, x(3), amp, Lc, sigma2, G(Nmodes), Gnorm, kx, ky kz, k0(3) 



    ! correlation length
    Lc = Lx
    k0(1) = TWOPI/Lx
    k0(2) = TWOPI/Ly
    k0(3) = TWOPI/Lz

    ! wave variance squared
    sigma2 = 1.0

    ! k-space grid spacing
    dk = TWOPI / Lc

    ! compute power spectrum envelope function
    Gnorm  = 0.0
    DO n = 1, Nmodes

        k = (n-0.5) * dk
        G(n) = (2.0*TWOPI*(k**2)*dk) / (1.0 + (k*Lc)**gamma)
        Gnorm = Gnorm + G(n)
    
    END DO    
    G = G / Gnorm
   
    ! sum over fourier modes to get the turbulent magnetic field
    B = 0.0
    B(:,:,:,3) = 1.0 ! mean magnetic field 
    
    DO n = 1, Nmodes
    
        !k = (n-0.5) * dk

        ! compute random wave-vector, polarization vector and phase
        CALL RANDOM_NUMBER(rn)
        
        kx = k0(1) * (1+INT(rn(1)*nx))
        ky = k0(2) * (1+INT(rn(2)*ny))
        kz = k0(3) * (1+INT(rn(3)*nz))
        k = SQRT(kx**2 + ky**2 + kz**2)
        
        !theta = 0.5 * TWOPI * rn(1)
        !phi   = TWOPI * rn(2)
    
        theta = ACOS(kz/k) 
        phi   = ASIN(ky/kx)
        
        alf   = TWOPI * rn(4)
        beta  = TWOPI * rn(5)
    
        khat(1) = SIN(theta) * COS(phi)
        khat(2) = SIN(theta) * SIN(phi)
        khat(3) = COS(theta)
        
        zeta(1) = -SIN(phi) * COS(alf) + COS(theta) * COS(phi) * SIN(alf)
        zeta(2) = COS(phi) * COS(alf) + COS(theta) * SIN(phi) * SIN(alf)
        zeta(3) = -SIN(theta) * SIN(alf)
    
        ! compute fourier mode amplitude
        amp = SQRT(sigma2*G(n))
    
        ! now sample the B-field on our finite spatial grid
        DO iz = 1-nb, nz+nb
            x(3) = (iz-0.5) * dx
            DO iy = 1-nb, ny+nb
                x(2) = (iy-0.5) * dx
                DO ix = 1-nb, nx+nb
                    x(1) = (ix-0.5) * dx
                    
                    kdotx = khat(1)*x(1) + khat(2)*x(2) + khat(3)*x(3)  
                    
                    B(ix,iy,iz,1) = B(ix,iy,iz,1) + amp * zeta(1) * COS(k*kdotx + beta)     
                    B(ix,iy,iz,2) = B(ix,iy,iz,2) + amp * zeta(2) * COS(k*kdotx + beta)     
                    B(ix,iy,iz,3) = B(ix,iy,iz,3) + amp * zeta(3) * COS(k*kdotx + beta)     
            
                END DO    
            END DO    
        END DO    
  
    END DO
    



END SUBROUTINE generate_bfield_grid_ver1


! generates B field on the spatial grid
SUBROUTINE generate_bfield_grid_ver2()

    INTEGER :: n, ix, iy, iz, knx, kny, knz 
    REAL(4) :: rn(4), theta, phi, alf, beta, zeta(3), k(3), k0(3), dk, khat(3), kdotx, x(3), amp, Lc, sigma2, G(Nmodes), Gnorm



    ! correlation length
    Lc = Lx
    k0(1) = TWOPI/Lx
    k0(2) = TWOPI/Ly
    k0(3) = TWOPI/Lz

    ! wave variance squared
    sigma2 = 1.0

    ! k-space grid spacing
    dk = TWOPI / Lc

    ! compute power spectrum envelope function
    Gnorm  = 0.0
    DO knz = 1, nz
    DO kny = 1, ny
    DO knx = 1, nx

        k = SQRT((knx*k0(1))**2 + (kny*k0(2)**2 + (knz*k0(3))**2 )
        G(n) = (2.0*TWOPI*(k**2)*dk) / (1.0 + (k*Lc)**gamma)
        Gnorm = Gnorm + G(n)
    
    END DO    
    END DO    
    END DO    
    G = G / Gnorm
   
    ! sum over fourier modes to get the turbulent magnetic field
    B = 0.0
    B(:,:,:,3) = 1.0 ! mean magnetic field 
    
    DO knz = 1, nz
    DO kny = 1, ny
    DO knx = 1, nx
    
        k = SQRT((knx*k0(1))**2 + (kny*k0(2)**2 + (knz*k0(3))**2)

        ! compute random wave-vector, polarization vector and phase
        CALL RANDOM_NUMBER(rn)
        theta = 0.5 * TWOPI * rn(1)
        phi   = TWOPI * rn(2)
        alf   = TWOPI * rn(3)
        beta  = TWOPI * rn(4)
    
        khat(1) = SIN(theta) * COS(phi)
        khat(2) = SIN(theta) * SIN(phi)
        khat(3) = COS(theta)
        
        zeta(1) = -SIN(phi) * COS(alf) + COS(theta) * COS(phi) * SIN(alf)
        zeta(2) = COS(phi) * COS(alf) + COS(theta) * SIN(phi) * SIN(alf)
        zeta(3) = -SIN(theta) * SIN(alf)
    
        ! compute fourier mode amplitude
        amp = SQRT(sigma2*G(n))
    
        ! now sample the B-field on our finite spatial grid
        DO iz = 1-nb, nz+nb
            x(3) = (iz-0.5) * dx
            DO iy = 1-nb, ny+nb
                x(2) = (iy-0.5) * dx
                DO ix = 1-nb, nx+nb
                    x(1) = (ix-0.5) * dx
                    
                    kdotx = khat(1)*x(1) + khat(2)*x(2) + khat(3)*x(3)  
                    
                    B(ix,iy,iz,1) = B(ix,iy,iz,1) + amp * zeta(1) * COS(k*kdotx + beta)     
                    B(ix,iy,iz,2) = B(ix,iy,iz,2) + amp * zeta(2) * COS(k*kdotx + beta)     
                    B(ix,iy,iz,3) = B(ix,iy,iz,3) + amp * zeta(3) * COS(k*kdotx + beta)     
            
                END DO    
            END DO    
        END DO    
  
    END DO
    END DO
    END DO
    



END SUBROUTINE generate_bfield_grid_ver2



! generates B field at a specific location
SUBROUTINE generate_bfield_point(x,Bout)

    REAL(4), INTENT(IN) :: x(3)
    REAL(4), INTENT(INOUT) :: Bout(3)
    INTEGER :: n, ix, iy, iz 
    REAL(4) :: rn(4), theta, phi, alf, beta, zeta(3), k, dk, khat(3), kdotx, amp, Lc, sigma2, G(Nmodes), Gnorm



    ! correlation length
    Lc = Lx

    ! wave variance squared
    sigma2 = 1.0

    ! k-space grid spacing
    dk = 1.0 / Lc

    ! compute power spectrum envelope function
    Gnorm  = 0.0
    DO n = 1, Nmodes

        k = (n-0.5) * dk
        G(n) = (2.0*TWOPI*(k**2)*dk) / (1.0 + (k*Lc)**gamma)
        Gnorm = Gnorm + G(n)
    
    END DO    
    G = G / Gnorm
   
    ! sum over fourier modes to get the turbulent magnetic field
    Bout = 0.0
    Bout(3) = 1.0 ! mean magnetic field 
    
    DO n = 1, Nmodes
    
        k = (n-0.5) * dk

        ! compute random wave-vector, polarization vector and phase
        CALL RANDOM_NUMBER(rn)
        theta = 0.5 * TWOPI * rn(1)
        phi   = TWOPI * rn(2)
        alf   = TWOPI * rn(3)
        beta  = TWOPI * rn(4)
    
        khat(1) = SIN(theta) * COS(phi)
        khat(2) = SIN(theta) * SIN(phi)
        khat(3) = COS(theta)
        
        zeta(1) = -SIN(phi) * COS(alf) + COS(theta) * COS(phi) * SIN(alf)
        zeta(2) = COS(phi) * COS(alf) + COS(theta) * SIN(phi) * SIN(alf)
        zeta(3) = -SIN(theta) * SIN(alf)
    
        ! compute fourier mode amplitude
        amp = SQRT(sigma2*G(n))
    
        ! now sample the B-field at the input location
        kdotx = khat(1)*x(1) + khat(2)*x(2) + khat(3)*x(3)  
                    
        Bout(1) = Bout(1) + amp * zeta(1) * COS(k*kdotx + beta)     
        Bout(2) = Bout(2) + amp * zeta(2) * COS(k*kdotx + beta)     
        Bout(3) = Bout(3) + amp * zeta(3) * COS(k*kdotx + beta)     
  
    END DO
    



END SUBROUTINE generate_bfield_point





END MODULE generate_bfield_mod