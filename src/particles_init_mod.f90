MODULE particles_init_mod

USE constants_mod
USE grid_arrays_mod

IMPLICIT NONE



CONTAINS


! routine for initializing particle states at the beginning of simulation
SUBROUTINE particles_init()

    INTEGER :: i, j, k, ilow, ihi, jlow, jhi, klow, khi, nd, nn
    REAL(4) :: p(3), x(3), ve(3), vp(3)


    ! specify region to populate with particles
    ilow = 0.5*nx
    ihi = ilow+10
    jlow = ilow
    jhi = ihi
    klow = nz/2
    khi = nz/2
    
    ! specify particle number density
    nd = 10

    ! place particles uniformly across specified region of grid with approximately Maxwellian velocity distribution
    DO k = klow, khi, 1 
        DO j = jlow, jhi, 1
            DO i = ilow, ihi, 1
                DO nn = 1, nd
            
                    !IF(Ne .GE. 1 .OR. Np .GE. 1) EXIT
            
                    CALL RANDOM_NUMBER(p)                                        
                    p = p - 0.5
                    x(1) = (i - 0.5 + p(1)) * dx
                    x(2) = (j - 0.5 + p(2)) * dx 
                    x(3) = (k - 0.5) * dx       !+ p(3) 
            
                    CALL RANDOM_NUMBER(p)                                        
                    p = p - 0.5
                    ve(1) = v0x + vth_e * (vth_e + SUM(p))
                                                
                    CALL RANDOM_NUMBER(p)                                       
                    p = p - 0.5
                    ve(2) = vth_e * SUM(p)

                    CALL RANDOM_NUMBER(p)                                       
                    p = p - 0.5                
                    ve(3) = vth_e * SUM(p)
                               
                    CALL create_particle(electron,Ne,x,ve)    
                           
                    CALL RANDOM_NUMBER(p)                                      

                    p = p - 0.5

                    vp(1) = v0x + vth_p*(vth_p + SUM(p))
                        
                    CALL RANDOM_NUMBER(p)                                      
                    p = p - 0.5
                    vp(2) = vth_p * SUM(p)

                    CALL RANDOM_NUMBER(p)                    
                    p = p - 0.5                  
                    vp(3) = vth_p * SUM(p)          
                   
                    CALL create_particle(proton,Np,x,vp)
           
                               
                END DO
            END DO
        END DO
	END DO


    PRINT*,''
    PRINT*,'Total number of electrons created = ',Ne
    PRINT*,'Total number of protons created = ',Np
    PRINT*,''


END SUBROUTINE particles_init



! routine for inserting a new particle into the particle array
SUBROUTINE create_particle(p, N, x_in, v_in)

    TYPE(particle), INTENT(INOUT) :: p(:)
    INTEGER, INTENT(INOUT) :: N
    REAL(4), INTENT(IN) :: x_in(3), v_in(3)
    
    ! incerement particle counter
    N = N + 1
    
    ! initilize particle state
    p(N)%x = x_in
    p(N)%v = v_in


END SUBROUTINE create_particle


! routine for removing an existing particle from the particle array
SUBROUTINE destroy_particle(p, N, x_in, v_in)

    TYPE(particle), INTENT(INOUT) :: p(:)
    INTEGER, INTENT(INOUT) :: N
    INTEGER, INTENT(IN) :: x_in(3), v_in(3)
    
    ! decerement particle counter
    N = N - 1
    
    ! no need to clear out the values stored in the empty slot since it'll just get overwritten if a new particle gets created later
    

END SUBROUTINE destroy_particle




END MODULE particles_init_mod 