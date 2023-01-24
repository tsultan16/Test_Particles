MODULE mover_mod

USE constants_mod
USE grid_arrays_mod
USE generate_bfield_mod

IMPLICIT NONE

! rk parameters
REAL(4), PARAMETER :: A1 = .0, A2 = .2, A3 = .3, A4 = .6, A5 = 1., A6 = .875, &
                      B21 = .2, &
                      B31 = 3./40., B32 = 9./40. , &
                      B41 = .3, B42 = -.9, B43 = 1.2, &
                      B51 = -11./54., B52 = 2.5, B53 = -70./27., B54 = 35./27., &
                      B61 = 1631./55296., B62 = 175./512., B63 = 575./13824., B64 = 44275./110592., B65 = 253./4096., &
                      C1 = 37./378., C2 = 0.0, C3 = 250./621., C4 = 125./594., C5 = 0.0, C6 = 512./1771., &
                      CE1 = 2825./27648., CE2 = C2, CE3 = 18575./48384., CE4 = 13525./55296., CE5 = 277./14336., CE6 = .25 
                      
REAL(4), PARAMETER :: eps = 0.0001, pshrink = 0.25, safety = 0.9
LOGICAL, PARAMETER :: adaptive_dt = .TRUE.   ! switch for enabling adaptive step-size control
  

! particle charge to mass ratio                
REAL(4) :: my_qm

CONTAINS


SUBROUTINE move_particles(t, dt, qm, p, N)

    REAL(4) :: t, dt, qm
    TYPE(particle), INTENT(INOUT) :: p(:)
    INTEGER, INTENT(IN) :: N
    INTEGER :: i
    REAL(4) :: x(6)
        

    my_qm = qm


    ! loop over particles and integrate their equations of motion
    !$OMP PARALLEL PRIVATE(i,x) SHARED(p,t,dt,qm)
    !$OMP DO
    DO i = 1, N
        
        ! copy particle phase-space coordinates into buffer
        x(1:3) = p(i)%x        
        x(4:6) = p(i)%v        
    
        CALL integrate_particle(x, t, dt)        
    
        ! copy back from buffer into particle array
        p(i)%x = x(1:3)        
        p(i)%v = x(4:6)         
    
    
    END DO
    !$OMP END PARALLEL    


END SUBROUTINE move_particles



SUBROUTINE integrate_particle(x,t0,dt0)

    REAL(4), INTENT(INOUT) :: x(6)
    REAL(4), INTENT(IN) :: t0, dt0
    REAL(4) :: x0(6), xerr(6), xscal(6), del0(6), del1(6), dt, dtfac, t, tf 
    
    dt = dt0  
    t  = t0
    tf = t0 + dt0
        
    ! this loop does sub-cycling when adaptive step-size control requires smaller step-size
    DO WHILE(t .LT. tf)

        ! copy initial state into buffer
        x0 = x
    
        ! 5th order Runge-Kutta time integration
        CALL rk4(x0, x, del1, dt, t)
        
        ! apply adaptive step-size control
        IF(adaptive_dt) THEN
            
            !****************************
            ! adaptive step-size control
            !****************************
            xscal = ABS(x0) + ABS(dt * f(t,x0))  ! scale relative to required error tolerance 
            del0 = eps * xscal                   ! required error 
            
            dtfac = MINVAL(ABS(del0/del1))
                        
            IF(dtfac .LT. 1.0) THEN
            
                dt = MIN(safety * (dtfac**pshrink) * dt, tf-t)
                                
                IF(dt .LT. 1.e-15) THEN
                    PRINT*,'**WARNING**'
                    PRINT*,'Step size underflow...dt = ',dt
                END IF
                
                ! re-compute RK integration using new time step size
                CALL rk4(x0, x, xerr, dt, t)
            
            END IF
        
        END IF
        
        t = t + dt
        dt = MIN(dt, tf-t)

    END DO


END SUBROUTINE integrate_particle



! advances solution over time step dt using 4th order (Cash-Karp) embedded Runge-Kutta and computes the truncation error 
SUBROUTINE rk4(x0, xout, xerr, dt, t0)

    REAL(4), INTENT(IN) ::  dt, t0     
    REAL(4), INTENT(INOUT) :: x0(6), xout(6), xerr(6) 
    INTEGER :: i
    REAL(4) :: dt2, dxt(6,6), xt(6), xt4(6)
    
    ! RK step 1
    dxt(1,1:6) = f(t0+A1*dt,x0)
    DO i = 1, 6
        xt(i) = x0(i) + dt*B21*dxt(1,i)
    END DO    
    CALL periodic_bc(xt(1:3))

    !RK step 2
    dxt(2,1:6) = f(t0+A2*dt,xt)
    DO i = 1, 6
        xt(i) = x0(i) + dt * (B31*dxt(1,i) + B32*dxt(2,i))
    END DO
    CALL periodic_bc(xt(1:3))
      
    !RK step 3
    dxt(3,1:6) = f(t0+A3*dt,xt)
    DO i = 1, 6
        xt(i) = x0(i) + dt * (B41*dxt(1,i) + B42*dxt(2,i) + B43*dxt(3,i))
    END DO
    CALL periodic_bc(xt(1:3))
   
    !RK step 4
    dxt(4,1:6) = f(t0+A4*dt,xt)
    DO i = 1, 6
        xt(i) = x0(i) + dt * (B51*dxt(1,i) + B52*dxt(2,i) + B53*dxt(3,i) + B54*dxt(4,i))
    END DO
    CALL periodic_bc(xt(1:3))
    
    !RK step 5
    dxt(5,1:6) = f(t0+A5*dt,xt)
    DO i = 1, 6
        xt(i) = x0(i) + dt * (B61*dxt(1,i) + B62*dxt(2,i) + B63*dxt(3,i) + B64*dxt(4,i) + B65*dxt(5,i))
    END DO
    CALL periodic_bc(xt(1:3))
    
    !RK step 6 (final accumulation)
    dxt(6,1:6) = f(t0+A6*dt,xt)
    DO i = 1, 6
        xout(i) = x0(i) + dt * (C1*dxt(1,i) + C2*dxt(2,i) + C3*dxt(3,i) + C4*dxt(4,i) + C5*dxt(5,i) + C6*dxt(6,i))       ! embedded 4th order RK
        xt4(i)  = x0(i) + dt * (CE1*dxt(1,i) + CE2*dxt(2,i) + CE3*dxt(3,i) + CE4*dxt(4,i) + CE5*dxt(5,i) + CE6*dxt(6,i)) ! 5th order RK
    END DO    
        
    ! compute trunctaion error
    xerr = xout - xt4
    
    CALL periodic_bc(xout(1:3))
    
   
END SUBROUTINE rk4



FUNCTION f(t,x) RESULT(ftx)

    REAL(4), INTENT(IN) :: t, x(6)
    REAL(4) :: ftx(6), Bint(3),a(3), xp(3), xc(3), Sint(-1:2,3)
    INTEGER :: ix, iy, iz, i, j, k


    xp(1:3) = x(1:3)
    CALL periodic_bc(xp)
    
    
    IF(fieldInterp) THEN
    
        !************************************************
        ! interpolate magnetic field to particle location
        !************************************************
        
        ix = xp(1)/dx + 0.5
        iy = xp(2)/dx + 0.5
        iz = xp(3)/dx + 0.5   
        
        ! 1st order interpolation co-efficients
        Sint = 0.0
        
        !DO i = 0,1
        
        !    xc(1) = ((ix+i)-0.5) - x(1)/dx
        !    xc(2) = ((iy+i)-0.5) - x(2)/dx
        !    xc(3) = ((iz+i)-0.5) - x(3)/dx
                   
        !    Sint(i,1) = S1_1D(xc(1))
        !    Sint(i,2) = S1_1D(xc(2))
        !    Sint(i,3) = S1_1D(xc(3))
            
        !END DO       
        
        !IF(ii .EQ. 8279) THEN ! .AND. tstep .EQ. 14) THEN
        !PRINT*,'my_qm = ',my_qm
        !PRINT*,'xp = ',xp(1),xp(2),xp(3)
        !PRINT*,'ix = ',ix,iy,iz
        !END IF
        !PRINT*,'Sint_x(0), Sint_x(1) = ',Sint(0,1), Sint(1,1)    
        
        ! 2nd order interpolation co-efficients
        DO i = -1,2
            xc(1) = ((ix+i)-0.5) - xp(1)/dx
            xc(2) = ((iy+i)-0.5) - xp(2)/dx
            xc(3) = ((iz+i)-0.5) - xp(3)/dx
            Sint(i,1) = S1_1D(xc(1))
            Sint(i,2) = S1_1D(xc(2))
            Sint(i,3) = S1_1D(xc(3))    
        END DO 
               
        ! now interpolate from grid to particle location
        Bint = 0.0
        DO k = -1,2
            DO j = -1,2
                DO i = -1,2
                    Bint(1) = Bint(1) + (Sint(i,1) * Sint(j,2) * Sint(k,3) * B(ix+i,iy+j,iz+k,1)) 
                    Bint(2) = Bint(2) + (Sint(i,1) * Sint(j,2) * Sint(k,3) * B(ix+i,iy+j,iz+k,2)) 
                    Bint(3) = Bint(3) + (Sint(i,1) * Sint(j,2) * Sint(k,3) * B(ix+i,iy+j,iz+k,3)) 
                END DO
            END DO
        END DO       

        !PRINT*,'Interpolated B = ',Bint

    ELSE

        CALL generate_bfield_point(xp,Bint)

    END IF    

        
    ! compute acceleration (Lorentz Force)
    a(1) = my_qm * (Bint(3)*x(5) - Bint(2)*x(6))
    a(2) = my_qm * (Bint(1)*x(6) - Bint(3)*x(4))
    a(3) = my_qm * (Bint(2)*x(4) - Bint(1)*x(5))

    ftx(:) = (/ x(4), x(5), x(6), a(1), a(2), a(3) /)


END FUNCTION f


! 1D particle 1st order shape function (piecewise linear)
FUNCTION S1_1D(x) RESULT (sx)

    REAL(4), INTENT(IN) :: x
    REAL(4) :: sx

    IF(ABS(x) .LT. 1.0) THEN
        sx = 1.d0 - ABS(x)     
    ELSE    
        sx = 0.d0
    END IF

END FUNCTION S1_1D


! 1D particle 2nd order shape function (quadratic spline)
FUNCTION S2_1D(x) RESULT (sx)

    REAL(4), INTENT(IN) :: x
    REAL(4) :: sx

    IF(ABS(x) .LT. 0.5d0) THEN
        sx = (3.d0/4.d0) - x**2     
    ELSE IF(ABS(x) .GE. 0.5d0 .AND. ABS(x) .LE. 1.5d0)THEN
        sx = 0.5d0 * ( 1.5d0 - ABS(x))**2
    ELSE    
        sx = 0.d0
    END IF

END FUNCTION S2_1D


! this routine applies periodic boundary conditions to particles
SUBROUTINE periodic_bc(x)

    REAL(4), INTENT(INOUT) :: x(3) 
        
    !x(1) = MOD(x(1),Lx)
    !x(2) = MOD(x(2),Ly)
    !x(3) = MOD(x(3),Lz)
    
    IF(x(1) .LT. 0.0) x(1) = x(1) + Lx
    IF(x(1) .GT. Lx) x(1) = x(1) - Lx

    IF(x(2) .LT. 0.0) x(2) = x(2) + Ly
    IF(x(2) .GT. Ly) x(2) = x(2) - Ly

    IF(x(3) .LT. 0.0) x(3) = x(3) + Lz 
    IF(x(3) .GT. Lz) x(3) = x(3) - Lz

    
END SUBROUTINE periodic_bc


END MODULE mover_mod 