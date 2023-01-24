! ODE integrator
! Reference: Numerical Recipes, Press et al.
!

PROGRAM integrator_mod

USE OMP_LIB 

IMPLICIT NONE


INTEGER, PARAMETER :: n = 6 ! number of dependent variables
INTEGER, PARAMETER :: nsteps = 100000
REAL(4), PARAMETER :: tend = 3000.0
REAL(4), PARAMETER :: eps = 0.05
INTEGER, PARAMETER :: KMAX = 8
REAL(4), PARAMETER :: a = 1.0, pgrow = 0.2, pshrink = 0.25, safety = 0.98
LOGICAL, PARAMETER :: adaptive_dt = .TRUE.
REAL(4) :: x(1:n), xscal(1:n), xerr(1:n), x1(1:n), x2(1:n), dt, dtnext, del1(1:n), del0(1:n), xexact(1:n)
REAL(4) :: x0(1:n), dxdt0(1:n),  t0, t, dtfac
INTEGER :: i
REAL(4) :: Amp, w, delta, t1, t2


! initial state
t0 = 0.0
x = (/ 0.0, 1.0, 0.0, 1.0, 0.0, 0.0 /)
x1 = x
dxdt0 = f(t0,x0)
dt = 0.15
t = t0
i = 1


! exact solution constant
!w = SQRT(16.0) 
delta = 2.0*ATAN(1.d0)
!Amp = x(1) / COS(delta) 


PRINT*,''
PRINT*,'x0 = ',x0
PRINT*,'dxdt0 = ',dxdt0
PRINT*,''

OPEN(UNIT = 11, FILE = 'rk4_out.dat', FORM = 'UNFORMATTED', ACCESS = 'STREAM')


t1 = OMP_GET_WTIME()

DO WHILE(t < tend .AND. i < nsteps)

    !WRITE(*,FMT='("Time step # ",i5,", dt = ",f9.2,", t = ",f9.2)') i, dt, t
    
    ! single full step
    !x0 = x1
    !CALL rk4(x0,x1,dt,t)
    !x1 = 0.0
    
    ! double half-steps
    !CALL rk4(x0,x2,0.5*dt,t)
    !x0 = x2 
    !CALL rk4(x0,x2,0.5*dt,t+0.5*dt)
    
    x0 = x
    CALL rk5(x0, x, del1, dt, t)
    
    IF(adaptive_dt .AND. i .GT. 2) THEN
        
        !****************************
        ! adaptive step-size control
        !****************************
        xscal = ABS(x0) + ABS(dt * f(t,x0)) ! scale relative to required error tolerance 
        del0 = eps * x0                      ! required error 
        
        dtfac = MINVAL(ABS(del0/del1))
        !PRINT*,'dt factor = ',dtfac
        
        IF(dtfac .GE. 1.0) THEN
            dt = safety * (dtfac**pgrow  ) * dt 
        ELSE
            dt = safety * (dtfac**pshrink) * dt
        END IF
        
        
        !PRINT*,'Desired error = ',del0
        !PRINT*,'Truncation Error = ',del1
        !PRINT*,'dtnext = ',dtnext
        
        IF(dt .LT. 1.e-15) THEN
            PRINT*,'**WARNING**'
            PRINT*,'Step size underflow...dtnext = ',dtnext
        END IF
        
        ! re-compute RK integration using new time step size
        CALL rk5(x0, x, xerr, dt, t)
    
    END IF
    
    CALL mmid(4,x0,x1,dt,t)
    
    t = t + dt
    i = i + 1
    
    ! exact solution
    xexact(1) = SIN(t+delta)
    xexact(2) = COS(t+delta)
    xexact(3) = 0.0
    xexact(4) = COS(t+delta)
    xexact(5) = -SIN(t+delta)
    xexact(6) = 0.0
    
    
    !PRINT*,'Exact truncation error (RK5) = ',SQRT(SUM((xexact(1:2)-x(1:2))**2))
    !PRINT*,'Exact truncation error (MMID) = ',SQRT(SUM((xexact(1:2)-x1(1:2))**2))
    
    !save to file
    WRITE(11) t,x(1),x(2),x(3),x(4),x(5),x(6)
   
   
END DO

PRINT*,'Done.'
PRINT*,'# of steps = ',i

t2 = OMP_GET_WTIME()

PRINT*,'t1, t2 = ',t1,t2
PRINT*,'Total time (sec) = ',t2-t1
PRINT*,'Time per step (sec) = ',(t2-t1)/i
PRINT*,''


CLOSE(UNIT = 11)




CONTAINS



! Given the initial values (t = t0) of the vector of dependent variables x(n) and their derivatives dxdt(n), this subroutine advances the solution
! over a step size dt using the Bulirsch-Stoer Method.

!############################## INCOMPLETE  #########################################

SUBROUTINE bulirsch_stoer(x0, xout, xerr, xscal, dt, t0)

    REAL(4), INTENT(IN) :: x0(n), xscal(n), t0     
    REAL(4), INTENT(INOUT) :: xout(n), xerr(n), dt 
    INTEGER :: i, k, km
    REAL(4) :: x_extrp(1:n), xn(1:n), errk(1:KMAX), h, errmax
    INTEGER :: nseq
    REAL(4), PARAMETER :: SAFE1 = 0.25
    REAL(4), PARAMETER :: SAFE2 = 0.7

    x_extrp = x0

    ! Loop over sequence of even numbers (nseq = 2,4,6,...,2*KMAX). For each, we compute the solution using the modified midpoint method
    ! using sub-step size h = dt/nseq. We use the sequence of solutions that have been computed so far and use polynomial extrapolation
    ! to estimate a zero-step size solution and the corresponding error and check that against tolerance to decide whether to stop or proceed to the next
    ! even number in the sequence.    
    DO k = 1, KMAX

            nseq = 2*k    

            ! compute modified midpoint solution  
            CALL mmid(nseq,x0,xn,dt,t0)

            h = (dt/SNGL(nseq))**2  ! squared since error series only has even powers of h    
            
            ! do polynomial extrapolation to approximate zero step-size solution and compute error estimate
            CALL pol_extrp(k,h,xn,xerr,x_extrp)
    
            IF(k .NE. 1) THEN
                errmax = 1.e-30 ! set this to a small positive number
                errmax = MAX(errmax, MAXVAL(ABS(xerr/xscal))) / eps
                km = k - 1
                errk(km) = (errmax/SAFE1)**(1.0/(2.0*km+1.0))
            END IF
    
            !IF((k .NE. 1) .AND. ((k .GE. kopt-1) .OR. first)) THEN
                
           
            !END IF
            
            xout = xn
    
    END DO



END SUBROUTINE bulirsch_stoer


! uses polynomial extrapolation to compute 'n' functions at 'h = 0' by fitting a polynomial through a sequence of 'iest'
! points 'xest' with progressivele smaller values 'h=hest'. The ezxtrapolated values are stored in 'x_extrp' and estimated error is stored in 'xerr'  
SUBROUTINE pol_extrp(iest,hest,xest,xerr,x_extrp)

    INTEGER, INTENT(IN) :: iest
    REAL(4), INTENT(IN) :: hest, xest(1:n)
    REAL(4), INTENT(OUT) :: xerr(1:n), x_extrp(1:n)
  
    INTEGER :: j, k1
    REAL(4) :: delta, f1, f2, q, d(1:n)
    REAL(4), SAVE :: qcol(1:n,1:KMAX), h(1:KMAX)


    ! save current values   
    h(iest) = hest
    xerr(1:n) = xest(1:n)
    x_extrp(1:n) = xest(1:n)
    
    IF(iest .EQ. 1) THEN
       qcol(1:n,iest) = xest(1:n)   ! store first estimate in first column
    ELSE
    
        d(1:n) = xest(1:n)
        
        DO k1 = 1, iest-1
        
            delta = 1.0/(h(iest-k1)-hest)
            f1 = hest*delta
            f2 = h(iest-k1)*delta
            
            ! propagate tableau 1 diagonal more
            DO j = 1, n
                q = qcol(j,k1)        
                qcol(j,k1) = xerr(j)
                delta = d(j) - q
                xerr(j) = f1*delta
                d(j) = f2*delta
                x_extrp(j) = x_extrp(j) + xerr(j)
            END DO            
        END DO
        
        qcol(1:n,iest) = xerr(1:n)
    
    END IF        
       

END SUBROUTINE pol_extrp


! Given the initial values (t = t0) of the vector of dependent variables x(n) and their derivatives dxdt(n), this subroutine advances the solution
! over a step size dt using the 4th order Runge-Kutta method.
SUBROUTINE rk4(x0, xout, dt, t0)

    REAL(4), INTENT(IN) :: x0(n), dt, t0     
    REAL(4), INTENT(OUT) :: xout(n) 
    INTEGER :: i
    REAL(4) :: dt2, dxt(4,n), xt(n)
    REAL(4) :: A1, A2, A3, A4, B21, B31, B32, B41, B42, B43, C1, C2, C3, C4
    PARAMETER ( A1 = .0, A2 = .5, A3 = .5, A4 = 1. , &
                B21 = .5, &
                B31 = .0, B32 = .5 , &
                B41 = .0, B42 = .0, B43 = 1., &
                C1 = 1./6., C2 = 1./3., C3 = 1./3., C4 = 1./6 )
   
    dt2 = dt/2.0
    
    ! RK step 1
    dxt(1,1:n) = f(t0+A1*dt,x0)
    DO i = 1, n
        xt(i) = x0(i) + dt*B21*dxt(1,i)
    END DO   

    !RK step 2
    dxt(2,1:n) = f(t0+A2*dt,xt)
    DO i = 1, n
        xt(i) = x0(i) + dt * (B31*dxt(1,i) + B32*dxt(2,i))
    END DO
    
    !RK step 3
    dxt(3,1:n) = f(t0+A3*dt,xt)
    DO i = 1, n
        xt(i) = x0(i) + dt * (B41*dxt(1,i) + B42*dxt(2,i) + B43*dxt(3,i))
    END DO
   
    !RK step 4 (final)
    dxt(4,1:n) = f(t0+A4*dt,xt)
    DO i = 1, n
        xout(i) = x0(i) + dt * (C1*dxt(1,i) + C2*dxt(2,i) + C3*dxt(3,i) + C4*dxt(4,i))  
    END DO   
    
    
END SUBROUTINE rk4



! advances solution by over time step dt using 5th order (cash-Karp) Runge Kutta and also computes the truncation error using the embedded 4th order Runge-Kutta
SUBROUTINE rk5(x0, xout, xerr, dt, t0)

    REAL(4), INTENT(IN) :: x0(n), dt, t0     
    REAL(4), INTENT(OUT) :: xout(n), xerr(n) 
    INTEGER :: i
    REAL(4) :: dt2, dxt(6,n), xt(n), xt4(n)
    REAL(4) :: A1, A2, A3, A4, A5, A6, B21, B31, B32, B41, B42, B43, B51,B52,B53,B54, B61, B62, B63, B64, B65, C1, C2, C3, C4, C5, C6, CE1, CE2, CE3, CE4, CE5, CE6
    PARAMETER ( A1 = .0, A2 = .2, A3 = .3, A4 = .6, A5 = 1., A6 = .875, &
                B21 = .2, &
                B31 = 3./40., B32 = 9./40. , &
                B41 = .3, B42 = -.9, B43 = 1.2, &
                B51 = -11./54., B52 = 2.5, B53 = -70./27., B54 = 35./27., &
                B61 = 1631./55296., B62 = 175./512., B63 = 575./13824., B64 = 44275./110592., B65 = 253./4096., &
                C1 = 37./378., C2 = .0, C3 = 250./621., C4 = 125./594.,C6 = 512./1771., &
                CE1 = 2825./27648., CE2 = C2, CE3 = 18575./48384., CE4 = 13525./55296., CE5 = 277./14336., CE6 = .25 )
    
    ! RK step 1
    dxt(1,1:n) = f(t0+A1*dt,x0)
    DO i = 1, n
        xt(i) = x0(i) + dt*B21*dxt(1,i)
    END DO
    

    !RK step 2
    dxt(2,1:n) = f(t0+A2*dt,xt)
    DO i = 1, n
        xt(i) = x0(i) + dt * (B31*dxt(1,i) + B32*dxt(2,i))
    END DO
    
    !RK step 3
    dxt(3,1:n) = f(t0+A3*dt,xt)
    DO i = 1, n
        xt(i) = x0(i) + dt * (B41*dxt(1,i) + B42*dxt(2,i) + B43*dxt(3,i))
    END DO
   
    !RK step 4
    dxt(4,1:n) = f(t0+A4*dt,xt)
    DO i = 1, n
        xt(i) = x0(i) + dt * (B51*dxt(1,i) + B52*dxt(2,i) + B53*dxt(3,i) + B54*dxt(4,i))
    END DO
    
    !RK step 5
    dxt(5,1:n) = f(t0+A5*dt,xt)
    DO i = 1, n
        xt(i) = x0(i) + dt * (B61*dxt(1,i) + B62*dxt(2,i) + B63*dxt(3,i) + B64*dxt(4,i) + B65*dxt(5,i))
    END DO
    
    
    !RK step 6 (final accumulation)
    dxt(6,1:n) = f(t0+A6*dt,xt)
    DO i = 1, n
        xout(i) = x0(i) + dt * (C1*dxt(1,i) + C2*dxt(2,i) + C3*dxt(3,i) + C4*dxt(4,i) + C5*dxt(5,i) + C6*dxt(6,i))       ! 5th order RK
        xt4(i)  = x0(i) + dt * (CE1*dxt(1,i) + CE2*dxt(2,i) + CE3*dxt(3,i) + CE4*dxt(4,i) + CE5*dxt(5,i) + CE6*dxt(6,i)) ! embedded 4th order RK
    END DO    
    
    
    ! compute trunctaion error
    xerr = xt - xt4
   
   
END SUBROUTINE rk5


! Modified midpoint method
SUBROUTINE mmid(nstep, x0, xout, dt, t0)

    REAL(4), INTENT(IN) :: x0(n), dt, t0     
    REAL(4), INTENT(OUT) :: xout(n) 
    INTEGER, INTENT(IN) :: nstep
    INTEGER :: i
    REAL(4) :: h, x1(n), x2(n), xn(n), dxdt(n)

   
    ! sub step size
    h = dt/ nstep
    
    ! initilaze the first two points
    x1(1:n) = x0 
    x2(1:n) = x0 + h * f(t0,x0)
    
    
    ! compute intermediate points
    DO i = 1, nstep-1
        dxdt(1:n) = f(t0+i*h,x2)  
        xn(1:n) = x1(1:n) + 2.0*h * dxdt(1:n)  
        x1 = x2
        x2 = xn
    END DO
    
    ! compute final approximation
    dxdt(1:n) = f(t0+dt,x2)
    xout(1:n) = 0.5 * (x2(1:n) + x1(1:n) + h*dxdt(1:n))
    
END SUBROUTINE mmid       


FUNCTION f(t,x) RESULT(ftx)

    REAL(4), INTENT(IN) :: t, x(n)
    REAL(4) :: ftx(n)
    REAL(4), PARAMETER :: k = 16


    ! uniform acceleration 1D
    !ftx(:) = (/ x(2), a /)

    ! uniform magnetic field (magnitude = 1) in z direction
    ftx(:) = (/ x(4), x(5), x(6), x(5), -x(4), 0.0 /)


END FUNCTION f



END PROGRAM integrator_mod

