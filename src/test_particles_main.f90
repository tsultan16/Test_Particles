PROGRAM test_particles_main

USE constants_mod
USE grid_arrays_mod
USE mover_mod
USE particles_init_mod
USE file_io_mod
USE generate_bfield_mod
USE OMP_LIB

IMPLICIT NONE


INTEGER, parameter :: nsteps = 500
REAL(4) :: dt, tsim, qm
REAL(8) :: t1, t2, t3, t4 
INTEGER :: i, j, k

! initialze grid arrays
PRINT*,'Initializing grid...'
CALL create_grid_arrays()

! initialze partiles
PRINT*,'Initializing particles...'
CALL particles_init()
CALL file_output(0)

! initialize magnetic field
PRINT*,'Initializing magnetic field...'
B = 0.0
!B(:,:,:,3) = 1.0  ! uniform B field along z-direction

! linearly varying field
!DO k = 1-nb, nz+nb
!    DO j = 1-nb, ny+nb
!        DO i = 1-nb, nx+nb

!            B(i,j,k,3) = 1.0 + 0.8*j*dx

!        END DO
!    END DO
!END DO


tsim = 0.0

t1 = OMP_GET_WTIME()

! begin simulation loop

PRINT*,'Beginning simulation loop...'

DO i = 1, nsteps

   tstep = i

   ! set the time step size 
   dt = 0.25

   !PRINT*,'Moving particles...'

   ! advance particles by time step.. electrons first, then protons
   qm = q_e/m_e
   CALL move_particles(tsim, dt, qm, electron, Ne)

   qm = q_p/m_p
   CALL move_particles(tsim, dt, qm, proton, Np)

   ! update simulation time
   tsim = tsim + dt

   ! save to file
   !PRINT*,'Saving to file...'
   CALL file_output(i)

   WRITE(*,'(" Time Step # = ",i6, ", dt = ", e9.2, ", Simulation clock time = ", f10.3)') i,dt,tsim

END DO

t2 = OMP_GET_WTIME()

CALL destroy_grid_arrays()

PRINT*,''
PRINT*,'Simulation Loop Wall Time (sec) = ',t2-t1
PRINT*,'Done.'




CONTAINS





END PROGRAM test_particles_main