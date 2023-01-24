MODULE constants_mod

IMPLICIT NONE

! particle parameters
INTEGER, PARAMETER :: Ne_max = 100000    ! max number of electrons 
INTEGER, PARAMETER :: Np_max = 100000    ! max number of protons


! grid parameters
INTEGER, PARAMETER :: nx = 256
INTEGER, PARAMETER :: ny = nx
INTEGER, PARAMETER :: nz = nx
INTEGER, PARAMETER :: nb = 2


! particle parameters
REAL(4), PARAMETER  :: vmin  = -3.0      ! lower cut-off for velocity space distribution
REAL(4), PARAMETER  :: vmax  = 3.0       ! upper cut-off for velocity space distribution
REAL(4), PARAMETER  :: v0x   = 0.1      ! intial drift velocity
REAL(4), PARAMETER  :: v0y   = 0.0 
REAL(4), PARAMETER  :: v0z   = 0.0 
REAL(4), PARAMETER  :: vth_e  = 0.0   !0.1d0     ! R.M.S. thermal velocity
REAL(4), PARAMETER  :: vth_p  = 0.0   ! R.M.S. thermal velocity
INTEGER, PARAMETER :: vbins = 1000        ! velocity space bins 


! physics parameters
REAL(4), PARAMETER :: eps0 = 1.0        ! vacuum dielectric permittivity
REAL(4), PARAMETER :: c    = 0.45       ! speed of light in vacuum (also serves as the Courant number, for stability set this to be slightly <0.5)
REAL(4), PARAMETER :: q_e  = -1.e-4     ! electron charge
REAL(4), PARAMETER :: q_p  = -q_e       ! proton charge
REAL(4), PARAMETER :: m_e = 1.0*q_p     ! electron mass
REAL(4), PARAMETER :: m_p = 1.0*m_e     ! proton mass
REAL(4), PARAMETER :: wc_e = 0.0  	    ! electron cyclotron frequency
REAL(4), PARAMETER :: wc_p = 0.0        ! proton cyclotron frequency
 

! field interpolation parameters
LOGICAL, PARAMETER :: fieldInterp = .TRUE.

! Misc
REAL(4), PARAMETER :: TWOPI = 8.d0*ATAN(1.d0)


END MODULE constants_mod