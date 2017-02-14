!---------------------------------------------------------------------!
!     CALCULATION STUFF
!---------------------------------------------------------------------!
      integer, parameter :: layers = 1400
      integer, parameter :: n = 25
      real*8, parameter :: Dz = 0.03 ! layer spacing [km]
!---------------------------------------------------------------------!
!     PLANET-SPECIFIC (currently, Jupiter)
!---------------------------------------------------------------------!
      real*8, parameter :: Tf = 1.24d2  ! effective temperature [K]
      real*8, parameter :: qbelow = 3.d-5 ! subcloud mole fraction (at 0.6 bars)
      real*8, parameter :: g = 25.d0 ! acceleration due to gravity (Jupiter) [m/s2]
      real*8, parameter :: mu = 2.2d0 ! atmospheric molecular weight [g/mol]
      real*8, parameter	:: cp = 7.*(Rgas*1E-7)/(2.*em) ! sp. heat @ constant P (ideal hydrogen) [J/K g]
      real*8, parameter :: rho_0 = 0.00016 ! rho_0 = atmospheric density at 1 bar [g/cm3]
!---------------------------------------------------------------------!
!     NOT SURE WHERE THESE BELONG
!---------------------------------------------------------------------!
      real*8, parameter :: eps = 59.7*kb ! depth of Lennard-Jones poten. well for atmo. [K]
      real*8, parameter :: d = 2.827d-8 ! molecular diameter of hydrogen [cm]
!---------------------------------------------------------------------!
!     ADJUSTABLE MODEL PARAMETERS
!---------------------------------------------------------------------!
      real*8, parameter :: frain = 3.d0
      real*8, parameter :: Kmin = 1.d5 ! [cm2/s]
      real*8, parameter :: A = 1.d-1
      real*8, parameter :: Sc = 0.d0
      real*8, parameter :: sig = 2.d0
!---------------------------------------------------------------------!
!     HOPEFULLY WILL APPLY TO ANY PLANET
!---------------------------------------------------------------------!
      real*8, parameter :: Tau_a = -g/cp ! [K/km]
      real*8, parameter :: rho_amm = 0.84 ! density of ammonia [g/cm3]
      real*8, parameter :: mu_a = 17.031 ! molecular wt of ammonia
