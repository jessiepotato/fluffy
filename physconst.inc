c--
c-- updated to 2002 CODATA values by PHH 30/Jul/2007
c-- updated to 2010 CODATA values by PHH 11/Sep/2012
c
c--
c-- mathematical and physical constants
c--------------------------------------
c--
      real*8, parameter :: ev2erg=1.602 176 565d-12
      real*8, parameter :: c=2.99792458D+10
      real*8, parameter :: e=4.80325d-10
      real*8, parameter :: h=6.626 069 57d-27
      real*8, parameter :: kb=1.380 6488d-16
      real*8, parameter :: kev=kb/ev2erg
      real*8, parameter :: me=9.109 382 91 d-28
      real*8, parameter :: wlyman=911.7535d+0
      real*8, parameter :: mh=1.67333d-24
      real*8, parameter :: abohr=0.529 177 210 92d-8
      real*8, parameter :: whe2lyman=227.88486d0
      real*8, parameter :: pi=3.1415926535897932385d0
      real*8, parameter :: pi2=2.d0*pi
      real*8, parameter :: pi4=4.d0*pi
      real*8, parameter :: km2cm=1d+05
      real*8, parameter :: a2cm=1d-08
      real*8, parameter :: c1=h*c/ev2erg
      real*8, parameter :: rydberg=1.0 973 731 568 539d5
      real*8, parameter :: a_m_u=1.660 538 921d-24
      real*8, parameter :: Navog=6.022 141 29d23
      real*8, parameter :: stefco=5.670 373d-05
      real*8, parameter :: sigbolz=stefco
      real*8, parameter :: clight=c
      real*8, parameter :: boltz=kb
      real*8, parameter :: gravconst=6.673 84d-8
      real*8, parameter :: rg=8.314 4621 d+07
      real*8, parameter :: year=31556925.9747d0
      real*8, parameter :: hsun=3.826d33/(pi4*pi4)
      real*8, parameter :: msun=1.989d33
      real*8, parameter :: mspyr=msun/year
      real*8, parameter :: Rgas=kb/a_m_u   ! universal gas constant
      real*8, parameter :: theta_const=5.0397788085937500D+03   ! this is ~5040.d0: lg(e)*(ev2erg/k)
      real*8, parameter :: em = 2.*1.00794D0 ! molec wt H2 [g mol-1]
      real*8, parameter :: mu_amm = 1.7031d1 ! molecular weight of ammonia [g/mol]
