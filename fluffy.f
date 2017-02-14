      module interp
        implicit none
        include 'physconst.inc' !Rgas is in erg K-1 mol-1 (convert to J by x10^-7)
        include 'cloudvars.inc'
        include 'declaration.inc'
        contains
!----------------------------------------------------------------------!
!       Find temp at given alt by interpolation                        !
!----------------------------------------------------------------------!
      function TatZ(Pz,Pin,Tin) result(Tz)
      real*8, intent(in) :: Pz ! input
      real*8,intent(in), dimension (1 : 21) :: Pin, Tin ! input
      real*8 :: Tz ! output
      integer :: l ! only a location
      l = minloc(abs(Pin-Pz),1)
      if ( abs(Pin(l+1)-Pz) <= abs(Pin(l-1)-Pz) ) then
        Tz = Tin(l) + (Pz-Pin(l))*(Tin(l+1)-Tin(l))/(Pin(l+1)-Pin(l))
      elseif ( abs(Pin(l+1)-Pz) > abs(Pin(l-1)-Pz) ) then
        Tz = Tin(l-1)+(Pz-Pin(l-1))*(Tin(l)-Tin(l-1))/(Pin(l)-Pin(l-1))
      end if
      end function TatZ
!----------------------------------------------------------------------!
!       Find local temperature gradient                                !
!----------------------------------------------------------------------!
      function GRADI(Pz,Pin,Tin) result(tau)
      real*8, intent(in) :: Pz ! input
      real*8,intent(in), dimension (1 : 21) :: Pin, Tin ! input
      real*8 :: tau ! output
      integer :: l ! only a location
      l = minloc(abs(Pin-Pz),1)
      if ( abs(Pin(l+1)-Pz) <= abs(Pin(l-1)-Pz) ) then
        tau = (Tin(l+1)-Tin(l))/(Pin(l+1)-Pin(l))
      elseif ( abs(Pin(l+1)-Pz) > abs(Pin(l-1)-Pz) ) then
        tau = (Tin(l)-Tin(l-1))/(Pin(l)-Pin(l-1))
      end if
      end function GRADI
!----------------------------------------------------------------------!
!    Find eddy diffusion coefficient at given alt                      !
!    This is Equation 5 from Ackerman & Marley 2001                    !
!----------------------------------------------------------------------!
      function EDDY(H, len, ro) 	  result(K)
      real*8, intent(in) :: H, len, ro
      real*8, parameter :: F = (stefco*1.D-3)*(Tf)**4.        ! Flux [J/m2 s]
      real*8, parameter :: g_r = Rgas*1.D-7
      real :: K
      K = H*1.D5/3.*((len/H)**(4./3.))*(1.D3*g_r*F/(mu*ro*cp))**(1./3.)
!      write (20,*) Pz(i+1),len
      end function EDDY
!----------------------------------------------------------------------!
!       A function to calculate the saturation vapour mixing ratio for !
!       a given temperature and ambient pressure                       !
!----------------------------------------------------------------------!
      function SatVap(T,P)  result(Qs)
        real*8, intent(in) :: T, P  ! input temperature, pressure    !
        real*8 :: Qs    ! output is saturation vapour MR !
        qs = (exp(10.53 - 2161./T - 86596./T**2.))/P
      end function SatVap
!----------------------------------------------------------------------!
!       find sedimentation velocities                                  !
!----------------------------------------------------------------------!
      subroutine VSED(w, T, P, r, vs)
        real*8, intent(in) :: w  ! input mixing velocity
        real*8, dimension(1 : n), intent(out) :: vs, r ! output sedimentation velocity
        real*8 :: x,y,T,P,dr,dl
        dl = 1.d-6*Rgas*T/(sqrt(2.d0)*pi*d**2*Navog*P) ! molec. MFP in [cm]
        dr = rho_amm-rho_a ! density contrast btwn condensate & atm [g/cm^3]
        rw = (-1.26*dl+sqrt((1.26*dl)**2+4*w*eta*9./(2.*g*dr)))/2.
        do j=1,n
          r(j) = rw/sig + rw*(1.-1./sig)*(j-1.)/(n-1.)
          x = log(32.d2*rho_a*g*dr*r(j)**3/(3.*eta**2))
          y = .8*x-.01*x**2
            if ( dexp(y) <= 1000 ) then
            vs(j) = dexp(y)*eta/(2.*r(j)*rho_a)
            else
            vs(j) = (1+1.26*dl/r(j))*sqrt(8.d2*g*r(j)*dr/(1.35*rho_a))
            endif
        enddo
      end subroutine VSED
!----------------------------------------------------------------------!
      function ALPH(x, y) result(a)
         real*8, dimension(1 : n), intent(in) :: x, y
         real*8 :: a
      	a = size(x)*sum(log(x)*log(y))-sum(log(x))*sum(log(y))
        a = a/(size(x)*sum((log(x))**2)-(sum(log(x)))**2)
      end function ALPH
!----------------------------------------------------------------------!
      end module interp

      program fluffy
        use interp
        implicit none

      open (unit=20,file='pressure.txt',action="read",status="old")
      read(20, *) P;
      close(UNIT=20)
      open (unit=20,file='temp.txt',action="read",status="old")
      read(20, *) T;
      close(UNIT=20)

!     SET SOME INITIAL VALUES
      z(1) = 0.
      Tz(1) = 165.
      Pz(1) = 1.
      Hs = (Rgas*1.D-7)*Tf/(mu*g)!*Tz(1)/(mu*g)
      Qv(1) = qbelow
      Qc(1) = 0.d0
      Qt(1) = Qv(1) + Qc(1)
      dt(1) = 0.d0
!     NOW CALCULATE THE LOT
      open (unit=20,file="output.txt",action="write",status="replace")
      do i=1, layers-1
        z(i+1) = i*Dz
!       ESTIMATING PRESSURE & DENSITY USING HYDROSTATIC EQUILIBRIUM
        Pz(i+1) = exp(-z(i+1)/Hs)
        rho_a = rho_0*exp(-z(i+1)/HS)
!       GETTING TEMPERATURE FROM MEASURED PROFILE BY INTERPOLATION
        Tz(i+1) = TatZ(Pz(i+1),P,T)
!       FIND LOCAL SCALE HEIGHT ... makes the match to A&M better
        Hs = (Rgas*1.D-7)*Tz(i+1)/(mu*g)
!       CALCULATE LOCAL TEMPERATURE LAPSE RATE dT/dz
!       USING
!       dT/dz = [dT/dP]*[dP/dz] = [dT/dP]*[-rho_a * g]
        Tau = GRADI(Pz(i+1),P,T)*(-rho_a*g)*1.d1
!        Tau = 60.85*(-rho_a*g)*1.d1
!       CALCULATE LOCAL MIXING LENGTH
        L = Hs*max(A,Tau/Tau_a)
!       CALCULATE EDDY DIFFUSION COEFFICIENT
        K = EDDY(Hs,L,rho_a)
         if ( K < Kmin ) then
	        K = Kmin ! have prescribed a minimum value for eddy diff. coeff.
       	 endif
!       FIND THE CONVECTIVE VELOCITY FROM MIXING LENGTH THEORY
        w = (K/L)*(1.D-7)
!
!       AT LAST, THE GOOD BIT
!       SATURATION MIXING RATIO
        Qs(i+1) = SatVap(Tz(i+1),Pz(i+1))
!       ACTUAL MIXING RATIO OF VAPOUR => EQ 2 OF ACKERMAN & MARLEY 2001
        Qv(i+1) = min( Qv(i), (1.D0+Sc)*Qs(i+1) )
!       FIND AMOUNT OF "NEW" CONDENSATE BORN IN THIS LAYER
        cond = max( 0.d0, Qv(i)-(1.D0+Sc)*Qs(i+1) )
!       FIND RATE OF CHANGE OF MIXING RATIO WITH HEIGHT
        Dqt = -(frain*w*(Qc(i)+cond)/(K))*1.d4
!       UPDATE THE TOTAL MIXING RATIO
        Qt(i+1) = Qt(i) + Dqt*Dz*1.d3
!       FINALLY, SET THE NEW MIXING RATIO OF CONDENSATE
        Qc(i+1) = Qt(i+1)-Qv(i+1)
!        write (20,*) Pz(i+1),Tz(i+1),Qc(i+1),Qv(i+1),Qt(i+1),Qs(i+1)
!       NOW THE SECOND PART
        if ( Qc(i+1) > 0.d0 ) then
        eta = ((kb*Tz(i+1)/eps)**0.16)/(pi*(d**2)*1.22)
        eta = (5./16.)*sqrt(pi*em*kb*Tz(i+1)/Navog)*eta
        call VSED(w,Tz(i+1),Pz(i+1),r,vf)
        alf = ALPH(r/rw,vf/(w*1.d2))
!       FIND GEOMETRIC MEAN RADIUS
        r_g = rw*frain**(1/alf)*exp(-(alf+6)*(log(sig))**2/2.)
!       FIND TOTAL NUMBER CONCENTRATION OF PARTICLES
        EN = (3.*(mu_a/mu)*rho_a*Qc(i+1)/(4.*pi*rho_amm*r_g**3))
        EN = EN*exp(-9.*(log(sig))**2/2.)
!       FIND EFFECTIVE RADIUS
        r_eff = rw*frain**(1/alf)*exp(-(alf+1)*(log(sig))**2/2.)
!       ESTIMATE OPACITY FOR GEOMETRIC SCATTERER
        dt(i+1) = 3.*(mu_a/mu)*rho_a*Qc(i+1)*Dz/(2.*rho_amm*r_eff)
        write (20,*) Pz(i+1), Qc(i+1)!r_g*1.d4
        endif
      enddo
!
      close (20)
!      print*, 1.d2*8.,8.d2
      end program fluffy
