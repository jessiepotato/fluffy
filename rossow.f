c***************************************************************************
c    Grain diffusion routine by Tristan Guillot <guillot@obs-nice.fr>
c                   From - Fri Apr 21 17:21:33 2000
c      This routine computes new abundances of condensates and their
c      averaged size. It returns right away when the initial abundance
c      of condensate is <=0.
c      Characteristic mixing time (tmix) must be between:
c                    1.d3 sec (nothing settles)
c               and  1.d4 sec (practically all grains settle)
c***************************************************************************


      subroutine rossow_grains(p,t,rho,grav,hp,tmix,qcond,rog,
     &                         qleft,amoy)
c     ------------------------------------------------------------
c     Version: 18/04/00
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Code de calcul de temps caracteristiques associes a la formation c
c     de nuages & grains dans les atmospheres planetaires et stellairesc
c     D'apres Rossow 1978                                              c
c     Note: On tient compte uniquement d'un equilibre entre            c
c      melange turbulent, condensation et sedimentation                c
c      (tps de coagulation & coalescence negliges)                     c
ccccccccccccccccccccccccsccccccccccccccccccccccccccccccccccccccccccccccc
c     (c) T. Guillot, Laboratoire Cassini, OCA, Nice
c         guillot@obs-nice.fr
c Entrees:
c     t,p,rho,grav,hp: temperature, pression, densite, gravite,
c                  echelle de hauteur de pression; tout en cgs
c     tmix: temps de melange turbulent (en secondes)
c     qcond: fraction molaire initiale totale des elements condenses
c     rog: densite des grains condenses (moyenne) (en g/cm3)
c Sorties:
c     qleft: fraction molaire des grains tel que
c            tmix <= (tfall=tcondensation)
c     amoy: taille moyenne des grains (tel que tcond=tmix) (en cm)

      implicit none
      real*8 t,p,rho,grav,hp,tmix,qcond,rog,qleft,amoy

c Variables locales
      integer i,itype,j,k,ier
      integer*4 leng
      real*8 minx,maxx,tcond,tfall,tcoag,tcoal,tgrow
      data minx/-8.d0/,maxx/0.d0/

c Fonctions externes
      real*8 ffall_mix,fgrow_mix,f
      external ffall_mix,fgrow_mix,f

c Common: pour appel de fonction uniquement (for function call only)
      real*8 rp,rt,rro,rgrav,rqleft,rrog,rtmix
      common /rossow/rp,rt,rro,rgrav,rqleft,rrog,rtmix

c----------------------------------------------------------------------

      if (qcond.le.0.d0) then
       qleft=0.d0
       amoy=0.d0
       return
      endif

      qleft=qcond

      rp=p
      rt=t
      rro=rho
      rgrav=grav
      rqleft=qleft
      rrog=rog
      rtmix=tmix

c On calcule la taille d'equilibre <a> melange/sedimentation
c Valeurs maximales de log10(size): minx(-8) et maxx (0) (en cm)
      call secante(ffall_mix,minx,maxx,amoy,1d-4,ier)
c      write(*,*)'Fall/Mix: amoy,ier',amoy,ier
      if ((ier.eq.1).and.(amoy.le.minx)) then
       qleft=0.d0
       amoy=0.d0
       return
      endif
c Tps de croissance pour cette taille de particules
      amoy=10.d0**amoy
      call rossow78_mod(p,t,rho,grav,qleft,rog,
     &      amoy,tcond,tfall,tcoag,tcoal)
      tgrow=1/(1.d0/tcond+1.d0/tcoag+1.d0/tcoal)
c      write(*,*)'tgrow,tcond:',tgrow,tcond,tcoag,tcoal,tfall
c On compare le tps de melange au temps de condensation pour la taille <a>
c Si la condensation est plus rapide, on recalcule la fraction de
c condenses pour obtenir tcon=tmix
      if (tgrow.lt.tmix) then
       qleft=qleft*tgrow/tmix
       do while (abs(qleft-rqleft).gt.1d-10)
        call rossow78_mod(p,t,rho,grav,qleft,rog,
     &       amoy,tcond,tfall,tcoag,tcoal)
        tgrow=1.d0/(1.d0/tcond+1.d0/tcoag+1.d0/tcoal)
c        write(*,*)qleft,rqleft,tcond,tgrow
        rqleft=qleft
        qleft=qleft*tgrow/tmix
       enddo
      else
c Si le melange est rapide, la fraction de condenses n'est pas affectee
c mais on recalcule la taille moyenne <a> d'un equilibre entre le
c melange et la condensation
       call secante(fgrow_mix,minx,maxx,amoy,1d-4,ier)
       amoy=10.d0**amoy
c       write(*,*)'Grow/Mix: amoy,ier',amoy,ier
      endif

      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function f(x)
      real*8 f,x
      f=x*x*x+1
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ffall_mix(x)
      implicit none
      real*8 ffall_mix,x

      real*8 size,tcond,tfall,tcoag,tcoal
      real*8 rp,rt,rro,rgrav,rqleft,rrog,rtmix
      common /rossow/rp,rt,rro,rgrav,rqleft,rrog,rtmix

      size=10.d0**x
      call rossow78_mod(rp,rt,rro,rgrav,rqleft,rrog,
     &     size,tcond,tfall,tcoag,tcoal)

      ffall_mix=log(tfall/rtmix)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function fgrow_mix(x)
      implicit none
      real*8 fgrow_mix,x

      real*8 size,tcond,tfall,tcoag,tcoal
      real*8 rp,rt,rro,rgrav,rqleft,rrog,rtmix
      common /rossow/rp,rt,rro,rgrav,rqleft,rrog,rtmix

      size=10.d0**x
      call rossow78_mod(rp,rt,rro,rgrav,rqleft,rrog,
     &     size,tcond,tfall,tcoag,tcoal)

      fgrow_mix=log((tcond)/rtmix)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine secante(f,a,b,x,dx,ier)
      implicit none
      real*8 f,a,b,x,dx
      integer ier
      real*8 f1,f2,fx,x1,x2
      real*8 rp,rt,rro,rgrav,rqleft,rrog,rtmix
      common /rossow/rp,rt,rro,rgrav,rqleft,rrog,rtmix

      x1=a
      x2=b
      f1=f(x1)
      f2=f(x2)
      if (f1*f2.gt.0.) then
       if (abs(f1).gt.abs(f2)) then
        x=x2
       else
        x=x1
       endif
       ier=1
       return
      endif
c On alterne dichotomie et secantes
      do while (abs(x2-x1).gt.dx)
       x=x1+(x2-x1)*f1/(f1-f2)
       fx=f(x)
       if (fx*f1.gt.0.d0) then
        x1=x
        f1=fx
       else
        x2=x
        f2=fx
       endif
       x=(x1+x2)/2
       fx=f(x)
       if (fx*f1.gt.0.d0) then
        x1=x
        f1=fx
       else
        x2=x
        f2=fx
       endif
      enddo

      ier=0
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rossow78_mod(p,t,ro,grav,xcond,rop,
     &     size,tcond,tfall,tcoag,tcoal)
c     Version: 16/04/00
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Calcul des temps caracteristiques de formation de nuages de      c
c     condensation (d'apres Rossow (1978), Icarus 30, 1-50)            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Auteur: T. Guillot
c Entrees:
c     c_mol: type de molecule (ex: 'H2O')
c     p,t,ro: pression, temperature, densite du milieu (cgs)
c     grav: acceleration gravifique (en cm/s2)
c     xcond: concentration totale des condenses
c Sorties:
c     lat: chaleur latente de transition
c     nb: dimension des tableaux suivants
c     size(nb): taille des "gouttes"
c     t...(nb): temps de condensation, sedimentation, coagulation,coalescence
c NOTE: les temps de coagulation et coalescence ne sont pas calcules

      implicit none

      character*8 c_mol
      real*8 p,t,ro,grav,xcond,lat,size,tcond,tfall,
     &     tcoag,tcoal

      integer i
      real*8 nu,f,smax,alpha,sig,mu,rn,lambda,a70,ros,rop,
     &     logamin,logamax,kn,re,nn
      logical phase

      real*8 pi,kbol,amu
      data pi/3.1415926539d0/,kbol/1.380658d-16/,amu/1.6605402d-24/

c-----------------------------------------------------------------------
 2000 format(7(1x,1pd10.3))

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Parametres fixes dans le programme:                              c
c     -Viscosite dynamique [nu]					       c
c     -Correction de masse finie [f]				       c
c     -Supersaturation maximale [smax]				       c
c     -Efficacite du collage [alpha]				       c
c     -Section efficace du gaz [sig]				       c
c     -Masse molaire moyenne du gaz [mu]			       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nu=ro/sqrt(t)
      f=5.d0
      smax=1d-2
      alpha=1.d0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Densite des condenses: ros                                       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ros=xcond*ro

      sig=pi*(1.d-8)**2
      mu=ro*kbol*t/(p*amu)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Autres quantites:                                                c
c     -Densite atmospherique en nombre de particules/cm3 [rn]	       c
c     -Libre parcours moyen du gaz [lambda]			       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      rn=ro/mu/amu
      lambda=sqrt(3.d0*pi)/(8.d0*rn*sig)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Taille des particules pour lesquelles Re=70                      c
c     Corrigee, d'apres Carlson et al. 1988                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      a70=(9.d0*270d0*nu*nu/(4.d0*rop*ro*grav))**(1.d0/3.d0)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     -Nombre de Knudsen [kn]                                          c
c     -Nombre de Reynolds [re]                                         c
c     -Densite de particules condensees pour une masse de condenses    c
c      supposee egale a la masse de vapeur saturante [nn]              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       kn=lambda/size
       re=70.d0*size/a70
       nn=ros/(4.d0/3.d0*pi*rop*size**3)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Calcul des temps de condensation, sedimentation, coagulation     c
c     et coalescence.						       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if (kn.lt.alpha) then
        tcond=ro*size*size/(2.d0*f*nu*ros*smax/rop)
        if (re.gt.1.d0) tcond=tcond/(1.d0+2.d-1*sqrt(re))
       else
        tcond=size/(3.d0*alpha*f*ros*smax/(2.d0*rop))/
     &       sqrt(2.d0*kbol*t/(pi*amu*mu))
       endif

       if (kn.gt.1.d0) then
        tfall=27d0*pi*ro*(2.d0*kbol*t/(pi*amu*mu))**(3.d0/2.d0)
     &       /(16.d0*rop*grav*grav*size)
       elseif (re.gt.70.d0) then
        tfall=sqrt(3.d0*ro*(kbol*t)**2/(40.d0*rop*(amu*mu)**2*
     &       grav**3*size))
       else
        tfall=9.d0*nu*kbol*t/(2.d0*rop*amu*mu*grav**2*size**2)
       endif

       if (kn.lt.1.d0) then
        tcoag=1.d0/((4.d0*kbol*t/(3.d0*nu))*nn)
       else
        tcoag=1.d0/(4.d0*sqrt(3.d0*size*kbol*t/rop)*nn)
       endif

       if (kn.gt.1.d0) then
        tcoal=1.d0/((4.d0*pi*rop*grav/(27.d0*ro))*
     &       sqrt(pi*amu*mu/(2.d0*kbol*t))*size**3*nn)
       elseif (re.gt.70.d0) then
        tcoal=1.d0/(sqrt(10.d0*pi*pi*rop*grav*size**5
     &        /(3.d0*ro))*nn)
       else
        tcoal=1.d0/((pi*rop*grav/(9.d0*nu))*size**4*nn)
       endif

c       write(*,'(1p,8d10.3)')p,t,size,tcond,tfall,tcoag,tcoal
      end
