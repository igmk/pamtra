!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       SUBROUTINE EPSI(EPSILO,F,T,CONSTI,SNOWDE)

!***** Computation of dielectric constant:
!***** Original programm (now called epsilon.bck) written by C. Prigent
!***** Modifications added Nov 2002: different del. const for ice, added 
!***** porous ice,
!***** rain water, as well as comments
!***** by Martina Wiedner (mwiedner@cfa.harvard.edu)

!***** epsilon: Complex number with EM1-iEM2
!***** F - frequency in GHz
!***** T - temperature in K
!***** CONSTI - number code for constituant: 1 - pure water, 2 - sea water,
!*****        3 - earth  4 - dry snow, 5 - ice, 6 - porous ice, 7 - rain water
!***** SEASAL - salt content in sea water
!***** LANDMV - soil wetness
!***** SNOWDE - density of snow or porous ice in [g/cm^3]

       IMPLICIT NONE

       REAL*8 t,f,C,TH,E0,ALPHA,&
       EM1,EM2,PI,SNOWDE,&
       BETA,theta,v,airde,icede
       INTEGER CONSTI                            
       COMPLEX*16 EPSILO,EA,EB,IMA,MAT,&
       epssnow,epssnowwet,epsice,epswater,epsair
   
       DATA pi/3.1415927/C/2.99792458E10/ 

       IMA=(0.,1.)
                      
!      * Pure water (Manabe,Liebe) = 1
!      * Sea (Klein,Swift) = 2
!      * Earth (solid+pure water) = 3
!      * Dry snow (ice+air, Ulaby et al.) = 4
!      * Ice (Warren) = 5
!      * Porous ice (treated like 4) = 6
!      * Rain (treated like 1) = 7

!      * For pure water, earth or rain using formulae from Manabe, Liebe (1987)
       IF ((CONSTI.EQ.31).OR.(CONSTI.EQ.35).OR.&
     (CONSTI.EQ.41).OR.(CONSTI.EQ.45)) THEN
         TH=300./T-1.
         E0=77.66+103.3*TH
         EA=(E0-5.48)*F/(F-IMA*(20.09-142.*TH+294.*TH*TH))
         EB=(5.48-3.51)*F/(F-IMA*(590.-1500.*TH))
         MAT=E0-EA-EB
         EM2=-dIMAG(MAT)
         EM1=DREAL(MAT)
       END IF

       IF ((CONSTI.EQ.33).OR.(CONSTI.EQ.32).OR.(CONSTI.EQ.34)&
      .OR.(CONSTI.EQ.43).OR.(CONSTI.EQ.42).OR.(CONSTI.EQ.44)&
      .OR.(CONSTI.EQ.36)) THEN
!         * first calculation of the ice permittivity

!         * ice also from Warren very old code
!         TS=T-273.15
!         LAMBDA=C*1.E2/(F*1.E9)
!         EMINF=3.168
!         ALPHA=.288+.0052*TS+.00023*TS*TS
!         SIG=1.26*EXP(-12500./(TS+273.)/1.9869)
!         EMS=203.168+2.5*TS+.15*TS*TS
!         LAMBS=9.990288E-4*EXP(13200./(TS+273.)/1.9869)
!         EM1=(EMS-EMINF)*(1+((LAMBS/LAMBDA)**(1-ALPHA))*
!     S   SIN(ALPHA*PI/2.))
!         EM1=EM1/(1.+2*((LAMBS/LAMBDA)**(1.-ALPHA))*
!     S   SIN(ALPHA*PI/2.)+((LAMBS/LAMBDA)**(2*(1.-ALPHA))))+EMINF
!         EM2=(EMS-EMINF)*((LAMBS/LAMBDA)**(1.-ALPHA))*COS(PI*ALPHA/2.)
!          EM2=EM2/(1.+2*((LAMBS/LAMBDA)**(1.-ALPHA))*SIN(ALPHA*PI/2)+
!     S    ((LAMBS/LAMBDA)**(2*(1.-ALPHA))))+SIG*LAMBDA/18.8496E10
!         EM2=5.E-3

!         * ice from Warren (1984) by Wiedner et al., 2004
!         According to Warren, EM1 is known quite accurately
!         Epsilon = Epsilon (re) - i Epsilon(im)
!         Epsilon (re) = n(re)^2 - n(im)^2.
!         For small n(im): Epsilon(re) = n(re)^2
!         There is a factor of 10 uncertainty in n(im), due to the lack
!         of accurate measurements. m(im) is about 0.002. On top there
!         seems to be a temperature and frequency dependency (about a factor
!         of 2 between 12 and 60GHz and another factor of 2 between -1 and
!         -60C). Given the large uncertainty in the absolute measure and the
!         anyhow very small index of absoprtion, we did not think it necessary
!         to add the temp and frequency dependency. However, test should be
!         performed to check the effect of a 10 times bigger and a 10times
!         smaller EM2 on the predicted brightness temperatures.
!         Epsilon (im) = 2*n(re)*n(im)
!          EMINF=3.093+0.72E-4*T+0.11E-5*T**2
!          EM1=EMINF
!          EM2=0.000001
!          IF (EM1.GT.0) THEN
!          EM2=2.0*0.002*SQRT(EM1)
!          END IF

!         * ice permittivity from Matzler, 2005
!         * we check that it is not that different from the others
          if (T.gt.243) then
             EM1=3.1884+9.1e-4*(T-273.)
          else
             EM1=3.1
          end if
          theta=(300./T)-1.
          alpha=(0.00504+0.0062*theta)*exp(-22.1*theta)
          beta=(0.0207*exp(335./T))/(T*(exp(335./T)-1.)**2)+1.16e-11*f*f
          beta=beta+exp(-9.963+0.0372*(T-273.16))
          EM2=alpha/F+beta*F+1.16e-11*F*F*F
          epsice=EM1-IMA*EM2

!        *  snow and graupel permittivity
!        Dry snow and porous ice after Ulbay, Moore & Fung, 1986, p.2060
!        (we decided to treat porous ice in the same way as dry snow)
         IF ((CONSTI.EQ.32).OR.(CONSTI.EQ.34).OR.(CONSTI.EQ.36)&
        .OR.(CONSTI.EQ.42).OR.(CONSTI.EQ.44)) THEN      
!           if (CONSTI.eq.4) SNOWDE = 0.100
!           if (CONSTI.eq.6) SNOWDE = 0.400
!            if (consti.eq.4) then
!               if (snowde.le.0.1) snowde=0.1
!            end if
!          * Used in Wiedner et al, 2004
!          * Ulaby, Moore & Fung, 1986, p.2060
!          (we decided to treat porous ice in the same way as dry snow)
!          EM1ICE=EM1
!          EM2ICE=EM2
!          EM1=(1.+.51*SNOWDE)**3        ! SNOWDE en g/cm^3
!          EM2=0.34*(SNOWDE/0.916)*EM2ICE/((1-0.417*(SNOWDE/0.916))**2)
!c          * previously used:
!c          Equation E.84 p.2063
!          EM2=3.*(SNOWDE/0.916)*EM2ICE*(EM1**2)*(2.*EM1+1.)/
!     S    (EM1ICE+2*EM1)/(EM1ICE+2*(EM1**2))
!c          According to figure E.28 p. 2065, eq. E88 provides a slightly
!c          better fit to the data, as it falls between the curve for
!c          -18 C and that for 0C

!           if (consti.eq.6) then
!c          * Maxwell Garnett formula (Skofronick et al., IEEE, 2002)
!c          * air inclusion; ice host
           icede=0.916
           airde=1.e-3
           epsair=(1.,0.)
           v=1.-(snowde-airde)/(icede-airde)
           epssnow=(epsair-epsice)/(epsair+2.*epsice)
           epssnow=1.+(3*v*epssnow)/(1.-v*epssnow)
           epssnow=epssnow*epsice
           EM1=dble(epssnow)
           EM2=-dimag(epssnow)
!          end if


!ccccccccccxxxxxxxxxxxxxx
!c          * Maxwell Garnett formula (Skofronick et al., IEEE, 2002)
!c           * air host; ice inclusion
!           icede=0.916
!           airde=1.e-3
!           epsair=(1.,0.)
!           v=(snowde-airde)/(icede-airde)
!           epssnow=(epsice-epsair)/(epsice+2.*epsair)
!           epssnow=1.+(3*v*epssnow)/(1.-v*epssnow)
!           epssnow=epssnow*epsair
!           EM1=dble(epssnow)
!           EM2=-dimag(epssnow)
 
!c          * if wet snow, Maxwell Garnett formula again (Skofronick et al)
          if (CONSTI.EQ.32.OR.CONSTI.EQ.42) then
!c            * first calculation of water (see above)
             TH=300./T-1.
             E0=77.66+103.3*TH
             EA=(E0-5.48)*F/(F-IMA*(20.09-142.*TH+294.*TH*TH))
             EB=(5.48-3.51)*F/(F-IMA*(590.-1500.*TH))
             epswater=E0-EA-EB
!            * wetness as a function of T
             if (T.le.258.15) v=0.        
             if (T.gt.258.15.and.T.lt.273.15) v=(T-258.15)/100.        
             if (T.ge.273.15) v=.15
             epssnowwet=(epswater-epssnow)/(epswater+2.*epssnow)
             epssnowwet=1.+(3*v*epssnowwet)/(1.-v*epssnowwet)
             epssnowwet=epssnowwet*epssnow
             EM1=dble(epssnowwet)
             EM2=-dimag(epssnowwet)
           end if

         END IF
       END IF   

       EPSILO=EM1-IMA*EM2

       return
       END
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


