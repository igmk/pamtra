      SUBROUTINE CAL_REFRACTIVE_INDEX(PHASE, TEMPR, FRQ, PARTICLE_SIZE,&
          AS_RATIO, PARTICLE_MASS, AXI, REF)

      IMPLICIT NONE
      REAL*8 TEMPR, FRQ, FREQ, PARTICLE_SIZE, AS_RATIO, AXI, MASS
      REAL*8 PARTICLE_MASS
      COMPLEX*16 REF
      REAL*8 SAL
! INPUT:
!  TEMPR: TEMPERAUTRE, IN K
!  FRQ: FREQUENCY, IN Hz,
!  FREQ: in REFRECTIVE INDEX CALCULATION, CONVERT FREQ TO GHz
!  PARTICLE_SIZE: MAXIMUM DIMENSION OF A SINGLE PARTICLE, IN METER
!  AS_RATIO: ASPECT RATIO, MINIMUN_DIMENSION/MAXIMUM_DIMENSION
!  MASS: MASS OF A SINGLE PARTICLE, IN gram
! OUTPUT:
!  AXI: EQUIVALENT VOLUME SPHERE RADIUS, IN METER
!  REF: REFRACTIVE INDEX OF SNOW/ICE
  
      CHARACTER*1 PHASE
      IF(PHASE.EQ.'S'.OR.PHASE.EQ.'s')THEN
	CALL REFSNOW(TEMPR, FRQ, PARTICLE_SIZE,&
          AS_RATIO, PARTICLE_MASS, AXI, REF)
      ELSEIF(PHASE.EQ.'L'.OR.PHASE.EQ.'l')THEN
        SAL = 0
        CALL REFWATER(SAL,FRQ,TEMPR,REF)
        AXI = PARTICLE_SIZE/2.
!	WRITE(*,*)'REFRACTIVE_INDEX.f','AXI = ',AXI
      ENDIF

      RETURN
      END

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     SNOW REFRACTIVE INDEX

      SUBROUTINE REFSNOW(TEMPR, FRQ, PARTICLE_SIZE,&
          AS_RATIO, PARTICLE_MASS, AXI, SNOW_REF)
      
      IMPLICIT NONE

      REAL*8 TEMPR, FRQ,FREQ, PARTICLE_SIZE, AS_RATIO, AXI, MASS
      REAL*8 PARTICLE_MASS
! INPUT:
!  TEMPR: TEMPERAUTRE, IN K
!  FRQ: FREQUENCY, IN Hz,
!  FREQ: in REFRECTIVE INDEX CALCULATION, CONVERT FREQ TO GHz
!  PARTICLE_SIZE: MAXIMUM DIMENSION OF A SINGLE PARTICLE, IN METER
!  AS_RATIO: ASPECT RATIO, MINIMUN_DIMENSION/MAXIMUM_DIMENSION
!  MASS: MASS OF A SINGLE PARTICLE, IN gram
! OUTPUT:
!  AXI: EQUIVALENT VOLUME SPHERE RADIUS, IN METER
!  SNOW_REF: REFRACTIVE INDEX OF SNOW
  
! INTERNAL VARIABLES
      REAL*8 MIN_DIM, MAX_DIM, VOLM
      REAL*8 AIR_DENS, ICE_DENS
      REAL*8 REALP, ALPHA1, B1, B, B2, BETAM, DELTA_BETA, BETA1, IMAGP
      COMPLEX*16 ICE_REF, SNOW_REF
      CHARACTER*1 PHASE
! MIN_DIM: MINIMUM DIMENSION, IN M
! MAX_DIM: MAXIMUM DIMENSION, IN M
! VOLM: PARTICLE VOLUME, IN CM3
! REALP,...IMAGP, TEMP VALUE FOR SNOW PERMITTIVITY CALCULATION
! AIR_DENS: THE DENSITY OF AIR, g/cm3 (MIXING ICE WITH AIR -- SNOW)
! ICE_DENS: THE DENSITY OF ICE, g/cm3
!      WRITE(*,*)PARTICLE_SIZE, AS_RATIO
      AIR_DENS = 1.25E-3
      ICE_DENS = 0.917
      FREQ = FRQ/1E9
      MASS = PARTICLE_MASS 

!      WRITE(*,*)MASS
! THE MINIMUM DIMENSION OF OBLATES, cm. PARTICLE_SIZE in meter
      MIN_DIM = PARTICLE_SIZE*AS_RATIO*1.0E2
! THE MAXIMUM DIMENSION OF OBLATES, cm. PARTICLE_SIZE in meter
      MAX_DIM = PARTICLE_SIZE*1.0E2
! PARTICLE VOLUME,
      VOLM = 4.0E0*3.1415926/3.0E0*(MAX_DIM/2.0)*(MAX_DIM/2.0)&
          *(MIN_DIM/2.0)
! EQUIVALENT-VOLUME SPHERE'S RADIUS, in meter
      AXI = 0.5*(PARTICLE_SIZE**3.0*AS_RATIO)**(1.0E0/3.0E0)
!====================================================================
!      WRITE(*,*)MAX_DIM, MIN_DIM,AS_RATIO

!=========== calculate refractive index of ice
!=== real part of ice dielectric
      REALP = 3.1884E0+9.1*1E-4*(TEMPR-273.0E0)
!=== imaginary part of ice dielectric%%%%%%%%
      alpha1 = (0.00504+0.0062*(300D0/TEMPR-1E0))&
     	*exp(-22.1*(300E0/TEMPR-1.0E0))
      B1 = 0.0207E0
      b = 335.0E0
      B2 = 1.16*1.0E-11
      betam = B1/TEMPR*exp(b/TEMPR)/(exp(b/TEMPR)-1E0)**2D0+B2*FREQ*FREQ
      delta_beta = exp(-9.963E0+0.0372E0*(TEMPR-273.16E0))
      beta1 = betam+delta_beta
      IMAGP = alpha1/FREQ+beta1*FREQ

      ICE_REF = REALP*(1.0,0)+(alpha1/FREQ+beta1*FREQ)*(0,1E0)
      ICE_REF = SQRT(ICE_REF)
!======================================================================
!        RELATIVE VOLUME OF ICE IN PARTICLE
      SNOW_REF = (2*(AIR_DENS-MASS/VOLM)/(AIR_DENS-ICE_DENS)&
     	     *ICE_REF*ICE_REF-2*(AIR_DENS-MASS/VOLM)/(AIR_DENS-ICE_DENS)&
            +ICE_REF*ICE_REF+2)/&
            (ICE_REF*ICE_REF+2-(AIR_DENS-MASS/VOLM)/(AIR_DENS-ICE_DENS)&
            *ICE_REF*ICE_REF+(AIR_DENS-MASS/VOLM)/(AIR_DENS-ICE_DENS))
      SNOW_REF = SQRT(SNOW_REF)
!      WRITE(*,*)AIR_DENS, MASS, VOLM, ICE_DENS, ICE_REF
!      WRITE(*,*)'ICE_REF',ICE_REF
!      WRITE(*,*)SNOW_REF
!**************************************************************************

      RETURN
      END

!     WATER REFRACTIVE INDEX
      SUBROUTINE REFWATER(SAL1,FREQ,TEMPR,REF)
      IMPLICIT NONE
      REAL*8 SAL1, FREQ, Temp, MRR, MRI, FR, SAL,TEMPR
      COMPLEX*16 first_term,second_term,third_term,EPS,NN,TERM1,TERM2
      REAL*8 a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8
      REAL*8 a_9,a_10,a_11,a_12,a_13,a_14,a_15,a_16,a_17,a_18
      REAL*8 EPS_S, EPS_1, tau_1, tau_2, EPS_INF
      REAL*8 c_alpha_0, d_alpha_0, alpha_0, d_P, P, sigma_35
      REAL*8 SIGMA,PI
      COMPLEX*16 REF
      REAL*8 Q,C_P,ALPHA_1
            
      SAL = SAL1 * 1E-3
!      write(*,*)freq,sal
!*** Check the input ranges:
! IF FREQ GT 1d12 AND verbose EQ 1 THEN print, '!!!Frequency range: 0-1000 GHz!!!, EXTRAPOLATION'
! IF T GT 303.15 OR T LT 273.15 AND verbose EQ 1 THEN print, '!!!Temperature range: 0-30 degC!!!, EXTRAPOLATION'
! IF SAL GT 40. OR SAL LT 0. AND verbose EQ 1 THEN print, '!!!Salinity range: 0-40 ppt!!!, EXTRAPOLATION'

!*** Convert Temperature from Kelvin to degree Celsius
      Temp = TEMPR - 273.15 
!;--------------------------------------------------------------------------------------------------------
!;COEFFS AND CALCULATION OF eps(FREQ, Temp, SAL) according to (5.21, p.445)
!;--------------------------------------------------------------------------------------------------------

!;*** Coefficients a_i (Table 5.5 or p. 454):

      a_1  =  0.46606917e-2
      a_2  = -0.26087876e-4
      a_3  = -0.63926782e-5
      a_4  =  0.63000075e1
      a_5  =  0.26242021e-2
      a_6  = -0.42984155e-2
      a_7  =  0.34414691e-4
      a_8  =  0.17667420e-3
      a_9  = -0.20491560e-6
      a_10 =  0.58366888e3
      a_11 =  0.12634992e3
      a_12 =  0.69227972e-4
      a_13 =  0.38957681e-6
      a_14 =  0.30742330e3
      a_15 =  0.12634992e3
      a_16 =  0.37245044e1
      a_17 =  0.92609781e-2
      a_18 = -0.26093754e-1


!;*** Calculate parameter functions (5.24)-(5.28), p.447

      EPS_S   = 87.85306 * EXP(-0.00456992 * Temp - a_1*SAL - &
      		a_2*SAL*SAL - a_3*SAL*Temp)
      EPS_1   = a_4 * EXP( -a_5*Temp - a_6*SAL - a_7*SAL*Temp)
      tau_1   = (a_8 + a_9*SAL) * EXP( a_10 / (Temp + a_11)) * 1e-9
      tau_2   = (a_12 + a_13*SAL) * EXP( a_14 / (Temp + a_15)) * 1e-9
      EPS_INF = a_16 + a_17*Temp + a_18*SAL
!      write(*,*)eps_s,eps_1,tau_1,tau_2,eps_inf

!;*** Calculate seawater conductivity (5.20), p.437

      IF (SAL.GT.0.) THEN 
	c_alpha_0 =  (6.9431 + 3.2841 * SAL - 0.099486 * SAL**2.)
	d_alpha_0 =  (84.85 + 69.024 * SAL + SAL**2.)
	alpha_0   =  c_alpha_0 / d_alpha_0 
	alpha_1   = 49.843 - 0.2276 * SAL + 0.00198 * SAL**2.
	Q = 1.000 + alpha_0*(Temp - 15.0) / (Temp + alpha_1)
	c_P = (37.5109 + 5.45216 * SAL + 0.014409 * SAL**2.)
	d_P = (1004.75 + 182.283 * SAL + SAL**2.)
	P = SAL * c_P / d_P
	sigma_35  = 2.903602 + 8.607d-2 * Temp+4.738817d-4*Temp**2. &
      		-2.991d-6 * Temp**3. + 4.3041d-9 * Temp**4.
        SIGMA = sigma_35 * P * Q
      ELSE
	SIGMA = 0
      ENDIF
!;just reduce PC time.... calculation would give the same!

      PI = ACOS(-1.0)
!;*** Finally apply the interpolation formula (5.21)
      TERM1 = 1.*(1.0D0,0) -2.*PI*FREQ*tau_1*(0,1.0D0)
      first_term  = (EPS_S - EPS_1) /  TERM1
      TERM2 = 1.*(1.0D0,0) -2.*PI*FREQ*tau_2*(0,1.0D0)
      second_term = (EPS_1 - EPS_INF) / TERM2
!;third_term  = DCOMPLEX(EPS_INF, (17.9751d * SIGMA / FREQ ))
      third_term = EPS_INF

      EPS = first_term + second_term + third_term
!      write(*,*)  term1,term2
!      write(*,*)  first_term,second_term,third_term
!;calculate refractivity and mass/volume absorption coefficient
!;frequency in Hz, lwc in g/m³
!      RE = (EPS-1)/(EPS+2)
!      MASS_ABSCOF = 6.d*3.1415926*IMAG(RE)*FREQ*1d-3/cl
!      VOL_ABSCOF = MASS_ABSCOF * LWC
!;*** Convert to refractive index

      NN  = SQRT(EPS)
      MRR = REAL(NN)
      MRI = IMAG(NN)
      REF = NN

      IF(MRI.LT.0)REF=CONJG(NN)      

      RETURN
      END

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5




      SUBROUTINE VECTOR_SUM(DIM1, VECTOR, SUM1)
            
      INTEGER DIM1
      REAL*8 SUM1, VECTOR(1:DIM1)
      INTEGER I
      SUM1 = 0    ! initialize the value
      DO I = 1, DIM1, 1
         SUM1 = SUM1 + VECTOR(I)
      ENDDO

      RETURN
      END

      SUBROUTINE GAUSS_DISTRIBUTION(BETA, CANTING_NUM, CANTING_STD, &
     	CANTING_MEAN, CANTING_WEIGHTS)
      REAL*8 CANTING_STD, CANTING_MEAN, CANTING_WEIGHTS
      REAL*8 SUM_P, CANTING_ANGLE, BETA
      INTEGER CANTING_NUM
      
      PI = ACOS(-1.0)
      SUM_P = 0.0
      
      DO 1300 I = 1, CANTING_NUM
      CANTING_ANGLE = 90.0/(REAL(CANTING_NUM)-1)*(REAL(I)-1)
      SUM_P = SUM_P + &
           EXP(-(CANTING_ANGLE-CANTING_MEAN)**2.0/(2.0*CANTING_STD**2))
1300  CONTINUE      

      IF (CANTING_STD.EQ.0) THEN
	  IF (BETA.EQ.0) THEN
              CANTING_WEIGHTS = 1
	  ELSE
	      CANTING_WEIGHTS = 0
	  ENDIF
      ELSE 
          CANTING_WEIGHTS = 1.0/SUM_P*&
           EXP(-((BETA-CANTING_MEAN)**2.0)/(2.0*CANTING_STD**2.0))
!      WRITE(*,*)CANTING_WEIGHTS
      ENDIF

      IF (CANTING_NUM.EQ.1) THEN
	CANTING_WEIGHTS = 1D0
      ENDIF



      RETURN
      END

!       SUBROUTINE GAUSS_LEGENDRE_QUADRATURE
!      .                          (NUM, ABSCISSAS, WEIGHTS)
! C        Generates the abscissas and weights for an even 2*NUM point
! C      Gauss-Legendre quadrature.  Only the NUM positive points are returned.
!       INTEGER  NUM
!       REAL*8   ABSCISSAS(1), WEIGHTS(1)
!       INTEGER  N, I, J, L
!       REAL*8   X, XP, PL, PL1, PL2, DPL, TINY
!       PARAMETER (TINY=3.0D-14)
! 
!       N = 2*NUM
!       DO J = 1, NUM
!         X = COS(3.141592654*(J-.25)/(N+.5))
!         I = 0
! 100     CONTINUE
!           PL1 = 1
!           PL = X
!           DO L = 2, N
!             PL2 = PL1
!             PL1 = PL
!             PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
!           ENDDO
!           DPL = N*(X*PL-PL1)/(X*X-1)
!           XP = X
!           X = XP - PL/DPL
!           I = I+1
!         IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
!         ABSCISSAS(NUM+1-J) = X
!         WEIGHTS(NUM+1-J) = 2/((1-X*X)*DPL*DPL)
!       ENDDO
! 
!       RETURN
!       END
! 
!       SUBROUTINE LOBATTO_QUADRATURE
!      .                          (NUM, ABSCISSAS, WEIGHTS)
! C        Generates the abscissas and weights for an even 2*NUM point
! C      Gauss-Legendre quadrature.  Only the NUM positive points are returned.
!       INTEGER  NUM
!       REAL*8   ABSCISSAS(*), WEIGHTS(*)
!       INTEGER  N, N1, I, J, L
!       REAL*8   X, XP, PL, PL1, PL2, DPL, D2PL, CI, TINY
!       PARAMETER (TINY=3.0D-14)
!       
! c	OPEN(UNIT=4,FILE='LOBATTO.DAT')
!       N = 2*NUM
!       N1 = N-1
!       CI = 0.50
!       IF (MOD(N,2) .EQ. 1) CI = 1.00
!       DO J = 1, NUM-1
!         X = SIN(3.141592654*(J-CI)/(N-.5))
!         I = 0
! 100     CONTINUE
!           PL1 = 1
!           PL = X
!           DO L = 2, N1
!             PL2 = PL1
!             PL1 = PL
!             PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
!           ENDDO
!           DPL = N1*(X*PL-PL1)/(X*X-1)
!           D2PL = (2.D0*X*DPL-N1*(N1+1)*PL) / (1D0-X*X)
!           XP = X
!           X = XP - DPL/D2PL
!           I = I+1
!         IF (ABS(X-XP).GT.TINY .AND. I.LT.10) GO TO 100
!         ABSCISSAS(J) = X
!         WEIGHTS(J) = 2.0D0/(N*N1*PL*PL)
! C        write(*,*) j,abscissas(j),ACOS(ABSCISSAS(J))*180./3.1415926,
! C     &      ACOS(-ABSCISSAS(J))*180./3.1415926,
! C     &	  weights(j)
!       ENDDO
!       ABSCISSAS(NUM) = 1.D0
!       WEIGHTS(NUM) = 2.D0/(N*N1)
! C      write(*,*)NUM,ABSCISSAS(NUM),
! C     &	ACOS(ABSCISSAS(NUM))*180./3.1415926,
! C     &	weights(j)
! 101   FORMAT(I3,4F15.8)
!       RETURN
!       END
