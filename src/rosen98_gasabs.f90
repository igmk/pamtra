!=======================================================
!  ROSENKRANZ (1998) model  -- reference "Water vapor microwave
!  continuum absorption: a comparison of measurements and results"
!  To appear in Radio Science
!
! Top level interface (GasabsR98) by G. Petty
!
!   This subroutine calculates the mass extinction coefficient of the
!   dry and vapor components of air.
!
! Input:
!     F  = frequency (GHz),
!     Tk = absolute temperature (K)
!     Rhowv = water vapor density (kg/m**3).
!     Pa = Total air pressure (Pascals).
! Output:
!     Absair = extinction by dry air  (Np/km)
!     Abswv  = extinction by water vapor (Np/km)
!
!  2/8/2001 - fixed division by zero when rhowv = 0  (G. Petty)
!

subroutine gasabsR98(f,Tk,Rhowv,Pa,absair,abswv)

  use kinds

  implicit none

  real(kind=dbl), intent(in) :: f        ! frequency [GHz]
  real(kind=dbl), intent(in) :: Tk       ! temperature [K]
  real(kind=dbl), intent(in) :: pa       ! pressure [Pa]
  real(kind=dbl), intent(in) :: rhowv    ! water vapor density [kg/m**3]

  real(kind=dbl), intent(out) :: absair  ! extinction by dry air [Np/km]
  real(kind=dbl), intent(out) :: abswv   ! extinction by water vapor [Np/km]

  real(kind=dbl) :: pmb         ! pressure [mb]
  real(kind=dbl) :: vapden      ! water vpor density [g/m**3]
  real(kind=dbl) :: e           ! water vapor pressure [Pa]
  real(kind=dbl) :: q           ! specific humidity
  real(kind=dbl) :: Tv          ! virtuel temperature [K]
  real(kind=dbl) :: rhoair      ! moist air density [kg/m**3]

  real(kind=dbl) :: absn2       ! function to calculate extinction by n_2
  real(kind=dbl) :: o2abs       ! function to calculate extinction by o2
  real(kind=dbl) :: abh2o       ! function to calculate extinction by h_2o

  ! check for "reasonable" input values

  if (f .le. 0.0 .or. F .gt. 800.0) then
     print*, 'Frequency out of range in GasabsR98:', F
     stop
  endif
  if (Tk .lt. 100.0) then
     print *, 'Temperature out of range in GasabsR98: ', Tk
     stop
  endif
  if (pa .lt. 10.0 .or. pa .gt. 1.2e5) then
     print *, 'Pressure out of range in GasabsR98: ', pa
     stop
  endif
  ! convert pressure from Pa to Mb
  pmb = pa/100.0

  ! convert vapor density from kg/m**3 to g/m**3

  vapden = rhowv*1000.0

  ! get volume extinction coefficients
  absair = absn2(Tk,pmb,f) + o2abs(Tk,pmb,vapden,f)
  abswv = abh2o(Tk,pmb,vapden,f)
  !       abs_n2 = absn2(Tk,pmb,f)
  !       abs_o2 = o2abs(Tk,pmb,vapden,f)

  ! convert vapor density to vapor pressure
  e = rhowv*(Tk*461.5) 
  ! calculate specific humidity
  q = 0.622*e/pa
  ! calculate virtual temperature
  Tv = (1. + 0.61*q)*Tk
  ! moist air density
  rhoair = pa/(Tv*287.06)

  ! convert above from Np/km to m**2/kg

  !      absair = 0.001*absair/rhoair
  !      if (rhowv .eq. 0.0) then
  !         abswv = 0.0
  !      else
  !         abswv = 0.001*abswv/rhowv
  !      endif

  return

end subroutine gasabsR98

!******************************************************************
!This file contains subroutines for computing atmospheric absorption at
!microwave wavelengths, as supplied by P. Rosenkranz (2/98) via his
!anonymous ftp server (mesa.mit.edu; login anonymous; go to phil/lpl_rt)AN
!Consolidated into one file by G. Petty
!
!**********************************************************
!  Begin summaries
!**********************************************************
!     FUNCTION ABSN2(T,P,F)
!     ABSN2 = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR
!             (NEPER/KM)
!     T = TEMPERATURE (K)
!     P = PRESSURE (MB)
!     F = FREQUENCY (GHZ)
!
!************************************************************************
!     FUNCTION O2ABS(TEMP,PRES,VAPDEN,FREQ)
!
!     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
!              IN NEPERS/KM
!
!     NAME    UNITS    DESCRIPTION        VALID RANGE
!
!     TEMP    KELVIN   TEMPERATURE        (UNCERTAIN)
!     PRES   MILLIBARS PRESSURE           (3 TO 1000)
!     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
!                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
!     FREQ    GHZ      FREQUENCY          (0 TO 900)
!*************************************************************
!     FUNCTION ABH2O(T,P,RHO,F)
!
! PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
! 
!      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
!      T       KELVIN    I   TEMPERATURE
!      P       MILLIBAR  I   PRESSURE              .1 TO 1000
!      RHO     G/M**3    I   WATER VAPOR DENSITY
!      F       GHZ       I   FREQUENCY             0 TO 800
!      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT
!
!****************************************************************
!     FUNCTION ABLIQ(WATER,FREQ,TEMP)
!     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
!     FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
!     (INT. J. IR & MM WAVES V.12(17) JULY 1991
!     WATER IN G/M**3
!     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
!     TEMP IN KELVIN
!        PWR 8/3/92
!
!*************************************************************
!
!  Begin actual function definitions
!
!**************************************************************
function ABSN2(T,P,F)

  use kinds

  !     ABSN2 = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR
  !             (NEPER/KM)
  !     T = TEMPERATURE (K)
  !     P = PRESSURE (MB)
  !     F = FREQUENCY (GHZ)
  !
  implicit none

  real(kind=dbl) :: t,th,absn2,p,f

  TH = 300./T
  ABSN2 = 6.4E-14*P*P*F*F*TH**3.55

  return

end function absn2
!
!
function O2ABS(TEMP,PRES,VAPDEN,FREQ)
  !
  !     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
  !              IN NEPERS/KM
  !
  !      5/1/95  P. Rosenkranz
  !
  !     ARGUMENTS:
  !
  !     NAME    UNITS    DESCRIPTION        VALID RANGE
  !
  !     TEMP    KELVIN   TEMPERATURE        (UNCERTAIN)
  !     PRES   MILLIBARS PRESSURE           (3 TO 1000)
  !     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
  !                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
  !     FREQ    GHZ      FREQUENCY          (0 TO 900)
  !
  !     REFERENCE FOR EQUATIONS AND COEFFICIENTS:
  !     P.W. ROSENKRANZ, CHAP. 2 AND APPENDIX, IN ATMOSPHERIC REMOTE SENSING
  !      BY MICROWAVE RADIOMETRY (M.A. JANSSEN, ED. 1993)
  !     AND H.J. LIEBE ET AL, JQSRT V.48, PP.629-643 (1992)
  !     (EXCEPT: SUBMILLIMETER LINE INTENSITIES FROM HITRAN92)
  !
  use kinds

  implicit none

  integer :: k

  real(kind=dbl) :: TEMP,PRES,VAPDEN,FREQ
  real(kind=dbl) :: th, th1, b, preswv, presda, den, dfnr, df
  real(kind=dbl) :: str, sf1, sf2, sum, y

  real(kind=dbl) :: o2abs
  real(kind=dbl) :: X,WB300,W300(40),F(40),Y300(40),S300(40),&
       V(40),BE(40)
  !      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
  data F/118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,&
       59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,&
       56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,&
       55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,&
       53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,&
       52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7631,&
       487.2494, 715.3932, 773.8397, 834.1453/
  data S300/.2936E-14,.8079E-15, .2480E-14,.2228E-14,&
       .3351E-14,.3292E-14, .3721E-14,.3891E-14,&
       .3640E-14,.4005E-14, .3227E-14,.3715E-14,&
       .2627E-14,.3156E-14, .1982E-14,.2477E-14,&
       .1391E-14,.1808E-14, .9124E-15,.1230E-14,&
       .5603E-15,.7842E-15, .3228E-15,.4689E-15,&
       .1748E-15,.2632E-15, .8898E-16,.1389E-15,&
       .4264E-16,.6899E-16, .1924E-16,.3229E-16,&
       .8191E-17,.1423E-16, .6460E-15, .7047E-14, .3011E-14,&
       .1826E-14, .1152E-13, .3971E-14/
  data BE/.009,.015, .083,.084, 2*.212, 2*.391, 2*.626,&
       2*.915, 2*1.260, 1.660,1.665, 2.119,2.115, 2.624,2.625,&
       2*3.194, 2*3.814, 2*4.484, 2*5.224, 2*6.004, 2*6.844,&
       2*7.744, .048, .044, .049, .145, .141, .145/
  !      WIDTHS IN MHZ/MB
  data WB300/.56/, X/.8/
  data W300/1.63, 1.646, 1.468, 1.449, 1.382, 1.360,&
       1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,&
       1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 2*1.05,&
       2*1.02,2*1.00,2*.97,2*.94,2*.92,2*.89, 3*1.92, 3*1.81/
  data Y300/  -0.0233,  0.2408, -0.3486,  0.5227,&
       -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,&
       0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,&
       0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,&
       0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,&
       0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, 6*0./
  data V/  0.0079, -0.0978,  0.0844, -0.1273,&
       0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,&
       0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,&
       0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,&
       0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,&
       0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545, 6*0./
  !
  TH = 300./TEMP
  TH1 = TH-1.
  B = TH**X
  PRESWV = VAPDEN*TEMP/217.
  PRESDA = PRES -PRESWV
  DEN = .001*(PRESDA*B + 1.1*PRESWV*TH)
  DFNR = WB300*DEN
  SUM = 1.6E-17*FREQ*FREQ*DFNR/(TH*(FREQ*FREQ + DFNR*DFNR))
  do K=1,40
     DF = W300(K)*DEN
     Y = .001*PRES*B*(Y300(K)+V(K)*TH1)
     STR = S300(K)*exp(-BE(K)*TH1)
     SF1 = (DF + (FREQ-F(K))*Y)/((FREQ-F(K))**2 + DF*DF)
     SF2 = (DF - (FREQ+F(K))*Y)/((FREQ+F(K))**2 + DF*DF)
     SUM = SUM + STR*(SF1+SF2)*(FREQ/F(K))**2
  end do
  O2ABS = .5034E12*SUM*PRESDA*TH**3/3.14159

  return

end function o2abs

!*************************************************************
function ABH2O(T,P,RHO,F)
  !
  !  NAME- ABH2O    LANGUAGE- FORTRAN 77
  !
  ! PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
  ! 
  !  CALLING SEQUENCE PARAMETERS-
  !    SPECIFICATIONS
  !      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
  !      T       KELVIN    I   TEMPERATURE
  !      P       MILLIBAR  I   PRESSURE              .1 TO 1000
  !      RHO     G/M**3    I   WATER VAPOR DENSITY
  !      F       GHZ       I   FREQUENCY             0 TO 800
  !      ABH2O   NEPERS/KM O   ABSORPTION COEFFICIENT
  !
  !   REFERENCES-
  !   LINE INTENSITIES FROM HITRAN92 (SELECTION THRESHOLD=
  !     HALF OF CONTINUUM ABSORPTION AT 1000 MB).
  !   WIDTHS MEASURED AT 22,183,380 GHZ, OTHERS CALCULATED:
  !     H.J.LIEBE AND T.A.DILLON, J.CHEM.PHYS. V.50, PP.727-732 (1969) &
  !     H.J.LIEBE ET AL., JQSRT V.9, PP. 31-47 (1969)  (22GHz);
  !     A.BAUER ET AL., JQSRT V.37, PP.531-539 (1987) & 
  !     ASA WORKSHOP (SEPT. 1989) (380GHz);
  !     AND A.BAUER ET AL., JQSRT V.41, PP.49-54 (1989) (OTHER LINES).
  !   AIR-BROADENED CONTINUUM BASED ON LIEBE & LAYTON, NTIA 
  !     REPORT 87-224 (1987); SELF-BROADENED CONTINUUM BASED ON 
  !     LIEBE ET AL, AGARD CONF. PROC. 542 (MAY 1993), 
  !     BUT READJUSTED FOR LINE SHAPE OF
  !     CLOUGH et al, ATMOS. RESEARCH V.23, PP.229-241 (1989).
  !
  !   REVISION HISTORY-
  !    DATE- OCT.6, 1988  P.W.ROSENKRANZ - EQS AS PUBL. IN 1993.
  !          OCT.4, 1995  PWR- USE CLOUGH'S DEFINITION OF LOCAL LINE
  !                   CONTRIBUTION,  HITRAN INTENSITIES, ADD 7 LINES.
  !          OCT. 24, 95  PWR -ADD 1 LINE.
  !          JULY 7, 97   PWR -SEPARATE COEFF. FOR SELF-BROADENING, 
  !                       REVISED CONTINUUM.
  !
  !   LOCAL VARIABLES:

  use kinds

  implicit none

  real(kind=dbl) :: T,P,RHO,F,ABH2O
  integer NLINES,I,J
  parameter (NLINES=15)
  real(kind=dbl) :: DF(2),S1(NLINES),B2(NLINES),W3(NLINES),FL(NLINES),X(NLINES),&
       WS(NLINES),XS(NLINES)
  real(kind=dbl) PVAP,PDA,DEN,TI,TI2,SUM,WIDTH,WSQ,S,BASE,RES,CON
  !     LINE FREQUENCIES:
  data FL/22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508,&
       443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,&
       620.7008, 752.0332, 916.1712/
  !     LINE INTENSITIES AT 300K:
  data S1/ .1310E-13, .2273E-11, .8036E-13, .2694E-11, .2438E-10,&
       .2179E-11, .4624E-12, .2562E-10, .8369E-12, .3263E-11, .6659E-12,&
       .1531E-08, .1707E-10, .1011E-08, .4227E-10/
  !     T COEFF. OF INTENSITIES:
  data B2/ 2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,&
       3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/
  !     AIR-BROADENED WIDTH PARAMETERS AT 300K:
  data W3/.00281, .00281, .0023, .00278, .00287, .0021, .00186,&
       .00263, .00215, .00236, .0026, .00321, .00244, .00306, .00267/
  !     T-EXPONENT OF AIR-BROADENING:
  data X/.69, .64, .67, .68, .54, .63, .60, .66, .66, .65, .69, .69,&
       .71, .68, .70/
  !     SELF-BROADENED WIDTH PARAMETERS AT 300K:
  data WS/.01349, .01491, .0108, .0135, .01541, .0090, .00788,&
       .01275, .00983, .01095, .01313, .01320, .01140, .01253, .01275/
  !     T-EXPONENT OF SELF-BROADENING:
  data XS/ .61, .85, .54, .74, .89, .52, .50, .67, .65, .64, .72,&
       1.0, .68, .84, .78/
  !
  if(RHO.le.0.) then
     ABH2O = 0.
     return
  endif
  PVAP = RHO*T/217.
  PDA = P -PVAP
  DEN = 3.335E16*RHO
  TI = 300./T
  TI2 = TI**2.5
  !
  !      CONTINUUM TERMS
  CON = (5.43E-10*PDA*TI**3 + 1.8E-8*PVAP*TI**7.5)*PVAP*F*F 
  !
  !      ADD RESONANCES
  SUM = 0.
  do I=1,NLINES
     WIDTH = W3(I)*PDA*TI**X(I) + WS(I)*PVAP*TI**XS(I)
     WSQ = WIDTH*WIDTH
     S = S1(I)*TI2*exp(B2(I)*(1.-TI))
     DF(1) = F - FL(I)
     DF(2) = F + FL(I)
     !  USE CLOUGH'S DEFINITION OF LOCAL LINE CONTRIBUTION
     BASE = WIDTH/(562500. + WSQ)
     !  DO FOR POSITIVE AND NEGATIVE RESONANCES
     RES = 0.
     do J=1,2
        if(abs(DF(J)).lt.750.) RES = RES + WIDTH/(DF(J)**2+WSQ) - BASE
     end do
     SUM = SUM + S*RES*(F/FL(I))**2
  end do
  ABH2O = .3183E-4*DEN*SUM + CON

  return

end function ABH2O
!*************************************************************
function ABLIQ(WATER,FREQ,TEMP)
  !     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
  !     FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
  !     (INT. J. IR & MM WAVES V.12(17) JULY 1991
  !     WATER IN G/M**3
  !     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
  !     TEMP IN KELVIN
  !        PWR 8/3/92
  !

  use kinds

  implicit none

  real(kind=dbl) :: abliq, water, freq, temp
  real(kind=dbl) :: theta1, eps0, eps1, eps2, fp, fs
  complex(kind=dbl) :: EPS,RE

  if(WATER.le.0.) then
     ABLIQ = 0.
     return
  endif
  THETA1 = 1.-300./TEMP
  EPS0 = 77.66 - 103.3*THETA1
  EPS1 = .0671*EPS0
  EPS2 = 3.52 + 7.52*THETA1
  FP = (316.*THETA1 + 146.4)*THETA1 +20.20
  FS = 39.8*FP
  EPS = (EPS0-EPS1)/cmplx(1.,FREQ/FP) +&
       (EPS1-EPS2)/cmplx(1.,FREQ/FS) +EPS2
  RE = (EPS-1.)/(EPS+2.)
  ABLIQ = -.06286*aimag(RE)*FREQ*WATER

  return

end function ABLIQ



function VAPOR(T)
  !
  !  COMPUTES SATURATION H2O VAPOR PRESSURE (OVER LIQUID)
  !  USING LIEBE'S APPROXIMATION (CORRECTED)
  !
  !  T IN KELVIN
  !  VAPOR IN MBAR
  !                   PWR 4/8/92
  !
  use kinds

  implicit none

  real(kind=dbl) :: T, TH, VAPOR

  TH = 300./T
  VAPOR = 35.3*exp(22.64*(1.-TH))*TH**5
  return
end function VAPOR
