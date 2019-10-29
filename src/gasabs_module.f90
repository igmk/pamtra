!+ functions to calculate absorption by various gases
! 
module gasabs_module 
    !
    ! Description:
    !   This module contains functions, that calculate the absortpion for
    !   different gases.
    !
    ! Current Code Owner: IGMK
    !
    ! History:
    !
    ! Version   Date       Comment
    ! -------   ----       -------
    ! 0.1       21/09/2009 Code adaption from G. Petty - M. Mech
    ! 0.2       27/02/2013 Application of European Standards for Writing and
    !                      Documenting Exchangeable Fortran 90 Code - M. Mech
    ! Code Description:
    !   Language:		Fortran 90.
    !   Software Standards: "European Standards for Writing and  
    !     Documenting Exchangeable Fortran 90 Code". 
    !
    ! Define procedures contained in this module:
    !
    !******************************************************************
    !This file contains subroutines for computing atmospheric absorption at
    !microwave wavelengths, as supplied by P. Rosenkranz (2/98) via his
    !anonymous ftp server (mesa.mit.edu; login anonymous; go to phil/lpl_rt)AN
    !Consolidated into one file by G. Petty
    !
    !**********************************************************
    !  Begin summaries
    !**********************************************************
    !     real function absn2(tempK,pres,freq)
    !     absn2 = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR
    !             (NEPER/KM)
    !     tempK = TEMPERATURE (K)
    !     pres = PRESSURE (MB)
    !     freq = FREQUENCY (GHZ)
    !
    !************************************************************************
    !     real function o2abs(tempK,pres,vapden,freq)
    !
    !     PURPOSE: RETURNS ABSORPTION COEFFICIENT DUE TO OXYGEN IN AIR,
    !              IN NEPERS/KM
    !
    !     NAME    UNITS    DESCRIPTION        VALID RANGE
    !
    !     tempK   KELVIN   TEMPERATURE        (UNCERTAIN)
    !     PRES   MILLIBARS PRESSURE           (3 TO 1000)
    !     VAPDEN  G/M**3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION
    !                      DUE TO GREATER BROADENING EFFICIENCY OF H2O)
    !     FREQ    GHZ      FREQUENCY          (0 TO 900)
    !*************************************************************
    !     real function abh2o(tempK,pres,rho,freq)
    !
    ! PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
    !
    !      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
    !      tempK   KELVIN    I   TEMPERATURE
    !      pres    MILLIBAR  I   PRESSURE              .1 TO 1000
    !      rho     G/M**3    I   WATER VAPOR DENSITY
    !      freq    GHZ       I   FREQUENCY             0 TO 800
    !      abh2o   NEPERS/KM O   ABSORPTION COEFFICIENT
    !
    !****************************************************************
    !     real function abliq(water,freq,tempK)
    !     COMPUTES ABSORPTION IN NEPERS/KM BY SUSPENDED WATER DROPLETS
    !     FROM EQUATIONS OF LIEBE, HUFFORD AND MANABE
    !     (INT. J. IR & MM WAVES V.12(17) JULY 1991
    !     WATER IN G/M**3
    !     FREQ IN GHZ     (VALID FROM 0 TO 1000 GHZ)
    !     tempK IN KELVIN
    !        PWR 8/3/92
    !
    !*************************************************************
contains 
    !*************************************************************
    !
    !  Begin actual function definitions
    !
    !**************************************************************
    function absn2(tempK,pres,freq)

        use kinds, only: dbl

        !     absn2 = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR (NEPER/KM)
        !     tempK = TEMPERATURE (K)
        !     pres = PRESSURE (MB)
        !     freq = FREQUENCY (GHZ)
        !
        implicit none

        real(kind=dbl), intent(in) :: tempK,pres,freq
        real(kind=dbl) :: absn2
        real(kind=dbl) :: th

        th = 300._dbl/tempK
        absn2 = 6.4d-14*pres*pres*freq*freq*th**3.55_dbl

        return

    end function absn2
    !
    !
    function o2abs(tempK,pres,vapden,freq)
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
        !     tempK   KELVIN   TEMPERATURE        (UNCERTAIN)
        !     pres    MILLIBARS PRESSURE           (3 TO 1000)
        !     vapden  G/M^3   WATER VAPOR DENSITY  (ENTERS LINEWIDTH CALCULATION DUE TO GREATER BROADENING EFFICIENCY OF H2O)
        !     freq    GHZ      FREQUENCY          (0 TO 900)
        !
        !     REFERENCE FOR EQUATIONS AND COEFFICIENTS:
        !     P.W. ROSENKRANZ, CHAP. 2 AND APPENDIX, IN ATMOSPHERIC REMOTE SENSING
        !      BY MICROWAVE RADIOMETRY (M.A. JANSSEN, ED. 1993)
        !     AND H.J. LIEBE ET AL, JQSRT V.48, PP.629-643 (1992)
        !     (EXCEPT: SUBMILLIMETER LINE INTENSITIES FROM HITRAN92)
        !
        use kinds, only: dbl, long
        use constants, only: pi
        
        implicit none

        integer(kind=long) :: k

        real(kind=dbl), intent(in) :: tempK,pres,vapden,freq
        real(kind=dbl) :: th, th1, b, preswv, presda, den, dfnr, df
        real(kind=dbl) :: str, sf1, sf2, loc_sum, y

        real(kind=dbl) :: o2abs
        !      WIDTHS IN MHZ/MB
        real(kind=dbl) :: x = 0.8, &
                          wb300 = 0.56
        real(kind=dbl), dimension(40) :: w300 = (/&
          1.63, 1.646, 1.468, 1.449, 1.382, 1.360,&
          1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171,&
          1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 1.05, 1.05,&
          1.02, 1.02, 1.0, 1.0, 0.97, 0.97, 0.94, 0.94, 0.92, 0.92, &
          0.89, 0.89, 1.92, 1.92, 1.92, 1.81, 1.81, 1.81/)
        real(kind=dbl), dimension(40) :: y300 = (/&
          -0.0233,  0.2408, -0.3486,  0.5227,&
          -0.5430,  0.5877, -0.3970,  0.3237, -0.1348,  0.0311,&
          0.0725, -0.1663,  0.2832, -0.3629,  0.3970, -0.4599,&
          0.4695, -0.5199,  0.5187, -0.5597,  0.5903, -0.6246,&
          0.6656, -0.6942,  0.7086, -0.7325,  0.7348, -0.7546,&
          0.7702, -0.7864,  0.8083, -0.8210,  0.8439, -0.8529, &
          0., 0., 0., 0., 0., 0./)
        real(kind=dbl), dimension(40) :: v = (/&
          0.0079, -0.0978,  0.0844, -0.1273,&
          0.0699, -0.0776,  0.2309, -0.2825,  0.0436, -0.0584,&
          0.6056, -0.6619,  0.6451, -0.6759,  0.6547, -0.6675,&
          0.6135, -0.6139,  0.2952, -0.2895,  0.2654, -0.2590,&
          0.3750, -0.3680,  0.5085, -0.5002,  0.6206, -0.6091,&
          0.6526, -0.6393,  0.6640, -0.6475,  0.6729, -0.6545,&
          0., 0., 0., 0., 0., 0./)
        !      LINES ARE ARRANGED 1-,1+,3-,3+,ETC. IN SPIN-ROTATION SPECTRUM
        real(kind=dbl), dimension(40) :: f = (/&
          118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910,&
          59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002,&
          56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685,&
          55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241,&
          53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368,&
          52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7631,&
          487.2494, 715.3932, 773.8397, 834.1453/)
        real(kind=dbl), dimension(40) :: s300 = (/&
          0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14,&
          0.3351E-14, 0.3292E-14, 0.3721E-14, 0.3891E-14,&
          0.3640E-14, 0.4005E-14, 0.3227E-14, 0.3715E-14,&
          0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14,&
          0.1391E-14, 0.1808E-14, 0.9124E-15, 0.1230E-14,&
          0.5603E-15, 0.7842E-15, 0.3228E-15, 0.4689E-15,&
          0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15,&
          0.4264E-16, 0.6899E-16, 0.1924E-16, 0.3229E-16,&
          0.8191E-17, 0.1423E-16, 0.6460E-15, 0.7047E-14,&
          0.3011E-14, 0.1826E-14, 0.1152E-13, 0.3971E-14/)
        real(kind=dbl), dimension(40) :: be = (/&
          0.009, 0.015, 0.083, 0.084, 0.212, 0.212, 0.391, 0.391, 0.626, 0.626,&
          0.915, 0.915, 1.26, 1.26, 1.66, 1.665, 2.119, 2.115, 2.624, 2.625,&
          3.194, 3.194, 3.814, 3.814, 4.484, 4.484, 5.224, 5.224, 6.004, 6.004, &
          6.844, 6.844,7.744, 7.744, 0.048, 0.044, 0.049, 0.145, 0.141, 0.145/)

        th = 300._dbl/tempK
        th1 = th-1._dbl
        b = th**x
        preswv = vapden*tempK/217._dbl
        presda = pres-preswv
        den = 0.001_dbl*(presda*b + 1.1_dbl*preswv*th)
        dfnr = wb300*den
        loc_sum = 1.6d-17*freq*freq*dfnr/(th*(freq*freq + dfnr*dfnr))
        do k=1,40
            df = w300(k)*den
            y = 0.001_dbl*pres*b*(y300(k)+v(k)*th1)
            str = s300(k)*exp(-be(k)*th1)
            sf1 = (df + (freq-f(k))*y)/((freq-f(k))**2 + df*df)
            sf2 = (df - (freq+f(k))*y)/((freq+f(k))**2 + df*df)
            loc_sum = loc_sum + str*(sf1+sf2)*(freq/f(k))**2
        end do
        o2abs = 0.5034d12*loc_sum*presda*th**3/pi

        return

    end function o2abs

    !*************************************************************
    function abh2o(tempK,pres,rho,freq)
        !
        !  NAME- ABH2O    LANGUAGE- FORTRAN 77
        !
        ! PURPOSE- COMPUTE ABSORPTION COEF IN ATMOSPHERE DUE TO WATER VAPOR
        !
        !  CALLING SEQUENCE PARAMETERS-
        !    SPECIFICATIONS
        !      NAME    UNITS    I/O  DESCRIPTON            VALID RANGE
        !      tempK   KELVIN    I   TEMPERATURE
        !      pres    MILLIBAR  I   PRESSURE              .1 TO 1000
        !      rho     G/M^3    I   WATER VAPOR DENSITY
        !      freq    GHZ       I   FREQUENCY             0 TO 800
        !      abh2o   NEPERS/KM O   ABSORPTION COEFFICIENT
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
        !    DATE  SEP.3, 2012  SELF- AND FOREIGN-BROADENING PARAMETER FOR H2O LINE AT 22.235 GHz
        !                           CHANGED ACCORDING TO Liljegren et al.,2005, IEEE TGaRS
        !                       SELF- AND FOREIGN-CONTINUUM CONTRIBUTION CHANGED
        !                           ACCORDING TO Turner et al., 2009, IEEE TGaRS
        !
        !   LOCAL VARIABLES:

        use kinds, only: dbl, long

        implicit none

        real(kind=dbl), intent(in) :: tempK,pres,rho,freq
        real(kind=dbl) :: abh2o
        integer(kind=long), parameter :: nlines = 15
        integer(kind=long) :: i,j
        real(kind=dbl) :: df(2)
        real(kind=dbl) :: pvap,pda,den,ti,ti2,loc_sum,width,wsq,s,base,res,con
        !     LINE FREQUENCIES:
        real(kind=dbl), dimension(nlines) :: fl = (/&
          22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508,&
          443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360,&
          620.7008, 752.0332, 916.1712/)
        !     LINE INTENSITIES AT 300K:
        real(kind=dbl), dimension(nlines) :: s1 = (/&
          0.1310E-13, 0.2273E-11, 0.8036E-13, 0.2694E-11, 0.2438E-10,&
          0.2179E-11, 0.4624E-12, 0.2562E-10, 0.8369E-12, 0.3263E-11,&
          0.6659E-12, 0.1531E-08, 0.1707E-10, 0.1011E-08, 0.4227E-10/)
        !     T COEFF. OF INTENSITIES:
        real(kind=dbl), dimension(nlines) :: b2 = (/&
          2.144, .668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405,&
          3.597, 2.379, 2.852, .159, 2.391, .396, 1.441/)
        !     AIR-BROADENED WIDTH PARAMETERS AT 300K:
        real(kind=dbl), dimension(nlines) :: w3 = (/&
          0.002656, 0.00281, 0.0023, 0.00278, 0.00287, 0.0021, 0.00186,&
          0.00263, 0.00215, 0.00236, 0.0026, 0.00321, 0.00244, 0.00306, 0.00267/)
        !     T-EXPONENT OF AIR-BROADENING:
        real(kind=dbl), dimension(nlines) :: x = (/&
          0.69, 0.64, 0.67, 0.68, 0.54, 0.63, 0.60, 0.66, 0.66,&
          0.65, 0.69, 0.69, 0.71, 0.68, 0.70/)
        !     SELF-BROADENED WIDTH PARAMETERS AT 300K:
        real(kind=dbl), dimension(nlines) :: ws = (/&
          0.0127488, 0.01491, 0.0108, 0.0135, 0.01541, 0.0090, 0.00788,&
          0.01275, 0.00983, 0.01095, 0.01313, 0.01320, 0.01140, 0.01253, 0.01275/)
        !     T-EXPONENT OF SELF-BROADENING:
        real(kind=dbl), dimension(nlines) :: xs = (/&
          0.61, 0.85, 0.54, 0.74, 0.89, 0.52, 0.50, 0.67, 0.65, 0.64, 0.72,&
          1.0, 0.68, 0.84, 0.78/)

        if(rho <= 0.) then
            abh2o = 0._dbl
            return
        endif
        pvap = rho * tempK / 217._dbl
        pda = pres - pvap
        den = 3.335d16 * rho
        ti = 300._dbl / tempK
        ti2 = ti**2.5
        !
        !      continuum terms
        con = (5.43d-10*1.105_dbl*pda*ti**3 + 1.8d-8*0.79_dbl*pvap*ti**7.5)*pvap*freq*freq
        !
        !      add resonances
        loc_sum = 0.
        do i=1,nlines
            width = w3(i)*pda*ti**x(i) + ws(i)*pvap*ti**xs(i)
            wsq = width*width
            s = s1(i)*ti2*exp(b2(i)*(1.-ti))
            df(1) = freq - fl(i)
            df(2) = freq + fl(i)
            !  use clough's definition of local line contribution
            base = width/(562500._dbl + wsq)
            !  do for positive and negative resonances
            res = 0._dbl
            do j=1,2
                if(abs(df(j)).lt.750.) res = res + width/(df(j)**2+wsq) - base
            end do
            loc_sum = loc_sum + s*res*(freq/fl(i))**2
        end do
        abh2o = 0.3183d-4*den*loc_sum + con

        return

    end function abh2o

end module gasabs_module 
 
!- End of module header
