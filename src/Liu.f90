!
! Liu Ocean Permittivity module.
!
! Module containing routines to compute the complex permittivities for
! sea water based on
!
!   Liu, Q. et al. (2010) An improved fast microwave water emissivity model.
!      IEEE Trans. Geosci. Remote Sensing, accepted June 25, 2010
!
!
! CREATION HISTORY:
!       Written by:     Quanhua (Mark) Liu, JCSDA 30-Jul-2009
!                       Quanhua.Liu@noaa.gov

MODULE Liu

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE kinds, ONLY: dbl
  use report_module
  USE constants, ONLY: PI, &
                                   eps0, &  ! Permittivity of vacuum (F/m)
                                   t_abs ! Temperature units conversion
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE

  ! ... Procedures
  PUBLIC :: Liu_Ocean_Permittivity

  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Liu.f90 29405 2013-06-20 20:19:52Z paul.vandelst@noaa.gov $'

  ! Literal constants
  ! -----------------
  REAL(dbl), PARAMETER :: ZERO   = 0.0_dbl
  REAL(dbl), PARAMETER :: ONE    = 1.0_dbl
  REAL(dbl), PARAMETER :: TWO    = 2.0_dbl
  REAL(dbl), PARAMETER :: THREE  = 3.0_dbl
  REAL(dbl), PARAMETER :: FOUR   = 4.0_dbl

  ! Scaling factors
  ! ---------------
  REAL(dbl), PARAMETER :: GHZ_TO_HZ = 1.0e+09_dbl ! Gigahertz   -> Hertz

  ! Fixed value for ionic conductivity denominator term, scaled
  ! for frequency. See the last term of eqn.(3) in reference.
  ! -----------------------------------------------------------
  REAL(dbl), PARAMETER :: TAU0 = TWO*PI*eps0*GHZ_TO_HZ
  
  ! Parameters for the Liu et al (2010) permittivity model
  ! ------------------------------------------------------
  ! The coefficients for the high-frequency permittivity temperature
  ! polynomial. Eqn.(4a) in reference.
  REAL(dbl), PARAMETER :: EINF_COEFF(0:1) = (/ 3.8_dbl, &
                                              2.48033e-02_dbl /)
  
  ! The coefficients for the static permittivity temperature
  ! and salinity polynomials. Eqn.(4b) in reference.
  REAL(dbl), PARAMETER :: ES_T_COEFF(0:3) = (/ 87.9181727_dbl      , &
                                              -4.031592248e-01_dbl, &
                                               9.493088010e-04_dbl, &
                                              -1.930858348E-06_dbl /)
  REAL(dbl), PARAMETER :: ES_S_COEFF(0:2) = (/-2.697e-03_dbl, &
                                             -7.3E-06_dbl  , &
                                             -8.9E-06_dbl   /)

  ! The coefficients for the intermediate frequency permittivity
  ! temperature and salinity polynomials. Eqn.(4c) in reference.
  REAL(dbl), PARAMETER :: E1_T_COEFF(0:2) = (/ 5.723_dbl     , &
                                              2.2379e-02_dbl, &
                                             -7.1237e-04_dbl /)
  REAL(dbl), PARAMETER :: E1_S_COEFF(0:2) = (/-6.28908E-03_dbl, &
                                              1.76032E-04_dbl, &
                                             -9.22144E-05_dbl /)

  ! The coefficients for the relaxation time temperature and
  ! salinity polynomials. Eqns.(4e) and (4f) in reference.
  REAL(dbl), PARAMETER :: TAU1_T_COEFF(0:3) = (/ 1.124465e-01_dbl , &
                                               -3.9815727e-03_dbl, &
                                                8.113381e-05_dbl , &
                                               -7.1824242e-07_dbl /)
  REAL(dbl), PARAMETER :: TAU1_S_COEFF(0:2) = (/-2.39357E-03_dbl, &  
                                                3.1353E-05_dbl , &  
                                               -2.52477E-07_dbl /)  

  REAL(dbl), PARAMETER :: TAU2_T_COEFF(0:3) = (/ 3.049979018e-03_dbl, &
                                               -3.010041629E-05_dbl, &
                                                4.811910733E-06_dbl, &
                                               -4.259775841E-08_dbl /)
  REAL(dbl), PARAMETER :: TAU2_S_COEFF(0:2) = (/ 1.49e-01_dbl, &
                                               -8.8E-04_dbl , &
                                               -1.05E-04_dbl /)

  ! The coefficients for the ionic conductivity exponential.
  ! Eqn.(4g) in reference.
  REAL(dbl), PARAMETER :: ALPHA_COEFF = -4.259775841E-08_dbl

  ! The coefficients for the ionic conductivity exponent term
  ! polynomial. Eqn.(4i) in reference.
  REAL(dbl), PARAMETER :: BETA_COEFF(0:5) = (/ 2.033E-02_dbl, &
                                              1.266E-04_dbl, &
                                              2.464E-06_dbl, &
                                             -1.849E-05_dbl, &
                                              2.551E-07_dbl, &
                                             -2.551E-08_dbl /)
     
  ! The coefficients for the ionic conductivity at 25C polynomial.
  ! Eqn.(4j) in reference.
  REAL(dbl), PARAMETER :: ALPHA25_COEFF(0:3) = (/ 1.82521e-01_dbl, &
                                                -1.46192E-03_dbl, &
                                                 2.09324E-05_dbl, &
                                                -1.28205E-07_dbl /)


CONTAINS


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Liu_Ocean_Permittivity
!
! PURPOSE:
!       Subroutine to compute ocean permittivity according to the reference,
!         Liu, Q. et al. (2010) An improved fast microwave water emissivity model.
!            IEEE Trans. Geosci. Remote Sensing, accepted June 25, 2010
!
! CALLING SEQUENCE:
!       CALL Liu_Ocean_Permittivity( Temperature , & ! Input
!                                    Salinity    , & ! Input
!                                    Frequency   , & ! Input
!                                    Permittivity, & ! Output
!                                    iVar          ) ! Internal variable output
!
! INPUTS:
!       Temperature:   Sea surface temperature
!                      UNITS:      Kelvin (K)
!                      TYPE:       REAL(dbl)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!       Salinity:      Water salinity
!                      UNITS:      ppt (parts per thousand)
!                      TYPE:       REAL(dbl)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!       Frequency:     Frequency
!                      UNITS:      GHz
!                      TYPE:       REAL(dbl)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Permittivity:  Ocean permittivity
!                      UNITS:      N/A
!                      TYPE:       COMPLEX(dbl)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!       iVar:          Structure containing internal variables required for
!                      subsequent tangent-linear or adjoint model calls.
!                      The contents of this structure are NOT accessible
!                      outside of this module.
!                      UNITS:      N/A
!                      TYPE:       TYPE(iVar_type)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE Liu_Ocean_Permittivity( &
    Temperature , & ! Input
    Salinity    , & ! Input
    Frequency   , & ! Input
    Permittivity & ! Output
    ) ! Internal variable output                   
    ! Arguments
    REAL(dbl),        INTENT(IN)     :: Temperature
    REAL(dbl),        INTENT(IN)     :: Salinity
    REAL(dbl),        INTENT(IN)     :: Frequency
    COMPLEX(dbl),     INTENT(OUT)    :: Permittivity
!    TYPE(iVar_type), INTENT(IN OUT) :: iVar
    ! Local variables
    REAL(dbl) :: einf
    REAL(dbl) :: tau1, tau2, es, e1
    REAL(dbl) :: re, ie

    REAL(dbl) :: t=ZERO                             ! Temperature in deg.C
    REAL(dbl) :: s=ZERO                             ! Salinity
    REAL(dbl) :: delta=ZERO, beta=ZERO              ! Ionic conductivity components
    REAL(dbl) :: alpha25=ZERO, alpha=ZERO           ! Ionic conductivity terms
    REAL(dbl) :: es_t=ZERO, es_s=ZERO               ! The temperature and salinity es terms
    REAL(dbl) :: e1_t=ZERO, e1_s=ZERO               ! The temperature and salinity e1 terms
    REAL(dbl) :: tau1_t=ZERO, tau1_s=ZERO           ! The temperature and salinity tau1 terms
    REAL(dbl) :: tau2_t=ZERO, tau2_s=ZERO           ! The temperature and salinity tau2 terms
    REAL(dbl) :: f1=ZERO, f2=ZERO                   ! The relaxation compound terms, f.tau
    REAL(dbl) :: del1=ZERO, del2=ZERO               ! The permittivity differences
    
    character(80) :: nameOfRoutine = 'Liu_Ocean_Permittivity'
    
  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Setup
    ! -----
    ! ...Initialise imaginary component of result
    ie = ZERO
    ! ...Save the inputs
    t = Temperature - t_abs
    s = Salinity


    ! Compute the TEMPERATURE polynomial parameterisations
    ! ----------------------------------------------------
    ! ...The high-frequency permittivity temperature polynomial (eqn.4a)
    einf = EINF_COEFF(0) + t*EINF_COEFF(1)
    ! ...The static permittivity temperature polynomial (eqn.4b)
    es = ES_T_COEFF(0) + t*(ES_T_COEFF(1) + &
                           t*(ES_T_COEFF(2) + &
                             t*ES_T_COEFF(3)))
    es_t = es  ! Save it
    ! ...The intermediate frequency permittivity temperature polynomial (eqn.4c)
    e1 = E1_T_COEFF(0) + t*(E1_T_COEFF(1) + &
                           t*E1_T_COEFF(2))
    e1_t = e1  ! Save it
    ! ...The Debye relaxation time constants temperature polynomials (eqns.4e & 4f)
    ! ...Units of tau: nanoseconds (for use with GHz frequencies)
    tau1 = TAU1_T_COEFF(0) + t*(TAU1_T_COEFF(1) + &
                               t*(TAU1_T_COEFF(2) + &
                                 t*TAU1_T_COEFF(3)))
    tau1_t = tau1  ! Save it
    tau2 = TAU2_T_COEFF(0) + t*(TAU2_T_COEFF(1) + &
                               t*(TAU2_T_COEFF(2) + &
                                 t*TAU2_T_COEFF(3)))
    tau2_t = tau2  ! Save it 


    ! Compute the SALINITY polynomial parameterisations
    ! -------------------------------------------------
    IF ( s > ZERO ) THEN
      ! ...The temperature difference from 25C (eqn.4h) used to compute ionic conductivity.
      delta = 25.0_dbl - t
      ! ...The beta term (eqn.4i) used to compute ionic conductivity
      beta = BETA_COEFF(0) + delta*(BETA_COEFF(1) + &
                                    delta*BETA_COEFF(2)) + &
                  (BETA_COEFF(3) + delta*(BETA_COEFF(4) + &
                                     delta*BETA_COEFF(5)))*S
      ! ...The ionic conductivity at 25C (eqn.4j)
      alpha25 = s*(ALPHA25_COEFF(0) + &
                       s*(ALPHA25_COEFF(1) + &
                         s*(ALPHA25_COEFF(2) + &
                           s*ALPHA25_COEFF(3))))
      ! ...The ionic conductivity (eqn.4g)
      alpha = alpha25*EXP(-delta*beta)
      !  ...The imaginary component dependent on ionic conductivity (eqn.3)
      ie = -alpha/(Frequency*TAU0)


      ! ...The static permittivity salinity polynomial (eqn.4b)
      es_s = ONE + s*(ES_S_COEFF(0) + s*ES_S_COEFF(1) + t*ES_S_COEFF(2))
      es = es * es_s 


      ! ...The intermediate frequency permittivity salinity polynomial (eqn.4c)
      e1_s = ONE + s*(E1_S_COEFF(0) + s*E1_S_COEFF(1) + t*E1_S_COEFF(2))
      e1 = e1 * e1_s 


      ! ...The Debye relaxation time constants salinity polynomials (eqns.4e & 4f)
      ! ...Units of tau: nanoseconds (for use with GHz frequencies)
      tau1_s = ONE + s*(TAU1_S_COEFF(0) + t*(TAU1_S_COEFF(1) + &
                                                      t*TAU1_S_COEFF(2)))
      tau1 = tau1 * tau1_s
      tau2_s = ONE + s*(TAU2_S_COEFF(0) + t*TAU2_S_COEFF(1) + &
                                                    (s**2)*TAU2_S_COEFF(2))
      tau2 = tau2 * tau2_s

    END IF
    
    
    ! Compute the complex permittivity
    ! --------------------------------
    ! ...The compound terms
    f1 = Frequency*tau1  ! Note there is no GHz->Hz conversion.
    f2 = Frequency*tau2  ! That is embedded in the tau values.
    del1 = es - e1
    del2 = e1 - einf
    ! ...The real part
    re = einf + del1/(ONE + f1**2) + &
                del2/(ONE + f2**2)
    ! ...The imaginary part
    ie = -ie + del1*f1/(ONE + f1**2) + &
               del2*f2/(ONE + f2**2)
    ! ...Combine them (e = re - j.ie)
    Permittivity = CMPLX(re,-ie,dbl)                                  

  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  END SUBROUTINE Liu_Ocean_Permittivity

END MODULE Liu
