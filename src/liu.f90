module liu
  !
  ! Liu Ocean Permittivity module.
  !
  ! Module contains routines to compute the complex permittivities for
  ! sea water based on
  !
  !   Liu, Q. et al. (2011): An improved fast microwave water emissivity model.
  !      IEEE Trans. Geosci. Remote Sensing, Vol. 49, No. 4, pp. 1238-1250, 2011.
  !
  !
  ! CREATION HISTORY:
  !       Written by:     Quanhua (Mark) Liu, JCSDA 30-Jul-2009
  !                       Quanhua.Liu@noaa.gov
  !
  !       Adapted to PAMTRA code: Mario Mech, February 2014 
  !                               mech@meteo.uni-koeln.de

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  use kinds, only: dbl
  use report_module
  use constants, only: pi, &
                                   eps0, &  ! Permittivity of vacuum (F/m)
                                   t_abs ! Temperature units conversion
  ! Disable implicit typing
  implicit none


  ! ------------
  ! Visibilities
  ! ------------
  private

  ! ... Procedures
  public :: liu_ocean_permittivity

  ! Literal constants
  ! -----------------
  real(dbl), parameter :: zero   = 0.0_dbl
  real(dbl), parameter :: one    = 1.0_dbl
  real(dbl), parameter :: two    = 2.0_dbl
  real(dbl), parameter :: three  = 3.0_dbl
  real(dbl), parameter :: four   = 4.0_dbl

  ! Scaling factors
  ! ---------------
  real(dbl), parameter :: ghz_to_hz = 1.0e09_dbl ! gigahertz   -> hertz

  ! Fixed value for ionic conductivity denominator term, scaled
  ! for frequency. See the last term of eqn.(3) in reference.
  ! -----------------------------------------------------------
  real(dbl), parameter :: tau0 = two*pi*eps0*ghz_to_hz
  
  ! Parameters for the Liu et al (2011) permittivity model
  ! ------------------------------------------------------
  ! The coefficients for the high-frequency permittivity temperature
  ! polynomial. eqn.(4a) in reference.
  real(dbl), parameter :: einf_coeff(0:1) = (/ 3.8_dbl, &
                                              2.48033e-02_dbl /)
  
  ! The coefficients for the static permittivity temperature
  ! and salinity polynomials. Eqn.(4b) in reference.
  real(dbl), parameter :: es_t_coeff(0:3) = (/ 87.9181727_dbl      , &
                                              -4.031592248e-01_dbl, &
                                               9.493088010e-04_dbl, &
                                              -1.930858348e-06_dbl /)
  real(dbl), parameter :: es_s_coeff(0:2) = (/-2.697e-03_dbl, &
                                             -7.3e-06_dbl  , &
                                             -8.9e-06_dbl   /)

  ! The coefficients for the intermediate frequency permittivity
  ! temperature and salinity polynomials. Eqn.(4c) in reference.
  real(dbl), parameter :: e1_t_coeff(0:2) = (/ 5.723_dbl     , &
                                              2.2379e-02_dbl, &
                                             -7.1237e-04_dbl /)
  real(dbl), parameter :: e1_s_coeff(0:2) = (/-6.28908e-03_dbl, &
                                              1.76032e-04_dbl, &
                                             -9.22144e-05_dbl /)

  ! The coefficients for the relaxation time temperature and
  ! salinity polynomials. Eqns.(4e) and (4f) in reference.
  real(dbl), parameter :: tau1_t_coeff(0:3) = (/ 1.124465e-01_dbl , &
                                               -3.9815727e-03_dbl, &
                                                8.113381e-05_dbl , &
                                               -7.1824242e-07_dbl /)
  real(dbl), parameter :: tau1_s_coeff(0:2) = (/-2.39357e-03_dbl, &  
                                                3.1353e-05_dbl , &  
                                               -2.52477e-07_dbl /)  

  real(dbl), parameter :: tau2_t_coeff(0:3) = (/ 3.049979018e-03_dbl, &
                                               -3.010041629e-05_dbl, &
                                                4.811910733e-06_dbl, &
                                               -4.259775841e-08_dbl /)
  real(dbl), parameter :: tau2_s_coeff(0:2) = (/ 1.49e-01_dbl, &
                                               -8.8e-04_dbl , &
                                               -1.05e-04_dbl /)

  ! The coefficients for the ionic conductivity exponential.
  ! Eqn.(4g) in reference.
  real(dbl), parameter :: alpha_coeff = -4.259775841e-08_dbl

  ! The coefficients for the ionic conductivity exponent term
  ! polynomial. Eqn.(4i) in reference.
  real(dbl), parameter :: beta_coeff(0:5) = (/ 2.033e-02_dbl, &
                                              1.266e-04_dbl, &
                                              2.464e-06_dbl, &
                                             -1.849e-05_dbl, &
                                              2.551e-07_dbl, &
                                             -2.551e-08_dbl /)
     
  ! The coefficients for the ionic conductivity at 25C polynomial.
  ! Eqn.(4j) in reference.
  real(dbl), parameter :: alpha25_coeff(0:3) = (/ 1.82521e-01_dbl, &
                                                -1.46192e-03_dbl, &
                                                 2.09324e-05_dbl, &
                                                -1.28205e-07_dbl /)
contains
!--------------------------------------------------------------------------------
!
! NAME:
!       liu_ocean_permittivity
!
! PURPOSE:
!       Subroutine to compute ocean permittivity according to the reference,
!         Liu, Q. et al. (2011) An improved fast microwave water emissivity model.
!            IEEE Trans. Geosci. Remote Sensing, Vol. 49, No. 4, pp. 1238-1250, 2011.
!
! CALLING SEQUENCE:
!       call liu_ocean_permittivity( temperature , & ! input
!                                    salinity    , & ! input
!                                    frequency   , & ! input
!                                    permittivity & ! output
!                                    )
!
! INPUTS:
!       temperature:   Sea surface temperature
!                      units:      kelvin (k)
!                      type:       real(dbl)
!                      dimension:  scalar
!                      attributes: intent(in)
!
!       salinity:      Water salinity
!                      units:      ppt (parts per thousand)
!                      type:       real(dbl)
!                      dimension:  scalar
!                      attributes: intent(in)
!
!       frequency:     Frequency
!                      units:      ghz
!                      type:       real(dbl)
!                      dimension:  scalar
!                      attributes: intent(in)
!
! OUTPUTS:
!       permittivity:  Ocean permittivity
!                      units:      n/a
!                      type:       complex(dbl)
!                      dimension:  scalar
!                      attributes: intent(out)
!
!--------------------------------------------------------------------------------

  subroutine liu_ocean_permittivity( &
    temperature , & ! input
    salinity    , & ! input
    frequency   , & ! input
    permittivity & ! output
    )                  
    ! Arguments
    real(dbl),        intent(in)     :: temperature
    real(dbl),        intent(in)     :: salinity
    real(dbl),        intent(in)     :: frequency
    complex(dbl),     intent(out)    :: permittivity

    ! Local variables
    real(dbl) :: einf
    real(dbl) :: tau1, tau2, es, e1
    real(dbl) :: re, ie

    real(dbl) :: t=zero                             ! temperature in deg.C
    real(dbl) :: s=zero                             ! salinity
    real(dbl) :: delta=zero, beta=zero              ! ionic conductivity components
    real(dbl) :: alpha25=zero, alpha=zero           ! ionic conductivity terms
    real(dbl) :: es_t=zero, es_s=zero               ! the temperature and salinity es terms
    real(dbl) :: e1_t=zero, e1_s=zero               ! the temperature and salinity e1 terms
    real(dbl) :: tau1_t=zero, tau1_s=zero           ! the temperature and salinity tau1 terms
    real(dbl) :: tau2_t=zero, tau2_s=zero           ! the temperature and salinity tau2 terms
    real(dbl) :: f1=zero, f2=zero                   ! the relaxation compound terms, f.tau
    real(dbl) :: del1=zero, del2=zero               ! the permittivity differences
    
    character(80) :: nameOfRoutine = 'liu_ocean_permittivity'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Setup
    ! -----
    ! ...Initialise imaginary component of result
    ie = zero
    ! ...Save the inputs
    t = temperature - t_abs
    s = salinity


    ! Compute the TEMPERATURE polynomial parameterisations
    ! ----------------------------------------------------
    ! ...The high-frequency permittivity temperature polynomial (eqn.4a)
    einf = einf_coeff(0) + t*einf_coeff(1)
    ! ...The static permittivity temperature polynomial (eqn.4b)
    es = es_t_coeff(0) + t*(es_t_coeff(1) + t*(es_t_coeff(2) + t*es_t_coeff(3)))
    es_t = es  ! Save it
    ! ...The intermediate frequency permittivity temperature polynomial (eqn.4c)
    e1 = e1_t_coeff(0) + t*(e1_t_coeff(1) + t*e1_t_coeff(2))
    e1_t = e1  ! Save it
    ! ...The Debye relaxation time constants temperature polynomials (eqns.4e & 4f)
    ! ...Units of tau: nanoseconds (for use with GHz frequencies)
    tau1 = tau1_t_coeff(0) + t*(tau1_t_coeff(1) + t*(tau1_t_coeff(2) + t*tau1_t_coeff(3)))
    tau1_t = tau1  ! Save it
    tau2 = tau2_t_coeff(0) + t*(tau2_t_coeff(1) + t*(tau2_t_coeff(2) + t*tau2_t_coeff(3)))
    tau2_t = tau2  ! Save it 


    ! Compute the SALINITY polynomial parameterisations
    ! -------------------------------------------------
    if ( s > zero ) then
      ! ...The temperature difference from 25C (eqn.4h) used to compute ionic conductivity.
      delta = 25.0_dbl - t
      ! ...The beta term (eqn.4i) used to compute ionic conductivity
      beta = beta_coeff(0) + delta*(beta_coeff(1) + delta*beta_coeff(2)) + &
                  (beta_coeff(3) + delta*(beta_coeff(4) + delta*beta_coeff(5)))*s
      ! ...The ionic conductivity at 25C (eqn.4j)
      alpha25 = s*(alpha25_coeff(0) + s*(alpha25_coeff(1) + s*(alpha25_coeff(2) + s*alpha25_coeff(3))))
      ! ...The ionic conductivity (eqn.4g)
      alpha = alpha25*exp(-delta*beta)
      !  ...The imaginary component dependent on ionic conductivity (eqn.3)
      ie = -alpha/(frequency*tau0)

      ! ...The static permittivity salinity polynomial (eqn.4b)
      es_s = one + s*(es_s_coeff(0) + s*es_s_coeff(1) + t*es_s_coeff(2))
      es = es * es_s 

      ! ...The intermediate frequency permittivity salinity polynomial (eqn.4c)
      e1_s = one + s*(e1_s_coeff(0) + s*e1_s_coeff(1) + t*e1_s_coeff(2))
      e1 = e1 * e1_s 

      ! ...The Debye relaxation time constants salinity polynomials (eqns.4e & 4f)
      ! ...Units of tau: nanoseconds (for use with GHz frequencies)
      tau1_s = one + s*(tau1_s_coeff(0) + t*(tau1_s_coeff(1) + t*tau1_s_coeff(2)))
      tau1 = tau1 * tau1_s
      tau2_s = one + s*(tau2_s_coeff(0) + t*tau2_s_coeff(1) + (s**2)*tau2_s_coeff(2))
      tau2 = tau2 * tau2_s

    end if
    
    ! Compute the complex permittivity
    ! --------------------------------
    ! ...The compound terms
    f1 = frequency*tau1  ! Note there is no GHz->Hz conversion.
    f2 = frequency*tau2  ! That is embedded in the tau values.
    del1 = es - e1
    del2 = e1 - einf
    ! ...The real part
    re = einf + del1/(one + f1**2) + del2/(one + f2**2)
    ! ...The imaginary part
    ie = -ie + del1*f1/(one + f1**2) + del2*f2/(one + f2**2)
    ! ...Combine them (e = re - j.ie)
    permittivity = cmplx(re,-ie,dbl)                                  

  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine liu_ocean_permittivity

END MODULE Liu
