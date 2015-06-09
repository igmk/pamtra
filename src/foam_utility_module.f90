module foam_utility_module
  !
  ! Helper module containing the foam-related utility routines for the
  ! CRTM implementation of FASTEM4 and FASTEM5
  !
  !
  ! CREATION HISTORY:
  !       Written by:     Original FASTEM1-5 authors
  !
  !       Refactored by:  Paul van Delst, November 2011
  !                       paul.vandelst@noaa.gov
  !
  !       Adopted to PAMTRA: Mario Mech, February 2015
  !                          mech@meteo.uni-koeln.de
  !
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  use kinds, only: dbl
  use report_module

  ! Disable implicit typing
  implicit none
  ! ------------
  ! Visibilities
  ! ------------
  private
  public :: foam_coverage
  public :: foam_reflectivity
  ! -----------------
  ! Module parameters
  ! -----------------
  ! Literal constants
  real(dbl), parameter :: zero = 0.0_dbl
  real(dbl), parameter :: one  = 1.0_dbl
  real(dbl), parameter :: two  = 2.0_dbl
  
  real(dbl), parameter, dimension(6) :: frcoeff = (/ &
    0.93_dbl, -1.748e-3_dbl, -7.336e-5_dbl,&
    1.044e-7_dbl,  0.40_dbl, -5.0e-2_dbl/)
  real(dbl), parameter, dimension(2) :: fccoeff = (/ &
    1.95e-5_dbl,   2.55_dbl/)

  contains
  
  ! ===================================================================
  ! Foam coverage.
  !
  !   Monahan, E.C., and O'Muircheartaigh, I.G., (1986)
  !     Whitecaps and the passive remote sensing of the ocean surface,
  !     International Journal of Remote Sensing, 7, pp627-642.
  !
  ! The neutral stability condition is used here (i.e. the difference
  ! between the skin and air temperature is assumed to be zero) so
  ! that the form of the foam coverage equation is the same as in
  ! Tang (1974) and Liu et al. (1998).
  !
  !   Liu, Q. et al. (1998) Monte Carlo simulations of the
  !     microwave emissivity of the sea surface.
  !     JGR, 103(C11), pp24983-24989
  !
  !   Tang, C. (1974) The effect of droplets in the air-sea
  !     transition zone on the sea brightness temperature.
  !     J. Phys. Oceanography, 4, pp579-593.
  !
  ! ===================================================================
  ! Forward model
  subroutine foam_coverage(wind_speed, coverage)
    real(dbl), intent(in)  :: wind_speed
    real(dbl), intent(out) :: coverage

    character(80) :: nameOfRoutine = 'foam_coverage'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
    
    if ( wind_speed < zero ) then
      coverage = zero
      return
    end if
    coverage = FCCoeff(1) * (wind_speed**fccoeff(2))
    
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine foam_coverage

  ! =============================================================
  ! Foam reflectivity
  !
  ! See section d in
  !
  !   Kazumori, M. et al. (2008) Impact Study of AMSR-E Radiances
  !     in the NCEP Global Data Assimilation System,
  !     Monthly Weather Review, 136, pp541-559
  !
  ! Function dependence is on zenith angle only.
  ! =============================================================
  subroutine foam_reflectivity(zenith_angle, frequency, Rv, Rh)
    ! Arguments

    real(dbl)              , intent(in)  :: zenith_angle
    real(dbl)              , intent(in)  :: frequency
    real(dbl)              , intent(out) :: Rv, Rh
    ! Local variables
    real(dbl) :: factor
    
    character(80) :: nameOfRoutine = 'Foam_Reflectivity'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! The vertical component is a fixed value
    Rv = one - frcoeff(1)  ! Fixed nadir emissivity
    
    ! The horizontal component uses a regression equation
    ! to compute a factor modifying the nadir emissivity
    factor = one + zenith_angle*(frcoeff(2) + &
                     zenith_angle*(frcoeff(3) + &
                       zenith_angle*frcoeff(4)  )  )
    Rh = one - factor*frcoeff(1)
    
    ! Frequency correction
    factor = frcoeff(5) * exp(frcoeff(6)*frequency)
    Rv = Rv * factor
    Rh = Rh * factor
    
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine foam_reflectivity  

end module foam_utility_module
