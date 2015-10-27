module slope_variance
  !
  ! Helper module containing the slope variance routines for the
  ! CRTM implementation of FASTEM4
  !
  !
  ! CREATION HISTORY:
  !       Written by:     Original FASTEM1/2/3 authors
  !
  !       Modified by:    Quanhua Liu, Quanhua.Liu@noaa.gov
  !                       Stephen English, Stephen.English@metoffice.gov.uk
  !                       July, 2009
  !
  !       Refactored by:  Paul van Delst, October 2010
  !                       paul.vandelst@noaa.gov
  !
  !       Adopted to PAMTRA: Mario Mech, February 2015
  !                          mech@meteo.uni-koeln.de
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  use kinds     , only: dbl
  use report_module
  use hyperbolic_step, only: step
  ! Disable implicit typing
  implicit none


  ! ------------
  ! Visibilities
  ! ------------
  private

  ! Science routines
  public :: compute_slope_variance
  ! Literal constants
  real(dbl), parameter :: zero = 0.0_dbl
  real(dbl), parameter :: one  = 1.0_dbl
  real(dbl), parameter :: two  = 2.0_dbl

  ! Wave slope variance parameters
  real(dbl), parameter :: var_coeffs(2) = (/0.0030_dbl, 0.00512_dbl/)  ! cox-munk coeffs
  real(dbl), parameter :: f_coeffs(3)   = (/1.0000_dbl, 0.02000_dbl, 0.30000_dbl/)

  ! Scale factor for x-input to hyperbolic step function
  real(dbl), parameter :: xscale = 10000.0_dbl

  CONTAINS

  ! =============================================
  ! Procedures to compute the wave slope variance
  ! =============================================
  ! Forward model
  subroutine compute_slope_variance( &
    frequency , &  ! input
    wind_speed, &  ! input
    variance    )  ! output
    ! Arguments
    real(dbl)       , intent(in)  :: frequency
    real(dbl)       , intent(in)  :: wind_speed
    real(dbl)       , intent(out) :: variance
    ! Local variables
    real(dbl) :: vara  = zero
    real(dbl) :: varb  = zero
    real(dbl) :: fterm = zero
    real(dbl) :: x1 = zero, g1 = zero, var1 = zero
    real(dbl) :: x2 = zero, g2 = zero

    real(dbl) :: c
    
    character(80) :: nameOfRoutine

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Compute the frequency term
    fterm   = (f_coeffs(2)*frequency + f_coeffs(3))
    ! Compute the slope variances
    c = (var_coeffs(1) + var_coeffs(2)*wind_speed) * f_coeffs(1)
    vara = c
    varb = c * fterm
    ! Return required value
    ! ...First IF statement
    ! ...IF ( varb >= vara ) THEN
    x1 = xscale*(varb - vara)
    call step(x1,g1)
    var1 = vara*g1 + varb*(one-g1)
    ! ...Second IF statement
    ! ...IF ( varb <= ZERO ) THEN
    x2 = xscale*varb
    call step(x2,g2)
    variance = var1*g2

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

    end subroutine compute_slope_variance

end module slope_variance
