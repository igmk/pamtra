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

MODULE Slope_Variance

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE kinds     , ONLY: dbl
  use report_module
  USE Hyperbolic_Step, ONLY: Step
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE

  ! Science routines
  PUBLIC :: Compute_Slope_Variance

  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Slope_Variance.f90 29405 2013-06-20 20:19:52Z paul.vandelst@noaa.gov $'

  ! Literal constants
  REAL(dbl), PARAMETER :: ZERO = 0.0_dbl
  REAL(dbl), PARAMETER :: ONE  = 1.0_dbl
  REAL(dbl), PARAMETER :: TWO  = 2.0_dbl

  ! Wave slope variance parameters
  REAL(dbl), PARAMETER :: VAR_COEFFS(2) = (/0.0030_dbl, 0.00512_dbl/)  ! Cox-Munk coeffs
  REAL(dbl), PARAMETER :: F_COEFFS(3)   = (/1.0000_dbl, 0.02000_dbl, 0.30000_dbl/)

  ! Scale factor for x-input to hyperbolic step function
  REAL(dbl), PARAMETER :: XSCALE = 10000.0_dbl


CONTAINS


  ! =============================================
  ! Procedures to compute the wave slope variance
  ! =============================================
  ! Forward model
  SUBROUTINE Compute_Slope_Variance( &
    Frequency , &  ! Input
    Wind_Speed, &  ! Input
    Variance    )  ! Output
    ! Arguments
    REAL(dbl)       , INTENT(IN)  :: Frequency
    REAL(dbl)       , INTENT(IN)  :: Wind_Speed
    REAL(dbl)       , INTENT(OUT) :: Variance
    ! Local variables
    REAL(dbl) :: vara  = ZERO
    REAL(dbl) :: varb  = ZERO
    REAL(dbl) :: Fterm = ZERO
    REAL(dbl) :: x1 = ZERO, g1 = ZERO, var1 = ZERO
    REAL(dbl) :: x2 = ZERO, g2 = ZERO

    REAL(dbl) :: c
    
    character(80) :: nameOfRoutine

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Compute the frequency term
    Fterm   = (F_COEFFS(2)*Frequency + F_COEFFS(3))
    ! Compute the slope variances
    c = (VAR_COEFFS(1) + VAR_COEFFS(2)*Wind_Speed) * F_COEFFS(1)
    vara = c
    varb = c * Fterm
    ! Return required value
    ! ...First IF statement
    ! ...IF ( varb >= vara ) THEN
    x1 = XSCALE*(varb - vara)
    CALL Step(x1,g1)
    var1 = vara*g1 + varb*(ONE-g1)
    ! ...Second IF statement
    ! ...IF ( varb <= ZERO ) THEN
    x2 = XSCALE*varb
    CALL Step(x2,g2)
    Variance = var1*g2

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

    END SUBROUTINE Compute_Slope_Variance

END MODULE Slope_Variance
