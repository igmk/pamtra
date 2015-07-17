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

MODULE Foam_Utility_Module

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE kinds     , ONLY: dbl
  use report_module

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: Foam_Coverage
  PUBLIC :: Foam_Reflectivity


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Foam_Utility_Module.f90 29405 2013-06-20 20:19:52Z paul.vandelst@noaa.gov $'

  ! Literal constants
  REAL(dbl), PARAMETER :: ZERO = 0.0_dbl
  REAL(dbl), PARAMETER :: ONE  = 1.0_dbl
  REAL(dbl), PARAMETER :: TWO  = 2.0_dbl
  
  real(dbl), parameter, dimension(6) :: frcoeff = (/ &
    0.93_dbl, -1.748e-3_dbl, -7.336e-5_dbl,&
    1.044e-7_dbl,  0.40_dbl, -5.0e-2_dbl/)
  real(dbl), parameter, dimension(2) :: fccoeff = (/ &
    1.95e-5_dbl,   2.55_dbl/)
CONTAINS


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
  ! Tang (1974) and Liu et al. (1998)..
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
  SUBROUTINE Foam_Coverage(wind_speed, coverage)
    REAL(dbl)              , INTENT(IN)  :: wind_speed
    REAL(dbl)              , INTENT(OUT) :: coverage

    character(80) :: nameOfRoutine = 'Foam_Coverage'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
    
    IF ( wind_speed < ZERO ) THEN
      coverage = ZERO
      RETURN
    END IF
    coverage = FCCoeff(1) * (wind_speed**FCCoeff(2))
    
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  END SUBROUTINE Foam_Coverage

  ! =============================================================
  ! Foam reflectivity
  !
  ! See section d in
  !
  !   Kazumori, M. et al. (2008) Impact Study of AMSR-E Radiances
  !     in the NCEP Global Data Assimilation System,
  !     Monthly Weather Review, 136, pp541-559
  !
  ! Function dependence is on zenith angle only so no TL
  ! or AD routine.
  ! =============================================================
  SUBROUTINE Foam_Reflectivity( &
    Zenith_Angle, &
    Frequency   , &
    Rv          , &
    Rh            )
    ! Arguments

    REAL(dbl)              , INTENT(IN)  :: Zenith_Angle
    REAL(dbl)              , INTENT(IN)  :: Frequency
    REAL(dbl)              , INTENT(OUT) :: Rv, Rh
    ! Local variables
    REAL(dbl) :: factor
    
    character(80) :: nameOfRoutine = 'Foam_Reflectivity'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! The vertical component is a fixed value
    Rv = ONE - FRCoeff(1)  ! Fixed nadir emissivity
    
    ! The horizontal component uses a regression equation
    ! to compute a factor modifying the nadir emissivity
    factor = ONE + Zenith_Angle*(FRCoeff(2) + &
                     Zenith_Angle*(FRCoeff(3) + &
                       Zenith_Angle*FRCoeff(4)  )  )
    Rh = ONE - factor*FRCoeff(1)
    
    ! Frequency correction
    factor = FRCoeff(5) * EXP(FRCoeff(6)*Frequency)
    Rv = Rv * factor
    Rh = Rh * factor
    
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  END SUBROUTINE Foam_Reflectivity  

END MODULE Foam_Utility_Module
