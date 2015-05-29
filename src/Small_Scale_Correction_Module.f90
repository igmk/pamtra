!
! Small_Scale_Correction_Module
!
! Module containing the small-scale correction procedures for the 
! CRTM implementations of FASTEM4 and FASTEM5
!
! Equation (A4) of
!
!   Liu, Q. et al. (2011) An Improved Fast Microwave Water
!     Emissivity Model, TGRSS, 49, pp1238-1250
!
! describes the fitting of the small-scale correction formulation
! given in equation (17a,b) of
!
!   Liu, Q. et al. (1998) Monte Carlo simulations of the microwave
!     emissivity of the sea surface, JGR, 103, pp24983-24989
!
! and originally in equation (30) of
!
!   Guissard,A. and P.Sobieski (1987) An approximate model
!     for the microwave brightness temperature of the sea,
!     Int.J.Rem.Sens., 8, pp1607-1627.
!
!
! CREATION HISTORY:
!       Written by:     Original FASTEM authors
!
!       Refactored by:  Paul van Delst, November 2011
!                       paul.vandelst@noaa.gov
!

MODULE Small_Scale_Correction_Module

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

  ! Science routines
  PUBLIC :: Small_Scale_Correction


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Small_Scale_Correction_Module.f90 29405 2013-06-20 20:19:52Z paul.vandelst@noaa.gov $'

  ! Literal constants
  REAL(dbl), PARAMETER :: ZERO = 0.0_dbl
  REAL(dbl), PARAMETER :: ONE  = 1.0_dbl
  REAL(dbl), PARAMETER :: TWO  = 2.0_dbl
  ! Minimum and maximum frequency
  REAL(dbl), PARAMETER :: MIN_FREQUENCY = 1.4_dbl
  REAL(dbl), PARAMETER :: MAX_FREQUENCY = 200.0_dbl
  ! Minimum and maximum wind speed
  REAL(dbl), PARAMETER :: MIN_WIND_SPEED = 0.3_dbl
  REAL(dbl), PARAMETER :: MAX_WIND_SPEED = 35.0_dbl

  real(dbl), parameter, dimension(8) :: SSCCoeff = (/&
     -5.020848E-006_dbl, 2.3297951E-008_dbl, 4.6625726E-008_dbl, -1.9765665E-009_dbl, &
     -7.0469823E-004_dbl, 7.5061193E-004_dbl, 9.8103876E-004_dbl, 1.54895E-004_dbl &
    /)
  ! --------------------------------------
  ! Structure definition to hold internal
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
!  TYPE :: iVar_type
!   PRIVATE
    ! Direct inputs
!    REAL(dbl) :: wind_speed = ZERO
!    REAL(dbl) :: frequency  = ZERO
    ! Flag to indicate if wind spped outside limits
    ! Final result (since equation is an exponential)
!    REAL(dbl) :: correction = ZERO
!  END TYPE iVar_type


CONTAINS


  ! =============================================================
  ! Procedures to compute the reflectivity small scale correction
  ! =============================================================
  ! Forward model
  SUBROUTINE Small_Scale_Correction( &
    Frequency , &  ! Input
    cos_Z     , &  ! Input
    Wind_Speed, &  ! Input
    Correction &  ! Output
     )  ! Internal variable output
    ! Arguments

    REAL(dbl),               INTENT(IN)     :: Frequency
    REAL(dbl),               INTENT(IN)     :: cos_Z
    REAL(dbl),               INTENT(IN)     :: Wind_Speed
    REAL(dbl),               INTENT(OUT)    :: Correction
    ! Local variables
    REAL(dbl) :: w2, freq_loc, wsp_loc
    LOGICAL  :: wind_speed_limited = .FALSE.
    ! Intermediate variables
    REAL(dbl) :: f2     = ZERO
    REAL(dbl) :: y      = ZERO
    REAL(dbl) :: cos2_z = ZERO
    
    character(80) :: nameOfRoutine = 'Small_Scale_Correction'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Check input
    freq_loc = frequency
    IF ( frequency < MIN_FREQUENCY ) freq_loc = MIN_FREQUENCY
    IF ( frequency > MAX_FREQUENCY ) freq_loc = MAX_FREQUENCY

    wind_speed_limited = .FALSE.
    wsp_loc = wind_speed
    IF ( wind_speed < MIN_WIND_SPEED ) THEN
      wsp_loc        = MIN_WIND_SPEED
      wind_speed_limited = .TRUE.
    END IF
    IF ( wind_speed > MAX_WIND_SPEED ) THEN
      wsp_loc         = MAX_WIND_SPEED
      wind_speed_limited = .TRUE.
    END IF

    ! Compute correction
    f2 = freq_loc**2
    w2 = wsp_loc**2
    ! ...Intermediate term regression equation   
    y = &
      (SSCCoeff(1) * wsp_loc    * freq_loc) + &
      (SSCCoeff(2) * wsp_loc    * f2       ) + &
      (SSCCoeff(3) * w2                 * freq_loc) + &
      (SSCCoeff(4) * w2                 * f2       ) + &
      (SSCCoeff(5) * w2                 / freq_loc) + &
      (SSCCoeff(6) * w2                 / f2       ) + &
      (SSCCoeff(7) * wsp_loc                    ) + &
      (SSCCoeff(8) * w2                                 )
    ! ... The correction
    cos2_z     = cos_Z**2  
    correction = EXP(-y*cos2_z)
    
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    
  END SUBROUTINE Small_Scale_Correction

END MODULE Small_Scale_Correction_Module
