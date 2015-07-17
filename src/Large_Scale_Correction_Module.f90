!
! Large_Scale_Correction_Module
!
! Module containing the large-scale correction procedures for the
! CRTM implementations of FASTEM4 and FASTEM5
!
! Equations (A5a) and (A5b) of
!
!   Liu, Q. et al. (2011) An Improved Fast Microwave Water
!     Emissivity Model, TGRSS, 49, pp1238-1250
!
! describes the fitting of the large-scale correction formulation.
! No explicit description of the data that was fitted is given.
!
!
! CREATION HISTORY:
!       Written by:     Original FASTEM authors
!
!       Refactored by:  Paul van Delst, November 2011
!                       paul.vandelst@noaa.gov
!

MODULE Large_Scale_Correction_Module

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE kinds,   ONLY: dbl
  use report_module
  use constants, only: deg2rad, pi

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE

  ! Science routines
  PUBLIC :: Large_Scale_Correction


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Large_Scale_Correction_Module.f90 29405 2013-06-20 20:19:52Z paul.vandelst@noaa.gov $'

  ! Literal constants
  REAL(dbl), PARAMETER :: ZERO = 0.0_dbl
  REAL(dbl), PARAMETER :: ONE  = 1.0_dbl
  REAL(dbl), PARAMETER :: TWO  = 2.0_dbl
  ! Number of "final" coefficients per polarisation
  ! Corresponds with middle dimension of LSCCOEFF data
  INTEGER, PARAMETER :: N_ZCOEFFS = 6
  ! Number of look-up table dimensions
  INTEGER, PARAMETER :: N_LUTDIMS = 4

  real(dbl), parameter, dimension(3,n_zcoeffs,2) :: lsccoeff = reshape ((/&
    -5.994667E-002_dbl,   9.341346E-004_dbl,  -9.56611E-007_dbl, &
     8.360313E-002_dbl,  -1.085991E-003_dbl,   6.735338E-007_dbl, &
    -2.617296E-002_dbl,   2.864495E-004_dbl,  -1.429979E-007_dbl, &
    -5.265879E-004_dbl,   6.880275E-005_dbl,  -2.916657E-007_dbl, &
    -1.671574E-005_dbl,   1.086405E-006_dbl,  -3.632227E-009_dbl, &
     1.16194E-004_dbl,  -6.349418E-005_dbl,   2.466556E-007_dbl, &
    -2.431811E-002_dbl,  -1.03181E-003_dbl,   4.519513E-006_dbl, &
     2.868236E-002_dbl,   1.186478E-003_dbl,  -5.257096E-006_dbl, &
    -7.93339E-003_dbl,  -2.422303E-004_dbl,   1.089605E-006_dbl, &
    -1.083452E-003_dbl,  -1.788509E-005_dbl,   5.464239E-009_dbl, &
    -3.855673E-005_dbl,   9.360072E-007_dbl,  -2.639362E-009_dbl, &
     1.101309E-003_dbl,   3.599147E-005_dbl,  -1.043146E-007_dbl &
    /) ,(/3,6,2/))
CONTAINS


  ! =============================================================
  ! Procedures to compute the reflectivity large scale correction
  ! =============================================================
  ! Forward model
  SUBROUTINE Large_Scale_Correction( &
    Frequency , &  ! Input
    cos_Z     , &  ! Input
    Wind_Speed, &  ! Input
    Rv_Large  , &  ! Output
    Rh_Large  &  ! Output
      )  ! Internal variable output
    ! Arguments

    REAL(dbl),               INTENT(IN)     :: Frequency
    REAL(dbl),               INTENT(IN)     :: cos_Z
    REAL(dbl),               INTENT(IN)     :: Wind_Speed
    REAL(dbl),               INTENT(OUT)    :: Rv_Large
    REAL(dbl),               INTENT(OUT)    :: Rh_Large
    LOGICAL :: zcoeff_invalid = .TRUE.
    REAL(dbl) :: sec_z  = ZERO
    real(dbl) :: wsp_loc
    REAL(dbl) :: zcoeff_v(N_ZCOEFFS) = ZERO
    REAL(dbl) :: zcoeff_h(N_ZCOEFFS) = ZERO

    character(80) :: nameOfRoutine = 'Large_Scale_Correction'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Setup
    zcoeff_invalid = .TRUE.
    wsp_loc     = Wind_Speed
    ! the large scale correction gives suspicious results for zenith angles
    ! larger than 75 deg
    if (cos_z < cos(deg2rad*75.)) then
	sec_z = ONE/cos(75./180.*pi)
    else
	sec_z = ONE/cos_z
    end if ! compute fitting coefficients for a given frequency
    
    ! Compute the frequency polynomial coefficients
    CALL Compute_ZCoeff(LSCCoeff(:,:,1), Frequency, zcoeff_v)
    CALL Compute_ZCoeff(LSCCoeff(:,:,2), Frequency, zcoeff_h)
    zcoeff_invalid = .FALSE.

    ! Compute the reflectivity corrections
    Rv_Large = zcoeff_v(1)                              + &
               zcoeff_v(2) * sec_z                 + &
               zcoeff_v(3) * sec_z**2              + &
               zcoeff_v(4) * wsp_loc            + &
               zcoeff_v(5) * wsp_loc**2         + &
               zcoeff_v(6) * wsp_loc*sec_z

    Rh_Large = zcoeff_h(1)                              + &
               zcoeff_h(2) * sec_z                 + &
               zcoeff_h(3) * sec_z**2              + &
               zcoeff_h(4) * wsp_loc            + &
               zcoeff_h(5) * wsp_loc**2         + &
               zcoeff_h(6) * wsp_loc*sec_z

    
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  CONTAINS

    SUBROUTINE Compute_ZCoeff(coeff,frequency,zcoeff)
      REAL(dbl), INTENT(IN)  :: coeff(:,:)
      REAL(dbl), INTENT(IN)  :: frequency
      REAL(dbl), INTENT(OUT) :: zcoeff(:)
      INTEGER :: i
      
      character(80) :: nameOfRoutine = 'Compute_ZCoeff'
      
    if (verbose >= 4) call report(info,'Start of ', nameOfRoutine)
      
      DO i = 1, SIZE(zcoeff)
        zcoeff(i) = coeff(1,i) + frequency*(coeff(2,i) + frequency*coeff(3,i))
      END DO
    if (verbose >= 4) call report(info,'End of ', nameOfRoutine)

    END SUBROUTINE Compute_ZCoeff

  END SUBROUTINE Large_Scale_Correction


END MODULE Large_Scale_Correction_Module
