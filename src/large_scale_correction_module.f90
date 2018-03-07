module large_scale_correction_module
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
  !       Adopted to PAMTRA: Mario Mech, February 2015
  !                          mech@meteo.uni-koeln.de
  !
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  use kinds,   only: dbl
  use report_module
  use constants, only: deg2rad, pi

  ! Disable implicit typing
  implicit none
  ! ------------
  ! Visibilities
  ! ------------
  private

  ! Science routines
  public :: large_scale_correction

  ! Literal constants
  real(dbl), parameter :: zero = 0.0_dbl
  real(dbl), parameter :: one  = 1.0_dbl
  real(dbl), parameter :: two  = 2.0_dbl
  ! Number of "final" coefficients per polarisation
  ! Corresponds with middle dimension of LSCCOEFF data
  integer, parameter :: n_zcoeffs = 6
  ! Number of look-up table dimensions
  integer, parameter :: n_lutdims = 4

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

  contains

  ! =============================================================
  ! Procedures to compute the reflectivity large scale correction
  ! =============================================================
  ! Forward model
  subroutine large_scale_correction( &
    frequency , &  ! input
    cos_z     , &  ! input
    wind_speed, &  ! input
    Rv_Large  , &  ! Output
    Rh_Large  &  ! Output
      )  ! Internal variable output
    ! Arguments

    real(dbl), intent(in) :: frequency
    real(dbl), intent(in) :: cos_z
    real(dbl), intent(in) :: wind_speed
    real(dbl), intent(out) :: Rv_Large
    real(dbl), intent(out) :: Rh_Large
    logical :: zcoeff_invalid = .true.
    real(dbl) :: sec_z  = zero
    real(dbl) :: wsp_loc
    real(dbl) :: zcoeff_v(n_zcoeffs) = zero
    real(dbl) :: zcoeff_h(n_zcoeffs) = zero

    character(80) :: nameOfRoutine = 'large_scale_correction'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Setup
    zcoeff_invalid = .true.
    wsp_loc     = wind_speed
    sec_z = one/cos_z
!     ! the large scale correction gives suspicious results for zenith angles
!     ! larger than 75 deg
!     if (cos_z < cos(deg2rad*75.)) then
! 	sec_z = one/cos(75./180.*pi)
!     else
! 	sec_z = one/cos_z
!     end if ! compute fitting coefficients for a given frequency
    
    ! Compute the frequency polynomial coefficients
    call compute_zcoeff(lsccoeff(:,:,1), frequency, zcoeff_v)
    call compute_zcoeff(lsccoeff(:,:,2), frequency, zcoeff_h)
    zcoeff_invalid = .false.

    ! Compute the reflectivity corrections (eqn. A5a and A5b)
    Rv_Large = zcoeff_v(1)                 + &
               zcoeff_v(2) * sec_z         + &
               zcoeff_v(3) * sec_z**2      + &
               zcoeff_v(4) * wsp_loc       + &
               zcoeff_v(5) * wsp_loc**2    + &
               zcoeff_v(6) * wsp_loc*sec_z

    Rh_Large = zcoeff_h(1)                 + &
               zcoeff_h(2) * sec_z         + &
               zcoeff_h(3) * sec_z**2      + &
               zcoeff_h(4) * wsp_loc       + &
               zcoeff_h(5) * wsp_loc**2    + &
               zcoeff_h(6) * wsp_loc*sec_z

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  contains

    subroutine compute_zcoeff(coeff,frequency,zcoeff)
      real(dbl), intent(in)  :: coeff(:,:)
      real(dbl), intent(in)  :: frequency
      real(dbl), intent(out) :: zcoeff(:)
      integer :: i
      
      character(80) :: nameOfRoutine = 'compute_zcoeff'
      
      if (verbose >= 4) call report(info,'Start of ', nameOfRoutine)
      
      do i = 1, size(zcoeff)
        zcoeff(i) = coeff(1,i) + frequency*(coeff(2,i) + frequency*coeff(3,i))
      end do

      if (verbose >= 4) call report(info,'End of ', nameOfRoutine)

    end subroutine compute_zcoeff

  end subroutine large_scale_correction

end module large_scale_correction_module
