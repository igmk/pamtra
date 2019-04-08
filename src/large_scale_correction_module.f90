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
  !                          mario.mech@.uni-koeln.de
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
  ! fastem version 4 or 5 (default 5)
  integer, parameter :: fastem_version = 5

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
  real(dbl), parameter, dimension(36) :: lcoef5 = (/&
   -5.994667E-02_dbl, 9.341346E-04_dbl,-9.566110E-07_dbl, 8.360313E-02_dbl,-1.085991E-03_dbl, &
   6.735338E-07_dbl,-2.617296E-02_dbl, 2.864495E-04_dbl,-1.429979E-07_dbl,-5.265879E-04_dbl, &
   6.880275E-05_dbl,-2.916657E-07_dbl,-1.671574E-05_dbl, 1.086405E-06_dbl,-3.632227E-09_dbl, &
   1.161940E-04_dbl,-6.349418E-05_dbl, 2.466556E-07_dbl,-2.431811E-02_dbl,-1.031810E-03_dbl, &
   4.519513E-06_dbl, 2.868236E-02_dbl, 1.186478E-03_dbl,-5.257096E-06_dbl,-7.933390E-03_dbl, &
  -2.422303E-04_dbl, 1.089605E-06_dbl,-1.083452E-03_dbl,-1.788509E-05_dbl, 5.464239E-09_dbl, &
  -3.855673E-05_dbl, 9.360072E-07_dbl,-2.639362E-09_dbl, 1.101309E-03_dbl, 3.599147E-05_dbl, &
  -1.043146E-07_dbl  /)
  real(dbl), parameter, dimension(36) :: lcoef4 = (/&
  -9.197134E-02_dbl, 8.310678E-04_dbl,-6.065411E-07_dbl, 1.350073E-01_dbl,-1.032096E-03_dbl, &
   4.259935E-07_dbl,-4.373322E-02_dbl, 2.545863E-04_dbl, 9.835554E-08_dbl,-1.199751E-03_dbl, &
   1.360423E-05_dbl,-2.088404E-08_dbl,-2.201640E-05_dbl, 1.951581E-07_dbl,-2.599185E-10_dbl, &
   4.477322E-04_dbl,-2.986217E-05_dbl, 9.406466E-08_dbl,-7.103127E-02_dbl,-4.713113E-05_dbl, &
   1.754742E-06_dbl, 9.720859E-02_dbl, 1.374668E-04_dbl,-2.591771E-06_dbl,-2.687455E-02_dbl, &
  -3.677779E-05_dbl, 7.548377E-07_dbl,-3.049506E-03_dbl,-5.412826E-05_dbl, 2.285387E-07_dbl, &
  -2.201640E-05_dbl, 1.951581E-07_dbl,-2.599185E-10_dbl, 2.297488E-03_dbl, 3.787032E-05_dbl, &
  -1.553581E-07_dbl/)

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
    integer :: j
    real(dbl), dimension(12) :: zc
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Setup
    zcoeff_invalid = .true.
    wsp_loc     = wind_speed
    sec_z = one/cos_z
    ! the large scale correction gives suspicious results for zenith angles
    ! larger than 75 deg. This is taken into account in fastemx routine
    ! if (cos_z < cos(deg2rad*75.)) then
    !   sec_z = one/cos(75./180.*pi)
    ! else
    !   sec_z = one/cos_z
    ! end if ! compute fitting coefficients for a given frequency

    if (fastem_version == 5) then
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

    elseif (fastem_version == 4) then
      Do j = 1, 12
        ! zc(j) = Lcoef5(j*3-2) + Lcoef5(j*3-1)*frequency + Lcoef5(j*3)*frequency**2
        zc(j) = Lcoef4(j*3-2) + Lcoef4(j*3-1)*frequency + Lcoef4(j*3)*frequency**2
      End do
      Rv_Large = zc(1) + zc(2)*sec_z + zc(3)*sec_z**2 + zc(4)*wsp_loc &
        + zc(5)*wsp_loc**2 + zc(6)*wsp_loc*sec_z
      Rh_Large = zc(7) + zc(8)*sec_z + zc(9)*sec_z**2 + zc(10)*wsp_loc &
        + zc(11)*wsp_loc**2 + zc(12)*wsp_loc*sec_z
    end if
    ! print*, wsp_loc, 180.0_dbl/pi*acos(cos_z), Rv_Large, Rh_Large
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
