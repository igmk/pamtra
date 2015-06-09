module small_scale_correction_module
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

  ! Science routines
  public :: small_scale_correction

  ! Literal constants
  real(dbl), parameter :: zero = 0.0_dbl
  real(dbl), parameter :: one  = 1.0_dbl
  real(dbl), parameter :: two  = 2.0_dbl
  ! minimum and maximum frequency
  real(dbl), parameter :: min_frequency = 1.4_dbl
  real(dbl), parameter :: max_frequency = 200.0_dbl
  ! minimum and maximum wind speed
  real(dbl), parameter :: min_wind_speed = 0.3_dbl
  real(dbl), parameter :: max_wind_speed = 35.0_dbl

  real(dbl), parameter, dimension(8) :: ssccoeff = (/&
     -5.020848e-006_dbl, 2.3297951e-008_dbl, 4.6625726e-008_dbl, -1.9765665e-009_dbl, &
     -7.0469823e-004_dbl, 7.5061193e-004_dbl, 9.8103876e-004_dbl, 1.54895e-004_dbl &
    /)

  contains

  ! =============================================================
  ! Procedures to compute the reflectivity small scale correction
  ! =============================================================
  ! Forward model
  subroutine small_scale_correction( &
    frequency , &  ! input
    cos_z     , &  ! input
    wind_speed, &  ! input
    correction &  ! output
     )  ! internal variable output
    ! arguments

    real(dbl), intent(in) :: frequency
    real(dbl), intent(in) :: cos_z
    real(dbl), intent(in) :: wind_speed
    real(dbl), intent(out) :: correction
    ! local variables
    real(dbl) :: w2, freq_loc, wsp_loc
    logical  :: wind_speed_limited = .false.
    ! intermediate variables
    real(dbl) :: f2     = zero
    real(dbl) :: y      = zero
    real(dbl) :: cos2_z = zero
    
    character(80) :: nameOfRoutine = 'Small_Scale_Correction'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Check input
    freq_loc = frequency
    if ( frequency < min_frequency ) freq_loc = min_frequency
    if ( frequency > max_frequency ) freq_loc = max_frequency

    wind_speed_limited = .false.
    wsp_loc = wind_speed
    if ( wind_speed < min_wind_speed ) then
      wsp_loc        = min_wind_speed
      wind_speed_limited = .true.
    end if
    if ( wind_speed > max_wind_speed ) then
      wsp_loc         = max_wind_speed
      wind_speed_limited = .true.
    end if

    ! Compute correction
    f2 = freq_loc**2
    w2 = wsp_loc**2
    ! ...Intermediate term regression equation  (eqn A4)
    y = &
      (ssccoeff(1) * wsp_loc * freq_loc) + &
      (ssccoeff(2) * wsp_loc * f2      ) + &
      (ssccoeff(3) * w2      * freq_loc) + &
      (ssccoeff(4) * w2      * f2      ) + &
      (ssccoeff(5) * w2      / freq_loc) + &
      (ssccoeff(6) * w2      / f2      ) + &
      (ssccoeff(7) * wsp_loc           ) + &
      (ssccoeff(8) * w2                )
    ! ... The correction
    cos2_z     = cos_Z**2  
    correction = exp(-y*cos2_z)
    
    if ((verbose >= 1) .and. (wind_speed_limited)) call report(warning, 'Wind speed limited in ', nameOfRoutine)
    
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    
  end subroutine small_scale_correction

end module small_scale_correction_module
