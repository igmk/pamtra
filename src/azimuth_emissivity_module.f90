module azimuth_emissivity_module
  !
  ! Helper module containing the azimuth emissivity routines for the
  ! CRTM implementation of FASTEM4 and FASTEM5
  !
  !
  ! CREATION HISTORY:
  !       Written by:     Original FASTEM1/2/3 authors
  !
  !       Modified by:    Quanhua Liu, Quanhua.Liu@noaa.gov
  !                       Stephen English, Stephen.English@metoffice.gov.uk
  !                       July, 2009
  !
  !       Refactored by:  Paul van Delst, December 2011
  !                       paul.vandelst@noaa.gov
  !
  !       Adopted to PAMTRA: Mario Mech, February 2015
  !                          mech@meteo.uni-koeln.de
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  use kinds, only: dbl
  use report_module
  use constants, only: pi, deg2rad

  ! Disable implicit typing
  implicit none

  ! ------------
  ! Visibilities
  ! ------------
  private

  ! Science routines
  public :: azimuth_emissivity

  real(dbl), parameter :: zero  = 0.0_dbl
  real(dbl), parameter :: one   = 1.0_dbl
  real(dbl), parameter :: two   = 2.0_dbl
  real(dbl), parameter :: three = 3.0_dbl

  ! Dimensions
  ! ...Number of component predictors for harmonic coefficients
  integer, parameter :: n_predictors = 10
  ! ...Number of Stokes parameters
  integer, parameter :: n_stokes = 4
  ! ...The number of harmonics considered in the trignometric parameterisation
  integer, parameter :: n_harmonics = 3

  real(dbl), parameter, dimension(n_predictors,n_stokes,n_harmonics) :: azcoeff = reshape ((/&
    1.318143E-2_dbl, -1.660586E-4_dbl, -7.102244E-3_dbl,  8.771616E-5_dbl, -3.418311E-3_dbl,&
    3.784895E-5_dbl, 5.763184E-5_dbl, -6.290578E-7_dbl, 1.839451E-3_dbl, -1.856317E-5_dbl, &
    6.459324E-3_dbl, -7.57005E-5_dbl, -3.777932E-3_dbl,  4.270676E-5_dbl, -1.247285E-3_dbl, &
    1.13624E-5_dbl, 2.123934E-5_dbl, -2.377368E-7_dbl, 7.070105E-4_dbl, -5.092876E-6_dbl, &
   -6.296038E-3_dbl,  3.835747E-5_dbl,  3.013694E-3_dbl, -9.366178E-6_dbl,  1.680703E-3_dbl, &
   -5.745778E-6_dbl, -2.942056E-5_dbl, 1.889216E-7_dbl, -9.058433E-4_dbl, -1.136992E-6_dbl, &
   -5.854263E-4_dbl,  5.546263E-6_dbl,  2.485058E-4_dbl, -1.531698E-6_dbl,  1.243394E-4_dbl, &
   -1.575561E-6_dbl, -2.437488E-6_dbl, 2.986237E-8_dbl, -5.5557E-5_dbl, 6.076001E-7_dbl, &
    4.605486E-3_dbl,  5.781246E-5_dbl, -2.746737E-3_dbl, -4.690045E-5_dbl,  1.512049E-4_dbl, &
    -7.411844E-9_dbl, -3.476559E-6_dbl, 1.466902E-7_dbl, -6.472364E-5_dbl, -1.776898E-6_dbl, &
   -1.863094E-2_dbl,  2.76866E-4_dbl,  7.62493E-3_dbl, -1.397481E-4_dbl,  3.550912E-3_dbl, &
   -5.533696E-5_dbl, -6.557083E-5_dbl, 9.948138E-7_dbl, -1.626538E-3_dbl, 2.307157E-5_dbl, &
   -2.880306E-2_dbl,  2.418851E-4_dbl,  1.290535E-2_dbl, -8.803702E-5_dbl,  5.057109E-6_dbl, &
   -2.715428E-5_dbl, -6.912266E-5_dbl, 7.852767E-7_dbl, 5.337096E-4_dbl, 6.585635E-6_dbl, &
    6.042016E-3_dbl, -1.135219E-4_dbl, -2.231061E-3_dbl,  5.729232E-5_dbl, -1.543391E-3_dbl, &
    2.288614E-5_dbl, 2.828443E-5_dbl, -4.384802E-7_dbl, 7.080137E-4_dbl, -9.827192E-6_dbl, &
    1.205735E-3_dbl, -1.748276E-5_dbl, -6.002919E-4_dbl,  1.174144E-5_dbl, -1.735732E-4_dbl, &
    2.148296E-6_dbl, 2.955853E-6_dbl, -3.609258E-8_dbl, 9.669164E-5_dbl, -1.282544E-6_dbl, &
   -7.610401E-4_dbl,  1.29312E-5_dbl,  3.796897E-4_dbl, -5.562741E-6_dbl,  8.865672E-5_dbl, &
   -1.313724E-6_dbl, 7.009076E-8_dbl, 2.426378E-8_dbl, -8.192732E-5_dbl, 5.333771E-7_dbl, &
   -1.834561E-3_dbl,  2.896784E-5_dbl,  7.613927E-4_dbl, -1.367783E-5_dbl,  4.887281E-4_dbl, &
   -5.81038E-6_dbl, -9.568319E-6_dbl, 1.207029E-7_dbl, -2.21079E-4_dbl, 2.159904E-6_dbl, &
   -2.054959E-4_dbl,  1.806305E-7_dbl,  1.144686E-4_dbl,  4.638982E-7_dbl,  3.581176E-5_dbl, &
   -3.870976E-7_dbl, -6.861957E-7_dbl, 6.98978E-9_dbl, -1.526136E-5_dbl, 1.887424E-7_dbl &
    /),(/n_predictors,n_stokes,n_harmonics/))
  
  contains

  ! ===========================================================
  ! Compute emissivity as a function of relative azimuth angle.
  ! ===========================================================
  
  ! Forward model
  subroutine azimuth_emissivity(&
    wind_speed   , &  ! input
    azimuth_angle, &  ! input
    frequency    , &  ! input
    cos_z        , &  ! input
    e_azimuth     &  ! output
      )  ! Internal variable output
    ! Arguments

    real(dbl), intent(in)     :: wind_speed   
    real(dbl), intent(in)     :: azimuth_angle
    real(dbl), intent(in)     :: frequency    
    real(dbl), intent(in)     :: cos_z        
    real(dbl), intent(out)    :: e_azimuth(:)
   
    ! Local variables
    integer :: i, m
    real(dbl) :: phi, angle
    real(dbl) :: predictor(n_predictors)

    real(dbl) :: sec_z = zero
    real(dbl) :: cos_angle(n_harmonics) = zero
    real(dbl) :: sin_angle(n_harmonics) = zero
    real(dbl) :: trig_coeff(n_stokes, n_harmonics) = zero
    
    character(80) :: nameOfRoutine
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
    
    ! Initialise output
    e_azimuth = zero

    sec_z      = one/cos_z
    
    ! Convert angle
    phi = azimuth_angle * deg2rad

    ! Compute the azimuth emissivity component predictors
    call compute_predictors( wind_speed, frequency, sec_z, predictor )

    ! Compute the azimuth emissivity vector
    harmonic_loop: do m = 1, n_harmonics

      ! Compute the angles
      angle = real(m,dbl) * phi
      cos_angle(m) = cos(angle)
      sin_angle(m) = sin(angle)

      ! Compute the coefficients
      do i = 1, n_stokes
        call compute_coefficient( &
               azcoeff(:,i,m), &
               predictor, &
               trig_coeff(i,m) )
      end do

      ! Compute the emissivities
      e_Azimuth(1) = e_Azimuth(1) + trig_coeff(1,m)*cos_angle(m) ! Vertical
      e_Azimuth(2) = e_Azimuth(2) + trig_coeff(2,m)*cos_angle(m) ! Horizontal
      e_Azimuth(3) = e_Azimuth(3) + trig_coeff(3,m)*sin_angle(m) ! +/- 45deg.
      e_Azimuth(4) = e_Azimuth(4) + trig_coeff(4,m)*sin_angle(m) ! Circular

    end do harmonic_loop

    ! Apply frequency correction 
    e_azimuth = e_azimuth * azimuth_freq_correction(frequency)

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine azimuth_emissivity



!################################################################################
!################################################################################
!##                                                                            ##
!##                        ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

  ! =============================================
  ! Compute predictors for the azimuth components
  ! =============================================
  
  ! Forward model
  subroutine compute_predictors( &
    wind_speed, &  ! input
    frequency , &  ! input
    sec_z     , &  ! input
    predictor   )  ! output
    ! Arguments
    real(dbl), intent(in)  :: wind_speed
    real(dbl), intent(in)  :: frequency 
    real(dbl), intent(in)  :: sec_z     
    real(dbl), intent(out) :: predictor(n_predictors) 

    character(80) :: nameOfRoutine
    
    if (verbose >= 5) call report(info,'Start of ', nameOfRoutine)
    
    ! Compute the predictors.
    predictor( 1) = one
    predictor( 2) = frequency
    predictor( 3) = sec_z
    predictor( 4) = sec_z * frequency
    predictor( 5) = wind_speed
    predictor( 6) = wind_speed * frequency
    predictor( 7) = wind_speed**2
    predictor( 8) = frequency * wind_speed**2
    predictor( 9) = wind_speed * sec_z
    predictor(10) = wind_speed * sec_z * frequency

    if (verbose >= 5) call report(info,'End of ', nameOfRoutine)

    end subroutine compute_predictors
    
  ! ==============================================================
  ! Compute the component coefficient from the regression equation
  ! ==============================================================
  
  ! Forward model
  subroutine compute_coefficient( &
    c           , &  ! input
    x           , &  ! input
    coefficient   )  ! output
    ! Arguments
    real(dbl), intent(in)  :: c(:)  ! regression coefficient
    real(dbl), intent(in)  :: x(:)  ! predictor
    real(dbl), intent(out) :: coefficient
    ! Local variables
    integer :: i

    character(80) :: nameOfRoutine
    
    if (verbose >= 5) call report(info,'Start of ', nameOfRoutine)

    ! Compute component coefficient
    coefficient = zero
    do i = 1, n_predictors
      coefficient = coefficient + c(i)*x(i)
    end do
    if (verbose >= 5) call report(info,'End of ', nameOfRoutine)

  end subroutine compute_coefficient

  pure function  azimuth_freq_correction( frequency ) result( fre_c )
    implicit none
    real( dbl ), intent(in) :: frequency
    real( dbl ) :: fre_c
    integer :: i
      ! Data for the frequency correction
    real(dbl), parameter :: x(9) = (/ 0.0_dbl, 1.4_dbl, 6.8_dbl, 10.7_dbl, 19.35_dbl, &
                                   37._dbl, 89._dbl, 150._dbl, 200._dbl/)
    real(dbl), parameter :: y(9) = (/ 0.0_dbl, 0.1_dbl, 0.6_dbl, 0.9_dbl, 1._dbl, &
                                   1.0_dbl, 0.4_dbl, 0.2_dbl, 0.0_dbl/)

    if( frequency <= zero .or. frequency >= 200.0_dbl ) then
      fre_c = zero
      return
    else
      do i = 1, 8
        if( frequency >= x(i) .and. frequency <= x(i+1) ) then
          fre_c = y(i) + (y(i+1)-y(i))/(x(i+1)-x(i))*(frequency-x(i))
          return
        end if
      end do
    end if
    fre_c = zero

  end function  azimuth_freq_correction

end module azimuth_emissivity_module
