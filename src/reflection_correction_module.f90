module reflection_correction_module
  !
  ! Helper module containing the reflection correction routines for the
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
  !       Modified by:    Quanhua Liu, Quanhua.Liu@noaa.gov
  !                       Stephen English, Stephen.English@metoffice.gov.uk
  !                       August 16, 2011
  !
  !       Refactored by:  Paul van Delst, November 2011
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
  use slope_variance , only: compute_slope_variance

  ! Disable implicit typing
  implicit none

  ! ------------
  ! Visibilities
  ! ------------
  private
  ! Science routines
  public :: reflection_correction

  ! Literal constants
  real(dbl), parameter :: zero = 0.0_dbl
  real(dbl), parameter :: one  = 1.0_dbl
  real(dbl), parameter :: two  = 2.0_dbl
    
  ! Rx_Rough regression equation parameters
  ! ...Number of predictors
  integer, parameter :: n_predictors = 7
  ! ...Number of summation loops and equation terms
  integer, parameter :: n_terms = 3

  real(dbl), parameter, dimension(3,n_predictors,2) :: rccoeff = reshape ((/&
      3.88242E-002_dbl,  0.195_dbl,      -4.25093E-002_dbl, &
   6.07698_dbl,       -3.13861_dbl,       -1.03383_dbl, &    
  -3.77867_dbl,        1.80284_dbl,       0.699556_dbl, &    
 -5.06455E-002_dbl, -0.262822_dbl,       7.03056E-002_dbl,&
   3.62055_dbl,       -1.20318_dbl,       -1.24971_dbl,    & 
  1.54014E-002_dbl,  7.59848E-002_dbl, -2.68604E-002_dbl,&
  -8.02073_dbl,        3.24658_dbl,        3.04165_dbl,   &  
  0.199277_dbl,       0.166155_dbl,       1.53272E-002_dbl,&
   3.99234_dbl,       -1.30968_dbl,      -0.874716_dbl,     &
  -1.69403_dbl,      -2.60998E-002_dbl,  0.540443_dbl,     &
 -0.282483_dbl,      -0.219994_dbl,      -2.03438E-002_dbl,&
  0.351731_dbl,        2.08641_dbl,      -0.693299_dbl,     &
  8.67861E-002_dbl,  6.1902E-002_dbl,  5.95251E-003_dbl,&
  -4.75191_dbl,      -4.30134E-002_dbl,   2.48524_dbl   &
    /),(/3,n_predictors,2/))
  
  contains

  ! ===============================================================
  ! Use the transmittance to compute anisotropic downward radiation 
  ! effect through a corrected surface reflection.
  ! ===============================================================
  
  ! Forward model
  subroutine reflection_correction( &
    frequency    , &  ! input
    cos_z        , &  ! input
    wind_speed   , &  ! input
    transmittance, &  ! input
    Rv_Mod       , &  ! Output
    Rh_Mod     &  ! Output
     )  ! Internal variable output
    ! Arguments

    real(dbl), intent(in)     :: frequency
    real(dbl), intent(in)     :: cos_z
    real(dbl), intent(in)     :: wind_speed
    real(dbl), intent(in)     :: transmittance
    real(dbl), intent(out)    :: Rv_Mod
    real(dbl), intent(out)    :: Rh_Mod

    ! Local variables
    real(dbl) :: variance
    integer :: i
    real(dbl) :: od = zero
    real(dbl) :: odx(n_terms-1)   = zero
    real(dbl) :: zx(n_predictors) = zero
    real(dbl) :: rv_rough = zero
    real(dbl) :: rh_rough = zero

    character(80) :: nameOfRoutine
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
    
    ! Compute the wave slope variance
    call compute_slope_variance( frequency, wind_speed,variance )

    ! Compute surface to space optical depth predictors
    od     = -log(transmittance) * cos_z
    odx(1) = log(od)
    odx(2) = odx(1)**2

    ! Compute effective angle predictors
    zx(1) = one
    zx(2) = variance
    zx(4) = one / cos_z
    zx(3) = zx(2) * zx(4)
    zx(5) = zx(3) * zx(3)
    zx(6) = zx(4) * zx(4)
    zx(7) = zx(2) * zx(2)

    ! Compute the rough surface reflectivity
    Rv_rough = one
    Rh_rough = one
    do i = 1, n_predictors
      rv_rough = rv_rough + zx(i) * (rccoeff(1,i,1) + &
                                                     odx(1)*rccoeff(2,i,1) + &
                                                     odx(2)*rccoeff(3,i,1)   )

      rh_rough = rh_rough + zx(i) * (rccoeff(1,i,2) + &
                                                     odx(1)*rccoeff(2,i,2) + &
                                                     odx(2)*rccoeff(3,i,2)   )
    end do

    
    ! Compute the reflectivity modifier
    rv_mod = (one - transmittance**rv_rough) / (one - transmittance)
    rh_mod = (one - transmittance**rh_rough) / (one - transmittance) 
    ! ...and save it
    rv_mod = rv_mod
    rh_mod = rh_mod 

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine reflection_correction

end module reflection_correction_module
