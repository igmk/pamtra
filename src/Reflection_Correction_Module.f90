!
! Helper module conmtaining the reflection correction routines for the
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

MODULE Reflection_Correction_Module

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE kinds     , ONLY: dbl
  use report_module
  USE Slope_Variance , ONLY: Compute_Slope_Variance

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! Science routines
  PUBLIC :: Reflection_Correction


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Reflection_Correction_Module.f90 29405 2013-06-20 20:19:52Z paul.vandelst@noaa.gov $'

  ! Literal constants
  REAL(dbl), PARAMETER :: ZERO = 0.0_dbl
  REAL(dbl), PARAMETER :: ONE  = 1.0_dbl
  REAL(dbl), PARAMETER :: TWO  = 2.0_dbl
    
  ! Rx_Rough regression equation parameters
  ! ...Number of predictors
  INTEGER, PARAMETER :: N_PREDICTORS = 7
  ! ...Number of summation loops and equation terms
  INTEGER, PARAMETER :: N_TERMS = 3

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
  
CONTAINS


  ! ===============================================================
  ! Use the transmittance to compute anisotropic downward radiation 
  ! effect through a corrected surface reflection.
  ! ===============================================================
  
  ! Forward model
  SUBROUTINE Reflection_Correction( &
    Frequency    , &  ! Input
    cos_z        , &  ! Input
    Wind_Speed   , &  ! Input
    Transmittance, &  ! Input
    Rv_Mod       , &  ! Output
    Rh_Mod     &  ! Output
     )  ! Internal variable output
    ! Arguments

    REAL(dbl)              , INTENT(IN)     :: Frequency
    REAL(dbl)              , INTENT(IN)     :: cos_Z
    REAL(dbl)              , INTENT(IN)     :: Wind_Speed
    REAL(dbl)              , INTENT(IN)     :: Transmittance
    REAL(dbl)              , INTENT(OUT)    :: Rv_Mod
    REAL(dbl)              , INTENT(OUT)    :: Rh_Mod

    ! Local variables
    REAL(dbl) :: variance
    INTEGER :: i
    REAL(dbl) :: od = ZERO
    REAL(dbl) :: odx(N_TERMS-1)   = ZERO
    REAL(dbl) :: zx(N_PREDICTORS) = ZERO
    REAL(dbl) :: Rv_Rough = ZERO
    REAL(dbl) :: Rh_Rough = ZERO

    character(80) :: nameOfRoutine
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
    
    ! Compute the wave slope variance
    CALL Compute_Slope_Variance( Frequency, Wind_Speed,Variance )


    ! Compute surface to space optical depth predictors
    od     = -LOG(Transmittance) * cos_z
    odx(1) = LOG(od)
    odx(2) = odx(1)**2


    ! Compute effective angle predictors
    zx(1) = ONE
    zx(2) = variance
    zx(4) = ONE / cos_z
    zx(3) = zx(2) * zx(4)
    zx(5) = zx(3) * zx(3)
    zx(6) = zx(4) * zx(4)
    zx(7) = zx(2) * zx(2)

    
    ! Compute the rough surface reflectivity
    Rv_Rough = ONE
    Rh_Rough = ONE
    DO i = 1, N_PREDICTORS
      Rv_Rough = Rv_Rough + zx(i) * (RCCoeff(1,i,1) + &
                                                     odx(1)*RCCoeff(2,i,1) + &
                                                     odx(2)*RCCoeff(3,i,1)   )

      Rh_Rough = Rh_Rough + zx(i) * (RCCoeff(1,i,2) + &
                                                     odx(1)*RCCoeff(2,i,2) + &
                                                     odx(2)*RCCoeff(3,i,2)   )
    END DO

    
    ! Compute the reflectivity modifier
    Rv_Mod = (ONE - Transmittance**Rv_Rough) / (ONE - Transmittance)
    Rh_Mod = (ONE - Transmittance**Rh_Rough) / (ONE - Transmittance) 
    ! ...and save it
    Rv_Mod = Rv_Mod
    Rh_Mod = Rh_Mod 

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  END SUBROUTINE Reflection_Correction

END MODULE Reflection_Correction_Module
