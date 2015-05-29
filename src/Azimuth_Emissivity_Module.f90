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

MODULE Azimuth_Emissivity_Module

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
  PUBLIC :: Azimuth_Emissivity


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Azimuth_Emissivity_Module.f90 29405 2013-06-20 20:19:52Z paul.vandelst@noaa.gov $'

  REAL(dbl), PARAMETER :: ZERO  = 0.0_dbl
  REAL(dbl), PARAMETER :: ONE   = 1.0_dbl
  REAL(dbl), PARAMETER :: TWO   = 2.0_dbl
  REAL(dbl), PARAMETER :: THREE = 3.0_dbl
  REAL(dbl), PARAMETER :: PI = 3.141592653589793238462643383279_dbl
  REAL(dbl), PARAMETER :: DEGREES_TO_RADIANS = PI / 180.0_dbl

  ! Dimensions
  ! ...Number of component predictors for harmonic coefficients
  INTEGER, PARAMETER :: N_PREDICTORS = 10
  ! ...Number of Stokes parameters
  INTEGER, PARAMETER :: N_STOKES = 4
  ! ...The number of harmonics considered in the trignometric parameterisation
  INTEGER, PARAMETER :: N_HARMONICS = 3

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
  
CONTAINS


  ! ===========================================================
  ! Compute emissivity as a function of relative azimuth angle.
  ! ===========================================================
  
  ! Forward model
  SUBROUTINE Azimuth_Emissivity(&
    Wind_Speed   , &  ! Input
    Azimuth_Angle, &  ! Input
    Frequency    , &  ! Input
    cos_z        , &  ! Input
    e_Azimuth     &  ! Output
      )  ! Internal variable output
    ! Arguments

    REAL(dbl)              , INTENT(IN)     :: Wind_Speed   
    REAL(dbl)              , INTENT(IN)     :: Azimuth_Angle
    REAL(dbl)              , INTENT(IN)     :: Frequency    
    REAL(dbl)              , INTENT(IN)     :: cos_z        
    REAL(dbl)              , INTENT(OUT)    :: e_Azimuth(:)
   
    ! Local variables
    INTEGER :: i, m
    REAL(dbl) :: phi, angle
    REAL(dbl) :: predictor(N_PREDICTORS)

    REAL(dbl) :: sec_z = ZERO
    REAL(dbl) :: cos_angle(N_HARMONICS) = ZERO
    REAL(dbl) :: sin_angle(N_HARMONICS) = ZERO
    REAL(dbl) :: trig_coeff(N_STOKES, N_HARMONICS) = ZERO
    
    character(80) :: nameOfRoutine
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
    
    ! Initialise output
    e_Azimuth = ZERO

    sec_z      = ONE/cos_z
    
    ! Convert angle
    phi = Azimuth_Angle * DEGREES_TO_RADIANS

    ! Compute the azimuth emissivity component predictors
    CALL Compute_Predictors( Wind_Speed, Frequency, sec_z, Predictor )

    ! Compute the azimuth emissivity vector
    Harmonic_Loop: DO m = 1, N_HARMONICS

      ! Compute the angles
      angle = REAL(m,dbl) * phi
      cos_angle(m) = COS(angle)
      sin_angle(m) = SIN(angle)

      ! Compute the coefficients
      DO i = 1, N_STOKES
        CALL Compute_Coefficient( &
               AZCoeff(:,i,m), &
               Predictor, &
               trig_coeff(i,m) )
      END DO

      ! Compute the emissivities
      e_Azimuth(1) = e_Azimuth(1) + trig_coeff(1,m)*cos_angle(m) ! Vertical
      e_Azimuth(2) = e_Azimuth(2) + trig_coeff(2,m)*cos_angle(m) ! Horizontal
      e_Azimuth(3) = e_Azimuth(3) + trig_coeff(3,m)*sin_angle(m) ! +/- 45deg.
      e_Azimuth(4) = e_Azimuth(4) + trig_coeff(4,m)*sin_angle(m) ! Circular

    END DO Harmonic_Loop

    ! Apply frequency correction 
    e_Azimuth = e_Azimuth * Azimuth_Freq_Correction(Frequency)

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  END SUBROUTINE Azimuth_Emissivity



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
  SUBROUTINE Compute_Predictors( &
    Wind_Speed, &  ! Input
    Frequency , &  ! Input
    sec_z     , &  ! Input
    Predictor   )  ! Output
    ! Arguments
    REAL(dbl), INTENT(IN)  :: Wind_Speed
    REAL(dbl), INTENT(IN)  :: Frequency 
    REAL(dbl), INTENT(IN)  :: sec_z     
    REAL(dbl), INTENT(OUT) :: Predictor(N_PREDICTORS) 

    character(80) :: nameOfRoutine
    
    if (verbose >= 5) call report(info,'Start of ', nameOfRoutine)
    
    ! Compute the predictors.
    Predictor( 1) = ONE
    Predictor( 2) = Frequency
    Predictor( 3) = sec_z
    Predictor( 4) = sec_z * Frequency
    Predictor( 5) = Wind_Speed
    Predictor( 6) = Wind_Speed * Frequency
    Predictor( 7) = Wind_Speed**2
    Predictor( 8) = Frequency * Wind_Speed**2
    Predictor( 9) = Wind_Speed * sec_z
    Predictor(10) = Wind_Speed * sec_z * Frequency

    if (verbose >= 5) call report(info,'End of ', nameOfRoutine)

    END SUBROUTINE Compute_Predictors

    
    
  ! ==============================================================
  ! Compute the component coefficient from the regression equation
  ! ==============================================================
  
  ! Forward model
  SUBROUTINE Compute_Coefficient( &
    c           , &  ! Input
    X           , &  ! Input
    Coefficient   )  ! Output
    ! Arguments
    REAL(dbl), INTENT(IN)  :: c(:)  ! regression coefficient
    REAL(dbl), INTENT(IN)  :: X(:)  ! predictor
    REAL(dbl), INTENT(OUT) :: Coefficient
    ! Local variables
    INTEGER :: i

    character(80) :: nameOfRoutine
    
    if (verbose >= 5) call report(info,'Start of ', nameOfRoutine)

    ! Compute component coefficient
    Coefficient = ZERO
    DO i = 1, N_PREDICTORS
      Coefficient = Coefficient + c(i)*X(i)
    END DO
    if (verbose >= 5) call report(info,'End of ', nameOfRoutine)

    END SUBROUTINE Compute_Coefficient

  PURE FUNCTION  Azimuth_Freq_Correction( Frequency ) RESULT( Fre_C )
    IMPLICIT NONE
    REAL( dbl ), INTENT(IN) :: Frequency
    REAL( dbl ) :: Fre_C
    INTEGER :: i
      ! Data for the frequency correction
    REAL(dbl), PARAMETER :: x(9) = (/ 0.0_dbl, 1.4_dbl, 6.8_dbl, 10.7_dbl, 19.35_dbl, &
                                   37._dbl, 89._dbl, 150._dbl, 200._dbl/)
    REAL(dbl), PARAMETER :: y(9) = (/ 0.0_dbl, 0.1_dbl, 0.6_dbl, 0.9_dbl, 1._dbl, &
                                   1.0_dbl, 0.4_dbl, 0.2_dbl, 0.0_dbl/)

    IF( Frequency <= ZERO .or. Frequency >= 200.0_dbl ) THEN
      Fre_C = ZERO
      RETURN
    ELSE
      DO i = 1, 8
        IF( Frequency >= x(i) .and. Frequency <= x(i+1) ) THEN
          Fre_C = y(i) + (y(i+1)-y(i))/(x(i+1)-x(i))*(Frequency-x(i))
          RETURN
        END IF
      END DO
    END IF
    Fre_C = ZERO

  END FUNCTION  Azimuth_Freq_Correction

END MODULE Azimuth_Emissivity_Module
