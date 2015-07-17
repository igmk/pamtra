!
! CRTM_FastemX
!
! Module containing the Fastem4/5 procedures. The difference between the Fastem4
! and Fastem5 models is realised purely through the coefficients read during CRTM
! initialisation.
!
! CREATION HISTORY:
!       Written by:     Quanhua (Mark) Liu, Quanhua.Liu@noaa.gov
!                       Stephen English, Stephen.English@metoffice.gov.uk
!                       Fuzhong Weng, Fuzhong.Weng@noaa.gov
!                       July, 2009
!
!       Refactored by:  Paul van Delst, paul.vandelst@noaa.gov
!                       October, 2010
!
!       Modified by:    Quanhua Liu, Quanhua.Liu@noaa.gov
!                       Stephen English, Stephen.English@metoffice.gov.uk
!                       August, 2011
!
!       Refactored by:  Paul van Delst, paul.vandelst@noaa.gov
!                       December, 2011
!
!       Renamed:        Separate Fastem4 and Fastem5 modules replaced with
!                       single module named "FastemX"
!                       Paul van Delst, paul.vandelst@noaa.gov
!                       March, 2012
!

MODULE CRTM_FastemX

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE kinds         , ONLY: dbl
  use report_module
  use settings, only: nstokes

  USE constants, &
    ONLY: PI,deg2rad

  USE Fresnel, &
    ONLY: Fresnel_Reflectivity

  USE Liu, &
    ONLY: Ocean_Permittivity    => Liu_Ocean_Permittivity

  USE Foam_Utility_Module, &
    ONLY: Foam_Coverage, &
          Foam_Reflectivity

  USE Small_Scale_Correction_Module, &
    ONLY: Small_Scale_Correction

  USE Large_Scale_Correction_Module, &
    ONLY: Large_Scale_Correction

  USE Reflection_Correction_Module, &
    ONLY: Reflection_Correction

  USE Azimuth_Emissivity_Module, &
    ONLY: Azimuth_Emissivity

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE

  ! Science routines
  PUBLIC :: Compute_FastemX


  ! -----------------
  ! Module parameters
  ! -----------------

  ! Stokes component information
  ! ...The vector indices
  INTEGER, PARAMETER :: Iv_IDX = 1 ! Describes vertical polarization
  INTEGER, PARAMETER :: Ih_IDX = 2 ! Describes horizontal polarization
  INTEGER, PARAMETER :: U_IDX  = 3 ! Describes plane of polarization
  INTEGER, PARAMETER :: V_IDX  = 4 ! Describes ellipticity of polarization

  ! Switch for using LUT for Fastem5 large scale correction
!  LOGICAL, PUBLIC, PARAMETER :: FASTEM5_LUT = .FALSE.

  ! Invalid indicators
  REAL(dbl), PARAMETER :: INVALID_AZIMUTH_ANGLE = -999.0_dbl
  REAL(dbl), PARAMETER :: INVALID_TRANSMITTANCE = -999.0_dbl  ! Disable non-specular correction
  real(dbl), parameter :: zero = 0.0_dbl,&
    one = 1.0_dbl,&
    two = 2.0_dbl, &
    three = 3.0_dbl, &
    point_5 = 0.5_dbl

    ! The zenith angle term
    REAL(dbl) :: cos_z = ONE
    ! The permittivity term
    COMPLEX(dbl) :: Permittivity = ZERO
    ! The Fresnel reflectivity terms
    REAL(dbl) :: Rv_Fresnel = ZERO
    REAL(dbl) :: Rh_Fresnel = ZERO
    ! Foam Terms
    REAL(dbl) :: Rv_Foam = ZERO
    REAL(dbl) :: Rh_Foam = ZERO
    REAL(dbl) :: Foam_Cover = ZERO
    ! Large scale correction reflectivities
    REAL(dbl) :: Rv_Large = ZERO
    REAL(dbl) :: Rh_Large = ZERO
    ! Small scale correction factor
    REAL(dbl) :: F_Small = ZERO
    ! Final reflectivities
    REAL(dbl) :: Rv = ZERO
    REAL(dbl) :: Rh = ZERO
    ! Azimuthal emissivity
    REAL(dbl) :: e_Azimuth(NSTOKES) = ZERO
    ! Anisotropic downward radiation correction
    REAL(dbl) :: Rv_Mod = ZERO
    REAL(dbl) :: Rh_Mod = ZERO


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_FastemX
!
! PURPOSE:
!       Subroutine to compute the Fastem4 or Fastem5 microwave sea surface
!       emissivity and reflectivity.
!
! CALLING SEQUENCE:
!       CALL Compute_FastemX( &
!              Frequency    , &  ! Input
!              Zenith_Angle , &  ! Input
!              Temperature  , &  ! Input
!              Salinity     , &  ! Input
!              Wind_Speed   , &  ! Input
!              iVar         , &  ! Internal variable output
!              Emissivity   , &  ! Output
!              Reflectivity , &  ! Output
!              Azimuth_Angle, &  ! Optional input
!              Transmittance  )  ! Optional input
!
!
! INPUTS:
!       Frequency:      Microwave frequency.
!                       UNITS:      GHz
!                       TYPE:       REAL(dbl)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Zenith_Angle:   Sensor zenith angle at the sea surface
!                       UNITS:      Degrees
!                       TYPE:       REAL(dbl)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Temperature:    Sea surface temperature
!                       UNITS:      Kelvin, K
!                       TYPE:       REAL(dbl)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Salinity:       Water salinity
!                       UNITS:      ppt (parts per thousand)
!                       TYPE:       REAL(dbl)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Wind_Speed:     Sea surface wind speed
!                       UNITS:      m/s
!                       TYPE:       REAL(dbl)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Emissivity:     The surface emissivity
!                       UNITS:      N/A
!                       TYPE:       REAL(dbl)
!                       DIMENSION:  Rank-1, 4-elements (n_Stokes)
!                       ATTRIBUTES: INTENT(OUT)
!
!       Reflectivity:   The surface reflectivity.
!                       UNITS:      N/A
!                       TYPE:       REAL(dbl)
!                       DIMENSION:  Rank-1, 4-elements (n_Stokes)
!                       ATTRIBUTES: INTENT(OUT)
!
! OPTIONAL INPUTS:
!       Azimuth_Angle:  Relative azimuth angle (wind direction - sensor azimuth)
!                       UNITS:      Degrees
!                       TYPE:       REAL(dbl)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Transmittance:  Total atmospheric transmittance
!                       UNITS:      N/A
!                       TYPE:       REAL(dbl)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE Compute_FastemX( &
    Frequency    , &  ! Input
    Zenith_Angle , &  ! Input
    Temperature  , &  ! Input
    Salinity     , &  ! Input
    Wind_Speed   , &  ! Input
!    iVar         , &  ! Internal variable output
    emissivity   , &  ! Output
    reflectivity , &  ! Output
    Azimuth_Angle, &  ! Optional Input
    Transmittance  )  ! Optional Input
    ! Arguments

    REAL(dbl),                INTENT(IN)  :: Frequency
    REAL(dbl),                INTENT(IN)  :: Zenith_Angle
    REAL(dbl),                INTENT(IN)  :: Temperature
    REAL(dbl),                INTENT(IN)  :: Salinity
    REAL(dbl),                INTENT(IN)  :: Wind_Speed
    REAL(dbl), intent(out) :: Emissivity(nstokes)
    REAL(dbl), intent(out) :: Reflectivity(nstokes)
    REAL(dbl),      OPTIONAL, INTENT(IN)  :: Azimuth_Angle
    REAL(dbl),      OPTIONAL, INTENT(IN)  :: Transmittance

    character(80) :: nameOfRoutine = 'Compute_FastemX'

    real(dbl) :: emis_rt(2), refl_rt(2)
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Setup
    ! ...Save forward input variables for TL and AD calculations

     cos_z = COS(Zenith_Angle*deg2rad)

    ! Permittivity calculation
    CALL Ocean_Permittivity( Temperature, Salinity, Frequency, &
                             Permittivity)


    ! Fresnel reflectivity calculation
    CALL Fresnel_Reflectivity( Permittivity, cos_z, &
                               Rv_Fresnel, Rh_Fresnel)

    ! Foam reflectivity calculation
    CALL Foam_Reflectivity( &
           Zenith_Angle, &
           Frequency   , &
           Rv_Foam, &
           Rh_Foam  )


    ! Foam coverage calculation
    CALL Foam_Coverage( &
           Wind_Speed, &
           Foam_Cover )


      CALL Large_Scale_Correction( &
             Frequency , &
             cos_z     , &
             Wind_Speed, &
             Rv_Large  , &
             Rh_Large &
             )
!    END IF

    ! Small scale correction calculation, Var%small_corr
    CALL Small_Scale_Correction( &
           Frequency , &
           cos_z     , &
           Wind_Speed, &
           F_Small    &
           )

    ! Compute the first two Stokes components of the emissivity
    Rv = Rv_Fresnel*F_Small -Rv_Large
    Rh = Rh_Fresnel*F_Small -Rh_Large
    Emissivity(Iv_IDX) = ONE - (ONE-Foam_Cover)*Rv - Foam_Cover*Rv_Foam
    Emissivity(Ih_IDX) = ONE - (ONE-Foam_Cover)*Rh - Foam_Cover*Rh_Foam

    ! Azimuthal component calculation

    e_Azimuth = ZERO
    IF ( PRESENT(Azimuth_Angle) ) THEN
      IF ( ABS(Azimuth_Angle) <= 360.0_dbl ) THEN
        CALL Azimuth_Emissivity( &
               Wind_Speed, &
               Azimuth_Angle  , &
               Frequency , &
               cos_z     , &
               e_Azimuth  &
               )
      END IF
    END IF

    ! Anisotropic downward radiation correction calculation
    Rv_Mod = ONE
    Rh_Mod = ONE
    IF ( PRESENT(Transmittance) ) THEN
      IF ( Transmittance > ZERO .AND. Transmittance < ONE ) THEN
        CALL Reflection_Correction( &
               Frequency , &
               cos_z     , &
               Wind_Speed, &
               Transmittance  , &
               Rv_Mod    , &
               Rh_Mod  &
                )
      END IF
    END IF

    ! Assemble the return...
    ! ...emissivities
    Emissivity(Iv_IDX) = Emissivity(Iv_IDX) + e_Azimuth(Iv_IDX)
    Emissivity(Ih_IDX) = Emissivity(Ih_IDX) + e_Azimuth(Ih_IDX)
    if (nstokes .gt. 2) then
      Emissivity(U_IDX)  = e_Azimuth(U_IDX)
      Emissivity(V_IDX)  = e_Azimuth(V_IDX)
    end if
    ! ...reflectivties

    Reflectivity(Iv_IDX)      = Rv_Mod * (ONE-Emissivity(Iv_IDX))
    Reflectivity(Ih_IDX)      = Rh_Mod * (ONE-Emissivity(Ih_IDX))
    if (nstokes .gt. 2) Reflectivity(U_IDX:V_IDX) = ZERO   ! 3rd, 4th Stokes from atmosphere are not included.

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  END SUBROUTINE Compute_FastemX
!
END MODULE CRTM_FastemX
!
