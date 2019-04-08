module fastemx
!
! fastemx
!
! Module containing the FASTEM4/5 procedures. The difference between the FASTEM4
! and FASTEM5 models is realised purely through the coefficients read during CRTM
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
!       Adopted to PAMTRA: Mario Mech, February 2015
!                          mech@meteo.uni-koeln.de
!

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  use kinds, only: dbl
  use report_module
  use settings, only: nstokes

  use constants, only: pi, deg2rad

  use fresnel, only: fresnel_reflectivity

  use liu, only: ocean_permittivity => liu_ocean_permittivity

  use foam_utility_module, only: foam_coverage, foam_reflectivity

  use small_scale_correction_module, only: small_scale_correction

  use large_scale_correction_module, only: large_scale_correction

  use reflection_correction_module, only: reflection_correction

  use azimuth_emissivity_module, only: azimuth_emissivity

  ! Disable implicit typing
  implicit none

  ! ------------
  ! Visibilities
  ! ------------
  private

  ! Science routines
  public :: compute_fastemx

  ! -----------------
  ! Module parameters
  ! -----------------

  ! Stokes component information
  ! ...The vector indices
  integer, parameter :: Iv_IDX = 1 ! Describes vertical polarization
  integer, parameter :: Ih_IDX = 2 ! Describes horizontal polarization
  integer, parameter :: U_IDX  = 3 ! Describes plane of polarization
  integer, parameter :: V_IDX  = 4 ! Describes ellipticity of polarization

  ! invalid indicators
  real(dbl), parameter :: invalid_azimuth_angle = -999.0_dbl
  real(dbl), parameter :: invalid_transmittance = -999.0_dbl  ! Disable non-specular correction
  real(dbl), parameter :: zero = 0.0_dbl,&
    one = 1.0_dbl,&
    two = 2.0_dbl, &
    three = 3.0_dbl, &
    point_5 = 0.5_dbl

  ! The zenith angle term
  real(dbl) :: cos_z = one
  ! The permittivity term
  complex(dbl) :: permittivity = zero
  ! The Fresnel reflectivity terms
  real(dbl) :: rv_fresnel = zero
  real(dbl) :: rh_fresnel = zero
  ! Foam Terms
  real(dbl) :: rv_foam = zero
  real(dbl) :: rh_foam = zero
  real(dbl) :: foam_cover = zero
  ! Large scale correction reflectivities
  real(dbl) :: rv_large = zero
  real(dbl) :: rh_large = zero
  ! Small scale correction factor
  real(dbl) :: f_small = zero
  ! Final reflectivities
  real(dbl) :: rv = zero
  real(dbl) :: rh = zero
  ! Azimuthal emissivity
  real(dbl) :: e_azimuth(nstokes) = zero
  ! Anisotropic downward radiation correction
  real(dbl) :: rv_mod = zero
  real(dbl) :: rh_mod = zero

  contains

  !################################################################################
  !################################################################################
  !##                                                                            ##
  !##                         ## PUBLIC MODULE ROUTINES ##                       ##
  !##                                                                            ##
  !################################################################################
  !################################################################################

  !--------------------------------------------------------------------------------
  ! NAME:
  !       compute_fastemx
  !
  ! PURPOSE:
  !       Subroutine to compute the FASTEM4 or FASTEM5 microwave sea surface
  !       emissivity and reflectivity.
  !
  ! CALLING SEQUENCE:
  !       call compute_fastemx( &
  !              frequency    , &  ! Input
  !              zenith_angle , &  ! Input
  !              temperature  , &  ! Input
  !              salinity     , &  ! Input
  !              wind_speed   , &  ! Input
  !              emissivity   , &  ! Output
  !              reflectivity , &  ! Output
  !              azimuth_angle, &  ! Optional input
  !              transmittance  )  ! Optional input
  !
  !
  ! INPUTS:
  !       frequency:      Microwave frequency.
  !                       units:      GHz
  !                       type:       real(dbl)
  !                       dimension:  scalar
  !                       attributes: intent(in)
  !
  !       zenith_angle:   Sensor zenith angle at the sea surface
  !                       units:      degrees
  !                       type:       real(dbl)
  !                       dimension:  scalar
  !                       attributes: intent(in)
  !
  !       temperature:    Sea surface temperature
  !                       units:      kelvin, k
  !                       type:       real(dbl)
  !                       dimension:  scalar
  !                       attributes: intent(in)
  !
  !       salinity:       Water salinity
  !                       units:      ppt (parts per thousand)
  !                       type:       real(dbl)
  !                       dimension:  scalar
  !                       attributes: intent(in)
  !
  !       wind_speed:     Sea surface wind speed (at 10 m)
  !                       units:      m/s
  !                       type:       real(dbl)
  !                       dimension:  scalar
  !                       attributes: intent(in)
  !
  ! OUTPUTS:
  !       emissivity:     The surface emissivity
  !                       units:      n/a
  !                       type:       real(dbl)
  !                       dimension:  rank-1, 4-elements (n_stokes)
  !                       attributes: intent(out)
  !
  !       reflectivity:   The surface reflectivity.
  !                       units:      n/a
  !                       type:       real(dbl)
  !                       dimension:  rank-1, 4-elements (n_stokes)
  !                       attributes: intent(out)
  !
  ! OPTIONAL INPUTS:
  !       azimuth_angle:  Relative azimuth angle (wind direction - sensor azimuth).
  !                       Definition of relative azimuth is decribed in Mätzler et al. 2006 p.254
  !                       units:      degrees
  !                       type:       real(dbl)
  !                       dimension:  scalar
  !                       attributes: intent(in), optional
  !
  !       transmittance:  Total atmospheric transmittance
  !                       units:      n/a
  !                       type:       real(dbl)
  !                       dimension:  scalar
  !                       attributes: intent(in), optional
  !
  !--------------------------------------------------------------------------------

  subroutine compute_fastemx( &
    frequency    , &  ! Input
    zenith_angle , &  ! Input
    temperature  , &  ! Input
    salinity     , &  ! Input
    wind_speed   , &  ! Input
    emissivity   , &  ! Output
    reflectivity , &  ! Output
    azimuth_angle, &  ! Optional Input
    transmittance  )  ! Optional Input
    ! Arguments

    real(dbl), intent(in)  :: frequency
    real(dbl), intent(in)  :: zenith_angle
    real(dbl), intent(in)  :: temperature
    real(dbl), intent(in)  :: salinity
    real(dbl), intent(in)  :: wind_speed
    real(dbl), intent(out) :: emissivity(nstokes)
    real(dbl), intent(out) :: reflectivity(nstokes)
    real(dbl), optional, intent(in)  :: azimuth_angle
    real(dbl), optional, intent(in)  :: transmittance

    character(80) :: nameOfRoutine = 'compute_fastemx'
    
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    cos_z = cos(zenith_angle*deg2rad)

    ! Permittivity calculation
    call ocean_permittivity(temperature, salinity, frequency, permittivity)

    ! Fresnel reflectivity calculation
    call fresnel_reflectivity(permittivity, cos_z, rv_fresnel, rh_fresnel)

    ! Foam reflectivity calculation
    call foam_reflectivity(zenith_angle, frequency, rv_foam, rh_foam)

    ! Foam coverage calculation
    call foam_coverage(wind_speed, foam_cover)

    call large_scale_correction(frequency, cos_z, wind_speed, rv_large, rh_large)

    ! Small scale correction calculation, Var%small_corr
    call small_scale_correction(frequency, cos_z, wind_speed, f_small)

    ! Compute the first two Stokes components of the emissivity
    ! For angles larger 70 we use only Fresnel + foam cover correction
    ! Large and small scale corrections are not parameterized for such 
    ! large angles

    if (zenith_angle < 70.) then
        Rv = rv_fresnel*f_small -rv_large
        Rh = rh_fresnel*f_small -rh_large
    else
        Rv = rv_fresnel
        Rh = rh_fresnel
    end if
    emissivity(iv_idx) = one - (one-foam_cover)*Rv - foam_cover*rv_foam
    emissivity(ih_idx) = one - (one-foam_cover)*Rh - foam_cover*rh_foam

    ! Azimuthal component calculation

    e_azimuth = zero
    if ( present(azimuth_angle) ) then
      if ( abs(azimuth_angle) <= 360.0_dbl ) then
        call azimuth_emissivity(wind_speed, azimuth_angle, frequency, cos_z, e_azimuth)
      end if
    end if

    ! Anisotropic downward radiation correction calculation
    rv_mod = one
    rh_mod = one
    if ( present(transmittance) ) then
      if ( transmittance > zero .and. transmittance < one ) then
        call reflection_correction(frequency, cos_z, wind_speed, transmittance, rv_mod, rh_mod)
      end if
    end if

    ! Assemble the return...
    ! ...emissivities
    emissivity(iv_idx) = emissivity(iv_idx) + e_azimuth(iv_idx)
    emissivity(ih_idx) = emissivity(ih_idx) + e_azimuth(ih_idx)
    if (nstokes .gt. 2) then
      emissivity(u_idx)  = e_azimuth(u_idx) ! knwon Warning: Array reference is out of bounds (3 > 2) in dimension 1
      emissivity(v_idx)  = e_azimuth(v_idx) ! knwon Warning: Array reference is out of bounds (4 > 2) in dimension 1
    end if
    ! ...reflectivties
    reflectivity(iv_idx)      = rv_mod * (one-emissivity(iv_idx))
    reflectivity(ih_idx)      = rh_mod * (one-emissivity(ih_idx))
    if (nstokes .gt. 2) reflectivity(u_idx:v_idx) = zero   ! 3rd, 4th stokes from atmosphere are not included. ! knwon Warning: Lower array reference is out of bounds (3 > 2) in dimension 1


    if (verbose >= 3) call report(info,'end of ', nameofroutine)

  end subroutine compute_fastemx

end module fastemx

