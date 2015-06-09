module fresnel
  !
  ! fresnel
  !
  ! Module containing routines to compute Fresnel reflectivities.
  !
  !
  ! CREATION HISTORY:
  !       Written by:     Masahiro Kazumori, JCSDA
  !                       Masahiro.Kazumori@noaa.gov
  !       Modified by:    Paul van Delst, CIMSS/SSEC 11-Apr-2007
  !                       paul.vandelst@noaa.gov
  !       Adopted to PAMTRA: Mario Mech, February 2015
  !                          mech@meteo.uni-koeln.de
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

  public :: fresnel_reflectivity
  ! -----------------
  ! Module parameters
  ! -----------------
  real(dbl), parameter :: zero   = 0.0_dbl
  real(dbl), parameter :: point5 = 0.5_dbl
  real(dbl), parameter :: one    = 1.0_dbl
  real(dbl), parameter :: two    = 2.0_dbl

contains
!--------------------------------------------------------------------------------
!
! NAME:
!       fresnel_reflectivity
!
! PURPOSE:
!       Subroutine to compute Fresnel reflectivities
!
! CALLING SEQUENCE:
!       call fresnel_reflectivity( permittivity, &  ! Input
!                                  cos_i       , &  ! Input
!                                  Rv          , &  ! Output
!                                  Rh            &  ! Output
!                                  )
!
! INPUT ARGUMENTS:
!       permittivity:  Permittivity of medium
!                      units:      farads per metre (f.m^-1)
!                      type:       complex(dbl)
!                      dimension:  scalar
!                      attributes: intent(in)
!
!       cos_i:         Cosine of incidence angle
!                      units:      n/a
!                      type:       real(dbl)
!                      dimension:  scalar
!                      attributes: intent(in)
!
! OUTPUT ARGUMENTS:
!       Rv:            Reflectivity for polarisation parallel to the 
!                      plane of incidence (i.e. vertical polarization)
!                      units:      n/a
!                      type:       real(dbl)
!                      dimension:  scalar
!                      attributes: intent(out)
!
!       Rh:            Reflectivity for polarisation perpendicular to the
!                      plane of incidence (i.e. horizontal polarization)
!                      units:      n/a
!                      type:       real(dbl)
!                      dimension:  scalar
!                      attributes: intent(out)
!
! CREATION HISTORY:
!       Written by:     Masahiro Kazumori, JCSDA
!                       Masahiro.Kazumori@noaa.gov
!       Modified by:    Paul van Delst, CIMSS/SSEC 11-Apr-2007
!                       paul.vandelst@noaa.gov
!       Adopted to PAMTRA: Mario Mech, February 2015
!                          mech@meteo.uni-koeln.de
!
!--------------------------------------------------------------------------------

  subroutine fresnel_reflectivity( permittivity, &  ! Input
                                   cos_i       , &  ! Input
                                   Rv          , &  ! Output
                                   Rh           &  ! Output
                                   )  ! Internal variable output
    ! Arguments
    complex(dbl),     intent(in)     :: permittivity
    real(dbl),        intent(in)     :: cos_i
    real(dbl),        intent(out)    :: rv
    real(dbl),        intent(out)    :: rh

    ! Local variables
    complex(dbl) :: zrv ! vertical
    complex(dbl) :: zrh ! horizontal
    complex(dbl) :: z1, z2
    ! The real and imaginary components
    real(dbl)    :: rzRv,izRv  ! Vertical
    real(dbl)    :: rzRh,izRh  ! Horizontal

    character(80) :: nameOfRoutine = 'fresnel_reflectivity'
    
  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Compute the complex reflectivity components
    z1 = sqrt(permittivity - one + (cos_i*cos_i))
    z2 = permittivity * cos_i
    zRh = (cos_i  -z1) / (cos_i  +z1)
    zRv = (z2-z1) / (z2+z1)

    ! The square of the vertical abs value
    rzRv = real(zRv,dbl)
    izRv = aimag(zRv)
    Rv = rzRv**2 + izRv**2

    ! The square of the horizontal abs value
    rzRh = real(zRh,dbl)
    izRh = aimag(zRh)
    Rh = rzRh**2 + izRh**2
    
  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine fresnel_reflectivity
  
end module fresnel
