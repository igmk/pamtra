!
! Fresnel
!
! Module containing routines to compute Fresnel reflectivities.
!
!
! CREATION HISTORY:
!       Written by:     Masahiro Kazumori, JCSDA
!                       Masahiro.Kazumori@noaa.gov
!       Modified by:    Paul van Delst, CIMSS/SSEC 11-Apr-2007
!                       paul.vandelst@noaa.gov

MODULE Fresnel

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE kinds, ONLY: dbl
  use report_module
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE

  PUBLIC :: Fresnel_Reflectivity
  ! -----------------
  ! Module parameters
  ! -----------------
  ! RCS Id for the module
  CHARACTER(*), PARAMETER :: MODULE_RCS_ID = &
  '$Id: Fresnel.f90 29405 2013-06-20 20:19:52Z paul.vandelst@noaa.gov $'
  REAL(dbl), PARAMETER :: ZERO   = 0.0_dbl
  REAL(dbl), PARAMETER :: POINT5 = 0.5_dbl
  REAL(dbl), PARAMETER :: ONE    = 1.0_dbl
  REAL(dbl), PARAMETER :: TWO    = 2.0_dbl


CONTAINS


!--------------------------------------------------------------------------------
!
! NAME:
!       Fresnel_Reflectivity
!
! PURPOSE:
!       Subroutine to compute Fresnel reflectivities
!
! CALLING SEQUENCE:
!       CALL Fresnel_Reflectivity( permittivity, &  ! Input
!                                  cos_i       , &  ! Input
!                                  Rv          , &  ! Output
!                                  Rh          , &  ! Output
!                                  iVar          )  ! Internal variable output
!
! INPUT ARGUMENTS:
!       permittivity:  Permittivity of medium
!                      UNITS:      Farads per metre (F.m^-1)
!                      TYPE:       COMPLEX(dbl)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!       cos_i:         Cosine of incidence angle
!                      UNITS:      N/A
!                      TYPE:       REAL(dbl)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       Rv:            Reflectivity for polarisation parallel to the 
!                      plane of incidence (i.e. vertical polarization)
!                      UNITS:      N/A
!                      TYPE:       REAL(dbl)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!       Rh:            Reflectivity for polarisation perpendicular to the
!                      plane of incidence (i.e. horizontal polarization)
!                      UNITS:      N/A
!                      TYPE:       REAL(dbl)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!       iVar:          Structure containing internal variables required for
!                      subsequent tangent-linear or adjoint model calls.
!                      The contents of this structure are NOT accessible
!                      outside of the Fresnel module.
!                      UNITS:      N/A
!                      TYPE:       TYPE(iVar_type)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
! CREATION HISTORY:
!       Written by:     Masahiro Kazumori, JCSDA
!                       Masahiro.Kazumori@noaa.gov
!       Modified by:    Paul van Delst, CIMSS/SSEC 11-Apr-2007
!                       paul.vandelst@noaa.gov
!
!--------------------------------------------------------------------------------

  SUBROUTINE Fresnel_Reflectivity( permittivity, &  ! Input
                                   cos_i       , &  ! Input
                                   Rv          , &  ! Output
                                   Rh           &  ! Output
                                   )  ! Internal variable output
    ! Arguments
    COMPLEX(dbl),     INTENT(IN)     :: permittivity
    REAL(dbl),        INTENT(IN)     :: cos_i
    REAL(dbl),        INTENT(OUT)    :: Rv
    REAL(dbl),        INTENT(OUT)    :: Rh

    ! Local variables
    COMPLEX(dbl) :: zRv ! Vertical
    COMPLEX(dbl) :: zRh ! Horizontal

    COMPLEX(dbl) :: z1, z2
    ! The real and imaginary components
    REAL(dbl)    :: rzRv,izRv  ! Vertical
    REAL(dbl)    :: rzRh,izRh  ! Horizontal

    character(80) :: nameOfRoutine = 'Fresnel_Reflectivity'
    
  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    ! Compute the complex reflectivity components
    z1 = SQRT(permittivity - ONE + (cos_i*cos_i))
    z2 = permittivity * cos_i
    zRh = (cos_i  -z1) / (cos_i  +z1)
    zRv = (z2-z1) / (z2+z1)

    ! The square of the vertical abs value
    rzRv = REAL(zRv,dbl)
    izRv = AIMAG(zRv)
    Rv = rzRv**2 + izRv**2

    ! The square of the horizontal abs value
    rzRh = REAL(zRh,dbl)
    izRh = AIMAG(zRh)
    Rh = rzRh**2 + izRh**2
    
  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  END SUBROUTINE Fresnel_Reflectivity
  
END MODULE Fresnel
