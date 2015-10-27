!
! Hyperbolic step functions for differentiable replacement of IF statements
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 09-Nov-2010
!                       paul.vandelst@noaa.gov
!

MODULE Hyperbolic_Step

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
  PUBLIC :: Step


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: Hyperbolic_Step.f90 29405 2013-06-20 20:19:52Z paul.vandelst@noaa.gov $'
  ! Literals
  REAL(dbl), PARAMETER :: ZERO   = 0.0_dbl
  REAL(dbl), PARAMETER :: POINT5 = 0.5_dbl
  REAL(dbl), PARAMETER :: ONE    = 1.0_dbl
  ! X-input maximum value
  REAL(dbl), PARAMETER :: XCUTOFF = 70.0_dbl

CONTAINS


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Step
!
! PURPOSE:
!       Subroutine to compute a hyperbolic, differentiable, step function:
!
!         g(x) = 0.5(1 + TANH(x))
!
!       NOTE: No input checking of the validity of the x-argument for use
!             with TANH() is done.
!
! CALLING SEQUENCE:
!       CALL Step( x, g )
!                             
!
! INPUTS:
!       x:      The function abscissa.
!               UNITS:      N/A
!               TYPE:       REAL(dbl)
!               DIMENSION:  Scalar
!               ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       g:      The hyperbolic step function value.
!               UNITS:      N/A
!               TYPE:       REAL(dbl)
!               DIMENSION:  Scalar
!               ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE Step( x, g )
    REAL(dbl), INTENT(IN)  :: x
    REAL(dbl), INTENT(OUT) :: g
    
    character(80) :: nameOfRoutine
    
    if (verbose >= 5) call report(info,'Start of ', nameOfRoutine)
    g = POINT5 * ( ONE + TANH(x) )
    if (verbose >= 5) call report(info,'End of ', nameOfRoutine)
  END SUBROUTINE Step


END MODULE Hyperbolic_Step
