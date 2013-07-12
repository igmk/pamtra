!+ mass-size relation parameters calculation
subroutine make_mass_size(errorstatus)

! Description:
!  Calculate the a and b parameters of the formula relating mass (m) 
!  and diameter (D) for liquid and ice particles:
!  m = a * D^b
!  When rho_ms is provided the as_ratio is used to calculate a_ms
!  and b_ms according to this formulas:
!  if as_ratio > 1  (oblate)  b=3 a=pi/6 * rho_ms * 1/as_ratio
!  if as_ratio <= 1 (prolate) b=3 a=pi/6 * rho_ms * (as_ratio)^2
! 
! Owner: IGMK
!
! History:
!
! Version Date         Comment
! ------- ----         -------
! 0.      21/05/2013   Creation - E. Orlandi
!
! Code Description:
!  Language: Fortran 90.
!  Software Standards: "European Standards for Writing and
!   Documenting Exchangeable Fortran 90 Code".
!
! Parent Module: 
!
! Declarations:

! Modules used:

  use kinds, only: dbl, &      ! integer parameter specifying double precision
                    long       ! integer parameter specifying long integer

  use report_module

  use constants, only: pi, rho_water

  use drop_size_dist, only: rho_ms, as_ratio, liq_ice, &    ! IN
                            a_ms, b_ms                      ! IN-OUT

! Imported Scalar Variables with intent (in):

  implicit none

!- End of header ---------------------------------------------------------------

! Error handling

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'make_mass_size'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

! check for "reasonable" input values

  if (b_ms > 3.5) then
     msg = 'b parameter of the mass-size relation cannot be greater than 3.5!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif (a_ms == 0.0) then
     msg = 'a parameter of mass-size relation cannot be equal to 0!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif (rho_ms == 0.0) then
     msg = 'Particle density rho cannot be equal to 0!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif (as_ratio == 0.0) then
     msg = 'Particle aspect ratio cannot be equal to 0!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif (liq_ice == -1 .and. rho_ms > 917.0) then
     msg = 'Density of ice cannot exceed 917.0 kg/m3!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif (liq_ice == 1 .and. (a_ms > 0.0 .or. b_ms > 0.0 .or. rho_ms > 0.0 .or. &
                             (as_ratio /= 1. .and. as_ratio > 0.) )) then
     msg = 'Water droplet: a, b, as_ratio replaced with pi/6 * 1000 kg/m3, 3, and 1!'
     err = warning
     call report(err, msg, nameOfRoutine)
  elseif (liq_ice == -1 .and. a_ms < 0.0 .and. b_ms < 0.0 .and. rho_ms < 0.0) then
     msg = 'Need to specify a & b and/or rho to calculate the mass-size relation parameters!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif (liq_ice == -1 .and. a_ms < 0.0 .and. b_ms < 0.0 .and. rho_ms < 0.0) then
     msg = 'Need to specify a & b or rho to calculate the mass-size relation parameters!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif ((a_ms < 0. .and. b_ms >= 0.0) .or. (a_ms > 0.0 .and. b_ms < 0.0)) then
     msg = 'When a is specified then also b need to be specified and vice versa'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  else
     err = success
  endif

! Liquid cloud droplet or rain --> use spherical particles with rho = 1000. kg/m3, b = 3 and a = pi/6 * rho
  if (liq_ice == 1) then
     b_ms = 3._dbl
     a_ms = pi / 6._dbl * rho_water
     errorstatus = err
     if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
     return
  endif

! a_ms and b_ms are specified --> nothing to do, exit the subroutine
  if (a_ms > 0.0 .and. b_ms >= 0 .and. liq_ice == -1) then
     errorstatus = err
     if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
     return
  endif

! Fixed density --> calculate a_ms, b_ms
! !  aspect ratio not defined or equal to 1 --> use spheres
  if (liq_ice == -1 .and. (as_ratio < 0.0 .or. as_ratio == 1.0)) then
     b_ms = 3._dbl
     a_ms = pi / 6._dbl * rho_ms
     errorstatus = err
     if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
     return
! ! as_ratio < 1.0 --> prolate spheroid
  elseif (liq_ice == -1 .and. as_ratio > 0.0 .and. as_ratio < 1.0) then
     b_ms = 3._dbl
     a_ms = pi / 6._dbl * rho_ms * (as_ratio) ** 2._dbl
     errorstatus = err
     if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
     return
! ! as_ratio > 1.0 --> oblate spheroid
  elseif (liq_ice == -1 .and. as_ratio > 1.0) then
     b_ms = 3._dbl
     a_ms = pi / 6._dbl * rho_ms / (as_ratio)
     errorstatus = err
     if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
     return
  endif

  msg = 'If you are here something wrong append...'
  errorstatus = fatal
  call report(errorstatus, msg, nameOfRoutine)
  return

end subroutine make_mass_size