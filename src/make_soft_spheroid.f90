!+ soft spheroid diameter and density calculation
subroutine make_soft_spheroid(errorstatus)

! Description:
!  The routine also calculate diameter and density (soft_d_eff, 
!  soft_rho_eff) of soft spheres for each of the size bins.
!  The relations used are:
!  m_i = a * D^b 
!  m_i = rho_i * volume_i
!  volume_i = pi/6 * rho_i * soft_d_eff^3 * (As_ratio)^2  for prolate (as_ratio < 1.)
!  volume_i = pi/6 * rho_i * soft_d_eff^3 / As_ratio      for oblate  (as_ratio > 1.)
!  From these relation we get soft_rho_eff, if rho_ms is set to -99, or 
!  soft_d_eff otherwise.
! 
! Owner: IGMK
!
! History:
!
! Version Date         Comment
! ------- ----         -------
! 0.      21/06/2013   Creation - E. Orlandi
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

  use constants, only: pi, rho_ice

  use drop_size_dist, only: rho_ms, as_ratio, a_ms, b_ms, d_ds,nbin, &    ! IN
                            soft_rho_eff, soft_d_eff                      ! OUT

  implicit none

!- End of header ---------------------------------------------------------------

! Local array:

  real(kind=dbl), dimension(nbin) :: mass         ! particle mass  [kg]

! Local scalar:

  integer(kind=long) :: i

! Error handling

  integer(kind=long), intent(out) :: errorstatus
!   integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'make_soft_spheroid'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  allocate(soft_rho_eff(nbin))
  allocate(soft_d_eff(nbin))

! Calculate particle mass
  do i=1,nbin
    mass(i) = a_ms * d_ds(i)**b_ms
  enddo

! density set to a fixed value by the user
! Calculate the diameter of the soft spheroids
  if (rho_ms > 0.) then
    soft_rho_eff(:) = rho_ms
! oblate spheroid or sphere
    if (as_ratio < 0. .or. as_ratio >= 1.) then
      do i=1,nbin
        if (as_ratio < 0.)  soft_d_eff(i) = ((6._dbl * mass(i))            / (pi *  rho_ms))**(1._dbl/3._dbl)
        if (as_ratio >= 1.) soft_d_eff(i) = ((6._dbl * mass(i) * as_ratio) / (pi *  rho_ms))**(1._dbl/3._dbl)
      enddo
    endif
! prolate spheroid or sphere
    if (as_ratio > 0. .and. as_ratio < 1.) then
      do i=1,nbin
        soft_d_eff(i) = ((6._dbl * mass(i)) / (pi * rho_ms * as_ratio**2._dbl))**(1._dbl/3._dbl)
      enddo
    endif
  endif

! Calculate the density of the soft spheroids
  if (rho_ms < 0.) then
    soft_d_eff = d_ds
    ! oblate spheroid or sphere
    if (as_ratio < 0. .or. as_ratio >= 1.) then
      do i=1,nbin
        if (as_ratio < 0.)  soft_rho_eff(i) = (6._dbl * mass(i))            / (pi *  d_ds(i)**3._dbl)
        if (as_ratio >= 1.) soft_rho_eff(i) = (6._dbl * mass(i) * as_ratio) / (pi *  d_ds(i)**3._dbl)
        if (soft_rho_eff(i) < 5._dbl) soft_rho_eff(i) = 5._dbl
        if (soft_rho_eff(i) > rho_ice) soft_rho_eff(i) = rho_ice
      enddo
    endif
! prolate spheroid or sphere
    if (as_ratio > 0. .and. as_ratio < 1.) then
      do i=1,nbin
        soft_rho_eff(i) = (6._dbl * mass(i)) / (pi *  d_ds(i)**3._dbl * as_ratio**2._dbl) 
        if (soft_rho_eff(i) < 5._dbl) soft_rho_eff(i) = 5._dbl
        if (soft_rho_eff(i) > rho_ice) soft_rho_eff(i) = rho_ice
      enddo
    endif
  endif

  if (minval(soft_rho_eff) <= 0. .or. minval(soft_d_eff) <= 0.) then
     msg = 'something wrong in make_soft_spheroid!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  endif

  errorstatus = success
  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

  if (verbose >= 4) then
    print*,'soft spheroids diameter [m]'
    print*,soft_d_eff
    print*,'soft spheroids density [kg/m3]'
    print*,soft_rho_eff
  endif

  return

end subroutine make_soft_spheroid