!+ drop-size distribution moments calculation
subroutine calc_moment(errorstatus)

! Description:
!
! Owner: IGMK
!
! History:
!
! Version Date         Comment
! ------- ----         -------
! 0.      28/05/2013   Creation - E. Orlandi
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
                   long        ! integer parameter specifying long integer

  use report_module

  use drop_size_dist, only: a_ms, b_ms, nbin, d_ds, d_bound_ds, n_ds, f_ds, &   ! IN
                             m_0, m_32, am_b                                    ! OUT

! Imported Scalar Variables with intent (in):

  implicit none

!- End of header ---------------------------------------------------------------

! Subroutine arguments

! Local scalars:

  real(kind=dbl) :: work2, work3
  integer(kind=long) :: i

! Error handling

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=16) :: nameOfRoutine = 'calc_moment'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

! check for "reasonable" input values
  if (maxval(d_ds) <= 0. .or. maxval(d_bound_ds) <= 0.) then
    msg = 'maximum particle diameter <= 0!!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  elseif (maxval(n_ds) <= 0.) then
    msg = 'maximum particle number concentration <= 0!!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  else
     err = success
  endif

! m_0 --> Total number concentration
  m_0 = 0._dbl
  do i=1,nbin
    m_0 = m_0 + (f_ds(i) + f_ds(i+1)) /2._dbl * (d_bound_ds(i+1) - d_bound_ds(i))
  enddo
  if (sum(n_ds) - m_0 > 1.d-12) then
    msg = 'Something wrong with total number concentration!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  endif

! m_3 / m_2 --> Effective radius
  work2 = 0._dbl
  do i=1,nbin
    work2 = work2 + (f_ds(i) + f_ds(i+1)) / 2._dbl * (d_bound_ds(i+1) - d_bound_ds(i)) * d_ds(i)**2
  enddo
  work3 = 0._dbl
  do i=1,nbin
    work3 = work3 + (f_ds(i) + f_ds(i+1)) / 2._dbl * (d_bound_ds(i+1) - d_bound_ds(i)) * d_ds(i)**3
  enddo
  m_32 = work3 / work2

! am_b  --> Total mass concentration
  am_b = 0._dbl
  do i=1,nbin
    am_b = am_b + a_ms * (f_ds(i) + f_ds(i+1)) / 2._dbl * (d_bound_ds(i+1) - d_bound_ds(i)) * d_ds(i)**b_ms
  enddo

  errorstatus = err
  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

  return
end subroutine calc_moment