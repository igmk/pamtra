!+ drop-size distribution calculation
subroutine make_dist(errorstatus)

! Description:
! Given the distribution parameters, calculate the particle number 
! concentration [#/m3] according to the following distributions:
! 
! LOG-NORMAL:      n(D) = n_t / (D * sig * sqrt(2. * pi)) * EXP[-( (ln(D) - d_ln)**2 )/(2. * sig**2) ]
! 
! MODIFIED-GAMMA:  n(D) = n_0 * D**mu * EXP(-lambda * D**gam)
! 
! The following distribution are special cases of the mod-gamma distribution
! MONODISPERSE:    n(D) = n_0                    --> mu=0.; gam=0.; lambda=0.
! EXPONENTIAL:     n(D) = n_0 * EXP[-lambda * D] --> mu=0.; gam=0.
! 
!
! Owner: IGMK
!
! History:
!
! Version Date         Comment
! ------- ----         -------
! 0.      23/05/2013   Creation - E. Orlandi
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

  use constants, only: pi, delta_d_mono

  use drop_size_dist, only: dist_name, d_1, d_2, nbin, n_0, lambda, mu, gam, n_t, d_ln, sig, d_mono, & ! IN
                            d_ds, d_bound_ds, f_ds, n_ds, delta_d_ds                                               ! OUT

  implicit none

!- End of header ---------------------------------------------------------------

! Local character:

  character(len=4) :: str_work

! Local scalar:

  real(kind=dbl) :: d_1_work, d_2_work, work1
  integer(kind=long) :: i

! Error handling

  integer(kind=long), intent(out) :: errorstatus
!   integer(kind=long) :: err = 0
!   character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'make_dist'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  d_1_work = d_1
  d_2_work = d_2

! If monodisperse distribution then set D_1, D_2
  str_work = trim(dist_name)
  if (str_work == 'mono') then
    d_1_work = d_mono - delta_d_mono*.5_dbl
    d_2_work = d_mono + delta_d_mono*.5_dbl
  endif

! Create particle diameter array
  work1 = (d_2_work - d_1_work) / nbin      ! Delta diameter
  do i=1,nbin+1 
    d_bound_ds(i) = d_1_work + work1 * (i-1)
    if (i <= nbin) d_ds(i) = d_bound_ds(i) + work1 * .5_dbl
  enddo

! calculate the particle number concentration
! ! LOG-NORMAL distribution
  if (trim(dist_name) == 'logn') then
    do i=1,nbin+1
      f_ds(i) = n_t / (d_bound_ds(i) * sig * sqrt(2._dbl * pi)) * &
                EXP(-( (log(d_bound_ds(i)) - d_ln)**2 )/(2. * sig**2) )
    enddo
  else
    do i=1,nbin+1
      f_ds(i) = n_0 * d_bound_ds(i)**mu * EXP(-lambda * d_bound_ds(i)**gam)
    enddo
  endif
  do i=1,nbin
    n_ds(i) = (f_ds(i) + f_ds(i+1)) / 2._dbl * (d_bound_ds(i+1)-d_bound_ds(i))  ! trapezoid approximation of the integral
    delta_d_ds(i) =  d_bound_ds(i+1) - d_bound_ds(i)
  enddo

! print*,'d_ds',d_ds
! print*,'f_ds',f_ds
! print*,'n_ds',n_ds

  errorstatus = success

  if (verbose >= 1) call report(info,'End of ', nameOfRoutine)

  return
end subroutine make_dist