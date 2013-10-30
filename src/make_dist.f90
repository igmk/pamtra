!+ drop-size distribution calculation
subroutine make_dist(errorstatus)

! Description:
! Given the distribution parameters, calculate the particle number 
! concentration [#/m3] according to the following distributions:
! 
! LOG-NORMAL:      n(D) = n_t / (D * sig * sqrt(2. * pi)) * EXP[-( (ln(D) - d_ln)**2 )/(2. * sig**2) ]
! 
! Normalized gamma:       
!       X =  D(i)/d_m
!       n(D) = n_0_star * (gamma(b_ms+1)/(b_ms+1)**(b_ms+1) * (b_ms+mu+1)**(b_ms+mu+1)/gamma(b_ms+mu+1)) * X**mu * exp(-(b_ms+mu+1)*X)
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
                            d_ds, d_bound_ds, f_ds, n_ds, delta_d_ds, b_ms, n_0_star, d_m
  implicit none

!- End of header ---------------------------------------------------------------

! Local character:

  character(len=4) :: str_work

! Local scalar:

  real(kind=dbl) :: d_1_work, d_2_work, work1, tmp1, tmp2, tmpX
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
! normalized modified gamma
!follows the definition from Testud et al, 
!but with Dmax as size descriptor and 
!WITHOUT fixed exponent of mass size relation! Instead, use, b of mass size relation
!d_m defined as M_(b+1)/M_b (with M_i the ith moment)
!N_0_star = (b_ms+1)**(b_ms+1)/gamma(b_ms+1) *M_b**(b_ms+2)/M_(b+1)**(b_ms+1)
  else if (trim(dist_name) == 'norm_mgamma') then
    do i=1,nbin+1
      tmpX =  d_bound_ds(i)/d_m
      tmp1 = gamma(b_ms+1)/(b_ms+1)**(b_ms+1) * (b_ms+mu+1)**(b_ms+mu+1)/gamma(b_ms+mu+1)
      tmp2 = exp(-(b_ms+mu+1)*tmpX)
      f_ds(i) = n_0_star * tmp1 * tmpX**mu * tmp2
    enddo


! modified gamma, gamma or exponetial distribution
  else
    do i=1,nbin+1
      f_ds(i) = n_0 * d_bound_ds(i)**mu * EXP(-lambda * d_bound_ds(i)**gam)
    enddo
  endif

  do i=1,nbin
    delta_d_ds(i) =  d_bound_ds(i+1) - d_bound_ds(i)
    n_ds(i) = (f_ds(i) + f_ds(i+1)) / 2._dbl * delta_d_ds(i)  ! trapezoid approximation of the integral
  enddo

! print*,'d_ds',d_ds
! print*,'f_ds',f_ds
! print*,'n_ds',n_ds

  errorstatus = success

  if (verbose >= 1) call report(info,'End of ', nameOfRoutine)

  return
end subroutine make_dist