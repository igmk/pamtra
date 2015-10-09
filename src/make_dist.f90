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
  
  use settings, only: hydro_adaptive_grid

  use drop_size_dist, only: dist_name, d_1, d_2, nbin, n_0, lambda, mu, gam, n_t, d_ln, sig, d_mono, & ! IN
                            d_ds, d_bound_ds, f_ds, n_ds, delta_d_ds, b_ms, n_0_star, d_m, a_ms
                            
  use descriptor_file
  implicit none

!- End of header ---------------------------------------------------------------

! Local character:

  character(len=80) :: msg

! Local scalar:

  real(kind=dbl) :: d_1_work, d_2_work, work1, tmp1, tmp2, tmpX, n_tot, am_b
  real(kind=dbl) :: d_1_new, d_2_new, min_lin, min_log, thres_n
  integer(kind=long) :: i, ibig, nbin_work, nbin_log, nbin_lin
  integer(kind=long) :: erroalloc

! Error handling

  integer(kind=long), intent(out) :: errorstatus
!   integer(kind=long) :: err = 0
!   character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'make_dist'
  logical :: skip
  
  real(kind=dbl), allocatable, dimension(:) :: f_ds_work, d_bound_ds_work, d_bound_ds_lin, d_bound_ds_log
  
  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  skip = .false.
  thres_n = .1
    
  ! Monodisperse distribution
  if ((trim(dist_name) == 'mono')  .or. (trim(dist_name) == 'mono_cosmo_ice')) then
    d_bound_ds(1) = d_mono - delta_d_mono*.5_dbl
    d_bound_ds(2) = d_mono
    d_bound_ds(3) = d_mono + delta_d_mono*.5_dbl
    do i=1,nbin+1
      f_ds(i) = n_0 * d_bound_ds(i)**mu * EXP(-lambda * d_bound_ds(i)**gam)
    enddo
    skip = .true.
  endif

! Fisrt loop:
! -uses a very large number of bin and wide range of diameters to found where f_ds_work(D) is smaller than thres_n
! Second loop:
! -uses the user defined nbin BUT the nre defined diameter range
  bigloop: do ibig=1,2 ! Loop 2 is done only if hydro_adaptive_grid = .true.
    if (skip) exit bigloop
    if ((ibig == 1) .and. (.not.hydro_adaptive_grid)) then 
      d_1_work = d_1
      d_2_work = d_2
      nbin_work = nbin
      skip = .true.
    endif
    if ((ibig == 1) .and. hydro_adaptive_grid) then 
      d_1_work = 1.d-8
      d_2_work = 2.d-2
      nbin_work = 5.d2
    endif
    if (ibig == 2) then 
      d_1_work = d_1_new
      d_2_work = d_2_new
      nbin_work = nbin
      skip = .true.
    endif
    
    if (allocated(f_ds_work)) deallocate(f_ds_work,STAT=erroalloc)
    if (allocated(d_bound_ds_work)) deallocate(d_bound_ds_work,STAT=erroalloc)
    allocate(f_ds_work(nbin_work+1),STAT=erroalloc)
    allocate(d_bound_ds_work(nbin_work+1),STAT=erroalloc)
    nbin_lin = int(nbin_work/2._dbl)
    nbin_log = int(nbin_work/2._dbl)
    if (nbin_lin + nbin_log == nbin_work) nbin_lin = nbin_lin - 1
    if (allocated(d_bound_ds_lin)) deallocate(d_bound_ds_lin,STAT=erroalloc)
    if (allocated(d_bound_ds_log)) deallocate(d_bound_ds_log,STAT=erroalloc)
    allocate(d_bound_ds_lin(nbin_lin+1),STAT=erroalloc)
    allocate(d_bound_ds_log(nbin_log+1),STAT=erroalloc)

!   ! Create particle diameter array
    work1 = (d_2_work - d_1_work) / (nbin_lin + 1)      ! Delta diameter, +1 is to avoid having twice the first and last values
    do i=1,nbin_lin+1 
      d_bound_ds_lin(i) = d_1_work + work1 * (i)        ! this have the largest diameter defined (d_2_work)
    enddo
! !   Create a logspace diameter array
    work1 = (dlog(d_2_work)-dlog(d_1_work)) / (nbin_log + 1)   ! Delta diameter, +1 is to avoid having twice the first and last values
    do i = 1, nbin_log+1
      d_bound_ds_log(i) = exp(work1 * (i-1)  + dlog(d_1_work)) ! this have the smallest diameter defined (d_1_work)
    enddo

! Put linear and log bin together and sort them
    do i=1,nbin_work+1
      min_lin = minval(d_bound_ds_lin)
      min_log = minval(d_bound_ds_log)
      if (min_log <= min_lin) then
        d_bound_ds_work(i) = min_log
        d_bound_ds_log(minloc(d_bound_ds_log)) = 1.e30
      endif
      if (min_log > min_lin) then
        d_bound_ds_work(i) = min_lin
        d_bound_ds_lin(minloc(d_bound_ds_lin)) = 1.e30
      endif
    enddo

  ! calculate the particle number concentration
  ! ! LOG-NORMAL distribution
    if (trim(dist_name) == 'logn') then
      do i=1,nbin_work+1
	f_ds_work(i) = n_t / (d_bound_ds_work(i) * sig * sqrt(2._dbl * pi)) * &
		  EXP(-( (log(d_bound_ds_work(i)) - d_ln)**2 )/(2. * sig**2) )
      enddo
  ! normalized modified gamma
  !follows the definition from Testud et al, 
  !but with Dmax as size descriptor and 
  !WITHOUT fixed exponent of mass size relation! Instead, use, b of mass size relation
  !d_m defined as M_(b+1)/M_b (with M_i the ith moment)
  !N_0_star = (b_ms+1)**(b_ms+1)/gamma(b_ms+1) *M_b**(b_ms+2)/M_(b+1)**(b_ms+1)
    else if (trim(dist_name) == 'norm_mgamma') then
      do i=1,nbin_work+1
	tmpX =  d_bound_ds_work(i)/d_m
	tmp1 = gamma(b_ms+1)/(b_ms+1)**(b_ms+1) * (b_ms+mu+1)**(b_ms+mu+1)/gamma(b_ms+mu+1)
	tmp2 = exp(-(b_ms+mu+1)*tmpX)
	f_ds_work(i) = n_0_star * tmp1 * tmpX**mu * tmp2
      enddo

    else if ((trim(dist_name) == 'mgamma')     .or. (trim(dist_name) == 'exp') .or. &
             (trim(dist_name) =='exp_field_t') .or. (trim(dist_name) == 'exp_cosmo_snow') .or. &
             (trim(dist_name) == 'exp_ryan') .or. (trim(dist_name) == 'mgamma_MNH') ) then
      do i=1,nbin_work+1
	f_ds_work(i) = n_0 * d_bound_ds_work(i)**mu * EXP(-lambda * d_bound_ds_work(i)**gam)
      enddo
      
    else if ((trim(dist_name) /= 'mono')  .and. (trim(dist_name) /= 'mono_cosmo_ice') .and. &
             (trim(dist_name) /= 'const') .and. (trim(dist_name) /= 'const_cosmo_ice')) then 
      msg = 'did not undestand drop size name'
      errorstatus = fatal
      call report(errorstatus, msg, nameOfRoutine)
      return
    endif

!  Find the d_1 and d_2 where f_ds_work(d) = thres_n
! STEP INTO this cycle ONLY if first loop (ibig == 1) and ONLY if adaptive grid
   if (hydro_adaptive_grid .and. (ibig == 1)) then
     d_1_new = -1._dbl
     d_2_new = -1._dbl
     
     search_loop1: do i=nbin_work,1,-1
       d_2_new = d_bound_ds_work(i+1)         ! smallest possible d_2 = d_bound_ds_work(2)
       if (f_ds_work(nbin_work+1) >= thres_n) then !right tail gt .1
         d_2_new = d_bound_ds_work(nbin_work+1)
         exit search_loop1
       endif
       if ((f_ds_work(i) >= thres_n) .and. (f_ds_work(i+1) < thres_n)) then
         exit search_loop1
       endif
     enddo search_loop1
     search_loop2: do i=1,nbin_work
       d_1_new = d_bound_ds_work(i)           ! largest possible d_1 =  d_bound_ds_work(nbin_work)
       if (((f_ds_work(i) < thres_n) .and. (f_ds_work(i+1) >= thres_n)) .or. (f_ds_work(i) > f_ds_work(i+1))) then
         exit search_loop2
       endif
     enddo search_loop2

     if ((d_1_new < -1.d-30) .or. (d_2_new < -1.d-30) .or. (d_1_new > d_2_new)) then
      msg = 'something went wrong with the adaptive grid'
      errorstatus = fatal
      print*, d_1_new, d_2_new
      call report(errorstatus, msg, nameOfRoutine)
      return
     endif
   endif

!Fill f_ds and d_bound_ds only at the last cycle
  if (skip) then 
    f_ds(:) = f_ds_work(:)
    d_bound_ds(:) = d_bound_ds_work(:)
  endif

enddo bigloop
    
  do i = 1, nbin
      d_ds(i) = (d_bound_ds(i) + d_bound_ds(i+1)) * .5_dbl
  enddo
    
  do i=1,nbin
    delta_d_ds(i) =  d_bound_ds(i+1) - d_bound_ds(i)
    n_ds(i) = (f_ds(i) + f_ds(i+1)) / 2._dbl * delta_d_ds(i)  ! Trapezoidal approximation of the integral
  enddo

  n_tot = SUM(n_ds)  
  
  do i=1,nbin
    am_b = am_b + a_ms * (f_ds(i) + f_ds(i+1)) / 2._dbl * (d_bound_ds(i+1) - d_bound_ds(i)) * d_ds(i)**b_ms
  enddo
  
  
  !remove numerical instabilities
  WHERE (n_ds < n_tot/nbin * 1d-60) n_ds = 0.d0

! print*, "lambda", lambda, "mu", mu, "n_0", n_0, "gam", gam
! print*, n_ds
! print*,'d_ds',d_ds
! print*,'d_bound_ds',d_bound_ds
! print*,'delta_d_ds',delta_d_ds
! print*,'f_ds',f_ds
! print*,'n_ds',n_ds
!    do i=1,nbin+1
!    print*,i,f_ds(i),d_bound_ds(i)!,n_ds(i)
!    enddo

  errorstatus = success

  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

  return
end subroutine make_dist