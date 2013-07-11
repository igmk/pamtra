!+ drop-size distribution parameters calculation
subroutine make_dist_params(errorstatus)

! Description:
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
                    long       ! integer parameter specifying long integer

  use report_module

  use constants, only: pi, rho_water, delta_d_mono

  use drop_size_dist, only: dist_name, p_1, p_2, p_3, p_4, a_ms, b_ms, d_1, moment_in, & ! IN
                            q_h, n_tot, r_eff, t,                                      & ! IN
                            n_0, lambda, mu, gam, n_t, sig, d_ln, d_mono                 ! OUT

! Imported Scalar Variables with intent (in):

  implicit none

!- End of header ---------------------------------------------------------------

! Local scalars:

  real(kind=dbl) :: work1, work2, work3
  real(kind=dbl) :: ztc, hlp, alf, bet, m2s, m3s
  integer(kind=long) :: nn

! Local arrays:

  real(kind=dbl), dimension(10) :: mma, mmb

! Error handling

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=16) :: nameOfRoutine = 'make_dist_params'

! Used functions

  real(kind=dbl) :: nan

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

! Initialize output variables
  n_0 = -99.
  lambda = -99.
  gam = nan()
  mu = nan()
  n_t = -99.
  sig = -99.
  d_ln = nan()
  d_mono = -99.

! check for "reasonable" input values
  if (moment_in == 3 .and. q_h <= 0.) then      
    msg = 'if moment_in eq 3 then input moment q_h must be greater than 0!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  elseif (moment_in == 1 .and. n_tot <= 0.) then      
    msg = 'if moment_in eq 2 then input moment n_tot must be greater than 0!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  elseif (moment_in == 2 .and. r_eff <= 0.) then      
    msg = 'if moment_in eq 3 then input moment r_eff must be greater than 0!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  elseif (moment_in == 13 .and. (q_h <= 0. .or. n_tot <= 0.)) then      
    msg = 'if moment_in eq 12 then input moment q_h and n_tot must be greater than 0!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  elseif (moment_in == 23 .and. (q_h <= 0. .or. r_eff <= 0.)) then      
    msg = 'if moment_in eq 13 then input moment q_h and r_eff must be greater than 0!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  elseif (moment_in == 12 .and. (n_tot <= 0. .or. r_eff <= 0.)) then      
    msg = 'if moment_in eq 23 then input moment n_tot and r_eff must be greater than 0!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  else
    err = success
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MONODISPERSE distribution   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (trim(dist_name) == 'mono' .or. trim(dist_name) == 'mono_cosmo_ice') then 
! Set parameter for the Gamma dist. to get a monodisperse dist.
    lambda = 0._dbl
    mu = 0._dbl
    gam = 0._dbl
    if (trim(dist_name) == 'mono') then
! ! fixed radius (via d_1)
      if (d_1 /= -99. .and. p_1 == -99.) then
        d_mono = d_1
        if (moment_in == 3)    n_0 = q_h / (delta_d_mono * a_ms * d_1**b_ms)
        if (moment_in == 1)    n_0 = n_tot / delta_d_mono
      endif
! ! fixed n_tot (via p_1)
      if (p_1 /= -99. .and. d_1 == -99.) then
        n_0 = p_1 / delta_d_mono
        if (moment_in == 3)    d_mono = (q_h / (p_1 * a_ms))**(1._dbl / b_ms)
        if (moment_in == 2)    d_mono = r_eff / 2._dbl
      endif
    endif
    if (trim(dist_name) == 'mono_cosmo_ice' .and. moment_in == 3) then
! ! Monodisperse size distribution coherent with COSMO-de 1-moment scheme
! Number_concentration of activated ice ctystal is temperature dependent
! from COSMO-de code src_gscp.f90 routine: hydci_pp_gr
! Radius is derived from mass-size relation m=aD^3
! a=130 kg/m^3 (hexagonal plates with aspect ratio of 0.2 -> thickness=0.2*Diameter)
      work1 = 1.d2 * exp(0.2_dbl * (273.15_dbl - t)) ! N_tot
      n_0 = work1 / delta_d_mono
      d_mono = (q_h / (work1 * a_ms))**(1._dbl / b_ms)
!  CHECK if dia1 > maxdiam=2.d-4 (maximum diameter for COSMO)
!  then recalculate the drop mass using 2.d-4 as particle diameter
      if (d_mono > 2.d-4) then 
        d_mono = 2.d-4
        n_0 = q_h / (a_ms * d_mono**b_ms) / delta_d_mono
      endif
    endif
! ! Check that the variables have been filled in
    if (lambda /= 0._dbl .or. mu /= 0._dbl .or. gam /= 0._dbl .or. &
       n_0 <= 0._dbl .or. d_mono <= 0._dbl ) then
      msg = 'Monodisperse case: something wrong or this parameters combination is not yet implemented...'
      errorstatus = fatal
      call report(errorstatus, msg, nameOfRoutine)
      return
    endif
    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EXPONENTIAL distribution   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (trim(dist_name) == 'exp' .or. trim(dist_name) == 'exp_field_t' .or.       &
      trim(dist_name) == 'exp_cosmo_snow' .or. trim(dist_name) == 'exp_ryan')     then 
! Set parameter for the Gamma dist. to get an exponential dist.
    gam = 1._dbl
    mu = 0._dbl
    if (trim(dist_name) == 'exp') then
! ! fixed n_tot (via p_1)
      if (p_1 /= -99. .and. p_2 == -99. .and. p_3 == -99.) then
        if (moment_in == 3)  lambda = (p_1 * a_ms * dgamma(b_ms+1._dbl) / q_h)**(1._dbl / b_ms)
        if (moment_in == 2)  lambda = 3._dbl / r_eff
        n_0    = p_1 * lambda
      endif
! ! fixed r_eff (via p_2)
      if (p_1 == -99. .and. p_2 /= -99. .and. p_3 == -99.) then
        lambda = 3._dbl / p_2
        if (moment_in == 3)  n_0 = (q_h * lambda**(b_ms+1._dbl)) / (a_ms * dgamma(b_ms+1._dbl))
        if (moment_in == 1)  n_0 = n_tot * lambda
      endif
! ! fixed N_0 (via p_3)
      if (p_1 == -99. .and. p_2 == -99. .and. p_3 /= -99.) then
        n_0 = p_3
        if (moment_in == 3)  lambda = (n_0 * a_ms * dgamma(b_ms+1._dbl) / q_h)**(1._dbl / (b_ms+1._dbl))
        if (moment_in == 1)  lambda = n_0 / n_tot
        if (moment_in == 2)  lambda = 3._dbl / r_eff
      endif
! ! both moments from input file
      if (p_1 == -99. .and. p_2 == -99.) then
        if (moment_in == 13)  then
          lambda = (n_tot * a_ms * dgamma(b_ms+1._dbl) / q_h)**(1._dbl / b_ms)
          n_0    = n_tot * lambda
        endif
        if (moment_in == 23)  then
          lambda = 3._dbl / r_eff
          n_0    = (q_h * lambda**(b_ms+1)) / (a_ms * dgamma(b_ms+1._dbl))
        endif
        if (moment_in == 12)  then
          lambda = 3._dbl / r_eff
          n_0    = n_tot * lambda
        endif
      endif
    endif
! ! Ryan (2000 JAS) Lambda = Lambda(t)
    if (trim(dist_name) == 'exp_ryan') then
      lambda = 1220._dbl * 10._dbl**(-0.0245_dbl * (273.15_dbl-t))
      if (moment_in == 3) n_0 = (q_h * lambda**(b_ms+1)) / (a_ms * dgamma(b_ms+1._dbl))
      if (moment_in == 1) n_0 = n_tot * lambda
    endif
! ! Field et al. (2005 QJRM, end of page 2008) n_0 = n_0(T)
    if (trim(dist_name) == 'exp_field_t') then
      n_0 = 5.65d5 * exp(0.107_dbl * (273.15_dbl - t))
      if (moment_in == 3) lambda = (a_ms * n_0 * dgamma(b_ms+1._dbl) / q_h)**(1._dbl /(b_ms+1._dbl))
      if (moment_in == 2) lambda = 3._dbl / r_eff
    endif
! ! Field et al. (2005 QJRM, page 2007) n_0 = n_0(T,q_h)
! ! AS IMPLEMENTED IN cosmo 2-MOMENT scheme
    if (trim(dist_name) == 'exp_cosmo_snow' .and. moment_in == 3) then
! taken from COSMO-de routine hydci_pp_gr in src_gscp.f90
      mma = (/   5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
             0.312550,  0.000204,  0.003199, 0.000000, -0.015952 /)
      mmb = (/   0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
             0.060366,  0.000079,  0.000594, 0.000000, -0.003577 /)
! Calculate n0s using the temperature-dependent moment
      ztc = t - 273.15                           ! temperature in C
      ztc = MAX(MIN(ztc,0.0),-40.0)              !limited to -40
! the moments of every order of the particle size distribution f(D) depend from the second one (proportional to mass) and temperature
      nn  = 3
      hlp = mma(1)      +mma(2)*ztc      +mma(3)*nn       +mma(4)*ztc*nn+mma(5)*ztc**2 &
          + mma(6)*nn**2+mma(7)*ztc**2*nn+mma(8)*ztc*nn**2+mma(9)*ztc**3+mma(10)*nn**3
      alf = 10.0d0**hlp
      bet = mmb(1)      +mmb(2)*ztc      +mmb(3)*nn       +mmb(4)*ztc*nn+mmb(5)*ztc**2 &
          + mmb(6)*nn**2+mmb(7)*ztc**2*nn+mmb(8)*ztc*nn**2+mmb(9)*ztc**3+mmb(10)*nn**3
      m2s = q_h / 0.038d0 ! 0.038 = Formfactor in the mass-size relation of snow particles [kg/m^2]
      m3s = alf*EXP(bet*LOG(m2s))
      hlp  =  5.65d5 * EXP(-0.107d0*ztc) ! N0snow as a function of solely T
      n_0 = 13.5 * m2s**4 / m3s**3       ! N0snow as a function of T and snow mixing ratio
      n_0 = MAX(n_0,0.5*hlp)
      n_0 = MIN(n_0,1e2*hlp)
      n_0 = MIN(n_0,1e9)
      n_0 = MAX(n_0,1e6)
      lambda = (a_ms * n_0 * dgamma(b_ms+1._dbl) / q_h)**(1._dbl /(b_ms+1._dbl))
    endif
! ! Check that the variables have been filled in
    if (gam /= 1._dbl .or. mu /= 0._dbl .or. n_0 <= 0._dbl .or. lambda <= 0._dbl) then
      print*, gam, mu, n_0, lambda
      msg = 'Exponential case: something wrong or this parameters combination is not yet implemented...'
      errorstatus = fatal
      call report(errorstatus, msg, nameOfRoutine)
      return
    endif
    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! LOG-NORMAL distribution   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (trim(dist_name) == 'logn') then
! ! fixed n_tot (via p_1) and sigma (via p_3)
    if (p_1 /= -99. .and. p_3 /= -99.) then
      n_t = p_1
      sig = p_3
      if (moment_in == 3) d_ln = (1._dbl/b_ms) * (log(q_h / (a_ms * n_t)) - (0.5_dbl * (b_ms * sig)**2._dbl))
      if (moment_in == 2) d_ln = log(2._dbl * r_eff) - 2.5_dbl * sig**2._dbl
    endif
! ! fixed n_tot (via p_1) and r_eff (via p_2)
    if (p_1 /= -99. .and. p_2 /= -99. .and. p_3 == -99. .and. moment_in == 3) then
      n_t = p_1
      work1 = q_h / (a_ms * n_t * (2._dbl * p_2)**2._dbl)
      work2 = 2._dbl / (b_ms*(b_ms-5._dbl))
      sig = sqrt(log(work1) * work2)
      d_ln = log(2._dbl * p_2) - 2.5_dbl * sig**2._dbl
    endif
! ! fixed r_eff (via p_2) and sigma (via p_3)
    if (p_1 == -99. .and. p_2 /= -99. .and. p_3 /= -99.) then
      sig = p_3
      d_ln = log(2._dbl * p_2) - 2.5_dbl * sig**2._dbl
      if (moment_in == 3) n_t = q_h / a_ms * exp(-1._dbl*((b_ms * d_ln) + (0.5_dbl *(b_ms * sig)**2._dbl)))
      if (moment_in == 1) n_t = n_tot
    endif
! ! fixed sigma (via p_3)
    if (p_1 == -99. .and. p_2 == -99. .and. p_3 /= -99.) then
      sig = p_3
      if (moment_in == 13) then
        n_t  = n_tot
        d_ln = (1._dbl/b_ms) * (log(q_h / (a_ms * n_t)) - (0.5_dbl * (b_ms * sig)**2._dbl))
      endif
      if (moment_in == 23) then
        d_ln = log(2._dbl * r_eff) - 2.5_dbl * sig**2._dbl
        n_t  = q_h / a_ms * exp(-1._dbl*((b_ms * d_ln) + (0.5_dbl *(b_ms * sig)**2._dbl)))
      endif
      if (moment_in == 12) then
        n_t  = n_tot
        d_ln = log(2._dbl * r_eff) - 2.5_dbl * sig**2._dbl
      endif
    endif
! ! Check that the variables have been filled in
    if (n_t <= 0._dbl .or. isnan(d_ln) .or. sig <= 0._dbl) then
      msg = 'Log-normal case: something wrong or this parameters combination is not yet implemented...'
      errorstatus = fatal
      call report(errorstatus, msg, nameOfRoutine)
      return
    endif
    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODIFIED GAMMA distribution   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (trim(dist_name) == 'mgamma') then
! ! The user MUST specify mu and gam parameters
    if (p_3 == -99. .or. p_4 == -99.) then 
      msg = 'Modified Gamma case: p_3 and p_4 parameters must be specified...'
      errorstatus = fatal
      call report(errorstatus, msg, nameOfRoutine)
      return
    endif
    mu  = p_3
    gam = p_4
! ! fixed n_tot (via p_1)
    if (p_1 /= -99. .and. p_2 == -99.) then
      if (moment_in == 3)  then
        work1 = a_ms * gam**(b_ms-1._dbl) * p_1**(b_ms) / q_h
        work2 = dgamma((mu + b_ms + 1._dbl) / gam)
        work3 = dgamma((mu + 1._dbl) / gam)
        n_0 = (work1 * work2 / work3**(b_ms))**(1._dbl / (b_ms - 1._dbl))
        lambda = (n_0 * work3 / (gam * p_1))**(gam / (mu + 1._dbl))
      endif
      if (moment_in == 2)  then
        work1 = dgamma((mu + 4._dbl) / gam)
        work2 = dgamma((mu + 3._dbl) / gam) * r_eff * 2._dbl
        lambda = (work1 / work2)**gam
        work1 = (mu + 1._dbl) / gam
        n_0 = p_1 * gam * lambda**work1 / dgamma(work1)
      endif
    endif
! ! fixed r_eff (via p_2)
    if (p_1 == -99. .and. p_2 /= -99.) then
      if (moment_in == 3)  then
        work1 = dgamma((mu + 4._dbl) / gam)
        work2 = dgamma((mu + 3._dbl) / gam) * p_2 * 2._dbl
        lambda = (work1 / work2)**gam
        work1 = q_h * gam / a_ms
        work2 = (mu + b_ms + 1._dbl) / gam
        n_0 = work1 * lambda**work2 / dgamma(work2)
      endif
      if (moment_in == 1)  then
        work1 = dgamma((mu + 4._dbl) / gam)
        work2 = dgamma((mu + 3._dbl) / gam) * p_2 * 2._dbl
        lambda = (work1 / work2)**gam
        work1 = (mu + 1._dbl) / gam
        n_0 = n_tot * gam * lambda**work1 / dgamma(work1)
      endif
    endif
! ! 2 moments from the input file
    if (p_1 == -99. .and. p_2 == -99.) then
      if (moment_in == 13)  then
        work1 = a_ms * gam**(b_ms-1._dbl) * n_tot**(b_ms) / q_h
        work2 = dgamma((mu + b_ms + 1._dbl) / gam)
        work3 = dgamma((mu + 1._dbl) / gam)
        n_0 = (work1 * work2 / work3**(b_ms))**(1._dbl / (b_ms - 1._dbl))
        lambda = (n_0 * work3 / (gam * n_tot))**(gam / (mu + 1._dbl))
      endif
      if (moment_in == 23)  then
        work1 = dgamma((mu + 4._dbl) / gam)
        work2 = dgamma((mu + 3._dbl) / gam) * r_eff * 2._dbl
        lambda = (work1 / work2)**gam
        work1 = q_h * gam / a_ms
        work2 = (mu + b_ms + 1._dbl) / gam
        n_0 = work1 * lambda**work2 / dgamma(work2)
      endif
      if (moment_in == 12)  then
        work1 = dgamma((mu + 4._dbl) / gam)
        work2 = dgamma((mu + 3._dbl) / gam) * r_eff * 2._dbl
        lambda = (work1 / work2)**gam
        work1 = (mu + 1._dbl) / gam
        n_0 = n_tot * gam * lambda**work1 / dgamma(work1)
      endif
    endif
! ! Check that the variables have been filled in
    if (lambda <= 0._dbl .or. n_0 <= 0._dbl .or. isnan(mu) .or. isnan(gam)) then
      msg = 'Modified-gamma case: something wrong or this parameters combination is not yet implemented...'
      errorstatus = fatal
      call report(errorstatus, msg, nameOfRoutine)
      return
    endif
    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return
  endif

!  If we are here something went wrong...
  msg = 'Distribution not yet implemented: '//trim(dist_name)
  errorstatus = fatal
  call report(errorstatus, msg, nameOfRoutine)
  return
  
end subroutine make_dist_params