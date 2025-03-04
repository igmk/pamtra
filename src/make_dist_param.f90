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

  use nan, only: nan_dbl ! function returning not-a-number

  use report_module

  use constants, only: pi, rho_water, delta_d_mono

  use drop_size_dist, only: dist_name, rho_ms, p_1, p_2, p_3, p_4, a_ms, b_ms, d_1, d_2, moment_in, & ! IN
       q_h, n_tot, r_eff, layer_t, nbin, liq_ice,                            & ! IN
       n_0, lambda, mu, gam, n_t, sig, d_ln, d_mono, d_m, n_0_star ! OUT

  ! variables used for surface and pressure dependent parameters (ECHAM/HIRHAM)
  use vars_index, only: i_x, i_y, i_z
  use vars_atmosphere, only: sfc_type, atmo_press, atmo_temp

  ! Imported Scalar Variables with intent (in):

  implicit none

  !- End of header ---------------------------------------------------------------

  ! Local scalars:

  real(kind=dbl) :: work1, work2, work3, delta_d_const, rho_dry, kappa
  real(kind=dbl) :: ztc, hlp, alf, bet, m2s, m3s
  integer(kind=long) :: nn

  ! Local arrays:

  real(kind=dbl), dimension(10) :: mma, mmb

  ! functions

  real(kind=dbl) :: rho_air

  ! Error handling

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=16) :: nameOfRoutine = 'make_dist_params'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  ! Initialize output variables
  n_0 = -99.
  lambda = -99.
  gam = nan_dbl()
  mu = nan_dbl()
  n_t = -99.
  sig = -99.
  d_ln = nan_dbl()
  d_mono = -99.
  d_m = -99.
  n_0_star = -99.

  ! check for "reasonable" input values
  if ((moment_in == 3) .and. (q_h <= 0.)) then      
     msg = 'if moment_in eq 3 then input moment q_h must be greater than 0!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif ((moment_in == 1) .and. (n_tot <= 0.)) then      
     msg = 'if moment_in eq 1 then input moment n_tot must be greater than 0!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif ((moment_in == 2) .and. (r_eff <= 0.)) then      
     msg = 'if moment_in eq 2 then input moment r_eff must be greater than 0!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif ((moment_in == 13) .and. ((q_h <= 0.) .or. (n_tot <= 0.))) then      
     msg = 'if moment_in eq 13 then input moment q_h and n_tot must be greater than 0!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif ((moment_in == 23) .and. ((q_h <= 0.) .or. (r_eff <= 0.))) then      
     msg = 'if moment_in eq 23 then input moment q_h and r_eff must be greater than 0!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif ((moment_in == 12) .and. ((n_tot <= 0.) .or. (r_eff <= 0.))) then      
     msg = 'if moment_in eq 12 then input moment n_tot and r_eff must be greater than 0!'
     errorstatus = fatal
     call report(errorstatus, msg, nameOfRoutine)
     return
  else
     err = success
  endif


  ! print*, dist_name, "ZZZ", len(dist_name)
  ! print*, trim(dist_name), trim(dist_name) == 'mono', dist_name == 'mono', dist_name(1:len(dist_name)) == 'mono',&
  !   dist_name(1:6) == 'mono', dist_name(1:len(dist_name))
  ! 
  ! print*, trim(dist_name) == "'mono'", dist_name(2:5),dist_name(2:5) == "mono"
  ! STOP

  if (verbose >= 5) then
      print*, 'dist:', dist_name, 'trim(dist):', trim(dist_name), 'moment_in', moment_in
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! MONODISPERSE distribution   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((trim(dist_name) == 'mono') .or. (trim(dist_name) == 'mono_cosmo_ice') .or. &
    (trim(dist_name) == 'mono_echam_cl') .or. (trim(dist_name) == 'mono_echam_ice')) then
     ! Set parameter for the Gamma dist. to get a monodisperse dist.
     if (verbose >= 7) print*, 'Distribution: Some kind of mono'
     lambda = 0._dbl
     mu = 0._dbl
     gam = 0._dbl
     if (trim(dist_name) == 'mono') then
        if (verbose >= 7) print*, 'Distribution: mono'
        ! ! fixed radius (via d_1)
        if ((d_1 /= -99.) .and. (p_1 == -99.)) then
           if (verbose >= 7) print*, 'Distribution: mono, fixed radius (via d_1)'
           d_mono = d_1
           if (moment_in == 3)    n_0 = q_h / (delta_d_mono * a_ms * d_1**b_ms)
           if (moment_in == 1)    n_0 = n_tot / delta_d_mono
        endif
        ! ! fixed n_tot (via p_1)
        if ((p_1 /= -99.) .and. (d_1 == -99.)) then
           if (verbose >= 7) print*, 'Distribution: mono, fixed n_tot (via p_1)'
           n_0 = p_1 / delta_d_mono
           if (moment_in == 3)    d_mono = (q_h / (p_1 * a_ms))**(1._dbl / b_ms)
           if (moment_in == 2)    d_mono = r_eff * 2._dbl
        endif
        ! ! mass_conc and Reff from input
        if (moment_in == 23) then
           if (verbose >= 7) print*, 'Distribution: mono, mass_conc and Reff from input'
           d_mono = r_eff * 2._dbl
           n_0 = q_h / (delta_d_mono * a_ms * d_mono**b_ms)
        endif
     endif
     if ((trim(dist_name) == 'mono_cosmo_ice') .and. (moment_in == 3 .or. moment_in == 13)) then
        if (verbose >= 7) print*, 'Distribution: mono_cosmo_ice'
        ! ! Monodisperse size distribution coherent with COSMO-de 1-moment scheme
        ! Number_concentration of activated ice crystal is temperature dependent
        ! from COSMO-de code src_gscp.f90 routine: hydci_pp_gr
        ! Radius is derived from mass-size relation m=aD^3
        ! a=130 kg/m^3 (hexagonal plates with aspect ratio of 0.2 -> thickness=0.2*Diameter)
        work1 = 1.d2 * exp(0.2_dbl * (273.15_dbl - layer_t)) ! N_tot
        n_0 = work1 / delta_d_mono
        d_mono = (q_h / (work1 * a_ms))**(1._dbl / b_ms)
        !  CHECK if dia1 > maxdiam=2.d-4 (maximum diameter for COSMO)
        !  then recalculate the drop mass using 2.d-4 as particle diameter
        if (d_mono > 2.d-4) then 
           d_mono = 2.d-4
           n_0 = q_h / (a_ms * d_mono**b_ms) / delta_d_mono
        endif
     endif
     if (trim(dist_name) == 'mono_echam_cl') then
        if (verbose >= 7) print*, 'Distribution: mono_echam_cl'
        ! mono disperse distribution coherent with ECHAM6 1 moment scheme for liquid clouds
        if (sfc_type(i_x,i_y) > 0.5) then ! land
          n_0 = 220.d6  ! per m^-3
          kappa = 1.143
        else
          n_0 = 80.d6 ! per m^-3
          kappa = 1.077
        end if
        if (atmo_press(i_x,i_y,i_z) < 8.d4) then
!        if (atmo_hgt(i_x,i_y,i_z) > 1.d3) then
         ! n_0 = n_0 * exp(log(50._dbl/n_0)/10.d3*atmo_hgt(i_x,i_y,i_z))
          n_0 = (50._dbl + (n_0/1.d6 - 50._dbl) * exp(1-(8.d4/max(1.d4,atmo_press(i_x,i_y,i_z)))**2))* 1.d6
        end if
        rho_dry = rho_air(atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z))
        d_mono = kappa*2._dbl*(q_h*rho_dry*3._dbl/(4._dbl*pi*n_0*rho_water))**(1._dbl/3._dbl)
        if (d_mono > 48.d-6) then
          d_mono = 48.d-6
          n_0 = 6.*q_h*rho_dry/(pi*d_mono**3)
        end if 
     endif
     if (trim(dist_name) == 'mono_echam_ice') then
        if (verbose >= 7) print*, 'Distribution: mono_echam_ice'
        ! mono disperse distribution coherent with ECHAM6 1 moment scheme for liquid clouds

        rho_dry = rho_air(atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z))
        d_mono = 2.d-6*(sqrt(2809._dbl*(83.8_dbl*(1.d3*q_h*rho_dry)**0.216_dbl)**3.+5113188._dbl)-2261._dbl)**(1._dbl/3._dbl) ! mean ice crystal volume diameter
        n_0 = q_h*rho_dry / (delta_d_mono * rho_ms * (pi/6._dbl) * d_mono**3._dbl)
     endif
     ! ! Check that the variables have been filled in
     if ((lambda /= 0._dbl) .or. (mu /= 0._dbl) .or. (gam /= 0._dbl) .or. &
          (n_0 <= 0._dbl) .or. (d_mono <= 0._dbl) ) then
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
  ! CONSTANT distribution   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if ((trim(dist_name) == 'const') .or. (trim(dist_name) == 'const_cosmo_ice')) then 
     ! Set parameter for the Gamma dist. to get a constant dist.
     lambda = 0._dbl
     mu = 0._dbl
     gam = 0._dbl
     delta_d_const = (d_2 - d_1)/(nbin-1)
     if (trim(dist_name) == 'const') then
        ! ! fixed radius (via d_1)
        if ((d_1 /= -99.) .and. (d_2 /= -99.) .and. (p_1 == -99.) .and. (p_2 == -99.) &
             .and. (p_1 == -99.) .and. (p_4 == -99.)) then
           if (moment_in == 3)    n_0 = q_h / (delta_d_const * a_ms * d_1**b_ms)
           if (moment_in == 1)    n_0 = n_tot / delta_d_const
           n_0 = n_0 /(nbin-1)
        endif
     endif
     if ((trim(dist_name) == 'const_cosmo_ice') .and. (moment_in == 3) &
          .and. (nbin == 2)) then
        ! ! Monodisperse size distribution coherent with COSMO-de 1-moment scheme
        ! Number_concentration of activated ice crystal is temperature dependent
        ! from COSMO-de code src_gscp.f90 routine: hydci_pp_gr
        ! Radius is derived from mass-size relation m=aD^3
        ! a=130 kg/m^3 (hexagonal plates with aspect ratio of 0.2 -> thickness=0.2*Diameter)
        work1 = 1.d2 * exp(0.2_dbl * (273.15_dbl - layer_t)) ! N_tot
        n_0 = work1 / delta_d_mono
        d_mono = (q_h / (work1 * a_ms))**(1._dbl / b_ms)
        !  CHECK if dia1 > maxdiam=2.d-4 (maximum diameter for COSMO)
        !  then recalculate the drop mass using 2.d-4 as particle diameter
        if (d_mono > 2.d-4) then 
           d_mono = 2.d-4
           n_0 = q_h / (a_ms * d_mono**b_ms) / delta_d_mono
        endif
        d_1 = d_mono - (0.5_dbl* delta_d_mono)
        d_2 = d_mono + (0.5_dbl* delta_d_mono)
     endif

     ! ! Check that the variables have been filled in
     if ((lambda /= 0._dbl) .or. (mu /= 0._dbl) .or. (gam /= 0._dbl) .or. &
          (n_0 <= 0._dbl) .or. ISNAN(n_0) .or. (n_0 >= HUGE(n_0))) then
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
  if ((trim(dist_name) == 'exp') .or. (trim(dist_name) == 'exp_field_t') .or.       &
       (trim(dist_name) == 'exp_cosmo_snow') .or. (trim(dist_name) == 'exp_ryan'))     then 
     ! Set parameter for the Gamma dist. to get an exponential dist.
     gam = 1._dbl
     mu = 0._dbl
     if (trim(dist_name) == 'exp') then
        ! !  everything fixed 
        if ((p_1 /= -99.) .and. (p_2 /= -99.) .and. (moment_in == 0)) then
           lambda = p_1
           n_0 = p_2
        end if
        ! !  inverse expenonetial like ecmwf/ifs: p_1 := n_ax, p_2 := n_bx 
        if ((p_1 /= -99.) .and. (p_2 /= -99.) .and. (moment_in == 3)) then
           lambda = (p_1 * a_ms * gamma(b_ms+1._dbl) / q_h)**(1._dbl / (b_ms + 1 - p_2))
           n_0 = p_1 * lambda ** p_2
        end if        
        ! ! fixed n_tot (via p_1)
        if ((p_1 /= -99.) .and. (p_2 == -99.) .and. (p_3 == -99.)) then
           if (moment_in == 3)  lambda = (p_1 * a_ms * gamma(b_ms+1._dbl) / q_h)**(1._dbl / b_ms)
           if (moment_in == 2)  lambda = 3._dbl / r_eff
           n_0    = p_1 * lambda
        endif
        ! ! fixed r_eff (via p_2)
        if ((p_1 == -99.) .and. (p_2 /= -99.) .and. (p_3 == -99.)) then
           lambda = 3._dbl / p_2
           if (moment_in == 3)  n_0 = (q_h * lambda**(b_ms+1._dbl)) / (a_ms * gamma(b_ms+1._dbl))
           if (moment_in == 1)  n_0 = n_tot * lambda
        endif
        ! ! fixed N_0 (via p_3)
        if ((p_1 == -99.) .and. (p_2 == -99.) .and. (p_3 /= -99.)) then
           n_0 = p_3
           if (moment_in == 3)  lambda = (n_0 * a_ms * gamma(b_ms+1._dbl) / q_h)**(1._dbl / (b_ms+1._dbl))
           if (moment_in == 1)  lambda = n_0 / n_tot
           if (moment_in == 2)  lambda = 3._dbl / r_eff
        endif
        ! ! both moments from input file
        if ((p_1 == -99.) .and. (p_2 == -99.)) then
           if (moment_in == 13)  then
              lambda = (n_tot * a_ms * gamma(b_ms+1._dbl) / q_h)**(1._dbl / b_ms)
              n_0    = n_tot * lambda
           endif
           if (moment_in == 23)  then
              lambda = 3._dbl / r_eff
              n_0    = (q_h * lambda**(b_ms+1)) / (a_ms * gamma(b_ms+1._dbl))
           endif
           if (moment_in == 12)  then
              lambda = 3._dbl / r_eff
              n_0    = n_tot * lambda
           endif
        endif
     endif
     ! ! Ryan (2000 JAS) Lambda = Lambda(layer_t)
     if (trim(dist_name) == 'exp_ryan') then
        lambda = 1220._dbl * 10._dbl**(-0.0245_dbl * (273.15_dbl-layer_t))
        if (moment_in == 3) n_0 = (q_h * lambda**(b_ms+1)) / (a_ms * gamma(b_ms+1._dbl))
        if (moment_in == 1) n_0 = n_tot * lambda
     endif
     ! ! Field et al. (2005 QJRM, end of page 2008 + end of page 2009 for the relation between N_0 and N_0,23) n_0 = n_0(T)
     if (trim(dist_name) == 'exp_field_t') then
        n_0 = 7.628d6 * exp(0.107_dbl * (273.15_dbl - layer_t))
        if (moment_in == 3) lambda = (a_ms * n_0 * gamma(b_ms+1._dbl) / q_h)**(1._dbl /(b_ms+1._dbl))
        if (moment_in == 2) lambda = 3._dbl / r_eff
     endif
     ! ! Field et al. (2005 QJRM, page 2007) n_0 = n_0(T,q_h)
     ! ! AS IMPLEMENTED IN cosmo 2-MOMENT scheme
     if ((trim(dist_name) == 'exp_cosmo_snow') .and. (moment_in == 3)) then
        ! taken from COSMO-de routine hydci_pp_gr in src_gscp.f90
        mma = (/   5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
             0.312550,  0.000204,  0.003199, 0.000000, -0.015952 /)
        mmb = (/   0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
             0.060366,  0.000079,  0.000594, 0.000000, -0.003577 /)
        ! Calculate n0s using the temperature-dependent moment
        ztc = layer_t - 273.15                     ! temperature in C
        ztc = MAX(MIN(ztc,0.0),-40.0)              !limited to -40
        ! the moments of every order of the particle size distribution f(D) depend from the second one (proportional to mass) and temperature
        nn  = 3
        hlp = mma(1)      +mma(2)*ztc      +mma(3)*nn       +mma(4)*ztc*nn+mma(5)*ztc**2 &
             + mma(6)*nn**2+mma(7)*ztc**2*nn+mma(8)*ztc*nn**2+mma(9)*ztc**3+mma(10)*nn**3
        alf = 10.0d0**hlp
        bet = mmb(1)      +mmb(2)*ztc      +mmb(3)*nn       +mmb(4)*ztc*nn+mmb(5)*ztc**2 &
             + mmb(6)*nn**2+mmb(7)*ztc**2*nn+mmb(8)*ztc*nn**2+mmb(9)*ztc**3+mmb(10)*nn**3
        m2s = q_h / a_ms ! Formfactor in the mass-size relation of snow particles [kg/m^2] (COSMO value =0.038)
        m3s = alf*EXP(bet*LOG(m2s))
        hlp  =  7.628d6 * EXP(-0.107d0*ztc) ! N0snow as a function of solely T
        n_0 = 13.5 * m2s**4 / m3s**3       ! N0snow as a function of T and snow mixing ratio
        n_0 = MAX(n_0,0.5*hlp)
        n_0 = MIN(n_0,1e2*hlp)
        n_0 = MIN(n_0,1e9)
        n_0 = MAX(n_0,1e6)
        lambda = (a_ms * n_0 * gamma(b_ms+1._dbl) / q_h)**(1._dbl /(b_ms+1._dbl))
     endif
     ! ! Check that the variables have been filled in
     if ((gam /= 1._dbl) .or. (mu /= 0._dbl) .or. (n_0 <= 0._dbl) .or. (lambda <= 0._dbl)) then
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
     if ((p_1 /= -99.) .and. (p_3 /= -99.)) then
        n_t = p_1
        sig = p_3
        if (moment_in == 3) d_ln = (1._dbl/b_ms) * (log(q_h / (a_ms * n_t)) - (0.5_dbl * (b_ms * sig)**2._dbl))
        if (moment_in == 2) d_ln = log(2._dbl * r_eff) - 2.5_dbl * sig**2._dbl
     endif
     ! ! fixed n_tot (via p_1) and r_eff (via p_2)
     if ((p_1 /= -99.) .and. (p_2 /= -99.) .and. (p_3 == -99.) .and. (moment_in == 3)) then
        n_t = p_1
        work1 = q_h / (a_ms * n_t * (2._dbl * p_2)**2._dbl)
        work2 = 2._dbl / (b_ms*(b_ms-5._dbl))
        sig = sqrt(log(work1) * work2)
        d_ln = log(2._dbl * p_2) - 2.5_dbl * sig**2._dbl
     endif
     ! ! fixed r_eff (via p_2) and sigma (via p_3)
     if ((p_1 == -99.) .and. (p_2 /= -99.) .and. (p_3 /= -99.)) then
        sig = p_3
        d_ln = log(2._dbl * p_2) - 2.5_dbl * sig**2._dbl
        if (moment_in == 3) n_t = q_h / a_ms * exp(-1._dbl*((b_ms * d_ln) + (0.5_dbl *(b_ms * sig)**2._dbl)))
        if (moment_in == 1) n_t = n_tot
     endif
     ! ! fixed sigma (via p_3)
     if ((p_1 == -99.) .and. (p_2 == -99.) .and. (p_3 /= -99.)) then
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
     if ((n_t <= 0._dbl) .or. isnan(d_ln) .or. (sig <= 0._dbl)) then
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
  if (trim(dist_name) == 'mgamma' .or. trim(dist_name) == 'mgamma_MNH') then
     ! ! The user MUST specify mu and gam parameters
     if ((p_3 == -99.) .or. (p_4 == -99.)) then 
        msg = 'Modified Gamma case: p_3 and p_4 parameters must be specified...'
        errorstatus = fatal
        call report(errorstatus, msg, nameOfRoutine)
        return
     endif
     mu  = p_3
     gam = p_4
     ! ! fixed n_tot (via p_1)
     if ((p_1 /= -99.) .and. (p_2 == -99.)) then
        if (moment_in == 3)  then
           work2 = gamma((mu + b_ms + 1._dbl) / gam)
           work3 = gamma((mu + 1._dbl) / gam)
           lambda = (a_ms / q_h * p_1 * work2 / work3)**(gam / b_ms)
           n_0 = gam * p_1 / work3 * lambda**((mu + 1._dbl) / gam)
        endif
        if (moment_in == 2)  then
           work1 = gamma((mu + 4._dbl) / gam)
           work2 = gamma((mu + 3._dbl) / gam) * r_eff * 2._dbl
           lambda = (work1 / work2)**gam
           work1 = (mu + 1._dbl) / gam
           n_0 = p_1 * gam * lambda**work1 / gamma(work1)
        endif
     endif
     ! ! fixed r_eff (via p_2)
     if ((p_1 == -99.) .and. (p_2 /= -99.)) then
        if (moment_in == 3)  then
           work1 = gamma((mu + 4._dbl) / gam)
           work2 = gamma((mu + 3._dbl) / gam) * p_2 * 2._dbl
           lambda = (work1 / work2)**gam
           work1 = q_h * gam / a_ms
           work2 = (mu + b_ms + 1._dbl) / gam
           n_0 = work1 * lambda**work2 / gamma(work2)
        endif
        if (moment_in == 1)  then
           work1 = gamma((mu + 4._dbl) / gam)
           work2 = gamma((mu + 3._dbl) / gam) * p_2 * 2._dbl
           lambda = (work1 / work2)**gam
           work1 = (mu + 1._dbl) / gam
           n_0 = n_tot * gam * lambda**work1 / gamma(work1)
        endif
     endif
     ! ! MESO-NH distribution
     ! Ntot = C * lambda^x
     ! C = p_1 and x = p_2
     if (trim(dist_name) == 'mgamma_MNH') then
        if ((p_1 < 0.) .or. (p_2 == -99.)) then
           msg = 'Modified-gamma case: something wrong or this parameters combination is not yet implemented...'
           errorstatus = fatal
           call report(errorstatus, msg, nameOfRoutine)
           return
        endif
        work2 = gamma((mu + b_ms + 1._dbl) / gam)
        work3 = gamma((mu + 1._dbl) / gam)
        lambda = (a_ms * p_1 / q_h * work2 / work3)**(gam / (b_ms - p_2))
        n_0 = gam * p_1 / work3 * lambda**((mu + 1._dbl + p_2) / gam)
     endif
     ! ! 2 moments from the input file
     if ((p_1 == -99.) .and. (p_2 == -99.)) then
        if (moment_in == 13)  then
           work2 = gamma((mu + b_ms + 1._dbl) / gam)
           work3 = gamma((mu + 1._dbl) / gam)
           lambda = (a_ms / q_h * n_tot * work2 / work3)**(gam / b_ms)
           n_0 = gam * n_tot / work3 * lambda**((mu + 1._dbl) / gam)
        endif
        if (moment_in == 23)  then
           work1 = gamma((mu + 4._dbl) / gam)
           work2 = gamma((mu + 3._dbl) / gam) * r_eff * 2._dbl
           lambda = (work1 / work2)**gam
           work1 = q_h * gam / a_ms
           work2 = (mu + b_ms + 1._dbl) / gam
           n_0 = work1 * lambda**work2 / gamma(work2)
        endif
        if (moment_in == 12)  then
           work1 = gamma((mu + 4._dbl) / gam)
           work2 = gamma((mu + 3._dbl) / gam) * r_eff * 2._dbl
           lambda = (work1 / work2)**gam
           work1 = (mu + 1._dbl) / gam
           n_0 = n_tot * gam * lambda**work1 / gamma(work1)
        endif
     endif
     ! ! Check that the variables have been filled in
     if ((lambda <= 0._dbl) .or. (n_0 <= 0._dbl) .or. isnan(mu) .or. isnan(gam)) then
        msg = 'Modified-gamma case: something wrong or this parameters combination is not yet implemented...'
        errorstatus = fatal
        call report(errorstatus, msg, nameOfRoutine)
        return
     endif
     errorstatus = err
     if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
     return
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NORMALIZED GAMMA distribution   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (trim(dist_name) == 'norm_gamma') then
     call assert_true(err,(moment_in == 0) .or. (moment_in == 23),&
          'Normalized Modified Gamma case: currently only implemented for moment_in = 0 and 23')
     if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
     end if
     if (moment_in == 0) then 
        ! ! The user MUST specify d_m, n_0_star and mu parameters
        call assert_false(err,((p_1 == -99.) .or. (p_2 == -99.)),&
             'Normalized Modified Gamma case: p_1 and p_2 parameters must be specified')
        call assert_true(err,((p_3 /= -99.) .NEQV. (p_4 /= -99.)),& ! NEQV = xor
             'Normalized Modified Gamma case: p_3 xor p_4 parameters must be specified...' )
        if (err > 0) then
           errorstatus = fatal
           msg = "assertation error"
           call report(errorstatus, msg, nameOfRoutine)
           return
        end if
        d_m  = p_1
        n_0_star = p_2
        if (p_3 /= -99.) mu = p_3
        if (p_4 /= -99.) mu = p_4 -(b_ms+1) !shifted mu value for better numerical handling in optimal estimation!
     elseif (moment_in == 23) then 
        ! attention, the normalized form requires the ration of the ratio of the (b_ms+1)nth to the (b_ms)nth moment, reff is usually 3rd to 2nd moment
        d_m  = r_eff * 2 

        !eq A14 of Maahn et al 2015:
        work1 = gamma(b_ms+1)
        work2 = work1/(b_ms+1)**(b_ms+1)
        !eq8 of Testud 2001 & eq A9 of Maahn et al 2015:
        n_0_star = q_h/(a_ms * d_m**(b_ms+1) * work2)


        if (verbose >= 5)  print*, q_h, d_m, n_0_star, a_ms, b_ms, work1, work2

        if (p_3 /= -99.) then
          mu = p_3
        else
         call assert_true(err,(liq_ice == 1),& 
             'Mu can be only estimated from d_m for rain' )
          !Mroz et al. 2023 eq 21, for rain only!
          mu = 10. * (d_m*1000.)**(-0.8) - 4.
        end if  


     end if

     if (verbose >= 5)  print*, d_m, n_0_star, mu

     if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
     end if
     errorstatus = err
     if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
     return
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! NORMALIZED MODIFIED GAMMA distribution   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (trim(dist_name) == 'norm_mgamma') then
     ! ! The user MUST specify d_m, n_0_star, mu and gamparameters
     call assert_false(err,((p_1 == -99.) .or. (p_2 == -99.) .or.&
                            (p_3 == -99.) .or. (p_4 == -99.)),&
          'Normalized Modified Gamma case: p_1 to p_4 parameters must be specified')
     call assert_true(err,(moment_in == 0),&
          'Normalized Modified Gamma case: currently only implemented for moment_in = 0')
     if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
     end if
     d_m  = p_1
     n_0_star = p_2
     mu = p_3
     gam = p_4
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
