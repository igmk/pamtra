! Subroutine for the setup of the parameters of the snow particle size distribution.
!
!
subroutine snow_ssp(f,swc,t,maxleg,kext, salb, back,  &
     nlegen, legen, legen2, legen3, legen4, qs, nc)

  use kinds
  use nml_params, only: verbose, lphase_flag, n_0snowDsnow, EM_snow, n_moments, isnow_n0, SD_snow
  use constants, only: pi, im
  use double_moments_module
  use conversions

  implicit none

  integer :: nbins, nlegen, nn

  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
       swc,&
       qs,&
       t,&
       f

  real(kind=dbl), optional, intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, ad, bd, alpha, gamma, b_snow, a_msnow

  real(kind=dbl), intent(out) :: &
       kext,&
       salb,&
       back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex, m_air

  real(kind=dbl) :: gammln

  real(kind=dbl), dimension(10) :: mma, mmb

  real(kind=dbl) :: ztc, hlp, alf, bet, m2s, m3s

  if (verbose .gt. 1) print*, 'Entering snow_ssp'

  call ref_ice(t, f, refre, refim)
  mindex = refre-Im*refim  ! mimicking a
  m_air = 1.0d0 - 0.0d0 * Im

  if (n_moments .eq. 1) then
  	if (SD_snow .eq. 'C') then
     !	b_snow = 2.0d0
     !	a_msnow = 0.038d0
     b_snow = 2.2850d0
     a_msnow = 0.2124d0

     dia1 = 0.51d-5 ! minimum maximum diameter [m] after kneifel
     dia2 = 1.d-2 ! maximum maximum diameter [m] after kneifel

     !option isnow_n0 as in COSMO-de 1 moment scheme
     !isnow_n0temp = 2 intercept parameter of snow depend on T and qs (snow mixing ratio) Field 2005
     !isnow_n0temp = 1 intercept parameter of snow depend on T  	Field 2005

     if (isnow_n0 .eq. 2) then
        !//taken from COSMO-de routine hydci_pp_gr in src_gscp.f90
        ! Coeffs for moment relation based on 2nd moment (Field 2005)
        mma = (/   5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
             0.312550,  0.000204,  0.003199, 0.000000, -0.015952 /)
        mmb = (/   0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
             0.060366,  0.000079,  0.000594, 0.000000, -0.003577 /)
        ! Calculate n0s using the temperature-dependent moment
        ! relations of Field et al. (2005)
        ztc = t - 273.15				! temperature in C
        ztc = MAX(MIN(ztc,0.0),-40.0) 	!limited to -40

        !formula from Field (2005)
        !the moments of every order of the particle size distribution f(D) depend from the second one (proportional to mass) and temperature
        !N0 for snow depends on both temperature and snow mixing ratio
        nn  = 3
        hlp = mma(1)      +mma(2)*ztc      +mma(3)*nn       +mma(4)*ztc*nn+mma(5)*ztc**2 &
             + mma(6)*nn**2+mma(7)*ztc**2*nn+mma(8)*ztc*nn**2+mma(9)*ztc**3+mma(10)*nn**3
    	alf = 10.0d0**hlp
    	bet = mmb(1)      +mmb(2)*ztc      +mmb(3)*nn       +mmb(4)*ztc*nn+mmb(5)*ztc**2 &
             + mmb(6)*nn**2+mmb(7)*ztc**2*nn+mmb(8)*ztc*nn**2+mmb(9)*ztc**3+mmb(10)*nn**3
    	m2s = qs / 0.038d0 ! 0.038 = Formfactor in the mass-size relation of snow particles [kg/m^3]
    	m3s = alf*EXP(bet*LOG(m2s))
    	hlp  = n_0snowDsnow * 1.d6*EXP(-0.107d0*ztc)	! N0snow as a function of solely T from Field (2005)				[1/m^4]
    	ad = 13.5 * m2s**4 / m3s**3						! N0snow as a function of T and snow mixing ratio from Field (2005)	[1/m^4]
        ad = MAX(ad,0.5*hlp)
        ad = MIN(ad,1e2*hlp)
        ad = MIN(ad,1e9)
        ad = MAX(ad,1e6)
        !//end of COSMO code
     else
        ! Field param. ! multiplied by 10^6 is 1/m^4
        ad = n_0snowDsnow * 1.d6 * exp(-0.107d0 * (t - 273.15))
     endif
     bd = (exp(gammln(b_snow + 1)) * a_msnow * ad/swc)**(1.0d0/(1.0d0 + b_snow))  ! [m**-1]
     nbins = 100
     alpha = 0.d0 ! exponential SD
     gamma = 1.d0
    else if (SD_snow .eq. 'M') then
     dia1 = 0.51d-5 ! minimum maximum diameter [m] after kneifel
     dia2 = 1.d-2 ! maximum maximum diameter [m] after kneifel

     a_msnow = 0.02d0
     b_snow = 1.9d0
     bd = (swc/(a_msnow*5.d0*exp(gammln(b_snow+1.))))**(1.d0/(1.d0-b_snow))  ! [m**-1]
	 ad = 5.*bd**2.
     nbins = 100
     alpha = 0.d0 ! exponential SD
     gamma = 1.d0
	end if
  else if (n_moments .eq. 2) then
     if (.not. present(nc)) stop 'STOP in routine snow_ssp'
     call double_moments(swc,nc,gamma_snow(1),gamma_snow(2),gamma_snow(3),gamma_snow(4), &
          ad,bd,alpha,gamma,a_msnow,b_snow)
     nbins = 100
     dia1 = 1.d-6 ! minimum maximum diameter [m] after kneifel
     dia2 = 2.d-2 ! maximum maximum diameter [m] after kneifel
  else
     stop'Number of moments is not specified'
  end if

  if (EM_snow .eq. 'icesf') then 
     call mie_densitysizedep_spheremasseq(f, mindex,      &
          a_msnow, b_snow, dia1, dia2, nbins, maxleg,   &
          ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
          back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
          LEGEN4, SD_snow)
  elseif (EM_snow .eq. 'surus') then 
    print*, "here"
     call mie_icefactor(f, t,mindex,      &
          a_msnow, b_snow, dia1, dia2, nbins, maxleg,   &
          ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
          back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
          LEGEN4, SD_snow,0.863*1.e-3*f+0.115,42)
  elseif (EM_snow .eq. 'liudb') then
     dia1 = 1.02d-4
     dia2 = 2.d-2
     call dda_db_liu(f,t,mindex, &
          dia1,dia2,nbins,maxleg,ad,&
          bd, alpha, gamma, lphase_flag,kext, salb,&
          back, nlegen, legen, legen2, legen3,&
          legen4, SD_snow)
  else 
     write (*, *) 'no em mod', EM_snow
     stop
  endif

  if (verbose .gt. 1) print*, 'Exiting snow_ssp'

  return

end subroutine snow_ssp
