! Subroutine for the setup of the parameters of the snow particle size distribution.
!
!
subroutine snow_ssp(f,qs,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4)
  
  use kinds
  use nml_params, only: verbose, lphase_flag, n_0snowDsnow, EM_snow
  use constants, only: pi, im

  implicit none

  integer :: numrad, nlegen, nn, isnow_n0temp
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qs,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: rad1, rad2, swc, ad, bd, alpha, gamma, b_snow, a_msnow

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex, m_air

  real(kind=dbl) :: spec2abs, gammln
  
  real(kind=dbl), dimension(10) :: mma, mmb

  real(kind=dbl) :: ztc, hlp, alf, bet, m2s, m3s

  if (verbose .gt. 1) print*, 'Entering snow_ssp'

	b_snow = 2.0d0     ! MKS system
	a_msnow = 0.038d0 


	call ref_ice(t, f, refre, refim)
	mindex = refre-Im*refim  ! mimicking a

	m_air = 1.0d0 - 0.0d0 * Im 

	rad1 = 1.d-6 ! minimum maximum diameter [m] after kneifel
	rad2 = 2.d-2 ! maximum maximum diameter [m] after kneifel

	!option isnow_n0temp as in COSMO-de 1 moment scheme
	!isnow_n0temp = 2 intercept parameter of snow depend on T and qs (snow mixing ratio) Field 2005
	!isnow_n0temp = 1 intercept parameter of snow depend on T  	Field 2005
	isnow_n0temp=2

if (isnow_n0temp .eq. 2) then
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
    m2s = qs / 0.038 ! 0.038 = Formfactor in the mass-size relation of snow particles [kg/m^3]
    m3s = alf*EXP(bet*LOG(m2s))
    hlp  = n_0snowDsnow * 1.d6*EXP(-0.107d0*ztc)	! N0snow as a function of solely T from Field (2005)				[1/m^4]
    ad = 13.5 * m2s**4 / m3s**3						! N0snow as a function of T and snow mixing ratio from Field (2005)	[1/m^4]
    ad = MAX(ad,0.5*hlp)
    ad = MIN(ad,1e2*hlp)
    ad = MIN(ad,1e9)
    ad = MAX(ad,1e6)
!//end of COSMO code
elseif (isnow_n0temp .eq. 1) then
	! Field param. ! multiplied by 10^6 is 1/m^4
	ad = n_0snowDsnow * 1.d6 * exp(-0.107d0 * (t - 273.15))
else
	print*, "isnow_n0temp not set"
	stop
endif

	swc =  spec2abs(qs,t,p,q) ! [kg/m^3]
	bd = (exp(gammln(b_snow + 1)) * a_msnow * ad/swc)**(1.0d0/(1.0d0 + b_snow))  ! [m**-1]   

	!formula 3.12 Mario Mech´s but for radii and units converted
	numrad = 100 

	if (EM_snow .eq. 'icesf') then 
	  call mie_densitysizedep_spheremasseq(f, mindex,      &
		a_msnow, b_snow, rad1/2., rad2/2., numrad, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
		LEGEN4, 'C')
	elseif (EM_snow .eq. 'surus') then 
	  call mie_icefactor(f, t,mindex,      &
		a_msnow, b_snow, rad1/2., rad2/2., numrad, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
		LEGEN4, 'C',0.863*1.e-3*f+0.115,42)
	elseif (EM_snow(1:3) .eq. 'liu') then
	    call dda_db_liu(f,t,9,mindex, &
		rad1/2.,rad2/2.,numrad,maxleg,ad,&
		bd, alpha, gamma, lphase_flag,kext, salb,&
		back, nlegen, legen, legen2, legen3,&
		legen4, 'C')
	else 
	    write (*, *) 'no em mod', EM_snow
	    stop
	endif

  if (verbose .gt. 1) print*, 'Exiting snow_ssp'

  return

end subroutine snow_ssp
