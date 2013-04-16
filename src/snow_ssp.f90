! Subroutine for the setup of the parameters of the snow particle size distribution.
!
!
subroutine snow_ssp(f,swc,t,press,maxleg,nc,kext, salb, back,  &
     nlegen, legen, legen2, legen3, legen4,&
     scatter_matrix,extinct_matrix, emis_vector,snow_spec)

  use kinds
  use settings, only: lphase_flag, n_0snowDsnow, EM_snow, &
	n_moments, isnow_n0, SD_snow, snow_density,liu_type,nstokes,&
        radar_nfft_aliased, radar_mode, active
  use constants, only: pi, im
  use double_moments_module
  use conversions
        use report_module

  implicit none

  integer :: nbins, nbins_spec, nlegen, nn,alloc_status

  integer, intent(in) :: maxleg
    integer, parameter ::  nquad = 16
  real(kind=dbl), intent(in) :: &
       swc,&
       t,&
       f,press

  real(kind=dbl), intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, ad, bd, alpha, gamma, b_snow, a_msnow

  real(kind=dbl), intent(out) :: &
       kext,&
       salb,&
       back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4
    real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
    real(kind=dbl), dimension(nstokes,nstokes,nquad,2), intent(out) :: extinct_matrix
    real(kind=dbl), dimension(nstokes,nquad,2), intent(out) :: emis_vector
  complex(kind=dbl) :: mindex, m_air
  real(kind=dbl), allocatable, dimension(:):: diameter_spec, back_spec
  real(kind=dbl) :: gammln

  real(kind=dbl), dimension(10) :: mma, mmb

  real(kind=dbl) :: ztc, hlp, alf, bet, m2s, m3s, a_as_snow, b_as_snow
  real(kind=dbl), intent(out), dimension(radar_nfft_aliased) :: snow_spec
  character(5) ::  particle_type

  if (verbose .gt. 1) print*, 'Entering snow_ssp'
  if ((n_moments .eq. 1) .and. (EM_snow .eq. "tmatr")) stop "1moment tmatr not tested yet for snow"

  call ref_ice(t, f, refre, refim)


  mindex = refre-Im*refim  ! mimicking a

  m_air = 1.0d0 - 0.0d0 * Im

  if (n_moments .eq. 1) then
  	if (SD_snow .eq. 'C') then
     	b_snow = 2.0d0
     	a_msnow = 0.038d0 ! Locatelli and Hobbs (1974)
	!area-size relation in SI
        a_as_snow = 0.3970874893692325d0 !0.2285 in CGS
        b_as_snow = 1.88d0 !from mitchell 1996 similar to a_msnow&b_snow
!     b_snow = 2.2850d0
!     a_msnow = 0.2124d0

     dia1 = 0.51d-10 ! minimum maximum diameter [m] after kneifel
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
    	m2s = swc / 0.038d0 ! 0.038 = Formfactor in the mass-size relation of snow particles [kg/m^2]
    	m3s = alf*EXP(bet*LOG(m2s))
    	hlp  = n_0snowDsnow * 1.d6*EXP(-0.107d0*ztc)	! N0snow as a function of solely T from Field (2005)				[1/m^4]
    	ad = 13.5 * m2s**4 / m3s**3						! N0snow as a function of T and snow mixing ratio from Field (2005)	[1/m^4]
        ad = MAX(ad,0.5*hlp)
        ad = MIN(ad,1e2*hlp)
        ad = MIN(ad,1e9)
        ad = MAX(ad,1e6)
        !//end of COSMO code
     else if (isnow_n0 .eq. 1) then
        ! Field param. ! multiplied by 10^6 is 1/m^4
        ad = n_0snowDsnow * 1.d6 * exp(-0.107d0 * (t - 273.15))

     else if (isnow_n0 .eq. 0) then
        !fixed N0
        ad = n_0snowDsnow * 1.d6
     else
        print *, "isnow_n0: ", isnow_n0
        stop "Wrong isnow_n0 value"
     endif
     bd = (exp(gammln(b_snow + 1)) * a_msnow * ad/swc)**(1.0d0/(1.0d0 + b_snow))  ! [m**-1]
     alpha = 0.d0 ! exponential SD
     gamma = 1.d0
    else if (SD_snow .eq. 'M') then
     dia1 = 0.51d-10 ! minimum maximum diameter [m] after kneifel
     dia2 = 1.d-2 ! maximum maximum diameter [m] after kneifel

     a_msnow = 0.02d0
     b_snow = 1.9d0
     !area-size relation in SI
     a_as_snow = 0.3970874893692325d0 !0.2285 in CGS
     b_as_snow = 1.88d0 !from mitchell 1996 similar to a_msnow&b_snow
     bd = (swc/(a_msnow*5.d0*exp(gammln(b_snow+1.))))**(1.d0/(1.d0-b_snow))  ! [m**-1]
	 ad = 5.d0*bd**2.d0
	 dia2 = log(ad)/bd
     if (dia2 .gt. 2.d-2) dia2 = 2.d-2
     alpha = 0.d0 ! exponential SD
     gamma = 1.d0
    else
      print*, "did not understand SD_snow: ", SD_snow
      stop
     end if
  else if (n_moments .eq. 2)  then
     if ((nc .eq. 0.d0) .or. (SD_snow .ne. "G"))stop 'STOP in routine snow_ssp'
     call double_moments(swc,nc,gamma_snow(1),gamma_snow(2),gamma_snow(3),gamma_snow(4), &
          ad,bd,alpha,gamma,a_msnow,b_snow)
     dia1 = 1.d-10
     dia2 = 2.d-2
     !to do: implement area-size in double_moments file!
     a_as_snow = 0.3970874893692325d0 !0.2285 in CGS
     b_as_snow = 1.88d0 !from mitchell 1996 similar to a_msnow&b_snow
  else
     stop'Number of moments is not specified'
  end if


  if ((EM_snow .eq. 'densi') .or. (EM_snow .eq. 'surus')) then
    nbins = 100
    nbins_spec = nbins+1 !Mie routine uses nbins+1 bins!
  elseif (EM_snow .eq. 'liudb') then
    nbins = 100
    nbins_spec = nbins
  elseif (EM_snow .eq. 'tmatr') then
    nbins = 50
    nbins_spec = nbins
  else
     write (*, *) 'no em mod', EM_snow
     stop
  end if

  allocate(diameter_spec(nbins_spec),stat=alloc_status)
  allocate(back_spec(nbins_spec),stat=alloc_status)

  if (EM_snow .eq. 'densi' .or. EM_snow .eq. 'surus') then
   if (EM_snow .eq. 'surus') snow_density = 0.863*f+115.d0
     call mie_densitydep_spheremasseq(f, t,mindex,      &
          a_msnow, b_snow, dia1, dia2, nbins, maxleg,   &
          ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
          back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
          LEGEN4, SD_snow,snow_density,swc,&
          diameter_spec, back_spec)
      scatter_matrix= 0.d0
      extinct_matrix= 0.d0
      emis_vector= 0.d0
  elseif (EM_snow .eq. 'liudb') then
     dia1 = 1.02d-4
     dia2 = 2.d-2
     call dda_db_liu(f,t,liu_type,mindex, &
          dia1,dia2,nbins,maxleg,ad,&
          bd, alpha, gamma, lphase_flag,kext, salb,&
          back, nlegen, legen, legen2, legen3,&
          legen4, SD_snow,&
          diameter_spec, back_spec)
      scatter_matrix= 0.d0
      extinct_matrix= 0.d0
      emis_vector= 0.d0
  elseif (EM_snow .eq. 'tmatr') then
    call tmatrix_snow(f, swc, t, nc, &
          ad, bd, alpha, gamma, a_msnow, b_snow, SD_snow, nbins, scatter_matrix,extinct_matrix, emis_vector,&
          diameter_spec, back_spec)
    back = scatter_matrix(1,16,1,16,2) !scatter_matrix(A,B;C;D;E) backscattering is M11 of Mueller or Scattering Matrix (A;C=1), in quadrature 2 (E) first 16 (B) is 180deg (upwelling), 2nd 16 (D) 0deg (downwelling). this definition is lokkiing from BELOW, scatter_matrix(1,16,1,16,3) would be from above!
    back = 4*pi*back!/k**2 !eq 4.82 Bohren&Huffman without k**2 (because of different definition of Mueller matrix according to Mishenko AO 2000). note that scatter_matrix contains already squard entries!
    kext = extinct_matrix(1,1,16,1) !11 of extinction matrix (=not polarized), at 0°, first quadrature. equal to extinct_matrix(1,1,16,2)

    !not needed by rt4
    salb = 0.d0
    nlegen = 0
    legen = 0.0d0
    legen2 = 0.0d0
    legen3 = 0.0d0
    legen4 = 0.0d0
  else 
     write (*, *) 'no em mod', EM_snow
     stop
  endif

  particle_type="snow" 

  if ((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments"))) then
    call radar_spectrum(nbins_spec,diameter_spec, back, back_spec,t,press,f,&
      particle_type,a_msnow,b_snow,a_as_snow,b_as_snow,snow_spec)
  else
    snow_spec(:)=0.d0
  end if

  deallocate(diameter_spec, back_spec)

  if (verbose .gt. 1) print*, 'Exiting snow_ssp'

  return

end subroutine snow_ssp
