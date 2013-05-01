subroutine ice_ssp(f,iwc,t,press,maxleg,nc, kext, salb, back,  &
     nlegen, legen, legen2, legen3, legen4,&
     scatter_matrix,extinct_matrix, emis_vector,ice_spec)

  use kinds

  use settings, only: lphase_flag, n_moments, EM_ice, SD_ice,&
      nstokes, radar_nfft_aliased, radar_mode, active, ad_ice, bd_ice,&
      alphad_ice, gammad_ice, liu_type_ice, diamin_ice, diamax_ice, &
      mass_size_ice_b, mass_size_ice_a, area_size_ice_b, area_size_ice_a,&
      as_ratio_ice

  use constants, only: pi, im
  use double_moments_module
  use conversions
        use report_module

  implicit none

  integer :: nbins, nbins_spec, nlegen,alloc_status

  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
       iwc,&
       t,&
       f,press

  real(kind=dbl), intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, del_d, den_ice, drop_mass, b_mice, a_mice
  real(kind=dbl) :: ad, bd, alpha, gamma, number_concentration,a_as_ice, b_as_ice

  real(kind=dbl), intent(out) :: &
       kext,&
       salb,&
       back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4
    integer, parameter ::  nquad = 16
    real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
    real(kind=dbl), dimension(nstokes,nstokes,nquad,2), intent(out) :: extinct_matrix
    real(kind=dbl), dimension(nstokes,nquad,2), intent(out) :: emis_vector
  complex(kind=dbl) :: mindex
  real(kind=dbl), allocatable, dimension(:):: diameter_spec, back_spec
  real(kind=dbl), intent(out), dimension(radar_nfft_aliased) :: ice_spec

  real(kind=dbl) :: Nt, re

  character(5) ::  particle_type
  if (verbose .gt. 1) print*, 'Entering ice_ssp'

  call ref_ice(t,f, refre, refim)
  mindex = refre-im*refim  ! mimicking a

  if (n_moments .eq. 1) then
	if (SD_ice .eq. 'C') then
     !	Monodisperse size distribution coherent with COSMO-de 1-moment scheme
     !	Number_concentration of activated ice ctystal is temperature dependent
     !	 from COSMO-de code src_gscp.f90 routine: hydci_pp_gr
     !	Radius is derived from mass-size relation m=aD^3
     !	 a=130 kg/m^3 (hexagonal plates with aspect ratio of 0.2 -> thickness=0.2*Diameter)
     number_concentration = 1.0d2*DEXP(0.2d0*(273.15d0-t)) 	! [1/m^3]
     drop_mass = iwc/number_concentration 					! [kg]
     del_d = 1.d-8	
del_d =1.d-9
     dia1 = (drop_mass/130.0d0)**(1.0d0/3.0d0)				! [m]
!    CHECK if dia1 > maxdiam=2.d-4 (maximum diameter for COSMO)
! 	 then recalculate the drop mass using 2.d-4 as particle diameter
	 if (dia1 .gt. 2.d-4) then
	 	dia1 = 2.d-4 										! [m] maximum allowed diameter
		drop_mass = 130.d0 * dia1**3						! [kg]
	 endif
	 dia2 = dia1 + del_d
     ad = iwc/(drop_mass*del_d) 	!intercept parameter [1/m^4]
     bd = 0.0d0
     nbins = 2
nbins=20
     alpha = 0.0d0     ! exponential SD
     gamma = 1.0d0
     den_ice=917.d0

     a_mice = 130d0
     b_mice = 3.d0
     !area-size relation in SI
     a_as_ice = 0.684 !0.684 also in CGS
     b_as_ice = 2.d0 !from mitchell 1996 similar to a_msnow&b_snow
    else if (SD_ice .eq. 'M') then
     del_d = 1.d-8    ! [m]
     dia1 = 6.d-5     ! [m] 60 micron diameter
     dia2 = dia1 + del_d
     den_ice = 917.d0  ! [kg/m^3]
     drop_mass = pi/6.d0 * dia1**3 * den_ice
     a_mice = 0.82d0
     b_mice = 2.5d0
     !area-size relation in SI
     a_as_ice = 0.12028493607054538 !0.24 in CGS
     b_as_ice = 1.85d0 !from mitchell 1996 similar to a_msnow&b_snow
     ad = iwc/(drop_mass*del_d) 	!intercept parameter [1/m^4]
     bd = 0.0d0
     nbins = 2
     alpha = 0.0d0     ! exponential SD
     gamma = 1.0d0 
	else if (SD_ice .eq. 'L') then !lognormal cloud distribution to test Pavlos Radar Spectrum
!teh hard and ugly way:

	    dia1 = 3.d-6 ! [m] 2 micron diameter 
	    dia2 = 16.d-6 ! [m] 50 micron diameter 
	    nbins = 48
	    alpha = 0.3 !S_x in Pavlos code
	    re=10.d-6
	    bd = re/exp(2.5d0*alpha**2); !ro in Pavlos Code
	    Nt = (3*iwc*1.d3)/(4*pi*1*(bd*1d2)**3*exp(4.5d0*alpha**2))
	    ad = Nt/(sqrt(2*pi) *alpha)
	    gamma = 1.d0
     den_ice=917.d0
     a_mice = 130d0
     b_mice = 3.d0
     a_as_ice = 0.12028493607054538 !0.24 in CGS
     b_as_ice = 1.85d0 !from mitchell 1996 similar to a_msnow&b_snow


	else if (SD_ice .eq. 'G') then !lMPACE
	    dia1 = diamin_ice ! [m] 
	    dia2 = diamax_ice! [m] 
	    nbins = 48
	    alpha = alphad_ice!from nml_params
	    bd = bd_ice !from nml_params
	    ad = ad_ice !from nml_params
	    gamma = gammad_ice!from nml_params
	    den_ice=917.d0
	    if (EM_ice .eq. 'liudb') then
	      if (liu_type_ice == 9) then
		b_mice = 1.511d0 !from Stefan Kneifel
		a_mice = 1.191d-3 !in SI
		b_as_ice = -0.377 + 2 !from liu 2004 area ratio
		a_as_ice = (0.261d0*pi/4.d0) * 10d0**(2d0*b_as_ice-4d0)
	      else if (liu_type_ice == 10) then
		b_mice = 1.82d0 !from Stefan Kneifel
		a_mice = 5.666d-3 !in SI
		b_as_ice = 0 + 2 !from liu 2004: I found from the den image a fixed area ratio 0f 0.4 (MX,2013) !random orientation is totally ignored!
		a_as_ice = (0.4*pi/4.d0) * 10d0**(2d0*b_as_ice-4d0)
	      else
		print*, "density for liu_type_ice not defined", liu_type_ice
		STOP
	      end if
	    else
! 	      b_mice = 1.7d0
! 	      a_mice = 1.07d-10 * 10**(6*b_mice - 3) !in SI
	      b_mice = mass_size_ice_b
	      a_mice = mass_size_ice_a
	      b_as_ice = area_size_ice_b !1.63d0 !from mitchell 1996 similar to a_msnow&b_snow
	      a_as_ice = area_size_ice_a !0.11d0 * 10**(2*b_as_ice-4)
	    end if
    else
      print*, "did not understand SD_ice: ", SD_ice
      stop
    end if
  else if (n_moments .eq. 2) then
     if (nc .eq. 0.d0) stop 'STOP in routine ice_ssp'
     call double_moments(iwc,nc,gamma_ice(1),gamma_ice(2),gamma_ice(3),gamma_ice(4), &
          ad,bd,alpha,gamma,a_mice,b_mice)
     nbins = 100
     den_ice=917.d0
     dia1 = 1.d-10	! minimum diameter [m]
     dia2 = 6.d-5	! maximum diameter [m]
     !to_do: implement area-size relation in two-moments!
     !area-size relation in SI
     a_as_ice = 0.684 !0.684 also in CGS
     b_as_ice = 2.d0 !from mitchell 1996 similar to a_msnow&b_snow
  else
     stop 'Number of moments is not specified'
  end if

  if ((EM_ice .eq. 'mieic')) then
    nbins_spec = nbins+1 !Mie routine uses nbins+1 bins!
  else
    nbins_spec = nbins
  end if
  allocate(diameter_spec(nbins_spec),stat=alloc_status)
  allocate(back_spec(nbins_spec),stat=alloc_status)

  if (EM_ice .eq. 'mieic') then
     call mie(f, mindex,      &
          dia1, dia2, nbins, maxleg,   &
          ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
          back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
          LEGEN4, SD_ice,den_ice,iwc,&
          diameter_spec, back_spec)
      scatter_matrix= 0.d0
      extinct_matrix= 0.d0
      emis_vector= 0.d0
  elseif (EM_ice .eq. 'densi' .or. EM_ice .eq. 'surus') then
   if (EM_ice .eq. 'surus') den_ice = 0.863*f+115.d0
     call mie_densitydep_spheremasseq(f, t,mindex,      &
          a_mice, b_mice, dia1, dia2, nbins, maxleg,   &
          ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
          back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
          LEGEN4, SD_ice,den_ice,iwc,&
          diameter_spec, back_spec)
      scatter_matrix= 0.d0
      extinct_matrix= 0.d0
      emis_vector= 0.d0

  elseif (EM_ice .eq. 'liudb') then
     call dda_db_liu(f,t,liu_type_ice,mindex, &
          dia1,dia2,nbins,maxleg,ad,&
          bd, alpha, gamma, lphase_flag,kext, salb,&
          back, nlegen, legen, legen2, legen3,&
          legen4, SD_ice,&
          diameter_spec, back_spec)
      scatter_matrix= 0.d0
      extinct_matrix= 0.d0
      emis_vector= 0.d0
  elseif (EM_ice .eq. 'tmatr') then


    call tmatrix_ice(f, iwc, t, nc, &
          ad, bd, alpha, gamma, a_mice, b_mice, SD_ice, dia1, dia2, nbins, &
          scatter_matrix,extinct_matrix, emis_vector,&
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

  
    elseif (EM_ice .eq. 'tmSQL') then

    call tmatrix_snow_sql(f, iwc, t, nc, &
          ad, bd, alpha, gamma, a_mice, b_mice, as_ratio_ice, SD_ice, dia1, dia2, nbins, &
          scatter_matrix,extinct_matrix, emis_vector,&
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

  
     write (*, *) 'no em mod ', EM_ice
     stop
  endif

  if ((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments"))) then
    particle_type ="ice"

    call radar_spectrum(nbins_spec,diameter_spec, back,  back_spec,t,press,f,&
      particle_type,a_mice,b_mice,a_as_ice,b_as_ice,ice_spec)

  else
    ice_spec(:)=0.d0
  end if
  deallocate(diameter_spec, back_spec)

  if (verbose .gt. 1) print*, 'Exiting ice_ssp'

  return

end subroutine ice_ssp
