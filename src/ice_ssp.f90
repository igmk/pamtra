subroutine ice_ssp(f,iwc,t,press,hgt,maxleg,nc, kext, salb, back,  &
     nlegen, legen, legen2, legen3, legen4,&
     scatter_matrix,extinct_matrix, emis_vector,ice_spec)

  use kinds
  use nml_params, only: verbose, lphase_flag, n_moments, EM_ice, SD_ice,&
      nstokes, radar_nfft, radar_spectrum
  use constants, only: pi, im
  use double_moments_module
  use conversions

  implicit none

  integer :: nbins, nbins_spec, nlegen,alloc_status

  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
       iwc,&
       t,&
       f,press,hgt

  real(kind=dbl), intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, del_d, den_ice, drop_mass, b_ice, a_mice
  real(kind=dbl) :: ad, bd, alpha, gamma, number_concentration

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
  real(kind=dbl), intent(out), dimension(radar_nfft) :: ice_spec

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
     alpha = 0.0d0     ! exponential SD
     gamma = 1.0d0
     den_ice=917.d0
    else if (SD_ice .eq. 'M') then
     del_d = 1.d-8    ! [m]
     dia1 = 6.d-5     ! [m] 60 micron diameter
     dia2 = dia1 + del_d
     den_ice = 917.d0  ! [kg/m^3]
     drop_mass = pi/6.d0 * dia1**3 * den_ice
     a_mice = 0.82d0
     b_ice = 2.5d0
     ad = iwc/(drop_mass*del_d) 	!intercept parameter [1/m^4]
     bd = 0.0d0
     nbins = 2
     alpha = 0.0d0     ! exponential SD
     gamma = 1.0d0 
    end if
  else if (n_moments .eq. 2) then
     if (nc .eq. 0.d0) stop 'STOP in routine ice_ssp'
     call double_moments(iwc,nc,gamma_ice(1),gamma_ice(2),gamma_ice(3),gamma_ice(4), &
          ad,bd,alpha,gamma,a_mice,b_ice)
     nbins = 100
     den_ice=917.d0
     dia1 = 1.d-10	! minimum diameter [m]
     dia2 = 6.d-5	! maximum diameter [m]
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
  elseif (EM_ice .eq. 'liudb') then
     call dda_db_liu(f,t,9,mindex, &
          dia1,dia2,nbins,maxleg,ad,&
          bd, alpha, gamma, lphase_flag,kext, salb,&
          back, nlegen, legen, legen2, legen3,&
          legen4, SD_ice,&
          diameter_spec, back_spec)
      scatter_matrix= 0.d0
      extinct_matrix= 0.d0
      emis_vector= 0.d0
  else
     write (*, *) 'no em mod', EM_ice
     stop
  endif

  if (radar_spectrum) then
    particle_type ="ice"

    call calc_radar_spectrum(nbins_spec,diameter_spec, back_spec,t,press,hgt,f,particle_type,ice_spec)
  else
    ice_spec(:)=0.d0
  end if
  deallocate(diameter_spec, back_spec)

  if (verbose .gt. 1) print*, 'Exiting ice_ssp'

  return

end subroutine ice_ssp
