subroutine ice_ssp(f,iwc,t,maxleg,kext, salb, back,  &
     nlegen, legen, legen2, legen3, legen4, nc)

  use kinds
  use nml_params, only: verbose, lphase_flag, n_moments, EM_ice, SD_ice
  use constants, only: pi, im
  use double_moments_module
  use conversions

  implicit none

  integer :: nbins, nlegen, ice_flag

  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
       iwc,&
       t,&
       f

  real(kind=dbl), optional, intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, del_d, den_ice, drop_mass, b_ice, a_mice
  real(kind=dbl) :: ad, bd, alpha, gamma, number_concentration

  real(kind=dbl), intent(out) :: &
       kext,&
       salb,&
       back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex

  if (verbose .gt. 1) print*, 'Entering ice_ssp'

  call ref_ice(t,f, refre, refim)
  mindex = refre-im*refim  ! mimicking a

  if (n_moments .eq. 1) then
	if (SD_ice .eq. 'C') then
     ice_flag=0
     if (ice_flag .eq. 1) then
        !	Monodisperse size distribution coherent with COSMO-de 1-moment scheme
        !	Number_concentration of activated ice ctystal is temperature dependent
        !	 from COSMO-de code src_gscp.f90 routine: hydci_pp_gr
        !	Radius is derived from mass-size relation m=aD^3
        !	 a=130 kg/m^3 (hexagonal plates with aspect ratio of 0.2 -> thickness=0.2*Diameter)
        number_concentration = 1.0d2*DEXP(0.2d0*(273.15d0-t)) 	! [1/m^3]
        drop_mass = iwc/number_concentration 					! [kg]
        del_d = 1.d-8											! [m]
        dia1 = (drop_mass/130.0d0)**(1.0d0/3.0d0)				! [m]
        dia2 = dia1 + del_d
     else
        ! monodisperse distribution
        ! Fixed diameter
        del_d = 1.d-8    ! [m]
        dia1 = 1.d-4     ! [m] 100 micron diameter
        dia2 = dia1 + del_d
        den_ice = 917.d0  ! [kg/m^3]
        drop_mass = pi/6.d0 * dia1**3 * den_ice
     endif
     ad = iwc/(drop_mass*del_d) 	!intercept parameter [1/m^4]
     bd = 0.0d0
     nbins = 2
     alpha = 0.0d0     ! exponential SD
     gamma = 1.0d0
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
     if (.not. present(nc)) stop 'STOP in routine ice_ssp'
     call double_moments(iwc,nc,gamma_ice(1),gamma_ice(2),gamma_ice(3),gamma_ice(4), &
          ad,bd,alpha,gamma,a_mice,b_ice)
     nbins = 100
     dia1 = 1.d-6	! minimum diameter [m]
     dia2 = 1.d-3	! maximum diameter [m]
  else
     stop 'Number of moments is not specified'
  end if

  if (EM_ice .eq. 'mieic') then
     call mie(f, mindex,      &
          dia1, dia2, nbins, maxleg,   &
          ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
          back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
          LEGEN4, SD_ice,den_ice,iwc)
  elseif (EM_ice .eq. 'liudb') then
     call dda_db_liu(f,t,9,mindex, &
          dia1,dia2,nbins,maxleg,ad,&
          bd, alpha, gamma, lphase_flag,kext, salb,&
          back, nlegen, legen, legen2, legen3,&
          legen4, SD_ice)
  else
     write (*, *) 'no em mod', EM_ice
     stop
  endif
  if (verbose .gt. 1) print*, 'Exiting ice_ssp'

  return

end subroutine ice_ssp
