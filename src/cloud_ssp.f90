subroutine cloud_ssp(f,cwc,t, press,hgt,maxleg,nc, kext, salb, back,  &
     nlegen, legen, legen2, legen3, legen4,&
     scatter_matrix,extinct_matrix, emis_vector,cloud_spec)

  ! This subroutine prepares the input parameters for the routines
  ! calculating the extinction, absorption, and backscattering coefficients.
  !
  ! Input:
  !   f    frequency [GHz]
  !   cwc   claud water content [kg/m^3]
  !   t    temperature [K]
  !
  ! Output:
  !  kext
  !  salb
  !  back
  !  legen[2-4] legendre coefficients for the phase function

  use kinds
  use nml_params, only: verbose, lphase_flag, n_moments, SD_cloud, &
      nstokes, EM_cloud, radar_nfft, radar_spectrum
  use constants, only: pi, im
  use double_moments_module
  use conversions

  implicit none

  integer :: nbins, nlegen, iautocon,alloc_status

  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
       cwc,&
       t, &
       f, press,hgt

  real(kind=dbl), intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: absind, abscof

  real(kind=dbl) :: dia1, dia2, del_d, den_liq, drop_mass, b_cloud, a_mcloud
  real(kind=dbl) :: ad, bd, alpha, gamma, number_density

  real(kind=dbl), intent(out) :: &
       kext,&
       salb,&
       back
  real(kind=dbl), intent(out), dimension(radar_nfft) :: cloud_spec

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4
    integer, parameter ::  nquad = 16
    real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
    real(kind=dbl), dimension(nstokes,nstokes,nquad,2), intent(out) :: extinct_matrix
    real(kind=dbl), dimension(nstokes,nquad,2), intent(out) :: emis_vector
  complex(kind=dbl) :: mindex
  real(kind=dbl), allocatable, dimension(:):: diameter_spec, qback_spec
  character(5) ::  particle_type

real(kind=dbl) :: re, Nt
  if (verbose .gt. 1) print*, 'Entering cloud_ssp'

  ! absind  absorptive index 
  ! abscof  absorption coefficient    [1/m]

  call ref_water(0.d0, t-273.15, f, refre, refim, absind, abscof)
  mindex = refre-im*refim

  den_liq = 1.d3 	! density of liquid water [kg/m^3]

  if (n_moments .eq. 1) then
	if (SD_cloud .eq. 'C') then
	    del_d = 1.d-8 	! delta_diameter for mie calculation [m]
	    dia1 = 2.d-5 	! [m] 20 micron diameter monodisperse
	    dia2 = dia1 + del_d
	    drop_mass = pi/6.d0 * den_liq * dia1**3 		    	! [kg]
	    ad = cwc / (drop_mass*del_d)	! intercept parameter [1/m^4]
	    bd = 0.d0
	    nbins = 2
	    alpha = 0.d0 ! exponential SD
	    gamma = 1.d0
	else if (SD_cloud .eq. 'M') then
	    del_d = 1.d-8 	!delta_diameter for mie calculation [m]
	    dia1 = 2.d-5 	! [m] 20 micron diameter monodisperse
	    dia2 = dia1 + del_d
	    drop_mass = pi/6.d0 * den_liq * dia1**3 		    	! [kg]
	    ad = cwc / (drop_mass*del_d)	! intercept parameter [1/m^4]
	    bd = 0.d0
	    nbins = 2
	    alpha = 0.d0 ! exponential SD
	    gamma = 1.d0
	else if (SD_cloud .eq. 'L') then !lognormal cloud distribution to test Pavlos Radar Spectrum
!teh hard and ugly way:

	    dia1 = 2.d-6 ! [m] 2 micron diameter 
	    dia2 = 50.d-6 ! [m] 50 micron diameter 
	    nbins = 48
	    alpha = 0.3 !S_x in Pavlos code
	    re=10.d-6
	    bd = re/exp(2.5d0*alpha**2); !ro in Pavlos Code
	    Nt = (3*CWC*1.d3)/(4*pi*1*(bd*1d2)**3*exp(4.5d0*alpha**2))
	    ad = Nt/(sqrt(2*pi) *alpha)

        else
	    stop "did not understand SD_cloud"
	end if
  else if (n_moments .eq. 2) then
     if (nc .eq. 0.d0) stop 'STOP in routine cloud_ssp'
     call double_moments(cwc,nc,gamma_cloud(1),gamma_cloud(2),gamma_cloud(3),gamma_cloud(4), &
          ad,bd,alpha,gamma,a_mcloud,b_cloud)
     nbins = 100
     dia1 = .5d-10	! minimum diameter [m]
     dia2 = 1.d-4	! maximum diameter [m]
  else
     stop 'Number of moments is not specified'
  end if

  allocate(diameter_spec(nbins+1),stat=alloc_status)
  allocate(qback_spec(nbins+1),stat=alloc_status)

  if (EM_cloud .eq. 'miecl') then
  call mie(f, mindex, dia1, dia2, nbins, maxleg, ad,       &
       bd, alpha, gamma, lphase_flag, kext, salb, back,  &
       nlegen, legen, legen2, legen3, legen4, SD_cloud,den_liq,cwc,&
       diameter_spec, qback_spec)
      scatter_matrix= 0.d0
      extinct_matrix= 0.d0
      emis_vector= 0.d0
  else
    stop "unknown EM_cloud"
  end if


  particle_type="cloud" 

  if (radar_spectrum) then
    call calc_radar_spectrum(nbins+1,diameter_spec, qback_spec,t,press,hgt,f,particle_type,cloud_spec)
  else
    cloud_spec(:)=0.d0
  end if

  deallocate(diameter_spec, qback_spec)

  if (verbose .gt. 1) print*, 'Exiting cloud_ssp'
  return

end subroutine cloud_ssp
