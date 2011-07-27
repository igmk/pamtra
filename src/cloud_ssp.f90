subroutine cloud_ssp(f,qc,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4, nc)
  
  ! This subroutine prepares the input parameters for the routines
  ! calculating the extinction, absorption, and backscattering coefficients.
  !
  ! Input:
  !   f    frequency [GHz]
  !   qc   claud water mass mixing ratio [kg/kg]
  !   t    temperature [K]
  !   p    pressure [Pa]
  !   q    specific humidity (mass of water vapor per moist air) [kg/kg]
  !
  ! Output:
  !  kext
  !  salb
  !  back
  !  legen[2-4] legendre coefficients for the phase function

  use kinds
  use nml_params, only: verbose, lphase_flag, n_moments
  use constants, only: pi, im
  use double_moments_module

  implicit none

  integer :: nbins, nlegen, iautocon

  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qc,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl), optional, intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: absind, abscof

  real(kind=dbl) :: dia1, dia2, del_d, den_liq, drop_mass, b_cloud, a_mcloud
  real(kind=dbl) :: lwc, ad, bd, alpha, gamma, number_density, nc_abs

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex

  real(kind=dbl) :: spec2abs

  character(1) :: dist_name

  if (verbose .gt. 1) print*, 'Entering cloud_ssp'

  ! absind  absorptive index 
  ! abscof  absorption coefficient    [1/m]

  call ref_water(0.d0, t-273.15, f, refre, refim, absind, abscof)
  mindex = refre-im*refim

  lwc = spec2abs(qc,t,p,q) 									! [kg/m^3]

  if (n_moments .eq. 1) then
  ! iautocon is defined in COSMO-de routine hydci_pp_gr in src_gscp.f90
  !		iautocon =  1 fixed number density eq to 1.0d8 [1/m^3] (eq to variable cloud_num in COSMO)
  !		iautocon /= 1 fixed diameter eq to 1.d-4 [m]
      iautocon = 0
      den_liq = 1.d3 	  									! density of liquid water [kg/m^3]
      del_d = 1.d-8     									! delta_diameter for mie calculation [m]
    if (iautocon .eq. 1) then
      number_density = 1.0d8								! fixed number density [1/m^3]
  	  drop_mass = lwc /  number_density						! [kg]
  	  dia1 = (6.d0 * drop_mass / (pi * den_liq))**(1./3.)	! monodisperse size distribution [m]
  	  dia2 = dia1 + del_d
  	  ad = lwc / (drop_mass*del_d)							! intercept parameter [1/m^4]
    else
	  dia1 = 2.d-5      									! [m] 20 micron diameter monodisperse
 	  dia2 = dia1 + del_d
  	  drop_mass = pi/6.d0 * den_liq * dia1**3 		    	! [kg]
  	  ad = lwc / (drop_mass*del_d)							! intercept parameter [1/m^4]
    end if
	  bd = 0.d0
      nbins = 2
      alpha = 0.d0 ! exponential SD
      gamma = 1.d0
      dist_name='C'
  else if (n_moments .eq. 2) then
    if (.not. present(nc)) stop 'STOP in routine cloud_ssp'
    nc_abs = spec2abs(nc,t,p,q) 							! [#/m^3]
    call double_moments(lwc,nc_abs,gamma_cloud(1),gamma_cloud(2),gamma_cloud(3),gamma_cloud(4), &
    	ad,bd,alpha,gamma,a_mcloud,b_cloud)
    nbins = 100
    dia1 = .5d-6	! minimum diameter [m]
    dia2 = 1.d-4	! maximum diameter [m]
    dist_name='G'
  else
    stop 'Number of moments is not specified'
  end if


  call mie(f, mindex, dia1, dia2, nbins, maxleg, ad,       &
	bd, alpha, gamma, lphase_flag, kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4, dist_name)

  if (verbose .gt. 1) print*, 'Exiting cloud_ssp'

  return

end subroutine cloud_ssp
