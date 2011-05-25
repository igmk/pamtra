subroutine cloud_ssp(f,qc,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4)
  
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
  use nml_params, only: verbose, lphase_flag
  use constants, only: pi, im

  implicit none

  integer :: numrad, nlegen, iautocon
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qc,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: absind, abscof

  real(kind=dbl) :: rad1, rad2, del_r, den_liq, drop_mass, lwc, ad, bd, alpha, gamma, number_density

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex

  real(kind=dbl) :: spec2abs

  if (verbose .gt. 1) print*, 'Entering cloud_ssp'

  ! absind  absorptive index 
  ! abscof  absorption coefficient    [1/m]

  call ref_water(0.d0, t-273.15, f, refre, refim, absind, abscof)
  mindex = refre-im*refim
  ! iautocon is defined in COSMO-de routine hydci_pp_gr in src_gscp.f90
  !		iautocon =  1 fixed number density eq to 1.0d8 [1/m^3] (eq to variable cloud_num in COSMO)
  !		iautocon /= 1 fixed radius eq to 1.d-5 [m]
  iautocon = 1

  den_liq = 1.d3 	  									! density of liquid water [kg/m^3]
  del_r = 1.d-8     									! delta_radius for mie calculation [m]
 if (iautocon .eq. 1) then
	number_density = 1.0d8								! fixed number density [1/m^3]
  	lwc = spec2abs(qc,t,p,q) 							! [kg/m^3]
  	drop_mass = lwc /  number_density					! [kg]
  	rad1 = (3./4. * drop_mass / pi / den_liq)**(1./3.)	! monodisperse size distribution [m]
  	rad2 = rad1 + del_r
  	ad = lwc / (drop_mass*del_r)						! intercept parameter [1/m^4]
!write(45,'(e13.6,x,e13.6,x,e13.6)')lwc,rad1,number_density
else
	rad1 = 1.d-5      									! [m] 10 micron radius monodisperse
 	rad2 = rad1 + del_r
  	drop_mass = 4./3. * pi * rad1**3 * den_liq 			! [kg]
!number_density=lwc / drop_mass
  	lwc = spec2abs(qc,t,p,q) 							! [kg/m^3]
  	ad = lwc / (drop_mass*del_r)						! intercept parameter [1/m^4]
!write(46,'(e13.6,x,e13.6,x,e13.6)')lwc,rad1,number_density
end if

  bd = 0.d0 
  alpha = 0.d0 ! exponential SD
  gamma = 1.d0 
  numrad = 2 

  call mie(f, mindex, rad1, rad2, numrad, maxleg, ad,       &
	bd, alpha, gamma, lphase_flag, kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4, 'C')

  if (verbose .gt. 1) print*, 'Exiting cloud_ssp'

  return

end subroutine cloud_ssp
