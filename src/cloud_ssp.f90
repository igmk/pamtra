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

  integer :: nbins, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qc,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: absind, abscof

  real(kind=dbl) :: dia1, dia2, del_d, den_liq, drop_mass, lwc, ad, bd, alpha, gamma

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
  del_d = 1.d-8     ! [m]
  dia1 = 2.d-5      ! [m] 20 micron diameter monodisperse
  dia2 = dia1 + del_d
  den_liq = 1.d3 ! density of liquid water [kg/m^3]
  drop_mass = 1./6. * pi * dia1**3 * den_liq ! [kg]
  lwc = spec2abs(qc,t,p,q) ! [kg/m^3]
  ad = lwc / (drop_mass*del_d)

  bd = 0.d0 
  alpha = 0.d0 ! exponential SD
  gamma = 1.d0 
  nbins = 2

  call mie(f, mindex, dia1, dia2, nbins, maxleg, ad,       &
	bd, alpha, gamma, lphase_flag, kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4, 'C')

  if (verbose .gt. 1) print*, 'Exiting cloud_ssp'

  return

end subroutine cloud_ssp
