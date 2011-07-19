subroutine ice_ssp(f,qi,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4)
  
  use kinds
  use nml_params, only: verbose, lphase_flag
  use constants, only: pi, im

  implicit none

  integer :: nbins, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qi,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, del_d, den_ice, drop_mass, iwc, ad, bd, alpha, gamma

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex

  real(kind=dbl) :: spec2abs

  if (verbose .gt. 1) print*, 'Entering ice_ssp'

    call ref_ice(t,f, refre, refim)
    mindex = refre-im*refim  ! mimicking a
    ! monodisperse distribution
    del_d = 1.d-8    ! [m]
    dia1 = 1.d-4     ! [m] 10 micron diameter
    dia2 = dia1 + del_d
    den_ice = 917.d0  ! [kg/m^3]
    drop_mass = pi/6. * dia1**3 * den_ice
    iwc = spec2abs(qi,t,p,q) ! [kg/m^3]

    ad = iwc/(drop_mass*del_d)
    bd = 0.0d0 
    alpha = 0.0d0     ! exponential SD
    gamma = 1.0d0 
    nbins = 2

    call mie(f, mindex, dia1, dia2, nbins, maxleg, ad,    &
	  bd, alpha, gamma, lphase_flag, kext, salb, back,     &
	  nlegen, legen, legen2, legen3, legen4, 'C')
  if (verbose .gt. 1) print*, 'Exiting ice_ssp'

  return

end subroutine ice_ssp
