subroutine rain_ssp(f,qr,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4)
  
  use kinds
  use nml_params, only: verbose, lphase_flag, n_0rainD, SD_rain
  use constants, only: pi, im

  implicit none

  integer :: nbins, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qr,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: absind, abscof

  real(kind=dbl) :: dia1, dia2, den_liq, rwc, ad, bd, alpha, gamma

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex

  real(kind=dbl) :: spec2abs

  if (verbose .gt. 1) print*, 'Entering rain_ssp'

    call ref_water(0.d0, t-273.15, f, refre, refim, absind, abscof)

    mindex = refre-im*refim
    dia1 = 1.d-4   ! minimum diameter [m]
    dia2 = 6.d-3   ! maximum diameter [m]
    den_liq = 1.d3  ! density of liquid water [kg/m^3]
    ! this is for integration over diameters
    rwc = spec2abs(qr,t,p,q) ! [kg/m^3]

    ad = n_0rainD*1.d6   ! [1/m^4]
    bd = (pi * den_liq * ad / rwc)**0.25

    nbins = 100
    alpha = 0.d0 ! exponential SD
    gamma = 1.d0 

    call mie(f, mindex, dia1, dia2, nbins, maxleg, ad,    &
	  bd, alpha, gamma, lphase_flag, kext, salb, back,     &
	  nlegen, legen, legen2, legen3, legen4, 'C')    

  if (verbose .gt. 1) print*, 'Exiting rain_ssp'

  return

end subroutine rain_ssp
