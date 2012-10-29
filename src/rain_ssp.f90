subroutine rain_ssp(f,rwc,cwc,t,maxleg,kext, salb, back,  &
     nlegen, legen, legen2, legen3, legen4, nc)

  use kinds
  use nml_params, only: verbose, lphase_flag, n_0rainD, SD_rain, n_moments
  use constants, only: pi, im
  use double_moments_module
  use conversions

  implicit none

  integer :: nbins, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
       rwc,&
       cwc,&
       t,&
       f

  real(kind=dbl), optional, intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: absind, abscof

  real(kind=dbl) :: dia1, dia2, den_liq, ad, bd, alpha, gamma, b_rain, a_mrain

  real(kind=dbl), intent(out) :: &
       kext,&
       salb,&
       back

  real(kind=dbl), dimension(200):: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex

  real(kind=dbl) :: gammln

  if (verbose .gt. 1) print*, 'Entering rain_ssp'

  call ref_water(0.d0, t-273.15, f, refre, refim, absind, abscof)
  mindex = refre-im*refim

  den_liq = 1.d3  ! density of liquid water [kg/m^3]

!  nbins = 100
  nbins = 50
!  dia1 = 1.d-10 ! minimum diameter [m]
!  dia2 = 6.d-3  ! maximum diameter [m]
  dia1 = 1.2d-4 ! minimum diameter [m]
  dia2 = 6.d-3  ! maximum diameter [m]

  if (n_moments .eq. 1) then

	if (SD_rain .eq. 'C') then

	  ! this is for integration over diameters

	  ad = n_0rainD*1.d6   ! [1/m^4]
   	  bd = (pi * den_liq * ad / rwc)**0.25
      alpha = 0.d0 ! exponential SD
      gamma = 1.d0
    else if (SD_rain .eq. 'M') then
      b_rain = 3.
      a_mrain = 524.
	  bd = (rwc/(a_mrain*1.d7*exp(gammln(b_rain+1.))))**(1./(-1.-b_rain))
	  ad = 1.d7
	  dia2 = log(ad)/bd
	  if (dia2 .gt. 6.d-3) dia2 = 6.d-3
      alpha = 0.d0 ! exponential SD
      gamma = 1.d0
    end if
  else if (n_moments .eq. 2) then
     if (.not. present(nc)) stop 'STOP in routine rain_ssp'
     call double_moments(rwc,nc,gamma_rain(1),gamma_rain(2),gamma_rain(3),gamma_rain(4), &
          ad,bd,alpha,gamma,a_mrain, b_rain,cwc)
  else
     stop 'Number of moments is not specified'
  end if
  call mie(f, mindex, dia1, dia2, nbins, maxleg, ad,    &
       bd, alpha, gamma, lphase_flag, kext, salb, back,     &
       nlegen, legen, legen2, legen3, legen4, SD_rain,den_liq,rwc)

  if (verbose .gt. 1) print*, 'Exiting rain_ssp'

  return

end subroutine rain_ssp
