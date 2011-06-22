subroutine rain_ssp(f,qr,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4, nr)
  
  use kinds
  use nml_params, only: verbose, lphase_flag, n_0rainD, SD_rain, n_moments
  use constants, only: pi, im
  use double_moments_module

  implicit none

  integer :: numrad, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qr,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl), optional, intent(in) :: nr

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: absind, abscof

  real(kind=dbl) :: rad1, rad2, den_liq, rwc, ad, bd, alpha, gamma, nr_abs, b_rain, a_mrain

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex

  real(kind=dbl) :: spec2abs

  character(1) :: dist_name

  if (verbose .gt. 1) print*, 'Entering rain_ssp'

    call ref_water(0.d0, t-273.15, f, refre, refim, absind, abscof)
    mindex = refre-im*refim
    rwc = spec2abs(qr,t,p,q) ! [kg/m^3]

  if (n_moments .eq. 1) then
    rad1 = 1.d-4   ! minimum diameter [m]
    rad2 = 6.d-3   ! maximum diameter [m]
    den_liq = 1.d3  ! density of liquid water [kg/m^3]
    ! this is for integration over diameters

    ad = n_0rainD*1.d6   ! [1/m^4]
    bd = (pi * den_liq * ad / rwc)**0.25

    numrad = 100 
    alpha = 0.d0 ! exponential SD
    gamma = 1.d0
    dist_name='C'
  else if (n_moments .eq. 2) then
    if (.not. present(nr)) stop 'STOP in routine rain_ssp'
    nr_abs = spec2abs(nr,t,p,q) 							! [#/m^3]
    call double_moments(rwc,nr_abs,gamma_rain(1),gamma_rain(2),gamma_rain(3),gamma_rain(4), &
    	ad,bd,alpha,gamma,a_mrain, b_rain)
    numrad = 100
    rad1 = 1.d-4	! minimum diameter [m]
    rad2 = 6.d-3	! maximum diameter [m]
    dist_name='G'
  else
    stop'Number of moments is not specified'
  end if

    call mie(f, mindex, rad1/2., rad2/2., numrad, maxleg, ad,    &
	  bd, alpha, gamma, lphase_flag, kext, salb, back,     &
	  nlegen, legen, legen2, legen3, legen4, dist_name)

  if (verbose .gt. 1) print*, 'Exiting rain_ssp'

  return

end subroutine rain_ssp
