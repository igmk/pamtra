function eps_water(s,T,f)

! This function calculates the complex permittivity
! of natural water including soluted salt.
! Valid parameter range: s 0-40, T 0-30, f 0-500 (1000)
!
! Input:
!	s  salinity [ppt]
!	t  temperature [°C]
!	f  frequency [GHz]
!
! Result:
!	eps_water complex permittivity of natural water
! 
! References:
!      Maetzler 2006: Thermal microwave radiation: Application for remote sensing
!      section 5.2.5.4

  use kinds 

  implicit none

  real(kind=dbl), intent(in) :: s,& ! salinity [0/00]
				T,& ! temperature [°C]
				f   ! frequency [GHz]

  real(kind=dbl), dimension(18), parameter :: a = (/&
	0.46606917e-02, & ! a1
       -0.26087876e-04, & ! a2
       -0.63926782e-05, & ! a3
	0.63000075e+01, & ! a4
	0.26242021e-02, & ! a5
       -0.42984155e-02, & ! a6
	0.34414691e-04, & ! a7
	0.17667420e-03, & ! a8
       -0.20491560e-06, & ! a9
	0.58366888e+03, & ! a10
	0.12634992e+03, & ! a11
	0.69227972e-04, & ! a12
	0.38957681e-06, & ! a13
	0.30742330e+03, & ! a14
	0.12634992e+03, & ! a15
	0.37245044e+01, & ! a16
	0.92609781e-02, & ! a17
       -0.26093754e-01  & ! a18
	/)

  real(kind=dbl), parameter :: pi = 3.141592653589793

  real(kind=dbl) :: tau_1, tau_2, q, sig

  real(kind=dbl) :: p, alpha_0, alpha_1
  real(kind=dbl) :: sig_s35

  complex(kind=dbl) :: eps_s, eps_1, eps_inf

  complex(kind=dbl) :: eps_water

  p = s*(37.5109+5.45216*s+0.014409*s**2)/(1004.75+182.283*s+s**2)
  alpha_0 = (6.9431+3.2841*s-0.099486*s**2)/(84.85+69.024*s+s**2)
  alpha_1 = 49.843-0.2276*s+0.00198*s**2
  sig_s35 = 2.903602+8.607e-02*T+4.738817e-04*T**2-2.991e-06*T**3+4.3041e-09*T**4

  eps_s = 87.85306*exp(-0.00456992*T-a(1)*s-a(2)*s**2-a(3)*s*T)
  eps_1 = a(4)*exp(-a(5)*T-a(6)*s-a(7)*s*T)
  tau_1 = (a(8)+a(9)*s)*exp(a(10)/(T+a(11)))
  tau_2 = (a(12)+a(13)*s)*exp(a(14)/(T+a(15)))
  eps_inf = a(16)+a(17)*T+a(18)*s
  q = 1+(alpha_0*(T-15))/(T+alpha_1)
  sig = sig_s35*p*q

  eps_water = (eps_s-eps_1)/(1-(0.,1)*2*pi*f*tau_1)+&
    (eps_1-eps_inf)/(1-(0.,1.)*2*pi*f*tau_2)+eps_inf+(0.,1.)*17.9751*sig/f

  return

end function eps_water

!! refractive index
!  ref_wat = sqrt(eps_water)
!   refre = real(ref_wat)
!   refim = aimag(ref_wat)

  ! absind  absorptive index 
  ! abscof  absorption coefficient    [1/m]

!   absind = refim/refre
!   abscoef = (4*pi*refim*f*1.e9/c)
