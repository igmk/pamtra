subroutine ref_ice(T,freq,refre,refim)

! This function calculates the components (Re,Im) of the complex index of refraction 
! of ice. No explicit freq. dependence for Re(RI), const. value for T < 240.
!
! Input:
!	T  temperature in Kelvin
!	f  frequency in GHz
!
! Output:
!	refre real part of complex index of refraction n_r
!	refim imaginary part of complex index of refraction n_i
! 
! References:
!      Maetzler 2006: Thermal microwave radiation: Application for remote sensing
!      Personal communication with Maetzler

  use kinds

  implicit none

  real(kind=dbl), intent(in) :: T, &
				freq

  real(kind=dbl), intent(out) :: refre,& ! real part of complex index of refraction n_r []
				 refim   ! imaginary part of complex index of refraction n_i []

  real(kind=dbl) :: eps_real,  &
		    mit,       &
		    alpha,     &
		    beta,     &
		    beta_m,    &
		    delta_beta,&
		    eps_imag

  complex(kind=dbl) :: eps

  if (T .ge. 240.) eps_real = 3.1884 + 9.1d-4*(T-273.)
  if (T .lt. 240.) eps_real = 3.1884 + 9.1d-4*(240. - 273.)

  ! "modified inverse temperature"
  mit = (300./T)-1

  alpha = (0.00504 + 0.0062*mit)*exp(-22.1*mit)

  beta_m = (0.0207/T)*(exp(335./T))/(exp(335./T)-1.)**2.
  beta_m = beta_m + 1.16d-11*(freq**2.)
  delta_beta = exp(-9.963 + 0.0372*(T-273.16))

  beta = beta_m + delta_beta
  eps_imag = (alpha/freq)+beta*freq

  eps = cmplx(eps_real, eps_imag)

  refre = real(sqrt(eps))
  refim = aimag(sqrt(eps))

  return

end subroutine ref_ice
