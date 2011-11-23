subroutine ref_ice(T,f,refre,refim)!,absind,abscoef)

  ! This function calculates the components (Re,Im) of the complex index of refraction 
  ! of ice. No explicit freq. dependence for Re(RI), const. value for T < 240.
  !
  ! Input:
  !  T  temperature in Kelvin
  !	 f  frequency in GHz
  !
  ! Output:
  !	 refre    real part of complex index of refraction n_r
  !	 refim    imaginary part of complex index of refraction n_i
  !  absind   absorption index
  !  abscoef  absorption coefficient
  !
  ! References:
  !      Maetzler 2006: Thermal microwave radiation: Application for remote sensing
  !      Personal communication with Maetzler

  use kinds
  use constants, only: pi, c

  implicit none

  real(kind=dbl), intent(in) :: T, &
       f

  real(kind=dbl), intent(out) :: refre,& ! real part of complex index of refraction n_r []
       refim!,& ! imaginary part of complex index of refraction n_i []
!       absind,&! absind absorptive index (n_i/n_r) []
!       abscoef ! abscoef absorption coefficient (4*pi*n_i*f/c) [1/m]

  complex(kind=dbl) :: epsi, eps_ice, ref_i

  epsi =  eps_ice(T,f)

  ref_i = sqrt(epsi)
  refre = real(ref_i)
  refim = aimag(ref_i)

!  absind = refim/refre
!  abscoef = (4*pi*refim*f*1.e9/c)

  return

end subroutine ref_ice
