subroutine ref_water(s,T,f,refre,refim,absind,abscoef)

! This function calculates the components (Re,Im) of the complex index of refraction 
! for natural water including soluted salt.
! Valid parameter range: s 0-40, t 0-30, f 0-500 (1000)
!
! Input:
!	s  salinity in ppt
!	t  temperature in Celsius degree
!	f  frequency in GHz
!
! Output:
!	refre real part of complex index of refraction n_r
!	refim imaginary part of complex index of refraction n_i
!	absind absorptive index (n_i/n_r)
!	abscoef absorption coefficient (4*pi*n_i*f/c) [1/m]
! 
! References:
!      Maetzler 2006: Thermal microwave radiation: Application for remote sensing

  use kinds 

  implicit none

  real(kind=dbl), intent(in) :: s,& ! salinity [0/00]
				T,& ! temperature [°C]
				f   ! frequency [GHz]

  real(kind=dbl), intent(out) :: refre,& ! real part of complex index of refraction n_r []
				 refim,& ! imaginary part of complex index of refraction n_i []
				 absind,&! absind absorptive index (n_i/n_r) []
				 abscoef ! abscoef absorption coefficient (4*pi*n_i*f/c) [1/m]


  complex(kind=dbl) :: eps_water, epsi, ref_wat

  real(kind=dbl), parameter :: pi = 3.141592653589793, &
			       c = 299792458
  !complex permittivity of natural water
  epsi =  eps_water(s,T,f)

  ref_wat = sqrt(epsi)
  refre = real(ref_wat)
  refim = aimag(ref_wat)

  absind = refim/refre
  abscoef = (4*pi*refim*f*1.e9/c)

  return

end subroutine ref_water