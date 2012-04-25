function eps_mix(eps0,eps1,mix_de)

  ! calculates the permittivity for a mixture of two materials
  ! with permittivities eps0 and eps1, respectively.
  !
  ! Input:
  !   eps0
  !   eps1
  !   mix_de
  !
  ! Result:
  !   eps_mix
  !
  ! References:
  !    A. H. Sihvola, 1989: "Self-consistency aspects of dielectric mixing theories.",
  !          IEEE Trans. on Geosci. and Remote Sens., 27, 403-415.
  use kinds

  implicit none

  real(kind=dbl) :: a, b, v, f

  real(kind=dbl) :: em_real, em_imag

  real(kind=dbl), intent(in) :: mix_de

  complex(kind=dbl), intent(in) :: eps0, eps1

  complex(kind=dbl) :: eps_mix

  v = 0.85d0

  f = mix_de/917.d0

  b = real(eps1) + 2.d0*eps0 - 2.d0*v*eps0 - f*(real(eps1) - eps0)*(1 + v)

  a = eps0*(real(eps1) + (2.d0 - v)*eps0 + f*(2.d0 - v)*(real(eps1 - eps0)))

  em_real = (sqrt(b**2+4.d0*v*a) - b)/(2.d0*v)

  em_imag = aimag(eps1)*(-1.d0*(em_real - eps0) + f*(em_real + 2.d0*eps0 + v*(em_real - eps0)))/&
  			((real(eps1) + 2.d0*eps0 + 2.d0*v*(em_real - eps0)) - f*(1.d0 + v)*(real(eps1) - eps0))

  eps_mix = cmplx(em_real, em_imag)

  return

end function eps_mix
