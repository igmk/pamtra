module double_moments_module

  use kinds

  implicit none

!Parameter of the drop size distribution when 2moments scheme is used
real (kind=dbl), dimension(4) :: 	gamma_cloud,&
  									gamma_rain, &
  									gamma_ice, &
  									gamma_snow, &
  									gamma_graupel, &
  									gamma_hail

end module double_moments_module
