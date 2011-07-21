subroutine dielec_water(s,T,f,K2)

! This function calculates the dielectric constant of warter based on the relative permettivity
! for natural water including soluted salt.
! Valid parameter range: s 0-40, t 0-30, f 0-500 (1000)
!
! Input:
!	s  salinity in ppt
!	t  temperature in Celsius degree
!	f  frequency in GHz
!
! Output:
!	K2 dielectric constant |K|^2
! 
! References:
!      Paul L Smith, “Equivalent Radar Reflectivity Factors for Snow and Ice Particles,” Journal of Climate and Applied Meteorology 23, no. 8 (1984): 1258-1260.


  use kinds 

  implicit none

  real(kind=dbl), intent(in) :: s,& ! salinity [0/00]
				T,& ! temperature [°C]
				f   ! frequency [GHz]

  real(kind=dbl), intent(out) :: K2 !dielectric constant |K|^2


  complex(kind=dbl) :: eps_water, epsi



  !complex permittivity of natural water
  epsi =  eps_water(s,T,f)

  K2 = abs((epsi-1)/(epsi+2))**2

  return

end subroutine dielec_water