function dielec_water(s,T,f)

  ! This function calculates the dielectric factor of water based on the relative permettivity
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
  !      Paul L Smith, “Equivalent Radar Reflectivity Factors for Snow and Ice Particles,”
  !                Journal of Climate and Applied Meteorology 23, no. 8 (1984): 1258-1260.

  use kinds 
  use eps_water, only: get_eps_water
  use report_module

  implicit none

  real(kind=dbl), intent(in) :: s,& ! salinity [0/00]
       T,& ! temperature [°C]
       f   ! frequency [GHz]

  real(kind=dbl) :: dielec_water ! dielectric factor |K|^2


  complex(kind=dbl) :: epsi
!   integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=30) :: nameOfRoutine = 'dielec_water'


  !complex permittivity of natural water
  call get_eps_water(err,s,T,f, epsi)
  if (err /= 0) then
      msg = 'error in get_eps_water!'
      call report(err, msg, nameOfRoutine)
!       errorstatus = err
!       return
      stop
  end if 
  dielec_water = abs((epsi-1)/(epsi+2))**2

  return

end function dielec_water
