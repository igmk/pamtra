subroutine get_atmosG0
  !
  !  Calculate layer averages of pressure, temperature, relative
  !  humidity, vapor pressure, absolute humidity, and  and vapor pressure in specified    
  !  layers, given the temperature, pressure, and specific humidity.    
  !
  !  Output:
  !	temp        ! temperature [K]
  !	relhum      ! rel. hum. [%]
  !	press       ! pressure [Pa]
  !	vapor_pressure    ! vapor pressure [Pa]
  !	rho_vap     ! vapor density/absolute humidity [kg/m^3]
  !	q_hum       ! specific humidity [kg/kg]

  use kinds
  use vars_atmosphere

  implicit none

  integer :: nz

  real(kind=dbl) :: e_sat_gg_water

  do nz = 1, nlyr 
     temp(nz) = 0.5 * (temp_lev(nz - 1) + temp_lev(nz)) 
     relhum(nz) = 0.5 * (relhum_lev(nz - 1) + relhum_lev(nz))
     press(nz) = (press_lev(nz) - press_lev(nz - 1))/log(press_lev(nz) / press_lev (nz - 1))  ! [Pa]
     hgt(nz) = (hgt_lev(nz-1)+hgt_lev(nz))*0.5d0
     delta_hgt_lev(nz) = hgt_lev(nz) - hgt_lev(nz-1)
     vapor_pressure(nz) = relhum(nz) * e_sat_gg_water(temp(nz)) ! Pa
     rho_vap(nz) = vapor_pressure(nz)/(temp(nz) * 461.5)  ! [kg/m3]
     q_hum(nz) = 0.622*vapor_pressure(nz)/(press(nz) - (1.- 0.622) * vapor_pressure(nz))  ! [kg/kg]
  end do

  return

end subroutine get_atmosG0
