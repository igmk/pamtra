subroutine get_atmosG0
  !
  !  Calculate layer averages of pressure, temperature, relative
  !  humidity, vapor pressure, absolute humidity, and  and vapor pressure in specified    
  !  layers, given the temperature, pressure, and specific humidity.    
  !
  !  Output:
  !	atmo_temp        ! temperature [K]
  !	atmo_relhum      ! rel. hum. [%]
  !	atmo_press       ! pressure [Pa]
  !	atmo_vapor_pressure    ! vapor pressure [Pa]
  !	atmo_rho_vap     ! vapor density/absolute humidity [kg/m^3]
  !	atmo_q_hum       ! specific humidity [kg/kg]

  use kinds
  use vars_atmosphere, only: &
    atmo_relhum_lev, &
    atmo_press_lev, &
    atmo_hgt_lev, &
    atmo_temp_lev, &
    atmo_q_hum, &
    atmo_temp, &
    atmo_relhum, &
    atmo_press, &
    atmo_hgt, &
    atmo_delta_hgt_lev, &
    atmo_vapor_pressure, &
    atmo_rho_vap, &
    atmo_q_hum, &
    atmo_nlyrs

  use vars_index, only: i_x, i_y

  implicit none

  integer :: nz

  real(kind=dbl) :: e_sat_gg_water

  do nz = 1, atmo_nlyrs(i_x,i_y) 
     atmo_temp(i_x,i_y,nz) = 0.5 * (atmo_temp_lev(i_x,i_y,nz) + atmo_temp_lev(i_x,i_y,nz+1)) 
     atmo_relhum(i_x,i_y,nz) = 0.5 * (atmo_relhum_lev(i_x,i_y,nz) + atmo_relhum_lev(i_x,i_y,nz+1))
     atmo_press(i_x,i_y,nz) = (atmo_press_lev(i_x,i_y,nz+1) - atmo_press_lev(i_x,i_y,nz))/ &
          log(atmo_press_lev(i_x,i_y,nz+1) / atmo_press_lev(i_x,i_y,nz))  ! [Pa]
     atmo_hgt(i_x,i_y,nz) = (atmo_hgt_lev(i_x,i_y,nz)+atmo_hgt_lev(i_x,i_y,nz+1))*0.5d0
     atmo_delta_hgt_lev(i_x,i_y,nz) = atmo_hgt_lev(i_x,i_y,nz+1) - atmo_hgt_lev(i_x,i_y,nz)
     atmo_vapor_pressure(i_x,i_y,nz) = atmo_relhum(i_x,i_y,nz) * &
          e_sat_gg_water(atmo_temp(i_x,i_y,nz)) ! Pa
     atmo_rho_vap(i_x,i_y,nz) = atmo_vapor_pressure(i_x,i_y,nz)/(atmo_temp(i_x,i_y,nz) * 461.5)  ! [kg/m3]
     atmo_q_hum(i_x,i_y,nz) = 0.622*atmo_vapor_pressure(i_x,i_y,nz)/&
          (atmo_press(i_x,i_y,nz) - (1.- 0.622) * atmo_vapor_pressure(i_x,i_y,nz))  ! [kg/kg]
  end do

  return

end subroutine get_atmosG0
