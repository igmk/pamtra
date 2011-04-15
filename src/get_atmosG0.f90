subroutine get_atmosG0
!
!  Calculate layer averages of pressure, temperature, relative
!  humidity, vapor pressure, absolute humidity, and  and vapor pressure in specified    
!  layers, given the temperature, pressure, and specific humidity.    
!
!  Output:
!	p(nx,ny)%temp              ! temperature [K]
!	p(nx,ny)%relhum            ! rel. hum. [%]
!	p(nx,ny)%press             ! pressure [Pa]
!	p(nx,ny)%vapor_pressure    ! vapor pressure [Pa]
!	p(nx,ny)%rho_vap           ! vapor density/absolute humidity [kg/m^3]
!	p(nx,ny)%q_hum             ! specific humidity [kg/kg]

  use kinds
  use vars_atmosphere

  implicit none

  integer :: nx, ny, nz

  real(kind=dbl) :: e_sat_gg_water

								    
  do nx = 1, ngridx 
    do ny = 1, ngridy 
      do nz = 1, nlyr 
									
      profiles(nx,ny)%temp(nz) = 0.5 * (profiles(nx,ny)%temp_lev(nz - 1) + profiles(nx,ny)%temp_lev(nz)) 
      profiles(nx,ny)%relhum(nz) = 0.5 * (profiles(nx,ny)%relhum_lev(nz - 1) + profiles(nx,ny)%relhum_lev(nz))
      profiles(nx,ny)%press(nz) = (profiles(nx,ny)%press_lev(nz) - profiles(nx,ny)%press_lev(nz - 1))/&
      log(profiles(nx,ny)%press_lev(nz) / profiles(nx,ny)%press_lev (nz - 1))  ! [Pa]   
								      
      profiles(nx,ny)%vapor_pressure(nz) = profiles(nx,ny)%relhum(nz) * e_sat_gg_water(profiles(nx,ny)%temp(nz)) ! Pa
      profiles(nx,ny)%rho_vap(nz) = profiles(nx,ny)%vapor_pressure(nz)/(profiles(nx,ny)%temp(nz) * 461.5)        ! [kg/m3]
      profiles(nx,ny)%q_hum(nz) = 0.622*profiles(nx,ny)%vapor_pressure(nz)/&
	(profiles(nx,ny)%press(nz) - (1.- 0.622) * profiles(nx,ny)%vapor_pressure(nz))        ! [kg/kg]
      end do 
    end do 
  end do 
								    
  return

end subroutine get_atmosG0
	      