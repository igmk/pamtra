subroutine get_atmosg(press, temp, vapor_pressure, rho_vap, nlyr, freq, gas_mod, abscoef)!, abscoef_o2,abscoef_h2o,abscoef_n2)
!                                                                        
!     Calculate average air pressure and vapor pressure in specified    
!     layers, given the temperature, pressure, and relative humidity    
!     from cloud-resolving model output.                                
!     vapor_pressure viene data in mb = hPa

  use kinds
  
  implicit none
      
  integer, intent(in) :: nlyr
  integer :: nz
  real(kind=dbl), dimension(nlyr), intent(in) :: press
  real(kind=dbl), dimension(nlyr), intent(in) :: temp
  real(kind=dbl), dimension(nlyr), intent(in) :: vapor_pressure
  real(kind=dbl), dimension(nlyr), intent(in) :: rho_vap
  real(kind=dbl), dimension(nlyr), intent(out) ::  abscoef
 ! real(kind=dbl), dimension(nlyr), intent(out) ::  abscoef_o2, abscoef_h2o, abscoef_n2
  real(kind=dbl) :: absair, abswv
!  real(kind=dbl) :: abs_o2, abs_n2                                                                        
  real(kind=dbl), intent(in) :: freq
  real(kind=dbl) :: tc
  
  character(3) :: gas_mod
      
    do nz = 1, nlyr                     
      tc = temp(nz) - 273.15 
      if (gas_mod .eq. 'L93') then
	    call mpm93(freq, press(nz)/1.d3, vapor_pressure/1.d3,tc, 0.0d0, abscoef(nz))
	    abscoef(nz) = abscoef(nz)/1000.
      else if (gas_mod .eq. 'R98') then
      ! Rosenkranz 1998 gas absorption model
      ! Input parameters:
      !    frequency                    GHz
      !    temperature                   K
      !    water vapor density         kg/m**3
      !    total air pressure            Pa
      !
      ! Output parameters:
      !    extinction by dry air       Np/km
      !    extinction by water vapor   Np/km      
	    call gasabsr98(freq,temp(nz),rho_vap(nz),press(nz),absair,abswv)!,abs_n2,abs_o2)

        abscoef(nz) = (absair + abswv)/1000.    ! [1/m]

      else
	    abscoef(nz) = 0
        print*, "No gas absorption model specified!"
	    stop
      end if
    end do
                                                                        
 return 

end subroutine get_atmosg                     
