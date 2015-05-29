module ocean_sfc_optics

  use kinds
  use report_module
  use crtm_fastemx, only: compute_fastemx
  use settings, only: nummu, mu_values, salinity
  use constants, only: rad2deg
  use vars_rt, only: rt_sfc_emissivity, rt_sfc_reflectivity
  use vars_index, only: i_x, i_y
  use vars_atmosphere, only: atmo_groundtemp, atmo_wind10u, atmo_wind10v
  
  implicit none
  
  contains

  subroutine ocean_sfc_optics_fastem5(errorstatus,freq)

  integer(long) :: i
  integer(long), parameter :: sl = 80
  real(dbl), intent(in) :: freq!, &
  
  REAL(dbl) :: Azimuth_Angle
  REAL(dbl) :: Transmittance 
  
  real(dbl) :: zenith_angle,&
    wind_speed
  character(sl) :: nameOfRoutine = 'ocean_sfc_optics_fastem5'

  ! Error handling

  integer(kind=long), intent(out) :: errorstatus
    
  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

  wind_speed = sqrt(atmo_wind10u(i_x,i_y)**2+atmo_wind10v(i_x,i_y)**2)
  
  do i = 1, nummu
    zenith_angle = acos(mu_values(i))*rad2deg
    call Compute_FastemX( &
      freq, &     ! Input
      zenith_angle, &  ! Input
      atmo_groundtemp(i_x,i_y), &   ! Input
      salinity, &      ! Input
      wind_speed, &    ! Input
      rt_sfc_emissivity(:,i), &    ! Output
      rt_sfc_reflectivity(:,i))!, & ! Output
!      Azimuth_Angle, & ! Optional input
!      Transmittance  ) ! Optional input
  end do
  
  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine ocean_sfc_optics_fastem5
  
end module ocean_sfc_optics