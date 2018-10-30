module ocean_sfc_optics

  use kinds
  use report_module
  use fastemx, only: compute_fastemx

  use settings, only: nummu, mu_values
  use constants, only: rad2deg
  use vars_rt, only: rt_sfc_emissivity, rt_sfc_reflectivity
  use vars_index, only: i_x, i_y
  use vars_atmosphere, only: atmo_groundtemp, atmo_wind10u, atmo_wind10v, sfc_salinity
  
  implicit none
  
  contains

  subroutine ocean_sfc_optics_tessem2(errorstatus,freq)

  integer(long) :: i
  real(dbl), intent(in) :: freq
  
  real(dbl) :: zenith_angle,&
    wind_speed
  
  ! Error handling

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=24) :: nameOfRoutine = 'ocean_sfc_optics_tessem2'
    
  errorstatus = 0
  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
  
  errorstatus = err

  wind_speed = sqrt(atmo_wind10u(i_x,i_y)**2+atmo_wind10v(i_x,i_y)**2)

  do i = 1, nummu
    zenith_angle = acos(mu_values(i))*rad2deg
  
  call tessem2(err,&
    freq,&
    zenith_angle,&
    wind_speed,&
    atmo_groundtemp(i_x,i_y),&
    sfc_salinity(i_x,i_y),&
    rt_sfc_emissivity(:,i),&
    rt_sfc_reflectivity(:,i))
        if (i == nummu) then
           rt_sfc_emissivity(:,i) = (rt_sfc_emissivity(1,i) + rt_sfc_emissivity(2,i)) / 2.
           rt_sfc_reflectivity(:,i) = (rt_sfc_reflectivity(1,i) + rt_sfc_reflectivity(2,i)) / 2.
       end if
  end do
  
  if (err /= 0) then
          msg = 'error in tessem2 calculation'
          call report(err,msg, nameOfRoutine)
          errorstatus = err
          return
  end if
  
  errorstatus = err
      
  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine ocean_sfc_optics_tessem2
  
  subroutine ocean_sfc_optics_fastemx(errorstatus,freq)

  integer(long) :: i
  real(dbl), intent(in) :: freq
  
  real(dbl) :: zenith_angle,&
    wind_speed
    

  ! Error handling

  integer(kind=long), intent(out) :: errorstatus
  character(len=24) :: nameOfRoutine = 'ocean_sfc_optics_fastemx'
    
  errorstatus = 0
  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

  wind_speed = sqrt(atmo_wind10u(i_x,i_y)**2+atmo_wind10v(i_x,i_y)**2)
  
  do i = 1, nummu
    zenith_angle = acos(mu_values(i))*rad2deg
    call compute_fastemx( &
      freq, &     ! Input
      zenith_angle, &  ! Input
      atmo_groundtemp(i_x,i_y), &   ! Input
      sfc_salinity(i_x,i_y), &      ! Input
      wind_speed, &    ! Input
      rt_sfc_emissivity(:,i), &    ! Output
      rt_sfc_reflectivity(:,i))!, & ! Output
!      Azimuth_Angle, & ! Optional input
!      Transmittance  ) ! Optional input
  end do
  
  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine ocean_sfc_optics_fastemx
  
end module ocean_sfc_optics
