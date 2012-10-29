subroutine calc_radar_spectrum(nbins1,diameter_spec, back_spec,&
  temp,press,hgt,frequency,particle_type,particle_spec)
!in
!nbins: No of bins + 1
!diameter_spec: Diameter Spectrum (SI)
!back_spec: backscattering cross section per volume in m²/m⁴ (includes number density)
!temp: temperature in K
!frequency in GHz
!particle_type: cl cloud|
!out
!particle_spec particle spectrum in dependence of radar Dopple velocity in m6m-3/ms-1[yet to check!]

  use kinds
  use nml_params
  use constants
  implicit none
  
  integer, parameter :: maxTurbTerms=500

  integer,intent(in) ::  nbins1 !nbins+1
  real(kind=dbl), dimension(nbins1),intent(in):: diameter_spec, back_spec
  character(5),intent(in) ::  particle_type
  real(kind=dbl), intent(in):: temp, frequency, press,hgt

  real(kind=dbl), dimension(nbins1):: vel_spec,dD_dU,back_vel_spec, back_spec_ref
  real(kind=dbl), intent(out), dimension(radar_nfft):: particle_spec
  real(kind=dbl), dimension(radar_nfft):: radar_velo
  real(kind=dbl):: del_v, K2, dielec_water, wavelength,press_correction
!   real(kind=dbl):: del_d, Ze
  integer :: ii, jj

  if (verbose .gt. 1) print*, 'Entering calc_radar_spectrum.f90'

  ! get |K|**2 and lambda
  K2 = dielec_water(0.D0,temp-t_abs,frequency)
  wavelength = c / (frequency*1.d9)   ! m



  !get the particle velocities depending on particle type!
  if (particle_type .eq. "cloud") then
    vel_spec = 30.d0*(diameter_spec*1.d3)**2 ! [m/s], diameter_spec in m!
  else if (particle_type .eq. "rain") then 
    press_correction = 0.d0 !todo
    vel_spec = (9.65 - 10.3*exp(-0.6*diameter_spec)) ! [m/s], diameter_spec in m!
  else if (particle_type .eq. "ice") then 
    press_correction = 0.d0 !todo
    vel_spec = 0.5d0 ! [m/s], diameter_spec in m!
  else if (particle_type .eq. "snow") then 
    press_correction = 0.d0 !todo
    vel_spec = 1.d0 ! [m/s], diameter_spec in m!
  else if (particle_type .eq. "graup") then 
    press_correction = 0.d0 !todo
    vel_spec = 2.d0 ! [m/s], diameter_spec in m!
  else if (particle_type .eq. "hail") then 
    press_correction = 0.d0 !todo
    vel_spec = 5.d0 ! [m/s], diameter_spec in m!
  else
    print*, particle_type, ": did not understand particle_type in calc_radar_spectrum"
    stop
  end if


  back_spec_ref = (1d0/ (K2*pi**5) ) * back_spec * (wavelength)**4 ![m⁶/m⁴]
  back_spec_ref =  back_spec_ref * 1d18 !now non-SI: [mm⁶/m⁴]
  back_spec_ref =  back_spec_ref * 1d-3 !now non-SI: [mm⁶/m³/mm]

!  del_d = (diameter_spec(2)-diameter_spec(1)) * 1.d3 ![mm]
! print*,"SUM(back_spec_ref)*del_d",SUM(back_spec_ref)*del_d*0.5
! print*,"log10 SUM(back_spec_ref)*del_d",10*log10(SUM(back_spec_ref)*del_d*0.5)

! call avint( back_spec_ref, diameter_spec, SIZE(diameter_spec), &
!  diameter_spec(0), diameter_spec(SIZE(diameter_spec)), Ze )
! print*, "Ze", Ze
! print*, "Ze log", 10*log10(Ze)
! print*, "Ze_PAV", Ze*0.5
! print*, "Ze_PAV log", 10*log10(Ze*0.5)

  !move from dimension to velocity!
  do jj=1,nbins1-1
    dD_dU(jj) = (diameter_spec(jj+1)-diameter_spec(jj))/(vel_spec(jj+1)-vel_spec(jj)) ![m/(m/s)]
  end do
  dD_dU(nbins1) = dD_dU(nbins1-1)
  back_vel_spec = back_spec_ref * 1d3 * ABS(dD_dU)  !non-SI: [mm⁶/m³/mm * mm/(m/s)]


  !get delta velocity
  del_v = (radar_max_V-radar_min_V)/radar_nfft ![m/s]
  !create array from min_v to max_v iwth del_v spacing -> velocity spectrum of radar
  radar_velo = (/(((ii*del_v)+radar_min_V),ii=0,radar_nfft)/) ! [m/s]
  
  !interpolate from N(D) bins to radar bins. Why not average?

  call interpolate_spectra(nbins1,radar_nfft,vel_spec,back_vel_spec,radar_velo,particle_spec) ! particle_spec in [mm⁶/m³/mm * mm/(m/s)]

  if (verbose .gt. 1) print*, 'Done calc_radar_spectrum.f90'

  return

end subroutine calc_radar_spectrum