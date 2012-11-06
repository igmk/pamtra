subroutine calc_radar_spectrum(nbins1,diameter_spec, back_spec,&
  temp,press,hgt,frequency,particle_type,particle_spec)
! this routine takes the backscattering spectrum depending on size and converts it
! into a spectrum depending on radar Doppler (=fall) velocity
! based on Spectra_simulator by P. Kollias
! converted from Matlab to Fortran by M. Maahn, IGMK (11/2012)
!in
!nbins1: No of bins + 1
!diameter_spec: Diameter Spectrum (SI)
!back_spec: backscattering cross section per volume in m²/m⁴ (includes number density)
!temp: temperature in K
!press: air pressure in Pa
!frequency in GHz
!particle_type: cloud|rain|snow|ice|hail|graupel
!out
!particle_spec particle spectrum in dependence of radar Dopple velocity in m6m-3/ms-1

  use kinds
  use nml_params
  use constants
  implicit none
  
  integer, parameter :: maxTurbTerms=500

  integer,intent(in) ::  nbins1 !nbins+1
  real(kind=dbl), dimension(nbins1),intent(in):: diameter_spec, back_spec
  character(5),intent(in) ::  particle_type
  real(kind=dbl), intent(in):: temp, frequency, press,hgt

  real(kind=dbl), dimension(nbins1):: vel_spec,dD_dU,back_vel_spec, back_spec_ref,del_v_model
  real(kind=dbl), intent(out), dimension(radar_nfft):: particle_spec
  real(kind=dbl), dimension(radar_nfft):: radar_velo
  real(kind=dbl):: del_v_radar, K2, dielec_water, wavelength, &
	press_correction
   real(kind=dbl):: del_d, Ze
  integer :: ii, jj

  if (verbose .gt. 1) print*, 'Entering calc_radar_spectrum.f90, particle type: ', particle_type

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
    vel_spec = 30.d0*(diameter_spec*1.d3)**2+6.d0 ! [m/s], diameter_spec in m!
  else if (particle_type .eq. "hail") then 
    press_correction = 0.d0 !todo
    vel_spec = 5.d0 ! [m/s], diameter_spec in m!
  else
    print*, particle_type, ": did not understand particle_type in calc_radar_spectrum"
    stop
  end if


  back_spec_ref = (1d0/ (K2*pi**5) ) * back_spec * (wavelength)**4 ![m⁶/m⁴]
  back_spec_ref =  back_spec_ref * 1d18 !now non-SI: [mm⁶/m³/m]


  del_d = (diameter_spec(2)-diameter_spec(1)) * 1.d3 ![mm]
! print*,"SUM(back_spec_ref * 1d-3)*del_d",SUM(back_spec_ref * 1d-3)*del_d*0.5
print*,particle_type," Ze SUM(back_spec_ref)*del_d",10*log10(SUM(back_spec_ref * 1d-3)*del_d*0.5)

! call avint( back_spec_ref * 1d-3, diameter_spec, SIZE(diameter_spec), &
!  diameter_spec(0), diameter_spec(SIZE(diameter_spec)), Ze )
! print*, "Ze", Ze
! print*, "Ze log", 10*log10(Ze)
! print*, "Ze_PAV", Ze*0.5
! print*, "Ze_PAV log", 10*log10(Ze*0.5)

  !move from dimension to velocity!
  do jj=1,nbins1-1
    dD_dU(jj) = (diameter_spec(jj+1)-diameter_spec(jj))/(vel_spec(jj+1)-vel_spec(jj)) ![m/(m/s)]
    !is all particles fall with the same velocity, dD_dU gets infinitive!
    if (abs(dD_dU(jj)) .ge. huge(dD_dU(jj))) stop "Stop in calc_radar_spectrum: dD_dU is infinitive"
    if (verbose .gt. 3) print*,"jj,dD_dU(jj)",jj,dD_dU(jj)
    del_v_model(jj) = ABS(vel_spec(jj+1)-vel_spec(jj))
print*, jj,(diameter_spec(jj+1)-diameter_spec(jj)),(vel_spec(jj+1)-vel_spec(jj)),dD_dU(jj)



  end do
  dD_dU(nbins1) = dD_dU(nbins1-1)
  del_v_model(nbins1) = del_v_model(nbins1-1)
  back_vel_spec = back_spec_ref * ABS(dD_dU)  !non-SI: [mm⁶/m³/m * m/(m/s)]



  !get delta velocity
  del_v_radar = (radar_max_V-radar_min_V)/radar_nfft ![m/s]
  !create array from min_v to max_v iwth del_v_radar spacing -> velocity spectrum of radar
  radar_velo = (/(((ii*del_v_radar)+radar_min_V),ii=0,radar_nfft)/) ! [m/s]


print*, "TODO make sure that nyquist range is not overspend!"
  
  !interpolate OR average (depending whos bins size is greater) from N(D) bins to radar bins. 

!   if (del_v_radar .lt. del_v_model) then
!   if (verbose .gt. 2) print*, "INTEGRATING particle spectrum"
!     call interpolate_spectra(nbins1,radar_nfft,vel_spec,back_vel_spec,radar_velo,particle_spec) ! particle_spec in [mm⁶/m³/m * m/(m/s)]
!   else
  if (verbose .gt. 2) print*, "AVERAGING particle spectrum"
    call average_spectra(nbins1,radar_nfft,vel_spec,back_vel_spec,radar_velo,particle_spec) ! particle_spec in [mm⁶/m³/m * m/(m/s)]
!   end if

print*,particle_type," Ze SUM(back_spec_ref)*del_d",      10*log10(SUM(back_spec_ref * 1d-3)*del_d*0.5)
print*,particle_type," Ze SUM(back_vel_spec)*del_v_model",10*log10(SUM(back_vel_spec*del_v_model)*0.5)
print*,particle_type," Ze SUM(particle_spec)*del_v_radar",10*log10(SUM(particle_spec)*del_v_radar*0.5)
  if (verbose .gt. 1) print*, 'Done calc_radar_spectrum.f90'

  return

end subroutine calc_radar_spectrum