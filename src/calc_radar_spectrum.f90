subroutine calc_radar_spectrum(nbins,diameter_spec, back_spec,&
  temp,press,hgt,frequency,particle_type,particle_spec)
! this routine takes the backscattering spectrum depending on size and converts it
! into a spectrum depending on radar Doppler (=fall) velocity
! based on Spectra_simulator by P. Kollias
! converted from Matlab to Fortran by M. Maahn, IGMK (11/2012)
!in
!nbins: No of bins + 1
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

  integer,intent(in) ::  nbins !for Mie it is actually nbins+1
  real(kind=dbl), dimension(nbins),intent(in):: diameter_spec, back_spec
  character(5),intent(in) ::  particle_type
  real(kind=dbl), intent(in):: temp, frequency, press,hgt

  real(kind=dbl), dimension(nbins):: vel_spec,dD_dU,back_vel_spec, back_spec_ref,&
				del_v_model, diameter_spec_cp
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

  diameter_spec_cp(:) = diameter_spec(:)

  !liu_db has certain maximum and minimum diameter, if this range is overspend, teh diameter is constant.
  !hence, dU_dV becomes infinitive. Thus, we find all bins iwth equal diameter and collect all back in one of them
  if (((particle_type .eq. "snow") .and. (EM_snow .eq. 'liudb')) &
      .or.&
      ((particle_type .eq. "ice") .and. (EM_ice .eq. 'liudb'))) then
    !do the loop to find d_max
    do ii=nbins,2,-1
       !check whetehr equal
       if (diameter_spec_cp(ii) .eq. diameter_spec_cp(ii-1)) then
        !add to next bin
	back_spec_ref(ii-1) = back_spec_ref(ii-1) + back_spec_ref(ii)
        !set original bin to zero
	back_spec_ref(ii) = 0.d0
	!some decreasing random mnumber for diameter...
	diameter_spec_cp(ii) = 0.01*ii
      else
       !only borders are affected, so exit here
       EXIT
      end if
    end do 
!     now d_min
    do ii=1,nbins-1
       if (diameter_spec_cp(ii) .eq. diameter_spec_cp(ii+1)) then
	  back_spec_ref(ii+1) = back_spec_ref(ii+1) + back_spec_ref(ii)
	  back_spec_ref(ii) = 0.d0
	  diameter_spec_cp(ii) = -0.01*(nbins-ii) !some increasing random mnumber...
      else
       EXIT
      end if
    end do 
  end if

print*, "EM_rain ",EM_rain

  !get the particle velocities depending on particle type!
  if (particle_type .eq. "cloud") then
!     vel_spec = 0.75d0+(9.65 - 10.3*exp(-0.6*diameter_spec_cp)) ! [m/s], diameter_spec_cp in m!
    vel_spec = 30.d0*(diameter_spec_cp*1.d3)**2 ! [m/s], diameter_spec_cp in m!
  else if (particle_type .eq. "rain") then 
    press_correction = 0.d0 !todo
    vel_spec = (9.65 - 10.3*exp(-0.6*diameter_spec_cp)) ! [m/s], diameter_spec_cp in m!
  else if (particle_type .eq. "ice") then 
    press_correction = 0.d0 !todo
    vel_spec = 30.d0*(diameter_spec_cp*1.d3)**2 ! [m/s], diameter_spec_cp in m!
  else if (particle_type .eq. "snow") then 
    press_correction = 0.d0 !todo
    vel_spec = 1.d-2*(diameter_spec_cp*1.d3)**2 ! [m/s], diameter_spec_cp in m!
  else if (particle_type .eq. "graup") then 
    press_correction = 0.d0 !todo
    vel_spec = (9.65 - 10.3*exp(-0.6*diameter_spec_cp)) ! [m/s], diameter_spec_cp in m!
  else if (particle_type .eq. "hail") then 
    press_correction = 0.d0 !todo
    vel_spec = (9.65 - 10.3*exp(-0.6*diameter_spec_cp)) ! [m/s], diameter_spec_cp in m!
  else
    print*, particle_type, ": did not understand particle_type in calc_radar_spectrum"
    stop
  end if


  back_spec_ref = (1d0/ (K2*pi**5) ) * back_spec * (wavelength)**4 ![m⁶/m⁴]
  back_spec_ref =  back_spec_ref * 1d18 !now non-SI: [mm⁶/m³/m]


  del_d = (diameter_spec_cp(2)-diameter_spec_cp(1)) * 1.d3 ![mm]
! print*,"SUM(back_spec_ref * 1d-3)*del_d",SUM(back_spec_ref * 1d-3)*del_d*0.5

print*,particle_type," Ze log SUM(back_spec_ref)*del_d",10*log10(SUM(back_spec_ref * 1d-3)*del_d),&
      "assumes equidistant d"



if (nbins > 2) then
 call avint( back_spec_ref * 1d-3, diameter_spec_cp * 1.d3, SIZE(diameter_spec_cp), &
  diameter_spec_cp(1) * 1.d3, diameter_spec_cp(SIZE(diameter_spec_cp)) * 1.d3, Ze )
 !print*, "Ze", Ze
 print*, "Ze log INT(back_spec_ref)", 10*log10(Ze)
else
 print*, "Ze log INT(back_spec_ref)", "too few bins to integrate"
end if
! print*, "Ze_PAV", Ze*0.5
! print*, "Ze_PAV log", 10*log10(Ze*0.5)

  !move from dimension to velocity!
  do jj=1,nbins-1
    dD_dU(jj) = (diameter_spec_cp(jj+1)-diameter_spec_cp(jj))/(vel_spec(jj+1)-vel_spec(jj)) ![m/(m/s)]
    !is all particles fall with the same velocity, dD_dU gets infinitive!
    if (abs(dD_dU(jj)) .ge. huge(dD_dU(jj))) then
      print*, "Stop in calc_radar_spectrum: dD_dU is infinitive", jj,&
	(diameter_spec_cp(jj+1)-diameter_spec_cp(jj)), (vel_spec(jj+1)-vel_spec(jj))
      stop
    end if
    if (verbose .gt. 3) print*,"jj,dD_dU(jj)",jj,dD_dU(jj)
    del_v_model(jj) = ABS(vel_spec(jj+1)-vel_spec(jj))



  end do
  dD_dU(nbins) = dD_dU(nbins-1)
  del_v_model(nbins) = del_v_model(nbins-1)
  back_vel_spec = back_spec_ref * ABS(dD_dU)  !non-SI: [mm⁶/m³/m * m/(m/s)]



  !get delta velocity
  del_v_radar = (radar_max_V-radar_min_V)/radar_nfft ![m/s]
  !create array from min_v to max_v iwth del_v_radar spacing -> velocity spectrum of radar
  radar_velo = (/(((ii*del_v_radar)+radar_min_V),ii=0,radar_nfft)/) ! [m/s]

  if ((MAXVAL(vel_spec) .gt. radar_max_V) .or.(MINVAL(vel_spec) .lt. radar_min_V)) then
    print*, "WARNING: Nyquist range is overspend for particle: ", particle_type
!     stop
  end if
  
  !interpolate OR average (depending whos bins size is greater) from N(D) bins to radar bins. 

!   if (del_v_radar .lt. del_v_model) then
!   if (verbose .gt. 2) print*, "INTEGRATING particle spectrum"
!     call interpolate_spectra(nbins,radar_nfft,vel_spec,back_vel_spec,radar_velo,particle_spec) ! particle_spec in [mm⁶/m³/m * m/(m/s)]
!   else
  if (verbose .gt. 2) print*, "AVERAGING particle spectrum"
    call rescale_spectra(nbins,radar_nfft,vel_spec,back_vel_spec,radar_velo,particle_spec) ! particle_spec in [mm⁶/m³/m * m/(m/s)]
!   end if



print*,particle_type," Ze SUM(back_vel_spec)*del_v_model",10*log10(SUM(back_vel_spec*del_v_model))
print*,particle_type," Ze SUM(particle_spec)*del_v_radar",10*log10(SUM(particle_spec)*del_v_radar)
  if (verbose .gt. 1) print*, 'Done calc_radar_spectrum.f90'

  return

end subroutine calc_radar_spectrum