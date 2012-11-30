subroutine radar_simulator(particle_spectrum,back,kexthydro,&
      frequency,temp,delta_h,nz,nx,ny,fi)
! This routine takes the backscattering spectrum depending on Doppler velocity, 
! adds noise and turbulence and simulates temporal averaging
!
! based on Spectra_simulator by P. Kollias
! converted from Matlab to Fortran by M. Maahn (2012)

! in
! particle_spectrum: backscattering particle spectrum per Doppler velocity [mm⁶/m³/(m/s)] NON-SI
! back: volumetric backscattering crossection in m²/m³
! kexthydro: hydrometeor absorption coefficient [Np/m]
! frequency: Frequency in GHz
! temp: temperature in K
! delta_h: heigth of layer in m
! nz,nx,ny,fi: level, grid x, y, frequency index
! 
! out is saved directly to vars_output module


  use kinds
  use nml_params
  use constants
  use vars_output, only: radar_spectra, radar_snr, radar_vel,&
	  radar_moments, radar_slope, radar_quality, Ze, Att_hydro !output of the radar simulator
  implicit none
  

  real(kind=dbl),intent(in) ::  frequency, temp, back, delta_h,kexthydro
  integer,intent(in) ::  nz,nx,ny,fi
  real(kind=dbl), dimension(radar_nfft_aliased),intent(in):: particle_spectrum
  real(kind=dbl), dimension(radar_nfft_aliased) :: spectra_velo_aliased
  real(kind=dbl), dimension(radar_maxTurbTerms):: turb
  real(kind=dbl), dimension(radar_nfft*radar_no_Ave):: x_noise
  real(kind=dbl), dimension(radar_no_Ave,radar_nfft):: noise_turb_spectra_tmp
  real(kind=dbl), dimension(radar_nfft):: noise_turb_spectra,&
      snr_turb_spectra,spectra_velo, turb_spectra_aliased, noise_removed_turb_spectra
  real(kind=dbl), dimension(:),allocatable:: turb_spectra
  real(kind=dbl), dimension(0:4):: moments  
  real(kind=dbl), dimension(2):: slope  
  real(kind=dbl):: SNR, del_v, ss, K2, wavelength, Ze_back, dielec_water, K, &
      min_V_aliased, max_V_aliased
  integer :: ii, tt, turbLen,alloc_status,ts_imin, ts_imax, startI, stopI

  if (verbose .gt. 1) print*, 'Entering radar_simulator.f90'

  ! get |K|**2 and lambda
  K2 = dielec_water(0.D0,temp-t_abs,frequency)
  wavelength = c / (frequency*1.d9)   ! [m]

  !first, calculate the attenuation for hydrometeors
  Att_hydro(nx,ny,nz,fi) = 10*log10(exp(kexthydro*delta_h))

  !transform backscattering in linear reflectivity units, 10*log10(back) would be in dBz
  Ze_back = 1.d18* (1.d0/ (K2*pi**5) ) * back * (wavelength)**4 ![mm⁶/m³] 

  if (radar_mode == "simple") then

    Ze(nx,ny,nz,fi) = 10*log10(Ze_back)

  else if ((radar_mode == "moments") .or. (radar_mode == "spectrum")) then


    !get delta velocity
    del_v = (radar_max_V-radar_min_V) / radar_nfft ![m/s]
    !create array from min_v to max_v iwth del_v spacing -> velocity spectrum of radar
    spectra_velo = (/(((ii*del_v)+radar_min_V),ii=0,radar_nfft)/) ! [m/s]

    !same for the extended spectrum
    min_V_aliased = radar_min_V - radar_aliasing_nyquist_interv*(radar_max_V-radar_min_V)
    max_V_aliased = radar_max_V + radar_aliasing_nyquist_interv*(radar_max_V-radar_min_V)
    spectra_velo_aliased = (/(((ii*del_v)+min_V_aliased),ii=0,radar_nfft_aliased)/) ! [m/s]

    !get turbulence
    if (radar_turbulence_st > 0.d0) then
      ss = radar_turbulence_st/del_v; !in array indices!

      turb(:) = 0.d0
      tt = 1
      do while (tt .le. 24.d0/del_v+1) 
	if (tt .gt. radar_maxTurbTerms) then
	  print*,radar_maxTurbTerms, INT(12.d0/del_v+1.d0),&
	  ": maximum of turbulence terms reached. increase radar_maxTurbTerms (nml_params.f90)"
	  stop
	end if
	turb(tt) = 1.d0/(sqrt(2.d0*pi)*ss) * exp(-(tt-(12.d0/del_v+1))**2.d0/(2.d0*ss**2.d0))
	tt = tt+1
      end do

      turbLen=tt-1

      if (SIZE(particle_spectrum)+turbLen-1 .lt. floor(12.d0/del_v+1)+radar_nfft-1) then
	print*, SIZE(particle_spectrum)+turbLen-1,floor(12.d0/del_v+1)+radar_nfft-1,&
	    "vector resulting from convolution to short!  (radar_simulator.f90)"
	stop
      end if

      allocate(turb_spectra(SIZE(particle_spectrum)+turbLen-1),stat=alloc_status)

      !convolute spectrum and noise
      call convolution(particle_spectrum,SIZE(particle_spectrum),turb(1:turbLen),turbLen,turb_spectra)

      !I don't like Nans here'
      where(ISNAN(turb_spectra)) turb_spectra = 0.d0
      ts_imin = floor(12.d0/del_v+1)
      ts_imax = floor(12.d0/del_v+1)+radar_nfft_aliased-1
    else
      !add no turbulence
      ts_imin = 1
      ts_imax = radar_nfft_aliased
      allocate(turb_spectra(radar_nfft_aliased),stat=alloc_status)
      turb_spectra = particle_spectrum
    end if


    if ((turb_spectra(ts_imin) .ne. 0.d0) .or.(turb_spectra(ts_imax) .ne. 0.d0)) then
      print*, "WARNING: radar_aliasing_nyquist_interv too small to handle aliasing effects, increase it!"
      stop
    end if
    
    !lets look for aliasing effects. if we calculated radar_aliasing_nyquist_interv for a broader spectrum than necessary, fold it again:
    if (radar_aliasing_nyquist_interv > 0) then
      turb_spectra_aliased(:) = 0.d0
      do, ii=1,1+2*radar_aliasing_nyquist_interv
	!get indices
	startI = ts_imin + (ii-1)*radar_nfft
	stopI =  ts_imax -  (1+2*radar_aliasing_nyquist_interv - ii)*radar_nfft
	!appy aliasing
	turb_spectra_aliased = turb_spectra_aliased + turb_spectra(startI:stopI)
      end do
    else
      turb_spectra_aliased = turb_spectra(ts_imin:ts_imax)
    end if

    !get the SNR
    SNR = 10.d0*log10(Ze_back/radar_Pnoise)

    !this here is for scaling, if we have now a wrong Ze due to all the turbulence, rescaling etc...
    K = (Ze_back/SUM(turb_spectra_aliased*del_v))
    snr_turb_spectra = (K* turb_spectra_aliased + radar_Pnoise/(radar_nfft*del_v))
  !   snr_turb_spectra =turb_spectra_aliased + radar_Pnoise/(radar_nfft*del_v)

    !init_random_seed works with system clock, so if called very often it creates the same random numbers. Thus, create a big one now
    call init_random_seed()
    call RANDOM_NUMBER(x_noise)

    do tt = 1, radar_no_Ave
      noise_turb_spectra_tmp(tt,:) = -log(x_noise((tt-1)*radar_nfft+1:tt*radar_nfft))*snr_turb_spectra
    end do

    
    if (radar_no_Ave .eq. 1) then
      noise_turb_spectra = noise_turb_spectra_tmp(1,:)
    else
      noise_turb_spectra = SUM(noise_turb_spectra_tmp,DIM=1)/radar_no_Ave
    end if

  !apply spectral resolution
  noise_turb_spectra = noise_turb_spectra * del_v

  if (verbose .gt. 3) then
    print*,"second K",K
    print*,"TOTAL"," Ze back",10*log10(Ze_back)
    print*,"TOTAL"," Ze SUM(particle_spectrum)*del_v",10*log10(SUM(particle_spectrum)*del_v)
    print*,"TOTAL"," Ze SUM(turb_spectra)*del_v",10*log10(SUM(turb_spectra)*del_v)
    print*,"TOTAL"," Ze SUM(turb_spectra_aliased)*del_v",10*log10(SUM(turb_spectra_aliased)*del_v)
    print*,"TOTAL"," Ze SUM(snr_turb_spectra)*del_v",10*log10(SUM(snr_turb_spectra)*del_v)
    print*,"TOTAL"," Ze SUM(noise_turb_spectra)*del_v",10*log10(SUM(noise_turb_spectra))
    print*,"TOTAL"," Ze SUM(snr_turb_spectra)*del_v-radar_Pnoise",10*log10(SUM(snr_turb_spectra)*del_v-radar_Pnoise)
    print*,"TOTAL"," Ze SUM(noise_turb_spectra)*del_v-radar_Pnoise",10*log10(SUM(noise_turb_spectra)-radar_Pnoise)
  end if


  call radar_calc_moments(noise_turb_spectra,noise_removed_turb_spectra,moments,slope)
  if (verbose .gt. 3) then
    print*,"TOTAL"," Ze moments",10*log10(moments(0))
    print*,"#####################"
  end if

      ! collect results for output
  !   radar_spectra(nx,ny,nz,fi,:) = 10*log10(particle_spectrum(513:1024))

    !if wanted, apply the noise correction to the spectrum to be saved.
    if (radar_save_noise_corrected_spectra) noise_turb_spectra = noise_removed_turb_spectra

    radar_spectra(nx,ny,nz,fi,:) = 10*log10(noise_turb_spectra)
    radar_snr(nx,ny,nz,fi) = SNR
    radar_vel(:) = spectra_velo(:)
    radar_moments(nx,ny,nz,fi,:) = moments(1:4)
    radar_slope(nx,ny,nz,fi,:) = slope(:)
    radar_quality(nx,ny,nz,fi) = 0
    Ze(nx,ny,nz,fi) = 10*log10(moments(0))


    deallocate(turb_spectra)
  else
    print*, "did not understand radar_mode", radar_mode
    stop

  end if



  if (verbose .gt. 1) print*, 'Exiting radar_simulator.f90'

  return
end subroutine radar_simulator