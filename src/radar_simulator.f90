subroutine radar_simulator(errorstatus,particle_spectrum,back,kexthydro,&
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
    use settings
    use constants
    use vars_output, only: radar_spectra, radar_snr, radar_vel,&
    radar_moments, radar_slope, radar_quality, Ze, Att_hydro !output of the radar simulator
    use report_module
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
    integer::quality_2ndPeak, quailty_aliasing
    real(kind=dbl), dimension(:),allocatable:: turb_spectra
    real(kind=dbl), dimension(0:4):: moments
    real(kind=dbl), dimension(2):: slope
    real(kind=dbl):: SNR, del_v, ss, K2, wavelength, Ze_back, dielec_water, K, &
    min_V_aliased, max_V_aliased
    integer :: ii, tt, turbLen,alloc_status,ts_imin, ts_imax, startI, stopI

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'radar_simulator'    
    
    interface
        subroutine convolution(errorstatus,X,M,A,N,Y)
            use kinds
            implicit none
	    integer(kind=long), intent(out) :: errorstatus
            INTEGER, intent(in) :: M  ! Size of input vector X
            INTEGER, intent(in) :: N   ! Size of convolution filter A
            REAL(kind=dbl), intent(in), DIMENSION(M) :: X
            REAL(kind=dbl), intent(in), DIMENSION(N) :: A
            REAL(kind=dbl), intent(out), DIMENSION(M+N-1) :: Y
        end subroutine convolution

        subroutine radar_calc_moments(errorstatus,radar_spectrum_in,radar_spectrum_out,moments,slope,quality)
            use kinds
            use settings, only: radar_nfft
            implicit none
	    integer(kind=long), intent(out) :: errorstatus
            real(kind=dbl), dimension(radar_nfft), intent(in):: radar_spectrum_in
            real(kind=dbl), dimension(radar_nfft), intent(out):: radar_spectrum_out
            real(kind=dbl), dimension(0:4), intent(out):: moments
            real(kind=dbl), dimension(2), intent(out):: slope
            integer, intent(out) :: quality
        end subroutine radar_calc_moments

        subroutine random(errorstatus,n, pseudo, x_noise)
            use kinds
            implicit none
	    integer(kind=long), intent(out) :: errorstatus
            integer, intent(in) :: n
            logical, intent(in) :: pseudo
            real(kind=dbl), intent(out), dimension(n) :: x_noise
        end subroutine random
  
    end interface

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

    if (ANY(ISNAN(particle_spectrum))) then
	print*,particle_spectrum
	errorstatus = fatal
	msg = "got nan in values in backscattering spectrum"
	call report(errorstatus, msg, nameOfRoutine)
	return
    end if

    if (ISNAN(back) .or. back < 0.d0) then
	print*,back
	errorstatus = fatal
	msg = "got nan or negative vaue in linear Ze"
	call report(errorstatus, msg, nameOfRoutine)
	return
    end if
      
    
    ! get |K|**2 and lambda
    K2 = dielec_water(0.D0,temp-t_abs,frequency)
    wavelength = c / (frequency*1.d9)   ! [m]

    !first, calculate the attenuation for hydrometeors
    Att_hydro(nx,ny,nz,fi) = 10*log10(exp(kexthydro*delta_h))

    !transform backscattering in linear reflectivity units, 10*log10(back) would be in dBz
    Ze_back = 1.d18* (1.d0/ (K2*pi**5) ) * back * (wavelength)**4 ![mm⁶/m³]

    if (radar_mode == "simple") then
	if (Ze_back .eq. 0.d0) then
	  Ze(nx,ny,nz,fi) = -9999.d0
        else 
	  Ze(nx,ny,nz,fi) = 10*log10(Ze_back)
	end if
      if (verbose >= -3) print*, "nx,ny,nz,fi,ze", nx,ny,nz,fi,Ze(nx,ny,nz,fi)
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
            ss = radar_turbulence_st/del_v;            !in array indices!

            turb(:) = 0.d0
            tt = 1
            do while (tt .le. 24.d0/del_v+1)
                if (tt .gt. radar_maxTurbTerms) then
                    print*,radar_maxTurbTerms, INT(12.d0/del_v+1.d0)
		    errorstatus = fatal
		    msg = "maximum of turbulence terms reached. increase radar_maxTurbTerms (settings.f90)"
		    call report(errorstatus, msg, nameOfRoutine)
		    return
                end if
                turb(tt) = 1.d0/(sqrt(2.d0*pi)*ss) * exp(-(tt-(12.d0/del_v+1))**2.d0/(2.d0*ss**2.d0))
                tt = tt+1
            end do

            turbLen=tt-1

            if (SIZE(particle_spectrum)+turbLen-1 .lt. floor(12.d0/del_v+1)+radar_nfft-1) then
	      print*, SIZE(particle_spectrum)+turbLen-1,floor(12.d0/del_v+1)+radar_nfft-1
	      errorstatus = fatal
	      msg =  "vector resulting from convolution to short!"
	      call report(errorstatus, msg, nameOfRoutine)
	      return
            end if

            allocate(turb_spectra(SIZE(particle_spectrum)+turbLen-1),stat=alloc_status)

            !convolute spectrum and noise
            call convolution(err,particle_spectrum,SIZE(particle_spectrum),turb(1:turbLen),turbLen,turb_spectra)
	    if (err /= 0) then
		msg = 'error in convolution!'
		call report(err, msg, nameOfRoutine)
		errorstatus = err
		return
	    end if   
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


        !     if ((turb_spectra(ts_imin) .ne. 0.d0) .or.(turb_spectra(ts_imax) .ne. 0.d0)) then
        !       print*, "WARNING: radar_aliasing_nyquist_interv too small to handle aliasing effects, increase it!"
        !       stop
        !     end if
    
    
        quailty_aliasing = 0
        !lets look for aliasing effects. if we calculated radar_aliasing_nyquist_interv for a broader spectrum than necessary, fold it again:
        if (radar_aliasing_nyquist_interv > 0) then
            turb_spectra_aliased(:) = 0.d0
            do, ii=1,1+2*radar_aliasing_nyquist_interv
                !get indices
                startI = ts_imin + (ii-1)*radar_nfft
                stopI =  ts_imax -  (1+2*radar_aliasing_nyquist_interv - ii)*radar_nfft
                if ((ii .ne. radar_aliasing_nyquist_interv + 1) &
                .and. (SUM(turb_spectra(startI:stopI)) .gt. 0.d0)) then
                    quailty_aliasing = 1
                end if
                !appy aliasing
                turb_spectra_aliased = turb_spectra_aliased + turb_spectra(startI:stopI)
            end do
        else
            turb_spectra_aliased = turb_spectra(ts_imin:ts_imax)
        end if

        if (verbose > 2) then
            if (quailty_aliasing .ne. 0) then
                print*, "radar quality: aliasing found"
            else
                print*, "radar quality: NO aliasing found"
            end if
        end if


        !get the SNR
        SNR = 10.d0*log10(Ze_back/radar_Pnoise)
        !this here is for scaling, if we have now a wrong Ze due to all the turbulence, rescaling etc...
        K = (Ze_back/SUM(turb_spectra_aliased*del_v))
        snr_turb_spectra = (K* turb_spectra_aliased + radar_Pnoise/(radar_nfft*del_v))
        !   snr_turb_spectra =turb_spectra_aliased + radar_Pnoise/(radar_nfft*del_v)



        if (radar_no_Ave .eq. 0) then !0 means infinity-> no noise
            noise_turb_spectra = snr_turb_spectra
        else

            !get noise. if jacobian_mode, random number generator is always initiated with the same number
            if (verbose > 2) print*, "get noise"
            call random(err,radar_no_Ave*radar_nfft,jacobian_mode,x_noise)
	    if (err /= 0) then
		msg = 'error in random!'
		call report(err, msg, nameOfRoutine)
		errorstatus = err
		return
	    end if   
            do tt = 1, radar_no_Ave
                noise_turb_spectra_tmp(tt,:) = -log(x_noise((tt-1)*radar_nfft+1:tt*radar_nfft))*snr_turb_spectra
            end do

            if (radar_no_Ave .eq. 1) then
                noise_turb_spectra = noise_turb_spectra_tmp(1,:)
            else
                noise_turb_spectra = SUM(noise_turb_spectra_tmp,DIM=1)/radar_no_Ave
            end if
        end if

        !apply spectral resolution
        noise_turb_spectra = noise_turb_spectra * del_v !now [mm⁶/m³]

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


        call radar_calc_moments(err,noise_turb_spectra,noise_removed_turb_spectra,moments,slope,quality_2ndPeak)
	if (err /= 0) then
	  msg = 'error in radar_calc_moments!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
      end if   
        if (verbose .gt. 3) then
            print*,"TOTAL"," Ze moments",10*log10(moments(0))
            print*,"#####################"
        end if

          ! collect results for output
        !   radar_spectra(nx,ny,nz,fi,:) = 10*log10(particle_spectrum(513:1024))

        !if wanted, apply the noise correction to the spectrum to be saved.
        if (radar_save_noise_corrected_spectra) noise_turb_spectra = noise_removed_turb_spectra

        WHERE (ISNAN(noise_turb_spectra)) noise_turb_spectra = -9999.d0
        IF (ISNAN(moments(0))) moments(0) = -9999.d0

        radar_spectra(nx,ny,nz,fi,:) = 10*log10(noise_turb_spectra)
        radar_snr(nx,ny,nz,fi) = SNR
        radar_vel(:) = spectra_velo(:)
        radar_moments(nx,ny,nz,fi,:) = moments(1:4)
        radar_slope(nx,ny,nz,fi,:) = slope(:)
        radar_quality(nx,ny,nz,fi) = quailty_aliasing + quality_2ndPeak
        Ze(nx,ny,nz,fi) = 10*log10(moments(0))


        deallocate(turb_spectra)
    else
      errorstatus = fatal
      msg =   "did not understand radar_mode"// radar_mode
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if

    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return
end subroutine radar_simulator
