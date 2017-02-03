subroutine radar_simulator(errorstatus,particle_spectrum,back,kexthydro,delta_h)
    ! This routine takes the backscattering spectrum depending on Doppler velocity,
    ! adds noise and turbulence and simulates temporal averaging
    !
    ! based on Spectra_simulator by P. Kollias
    ! converted from Matlab to Fortran by M. Maahn (2012)
    !
    ! out is saved directly to vars_output module


    use kinds
    use settings
    use constants
    use radar_moments, only: radar_calc_moments
    use vars_atmosphere, only: atmo_airturb, &
      atmo_radar_prop, &
      atmo_lat, &
      atmo_lon, &
      atmo_nlyrs, &
      atmo_wind_uv, &
      atmo_turb_edr
    use vars_output, only: out_radar_spectra, out_radar_snr, out_radar_vel,out_radar_hgt, &
    out_radar_moments, out_radar_slopes, out_radar_edges, out_radar_quality, out_ze, out_att_hydro, & !output of the radar simulator
      out_att_atmo, &
      out_debug_radarvel, &
      out_debug_radarback_wturb, &
      out_debug_radarback_wturb_wnoise
    use radar_spectral_broadening, only: estimate_spectralBroadening  
    use report_module
    use vars_index, only: i_x,i_y, i_z, i_f, i_p, i_n

    implicit none

    real(kind=dbl), dimension(radar_npol),intent(in) ::  back !volumetric backscattering crossection in m²/m³
    real(kind=dbl),intent(in) ::  delta_h !heigth of layer in m
    real(kind=dbl),intent(in) ::  kexthydro !hydrometeor absorption coefficient [Np/m]
    real(kind=dbl), dimension(radar_npol,radar_nfft_aliased),intent(in):: particle_spectrum !backscattering particle spectrum per Doppler velocity [mm⁶/m³/(m/s)] NON-SI

    real(kind=dbl), dimension(radar_nfft_aliased) :: particle_spectrum_att
    real(kind=dbl), dimension(radar_nfft_aliased) :: spectra_velo_aliased
    real(kind=dbl), dimension(radar_nfft_aliased):: turb
    real(kind=dbl), dimension(radar_nfft*radar_no_Ave(i_f)):: x_noise
    real(kind=dbl), dimension(radar_no_Ave(i_f),radar_nfft):: noise_turb_spectra_tmp
    real(kind=dbl), dimension(radar_nfft):: noise_turb_spectra,&
    snr_turb_spectra,spectra_velo, turb_spectra_aliased, noise_removed_turb_spectra
    integer::quality_moments, quailty_aliasing
    real(kind=dbl), dimension(2) :: rand_number
    real(kind=dbl), dimension(2*radar_nfft_aliased-1):: turb_spectra
    real(kind=dbl), dimension(0:4,radar_nPeaks):: moments
    real(kind=dbl), dimension(2,radar_nPeaks):: slope
    real(kind=dbl), dimension(2,radar_nPeaks):: edge
    real(kind=dbl) :: noise_out, PIA
    real(kind=dbl) :: noise, noise_max
    real(kind=dbl):: SNR, del_v, ss, K2, wavelength, Ze_back, K, &
    min_V_aliased, max_V_aliased, receiver_uncertainty, radar_Pnoise, frequency
    integer(kind=long) :: ii, tt,ts_imin, ts_imax, startI, stopI, seed
    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=15) :: nameOfRoutine = 'radar_simulator'

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

        subroutine random(errorstatus,n, seedval, x_noise)
            use kinds
            implicit none
	    integer(kind=long), intent(out) :: errorstatus
            integer, intent(in) :: n
            integer, intent(in) :: seedval
            real(kind=dbl), intent(out), dimension(n) :: x_noise
        end subroutine random

        subroutine radar_hildebrand_sekhon(errorstatus,spectrum,n_ave,n_ffts,&
          noise_mean,noise_max)
            use kinds
            implicit none
            integer(kind=long), intent(out) :: errorstatus
            integer, intent(in) :: n_ave, n_ffts
            real(kind=dbl), dimension(n_ffts), intent(in) :: spectrum
            real(kind=dbl), intent(out) :: noise_mean
            real(kind=dbl), intent(out) :: noise_max
        end subroutine radar_hildebrand_sekhon
    end interface
    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    err = 0

!     call assert_true(err,(radar_nPeaks == 1),&
!         "only radar_nPeaks=1 allowed as of today")
    call assert_false(err,(ANY(ISNAN(particle_spectrum))),&
        "got nan in values in backscattering spectrum")
    call assert_false(err,(ANY(ISNAN(back)) .or. ANY(back < 0.d0)),&
        "got nan or negative value in linear Ze")
    call assert_false(err,(SUM(particle_spectrum) < 0.d0),&
        "sum particle_spectrum < 0")
    if (err > 0) then
      errorstatus = fatal
      msg = "assertation error"
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if

    if (verbose >= 10)print*, "particle_spectrum"
    if (verbose >= 10)print*, particle_spectrum

    frequency = freqs(i_f)
    ! get |K|**2 and lambda

    K2 = radar_K2(i_f)!dielec_water(0.D0,radar_K2_temp-t_abs,frequency)
    wavelength = c / (frequency*1.d9)   ! [m]

    !first, calculate the attenuation for hydrometeors
    out_att_hydro(i_x,i_y,i_z,i_f,1) = 10*log10(exp(kexthydro*delta_h))

    do i_p = 1  , radar_npol
      if (verbose >= 4) print*, "do polarisation", i_p, radar_pol(i_p)

      ! if (back(i_p) == 0.0) then
      !   if (verbose >= 3) print*, "Skipping polarisation", i_p, radar_pol(i_p), "backscattering is zero"
      !   CYCLE
      ! end if

      !transform backscattering in linear reflectivity units, 10*log10(back) would be in dBz
      Ze_back = 1.d18* (1.d0/ (K2*pi**5) ) * back(i_p) * (wavelength)**4 ![mm⁶/m³]

      !take care of path integrated attenuation
      PIA = 0.d0
      if (TRIM(radar_attenuation) == "top-down") then
        if (i_z < atmo_nlyrs(i_x,i_y)) then
          PIA = 2 * (SUM(out_att_hydro(i_x,i_y,atmo_nlyrs(i_x,i_y):i_z+1:-1,i_f,1)) + &
              SUM(out_att_atmo(i_x,i_y,atmo_nlyrs(i_x,i_y):i_z+1:-1,i_f)))
        end if
        PIA = PIA + out_att_hydro(i_x,i_y,i_z,i_f,1) + out_att_atmo(i_x,i_y,i_z,i_f)
      else if (TRIM(radar_attenuation) == "bottom-up") then
        if (i_z > 1) then
          PIA = 2 * (SUM(out_att_hydro(i_x,i_y,1:i_z-1,i_f,1)) + SUM(out_att_atmo(i_x,i_y,1:i_z-1,i_f)))
        end if
        PIA = PIA + out_att_hydro(i_x,i_y,i_z,i_f,1) + out_att_atmo(i_x,i_y,i_z,i_f)
      else if (TRIM(radar_attenuation) /= "disabled") then
        errorstatus = fatal
        msg = "do not understand radar_attenuation: "//radar_attenuation
        call report(errorstatus, msg, nameOfRoutine)
        return
      end if

      PIA = 10d0**(0.1d0*PIA) !linearize
      Ze_back = Ze_back/PIA
      particle_spectrum_att = particle_spectrum(i_p,:)/PIA
      if (verbose >= 10)print*, "particle_spectrum_att"
      if (verbose >= 2)print*, particle_spectrum_att

      if (radar_mode == "simple") then
        call assert_true(err,(radar_nPeaks == 1),&
            "radar_nPeaks=1 required for 'simple' radar_mode")
        if (err > 0) then
          errorstatus = fatal
          msg = "assertation error"
          call report(errorstatus, msg, nameOfRoutine)
          return
        end if
          if (Ze_back .eq. 0.d0) then
            out_Ze(i_x,i_y,i_z,i_f,i_p,1) = -9999.d0
          else
            out_Ze(i_x,i_y,i_z,i_f,i_p,1) = 10*log10(Ze_back)
          end if
        if (verbose >= 3) print*, "i_x,i_y,i_z,i_f,out_Ze", i_x,i_y,i_z,i_f,out_Ze(i_x,i_y,i_z,i_f,i_p,1)

      else if ((radar_mode == "moments") .or. (radar_mode == "spectrum")) then
          !calculate the noise level depending on range:
          ! did not find any value in the atmo arrays, take the one from namelist file!
          if (ISNAN(atmo_radar_prop(i_x,i_y,1)) .or. (atmo_radar_prop(i_x,i_y,1) == -9999.)) then
            radar_Pnoise = 10**(0.1*radar_Pnoise0(i_f)) * &
              (out_radar_hgt(i_x,i_y,i_z)/1000.)**2
            if (verbose >= 3) print*, "took radar noise from nml file", 10*log10(radar_Pnoise), &
                radar_Pnoise0(i_f), out_radar_hgt(i_x,i_y,i_z)
          else
            ! take the one from the atmo files
            radar_Pnoise = 10**(0.1*atmo_radar_prop(i_x,i_y,1)) * &
              (out_radar_hgt(i_x,i_y,i_z)/1000.)**2
            if (verbose >= 3) print*, "took radar noise from atmo array", 10*log10(radar_Pnoise), &
                    atmo_radar_prop(i_x,i_y,1), out_radar_hgt(i_x,i_y,i_z)
          end if
          call assert_true(err,(radar_Pnoise > 0),&
              "nan or negative radar_Pnoise")
          if (err > 0) then
            errorstatus = fatal
            msg = "assertation error"
            call report(errorstatus, msg, nameOfRoutine)
            return
          end if


          !get delta velocity
          del_v = (radar_max_V(i_f)-radar_min_V(i_f)) / radar_nfft ![m/s]
          !create array from min_v to max_v iwth del_v spacing -> velocity spectrum of radar
          spectra_velo = (/(((ii*del_v)+radar_min_V(i_f)),ii=0,radar_nfft-1)/) ! [m/s]

          !same for the extended spectrum
          min_V_aliased = radar_min_V(i_f) - radar_aliasing_nyquist_interv*(radar_max_V(i_f)-radar_min_V(i_f))
          max_V_aliased = radar_max_V(i_f) + radar_aliasing_nyquist_interv*(radar_max_V(i_f)-radar_min_V(i_f))
          spectra_velo_aliased = (/(((ii*del_v)+min_V_aliased),ii=0,radar_nfft_aliased-1)/) ! [m/s]


          if (isnan(atmo_airturb(i_x,i_y,i_z)) .and. (.not. isnan(atmo_turb_edr(i_x,i_y,i_z)*atmo_wind_uv(i_x,i_y,i_z))) ) then

            call estimate_spectralBroadening(err,atmo_turb_edr(i_x,i_y,i_z),atmo_wind_uv(i_x,i_y,i_z),out_radar_hgt(i_x,i_y,i_z),&
                radar_fwhr_beamwidth_deg(i_f),radar_integration_time(i_f),wavelength,radar_kolmogorov_constant,ss)
             if (err > 0) then
                errorstatus = fatal
                msg = 'Error in atmo_airturb'
                call report(errorstatus, msg, nameOfRoutine)
                return
             end if

            if (verbose > 10) print*, i_x,i_y,i_z,i_f, ss


          else if (.not.  isnan(atmo_airturb(i_x,i_y,i_z)) .and. &
                  (isnan(atmo_turb_edr(i_x,i_y,i_z)*atmo_wind_uv(i_x,i_y,i_z)))) then
            ss = atmo_airturb(i_x,i_y,i_z)
          else if ((isnan(atmo_airturb(i_x,i_y,i_z))) .and. (isnan(atmo_turb_edr(i_x,i_y,i_z))) .and. &
                 (isnan(atmo_wind_uv(i_x,i_y,i_z))) ) then! no value provided.. take zero
            ss = 0.d0
          else
            errorstatus = fatal
            msg = "Didn't get valid broadening value"
            call report(errorstatus, msg, nameOfRoutine)
            return

          end if  


          ! The convolution function will return a vector of length 2*radar_nfft_aliased-1. 
          ! Get the indices to cut out the required part.
          ! Error here result in a shifted spectrum
          ts_imin = (radar_nfft_aliased/2.) + 1
          ts_imax = 2*(radar_nfft_aliased/2.)

         !get turbulence (no turbulence in clear sky...)
         turb(:) = 0.d0
         if ((ss > 0.d0) .and. (back(i_p) > 0)) then
              do tt = 1  , radar_nfft_aliased
                  ! gaussian function with same length as radar spectrum, centered around zero
                  turb(tt) = exp(-(spectra_velo_aliased(tt)-0)**2.d0/(2.d0*ss**2.d0))
              end do
              turb(:) = turb/ sum(turb) !normalize to unity area   
        end if

        ! skip in case ss was so little that sum(turb) is zero)
        if ((SUM(turb) > 0.d0) .and. (back(i_p) > 0)) then
          call assert_false(err,(ANY(ISNAN(turb)) .or. (SUM(turb) <= 0.d0)),&
              "got nan or negative value in linear turb")
          call assert_true(err,((ABS(turb(1)) < almostZero) .and.(ABS(turb(radar_nfft_aliased))  < almostZero)),&
              "increase radar_aliasing_nyquist_interv, turbulence too large")
          if (err > 0) then
            errorstatus = fatal
            msg = "assertation error"
            call report(errorstatus, msg, nameOfRoutine)
            return
          end if


              !convolute spectrum and noise
              call convolution(err,particle_spectrum_att,radar_nfft_aliased,turb,radar_nfft_aliased,turb_spectra)
              if (err /= 0) then
                  msg = 'error in convolution!'
                  call report(err, msg, nameOfRoutine)
                  errorstatus = err
                  return
              end if

              !I don't like Nans values here
              where(ISNAN(turb_spectra)) turb_spectra = 0.d0
              ! negative number are resulting from numerical effects
              where(turb_spectra<0) turb_spectra = 0.d0
          else
              turb_spectra(:) = 0.d0
              turb_spectra(ts_imin:ts_imax) = particle_spectrum_att

          end if
          if (verbose >= 10)print*, "turb_spectra"
          if (verbose >= 10)print*, SHAPE(turb_spectra)
          if (verbose >= 10)print*, turb_spectra

          quailty_aliasing = 0
          !lets look for aliasing effects. if we calculated radar_aliasing_nyquist_interv for a broader spectrum than necessary, fold it again:
          if (radar_aliasing_nyquist_interv > 0) then
              turb_spectra_aliased(:) = 0.d0
              do, ii=1,1+2*radar_aliasing_nyquist_interv
                  !get indices
                  startI = ts_imin + (ii-1)*radar_nfft
                  stopI =  ts_imax -  (1+2*radar_aliasing_nyquist_interv - ii)*radar_nfft
                  if ((ii .ne. radar_aliasing_nyquist_interv + 1) &
                  .and. (SUM(turb_spectra(startI:stopI)) .gt. almostZero)) then
                      quailty_aliasing = 1
                  end if
                  !appy aliasing
                  turb_spectra_aliased = turb_spectra_aliased + turb_spectra(startI:stopI)
              end do
          else ! aliasing effects not considered
              turb_spectra_aliased = turb_spectra(ts_imin:ts_imax)
          end if

          if (verbose > 2) then
              if (quailty_aliasing .ne. 0) then
                  print*, "radar quality: aliasing found"
              else
                  print*, "radar quality: NO aliasing found"
              end if
          end if

          call assert_false(err,(ANY(ISNAN(turb_spectra_aliased)) .or. ANY(turb_spectra_aliased < 0.d0)),&
              "got nan or negative value in linear turb_spectra_aliased")
          ! call assert_false(err,(ALL(turb_spectra_aliased==0)),&
          !     "all values of turb_spectra_aliased == 0")
          if (err > 0) then
            errorstatus = fatal
            msg = "assertation error"
            call report(errorstatus, msg, nameOfRoutine)
            return
          end if

          if (verbose == 666) then
            print*, "##########################################"
            print*, "velocity (m/s)"
            print*, spectra_velo
            print*, "##########################################"
            print*, "particle_spec with turbulence (v) [mm⁶/m³/(m/s)] ", MAXVAL(turb_spectra_aliased)
            print*, SHAPE(turb_spectra_aliased)
            print*, turb_spectra_aliased
            print*, "##########################################"
          end if
          !spetial output for testing the radar simulator
        out_debug_radarvel(:) = spectra_velo(:)
        out_debug_radarback_wturb(:) = turb_spectra_aliased


          !get the SNR
          SNR = 10.d0*log10(Ze_back/radar_Pnoise)
          !this here is for scaling, if we have now a wrong Ze due to all the turbulence, rescaling etc. 
          !Can happen e.g. due to numeric issues when applying very large or very small turbulence. Skip 
          !if Ze_back ==0.
          if (Ze_back >0) then
            K = (Ze_back/SUM(turb_spectra_aliased*del_v))
          else
            K = 1.d0
          end if

          if (verbose >= 4) print*, "first K", K
          snr_turb_spectra = (K* turb_spectra_aliased + radar_Pnoise/(radar_nfft*del_v))
          !   snr_turb_spectra =turb_spectra_aliased + radar_Pnoise/(radar_nfft*del_v)

          call assert_false(err,(ANY(ISNAN(snr_turb_spectra)) .or. ANY(snr_turb_spectra < 0.d0)),&
              "got nan or negative value in linear snr_turb_spectra")
          call assert_false(err,(ISNAN(K)) .or. (K >= HUGE(K)),&
              "K is nan or infinitive")
          if (err > 0) then
            errorstatus = fatal
            msg = "assertation error"
            call report(errorstatus, msg, nameOfRoutine)
            return
          end if

          if (radar_no_Ave(i_f) .eq. 0) then !0 means infinity-> no noise
              noise_turb_spectra = snr_turb_spectra
          else

              !get noise.
              if (verbose > 2) print*, "get noise"
              if (randomseed == -1) then
                !get it from lat lon
                seed = INT(ABS((atmo_lat(i_x,i_y) * atmo_lon(i_x,i_y) * 1000)))
              else
                seed = randomseed
              end if
              call random(err,radar_no_Ave(i_f)*radar_nfft,seed,x_noise)
              if (err /= 0) then
                  msg = 'error in random!'
                  call report(err, msg, nameOfRoutine)
                  errorstatus = err
                  return
              end if
              do tt = 1, radar_no_Ave(i_f)
                  noise_turb_spectra_tmp(tt,:) = -log(x_noise((tt-1)*radar_nfft+1:tt*radar_nfft))*snr_turb_spectra
              end do

              if (radar_no_Ave(i_f) .eq. 1) then
                  noise_turb_spectra = noise_turb_spectra_tmp(1,:)
              else
                  noise_turb_spectra = SUM(noise_turb_spectra_tmp,DIM=1)/radar_no_Ave(i_f)
              end if
          end if


          !spetial output for testing the radar simulator
          if (verbose == 666) then
            print*, "##########################################"
            print*, "particle_spec with turbulence and noise (v) [mm⁶/m³/(m/s)] ", MAXVAL(noise_turb_spectra)
            print*, noise_turb_spectra
            print*, "##########################################"
          end if
          out_debug_radarback_wturb_wnoise(:) = noise_turb_spectra(:)

          !apply spectral resolution
          noise_turb_spectra = noise_turb_spectra * del_v !now [mm⁶/m³]

          if (verbose >= 4) then
              print*,"first K",K
              print*,"TOTAL"," Ze back",10*log10(Ze_back)
              print*,"TOTAL"," Ze SUM(particle_spectrum)*del_v",10*log10(SUM(particle_spectrum)*del_v)
              print*,"TOTAL"," Ze SUM(particle_spectrum_att)*del_v",10*log10(SUM(particle_spectrum_att)*del_v)
              print*,"TOTAL"," Ze SUM(turb_spectra)*del_v",10*log10(SUM(turb_spectra)*del_v)
              print*,"TOTAL"," Ze SUM(turb_spectra_aliased)*del_v",10*log10(SUM(turb_spectra_aliased)*del_v)
              print*,"TOTAL"," Ze SUM(snr_turb_spectra)*del_v",10*log10(SUM(snr_turb_spectra)*del_v)
              print*,"TOTAL"," Ze SUM(noise_turb_spectra)*del_v",10*log10(SUM(noise_turb_spectra))
              print*,"TOTAL"," Ze SUM(snr_turb_spectra)*del_v-radar_Pnoise",10*log10(SUM(snr_turb_spectra)*del_v-radar_Pnoise)
              print*,"TOTAL"," Ze SUM(noise_turb_spectra)*del_v-radar_Pnoise",10*log10(SUM(noise_turb_spectra)-radar_Pnoise)
          end if

          if (verbose >= 5) then
              print*," linear spectrum befor receiver uncertainty"
              print*,noise_turb_spectra
              print*,"#####################"
          end if

            !apply a receiver uncertainty:
          if (radar_receiver_uncertainty_std(i_f) /= 0) then
            !get random
            call random(err,2,seed,rand_number)
            if (err /= 0) then
                msg = 'error in random!'
                call report(err, msg, nameOfRoutine)
                errorstatus = err
                return
            end if
            !apply a gaussian distribution to random numbers
            receiver_uncertainty =  radar_receiver_uncertainty_std(i_f) * sqrt( -2.0d0 * log ( rand_number(1))) &
                          * cos(2.0d0 * pi * rand_number(2))
            !make linear
            receiver_uncertainty = 10**(0.1*receiver_uncertainty)
            !apply to spectrum
            noise_turb_spectra = noise_turb_spectra * receiver_uncertainty
            radar_Pnoise = radar_Pnoise * receiver_uncertainty
          end if

          !apply a receiver miscalibration:
          noise_turb_spectra = noise_turb_spectra * 10**(0.1*radar_receiver_miscalibration(i_f))
          radar_Pnoise = radar_Pnoise * 10**(0.1*radar_receiver_miscalibration(i_f))
          

          !calculate noise level (actually we know already the result which is noise_model)
          if (radar_use_hildebrand) then
              call radar_hildebrand_sekhon(err,noise_turb_spectra,radar_no_Ave(i_f),radar_nfft,&
                noise,noise_max)
              if (err /= 0) then
                  msg = 'error in radar_hildebrand_sekhon!'
                  call report(err, msg, nameOfRoutine)
                  errorstatus = err
                  return
              end if   
              if (verbose .ge. 3) print*, i_f, 'calculated noise, noise_max:', noise, noise_max
              if (radar_noise_distance_factor(i_f) > 0) &
                noise_max = radar_noise_distance_factor(i_f)*noise
          else
              noise = radar_Pnoise/radar_nfft !no devison by del_v neccessary!
              noise_max = radar_noise_distance_factor(i_f)*noise
          end if          

          !noise for the output
          noise_out = noise * radar_nfft

          call radar_calc_moments(err,radar_nfft,radar_nPeaks,&
            noise_turb_spectra,noise,noise_max, &
            noise_removed_turb_spectra,&
            moments,slope,edge,quality_moments)
          if (err /= 0) then
            msg = 'error in radar_calc_moments!'
            call report(err, msg, nameOfRoutine)
            errorstatus = err
            return
        end if
          if (verbose >= 4) then
            do i_n = 1  , radar_nPeaks
              print*,"TOTAL#",i_n," Ze moments log ",10*log10(moments(0,i_n)), "lin ", moments(0,i_n)
              print*,"#####################"
            end do
          end if

            ! collect results for output
          !   out_radar_spectra(i_x,i_y,i_z,fi,:) = 10*log10(particle_spectrum_att(513:1024))

          !if wanted, apply the noise correction to the spectrum to be saved.
          if (radar_save_noise_corrected_spectra) noise_turb_spectra = noise_removed_turb_spectra

          if (verbose >= 5) then
              print*,"final linear spectrum"
              print*,noise_turb_spectra
              print*,"#####################"
          end if


          out_radar_spectra(i_x,i_y,i_z,i_f,i_p,:) = 10*log10(noise_turb_spectra)
          WHERE (ISNAN(noise_turb_spectra)) noise_turb_spectra = -9999.d0
          if (verbose >= 5 ) then
              print*,"final log spectrum"
              print*, "i_x,i_y,i_z,i_f,i_p,", i_x,i_y,i_z,i_f,i_p
              print*,out_radar_spectra(i_x,i_y,i_z,i_f,i_p,:)
              print*,"#####################"
          end if
          do i_n = 1  , radar_nPeaks
            out_radar_snr(i_x,i_y,i_z,i_f,i_p,i_n) = SNR !same SNR for all values as of now...
            out_radar_vel(i_f,:) = spectra_velo(:)
            out_radar_moments(i_x,i_y,i_z,i_f,i_p,i_n,:) = moments(1:4,i_n)
            out_radar_slopes(i_x,i_y,i_z,i_f,i_p,i_n,:) = slope(:,i_n)
            out_radar_edges(i_x,i_y,i_z,i_f,i_p,i_n,:) = edge(:,i_n)
            out_radar_quality(i_x,i_y,i_z,i_f,i_p,i_n) = quailty_aliasing + quality_moments !same quailty for all values as of now...

            moments(0,i_n) = 10*log10(moments(0,i_n))
            IF (ISNAN(moments(0,i_n))) moments(0,i_n) = -9999.d0
            out_Ze(i_x,i_y,i_z,i_f,i_p,i_n) = moments(0,i_n)
            if (verbose >= 5 ) then
              print*, "i_x,i_y,i_z,i_f,i_p,i_n ",i_x,i_y,i_z,i_f,i_p,i_n
              print*, "out_radar_snr",out_radar_snr(i_x,i_y,i_z,i_f,i_p,i_n)
              print*, "out_radar_moments",out_radar_moments(i_x,i_y,i_z,i_f,i_p,i_n,:)
              print*, "out_radar_slopes",out_radar_slopes(i_x,i_y,i_z,i_f,i_p,i_n,:)
              print*, "out_radar_edges",out_radar_edges(i_x,i_y,i_z,i_f,i_p,i_n,:)
              print*, "out_radar_quality",out_radar_quality(i_x,i_y,i_z,i_f,i_p,i_n)
              print*, "out_Ze",out_Ze(i_x,i_y,i_z,i_f,i_p,i_n)
            end if
          end do
      else
        errorstatus = fatal
        msg =   "did not understand radar_mode"// radar_mode
        call report(errorstatus, msg, nameOfRoutine)
        return
      end if

    end do !radar_npol

    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return
end subroutine radar_simulator







