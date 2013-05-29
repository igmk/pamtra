subroutine radar_calc_moments(errorstatus,radar_spectrum_in,radar_spectrum_out,moments,slope,quality)

    ! written by P. Kollias, tranlated to Fortran by M. Maahn (12.2012)
    ! calculate the 0th -4th moment and the slopes of the peak of a radar spectrum!
    !
    ! in
    ! radar_spectrum_in: radar spectrum with noise [mm⁶/m³]
    ! out
    ! radar_spectrum_out: radar spectrum with noise removed [mm⁶/m³]
    ! moments, dimension(0:4):0th - 4th moment [mm⁶/m³, m/s, m/s,-,-]
    ! slope, dimension(2): left(0) and right(1) slope of the peak [dB/(m/s)]


    use kinds
    use settings
    use constants
    use report_module
    implicit none


    real(kind=dbl), dimension(radar_nfft), intent(in):: radar_spectrum_in
    real(kind=dbl), dimension(radar_nfft), intent(out):: radar_spectrum_out
    real(kind=dbl), dimension(0:4), intent(out):: moments
    real(kind=dbl), dimension(2), intent(out):: slope
    integer, intent(out) :: quality
    real(kind=dbl) :: del_v, noise
    integer :: spec_max, spec_max_a(1), right_edge,left_edge, &
    spec_max_cp, right_edge_cp, left_edge_cp, ii, jj
    real(kind=dbl), dimension(radar_nfft) :: spectra_velo, specLog,&
    radar_spectrum_cp, radar_spectrum_smooth

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'convolution'     
    
    interface
        SUBROUTINE SMOOTH_SAVITZKY_GOLAY(errorstatus,dataIn,length,dataOut)
            use kinds
            implicit none
	    integer(kind=long), intent(out) :: errorstatus
            integer, intent(in) :: length
            real(kind=dbl), intent(in), dimension(length) :: dataIn
            real(kind=dbl), intent(out), dimension(length) :: dataOut
        END SUBROUTINE SMOOTH_SAVITZKY_GOLAY

        subroutine radar_hildebrand_sekhon(errorstatus,spectrum,n_ave,n_ffts,noise_mean)
            use kinds
            implicit none
	    integer(kind=long), intent(out) :: errorstatus
            integer, intent(in) :: n_ave, n_ffts
            real(kind=dbl), dimension(n_ffts), intent(in) :: spectrum
            real(kind=dbl), intent(out) :: noise_mean
        end subroutine radar_hildebrand_sekhon
    end interface

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    
    del_v = (radar_max_V-radar_min_V) / radar_nfft

    spectra_velo = (/(((ii*del_v)+radar_min_V),ii=0,radar_nfft-1)/) ! [m/s]

    !find maximum of spectrum -> most significant peak
    spec_max_a = MAXLOC(radar_spectrum_in)
    if ((spec_max_a(1) == 1) .or. (spec_max_a(1) == radar_nfft)) then
        spec_max_a = MAXLOC(radar_spectrum_in(2:radar_nfft-1))
    end if
    spec_max = spec_max_a(1)

    !calculate noise level (actually we know already the result which is radar_pnoise)
    if (radar_use_hildebrand) then
        call radar_hildebrand_sekhon(err,radar_spectrum_in,radar_no_Ave,radar_nfft,noise)
	if (err /= 0) then
	    msg = 'error in radar_hildebrand_sekhon!'
	    call report(err, msg, nameOfRoutine)
	    errorstatus = err
	    return
	end if   
        if (verbose .gt. 2) print*, 'calculated noise:', noise
    else
        noise = radar_pnoise/radar_nfft !no devison by del_v neccessary!
    end if

    !remove noise
    radar_spectrum_out = radar_spectrum_in - noise ! mm⁶/m³

    !!get the borders of the most significant peak
    do ii = spec_max+1, radar_nfft
        if (radar_spectrum_out(ii) <= 0.25*noise ) EXIT
    end do
    right_edge = ii
    do jj = spec_max-1, 1, -1
        if (radar_spectrum_out(jj) <= 0.25*noise ) EXIT
    end do
    left_edge = jj

    !check whether peak is present:
    if (SUM(radar_spectrum_in(left_edge+1:right_edge-1))/(right_edge-left_edge-1)/noise <radar_min_spectral_snr) then
        !no peak
        radar_spectrum_out = -9999
        moments = -9999
        slope = -9999
        quality = 64
  
    else
        !!look for a second peak for quality array
        radar_spectrum_cp = radar_spectrum_out
        !set peak to zero
        radar_spectrum_cp(left_edge+1:right_edge-1) = 0.d0
        !find maximum of spectrum -> most significant peak
        spec_max_a = MAXLOC(radar_spectrum_cp)
        if ((spec_max_a(1) == 1) .or. (spec_max_a(1) == radar_nfft)) then
            spec_max_a = MAXLOC(radar_spectrum_cp(2:radar_nfft-1))
        end if
        spec_max_cp = spec_max_a(1)
        !find second peak
        do ii = spec_max_cp+1, radar_nfft
            if (radar_spectrum_cp(ii) <= 0 ) EXIT
        end do
        right_edge_cp = ii
        do jj = spec_max_cp-1, 1, -1
            if (radar_spectrum_cp(jj) <= 0 ) EXIT
        end do
        left_edge_cp = jj

        !     fstPeak = SUM(radar_spectrum_out(left_edge:right_edge))
        !     secPeak = SUM(radar_spectrum_out(left_edge_cp:right_edge_cp))

        !second peak needs to fullfill: width min 0.01*radar_nfft+3, power min 1% of primary
        !     if ((right_edge_cp - left_edge_cp > 0.01*radar_nfft+3) .and. (secPeak > 0.01d0*fstPeak)) then
        if (SUM(radar_spectrum_in(left_edge_cp+1:right_edge_cp-1))/&
        (right_edge_cp-left_edge_cp-1)/noise >radar_min_spectral_snr) then
            quality = 2
            if (verbose > 2) print*, "second peak found", &
            SUM(radar_spectrum_in(left_edge_cp+1:right_edge_cp-1))/(right_edge_cp-left_edge_cp-1)/noise
        else
            quality = 0
            if (verbose > 2) print*, "NO second peak found", &
            SUM(radar_spectrum_in(left_edge_cp+1:right_edge_cp-1))/(right_edge_cp-left_edge_cp-1)/noise
        end if

        !make the spectrum smooth
        call smooth_savitzky_golay(err,radar_spectrum_in, radar_nfft, radar_spectrum_smooth)
	if (err /= 0) then
	    msg = 'error in smooth_savitzky_golay!'
	    call report(err, msg, nameOfRoutine)
	    errorstatus = err
	    return
	end if  
        !set remaining sectrum to zero

        radar_spectrum_out(1:left_edge) = 0.d0
        radar_spectrum_out(right_edge:radar_nfft) = 0.d0
        radar_spectrum_smooth(1:left_edge) = 0.d0
        radar_spectrum_smooth(right_edge:radar_nfft) = 0.d0

        !get the (log)slope of the peak
        specLog = 10*log10(radar_spectrum_smooth)
        slope(:) = 0.d0 ! dB/(m/s)
        slope(1) = (specLog(spec_max)-specLog(left_edge+1))/(spectra_velo(spec_max)-spectra_velo(left_edge+1))
        slope(2) = (specLog(right_edge-1)-specLog(spec_max))/(spectra_velo(right_edge-1)-spectra_velo(spec_max))

        WHERE (ISNAN(slope)) slope = -9999.d0

        !calculate the moments
        moments(0) = SUM(radar_spectrum_smooth) ! mm⁶/m³
        moments(1) = SUM(radar_spectrum_smooth*spectra_velo)/moments(0) ! m/s
        moments(2) = SQRT(SUM(radar_spectrum_smooth * (spectra_velo-moments(1))**2)/moments(0)) ! m/s
        moments(3) = SUM(radar_spectrum_smooth * (spectra_velo-moments(1))**3)/(moments(0)*moments(2)**3) ![-]
        moments(4) = SUM(radar_spectrum_smooth * (spectra_velo-moments(1))**4)/(moments(0)*moments(2)**4) ![-]

    end if

    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return

end subroutine radar_calc_moments
