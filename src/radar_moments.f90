module radar_moments

  use kinds
  use settings, only: &
    radar_max_V,&
    radar_min_V,&
    radar_no_Ave,&
    radar_smooth_spectrum, &
    radar_use_wider_peak, &
    radar_peak_min_bins, &
    radar_peak_min_snr,&
    radar_peak_snr_definition
  use constants
  use report_module
  implicit none

  contains 
subroutine radar_calc_moments(errorstatus,radar_nfft,radar_nPeaks,radar_spectrum_in,&
  noise,noise_max, &
  radar_spectrum_out,moments,slope,edge,quality)

    ! written by P. Kollias, translated to Fortran by M. Maahn (12.2012)
    ! calculate the 0th -4th moment and the slopes of the peak of a radar spectrum!
    !
    ! in
    ! radar_spectrum_in: radar spectrum with noise [mm⁶/m³]
    ! noise: mean spectral noise in [mm⁶/m³]
    ! noise_max: max spectral noise in [mm⁶/m³]
    ! out
    ! radar_spectrum_out: radar spectrum with noise removed [mm⁶/m³]
    ! moments, dimension(0:4):0th - 4th moment [mm⁶/m³, m/s, m/s,-,-]
    ! slope, dimension(2): left(0) and right(1) slope of the peak [dB/(m/s)]
    ! edge, dimension(2): left(0) and right(1) edge the peak [m/s]
    ! quality

    use vars_index, only: i_f
    implicit none

    integer, intent(in) :: radar_nfft, radar_nPeaks
    real(kind=dbl), dimension(radar_nfft), intent(in):: radar_spectrum_in
    real(kind=dbl), intent(in):: noise, noise_max
    real(kind=dbl), dimension(radar_nfft), intent(out):: radar_spectrum_out
    real(kind=dbl), dimension(0:4,radar_nPeaks), intent(out):: moments
    real(kind=dbl), dimension(2,radar_nPeaks), intent(out):: slope
    real(kind=dbl), dimension(2,radar_nPeaks), intent(out):: edge
    integer, intent(out) :: quality
    real(kind=dbl), dimension(radar_nPeaks+1,radar_nfft):: radar_spectrum_arr
    real(kind=dbl), dimension(radar_nfft):: radar_spectrum_only_noise
    real(kind=dbl), dimension(radar_nfft):: radar_spectrum_4mom
    real(kind=dbl) :: del_v, specMax, noiselog, specSNR
    integer :: spec_max_ii, spec_max_ii_a(1), right_edge,left_edge, &
    ii, jj, right_edge4slope, &
    left_edge4slope
    real(kind=dbl), dimension(radar_nfft) :: spectra_velo, specLog
    logical :: additionalPeaks
    integer :: nn, kk

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'radar_calc_moments'     
    
    interface
        SUBROUTINE SMOOTH_SAVITZKY_GOLAY(errorstatus,dataIn,length,dataOut)
            use kinds
            implicit none
            integer(kind=long), intent(out) :: errorstatus
            integer, intent(in) :: length
            real(kind=dbl), intent(in), dimension(length) :: dataIn
            real(kind=dbl), intent(out), dimension(length) :: dataOut
        END SUBROUTINE SMOOTH_SAVITZKY_GOLAY

    end interface

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    
    if (verbose >= 10) print*, "radar_nfft,radar_nPeaks,noise, noise_max,radar_spectrum_in"
    if (verbose >= 10) print*, radar_nfft,radar_nPeaks,noise,noise_max,"spec:",radar_spectrum_in
 
    err = 0

      call assert_false(err,any(ISNAN(radar_spectrum_in)),&
          "found nan in radar_spectrum")
      if (err /= 0) then
        msg = 'assertation error'
        call report(err, msg, nameOfRoutine)
        errorstatus = err
        return
      end if


    !initilaize
    moments(:,:) = -9999
    edge(:,:) = -9999
    slope(:,:) = -9999


    del_v = (radar_max_V(i_f)-radar_min_V(i_f)) / radar_nfft
    spectra_velo = (/(((ii*del_v)+radar_min_V(i_f)),ii=0,radar_nfft-1)/) ! [m/s]
    quality = 0
    additionalPeaks = .false.


  do nn = 1, radar_nPeaks + 1
    radar_spectrum_arr(nn,:) = radar_spectrum_in
  end do

  radar_spectrum_only_noise = radar_spectrum_in


  do nn = 1, radar_nPeaks + 1

    !find maximum of spectrum (which must not be at the borders!) -> most significant peak
    spec_max_ii_a = MAXLOC(radar_spectrum_arr(nn,2:radar_nfft-1))
    spec_max_ii = spec_max_ii_a(1) +1 
    if (verbose >= 6) print*, "found maximum at ", spec_max_ii, radar_spectrum_arr(nn,spec_max_ii)

    if (radar_spectrum_arr(nn,spec_max_ii) < noise_max) then
      if (verbose >= 3) print*, "Skipped peak ", nn, " because of:", &
        "spec_max_ii) < noise_max",radar_spectrum_arr(nn,spec_max_ii) < noise_max, &
        radar_spectrum_arr(nn,spec_max_ii), noise_max
      if (nn == 1) then 
        radar_spectrum_out(:) = -9999 !right now, only primary peak is processed for radar_spectrum_out...
        quality = quality + 64 !no peak found at all
      end if
      EXIT !loop. no more peaks
    end if


    !!get the borders of the most significant peak
    do ii = spec_max_ii+1, radar_nfft
        if (verbose >= 6) print*, "to the right:",nn,ii, radar_spectrum_arr(nn,ii), noise_max, spectra_velo(ii), &
            radar_spectrum_arr(nn,ii) <= noise_max
        if (radar_spectrum_arr(nn,ii) <= noise_max ) EXIT
    end do
    ! if (ii > radar_nfft) ii = radar_nfft !Fortran tends to go one step too far if EXIT does not happen
    right_edge = ii
    right_edge4slope = right_edge
    if ((radar_use_wider_peak) .and. &
        (radar_spectrum_arr(nn,right_edge) >= noise ) .and. &
         (right_edge < radar_nfft) ) then
      right_edge = right_edge +1
      if (verbose >= 6) print*, nn, "extended right edge to ", right_edge
    end if
    do jj = spec_max_ii-1, 1, -1
        if (verbose >= 6) print*, "to the left:",nn,jj, radar_spectrum_arr(nn,jj), noise_max, spectra_velo(jj), &
            radar_spectrum_arr(nn,jj) <= noise_max
        if (radar_spectrum_arr(nn,jj) <= noise_max ) EXIT
    end do
    ! if (jj < 1) jj = 1 !Fortran tends to go one step too far if EXIT does not happen
    left_edge = jj
    left_edge4slope = left_edge
    if ((radar_use_wider_peak) .and. &
        (radar_spectrum_arr(nn,left_edge) >= noise )  .and. &
         (left_edge > 1) ) then

      left_edge = left_edge -1
      if (verbose >= 6) print*, nn, "extended left edge to ", left_edge
    end if

    if (verbose >= 5)  print*, nn, "found peak from", left_edge, radar_spectrum_arr(nn,left_edge+1),&
      "to", right_edge-1, radar_spectrum_arr(nn,right_edge-1), "noiseMax", noise_max

    !we don't want to find this peak again, so set it to zero for the other spectra
    do kk=nn+1, radar_nPeaks+1
      radar_spectrum_arr(kk,left_edge+1:right_edge-1) = 0.d0
    end do 
    
    if (radar_peak_snr_definition == 'log') then
      specSNR = 10 * log10 (SUM(radar_spectrum_arr(nn,left_edge+1:right_edge-1))/(noise*radar_nfft) )
    else if (radar_peak_snr_definition == 'specLin') then !old defintiion kept for backwards compatibility
      specSNR = SUM(radar_spectrum_arr(nn,left_edge+1:right_edge-1))/(right_edge-left_edge-1)/noise
    else
            msg = 'Do not understand radar_peak_snr_definition'
            call report(fatal, msg, nameOfRoutine)
            errorstatus = fatal
            return
    end if
    !check whether peak is NOT present:
    if ((specSNR <radar_peak_min_snr(i_f)) &
      .or. (right_edge-left_edge <= radar_peak_min_bins(i_f))  .or. (right_edge-left_edge ==  radar_nfft+ 1) ) then

      !no or too thin peak or too wide peak
      if (verbose >= 3) print*, "Skipped peak ", nn, " because of:", &
        "radar_peak_min_snr(i_f)", specSNR <radar_peak_min_snr(i_f), specSNR, radar_peak_min_snr(i_f), &
        "too thin", (right_edge-left_edge <= radar_peak_min_bins(i_f)), &
        "too wide", (right_edge-left_edge ==  radar_nfft+ 1)

      if (nn == 1) then 
        radar_spectrum_out(:) = -9999 !right now, only primary peak is processed for radar_spectrum_out...
        quality = quality + 64 !no peak found at all
      end if
      EXIT !loop. no more peaks

    else !peak is present!
      if (verbose >= 3) print*, nn, "peak ", nn, " confirmed with spec SNR of", specSNR

      radar_spectrum_only_noise(left_edge+1:right_edge-1) = -9999 ! in this spectrum we want ALL peaks removed

      if (nn > radar_nPeaks) then 
        additionalPeaks = .true.
        !additional peaks are as of now NOT processed...
        EXIT !loop. no more peaks
      end if

      if (radar_smooth_spectrum) then
        !make the spectrum smooth and remove noise
        call smooth_savitzky_golay(err,radar_spectrum_in, radar_nfft, radar_spectrum_4mom)
        if (err /= 0) then
            msg = 'error in smooth_savitzky_golay!'
            call report(err, msg, nameOfRoutine)
            errorstatus = err
            return
        end if  
      else
        radar_spectrum_4mom(:) = radar_spectrum_in(:)
      end if

!     remove noise for moment estimation
      radar_spectrum_4mom = radar_spectrum_4mom - noise

!     for moments estimation, set remaining sectrum to zero
      if (left_edge  >= 1) radar_spectrum_4mom(1:left_edge) = 0.d0
      if (right_edge <= radar_nfft) radar_spectrum_4mom(right_edge:radar_nfft) = 0.d0

      if (verbose >= 5) print*, "radar_spectrum_smooth, peak only", SHAPE(radar_spectrum_4mom),  radar_spectrum_4mom
      if (verbose >= 5) print*, "spectra_velo", SHAPE(spectra_velo), spectra_velo

!     calculate the moments
      moments(0,nn) = SUM(radar_spectrum_4mom) ! mm⁶/m³
      moments(1,nn) = SUM(radar_spectrum_4mom*spectra_velo)/moments(0,nn) ! m/s
      moments(2,nn) = SQRT(SUM(radar_spectrum_4mom * (spectra_velo-moments(1,nn))**2)/moments(0,nn)) ! m/s
      moments(3,nn) = SUM(radar_spectrum_4mom * (spectra_velo-moments(1,nn))**3)/(moments(0,nn)*moments(2,nn)**3) ![-]
      moments(4,nn) = SUM(radar_spectrum_4mom * (spectra_velo-moments(1,nn))**4)/(moments(0,nn)*moments(2,nn)**4) ![-]

      edge(1,nn) = spectra_velo(left_edge+1)
      edge(2,nn) = spectra_velo(right_edge-1)      

      if (nn == 1) then 
  !     output spectrum is with main peak only and noise removed but without the smoothing
        radar_spectrum_out = radar_spectrum_in - noise
        if (left_edge  >= 1) radar_spectrum_out(1:left_edge) = 0.d0
        if (right_edge <= radar_nfft) radar_spectrum_out(right_edge:radar_nfft) = 0.d0
      end if

    if (moments(1,nn) == -9999) then
      slope(:,nn) = -9999
    !sometimes the estimation of noise goes wrong
    else if ((MAXVAL(radar_spectrum_only_noise) - noise) <= 0) then
      slope(:,nn) = -9999
      moments(:,nn) = -9999
      edge(:,nn) = -9999
      quality = quality + 128 !error in noise estiamtion
    else 
      if (verbose >= 5) print*,MAXVAL(radar_spectrum_only_noise), noise
    end if !moments(1,1) == -9999


    noiselog = 10*log10(noise)
    specMax = 10*log10(MAXVAL(radar_spectrum_4mom)) !without any noise removed!

      call assert_false(err,nn>radar_nPeaks,&
          "nn>radar_nPeaks")
      call assert_false(err,spec_max_ii == left_edge,&
          "spectra_velo(spec_max_ii) == spectra_velo(left_edge)")
      call assert_false(err,right_edge == spec_max_ii,&
          "spectra_velo(right_edge) == spectra_velo(spec_max_ii)")
      if (err /= 0) then
        msg = 'assertation error'
        call report(err, msg, nameOfRoutine)
        errorstatus = err
        return
      end if

      slope(:,nn) = 0.d0 ! dB/(m/s)
      if (left_edge4slope >= 1) then
        slope(1,nn) = (specMax-noiselog)/&
                    (spectra_velo(spec_max_ii)-spectra_velo(left_edge4slope))
      else
        slope(1,nn) = -9999
      end if
      if (right_edge4slope <= radar_nfft) then
        slope(2,nn) = (noiselog-specMax)/&
                    (spectra_velo(right_edge4slope)-spectra_velo(spec_max_ii))
      else
        slope(2,nn) = -9999
      end if


    if (verbose >= 5) print*, "nn", nn
    if (verbose >= 5) print*, "left slope",specMax, noiselog, spectra_velo(spec_max_ii),&
      spectra_velo(left_edge4slope), slope(1,nn)
    if (verbose >= 5) print*, "right slope",noiselog,specMax, spectra_velo(right_edge4slope),&
      spectra_velo(spec_max_ii), slope(2,nn)
    if (verbose >= 5) print*, "left_edge", left_edge
    if (verbose >= 5) print*, "right_edge", right_edge
    if (verbose >= 5) print*, "spec_max_ii", spec_max_ii
    if (verbose >= 5) print*, "quality", quality
    if (verbose >= 5) print*, "edge", edge(:,nn)
    if (verbose >= 5) print*, "moments", moments

    call assert_false(err,slope(1,nn) >= HUGE(slope(1,nn)),&
        "inf in left slope")
    call assert_false(err,slope(2,nn) >= HUGE(slope(2,nn)),&
        "inf in right slope")

    end if !skip peak
  end do !no peaks

    if (additionalPeaks) quality = quality + 2

    if (verbose >= 5) print*, "radar_nfft", radar_nfft
    if (verbose >= 5) print*, "radar_spectrum_out", radar_spectrum_out



    call assert_false(err,any(ISNAN(slope)),&
        "nan in slope")

    call assert_false(err,any(ISNAN(moments)),&
        "nan in moments")
    call assert_false(err,any(ISNAN(radar_spectrum_out)),&
        "nan in radar_spectrum_out")   
    if (err /= 0) then
      msg = 'assertation error'
      call report(err, msg, nameOfRoutine)
      errorstatus = err

! print*, noiselog, specMax, spectra_velo(left_edge4slope), spectra_velo(right_edge4slope), spectra_velo(spec_max_ii)
! print*, 'radar_spectrum_4mom', radar_spectrum_4mom

      return
    end if


    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return

end subroutine radar_calc_moments
end module radar_moments
