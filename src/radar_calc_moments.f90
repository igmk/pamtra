subroutine radar_calc_moments(radar_spectrum_in,radar_spectrum_out,moments,slope)

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
  use nml_params
  use constants
  implicit none


  real(kind=dbl), dimension(radar_nfft), intent(in):: radar_spectrum_in
  real(kind=dbl), dimension(radar_nfft), intent(out):: radar_spectrum_out
  real(kind=dbl), dimension(0:4), intent(out):: moments  
  real(kind=dbl), dimension(2), intent(out):: slope  
  
  real(kind=dbl) :: del_v, noise
  integer :: spec_max, spec_max_a(1), right_edge, left_edge, ii, jj
  real(kind=dbl), dimension(radar_nfft) :: spectra_velo, specLog

  if (verbose .gt. 1) print*, 'Entering radar_moments.f90'

  del_v = (radar_max_V-radar_min_V) / radar_nfft

  spectra_velo = (/(((ii*del_v)+radar_min_V),ii=0,radar_nfft)/) ! [m/s]

  !make the spectrum smooth
  call smooth_savitzky_golay(radar_spectrum_in, radar_nfft, radar_spectrum_out)

  !find maximum of spectrum -> most significant peak
  spec_max_a = MAXLOC(radar_spectrum_out)
  if ((spec_max_a(1) == 1) .or. (spec_max_a(1) == radar_nfft)) then
    spec_max_a = MAXLOC(radar_spectrum_out(2:radar_nfft-1))
  end if
  spec_max = spec_max_a(1)

  !calculate noise level (actually we know already the result which is radar_pnoise)
  if (radar_use_hildebrand) then
    call radar_hildebrand_sekhon(radar_spectrum_out,radar_no_Ave,radar_nfft,noise)
      if (verbose .gt. 2) print*, 'calculated noise:', noise
  else
    noise = radar_pnoise/radar_nfft !no devison by del_v neccessary!
  end if

  !!get the borders of the most significant peak
  do ii = spec_max+1, radar_nfft
    if (radar_spectrum_out(ii) - noise <= 0 ) EXIT
  end do
  right_edge = ii
  do jj = spec_max-1, 1, -1
    if (radar_spectrum_out(jj) - noise <= 0 ) EXIT
  end do
  left_edge = jj


  !remove noise and set remaining sectrum to zero
  radar_spectrum_out = radar_spectrum_out - noise ! mm⁶/m³
  radar_spectrum_out(1:left_edge) = 0.d0
  radar_spectrum_out(right_edge:radar_nfft) = 0.d0

  !get the (log)slope of the peak
  specLog = 10*log10(radar_spectrum_out)
  slope(:) = 0.d0 ! dB/(m/s)
  slope(1) = (specLog(spec_max)-specLog(left_edge+1))/(spectra_velo(spec_max)-spectra_velo(left_edge+1))
  slope(2) = (specLog(right_edge-1)-specLog(spec_max))/(spectra_velo(right_edge-1)-spectra_velo(spec_max))


  !calculate the moments
  moments(0) = SUM(radar_spectrum_out) ! mm⁶/m³
  moments(1) = SUM(radar_spectrum_out*spectra_velo)/moments(0) ! m/s
  moments(2) = SQRT(SUM(radar_spectrum_out * (spectra_velo-moments(1))**2)/moments(0)) ! m/s
  moments(3) = SUM(radar_spectrum_out * (spectra_velo-moments(1))**3)/(moments(0)*moments(2)**3) ![-]
  moments(4) = SUM(radar_spectrum_out * (spectra_velo-moments(1))**4)/(moments(0)*moments(2)**4) ![-]

  if (verbose .gt. 1) print*, 'Exiting radar_moments.f90'


end subroutine radar_calc_moments