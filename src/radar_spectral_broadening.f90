module radar_spectral_broadening

  contains

  subroutine estimate_spectralBroadening(EDR,wind_uv,height,beamwidth_deg,integration_time,wavelength,kolmogorov,specBroad)
   
  ! everything in SI units
  ! estimate spectral broadening of radar Doppler spectrum due to turbulence and horizontal wind. Wind gradient is neglected

    use kinds
    use constants, only: pi
    use report_module

    implicit none

    real(kind=dbl), intent(in) :: EDR
    real(kind=dbl), intent(in) :: wind_uv
    real(kind=dbl), intent(in) :: height
    real(kind=dbl), intent(in) :: beamwidth_deg !full width half radiation
    real(kind=dbl), intent(in) :: integration_time
    real(kind=dbl), intent(in) :: wavelength
    real(kind=dbl), intent(in) :: kolmogorov
    real(kind=dbl), intent(out) :: specBroad

    real(kind=dbl) :: L_s, L_lambda, sig_B, sig_T, beamwidth

      beamwidth = beamwidth_deg/2./180.*pi 

      L_s = (wind_uv * integration_time) + 2*height*sin(beamwidth)!CORRECT, formular of shupe or oconnor is for full width beam width
      L_lambda = wavelength / 2.
      sig_B = sqrt(wind_uv**2 * (beamwidth)**2 / 2.76)
      sig_T = sqrt(3*kolmogorov/2. * (EDR/(2.*pi))**(2./3.) * (L_s**(2./3.) - L_lambda**(2./3.)))
      specBroad = sqrt(sig_B**2 + sig_T**2)


  end subroutine estimate_spectralBroadening
end module radar_spectral_broadening
