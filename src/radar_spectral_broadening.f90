module radar_spectral_broadening

  contains

  subroutine estimate_spectralBroadening(errorstatus,&
    EDR,&
    wind_uv,&
    height,&
    beamwidth_deg,&
    integration_time,&
    wavelength,&
    kolmogorov,&
    specBroad)
   
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

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err
    character(len=80) :: msg
    character(len=15) :: nameOfRoutine = 'estimate_spectralBroadening'

    real(kind=dbl) :: L_s, L_lambda, sig_B, sig_T, beamwidth

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
    err = 0

    call assert_false(err,(ISNAN(EDR)),&
        "got nan in values in EDR")
    call assert_false(err,(ISNAN(wind_uv)),&
        "got nan in values in wind_uv")
    call assert_false(err,(ISNAN(height)),&
        "got nan in values in height")
    call assert_false(err,(ISNAN(beamwidth_deg)),&
        "got nan in values in beamwidth_deg")
    call assert_false(err,(ISNAN(integration_time)),&
        "got nan in values in integration_time")
    call assert_false(err,(ISNAN(wavelength)),&
        "got nan in values in wavelength")
    call assert_false(err,(ISNAN(kolmogorov)),&
        "got nan in values in kolmogorov")
    if (err > 0) then
      errorstatus = fatal
      msg = "assertation error"
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if


      beamwidth = beamwidth_deg/2./180.*pi 

      L_s = (wind_uv * integration_time) + 2*height*sin(beamwidth)!CORRECT, formular of shupe or oconnor is for full width beam width
      L_lambda = wavelength / 2.
      sig_B = sqrt(wind_uv**2 * (beamwidth)**2 / 2.76)
      sig_T = sqrt(3*kolmogorov/2. * (EDR/(2.*pi))**(2./3.) * (L_s**(2./3.) - L_lambda**(2./3.)))
      specBroad = sqrt(sig_B**2 + sig_T**2)

    errorstatus = err
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    return

  end subroutine estimate_spectralBroadening
end module radar_spectral_broadening
