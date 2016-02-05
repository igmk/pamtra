function eps_ice(T,f)

    ! This function calculates the complex permittivity of ice.
    ! No explicit freq. dependence for Re(eps), const. value for T < 240.
    !
    ! Input:
    !	 T  temperature in Kelvin
    !	 f  frequency in GHz
    !
    ! Result:
    !	 eps_ice complex permittivity of natural ice
    !
    ! References:
    !      Maetzler 2006: Thermal microwave radiation: Application for remote sensing
    !      Personal communication with Maetzler

    use kinds, only: dbl

    implicit none

    real(kind=dbl), intent(in) :: T, &
    f

    real(kind=dbl) :: eps_real,  &
    mit,       &
    alpha,     &
    beta,     &
    beta_m,    &
    delta_beta,&
    eps_imag

    complex(kind=dbl) :: eps_ice

    eps_real = 0._dbl

    if (T >= 240._dbl) eps_real = 3.1884_dbl + 9.1d-4*(T - 273._dbl)
    if (T < 240._dbl) eps_real = 3.1884_dbl + 9.1d-4*(240._dbl - 273._dbl)

    ! "modified inverse temperature"
    mit = (300._dbl/T)-1._dbl

    alpha = (0.00504_dbl + 0.0062_dbl*mit)*exp(-22.1_dbl*mit)

    beta_m = (0.0207_dbl/T)*(exp(335._dbl/T))/(exp(335._dbl/T)-1._dbl)**2.
    beta_m = beta_m + 1.16d-11*(f**2.)
    delta_beta = exp(-9.963_dbl + 0.0372_dbl*(T-273.16_dbl))

    beta = beta_m + delta_beta
    eps_imag = (alpha/f)+beta*f

    eps_ice = dcmplx(eps_real, eps_imag)

    return

end function eps_ice
