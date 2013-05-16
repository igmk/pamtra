module double_moments_module

    use kinds

    implicit none

    !Parameter of the drop size distribution when 2moments scheme is used
    real (kind=dbl), dimension(5) :: gamma_cloud,&
    gamma_rain, &
    gamma_ice, &
    gamma_snow, &
    gamma_graupel, &
    gamma_hail

contains

    subroutine double_moments( &
    rohydro,ntot,nu_mass_const,mu_mass,alpha_dm,beta_dm, &     !IN
    a_dia,lambda_dia,nu_dia,mu_dia,alpha_md,beta_md,&      !OUT
    cwc)

        !This routine convert the parameters of a drop size distribution as a function of MASS F(m)
        !     F(m)=A_mass * m ^ Nu_mass * exp(- Lambda_mass * m ^ Mu_mass)
        !to the parameters of a drop size distribution as a function of diameter f(d)
        !     F(r)=A_dia * r ^ Nu_dia * exp(- Lambda_dia * d ^ Mu_dia)
        !
        !As a reference for conversion from MASS to DIAMETER see table.2 and eq.51 in
        !GW Petty and W Huang JAS, 2011
        !
        !As a reference for the calculation of a and lambda for F(m) see appendix.a Seifert and Beheng, 2006, MAP
        !
        !INPUT:
        !rohydro      : density of hydrometeor                        [kg/m^3]
        !ntot         : total number concentration                    [#/m^3]
        !nu_mass_const: Nu parameter of F(m)
        !mu_mass      : Mu parameter of F(m)
        !alpha_dm     : alpha parameter for diameter-mass relation    [m*kb^(-beta_dm)]
        !beta_dm      : beta parameter for diameter-mass relation
        !cwc          : in case of rain outside the cloud nu_mass is changed
        !
        !OUTPUT
        !a_dia        : A parameter of F(d)                           [#/m^(nu_dia+4)]
        !lambda_dia   : Lambda parameter of F(d)                      [m^(-mu_dia)]
        !nu_dia       : Nu parameter of F(d)
        !mu_dia       : Mu parameter of F(d)
        !alpha_md     : alpha parameter for mass-diameter relation    [m*kb^(-beta_md)]
        !beta_md      : beta parameter for mass-diameter relation

        use constants, only: pi
        use report_module

        implicit none

        real (kind=dbl), intent(in)           ::      rohydro,ntot,nu_mass_const,mu_mass,alpha_dm,beta_dm
        real (kind=dbl), intent(out)          ::      a_dia,lambda_dia,nu_dia,mu_dia, alpha_md, beta_md
        real (kind=dbl), optional, intent(in) ::      cwc

        real (kind=dbl)           ::      a_temp, lambda_temp
        real (kind=dbl)           ::      gammln
        real (kind=dbl)           ::      arg1, arg2, gamma1, gamma2, rho_w
        real (kind=dbl)           ::      rain_cmu0, rain_cmu1, rain_cmu2, rain_cmu3, rain_cmu4
        real (kind=dbl)           ::      x_r, D_m, mue, nu_mass

        if (verbose .gt. 1) print*, 'Entering double_moments'

        !..mu-Dm-relation for raindrops based on 1d-bin model
        rain_cmu0 = 6.0             ! Axel's 2008 relation
        rain_cmu1 = 30.0            !
        rain_cmu2 = 1.00d+3         !
        rain_cmu3 = 1.10d-3         ! D_eq,break
        rain_cmu4 = 1.0             !
        rho_w = 1000.               !density of water

        ! IF BELOW CLOUD then use the mu-D_m relation
        IF (present(cwc) .and. cwc .lt. 3.d-7) THEN
            x_r = rohydro / ntot                           !MEAN mass in SI
            D_m = ( 6. / (rho_w*pi) * x_r )**(1./3.)      !corresponding mean diameter assuming a spherical shape
            IF (D_m.LE.rain_cmu3) THEN
                mue = rain_cmu0*TANH((4.*rain_cmu2*(D_m-rain_cmu3))**2) &
                + rain_cmu4
                nu_mass = (mue+1)*beta_dm - 1
            ELSE
                mue = rain_cmu1*TANH((1.*rain_cmu2*(D_m-rain_cmu3))**2) &
                + rain_cmu4
                nu_mass = (mue+1)*beta_dm - 1
            ENDIF
        ELSE
            nu_mass = nu_mass_const
        end if


        arg1=(nu_mass+1.d0)/mu_mass
        arg2=(nu_mass+2.d0)/mu_mass
        gamma1=exp(gammln(arg1))
        gamma2=exp(gammln(arg2))
        !     a and lambda of the drop size distribution as a function of MASS
        lambda_temp=(rohydro*gamma1/ntot/gamma2)**(-mu_mass)
        a_temp=ntot*mu_mass*lambda_temp**(arg1)/gamma1

        !From diameter-mass relation      to      mass-diameter
        !     D=alpha_dm * m^beta_dm      -->     m=alpha_md * D^beta_md
        beta_md=1.d0/beta_dm
        alpha_md=(1.d0/alpha_dm)**beta_md

        !a, lambda, nu, mu from the drop size distribution as a function of mass F(m) to DIAMETER F(D)
        a_dia=a_temp*alpha_md**(nu_mass+1.)*beta_md
        lambda_dia=lambda_temp*alpha_md**(mu_mass)
        mu_dia=mu_mass*beta_md
        nu_dia=beta_md*(nu_mass+1.d0)-1.d0

        if (verbose .gt. 1) print*, 'Exiting double_moments'

        return
    end subroutine double_moments

    subroutine double_moments_module_read(moments_file)
        character(20) :: dummy
        character(20) :: moments_file

        open(118,file=moments_file)
        !read NU & Mu parameter of the drop size distr. as a function of MASS
        !, Alpha & Beta parameter of Diameter-Mass function, and x_max as maximum mass of particle
        !gamma_xxx(1)=nu    gamma_xxx(2)=mu     gamma_xxx(3)=alpha      gamma_xxx(4)=beta gamma_xxx(5)=x_max
        read(118,*) dummy,gamma_cloud
        read(118,*) dummy,gamma_rain
        read(118,*) dummy,gamma_ice
        read(118,*) dummy,gamma_snow
        read(118,*) dummy,gamma_graupel
        read(118,*) dummy,gamma_hail
        close(118)

        return
    end subroutine double_moments_module_read

end module double_moments_module
