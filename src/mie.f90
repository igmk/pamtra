subroutine mie(f, mindex, dia1, dia2, nbins, maxleg,   &
ad, bd, alpha, gamma, lphase_flag, extinction, albedo, back_scatt,  &
nlegen, legen, legen2, legen3, legen4, aerodist,density,wc,&
diameter, back_spec)

    ! subroutine mie_densitydep_spheremasseq(f, t, m_ice,    &
    !      a_mtox, bcoeff, dia1, dia2, nbins, maxleg, ad, bd, alpha, &
    !      gamma, lphase_flag, extinction, albedo, back_scatt, nlegen, legen,  &
    !      legen2, legen3, legen4, aerodist,density,wc)

    !out
    !diameter: diameter spectrum [m]
    !back_spec: backscattering cross section per volume per del_d [m²/m⁴]

    ! note that mindex has the convention with negative imaginary part
    !
    ! computes the mie scattering properties for a gamma or lognormal
    ! distribution of spheres.
    ! in addition to mie routine spectra of diamter and qback are returned
                                

    use kinds
    !  use settings, only: verbose
    use constants, only: pi,c
    use report_module

    implicit none

    integer :: maxleg, nlegen, nbins
    logical :: lphase_flag

    real(kind=dbl), intent(in) :: f  ! frequency [GHz]

    real(kind=dbl) :: wavelength, dia1, dia2
    real(kind=dbl) :: ad, bd, alpha, gamma
    complex(kind=dbl) :: mindex
    real(kind=dbl) :: extinction, albedo, back_scatt
    real(kind=dbl), dimension(200) :: legen, legen2, legen3, legen4
    integer, parameter :: maxn = 5000
    integer :: nterms, nquad, nmie, nleg
    integer :: i, l, m, ir
    real(kind=dbl) :: x, del_d, tmp, tot_mass, wc, density
    real(kind=dbl) :: qext, qscat, qback, scatter
    real(kind=dbl) :: distribution
    real(kind=dbl) :: mu(maxn), wts(maxn)
    real(kind=dbl) :: p1, pl, pl1, pl2, p2, p3, p4
    real(kind=dbl) :: sumqe, sumqs, sumqback
    real(kind=dbl) :: sump1(maxn), coef1(maxn), sump2(maxn), coef2(maxn),   &
    sump3(maxn), coef3(maxn), sump4(maxn), coef4(maxn)
    complex(kind=dbl) :: a(maxn), b(maxn), msphere

    real(kind=dbl), intent(out) :: diameter(nbins+1)
    real(kind=dbl), intent(out) :: back_spec(nbins+1)
    real(kind=dbl):: ndens(nbins+1)
    real(kind=dbl) :: K2, dielec_water !for debuging only

    character :: aerodist*1



    if (verbose .gt. 1) print*, 'Entering mie'

    legen4 = 0._dbl
    wavelength = c/(f*1.e9) !

    ! find the maximum number of terms required in the mie series,
    msphere = mindex
    x = pi * dia2 / wavelength
    nterms = 0
    call miecalc(nterms, x, msphere, a, b) ! miecalc returns nterms
    nlegen = 2 * nterms
    nlegen = min(maxleg, nlegen)
    nquad = (nlegen + 2 * nterms + 2) / 2
    if (nquad .gt. maxn) stop 'mie: maxn exceeded'

    !  get the gauss-legendre quadrature abscissas and weights
    call gausquad(nquad, mu, wts)

    sumqe = 0.0d0
    sumqs = 0.0d0
    sumqback = 0.0d0
    do i = 1, nquad
        sump1(i) = 0.0d0
        sump2(i) = 0.0d0
        sump3(i) = 0.0d0
        sump4(i) = 0.0d0
    end do

    !   integration loop over diameter of spheres
    if (nbins .gt. 0) then
        del_d = (dia2 - dia1) / nbins
    else
        del_d = 1.d0
    end if
  
    tot_mass = 0.
    do ir = 1, nbins+1
        diameter(ir) = dia1 + (ir - 1) * del_d
        ndens(ir) = distribution(ad, bd, alpha, gamma, diameter(ir), aerodist)  ! number density [1/m⁴]
        !first and last bin get only have of the particles (Why actually?)
        if ((ir .eq. 1 .or. ir .eq. nbins+1) .and. nbins .gt. 0) then
            ndens(ir) = 0.5 * ndens(ir)
        end if
        tot_mass = tot_mass + ndens(ir)*del_d*pi/6.d0*density*diameter(ir)**3.d0
        !if mass is missing put into last bin
        if ((ir .eq. nbins+1) .and. (tot_mass/wc*100. .lt. 99.9d0)) then
            ndens(ir) = ndens(ir) + (wc-tot_mass)/(del_d*pi/6.d0*density*diameter(ir)**3.d0)
            tot_mass = wc
        end if

        x = pi * diameter(ir) / wavelength ! size parameter
        nmie = 0
        call miecalc(nmie, x, msphere, a, b) ! calculate a and b
        !calculate the efficencies
        call miecross(nmie, x, a, b, qext, qscat, qback)
        ! sum up extinction, scattering, and backscattering as cross-sections/pi .pi is added in a later step
        qext =   qext  * ndens(ir) * (diameter(ir)/2.d0)**2         ! [m²/m⁴]!
        qscat =  qscat * ndens(ir) * (diameter(ir)/2.d0)**2        ! [m²/m⁴]!
        qback =  qback * ndens(ir) * (diameter(ir)/2.d0)**2        !  [m²/m⁴]! cross section per volume per del_d
 
        !integrate=sum up . del_d is added at a later step!
        sumqe = sumqe + qext
        sumqs = sumqs + qscat
        sumqback = sumqback + qback

        back_spec(ir) =  qback * pi  ! volumetric backscattering corss section for radar simulator per del_d in [m²/m⁴]

        if (lphase_flag) then
            nmie = min(nmie, nterms)
            do i = 1, nquad
                call mieangle(nmie, a, b, mu(i), p1, p2, p3, p4)
                sump1(i) = sump1(i) + p1 * ndens(ir)
                sump2(i) = sump2(i) + p2 * ndens(ir)
                sump3(i) = sump3(i) + p3 * ndens(ir)
                sump4(i) = sump4(i) + p4 * ndens(ir)
            end do
        end if
    end do

    !   multiply the sums by the integration delta and other constants
    !   put quadrature weights in angular array for later usage

    extinction = pi * sumqe * del_d      ! extinction  [m*m²/m⁴]!
    scatter = pi * sumqs * del_d         ! scattering [m*m²/m⁴]!
    back_scatt = pi * sumqback * del_d   ! back scattering [m*m²/m⁴]!
    albedo = scatter / extinction        ! single scattering albedo


    ! if the phase function is not desired then leave now
    if ( .not. lphase_flag) return

    tmp = (wavelength**2 / (pi * scatter)) * del_d
    do i = 1, nquad
        sump1(i) = tmp * sump1(i) * wts(i)
        sump2(i) = tmp * sump2(i) * wts(i)
        sump3(i) = tmp * sump3(i) * wts(i)
        sump4(i) = tmp * sump4(i) * wts(i)
    end do

    !  integrate the angular scattering functions times legendre
    !  polynomials to find the legendre coefficients
    do m = 1, nlegen + 1
        coef1(m) = 0.0d0
        coef2(m) = 0.0d0
        coef3(m) = 0.0d0
        coef4(m) = 0.0d0
    end do

    !  use upward recurrence to find legendre polynomials
    do i = 1, nquad
        pl1 = 1.0d0
        pl = 1.0d0
        do l = 0, nlegen
            m = l + 1
            if (l .gt. 0) pl = (2*l-1)*mu(i)*pl1/l-(l-1)*pl2/l
            coef1(m) = coef1(m) + sump1(i) * pl
            coef2(m) = coef2(m) + sump2(i) * pl
            coef3(m) = coef3(m) + sump3(i) * pl
            coef4(m) = coef4(m) + sump4(i) * pl
            pl2 = pl1
            pl1 = pl
        end do
    end do

    nleg = nlegen
    do l = 0, nleg
        m = l + 1
        legen(m) = (2 * l + 1) / 2.0 * coef1(m)
        legen2(m) = (2 * l + 1) / 2.0 * coef2(m)
        legen3(m) = (2 * l + 1) / 2.0 * coef3(m)
        legen4(m) = (2 * l + 1) / 2.0 * coef4(m)
        if (legen(m) .gt. 1.0e-7) nlegen = l
    end do


    !for debuging puposes, calculate the Rayleigh Backscattering! makes only sence for SMALL spheres
    if (verbose .gt. 2) then
        K2=dielec_water(0.D0,10.d0,f)
        print*,"RAYLEIGH BACKSCATTERING with K2 at 10deg"
        print*,"Z=",10.d0*log10(K2*SUM((1d3*diameter)**6*(ndens*1d-3)*&
        1d3*(diameter(2)-diameter(1)))),"dBz"
    end if

    if (verbose .gt. 1) print*, 'finished with mie'
    return
end subroutine mie
