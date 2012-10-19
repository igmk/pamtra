subroutine mie_spec(f, mindex, dia1, dia2, nbins, maxleg,   &
     ad, bd, alpha, gamma, lphase_flag, extinction, albedo, back_scatt,  &
     nlegen, legen, legen2, legen3, legen4, aerodist,density,wc,&
     diameter, qback_spec)
  ! note that mindex has the convention with negative imaginary part
  !     
  ! computes the mie scattering properties for a gamma or lognormal 
  ! distribution of spheres.          
  ! in addition to mie routine spectra of diamter, ndens and qback are returned
                                

  use kinds
  use nml_params, only: verbose

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
  real(kind=dbl), parameter :: pi = 3.14159265358979d0
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
  real(kind=dbl), intent(out) :: qback_spec(nbins+1)
  real(kind=dbl):: ndens

  character :: aerodist*1 

  real(kind=dbl), parameter :: c = 299792458.d0


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
  if (nbins .gt. 0) del_d = (dia2 - dia1) / nbins
  tot_mass = 0.
  do ir = 1, nbins+1
     diameter(ir) = dia1 + (ir - 1) * del_d
     ndens = distribution(ad, bd, alpha, gamma, diameter(ir), aerodist)  ! number density [1/m⁴]
!removed for TESTING only by MAX
!      if ((ir .eq. 1 .or. ir .eq. nbins+1) .and. nbins .gt. 0) then
!         ndens = 0.5 * ndens
!      end if
     tot_mass = tot_mass + ndens*del_d*pi/6.d0*density*diameter(ir)**3.d0
!removed for TESTING only by MAX
!      if ((ir .eq. nbins+1) .and. (tot_mass/wc*100. .lt. 99.9d0)) then
!       ndens = ndens + (wc-tot_mass)/(del_d*pi/6.d0*density*diameter(ir)**3.d0)
!       tot_mass = wc
!      end if



     x = pi * diameter(ir) / wavelength ! size parameter
     nmie = 0 
     call miecalc(nmie, x, msphere, a, b) ! calculate a and b
     !calculate the efficencies
     call miecross(nmie, x, a, b, qext, qscat, qback)
     ! sum up extinction, scattering, and backscattering as cross-sections/pi .pi is added in a later step
     qext =   qext  * ndens * (diameter(ir)/2.d0)**2         ! [m²/m⁴]!
     qscat =  qscat * ndens * (diameter(ir)/2.d0)**2        ! [m²/m⁴]!
     qback =  qback * ndens * (diameter(ir)/2.d0)**2        !  [m²/m⁴]!
 
     !integrate=sum up . del_d is added at a later step!
     sumqe = sumqe + qext 
     sumqs = sumqs + qscat
     sumqback = sumqback + qback 

     qback_spec(ir) =  qback * pi  ! volumetric backscattering corss section for radar simulator in [m²/m³]

     if (lphase_flag) then 
        nmie = min(nmie, nterms) 
        do i = 1, nquad 
           call mieangle(nmie, a, b, mu(i), p1, p2, p3, p4)
           sump1(i) = sump1(i) + p1 * ndens
           sump2(i) = sump2(i) + p2 * ndens
           sump3(i) = sump3(i) + p3 * ndens
           sump4(i) = sump4(i) + p4 * ndens
        end do
     end if
  end do

  !   multiply the sums by the integration delta and other constants
  !   put quadrature weights in angular array for later usage    
  if (nbins .eq. 0) del_d = 1.0d0

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

  if (verbose .gt. 1) print*, 'finished with mie'
  return 
end subroutine mie_spec
