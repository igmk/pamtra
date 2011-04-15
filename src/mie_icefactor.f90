subroutine mie_icefactor(f, t, m_ice,    &
  a_mtox, bcoeff, rad1, rad2, numrad, maxleg, ad, bd, alpha, &
  gamma, lphase_flag, extinction, albedo, back_scatt, nlegen, legen,  &
  legen2, legen3, legen4, aerodist,density,ice_type)                                 
!    computing the scattering properties according to                  
!    ice sphere model, i.e. the electromagnetic properties of the      
!     particle are computed by assuming that they are the same          
!     as the equivalent mass sphere
!                                     
! note that mindex has the convention with negative imaginary part      

use kinds

implicit none

real(kind=dbl), intent(in) :: f,  &! frequency [GHz]
 				t    ! temperature [K]

integer :: maxleg, nlegen, numrad, ice_type 
logical :: lphase_flag 
real(kind=dbl) :: wavelength, rad1, rad2 
real(kind=dbl), intent(in) :: ad, bd, alpha, gamma 
complex(kind=dbl) :: m_ice
real(kind=dbl) :: a_mtox, bcoeff 
real(kind=dbl) :: extinction, albedo, back_scatt, legen (200), legen2 (200),&
  legen3 (200), legen4 (200)                                        
integer, parameter :: maxn  = 5000
real(kind=dbl), parameter :: pi = 3.14159265358979d0
integer :: nterms, nquad, nmie, nleg 
integer :: i, l, m, ir
real(kind=dbl) :: x, delrad, radius, ndens, tmp, radius_ice, density
real(kind=dbl) :: qext, qscat, qback, scatter 
real(kind=dbl) :: distribution 
real(kind=dbl) :: mu (maxn), wts (maxn) 
real(kind=dbl) :: p1, pl, pl1, pl2, p2, p3, p4 
real(kind=dbl) :: sumqe, sumqs, sumqback
real(kind=dbl), dimension(maxn) :: sump1, coef1, sump2, coef2,   &
sump3, coef3, sump4, coef4            
complex(kind=dbl), dimension(maxn) :: a, b
complex(kind=dbl) :: msphere 
character :: aerodist * 1

real(kind=dbl), parameter :: c = 299792458.d0
							      
wavelength = c/(f*1.e9) !
																  
!           find the maximum number of terms required in the mie series,
!       call density_ice(a_mtox, bcoeff, rad2, dens_graup) 
!       rad2_ice = (dens_graup / 917.) **0.33333333 * rad2 
radius_ice = 0.5*(6.*a_mtox*(2.*rad2)**bcoeff/(pi*density*1.e3))**(1./3.)
msphere = m_ice !conjg(m_ice)

call epsi(msphere,f,t,ice_type,density) 

x = 2.0d0 * pi * radius_ice / wavelength 
nterms = 0 
call miecalc (nterms, x, msphere, a, b) 
nlegen = 2 * nterms 
nlegen = min(maxleg, nlegen) 
nquad = (nlegen + 2 * nterms + 2) / 2 
if (nquad.gt.maxn) stop 'mie: maxn exceeded' 
								  
!           get the gauss-legendre quadrature abscissas and weights     
call gausquad(nquad, mu, wts) 
								  
sumqe = 0.0d0 
sumqs = 0.0d0 
sumqback = 0.0d0 
do i = 1, nquad 
  sump1 (i) = 0.0d0 
  sump2 (i) = 0.0d0 
  sump3 (i) = 0.0d0 
  sump4 (i) = 0.0d0 
end do 
								  
!               integration loop over radius of spheres                 
if (numrad .gt. 0) delrad = (rad2 - rad1) / numrad 
do ir = 1, numrad+1 
  radius = rad1 + (ir - 1) * delrad 
  ndens = distribution (ad, bd, alpha, gamma, radius, aerodist)                                         
  if ( (ir .eq. 1 .or. ir .eq. numrad+1) .and. numrad .gt. 0) then 
    ndens = 0.5 * ndens 
  end if 
								    
  nmie = 0 
								    
  !       call density_ice(a_mtox, bcoeff, radius, dens_graup) 
  ! !         write(18,*)'dens',dens_graup                                  
  !       radius_ice = (dens_graup / 917.) **0.33333333 * radius 

  radius_ice = 0.5*(6.*a_mtox*(2.*radius)**bcoeff/(pi*density*1.e3))**(1./3.)
  x = 2.0d0 * pi * radius_ice / wavelength 
  call epsi(msphere,f,t,ice_type,density) 
								    
  call miecalc (nmie, x, msphere, a, b) 
  call miecross (nmie, x, a, b, qext, qscat, qback) 
  sumqe = sumqe+qext * ndens * radius_ice**2 
  sumqs = sumqs + qscat * ndens * radius_ice**2 
  sumqback = sumqback + qback * ndens * radius_ice**2 
  !          write(*,*)'mie in',qext,qscat,ndens,radius,x,nmie            
  if (lphase_flag) then 
    nmie = min0(nmie, nterms) 
    do i = 1, nquad 
      call mieangle (nmie, a, b, mu (i), p1, p2, p3, p4) 
      sump1 (i) = sump1 (i) + p1 * ndens 
      sump2 (i) = sump2 (i) + p2 * ndens 
      sump3 (i) = sump3 (i) + p3 * ndens 
      sump4 (i) = sump4 (i) + p4 * ndens 
    end do 
  end if 
end do 
								  
								  
!           multiply the sums by the integration delta and other constan
!             put quadrature weights in angular array for later         
if (numrad.eq.0) delrad = 1.0d0 
								  
extinction = pi * sumqe * delrad 
scatter = pi * sumqs * delrad 
back_scatt = pi * sumqback * delrad 
albedo = scatter / extinction 
								  
!         if the phase function is not desired then leave now           
if ( .not. lphase_flag) return 
								  
tmp = (wavelength**2 / (pi * scatter) ) * delrad 
do i = 1, nquad 
  sump1 (i) = tmp * sump1 (i) * wts (i) 
  sump2 (i) = tmp * sump2 (i) * wts (i) 
  sump3 (i) = tmp * sump3 (i) * wts (i) 
  sump4 (i) = tmp * sump4 (i) * wts (i) 
end do 
								  
!           integrate the angular scattering functions times legendre   
!             polynomials to find the legendre coefficients             
do m = 1, nlegen + 1 
  coef1 (m) = 0.0d0 
  coef2 (m) = 0.0d0 
  coef3 (m) = 0.0d0 
  coef4 (m) = 0.0d0 
end do 
!           use upward recurrence to find legendre polynomials          
do i = 1, nquad 
  pl1 = 1.0d0 
  pl = 1.0d0 
  do l = 0, nlegen 
    m = l + 1 
    if (l .gt. 0) pl = (2*l-1)*mu(i)*pl1/l-(l-1)*pl2/l                                                                 
    coef1 (m) = coef1 (m) + sump1 (i) * pl 
    coef2 (m) = coef2 (m) + sump2 (i) * pl 
    coef3 (m) = coef3 (m) + sump3 (i) * pl 
    coef4 (m) = coef4 (m) + sump4 (i) * pl 
    pl2 = pl1 
    pl1 = pl 
  end do 
end do 
nleg = nlegen 
do l = 0, nleg 
  m = l + 1 
  legen (m) = (2 * l + 1) / 2.0 * coef1 (m) 
  legen2 (m) = (2 * l + 1) / 2.0 * coef2 (m) 
  legen3 (m) = (2 * l + 1) / 2.0 * coef3 (m) 
  legen4 (m) = (2 * l + 1) / 2.0 * coef4 (m) 
  if (legen (m) .gt. 1.0e-7) nlegen = l 
end do 
								  
return 

end subroutine mie_icefactor                      
