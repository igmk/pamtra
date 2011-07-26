subroutine dda_db_liu(f, t, npt, mindex, rad1, rad2, numrad, maxleg,   &
  ad, bd, alpha, gamma, lphase_flag, extinction, albedo, back_scatt,  &
  nlegen, legen, legen2, legen3, legen4, aerodist)
                  
! note that mindex has the convention with negative imaginary part      
! computes the mie scattering properties for a gamma or lognormal 
! distribution of spheres.
                                       
  use kinds

  implicit none

  integer :: maxleg, nlegen, numrad

  integer :: iret, is_loaded

  logical :: lphase_flag

  integer, intent(in) :: npt    ! id of particle in db

  real(kind=dbl), intent(in) :: f,  &! frequency [GHz]
 				t    ! temperature [K]

  real(kind=dbl) :: wavelength, rad1, rad2, r_ice_eq
  real(kind=dbl) :: ad, bd, alpha, gamma 
  complex(kind=dbl) :: mindex 
  real(kind=dbl) :: extinction, albedo, back_scatt, &
	    abs_liu,sca_liu,bsc_liu,g,&
	    legen(200), legen2(200), legen3(200), legen4(200)
  real(kind=dbl), dimension(37) :: p_liu
  integer, parameter :: maxn = 5000
  real(kind=dbl), parameter :: pi = 3.14159265358979d0
  integer :: nterms, nquad, nmie, nleg 
  integer :: i, l, m, ir
  real(kind=dbl) :: x, delrad, radius, ndens, tmp, tot_mass
  real(kind=dbl) :: qext, qscat, qback, scatter 
  real(kind=dbl) :: distribution 
  real(kind=dbl) :: mu(maxn), wts(maxn) 
  real(kind=dbl) :: p1, pl, pl1, pl2, p2, p3, p4 
  real(kind=dbl) :: sumqe, sumqs, sumqback 
  real(kind=dbl) :: sump1(maxn), coef1(maxn), sump2(maxn), coef2(maxn),   &
    sump3(maxn), coef3(maxn), sump4(maxn), coef4(maxn)            
  complex(kind=dbl) :: a(maxn), b(maxn), msphere 
  character :: aerodist*1 

  real(kind=dbl), parameter :: c = 299792458.d0

  is_loaded = 0

  qext = 0.
  qscat = 0.
  qback = 0.

  wavelength = c/(f*1.e9) !
								    
  ! find the maximum number of terms required in the mie series,
  msphere = mindex
  x = 2.0d0 * pi * rad2 / wavelength
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
								    
  !   integration loop over radius of spheres 
  if (numrad .gt. 0) delrad = (rad2 - rad1) / numrad
  tot_mass = 0.
  do ir = 1, numrad+1 
    radius = rad1 + (ir - 1) * delrad
    ndens = distribution(ad, bd, alpha, gamma, radius, aerodist)  ! number density
    if ((ir .eq. 1 .or. ir .eq. numrad+1) .and. numrad .gt. 0) then 
	ndens = 0.5 * ndens 
    end if 
!    tot_mass = tot_mass +4./3.*pi*(radius)**3*ndens*2.*delrad*0.917d3
!    tot_mass = tot_mass +0.038*(2.*radius)**2*ndens*2.*delrad
!    print*, radius, ndens, tot_mass
    x = 2.0d0 * pi * radius / wavelength ! size parameter
    nmie = 0 
!    call miecalc(nmie, x, msphere, a, b) ! calculate a and b
!    call miecross(nmie, x, a, b, qext, qscat, qback)
    print*, f,t,npt,2*radius*1.e6
    call scatdb(f,t,npt,2*radius*1.e6,abs_liu,sca_liu,bsc_liu,g,p_liu,r_ice_eq,iret,is_loaded)
    print*, iret,abs_liu,sca_liu,bsc_liu,g,p_liu,r_ice_eq

    ! sum up extinction, scattering, and backscattering as cross-sections/pi
    sumqe = sumqe + qext * ndens * radius**2 
    sumqs = sumqs + qscat * ndens * radius**2 
    sumqback = sumqback + qback * ndens * radius**2 
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
  if (numrad .eq. 0) delrad = 1.0d0 

  extinction = pi * sumqe * delrad      ! [1/m]
  scatter = pi * sumqs * delrad         ! [1/m]
  back_scatt = pi * sumqback * delrad   ! [1/m]
  albedo = scatter / extinction 
								    
  ! if the phase function is not desired then leave now           
  if ( .not. lphase_flag) return 
								    
  tmp = (wavelength**2 / (pi * scatter)) * delrad 
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
								    
  return
end subroutine dda_db_liu                            

!   use kinds
! 
!   implicit none
! 
!   integer :: maxleg, nlegen, numrad 
!   logical :: lphase_flag 
! 
!   integer, intent(in) :: npt    ! id of particle in db
! 
!   real(kind=dbl), intent(in) :: f,  &! frequency [GHz]
! 				t    ! temperature [K]
! 
!   real(kind=dbl) :: wavelength, rad1, rad2 
!   real(kind=dbl) :: ad, bd, alpha, gamma 
!   complex(kind=dbl) :: mindex 
!   real(kind=dbl) :: extinction, albedo, back_scatt, &
! 	    legen(200), legen2(200), legen3(200), legen4(200)                                        
!   integer, parameter :: maxn = 5000
!   real(kind=dbl), parameter :: pi = 3.14159265358979d0
!   integer :: nterms, nquad, nmie, nleg 
!   integer :: i, l, m, ir
!   real(kind=dbl) :: x, delrad, radius, ndens, tmp 
!   real(kind=dbl) :: qext, qscat, qback, scatter 
!   real(kind=dbl) :: distribution 
!   real(kind=dbl) :: mu(maxn), wts(maxn) 
!   real(kind=dbl) :: p1, pl, pl1, pl2, p2, p3, p4 
!   real(kind=dbl) :: sumqe, sumqs, sumqback, one 
!   real(kind=dbl) :: sump1(maxn), coef1(maxn), sump2(maxn), coef2(maxn),   &
!     sump3(maxn), coef3(maxn), sump4(maxn), coef4(maxn)            
!   complex(kind=dbl) :: a(maxn), b(maxn), msphere 
!   character :: aerodist*1 
! 
!   real(kind=dbl), parameter :: c = 299792458.d0
! 								    
!   wavelength = c/(f*1.e9) !
! 								    
!   !           find the maximum number of terms required in the mie series,
!   msphere = mindex 
!   x = 2.0d0 * pi * rad2 / wavelength 
!   nterms = 0 
!   call miecalc(nterms, x, msphere, a, b) ! miecalc returns nterms
!   nlegen = 2 * nterms 
!   nlegen = min(maxleg, nlegen) 
!   nquad = (nlegen + 2 * nterms + 2) / 2 
!   if (nquad .gt. maxn) stop 'mie: maxn exceeded' 
! 								    
!   !  get the gauss-legendre quadrature abscissas and weights     
!   call gausquad(nquad, mu, wts) 
! 								    
!   sumqe = 0.0d0 
!   sumqs = 0.0d0 
!   sumqback = 0.0d0 
!   do i = 1, nquad 
!     sump1(i) = 0.0d0 
!     sump2(i) = 0.0d0 
!     sump3(i) = 0.0d0 
!     sump4(i) = 0.0d0 
!   end do 
! 								    
!   !   integration loop over radius of spheres                 
!   if (numrad .gt. 0) delrad = (rad2 - rad1) / numrad 
!   do ir = 1, numrad+1 
!     radius = rad1 + (ir - 1) * delrad 
!     ndens = distribution(ad, bd, alpha, gamma, radius, aerodist)  ! number density                                     
!     if ((ir .eq. 1 .or. ir .eq. numrad+1) .and. numrad .gt. 0) then 
! 	ndens = 0.5 * ndens 
!     end if 
!     x = 2.0d0 * pi * radius / wavelength ! size parameter
!     nmie = 0 
!     call miecalc(nmie, x, msphere, a, b) ! calculate a and b
!     call miecross(nmie, x, a, b, qext, qscat, qback)
!     ! sum up extinction, scattering, and backscattering as cross-sections/pi
!     sumqe = sumqe + qext * ndens * radius**2 
!     sumqs = sumqs + qscat * ndens * radius**2 
!     sumqback = sumqback + qback * ndens * radius**2 
!     if (lphase_flag) then 
!       nmie = min(nmie, nterms) 
!       do i = 1, nquad 
! 	call mieangle(nmie, a, b, mu(i), p1, p2, p3, p4) 
! 	sump1(i) = sump1(i) + p1 * ndens 
! 	sump2(i) = sump2(i) + p2 * ndens 
! 	sump3(i) = sump3(i) + p3 * ndens 
! 	sump4(i) = sump4(i) + p4 * ndens 
!       end do 
!     end if 
!   end do 
! 								    
! 								    
!   !   multiply the sums by the integration delta and other constants
!   !   put quadrature weights in angular array for later usage    
!   if (numrad .eq. 0) delrad = 1.0d0 
! 
!   extinction = pi * sumqe * delrad      ! [1/m]
!   scatter = pi * sumqs * delrad         ! [1/m]
!   back_scatt = pi * sumqback * delrad   ! [1/m]
!   albedo = scatter / extinction 
! 								    
!   ! if the phase function is not desired then leave now           
!   if ( .not. lphase_flag) return 
! 								    
!   tmp = (wavelength**2 / (pi * scatter)) * delrad 
!   do i = 1, nquad 
!     sump1(i) = tmp * sump1(i) * wts(i) 
!     sump2(i) = tmp * sump2(i) * wts(i) 
!     sump3(i) = tmp * sump3(i) * wts(i) 
!     sump4(i) = tmp * sump4(i) * wts(i) 
!   end do 
! 								    
!   !  integrate the angular scattering functions times legendre   
!   !  polynomials to find the legendre coefficients             
!   do m = 1, nlegen + 1 
!     coef1(m) = 0.0d0 
!     coef2(m) = 0.0d0 
!     coef3(m) = 0.0d0 
!     coef4(m) = 0.0d0 
!   end do 
!   !  use upward recurrence to find legendre polynomials          
!   do i = 1, nquad 
!     pl1 = 1.0d0 
!     pl = 1.0d0 
!     do l = 0, nlegen 
!       m = l + 1 
!       if (l .gt. 0) pl = (2*l-1)*mu(i)*pl1/l-(l-1)*pl2/l                                                              
!       coef1(m) = coef1(m) + sump1(i) * pl 
!       coef2(m) = coef2(m) + sump2(i) * pl 
!       coef3(m) = coef3(m) + sump3(i) * pl 
!       coef4(m) = coef4(m) + sump4(i) * pl 
!       pl2 = pl1 
!       pl1 = pl 
!     end do 
!   end do 
!   nleg = nlegen 
!   do l = 0, nleg 
!     m = l + 1 
!     legen(m) = (2 * l + 1) / 2.0 * coef1(m) 
!     legen2(m) = (2 * l + 1) / 2.0 * coef2(m) 
!     legen3(m) = (2 * l + 1) / 2.0 * coef3(m) 
!     legen4(m) = (2 * l + 1) / 2.0 * coef4(m) 
!     if (legen(m) .gt. 1.0e-7) nlegen = l 
!   end do 
! 								    
!   return 

! 
! subroutine dda_db_liu(f, t, nshp, m_ice,    &
!   m_air, a_mtox, bcoeff, rad1, rad2, numrad, maxleg, ad, bd, alpha, &
!   gamma, lphase_flag, extinction, albedo, back_scatt, nlegen, legen,  &
!   legen2, legen3, legen4, aerodist)                                 
! !c    computing the scattering properties according to                  
! !c    ice sphere model, i.e. the electromegnetic properties of the      
! !     particle are computed by assuming that they are the same          
! !     as the equivalent mass sphere                                     
! ! note that mindex has the convention with negative imaginary part      
! 								  
! 								  
! !       computes the mie scattering properties for a gamma or lognormal 
! !     distribution of spheres.   
!   use kinds
! 
!   implicit none
! 
!   integer :: maxleg, nlegen, numrad 
!   logical :: lphase_flag 
!   real(kind=dbl) :: wavelength, rad1, rad2 
!   real(kind=dbl) :: ad, bd, alpha, gamma 
!   complex(kind=dbl) :: m_ice, m_air, m_mg 
!   real(kind=dbl) :: a_mtox, bcoeff, dens_graup, fvol_ice 
!   real(kind=dbl) :: extinction, albedo, back_scatt, legen(200), legen2(200),&
!   legen3(200), legen4(200)                                        
!   integer, parameter :: maxn = 5000
!   real(kind=dbl), parameter :: pi = 3.14159265358979d0
!   integer :: nterms, nquad, nmie, nleg 
!   integer :: i, l, m, ir 
!   real(kind=dbl) :: x, delrad, radius, ndens, tmp, radius_ice, rad2_ice 
!   real(kind=dbl) :: qext, qscat, qback, scatter 
!   real(kind=dbl) :: distribution 
!   real(kind=dbl) :: mu(maxn), wts(maxn) 
!   real(kind=dbl) :: p1, pl, pl1, pl2, p2, p3, p4 
!   real(kind=dbl) :: sumqe, sumqs, sumqback, one 
!   real(kind=dbl) :: sump1(maxn), coef1(maxn), sump2(maxn), coef2(maxn),   &
!   sump3(maxn), coef3(maxn), sump4(maxn), coef4(maxn)            
!   complex(kind=dbl) :: a(maxn), b(maxn), msphere 
!   character :: aerodist * 1 
! 
!   integer :: is_loaded, iret, nshp
! 
!   real(kind=dbl) :: f, t
!   real :: cabs, csca, cbsc, g, re
!   real, dimension(37) :: ph
! 
!   real(kind=dbl), parameter :: c = 299792458
! 
!   is_loaded = 0
! 
!   wavelength = c/(f*1.e9) !
! 
!   ! find the maximum number of terms required in the mie series,
!   call density_ice(a_mtox, bcoeff, rad2, dens_graup) 
!   rad2_ice = (dens_graup / 0.917)**0.33333333 * rad2 
!   msphere = conjg(m_ice) 
!   x = 2.0d0 * pi * rad2_ice / wavelength  ! size parameter
!   print*, wavelength, rad2_ice 
!   nterms = 0 
!   call miecalc(nterms, x, msphere, a, b) 
!   nlegen = 2 * nterms 
!   nlegen = min0 (maxleg, nlegen) 
!   nquad = (nlegen + 2 * nterms + 2) / 2 
!   if (nquad .gt. maxn) stop 'mie: maxn exceeded' 
! 								    
!   ! get the gauss-legendre quadrature abscissas and weights     
!   call gausquad(nquad, mu, wts) 
! 								    
!   sumqe = 0.0d0 
!   sumqs = 0.0d0 
!   sumqback = 0.0d0 
!   do i = 1, nquad 
!     sump1 (i) = 0.0d0 
!     sump2 (i) = 0.0d0 
!     sump3 (i) = 0.0d0 
!     sump4 (i) = 0.0d0 
!   end do 
! 								    
!   ! integration loop over radius of spheres                 
!   if (numrad .gt. 0) delrad = (rad2 - rad1) / numrad 
!   do ir = 1, numrad+1 
!     radius = rad1 + (ir - 1) * delrad 
!     ndens = distribution(ad, bd, alpha, gamma, radius, aerodist)                                         
!     if ((ir .eq. 1 .or. ir .eq. numrad+1) .and. numrad .gt. 0) then 
!       ndens = 0.5 * ndens 
!     endif 
! 								      
!     nmie = 0 
! 								      
!     call density_ice(a_mtox, bcoeff, radius, dens_graup) 
!     ! write(18,*)'dens',dens_graup                                  
!     radius_ice = (dens_graup / 0.917) **0.33333333 * radius 
! 								      
!     x = 2.0d0 * pi * radius_ice / wavelength 
! 								      
!     call miecalc(nmie, x, msphere, a, b) 
!     call miecross(nmie, x, a, b, qext, qscat, qback) 
!     print*, f, t, nshp, 2.0*radius*1.0e3
!     call scatdb(real(f),real(t),9,real(2.0*radius*1.0e3),cabs,csca,cbsc,g,ph,re,iret,is_loaded)
!     print*, iret
!   !     if(iret .eq. 0) then
!   !       print*, cabs,csca,cbsc,g,re
!   !       do i = 1,37
!   ! 	print*, ph(i)
!   !       end do
!   !     end if
!   !    call scatdb(f,t,nshp,dmax,cabs,csca,cbsc,g,ph,re,iret,is_loaded)
!     sumqe = sumqe+qext * ndens * radius_ice**2 
!     sumqs = sumqs + qscat * ndens * radius_ice**2 
!     sumqback = sumqback + qback * ndens * radius_ice**2 
!     !          write(*,*)'mie in',qext,qscat,ndens,radius,x,nmie            
!     if (lphase_flag) then 
!       nmie = min(nmie, nterms) 
!       do i = 1, nquad 
! 	call mieangle(nmie, a, b, mu(i), p1, p2, p3, p4) 
! 	sump1(i) = sump1(i) + p1 * ndens 
! 	sump2(i) = sump2(i) + p2 * ndens 
! 	sump3(i) = sump3(i) + p3 * ndens 
! 	sump4(i) = sump4(i) + p4 * ndens 
!       end do 
!     end if 
!   end do
!   !  multiply the sums by the integration delta and other constan
!   !  put quadrature weights in angular array for later         
!   if (numrad .eq. 0) delrad = 1.0d0 
! 								    
!   extinction = pi * sumqe * delrad 
!   scatter = pi * sumqs * delrad 
!   back_scatt = pi * sumqback * delrad 
!   albedo = scatter / extinction 
! 								    
!   !  if the phase function is not desired then leave now           
!   if (.not. lphase_flag) return 
! 								    
!   tmp = (wavelength**2 / (pi * scatter)) * delrad 
!   do i = 1, nquad 
!     sump1(i) = tmp * sump1(i) * wts(i) 
!     sump2(i) = tmp * sump2(i) * wts(i) 
!     sump3(i) = tmp * sump3(i) * wts(i) 
!     sump4(i) = tmp * sump4(i) * wts(i) 
!   end do 
! 								    
!   !  integrate the angular scattering functions times legendre   
!   !  polynomials to find the legendre coefficients             
!   do m = 1, nlegen + 1 
!     coef1(m) = 0.0d0 
!     coef2(m) = 0.0d0 
!     coef3(m) = 0.0d0 
!     coef4(m) = 0.0d0 
!   end do 
!   !  use upward recurrence to find legendre polynomials          
!   do i = 1, nquad 
!     pl1 = 1.0d0 
!     pl = 1.0d0 
!     do l = 0, nlegen 
!       m = l + 1 
!       if (l .gt. 0) pl = (2 * l - 1) * mu (i) * pl1 / l - (l - 1) * pl2 / l                                                                 
!       coef1(m) = coef1(m) + sump1(i) * pl 
!       coef2(m) = coef2(m) + sump2(i) * pl 
!       coef3(m) = coef3(m) + sump3(i) * pl 
!       coef4(m) = coef4(m) + sump4(i) * pl 
!       pl2 = pl1 
!       pl1 = pl 
!     end do 
!   end do 
!   nleg = nlegen 
!   do l = 0, nleg 
!     m = l + 1 
!     legen(m) = (2 * l + 1) / 2.0 * coef1(m) 
!     legen2(m) = (2 * l + 1) / 2.0 * coef2(m) 
!     legen3(m) = (2 * l + 1) / 2.0 * coef3(m) 
!     legen4(m) = (2 * l + 1) / 2.0 * coef4(m) 
!     if (legen (m) .gt. 1.0e-7) nlegen = l 
!   end do 
! 								    
!   return 
! end subroutine dda_db_liu                             

! subroutine dda_db_liu(f, t, nshp, mindex, rad1, rad2, numrad, maxleg,   &
! ad, bd, alpha, gamma, lphase_flag, extinction, albedo, back_scatt,  &
! nlegen, legen, legen2, legen3, legen4, aerodist)                  
! ! note that mindex has the convention with negative imaginary part      
! ! computes the mie scattering properties for a gamma or lognormal 
! ! distribution of spheres.                                          
! 
!   use kinds
! 
!   implicit none
! 
!   integer :: maxleg, nlegen, numrad
!   logical :: lphase_flag 
!   real(kind=dbl) :: f, rad1, rad2, wavelength
!   real(kind=dbl) :: ad, bd, alpha, gamma 
!   complex(kind=dbl) :: mindex 
!   real(kind=dbl) :: extinction, albedo, back_scatt, &
! 	    legen(200), legen2(200), legen3(200), legen4(200)                                        
!   integer, parameter :: maxn = 5000
!   real(kind=dbl), parameter :: pi = 3.14159265358979d0
!   real(kind=dbl), parameter :: c = 299792458
!   integer :: nterms, nquad, nmie, nleg 
!   integer :: i, l, m, ir
!   real(kind=dbl) :: x, delrad, radius, ndens, tmp 
!   real(kind=dbl) :: qext, qscat, qback, scatter 
!   real(kind=dbl) :: distribution 
!   real(kind=dbl) :: mu(maxn), wts(maxn) 
!   real(kind=dbl) :: p1, pl, pl1, pl2, p2, p3, p4 
!   real(kind=dbl) :: sumqe, sumqs, sumqback, one 
!   real(kind=dbl) :: sump1(maxn), coef1(maxn), sump2(maxn), coef2(maxn),   &
!     sump3(maxn), coef3(maxn), sump4(maxn), coef4(maxn)            
!   real(kind=dbl) :: t
!   complex(kind=dbl) :: a(maxn), b(maxn), msphere 
!   
!   character :: aerodist*1
! 
!   integer :: is_loaded, iret, nshp
! 
!   real :: cabs, csca, cbsc, g, re
!   real, dimension(37) :: ph
! 
!   is_loaded = 0 
!   
!   !  find the maximum number of terms required in the mie series,
!   msphere = mindex
!   wavelength = c/(f*1.e6)
!   x = 2.0d0 * pi * rad2 / wavelength
!   nterms = 0 
!   call miecalc(nterms, x, msphere, a, b) ! this returns nterms
!   nlegen = 2 * nterms 
!   nlegen = min(maxleg, nlegen) 
!   nquad = (nlegen + 2 * nterms + 2) / 2 
!   if (nquad .gt. maxn) stop 'mie: maxn exceeded' 
! 								    
!   !  get the gauss-legendre quadrature abscissas and weights     
!   call gausquad(nquad, mu, wts) 
! 								    
!   sumqe = 0.0d0 
!   sumqs = 0.0d0 
!   sumqback = 0.0d0 
!   do i = 1, nquad 
!     sump1(i) = 0.0d0 
!     sump2(i) = 0.0d0 
!     sump3(i) = 0.0d0 
!     sump4(i) = 0.0d0 
!   end do 
! 								    
!   !   integration loop over radius of spheres                 
!   if (numrad .gt. 0) delrad = (rad2 - rad1) / numrad 
!   do ir = 1, numrad+1 
!     radius = rad1 + (ir - 1) * delrad 
!     ndens = distribution(ad, bd, alpha, gamma, radius, aerodist)  ! number density                                     
!     if ((ir .eq. 1 .or. ir .eq. numrad+1) .and. numrad .gt. 0) then 
! 	ndens = 0.5 * ndens 
!     end if 
!     x = 2.0d0 * pi * radius / wavelength ! size parameter
!     nmie = 0 
!     call miecalc(nmie, x, msphere, a, b) ! calculate a and b
!     call miecross(nmie, x, a, b, qext, qscat, qback)
!     print*, f, t , nshp, 2.0*radius*1.0e3
!     call scatdb(real(f),real(t),9,real(2.0*radius*1.0e3),cabs,csca,cbsc,g,ph,re,iret,is_loaded)
!     print*, iret
! !     if(iret .eq. 0) then
! !       print*, cabs,csca,cbsc,g,re
! !       do i = 1,37
! ! 	print*, ph(i)
! !       end do
! !     end if
! !    call scatdb(f,t,nshp,dmax,cabs,csca,cbsc,g,ph,re,iret,is_loaded)
! 
!     ! sum up extinction, scattering, and backscattering as cross-sections/pi
!     sumqe = sumqe + qext * ndens * radius**2 
!     sumqs = sumqs + qscat * ndens * radius**2 
!     sumqback = sumqback + qback * ndens * radius**2 
!     !          write(*,*)'mie in',qext,qscat,ndens,radius,x,nmie            
!     if (lphase_flag) then 
!       nmie = min(nmie, nterms) 
!       do i = 1, nquad 
! 	call mieangle(nmie, a, b, mu(i), p1, p2, p3, p4) 
! 	sump1(i) = sump1(i) + p1 * ndens 
! 	sump2(i) = sump2(i) + p2 * ndens 
! 	sump3(i) = sump3(i) + p3 * ndens 
! 	sump4(i) = sump4(i) + p4 * ndens 
!       end do 
!     end if 
!   end do 
! 								    
! 								    
!   !   multiply the sums by the integration delta and other constants
!   !   put quadrature weights in angular array for later usage    
!   if (numrad .eq. 0) delrad = 1.0d0 
! 								    
!   extinction = pi * sumqe * delrad 
!   scatter = pi * sumqs * delrad 
!   back_scatt = pi * sumqback * delrad 
!   albedo = scatter / extinction 
! 								    
!   ! if the phase function is not desired then leave now           
!   if ( .not. lphase_flag) return 
! 								    
!   tmp = (wavelength**2 / (pi * scatter)) * delrad 
!   do i = 1, nquad 
!     sump1(i) = tmp * sump1(i) * wts(i) 
!     sump2(i) = tmp * sump2(i) * wts(i) 
!     sump3(i) = tmp * sump3(i) * wts(i) 
!     sump4(i) = tmp * sump4(i) * wts(i) 
!   end do 
! 								    
!   !  integrate the angular scattering functions times legendre   
!   !  polynomials to find the legendre coefficients             
!   do m = 1, nlegen + 1 
!     coef1(m) = 0.0d0 
!     coef2(m) = 0.0d0 
!     coef3(m) = 0.0d0 
!     coef4(m) = 0.0d0 
!   end do 
!   !  use upward recurrence to find legendre polynomials          
!   do i = 1, nquad 
!     pl1 = 1.0d0 
!     pl = 1.0d0 
!     do l = 0, nlegen 
!       m = l + 1 
!       if (l .gt. 0) pl = (2 * l - 1) * mu (i) * pl1 / l - (l - 1) * pl2 / l                                                              
!       coef1(m) = coef1(m) + sump1(i) * pl 
!       coef2(m) = coef2(m) + sump2(i) * pl 
!       coef3(m) = coef3(m) + sump3(i) * pl 
!       coef4(m) = coef4(m) + sump4(i) * pl 
!       pl2 = pl1 
!       pl1 = pl 
!     end do 
!   end do 
!   nleg = nlegen 
!   do l = 0, nleg 
!     m = l + 1 
!     legen(m) = (2 * l + 1) / 2.0 * coef1(m) 
!     legen2(m) = (2 * l + 1) / 2.0 * coef2(m) 
!     legen3(m) = (2 * l + 1) / 2.0 * coef3(m) 
!     legen4(m) = (2 * l + 1) / 2.0 * coef4(m) 
!     if (legen(m) .gt. 1.0e-7) nlegen = l 
!   end do 
! 								    
!   return 
! end subroutine dda_db_liu
