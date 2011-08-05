subroutine dda_db_liu(f, t, npt, mindex, dia1, dia2, nbins, maxleg,   &
  ad, bd, alpha, gamma, lphase_flag, extinction, albedo, back_scatt,  &
  nlegen, legen, legen2, legen3, legen4, aerodist)
                  
! note that mindex has the convention with negative imaginary part      
! computes the mie scattering properties for a gamma or lognormal 
! distribution of spheres.
                                       
  use kinds
  use constants, only: pi,c
  use nml_params, only: verbose

  implicit none

  integer :: maxleg, nlegen, nbins

  logical :: lphase_flag

  integer, intent(in) :: npt    ! id of particle in db

  real(kind=dbl), intent(in) :: f,  &! frequency [GHz]
 				t    ! temperature [K]

  real(kind=dbl) :: wavelength, dia1, dia2
  real(kind=dbl) :: ad, bd, alpha, gamma 
  complex(kind=dbl) :: mindex 
  real(kind=dbl) :: extinction, albedo, back_scatt
  real(kind=dbl) :: legen(200), legen2(200), legen3(200), legen4(200)
  integer, parameter :: maxn = 5000
  integer :: nterms, nquad, nmie, nleg 
  integer :: i, l, m, ir
  real(kind=dbl) :: x, del_d,  ndens, tmp, fac
  real(kind=dbl) :: qext, qscat, qback, scatter 
  real(kind=dbl) :: distribution 
!  real(kind=dbl) :: mu(maxn), wts(maxn)
  real(kind=dbl) :: p1, pl, pl1, pl2, p2, p3, p4 
  real(kind=dbl) :: sumqe, sumqs, sumqback 
!  real(kind=dbl) :: sump1(maxn), coef1(maxn), sump2(maxn), coef2(maxn),   &
!    sump3(maxn), coef3(maxn), sump4(maxn), coef4(maxn)
  real(kind=dbl), allocatable, dimension(:) :: sump1, coef1, sump2, coef2,   &
    sump3, coef3, sump4, coef4
  complex(kind=dbl), dimension(maxn) :: a, b
  complex(kind=dbl) :: msphere
  character :: aerodist*1 

  real(kind=dbl), dimension(nbins) :: xi, weights

  real(kind=dbl) :: ntot

! variables communicating with the database. need to be single precission !!!

  integer :: iret, is_loaded
  integer, parameter :: nang_db = 37
  real(4) :: t_liu, f_liu
  real(4) :: diameter
  real(4) :: abs_liu,sca_liu,bsc_liu,g, r_ice_eq
  real(4), dimension(nang_db) :: p_liu
  real(4), dimension(nang_db) :: ang_db

  real(kind=dbl), allocatable, dimension(:) :: ang_quad, mu, wts, P1_quad

  if (verbose .gt. 1) print*, 'Entering dda_db_liu'

  do i = 0,36
   ang_db(i+1) = i*5.
  end do

  t_liu = real(t)   ! convert to 4byte real (database requirement!)
  f_liu = real(f)   ! convert to 4byte real (database requirement!)
  is_loaded = 0

  qext = 0.d0;  qscat = 0.d0;  qback = 0.d0
  nlegen = 0
  legen = 0.d0;  legen2 = 0.d0;  legen3 = 0.d0;  legen4 = 0.d0

  wavelength = c/(f*1.e9)

  ntot = 0.d0

  if (t_liu .lt. 234.) then
    print*, "temperature to low for database"
  	t_liu = 234. ! lowest value in database
  end if

  if (t_liu .gt. 273.15) then
    print*, "temperature to high for database"
  	t_liu = 273.15 ! highest value in database
  end if
								    
  ! find the maximum number of terms required in the mie series
  ! probably this should be the mass equivalent sphere diameter
  if (verbose .gt. 2) print*, 'find maximum number for the mie series'
  msphere = mindex
  x = pi * dia2 / wavelength
  nterms = 0 
  call miecalc(nterms, x, msphere, a, b) ! miecalc returns nterms
  nlegen = 2 * nterms 
  nlegen = min(maxleg, nlegen) 
  nquad = (nlegen + 2 * nterms + 2) / 2 
  if (nquad .gt. maxn) stop 'mie: maxn exceeded' 

  if (verbose .gt. 2) print*, 'allocating phase function variables'

  allocate(ang_quad(nquad),mu(nquad),wts(nquad),P1_quad(nquad))
  allocate(sump1(nquad), coef1(nquad), sump2(nquad), coef2(nquad))
  allocate(sump3(nquad), coef3(nquad), sump4(nquad), coef4(nquad))

  sumqe = 0.0d0
  sumqs = 0.0d0
  sumqback = 0.0d0
  sump1 = 0.d0; sump2 = 0.d0; sump3 = 0.d0; sump4 = 0.d0
  coef1 = 0.d0;  coef2 = 0.d0;  coef3 = 0.d0;  coef4 = 0.d0

  if (verbose .gt. 2) print*, '... done'

  !  get the Gauss-Legendre quadrature abscissas and weights for the phase function integration
  call gausquad(nquad, mu, wts)

  !  get the Gauss-Legendre quadrature abscissas and weights for abs, bsc, and ext
  call gausquad(nbins, xi , weights)

  if (verbose .gt. 2) print*, 'done calling quadratures'

								    
  !   integration loop over radius of spheres 
  if (nbins .gt. 0) del_d = (dia2 - dia1) / nbins

  if (verbose .gt. 3) print*, 'Doing database search for diameter intervall: '
  do ir = 1, nbins+1
    diameter = dia1 + (ir - 1) * del_d
!	diameter = (dia2-dia1)/2.d0*xi(ir)+(dia1+dia2)/2.d0
    ndens = distribution(ad, bd, alpha, gamma, dble(diameter), aerodist)  ! number density
    if ((ir .eq. 1 .or. ir .eq. nbins+1) .and. nbins .gt. 0) then
		ndens = 0.5d0 * ndens
    end if 
    ntot = ntot + ndens*del_d
  if (verbose .gt. 2) print*, ir, ' with: ',f_liu,t_liu,npt,diameter*1.e6
    call scatdb(f_liu,t_liu,npt,diameter*1.e6,abs_liu,sca_liu,bsc_liu,g,p_liu,r_ice_eq,iret,is_loaded)
  if (verbose .gt. 2) print*, 'got: ',iret, abs_liu,sca_liu,bsc_liu,g

	if (iret .ne. 0) print*, iret,f_liu,t_liu,npt,diameter*1.e6, abs_liu,sca_liu,bsc_liu

!    qext = (abs_liu+sca_liu)/(pi*(diameter/.2)**2)
!    qscat = sca_liu/(pi*(diameter/.2)**2)
!    qback = bsc_liu/(pi*(diameter/.2)**2)
!  write(36,*) pi*diameter/wavelength, abs_liu/(pi*(diameter/.2)**2),qscat,qback
    qext = (abs_liu+sca_liu)/(pi*(r_ice_eq*1.e-6)**2)
    qscat = sca_liu/(pi*(r_ice_eq*1.e-6)**2)
    qback = bsc_liu/(pi*(r_ice_eq*1.e-6)**2)
    ! sum up extinction, scattering, and backscattering as cross-sections/pi
!    sumqe = sumqe + qext * ndens * (diameter/.2)**2
!    sumqs = sumqs + qscat * ndens * (diameter/.2)**2
!    sumqback = sumqback + qback * ndens * (diameter/.2)**2
    sumqe = sumqe + qext * ndens * (r_ice_eq*1.e-6)**2
    sumqs = sumqs + qscat * ndens * (r_ice_eq*1.e-6)**2
    sumqback = sumqback + qback * ndens * (r_ice_eq*1.e-6)**2
    if (lphase_flag) then
      ang_quad = acos(mu(nquad:1:-1))*180.d0/pi
	  call interpolation(nang_db,nquad,dble(ang_db),dble(p_liu),ang_quad,P1_quad)
	  fac = sum(P1_quad*wts)
	  P1_quad = P1_quad*(2.d0/fac)
      do i = 1, nquad 
		sump1(i) = sump1(i) + P1_quad(i) * ndens * del_d
!		sump2(i) = sump2(i) + p2 * ndens
!		sump3(i) = sump3(i) + p3 * ndens
!		sump4(i) = sump4(i) + p4 * ndens
      end do 
    end if 
  end do 
 
  !   multiply the sums by the integration delta and other constants
  !   put quadrature weights in angular array for later usage    
  if (nbins .eq. 0) del_d = 1.0d0

  extinction = pi * sumqe * del_d      ! [1/m]
  scatter = pi * sumqs * del_d         ! [1/m]
  back_scatt = pi * sumqback * del_d   ! [1/m]
  albedo = scatter / extinction 

  ! if the phase function is not desired then leave now           
  if ( .not. lphase_flag) return 
								    
  !tmp = (wavelength**2 / (pi * scatter)) * del_d
  tmp = 1.d0/ntot
  do i = 1, nquad 
    sump1(i) = tmp * sump1(i) * wts(i) 
!    sump2(i) = tmp * sump2(i) * wts(i)
!    sump3(i) = tmp * sump3(i) * wts(i)
!    sump4(i) = tmp * sump4(i) * wts(i)
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
!      coef2(m) = coef2(m) + sump2(i) * pl
!      coef3(m) = coef3(m) + sump3(i) * pl
!      coef4(m) = coef4(m) + sump4(i) * pl
      pl2 = pl1 
      pl1 = pl 
    end do 
  end do 
  nleg = nlegen 
  do l = 0, nleg 
    m = l + 1 
    legen(m) = (2 * l + 1) / 2.0 * coef1(m) 
!    legen2(m) = (2 * l + 1) / 2.0 * coef2(m)
!    legen3(m) = (2 * l + 1) / 2.0 * coef3(m)
!    legen4(m) = (2 * l + 1) / 2.0 * coef4(m)
    if (legen(m) .gt. 1.0e-7) nlegen = l 
  end do 

  deallocate(ang_quad,mu,wts,P1_quad)
  deallocate(sump1, coef1, sump2, coef2)
  deallocate(sump3, coef3, sump4, coef4)

  return
end subroutine dda_db_liu
