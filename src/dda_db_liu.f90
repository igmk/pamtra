subroutine dda_db_liu(f, t, mindex, dia1, dia2, nbins, maxleg,   &
     ad, bd, alpha, gamma, lphase_flag, extinction, albedo, back_scatt,  &
     nlegen, legen, legen2, legen3, legen4, aerodist)

  ! note that mindex has the convention with negative imaginary part      
  ! computes the mie scattering properties for a gamma or lognormal 
  ! distribution of spheres.

  use kinds
  use constants, only: pi,c
  use nml_params, only: verbose, liu_type, data_path

  implicit none

  integer, intent(in) :: maxleg, nbins
  integer, intent(out) :: nlegen
  logical, intent(in) :: lphase_flag

  real(kind=dbl), intent(in) :: f,  &! frequency [GHz]
       t    ! temperature [K]

  real(kind=dbl), intent(in) :: dia1, dia2
  real(kind=dbl), intent(in) :: ad, bd, alpha, gamma
  real(kind=dbl) :: wavelength,mass_eq_dia
  complex(kind=dbl), intent(in) :: mindex 
  real(kind=dbl), intent(out) :: extinction, albedo, back_scatt
  real(kind=dbl), intent(out) :: legen(200), legen2(200), legen3(200), legen4(200)
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
  character, intent(out) :: aerodist*1 

  real(kind=dbl), dimension(nbins) :: xi, weights

  real(kind=dbl) :: ntot, swc

  ! variables communicating with the database. need to be single precission !!!

  integer :: iret, is_loaded
  integer, parameter :: nang_db = 37
  real(kind=sgl) :: t_liu, f_liu
  real(kind=sgl) :: diameter
  real(kind=sgl) :: abs_liu,sca_liu,bsc_liu,g, r_ice_eq
  real(kind=sgl), dimension(nang_db) :: p_liu
  real(kind=sgl), dimension(nang_db) :: ang_db
  real(kind=sgl), parameter, dimension(0:10) :: dmax = (/4835.,3304.,2632.,3246.,5059.,&
       10000.,10000.,10000.,10000.,10000.,12454./)*1.e-6
  real(kind=sgl), parameter, dimension(0:10) :: dmin = (/121.,83.,66.,81.,127.,&
       50.,50.,50.,50.,50.,75./)*1.e-6
  real(kind=sgl), parameter, dimension(0:10) :: a_m = (/33.996,106.526,210.529,112.290,29.680,&
       0.183,0.1287,0.1680,0.2124,1.191e-3,5.666e-3/)
  real(kind=sgl), parameter, dimension(0:10) :: b_m = (/3.,3.,3.,3.,3.,&
       2.274,2.264,2.274,2.285,1.511,1.82/)

  real(kind=dbl), allocatable, dimension(:) :: ang_quad, mu, wts, P1_quad

  if (verbose .gt. 1) print*, 'Entering dda_db_liu'

  ! initialize database phase function angles
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
     if (verbose .gt. 1) print*, "temperature to low for database"
     t_liu = 234. ! lowest value in database
  end if

  if (t_liu .gt. 273.15) then
     if (verbose .gt. 1) print*, "temperature to high for database"
     t_liu = 273.15 ! highest value in database
  end if

  ! find the maximum number of terms required in the mie series
  ! probably this should be the mass equivalent sphere diameter
  if (verbose .gt. 2) print*, 'find maximum number for the mie series'
  msphere = mindex
  mass_eq_dia = (6.*a_m(liu_type)*dia2**b_m(liu_type)/(pi*917.))**(1./3.)
 print *,liu_type,mass_eq_dia
! mass_eq_dia =
  !  x = pi * dia2 / wavelength
  x = pi * mass_eq_dia / wavelength
  nterms = 0 

  call miecalc(nterms, x, msphere, a, b) ! miecalc returns nterms
! print *,nterms
! if (nterms .lt. 10) nterms = 10
! print *,nterms
! print*,"nterms",nterms
! print*,"x",x
! print*,"msphere",msphere
! print*,"a",a(1:nterms+1)
! print*,"b",b(1:nterms+1)

!    215.05250000000001     
!  nterms           7
!  x   1.3706891306673277     
!  msphere (  1.7771803954058019     ,-3.73910886490081188E-004)
!  a ( 0.33399466247654391     ,-0.47135633586867320     ) ( 4.84774208805761019E-003,-6.92517345487174696E-002) ( 1.15153713481139045E-005,-3.22475391848596905E-003) ( 3.84789862938517094E-008,-9.19563066392508683E-005) ( 5.48400644088990285E-010,-1.71554630228806083E-006) ( 7.00219050838030447E-012,-2.23577460569478837E-008) ( 6.65132936506777237E-014,-2.14488189800856400E-010) (  0.0000000000000000     ,  0.0000000000000000     )
!  b ( 0.12482916528744245     ,-0.32986235763903937     ) ( 1.80592715518353099E-004,-1.30439488340791648E-002) ( 3.74357707100429803E-007,-3.54988066902067962E-004) ( 4.43089698259088587E-009,-6.58945399625693957E-006) ( 5.59674417804745679E-011,-8.61476221430638550E-008) ( 5.31775508880886119E-013,-8.30514194861544992E-010) ( 3.89219326418746779E-015,-6.13480922417889293E-012) (  0.0000000000000000     ,  0.0000000000000000     )

!    215.05250000000001     
! nterms   =       8
!  x   2.9004083825932123     
!  msphere (  1.7771803954058019     ,-3.73910886490081188E-004)
! a(1:11)=(/CMPLX(0.51689001598900541,0.49858772416940417),CMPLX(0.99822832846958609,2.42751773546725294E-002),&
! CMPLX(0.29194416867455830,-0.45412626571685111),CMPLX(3.32891096256099838E-003,-5.73564266202203077E-002),&
! CMPLX(2.18389679251775356E-005,-4.48707013838365188E-003),CMPLX(1.61152560358732911E-007,-2.64463284347134209E-004),&
! CMPLX(3.97775003444399327E-009,-1.16616121824232274E-005),CMPLX(1.26749209754190923E-010,-3.94337841882554094E-007),&
! CMPLX(3.32036221889263258E-012,-1.05126699906402068E-008),CMPLX(7.06982997341716537E-014,-2.26347837158777329E-010),&
! CMPLX(0.0000000000000000,0.0000000000000000)/)
! b(1:11)=(/CMPLX(0.92729840594756863,0.25840379142207026),CMPLX(0.74148769509971280,0.43709128097481809),&
! CMPLX(0.48798170365440957,-0.49764281437443714),CMPLX(6.88705776535616151E-004,-2.56575535806836848E-002),&
! CMPLX(2.65242116661926677E-006,-1.25037531402793882E-003),CMPLX(4.24037023424217557E-008,-5.14812235979604160E-005),&
! CMPLX(1.22302719003843511E-009,-1.68482277526771416E-006),CMPLX(3.07206158788031366E-011,-4.41133370762673934E-008),&
! CMPLX(6.37475175976837969E-013,-9.39635036180452891E-010),CMPLX(1.10312812806734579E-014,-1.65613659102433772E-011),&
! CMPLX(0.0000000000000000,0.0000000000000000)/)


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
  !  do ir = 1, nbins+1
  do ir = 1, nbins
     !    diameter = dia1 + (ir - 1) * del_d
     diameter = (dia2-dia1)/2.d0*xi(ir)+(dia1+dia2)/2.d0
     ndens = distribution(ad, bd, alpha, gamma, dble(diameter), aerodist)  ! number density
     !    if ((ir .eq. 1 .or. ir .eq. nbins+1) .and. nbins .gt. 0) then
     !		ndens = 0.5d0 * ndens
     !    end if
     ntot = ntot + ndens*weights(ir)

     if (diameter .gt. dmax(liu_type)) diameter = dmax(liu_type)
     if (diameter .lt. dmin(liu_type)) diameter = dmin(liu_type)
     !Why does this line cause make test2 and/or three to fail, if verbose>2???
!      if (verbose .gt. 2) print*, ir, ' with: ',f_liu,t_liu,liu_type,diameter*1.e6
print*,"HH",trim(data_path)
     call scatdb(f_liu,t_liu,liu_type,diameter*1.e6,abs_liu,sca_liu,bsc_liu,g,p_liu,r_ice_eq,iret,is_loaded,&
          trim(data_path))

     !   if (verbose .gt. 2) print*, 'got: ',iret, abs_liu,sca_liu,bsc_liu,g
     if (verbose .gt. 1) print*, iret,f_liu,t_liu,liu_type,diameter*1.e6, abs_liu,sca_liu,bsc_liu

     qext = (abs_liu+sca_liu)
     qscat = sca_liu
     qback = bsc_liu
     ! sum up extinction, scattering, and backscattering as cross-sections/pi

     sumqe = sumqe + qext * ndens * weights(ir)
     sumqs = sumqs + qscat * ndens * weights(ir)
     sumqback = sumqback + qback * ndens * weights(ir)
     if (lphase_flag) then
        ang_quad = acos(mu(nquad:1:-1))*180.d0/pi
        call interpolation(nang_db,nquad,dble(ang_db),dble(p_liu),ang_quad,P1_quad)
        fac = sum(P1_quad*wts)
        P1_quad = P1_quad*(2.d0/fac)
        do i = 1, nquad 
           sump1(i) = sump1(i) + P1_quad(i) * ndens  * weights(ir)!* del_d
           !		sump2(i) = sump2(i) + p2 * ndens
           !		sump3(i) = sump3(i) + p3 * ndens
           !		sump4(i) = sump4(i) + p4 * ndens
        end do
     end if
  end do
  sumqe = sumqe*(dia2-dia1)/2.d0
  sumqs = sumqs*(dia2-dia1)/2.d0
  sumqback = sumqback*(dia2-dia1)/2.d0
  sump1 = sump1*(dia2-dia1)/2.d0
  ntot = ntot*(dia2-dia1)/2.d0
  !   multiply the sums by the integration delta and other constants
  !   put quadrature weights in angular array for later usage    
  if (nbins .eq. 0) del_d = 1.0d0

  extinction = sumqe !* del_d      ! [1/m]
  scatter = sumqs !* del_d       ! [1/m]
  back_scatt = sumqback !* del_d   ! [1/m]
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
