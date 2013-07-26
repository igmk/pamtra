module mie_spheres

  use kinds
  use constants, only: pi,c, Im
  use settings, only: softsphere_adjust, lphase_flag, maxnleg
  use report_module
  use mie_scat_utlities  
  implicit none


  contains
  

  subroutine mie_spheres_wrapper(f, t,liq_ice,&
      a_mtox, bcoeff, dia1, dia2, nbins, maxleg, ad, bd, alpha, &
      gamma, extinction, albedo, back_scatt, nlegen, legen,  &
      legen2, legen3, legen4, aerodist,density,wc,&
      diameter, back_spec)
    !    computing the scattering properties according to                  
    !    ice sphere model, i.e. the electromagnetic properties of the      
    !     particle are computed by assuming that they are the same          
    !     as the equivalent mass sphere
    !                                     
    ! note that mindex has the convention with negative imaginary part      
    !out
    !diameter: diameter spectrum [m]
    !back_spec: backscattering cross section per volume per del_d [m²/m⁴]

    implicit none

    real(kind=dbl), intent(in) :: f,  &! frequency [GHz]
	t    ! temperature [K]
    integer, intent(in) :: liq_ice

    integer :: maxleg, nlegen, nbins
    real(kind=dbl) :: dia1, dia2
    real(kind=dbl), intent(in) :: ad, bd, alpha, gamma
    real(kind=dbl) :: a_mtox, bcoeff,wc
    real(kind=dbl) :: extinction, albedo, back_scatt, legen (200), legen2 (200),&
	legen3 (200), legen4 (200)                                        
    integer, parameter :: maxn  = 5000
    integer :: ir
    real(kind=dbl) :: del_d, density
    
    real(kind=dbl) :: distribution 
    real(kind=dbl), intent(out) :: back_spec(nbins)
    character :: aerodist * 1
    
    real(kind=dbl), intent(out) :: diameter(nbins)
    real(kind=dbl), dimension(nbins) :: particle_mass, density_vec
    real(kind=dbl), dimension(nbins) ::  ndens
    
      integer(kind=long) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=40) :: nameOfRoutine = 'mie_spheres_wrapper'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

    if (nbins .gt. 0) then
      del_d = (dia2 - dia1) / nbins
    else
      del_d = 1.d0
    end if

    do ir = 1, nbins
      !diameter(ir) is here the maximum extebd of the particle
      diameter(ir) = dia1 + (ir - 1) * del_d
      ndens(ir) = distribution(ad, bd, alpha, gamma, diameter(ir), aerodist)
      particle_mass(ir) = a_mtox*diameter(ir)**bcoeff ! particle mass [kg]
    density_vec(ir) = density
	
    end do

!     call calc_mie_spheres(err, f*1d9, t, liq_ice, nbins, diameter, ndens, density_vec, &
!       extinction, albedo, back_scatt, nlegen, legen,  &
!       legen2, legen3, legen4, back_spec)    
   stop "TODO: add del_d to wrapper"     
    if (err /= 0) then
	msg = 'error in calc_mie_spheres!'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	stop !return
    end if        
        
    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return 

  end subroutine mie_spheres_wrapper
  

  subroutine calc_mie_spheres(&
      errorstatus, &
      freq, & ! frequency [Hz]
      t, &
      liq_ice, &
      nbins, &
      diameter, &
      del_d, &
      ndens, & !ndens NOT NORMED PER del_d!
      density, &
      refre, &
      refim, &
      extinction, &
      albedo, &
      back_scatt, &
      nlegen, &
      legen, &
      legen2, &
      legen3, &
      legen4, &
      back_spec )    
    !    computing the scattering properties according to                  
    !    ice sphere model, i.e. the electromagnetic properties of the      
    !     particle are computed by assuming that they are the same          
    !     as the equivalent mass sphere
    !                                     
    ! note that mindex has the convention with negative imaginary part      
    !out
    !diameter: diameter spectrum [m]
    !back_spec: backscattering cross section per volume per del_d [m²/m⁴]

    implicit none

    real(kind=dbl), intent(in) :: freq  ! frequency [Hz]
    real(kind=dbl), intent(in) :: t    ! temperature [K]
    integer, intent(in) :: liq_ice
    integer, intent(in) :: nbins
    real(kind=dbl), intent(in), dimension(nbins+1) :: diameter
    real(kind=dbl), intent(in), dimension(nbins) :: del_d    
    real(kind=dbl), intent(in), dimension(nbins+1) ::  ndens
    real(kind=dbl), intent(in), dimension(nbins+1) :: density
    real(kind=dbl), intent(in) :: refre
    real(kind=dbl), intent(in) :: refim !positive(?)

    real(kind=dbl), intent(out) :: extinction
    real(kind=dbl), intent(out) :: albedo
    real(kind=dbl), intent(out) :: back_scatt
    real(kind=dbl), intent(out), dimension(200) :: legen, legen2, legen3, legen4 
    real(kind=dbl), intent(out), dimension(nbins+1) :: back_spec
    integer, intent(out) :: nlegen

    real(kind=dbl) :: wavelength
    complex(kind=dbl) :: m_ice
    real(kind=dbl) :: del_d_eff, ndens_eff    
    integer, parameter :: maxn  = 5000
    integer :: nterms, nquad, nmie, nleg 
    integer :: i, l, m, ir
    real(kind=dbl) :: x, tmp,diameter_eff
    real(kind=dbl) :: qext, qscat, qback, scatter 
    real(kind=dbl) :: mu(maxn), wts(maxn)
    real(kind=dbl) :: p1, pl, pl1, pl2, p2, p3, p4 
    real(kind=dbl) :: sumqe, sumqs, sumqback
    real(kind=dbl), dimension(maxn) :: sump1, coef1, sump2, coef2,   &
	sump3, coef3, sump4, coef4          
    real(kind=dbl) :: absind, abscof
    complex(kind=dbl), dimension(maxn) :: a, b
    complex(kind=dbl) :: msphere, eps_mix

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'calc_mie_spheres'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

    if (verbose >= 4) print*, "calc_mie_spheres(",&
      errorstatus, &
      freq, & ! frequency [Hz]
      t, &
      liq_ice, &
      nbins, &
       "diameter", diameter, &
      "ndens", ndens, &
      "density", density
   

      call assert_true(err,all(density>=0),&
          "density must be positive")  
      call assert_true(err,all(ndens>=0),&
          "ndens must be positive")   
      call assert_true(err,all(diameter>0),&
          "diameter must be positive")   
      call assert_true(err,all(del_d>0),&
          "del_d must be positive")   
      call assert_true(err,(nbins>0),&
          "nbins must be positive")   
      call assert_true(err,(freq>0),&
          "freq must be positive")   
      call assert_true(err,(t>0),&
          "t must be positive")   
      if (err > 0) then
          errorstatus = fatal
          msg = "assertation error"
          call report(errorstatus, msg, nameOfRoutine)
          return
      end if    

     
      wavelength = c/(freq) !

     
      x = pi * diameter(1) / wavelength
      nterms = 0 
      call miecalc (err,nterms, x, msphere, a, b) 
      if (err /= 0) then
	  msg = 'error in mieclac!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
      end if         
      nlegen = 2 * nterms 
      nlegen = min(maxnleg, nlegen) 
      nquad = (nlegen + 2 * nterms + 2) / 2 
      if (nquad.gt.maxn) then
	  errorstatus = fatal
	  msg = 'mie: maxn exceeded' 
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      end if
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

    do ir = 1, nbins+1

      if (ir == 1) then
        ndens_eff = ndens(1)/2.d0
        del_d_eff = del_d(1)
      else if (ir == nbins+1) then
        ndens_eff = ndens(nbins+1)/2.d0
        del_d_eff = del_d(nbins)
      else
        ndens_eff = ndens(ir)
        del_d_eff = del_d(ir)
      end if

      if (liq_ice == -1 .and. density(ir) /= 917.d0) then
	  m_ice = refre-Im*refim  ! mimicking a
	  msphere = eps_mix((1.d0,0.d0),m_ice,density(ir))
      else
		msphere = refre-im*refim
      end if
  
      nmie = 0 

      x = pi * diameter(ir) / wavelength

      call miecalc (err,nmie, x, msphere, a, b) 
      if (err /= 0) then
	  msg = 'error in mieclac!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
      end if         
      
      if (verbose >= 0) print*, "density(ir), diameter(ir), ndens_eff, msphere, x"
      if (verbose >= 0) print*, density(ir), diameter(ir), ndens_eff, msphere, x 
      
      call miecross (nmie, x, a, b, qext, qscat, qback)
      
      if (verbose >= 4) print*, "qext, qscat, qback"
      if (verbose >= 4) print*, qext, qscat, qback
      !from efficiencies cross sections
      qext =   qext   * (diameter(ir)/2.d0)**2 *pi        ! [m²]!
      qscat =  qscat  * (diameter(ir)/2.d0)**2 *pi       ! [m²]!
      qback =  qback  * (diameter(ir)/2.d0)**2 *pi       !  [m²]! cross section
   
      ! apply bin weights
      qext =   qext  * ndens_eff      ! [m²/m⁴]!
      qscat =  qscat * ndens_eff      ! [m²/m⁴]!
      qback =  qback * ndens_eff      !  [m²/m⁴]! cross section per volume
  
      if (verbose >= 4) print*, "qback * ndens_eff * (diameter(ir)/2.d0), pi, del_d_eff"
      if (verbose >= 4) print*, qback , ndens_eff ,(diameter(ir)/2.d0), pi, del_d_eff
  
      !integrate=sum up . del_d is already included in ndens, since ndens is not normed!
      sumqe = sumqe + ( qext * del_d_eff)
      sumqs = sumqs + ( qscat * del_d_eff)
      sumqback = sumqback + ( qback * del_d_eff)

      back_spec(ir) =  qback   ! volumetric backscattering corss section for radar simulator in backscat per volume per del_d[m²/m⁴]

      if (lphase_flag) then 
	  nmie = min0(nmie, nterms) 
	  do i = 1, nquad 
	    call mieangle (nmie, a, b, mu (i), p1, p2, p3, p4) 
	    sump1 (i) = sump1 (i) + p1 * ndens_eff * del_d_eff
	    sump2 (i) = sump2 (i) + p2 * ndens_eff * del_d_eff
	    sump3 (i) = sump3 (i) + p3 * ndens_eff * del_d_eff
	    sump4 (i) = sump4 (i) + p4 * ndens_eff * del_d_eff
	  end do
      end if
    end do

    !           multiply the sums by the integration delta and other constan
    !             put quadrature weights in angular array for later         

    extinction = sumqe 
    scatter = sumqs 
    back_scatt = sumqback 
    albedo = scatter / extinction 

      if (verbose >= 4) print*,"ir, extinction, scatter, back_scatt, albedo"
      if (verbose >= 4) print*,ir, extinction, scatter, back_scatt, albedo
       
    
    ! if the liq_ice function is not desired then leave now           
    if ( .not. lphase_flag) return 

    tmp = (wavelength**2 / (pi * scatter) ) 
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

      call assert_false(err,any(isnan(legen)),&
          "nan in legen")
      call assert_false(err,any(isnan(legen2)),&
          "nan in legen2")
      call assert_false(err,any(isnan(legen3)),&
          "nan in legen3")
      call assert_false(err,any(isnan(legen4)),&
          "nan in legen4")
      call assert_false(err,any(isnan(back_spec)),&
          "nan in back_spec")   
      call assert_true(err,(extinction>0),&
          "extinction must be positive")   
      call assert_true(err,(scatter>0),&
          "scatter must be positive") 
      call assert_true(err,(back_scatt>0),&
          "back_scatt must be positive") 
      if (err > 0) then
          errorstatus = fatal
          msg = "assertation error"
          call report(errorstatus, msg, nameOfRoutine)
          return
      end if    

    errorstatus = err    
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return 

  end subroutine calc_mie_spheres  
  
end module mie_spheres