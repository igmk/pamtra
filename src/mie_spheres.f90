module mie_spheres

  use kinds
  use constants, only: pi,c, Im
  use settings, only: lphase_flag, maxnleg, lhyd_absorption
  use report_module
  use mie_scat_utilities  
  use vars_index, only: i_x,i_y, i_z, i_h
  implicit none


  contains

  subroutine calc_mie_spheres(&
      errorstatus, &
      freq, & ! frequency [Hz]
      t, &
      liq_ice, &
      nbins, &
      diameter, &
      del_d, &
      ndens, & !normed
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
    real(kind=dbl), intent(in), dimension(nbins) :: diameter
    real(kind=dbl), intent(in), dimension(nbins) :: del_d    
    real(kind=dbl), intent(in), dimension(nbins) ::  ndens
    real(kind=dbl), intent(in), dimension(nbins) :: density
    real(kind=dbl), intent(in) :: refre
    real(kind=dbl), intent(in) :: refim !positive(?)

    real(kind=dbl), intent(out) :: extinction
    real(kind=dbl), intent(out) :: albedo
    real(kind=dbl), intent(out) :: back_scatt
    real(kind=dbl), intent(out), dimension(200) :: legen, legen2, legen3, legen4 
    real(kind=dbl), intent(out), dimension(nbins) :: back_spec
    integer, intent(out) :: nlegen

    real(kind=dbl) :: wavelength
    complex(kind=dbl) :: m_ice
    real(kind=dbl) :: del_d_eff, ndens_eff, n_tot
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
    character(len=30) :: nameOfRoutine = 'calc_mie_spheres'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    err = 0

    ! initialize intent(out) variables with 0
    extinction = 0.0d0
    albedo = 0.0d0
    back_scatt = 0.0d0
    legen(:) = 0.0d0
    legen2(:) = 0.0d0
    legen3(:) = 0.0d0
    legen4(:) = 0.0d0
    back_spec(:) = 0.0d0
    nlegen = 0
    if (verbose >= 4) print*, "calc_mie_spheres",&
!     if ((liq_ice == -1) .and. (nbins == 50)) print*, "calc_mie_spheres",&
      i_x,i_y, i_z, i_h, &
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
      call assert_true(err,SUM(ndens)>0,&
          "sum(ndens) must be greater zero")    
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

      !initialize
      back_spec(:) = 0.d0

      wavelength = c/(freq) !

      if ((liq_ice == -1) .and. (density(1) /= 917.d0)) then
          m_ice = refre-Im*refim  ! mimicking a
          msphere = eps_mix((1.d0,0.d0),m_ice,density(nbins))
      else
                msphere = refre-im*refim
      end if
     
      x = pi * diameter(nbins) / wavelength
      nterms = 0 
      n_tot = 0.d0
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
      sump1 (:) = 0.0d0 
      sump2 (:) = 0.0d0 
      sump3 (:) = 0.0d0 
      sump4 (:) = 0.0d0 


    do ir = 1, nbins
      !Do not process if no particles present
      if (ndens(ir) <= 0) CYCLE

      ndens_eff = ndens(ir)
      del_d_eff = del_d(ir)

      n_tot = n_tot + (ndens_eff * del_d_eff)

      if ((liq_ice == -1) .and. (density(ir) /= 917.d0)) then
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
      
      if (verbose >= 4) print*, "ir, density(ir), diameter(ir), ndens_eff, del_d_eff, msphere, x"
      if (verbose >= 4) print*, ir, density(ir), diameter(ir), ndens_eff, del_d_eff, msphere, x 
      
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
  
      if (verbose >= 4) print*, "qback* del_d_eff, ndens_eff , (diameter(ir)/2.d0), pi, del_d_eff"
      if (verbose >= 4) print*, qback * del_d_eff, ndens_eff ,(diameter(ir)/2.d0), pi, del_d_eff
  
      !integrate=sum up . del_d is already included in ndens, since ndens is not normed!
      sumqe = sumqe + ( qext * del_d_eff)
      sumqs = sumqs + ( qscat * del_d_eff)
      sumqback = sumqback + ( qback * del_d_eff)
      if (verbose >= 4) print*, diameter(ir), ndens_eff, del_d_eff, n_tot, sumqe, sumqs, sumqback

      if (verbose >= 4) print*, "NEW: sumqback, sumqs, sumqe"
      if (verbose >= 4) print*,  sumqback , sumqs, sumqe
      back_spec(ir) =  qback   ! volumetric backscattering cross section for radar simulator in backscat per volume per del_d[m²/m⁴]

      if (lphase_flag) then 
	  nmie = min0(nmie, nterms) 
	  do i = 1, nquad 
	    call mieangle (nmie, a, b, mu (i), p1, p2, p3, p4) 
	    sump1 (i) = sump1 (i) + p1 * ndens_eff * del_d_eff ! = M11 = M22
	    sump2 (i) = sump2 (i) + p2 * ndens_eff * del_d_eff ! = M12 = M21
	    sump3 (i) = sump3 (i) + p3 * ndens_eff * del_d_eff ! = M33 = M44
	    sump4 (i) = sump4 (i) + p4 * ndens_eff * del_d_eff ! = M34 = -M43
	  end do
      end if
    end do

    !           multiply the sums by the integration delta and other constan
    !             put quadrature weights in angular array for later         

    if (verbose >= 4) print*, "ntot", n_tot

        
        
    if (lhyd_absorption) then
        extinction = sumqe 
    else
        extinction = sumqe - sumqs !remove scattering from extinction
    end if
    scatter = sumqs 
    back_scatt = sumqback 
    albedo = scatter / extinction 
      if (verbose >= 4) print*, sumqe, sumqs, sumqback, n_tot
      if (verbose >= 4) print*, "extinction, scatter, back_scatt, albedo"
      if (verbose >= 4) print*,  extinction, scatter, back_scatt, albedo
       
    
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