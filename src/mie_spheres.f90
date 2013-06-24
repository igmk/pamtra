module mie_spheres

  use kinds
  use constants, only: pi,c, Im
  use settings, only: softsphere_adjust, lphase_flag, maxnleg
  use report_module
  use mie_scat_utlities  
  implicit none


  contains
  

  subroutine mie_spheres_wrapper(f, t,phase,&
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
    integer, intent(in) :: phase

    integer :: maxleg, nlegen, nbins
    real(kind=dbl) :: dia1, dia2
    real(kind=dbl), intent(in) :: ad, bd, alpha, gamma
    real(kind=dbl) :: a_mtox, bcoeff,tot_mass,wc
    real(kind=dbl) :: extinction, albedo, back_scatt, legen (200), legen2 (200),&
	legen3 (200), legen4 (200)                                        
    integer, parameter :: maxn  = 5000
    integer :: ir
    real(kind=dbl) :: del_d, density

    real(kind=dbl) :: distribution 
    real(kind=dbl), intent(out) :: back_spec(nbins+1)
    character :: aerodist * 1
    
    real(kind=dbl), intent(out) :: diameter(nbins+1)
    real(kind=dbl), dimension(nbins+1) :: particle_mass, density_vec
    real(kind=dbl), dimension(nbins+1) ::  ndens
    
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

    tot_mass = 0.
    do ir = 1, nbins+1
      !diameter(ir) is here the maximum extebd of the particle
      diameter(ir) = dia1 + (ir - 1) * del_d
      ndens(ir) = distribution(ad, bd, alpha, gamma, diameter(ir), aerodist)
      if ( (ir .eq. 1 .or. ir .eq. nbins+1) .and. nbins .gt. 0) then
	  ndens(ir) = 0.5 * ndens(ir) 
      end if
      particle_mass(ir) = a_mtox*diameter(ir)**bcoeff ! particle mass [kg]
      tot_mass = tot_mass + ndens(ir)*del_d*particle_mass(ir)
      if ((aerodist == "C") .or. (aerodist == "M")) then !MPACE!!
	if ((ir .eq. nbins+1) .and. (tot_mass/wc*100. .lt. 99.9d0)) then
	  ndens(ir) = ndens(ir) + (wc-tot_mass)/(del_d*particle_mass(ir))
	  tot_mass = wc
	end if
    density_vec(ir) = density
      end if	
	
	
    end do

    call calc_mie_spheres(err, f, t, phase, nbins, diameter, ndens, density_vec, &
      extinction, albedo, back_scatt, nlegen, legen,  &
      legen2, legen3, legen4, back_spec)    
        
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
      f, & ! frequency [GHz]
      t, &
      phase, &
      nbins, &
      diameter, &
      ndens, &
      density, &
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

    real(kind=dbl), intent(in) :: f  ! frequency [GHz]
    real(kind=dbl), intent(in) :: t    ! temperature [K]
    integer, intent(in) :: phase
    integer, intent(in) :: nbins
    real(kind=dbl), intent(in), dimension(nbins+1) :: diameter
    real(kind=dbl), intent(in), dimension(nbins+1) ::  ndens
    real(kind=dbl), intent(in), dimension(nbins+1) :: density
    
    real(kind=dbl), intent(out) :: extinction
    real(kind=dbl), intent(out) :: albedo
    real(kind=dbl), intent(out) :: back_scatt
    real(kind=dbl), intent(out), dimension(200) :: legen, legen2, legen3, legen4 
    real(kind=dbl), intent(out), dimension(nbins+1) :: back_spec
    integer, intent(out) :: nlegen

    real(kind=dbl) :: refre, refim
    real(kind=dbl) :: del_d
    real(kind=dbl) :: wavelength
    complex(kind=dbl) :: m_ice
                                       
    integer, parameter :: maxn  = 5000
    integer :: nterms, nquad, nmie, nleg 
    integer :: i, l, m, ir
    real(kind=dbl) :: x, tmp,diameter_eff, density_eff(nbins+1)
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
      f, & ! frequency [GHz]
      t, &
      phase, &
      nbins, &
       "diameter", diameter, &
      "ndens", ndens, &
      "density", density
   

      wavelength = c/(f*1.d9) !


	diameter_eff= diameter(nbins+1)
	density_eff(:) = density(:)


      
      x = pi * diameter_eff / wavelength
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

print*, "TODO: take real del_d del_d, remove diameter_eff and density_eff(ir) part"
    do ir = 1, nbins+1

      if (phase == 1) then
	call ref_water(0.d0, t-273.15, f, refre, refim, absind, abscof)
	msphere = refre-im*refim
      else if (phase == -1) then
	call ref_ice(t, f, refre, refim)
	if (density_eff(ir) == 917.d0) then
	  !ice sphere
	  msphere = refre-Im*refim 
	else
	  !softsphere
	  m_ice = refre-Im*refim  ! mimicking a
	  msphere = eps_mix((1.d0,0.d0),m_ice,density_eff(ir))
	end if
      else
	errorstatus = fatal
	print*,"phase=", phase
	msg = 'Did not understand variable phase'
	call report(errorstatus, msg, nameOfRoutine)
	return
      end if



      !diameter is here the maximum extebd of the particle
 
      if (ir <= nbins) then
	del_d = diameter(ir+1) - diameter(ir)
      else
	del_d = diameter(ir) - diameter(ir-1)
      end if
  
      nmie = 0 

      x = pi * diameter_eff / wavelength

	  

      if (verbose >= 4) print*, "density_eff(ir), diameter(ir), ndens(ir), msphere"
      if (verbose >= 4) print*, density_eff(ir), diameter(ir), ndens(ir), msphere

      call miecalc (err,nmie, x, msphere, a, b) 
      if (err /= 0) then
	  msg = 'error in mieclac!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
      end if         
      call miecross (nmie, x, a, b, qext, qscat, qback)

      ! sum up extinction, scattering, and backscattering as cross-sections/pi .pi is added in a later step
      qext =   qext  * ndens(ir) * (diameter_eff/2.d0)**2         ! [m²/m⁴]!
      qscat =  qscat * ndens(ir) * (diameter_eff/2.d0)**2        ! [m²/m⁴]!
      qback =  qback * ndens(ir) * (diameter_eff/2.d0)**2        !  [m²/m⁴]! cross section per volume per del_d
  
      !integrate=sum up . del_d is added at a later step!
      sumqe = sumqe + qext 
      sumqs = sumqs + qscat
      sumqback = sumqback + qback 

      back_spec(ir) =  qback * pi  ! volumetric backscattering corss section for radar simulator in [m²/m⁴]

      if (lphase_flag) then 
	  nmie = min0(nmie, nterms) 
	  do i = 1, nquad 
	    call mieangle (nmie, a, b, mu (i), p1, p2, p3, p4) 
	    sump1 (i) = sump1 (i) + p1 * ndens(ir) 
	    sump2 (i) = sump2 (i) + p2 * ndens(ir) 
	    sump3 (i) = sump3 (i) + p3 * ndens(ir) 
	    sump4 (i) = sump4 (i) + p4 * ndens(ir) 
	  end do
      end if
    end do

    !           multiply the sums by the integration delta and other constan
    !             put quadrature weights in angular array for later         

    extinction = pi * sumqe * del_d
    scatter = pi * sumqs * del_d
    back_scatt = pi * sumqback * del_d
    albedo = scatter / extinction 

    !         if the phase function is not desired then leave now           
    if ( .not. lphase_flag) return 

    tmp = (wavelength**2 / (pi * scatter) ) * del_d
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

    errorstatus = err    
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return 

  end subroutine calc_mie_spheres  
  
end module mie_spheres