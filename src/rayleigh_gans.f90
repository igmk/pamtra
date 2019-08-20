module rayleigh_gans

  use kinds
  use constants, only: pi,c, Im, rho_ice
  use settings, only: active,passive, nummu, mu_values, nstokes, quad_weights, maxnleg, lhyd_absorption, lphase_flag
  use report_module
  use mie_scat_utilities  
  use vars_index, only: i_x,i_y, i_z, i_h
  implicit none


  contains

  subroutine calc_self_similar_rayleigh_gans_rt3(&
      errorstatus, &
      freq, & 
      liq_ice, & ! again... not sure
      nbins, &
      diameter, &
      del_d, &
      ndens, & ! in calc_mie_spheres this is apparently normed
      mass, &  ! in mie I have density
      as_ratio,& 
      canting,&
      refre, &
      refim, &
      rg_kappa, &
      rg_beta, &
      rg_gamma, &
      rg_zeta, &
      extinction, &
      albedo, &
      back_scatt, &
      nlegen, &
      legen1, &
      legen2, &
      legen3, &
      legen4, &
      back_spec)

!SELF_SIMILAR_RAYLEIGH_GANS_PASSIVE  
!   Compute the scattering matrix, extinction matrix, and emission vector of an aggregate particle using the
!   Self-Similar Rayleigh-Gans Approximation
!
!   from here on original documentation (MM)
!
!   bcs = self_similar_rayleigh_gans(wavelength, dielectric_const, D, ...
!                                    volume, rg_kappa, rg_beta, rg_gamma)
!   where the arguments are:
!     bcs                  Backscatter cross-section (m2)
!     wavelength           Wavelength of radiation (m)
!     dielectric_const     Dielectric constant of material (complex)
!     D                    Particle size in direction of propagation (m)
!     volume               Volume of material (m3)
!     rg_kappa                Kurtosis parameter describing mean structure
!     rg_beta                 Prefactor of power-law describing fluctuations
!     rg_gamma                Exponent of power-law describing fluctuations
!
!   The input variables rg_beta, rg_gamma, rg_kappa and rg_zeta must be scalars, but any of the others
!   may be vectors or scalars, provided that any vectors are the same
!   length, and D is also a vector. The output backscatter cross section
!   will have the same size as D.
!
!   This function is particularly suitable for computing the millimetre-wave
!   backscatter cross-section of ice and snow aggregates. Note that the
!   input dielectric constant is the value for solid ice and the input
!   volume is the volume of solid ice within the particle.
!
!   Note that the output is sometimes referred to as radar cross-section and
!   has units of m2, rather than in some conventions where backscatter
!   cross-section has units m2 sr-1 and is a factor of 4pi smaller.
!
!   For further information on the meaning of the arguments, see:
!     Hogan, R. J., and C. D. Westbrook, 2014: Equation for the microwave
!     backscatter cross section of aggregate snowflakes using the
!     Self-Similar Rayleigh-Gans Approximation. J. Atmos. Sci., in press.
!   This function essentially implements Eq. 12.
!
!   This file is NOT the original written by Robin Hogan, but it is based to that
!   It implements SSRG in the natural way and uses rt3-style output to interface with
!   the parametrization needed for RT4

    implicit none

    real(kind=dbl), intent(in) :: freq  ! frequency [Hz]
    !real(kind=dbl), intent(in) :: t    ! temperature [K]
    integer, intent(in) :: liq_ice
    integer, intent(in) :: nbins
    real(kind=dbl), intent(in), dimension(nbins) :: diameter
    real(kind=dbl), intent(in), dimension(nbins) :: del_d    
    real(kind=dbl), intent(in), dimension(nbins) :: ndens
    real(kind=dbl), intent(in), dimension(nbins) :: mass ! in mie I have density
    real(kind=dbl), intent(in), dimension(nbins) :: as_ratio
    real(kind=dbl), intent(in), dimension(nbins) :: canting
    real(kind=dbl), intent(in) :: refre
    real(kind=dbl), intent(in) :: refim !positive(?)
    real(kind=dbl), intent(in) :: rg_kappa
    real(kind=dbl), intent(in) :: rg_beta
    real(kind=dbl), intent(in) :: rg_gamma
    real(kind=dbl), intent(in) :: rg_zeta

    real(kind=dbl), intent(out) :: extinction
    real(kind=dbl), intent(out) :: albedo
    real(kind=dbl), intent(out) :: back_scatt
    integer, parameter :: maxn  = 5000
    real(kind=dbl), intent(out), dimension(maxnleg) :: legen1, legen2, legen3, legen4
    real(kind=dbl), dimension(maxn) :: sump1, coef1, sump2, coef2,   &
                                       sump3, coef3, sump4, coef4
    complex(kind=dbl), dimension(maxn) :: a, b
    real(kind=dbl), intent(out), dimension(nbins) :: back_spec
    integer, intent(out) :: nlegen

    real(kind=dbl) :: wavelength
    real(kind=dbl) :: volume !volume of solid ice with same mass
    real(kind=dbl) :: K2 !dielectric variable square norm
    real(kind=dbl) :: wave_num !wavenumber
    real(kind=dbl) :: x !electrical size
    real(kind=dbl) :: term1, term2 ! partial terms in the ssrg formula
    real(kind=dbl) :: d_wave !size along beam propagation
    real(kind=dbl) :: prefactor, shape_fact, phas_func
    real(kind=dbl) :: pl, pl1, pl2
    real(kind=dbl) :: qext, qscat, qback, cabs, scatter, n_tot, del_d_eff, ndens_eff
    real(kind=dbl) :: sumqe, sumqs, sumqback, tmp
    real(kind=dbl) :: mu(maxn), wts(maxn)

    real(kind=dbl) :: angular_dep, scat_angle_rad
    !real(kind=dbl), dimension(2*nummu) :: angles_rad, phas_func

    integer(kind=long) :: ii, i, l, m
    integer(kind=long) :: jj
    integer(kind=long) :: jmax
    integer(kind=long) :: ia
    integer(kind=long) :: nterms, nleg, nquad
    !integer(kind=long) :: ia1
    !integer(kind=long) :: ia2

    complex(kind=dbl) :: dielectric_const, K
    complex(kind=dbl) :: s11,s22 ! elements of amplitude matrix

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'calc_self_similar_rayleigh_gans_rt3'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    err = 0

      call assert_true(err,(liq_ice==-1),&
          "only for ice, liq_ice must be -1")  
      call assert_true(err,all(mass>=0),&
          "mass must be positive")  
      call assert_true(err,all(ndens>=0),&
          "ndens must be positive")
      call assert_true(err,SUM(ndens)>=0,&
          "sum(ndens) must be greater equal zero")    
      call assert_true(err,all(diameter>0),&
          "diameter must be positive")   
      call assert_true(err,all(del_d>0),&
          "del_d must be positive")   
      call assert_true(err,(nbins>0),&
          "nbins must be positive")   
      call assert_true(err,(freq>0),&
          "freq must be positive")   
      call assert_true(err,all((canting == 0) .or. (canting == 90)),&
          "Not yet implemented: canting must be zero or 90deg")   
      call assert_true(err,all(as_ratio > 0.d0),&
          "nan or negative as_ratio")
!      call assert_false(err,active,&
!          "'active' must be turned off")
      if (err > 0) then
          errorstatus = fatal
          msg = "assertation error"
          call report(errorstatus, msg, nameOfRoutine)
          return
      end if    

    ! initialize intent(out) variables with 0
    extinction = 0.0d0
    albedo = 0.0d0
    back_scatt = 0.0d0
    legen1(:) = 0.0d0
    legen2(:) = 0.0d0
    legen3(:) = 0.0d0
    legen4(:) = 0.0d0
    back_spec(:) = 0.0d0
    nlegen = 0

    ! initialize partial inner results with 0
    sumqe = 0.0d0
    sumqs = 0.0d0
    sumqback = 0.0d0
    sump1 (:) = 0.0d0
    sump2 (:) = 0.0d0
    sump3 (:) = 0.0d0
    sump4 (:) = 0.0d0
    n_tot = 0.d0 ! total number density of particles

    wavelength = c/(freq) ! meters
    wave_num = 2.d0*pi/wavelength

    ! Complex dielectric factor
    dielectric_const = (refre+im*refim)**2 !solid ice ! TODO this sign + is suspicious!!!
    K = (dielectric_const-1.0d0)/(dielectric_const+2.0d0)
    K2 = abs(K)**2

    !!! build angles - here I copy from the Mie routine
    ! I am quite sure there is a generic method to derive the
    ! the number of angles to be used to sample scattering properties
    ! in order to get the Legendre series - Davide
    x = pi * diameter(nbins) / wavelength
    nterms = 0
    call miecalc (err,nterms, x, refre+im*refim, a, b)
    if (err /= 0) then
      msg = 'error in miecalc 4 ssrg!'
      call report(err, msg, nameOfRoutine)
      errorstatus = err
      return
    end if         
    nlegen = 2 * nterms 
    nlegen = min(maxnleg-1, nlegen) 
    nquad = (nlegen + 2 * nterms + 2) / 2 
    if (nquad.gt.maxn) then
      errorstatus = fatal
      msg = 'ssrg: maxn exceeded' 
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if
    call gausquad(nquad, mu, wts) ! got angles and weights    

    ! Begin loop over sizes
    sumqs = 0.0d0
    do ii = 1, nbins
      
      !Do not process if no particles present
      if (ndens(ii) <= 0) CYCLE

      ndens_eff = ndens(ii)
      del_d_eff = del_d(ii)
      n_tot = n_tot + (ndens_eff * del_d_eff) ! sum up particle number concentration

      if (canting(ii) == 0) then
        d_wave = diameter(ii) * as_ratio(ii)
      else if (canting(ii) == 90) then
        d_wave = diameter(ii)
      end if
      volume = mass(ii) / rho_ice       !volume of pure ice in particle
      prefactor = 9.d0*wave_num**4*K2*volume**2/(4.d0*pi)
      Cabs = 3.d0*wave_num*volume*dimag(K)
      !print*,"diam ",diameter(ii)," volume ",volume," freq ",freq," n",refre," j",refim
      !print*,"Cabs ",Cabs
      ! Integrate phase function over scattering angles to get Csca
      jmax=360
      qscat = 0.0d0
      do ia = 2, jmax ! every half degree, 0 and pi gives 0 contribution, but I need pi for back
        scat_angle_rad = dble(ia-1)*pi/dble(jmax-1)
        ! Electrical size
        x = wave_num * d_wave * sin(scat_angle_rad*0.5d0)
        call calc_shape_factor(x, rg_kappa, rg_beta, rg_gamma, rg_zeta, shape_fact)
        phas_func = prefactor*shape_fact*(1.d0 + cos(scat_angle_rad)**2)*0.5d0
        qscat = qscat + phas_func*sin(scat_angle_rad)
      end do
      qback = phas_func ! last phase function corresponds to backscattering
      qscat = 0.5*pi*qscat/dble(jmax) ! divide by domain after summation in the integral and coefficient 0.5

      !print*,"ssrg qscat ", qscat
      qext = qscat + Cabs
      !print*,"ssrg qback", qback

      !integrate=sum up . del_d is already included in ndens, since ndens is not normed!
      sumqe = sumqe + ( qext * ndens_eff * del_d_eff)
      sumqs = sumqs + ( qscat *ndens_eff * del_d_eff)
      sumqback = sumqback + ( qback * ndens_eff * del_d_eff)

      if (verbose >= 4) print*, "NEW: diameter(ii), ndens_eff, del_d_eff, n_tot, sumqback, sumqs, sumqe"
      if (verbose >= 4) print*, diameter(ii), ndens_eff, del_d_eff, n_tot, sumqback , sumqs, sumqe
      back_spec(ii) =  qback * ndens_eff  ! volumetric backscattering cross section for radar simulator in backscat per volume per del_d[m²/m⁴]


      do ia = 1, nquad
        scat_angle_rad = acos(mu(ia))
        x = wave_num * d_wave * sin(scat_angle_rad*0.5d0)
        call calc_shape_factor(x, rg_kappa, rg_beta, rg_gamma, rg_zeta, shape_fact)
        s22 = -im * 3. * wave_num**3 * K * volume * shape_fact**0.5d0 / (4.d0*pi) ! here I have multiplied by -j*wave_num because of mie_sphere convention
        s11 = s22*cos(scat_angle_rad)
        !print*,"s11 ",s11,"    s22 ",s22
        sump1(ia) = sump1(ia) + 0.5*(abs(s11)**2 + abs(s22)**2)*ndens_eff*del_d_eff
        sump2(ia) = sump2(ia) + 0.5*(abs(s11)**2 - abs(s22)**2)*ndens_eff*del_d_eff
        sump3(ia) = sump3(ia) + dreal(dconjg(s11)*s22)*ndens_eff*del_d_eff
        sump4(ia) = sump4(ia) + dimag(dconjg(s22)*s11)*ndens_eff*del_d_eff
        !print*,"sump1", sump1(ia)
      end do
    end do ! end loop sizes
    
!           multiply the sums by the integration delta and other constan
    !             put quadrature weights in angular array for later         

    if (verbose >= 4) print*, "ntot", n_tot        
        
    if (lhyd_absorption) then ! this might be the problem with negative extinction
        extinction = sumqe
    else
        extinction = sumqe - sumqs !remove scattering from extinction
    end if
    scatter = sumqs
    !print*,"ssrg scatter ", scatter 
    back_scatt = sumqback 
    !print*,"backscatt ", back_scatt
    albedo = scatter / extinction 
    !print*,"ssrg albedo ", albedo

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

    ! integrate the angular scattering functions times legendre   
    ! polynomials to find the legendre coefficients             
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
      !print*, m, coef1(m), coef2(m), coef3(m), coef4(m)
      legen1 (m) = (2 * l + 1) / 2.0 * coef1 (m) 
      legen2 (m) = (2 * l + 1) / 2.0 * coef2 (m) 
      legen3 (m) = (2 * l + 1) / 2.0 * coef3 (m) 
      legen4 (m) = (2 * l + 1) / 2.0 * coef4 (m) 
      if (legen1 (m) .gt. 1.0e-7) nlegen = l 
    end do
    !print*,legen1

    call assert_false(err,any(isnan(legen1)),&
        "nan in legen1")
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

  end subroutine calc_self_similar_rayleigh_gans_rt3

  subroutine calc_shape_factor(x, kappa, beta, gamma, zeta, shape)
    real(kind=dbl), intent(in) :: x, kappa, beta, gamma, zeta
    real(kind=dbl), intent(out) :: shape
    real(kind=dbl) :: scale, summ, xang
    integer :: jmax, jj

    !xang = x*sin(theta*0.5d0)

    ! Compute the first term in braces Eq.4 Hogan et al. (2017)
    scale = (cos(x)*((1.0d0+kappa/3.0d0)*(1.0d0/(2.0d0*x+pi)-1.0d0/(2.0d0*x-pi)) &
            - kappa*(1.0d0/(2.0d0*x+3.0d0*pi)-1.0d0/(2.0d0*x-3.0d0*pi))))**2
    
    ! Compute the summation component of the phi shape function for ssrga
    ! Decide how many terms are needed
    jmax = floor(5.d0*x/pi + 1.d0)
    ! Evaluate summation
    ! summ = 0.0d0 ! not needed anymore, initialized by the first term with rg_zeta
    summ = zeta*(2.d0)**(-gamma)*(1.d0/(2.d0*(x+pi))**2 + 1.d0/(2.d0*(x-pi))**2)
    do jj = 2, jmax
      summ = summ + (2.d0*jj)**(-gamma) &!* sin(x)**2 &
              *(1.d0/(2.d0*(x+pi*jj))**2 + 1.d0/(2.d0*(x-pi*jj))**2)
    end do
    summ = summ*beta*sin(x)**2
    shape = pi**2*0.25d0*(scale+summ)
  end subroutine calc_shape_factor


  subroutine calc_self_similar_rayleigh_gans_passive(&
      errorstatus, &
      freq, & ! frequency [Hz]
      liq_ice, &
      nbins, &
      diameter, &
      del_d, &
      ndens, & 
      mass, &
      as_ratio,&
      canting,&
      refre, &
      refim, &
      rg_kappa, &
      rg_beta, &
      rg_gamma, &
      scatter_matrix, &
      extinct_matrix, &
      emis_vector )

!SELF_SIMILAR_RAYLEIGH_GANS_PASSIVE  
!   Compute the scattering matrix, extinction matrix, and emission vector of an aggregate particle using the
!   Self-Similar Rayleigh-Gans Approximation
!
!   from here on original documentation (MM)
!
!   bcs = self_similar_rayleigh_gans(wavelength, dielectric_const, D, ...
!                                    volume, rg_kappa, rg_beta, rg_gamma)
!   where the arguments are:
!     bcs                  Backscatter cross-section (m2)
!     wavelength           Wavelength of radiation (m)
!     dielectric_const     Dielectric constant of material (complex)
!     D                    Particle size in direction of propagation (m)
!     volume               Volume of material (m3)
!     rg_kappa                Kurtosis parameter describing mean structure
!     rg_beta                 Prefactor of power-law describing fluctuations
!     rg_gamma                Exponent of power-law describing fluctuations
!
!   The input variables rg_beta and rg_gamma must be scalars, but any of the others
!   may be vectors or scalars, provided that any vectors are the same
!   length, and D is also a vector. The output backscatter cross section
!   will have the same size as D.
!
!   This function is particularly suitable for computing the millimetre-wave
!   backscatter cross-section of ice and snow aggregates. Note that the
!   input dielectric constant is the value for solid ice and the input
!   volume is the volume of solid ice within the particle.
!
!   Note that the output is sometimes referred to as radar cross-section and
!   has units of m2, rather than in some conventions where backscatter
!   cross-section has units m2 sr-1 and is a factor of 4pi smaller.
!
!   For further information on the meaning of the arguments, see:
!     Hogan, R. J., and C. D. Westbrook, 2014: Equation for the microwave
!     backscatter cross section of aggregate snowflakes using the
!     Self-Similar Rayleigh-Gans Approximation. J. Atmos. Sci., in press.
!   This function essentially implements Eq. 12.
!
!   This file was written by Robin Hogan, but no copyright is asserted: this
!   file is in the public domain.  Therefore, copying and distribution of
!   this file, with or without modification, and merging all or part of this
!   file into other works, are permitted in any medium without royalty.
!   This file is offered as-is, without any warranty.

    implicit none

    real(kind=dbl), intent(in) :: freq  ! frequency [Hz]
    integer, intent(in) :: liq_ice
    integer, intent(in) :: nbins
    real(kind=dbl), intent(in), dimension(nbins) :: diameter
    real(kind=dbl), intent(in), dimension(nbins) :: del_d    
    real(kind=dbl), intent(in), dimension(nbins) ::  ndens
    real(kind=dbl), intent(in), dimension(nbins) :: mass
    real(kind=dbl), intent(in), dimension(nbins) :: as_ratio
    real(kind=dbl), intent(in), dimension(nbins) :: canting
    real(kind=dbl), intent(in) :: refre
    real(kind=dbl), intent(in) :: refim !positive(?)
    real(kind=dbl), intent(in) :: rg_kappa
    real(kind=dbl), intent(in) :: rg_beta
    real(kind=dbl), intent(in) :: rg_gamma

    real(kind=dbl) :: wavelength
    real(kind=dbl) :: volume !volume of solid ice with same mass
    real(kind=dbl) :: K2 !dielectric variable square norm
    real(kind=dbl) :: wave_num !wavenumber
    real(kind=dbl) :: x !electrical size
    real(kind=dbl) :: term1, term2
    real(kind=dbl) :: d_wave !size along beam propagation

    real(kind=dbl) :: angular_dep, scat_angle_rad
    real(kind=dbl), dimension(2*nummu) :: angles_rad, phas_func

    integer(kind=long) :: ii
    integer(kind=long) :: jj
    integer(kind=long) :: jmax
    integer(kind=long) :: ia
    integer(kind=long) :: ia1
    integer(kind=long) :: ia2

    complex(kind=dbl) :: dielectric_const, K

    real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
    real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
    real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector

    real(kind=dbl), dimension(nstokes,2*nummu,nstokes,nummu) :: scatter_matrix_tmp
    real(kind=dbl), dimension(nstokes,nstokes,nummu) :: extinct_matrix_tmp
    real(kind=dbl), dimension(nstokes,nummu) :: emis_vector_tmp

    real(kind=dbl) :: scattering_integral_11, scattering_integral_12

    real(kind=dbl) :: fact_sca
    complex(kind=dbl) :: fact_ext
    complex(kind=dbl) :: s11,s22 ! elements of amplitude matrix

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'calc_self_similar_rayleigh_gans_passive'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    err = 0

      call assert_true(err,(liq_ice==-1),&
          "only for ice, liq_ice must be -1")  
      call assert_true(err,all(mass>=0),&
          "mass must be positive")  
      call assert_true(err,all(ndens>=0),&
          "ndens must be positive")
      call assert_true(err,SUM(ndens)>=0,&
          "sum(ndens) must be greater equal zero")    
      call assert_true(err,all(diameter>0),&
          "diameter must be positive")   
      call assert_true(err,all(del_d>0),&
          "del_d must be positive")   
      call assert_true(err,(nbins>0),&
          "nbins must be positive")   
      call assert_true(err,(freq>0),&
          "freq must be positive")   
      call assert_true(err,all((canting == 0) .or. (canting == 90)),&
          "Not yet implemented: canting must be zero or 90deg")   
      call assert_true(err,all(as_ratio > 0.d0),&
          "nan or negative as_ratio")
      call assert_false(err,active,&
          "'active' must be turned off")
      call assert_false(err,passive,&
          "This routine is not implemented yet. use ssrg-rt3 for active/passive calculations using self similar rayleigh gans") 
      if (err > 0) then
          errorstatus = fatal
          msg = "assertation error"
          call report(errorstatus, msg, nameOfRoutine)
          return
      end if    

    emis_vector = 0.d0
    extinct_matrix = 0.d0
    scatter_matrix = 0.d0
    scatter_matrix_tmp = 0.d0

    wavelength = c/(freq) !m
    wave_num = 2.d0*pi/wavelength

    ! build angles

    do ia = 1,nummu
      angles_rad(ia) = acos(mu_values(ia))
      angles_rad(ia+nummu) = acos(-1.*mu_values(ia))
    end do 
    !print*,refre,' +j ',refim
    dielectric_const = (refre+im*refim)**2 !solid ice
    K = (dielectric_const-1.0d0)/(dielectric_const+2.0d0)
    ! K2 = abs(((dielectric_const-1.0d0)/(dielectric_const+2.0d0)))**2
    fact_sca = 0.5e0/(wave_num**2)
    fact_ext = 2.0e0*pi*im/wave_num**2.e0 
    ! back_spec(:) = 0.d0
    ! sumqback = 0.d0
    
    do ii = 1, nbins

      if (canting(ii) == 0) then
        d_wave = diameter(ii) * as_ratio(ii)
      else if (canting(ii) == 90) then
        d_wave = diameter(ii)
      end if

      !volume of pure ice in particle
      volume = mass(ii) / rho_ice 

      ! Loop over scattering angles
      do ia1 = 1, nummu
        scattering_integral_11 = 0.
        scattering_integral_12 = 0.
        do ia2 = 1, 2*nummu
          ! scattering angle [rad]
          scat_angle_rad = abs(angles_rad(ia1) - angles_rad(ia2))

          ! Electrical size
          x = k*d_wave * sin(scat_angle_rad/2.)

          term1=(cos(x)*((1.0d0+rg_kappa/3.0d0)*(1.0d0/(2.0d0*x+pi)-1.0d0/(2.0d0*x-pi)) &
                  - rg_kappa*(1.0d0/(2.0d0*x+3.0d0*pi)-1.0d0/(2.0d0*x-3.0d0*pi))))**2

          ! Initialize the summation in the second term in the braces of Eq. 12
          term2 = 0.0d0

          ! Decide how many terms are needed
          jmax = floor(5.d0*x/pi + 1.d0)

          ! Evaluate summation
          do jj = 1, jmax
            term2 = term2 + (2.d0*jj)**(-rg_gamma) * sin(x)**2 &
                  *(1.d0/(2.d0*(x+pi*jj))**2 + 1.d0/(2.d0*(x-pi*jj))**2)
          end do
          s22 = 3./8. * wave_num**2 * K * volume * (term1 + rg_beta*term2)**0.5
          s22 = s22*wave_num
          s11 = s22*cos(scat_angle_rad) ! 
          !print*,scat_angle_rad,'  ',s11,' ssrg ',s22
          ! Put the terms together
          scatter_matrix_tmp(1,ia2,1,ia1) = (s11*dconjg(s11)+s22*dconjg(s22))*2.0d0*pi*fact_sca
          scatter_matrix_tmp(1,ia2,2,ia1) = (s11*dconjg(s11)-s22*dconjg(s22))*2.0d0*pi*fact_sca
          scatter_matrix_tmp(2,ia2,1,ia1) = (s11*dconjg(s11)-s22*dconjg(s22))*2.0d0*pi*fact_sca
          scatter_matrix_tmp(2,ia2,2,ia1) = (s11*dconjg(s11)+s22*dconjg(s22))*2.0d0*pi*fact_sca

          if (scat_angle_rad .eq. 0.) then
            extinct_matrix_tmp(1,1,ia1) = -real((s11 + s22)*fact_ext)
            extinct_matrix_tmp(1,2,ia1) = -real((s11 - s22)*fact_ext)
            extinct_matrix_tmp(2,1,ia1) = -real((s11 - s22)*fact_ext)
            extinct_matrix_tmp(2,2,ia1) = -real((s11 + s22)*fact_ext)
            !print*,extinct_matrix_tmp(1,1,ia1)
          end if
        ! End loop scattering angles
          scattering_integral_11 = scattering_integral_11 + scatter_matrix_tmp(1,ia2,1,ia1)*2.*pi*quad_weights(ia1)
          scattering_integral_12 = scattering_integral_12 + scatter_matrix_tmp(1,ia2,2,ia1)*2.*pi*quad_weights(ia1)
        end do
        emis_vector_tmp(1,ia1) = extinct_matrix_tmp(1,1,ia1) - scattering_integral_11
        emis_vector_tmp(2,ia1) = extinct_matrix_tmp(1,2,ia1) - scattering_integral_12
      end do

      !Apply psd

      scatter_matrix(1,:,1,:,1) = scatter_matrix(1,:,1,:,1) &
        + scatter_matrix_tmp(1,1:nummu,1,1:nummu)*ndens(ii) * del_d(ii)
      scatter_matrix(1,:,2,:,1) = scatter_matrix(1,:,2,:,1) &
        + scatter_matrix_tmp(1,1:nummu,2,1:nummu)*ndens(ii) * del_d(ii)
      scatter_matrix(2,:,1,:,1) = scatter_matrix(2,:,1,:,1) &
        + scatter_matrix_tmp(2,1:nummu,1,1:nummu)*ndens(ii) * del_d(ii)
      scatter_matrix(2,:,2,:,1) = scatter_matrix(2,:,2,:,1) &
        + scatter_matrix_tmp(2,1:nummu,2,1:nummu)*ndens(ii) * del_d(ii)
      scatter_matrix(1,:,1,:,2) = scatter_matrix(1,:,1,:,2) &
        + scatter_matrix_tmp(1,nummu+1:2*nummu,1,1:nummu)*ndens(ii) * del_d(ii)
      scatter_matrix(1,:,2,:,2) = scatter_matrix(1,:,2,:,2) &
        + scatter_matrix_tmp(1,nummu+1:2*nummu,2,1:nummu)*ndens(ii) * del_d(ii)
      scatter_matrix(2,:,1,:,2) = scatter_matrix(2,:,1,:,2) &
        + scatter_matrix_tmp(2,nummu+1:2*nummu,1,1:nummu)*ndens(ii) * del_d(ii)
      scatter_matrix(2,:,2,:,2) = scatter_matrix(2,:,2,:,2) &
        + scatter_matrix_tmp(2,nummu+1:2*nummu,2,1:nummu)*ndens(ii) * del_d(ii)

      extinct_matrix(:,:,:) = extinct_matrix(:,:,:) + extinct_matrix_tmp(:,:,:) * ndens(ii) * del_d(ii)
      emis_vector(:,:) = emis_vector(:,:) + emis_vector_tmp(:,:) * ndens(ii) * del_d(ii)
    end do


    errorstatus = err    
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return 

  end subroutine calc_self_similar_rayleigh_gans_passive

  ! subroutine ssrg_phas_function(mu)

  !       implicit none

  !       real(dbl), intent(in) :: mu

  !       angle_rad = acos(mu_values(ia))
  !       angular_dep = (1+mu_values(ia)**2)/2.
  !       ! Electrical size
  !       x = k*d_wave * sin(angle_rad/2.)

  !       term1=(cos(x)*((1.0d0+rg_kappa/3.0d0)*(1.0d0/(2.0d0*x+pi)-1.0d0/(2.0d0*x-pi)) &
  !               - rg_kappa*(1.0d0/(2.0d0*x+3.0d0*pi)-1.0d0/(2.0d0*x-3.0d0*pi))))**2

  !       ! Initialize the summation in the second term in the braces of Eq. 12
  !       term2 = 0.0d0

  !       ! Decide how many terms are needed
  !       jmax = floor(5.d0*x/pi + 1.d0)

  !       ! Evaluate summation
  !       do jj = 1, jmax
  !         term2 = term2 + (2.d0*jj)**(-rg_gamma) * sin(x)**2 &
  !               *(1.d0/(2.d0*(x+pi*jj))**2 + 1.d0/(2.d0*(x-pi*jj))**2)
  !       end do

  !       ! Put the terms together
  !       phas_func(ia) = prefactor*(term1 + rg_beta*term2) * angular_dep

  ! end subroutine ssrg_phas_function

  subroutine calc_self_similar_rayleigh_gans(&
      errorstatus, &
      freq, & ! frequency [Hz]
      liq_ice, &
      nbins, &
      diameter, &
      del_d, &
      ndens, & 
      mass, &
      as_ratio,&
      canting,&
      refre, &
      refim, &
      rg_kappa, &
      rg_beta, &
      rg_gamma, &
      back_spec, &
      sumqback )

!SELF_SIMILAR_RAYLEIGH_GANS  
!   Compute the backscatter cross section of an aggregate particle using the
!   Self-Similar Rayleigh-Gans Approximation
!
!   bcs = self_similar_rayleigh_gans(wavelength, dielectric_const, D, ...
!                                    volume, rg_kappa, rg_beta, rg_gamma)
!   where the arguments are:
!     bcs                  Backscatter cross-section (m2)
!     wavelength           Wavelength of radiation (m)
!     dielectric_const     Dielectric constant of material (complex)
!     D                    Particle size in direction of propagation (m)
!     volume               Volume of material (m3)
!     rg_kappa                Kurtosis parameter describing mean structure
!     rg_beta                 Prefactor of power-law describing fluctuations
!     rg_gamma                Exponent of power-law describing fluctuations
!
!   The input variables rg_beta and rg_gamma must be scalars, but any of the others
!   may be vectors or scalars, provided that any vectors are the same
!   length, and D is also a vector. The output backscatter cross section
!   will have the same size as D.
!
!   This function is particularly suitable for computing the millimetre-wave
!   backscatter cross-section of ice and snow aggregates. Note that the
!   input dielectric constant is the value for solid ice and the input
!   volume is the volume of solid ice within the particle.
!
!   Note that the output is sometimes referred to as radar cross-section and
!   has units of m2, rather than in some conventions where backscatter
!   cross-section has units m2 sr-1 and is a factor of 4pi smaller.
!
!   For further information on the meaning of the arguments, see:
!     Hogan, R. J., and C. D. Westbrook, 2014: Equation for the microwave
!     backscatter cross section of aggregate snowflakes using the
!     Self-Similar Rayleigh-Gans Approximation. J. Atmos. Sci., in press.
!   This function essentially implements Eq. 12.
!
!   This file was written by Robin Hogan, but no copyright is asserted: this
!   file is in the public domain.  Therefore, copying and distribution of
!   this file, with or without modification, and merging all or part of this
!   file into other works, are permitted in any medium without royalty.
!   This file is offered as-is, without any warranty.

    implicit none

    real(kind=dbl), intent(in) :: freq  ! frequency [Hz]
    integer, intent(in) :: liq_ice
    integer, intent(in) :: nbins
    real(kind=dbl), intent(in), dimension(nbins) :: diameter
    real(kind=dbl), intent(in), dimension(nbins) :: del_d    
    real(kind=dbl), intent(in), dimension(nbins) ::  ndens
    real(kind=dbl), intent(in), dimension(nbins) :: mass
    real(kind=dbl), intent(in), dimension(nbins) :: as_ratio
    real(kind=dbl), intent(in), dimension(nbins) :: canting
    real(kind=dbl), intent(in) :: refre
    real(kind=dbl), intent(in) :: refim !positive(?)
    real(kind=dbl), intent(in) :: rg_kappa
    real(kind=dbl), intent(in) :: rg_beta
    real(kind=dbl), intent(in) :: rg_gamma

    real(kind=dbl), intent(out), dimension(nbins) :: back_spec
    real(kind=dbl), intent(out) :: sumqback

    real(kind=dbl) :: wavelength
    real(kind=dbl) :: volume !volume of solid ice with same mass
    real(kind=dbl) :: K2 !dielectric variable
    real(kind=dbl) :: k !wavenumber
    real(kind=dbl) :: x !electrical size
    real(kind=dbl) :: prefactor
    real(kind=dbl) :: term1, term2
    real(kind=dbl) :: d_wave !size along beam propagation

    integer(kind=long) :: ii
    integer(kind=long) :: jj
    integer(kind=long) :: jmax

    complex(kind=dbl) :: dielectric_const

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'calc_self_similar_rayleigh_gans'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    err = 0

      call assert_true(err,(liq_ice==-1),&
          "only for ice, liq_ice must be -1")  
      call assert_true(err,all(mass>=0),&
          "mass must be positive")  
      call assert_true(err,all(ndens>=0),&
          "ndens must be positive")
      call assert_true(err,SUM(ndens)>=0,&
          "sum(ndens) must be greater equal zero")    
      call assert_true(err,all(diameter>0),&
          "diameter must be positive")   
      call assert_true(err,all(del_d>0),&
          "del_d must be positive")   
      call assert_true(err,(nbins>0),&
          "nbins must be positive")   
      call assert_true(err,(freq>0),&
          "freq must be positive")   
      call assert_true(err,all((canting == 0) .or. (canting == 90)),&
          "Not yet implemented: canting must be zero or 90deg")   
      call assert_true(err,all(as_ratio > 0.d0),&
          "nan or negative as_ratio")
      call assert_false(err,passive,&
          "'passive' must be turned off")   
      if (err > 0) then
          errorstatus = fatal
          msg = "assertation error"
          call report(errorstatus, msg, nameOfRoutine)
          return
      end if    

    wavelength = c/(freq) !m
    k = 2.d0*pi/wavelength

    dielectric_const = (refre-im*refim)**2 !solid ice
    K2 = abs(((dielectric_const-1.0d0)/(dielectric_const+2.0d0)))**2

    back_spec(:) = 0.d0
    sumqback = 0.d0
    do ii = 1, nbins

      if (canting(ii) == 0) then
        d_wave = diameter(ii) * as_ratio(ii)
      else if (canting(ii) == 90) then
        d_wave = diameter(ii)
      end if

      ! Electrical size
      x = k*d_wave
      !volume of pure ice in particle
      volume = mass(ii) / rho_ice 

      ! Factor outside the braces in Eq. 12 of Hogan and Westbrook
      prefactor = 9.0d0 * pi * k**4 *K2 * volume**2 /16

      term1=(cos(x)*((1.0d0+rg_kappa/3.0d0)*(1.0d0/(2.0d0*x+pi)-1.0d0/(2.0d0*x-pi)) &
              - rg_kappa*(1.0d0/(2.0d0*x+3.0d0*pi)-1.0d0/(2.0d0*x-3.0d0*pi))))**2

      ! Initialize the summation in the second term in the braces of Eq. 12
      term2 = 0.0d0
      ! Decide how many terms are needed
      jmax = floor(5.d0*x/pi + 1.d0)

      ! Evaluate summation
      do jj = 1, jmax
        term2 = term2 + (2.d0*jj)**(-rg_gamma) * sin(x)**2 &
              *(1.d0/(2.d0*(x+pi*jj))**2 + 1.d0/(2.d0*(x-pi*jj))**2)
      end do

      ! Put the terms together
      back_spec(ii) = prefactor*(term1 + rg_beta*term2)
      
      !Apply psd
      back_spec(ii) =  back_spec(ii) *ndens(ii)

      sumqback = sumqback + ( back_spec(ii) * del_d(ii))

    end do


    errorstatus = err    
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return 

  end subroutine calc_self_similar_rayleigh_gans

  subroutine calc_rayleigh_gans(&
      errorstatus, &
      freq, & ! frequency [Hz]
      liq_ice, &
      nbins, &
      diameter, &
      del_d, &
      ndens, & 
      density, &
      as_ratio,&
      canting,&
      refre, &
      refim, &
      back_spec, &
      sumqback )

!RAYLEIGH_GANS  
!   Compute the backscatter cross section of an aggregate particle using the
!   Self-Similar Rayleigh-Gans Approximation
!
!   bcs = self_similar_rayleigh_gans(wavelength, dielectric_const, D, ...
!                                    volume, rg_kappa, rg_beta, rg_gamma)
!   where the arguments are:
!     bcs                  Backscatter cross-section (m2)
!     wavelength           Wavelength of radiation (m)
!     dielectric_const     Dielectric constant of material (complex)
!     D                    Particle size in direction of propagation (m)
!     volume               Volume of material (m3)
!     rg_kappa                Kurtosis parameter describing mean structure
!     rg_beta                 Prefactor of power-law describing fluctuations
!     rg_gamma                Exponent of power-law describing fluctuations
!
!   The input variables rg_beta and rg_gamma must be scalars, but any of the others
!   may be vectors or scalars, provided that any vectors are the same
!   length, and D is also a vector. The output backscatter cross section
!   will have the same size as D.
!
!   This function is particularly suitable for computing the millimetre-wave
!   backscatter cross-section of ice and snow aggregates. Note that the
!   input dielectric constant is the value for solid ice and the input
!   volume is the volume of solid ice within the particle.
!
!   Note that the output is sometimes referred to as radar cross-section and
!   has units of m2, rather than in some conventions where backscatter
!   cross-section has units m2 sr-1 and is a factor of 4pi smaller.
!
!   For further information on the meaning of the arguments, see:
!     Hogan, R. J., and C. D. Westbrook, 2014: Equation for the microwave
!     backscatter cross section of aggregate snowflakes using the
!     Self-Similar Rayleigh-Gans Approximation. J. Atmos. Sci., in press.
!   This function essentially implements Eq. 12.
!
!   This file was written by Robin Hogan, but no copyright is asserted: this
!   file is in the public domain.  Therefore, copying and distribution of
!   this file, with or without modification, and merging all or part of this
!   file into other works, are permitted in any medium without royalty.
!   This file is offered as-is, without any warranty.

    implicit none

    real(kind=dbl), intent(in) :: freq  ! frequency [Hz]
    integer, intent(in) :: liq_ice
    integer, intent(in) :: nbins
    real(kind=dbl), intent(in), dimension(nbins) :: diameter
    real(kind=dbl), intent(in), dimension(nbins) :: del_d    
    real(kind=dbl), intent(in), dimension(nbins) ::  ndens
    real(kind=dbl), intent(in), dimension(nbins) :: density
    real(kind=dbl), intent(in), dimension(nbins) :: as_ratio
    real(kind=dbl), intent(in), dimension(nbins) :: canting
    real(kind=dbl), intent(in) :: refre
    real(kind=dbl), intent(in) :: refim !positive(?)

    real(kind=dbl), intent(out), dimension(nbins) :: back_spec
    real(kind=dbl), intent(out) :: sumqback

    real(kind=dbl) :: wavelength
    real(kind=dbl) :: K2 !dielectric variable
    real(kind=dbl) :: k !wavenumber
    real(kind=dbl) :: x !electrical size
    real(kind=dbl) :: prefactor
    real(kind=dbl) :: frac
    real(kind=dbl) :: d_wave !size along beam propagation

    integer(kind=long) :: ii

    complex(kind=dbl) :: dielectric_const

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'calc_rayleigh_gans'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    err = 0

      call assert_true(err,all(ndens>=0),&
          "ndens must be positive")
      call assert_true(err,SUM(ndens)>=0,&
          "sum(ndens) must be greater equal zero")    
      call assert_true(err,all(diameter>0),&
          "diameter must be positive")   
      call assert_true(err,all(del_d>0),&
          "del_d must be positive")   
      call assert_true(err,(nbins>0),&
          "nbins must be positive")   
      call assert_true(err,(freq>0),&
          "freq must be positive")   
      call assert_true(err,all((canting == 0) .or. (canting == 90)),&
          "canting must be zero or 90deg")   
      call assert_true(err,all(as_ratio > 0.d0),&
          "nan or negative as_ratio")
      call assert_false(err,passive,&
          "'passive' must be turned off")   
      if (err > 0) then
          errorstatus = fatal
          msg = "assertation error"
          call report(errorstatus, msg, nameOfRoutine)
          return
      end if    

    wavelength = c/(freq) !m
    k = 2.d0*pi/wavelength

    dielectric_const = (refre-im*refim)**2 !solid ice
    K2 = abs(((dielectric_const-1.0d0)/(dielectric_const+2.0d0)))**2

    back_spec(:) = 0.d0
    sumqback = 0.d0
    do ii = 1, nbins

      if (canting(ii) == 0) then
        d_wave = diameter(ii) * as_ratio(ii)
      else if (canting(ii) == 90) then
        d_wave = diameter(ii)
      end if

      ! Electrical size
      x = k*d_wave
      
      if (liq_ice == -1) then
        frac = density(ii) / rho_ice     
      else if (liq_ice == 1) then
        frac = 1.d0
      else
        errorstatus = fatal
        msg = "liq_ice must be 1 or -1"
        call report(errorstatus, msg, nameOfRoutine)
        return
      end if



      ! Factor outside the braces in Eq. 1 of Hogan and Westbrook
      prefactor = 9.0d0 * pi *K2 * frac**2 / (16.d0 * k**2 * as_ratio(ii)**4)
      
      back_spec(ii) = prefactor * (sin(x) - x * cos(x))**2
      
      !Apply psd
      back_spec(ii) =  back_spec(ii) *ndens(ii)

      sumqback = sumqback + ( back_spec(ii) * del_d(ii))

    end do

    errorstatus = err    
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return 

  end subroutine calc_rayleigh_gans


  subroutine calc_rayleigh(&
      errorstatus, &
      freq, & ! frequency [Hz]
      liq_ice, &
      nbins, &
      diameter, &
      del_d, &
      ndens, & 
      density, &
      refre, &
      refim, &
      back_spec, &
      sumqback )

!pure RAYLEIGH

use settings, only: radar_K2

    implicit none

    real(kind=dbl), intent(in) :: freq  ! frequency [Hz]
    integer, intent(in) :: liq_ice
    integer, intent(in) :: nbins
    real(kind=dbl), intent(in), dimension(nbins) :: diameter
    real(kind=dbl), intent(in), dimension(nbins) :: del_d    
    real(kind=dbl), intent(in), dimension(nbins) ::  ndens
    real(kind=dbl), intent(in), dimension(nbins) :: density
    real(kind=dbl), intent(in) :: refre
    real(kind=dbl), intent(in) :: refim !positive(?)

    real(kind=dbl), intent(out), dimension(nbins) :: back_spec
    real(kind=dbl), intent(out) :: sumqback

    real(kind=dbl) :: wavelength
    real(kind=dbl) :: K2 !dielectric variable
    real(kind=dbl) :: prefactor

    integer(kind=long) :: ii

    complex(kind=dbl) :: dielectric_const
    complex(kind=dbl) :: eps_mix, m_ice

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'calc_rayleigh_gans'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    err = 0

      call assert_true(err,all(ndens>=0),&
          "ndens must be positive")
      call assert_true(err,SUM(ndens)>=0,&
          "sum(ndens) must be greater equal zero")    
      call assert_true(err,all(diameter>0),&
          "diameter must be positive")   
      call assert_true(err,all(del_d>0),&
          "del_d must be positive")   
      call assert_true(err,(nbins>0),&
          "nbins must be positive")   
      call assert_true(err,(freq>0),&
          "freq must be positive")   
      call assert_false(err,passive,&
          "'passive' must be turned off")   
      if (err > 0) then
          errorstatus = fatal
          msg = "assertation error"
          call report(errorstatus, msg, nameOfRoutine)
          return
      end if    

    wavelength = c/(freq) !m

    if ((liq_ice == -1) .and. (density(1) /= 917.d0)) then
        m_ice = refre-Im*refim  ! mimicking a
        dielectric_const = eps_mix((1.d0,0.d0),m_ice,density(nbins))
        dielectric_const = dielectric_const**2
    else
              dielectric_const = (refre-im*refim)**2
    end if

    K2 = abs(((dielectric_const-1.0d0)/(dielectric_const+2.0d0)))**2
    back_spec(:) = 0.d0
    sumqback = 0.d0
    do ii = 1, nbins

      prefactor = pi**5 * K2 / wavelength**4
      
      back_spec(ii) = prefactor * diameter(ii)**6
 
     
      !Apply psd
      back_spec(ii) =  back_spec(ii) *ndens(ii)

      sumqback = sumqback + ( back_spec(ii) * del_d(ii))

    end do


    errorstatus = err    
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return 

  end subroutine calc_rayleigh

end module rayleigh_gans



