module rayleigh_gans

  use kinds
  use constants, only: pi,c, Im, rho_ice
  use settings, only: passive
  use report_module
  use mie_scat_utilities  
  use vars_index, only: i_x,i_y, i_z, i_h
  implicit none


  contains

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
    real(kind=dbl) :: d_eq !mass equivalent diameter
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

      d_eq = (mass(ii) * 6.d0 / (rho_ice * pi))**(1.d0/3.d0) !better estimate it idrectly from density...
      volume = d_eq**3 * pi /6.d0

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



