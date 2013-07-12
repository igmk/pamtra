
module dia2vel
    !
    ! Description:
    !collection of subroutines to convert hydrometeor diameter to fall velocity
    !
    ! Current Code Owner: IGMK
    !
    ! History:
    !
    ! Version   Date       Comment
    ! -------   ----       -------
    ! 0.1       28/11/2012 created M.Maahn
    ! 0.2       15/04/2013 Application of European Standards for Writing and
    !                      Documenting Exchangeable Fortran 90 Code - M. Maahn

    ! Code Description:
    !   Language:		Fortran 90.
    !   Software Standards: "European Standards for Writing and  
    !     Documenting Exchangeable Fortran 90 Code". 
    !


  contains

  subroutine dia2vel_heymsfield10_particles_ms_as &
    (errorstatus,&      ! out
    nDia,&              !in
    diaSpec_SI,&        !in
    rho_air_SI,&        !in
    my_SI,&             !in
    mass_size_a_SI,&    !in
    mass_size_b,&       !in
    area_size_a_SI,&    !in
    area_size_b,&       !in
    velSpec)            !out
  !to be deleted after transition to pamtra v1.0
      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none 
      
      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec_SI
      real(kind=dbl), intent(in) :: rho_air_SI, my_SI,mass_size_a_SI,mass_size_b,&
      area_size_a_SI,area_size_b
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
      real(kind=dbl), dimension(ndia)::mass,area
  
      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=63) :: nameOfRoutine = 'dia2vel_heymsfield10_particles_ms_as'  
      
      mass = mass_size_a_SI * diaSpec_SI**mass_size_b
      area = area_size_a_SI * diaSpec_SI**area_size_b
  
      call dia2vel_heymsfield10_particles &
        (err,&      ! out
        nDia,&              !in
        diaSpec_SI,&        !in
        rho_air_SI,&        !in
        my_SI,&             !in
        mass,&              !in
        area,&              !in
        velSpec)            !out
  
  errorstatus = err
  
  end subroutine  dia2vel_heymsfield10_particles_ms_as
  
  subroutine dia2vel_heymsfield10_particles &
    (errorstatus,&     	! out
    nDia,& 		!in
    diaSpec_SI,&	!in
    rho_air_SI,&	!in
    my_SI,&		!in
    mass,&              !in
    area,&              !in
    velSpec)		!out

      !in
      !nDia: no of diameters
      !diaSpec_SI = diameter spectrum [m]
      !rho_air_SI = density of air [kg/m³]
      !my_SI = kinematic viscosity of air [m²/s]
      !mass = mass of the particle [SI]
      !area = cross section area [SI]
      !out
      !velSpec: velocity spectrum [m/s]

      !Heymsfield, A. J. & Westbrook, C. D. Advances in the Estimation of Ice Particle Fall Speeds
      !Using Laboratory and Field Measurements. Journal of the Atmospheric Sciences 67, 2469–2482 (2010).

      ! equations 9-11


      ! Modules used:
      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none

      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec_SI, mass,area
      real(kind=dbl), intent(in) :: rho_air_SI, my_SI
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
      real(kind=dbl), dimension(ndia)::area_proj, Xstar, Re
      real(kind=dbl) ::k, delta_0, C_0, eta


      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=33) :: nameOfRoutine = 'dia2vel_heymsfield10_particles'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

      ! no check for baundaries due to different possible particle types
      err = success

      k = 0.5d0 !defined in the paper
      delta_0 = 8.0d0
      C_0 = 0.35d0

      area_proj = area/((pi/4.d0)*diaSpec_SI**2)
      eta = my_SI * rho_air_SI !now dynamic viscosity

      Xstar = 8.d0 * rho_air_SI * mass * g / (pi * area_proj**(1.d0-k) * eta**2) !eq 9
      Re = delta_0**2/4.d0 * ((1.d0+((4.d0*sqrt(Xstar))/(delta_0**2*sqrt(C_0))))**0.5 - 1 )**2 !eq10
     
      velSpec = eta * Re /(rho_air_SI*diaSpec_SI)

      if (verbose .gt. 1) print*, 'Exiting dia2vel_khvorostyanov01_particles in dia2vel.f90'

      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

      return

  end subroutine dia2vel_heymsfield10_particles


  subroutine dia2vel_khvorostyanov01_particles &
    (errorstatus,&     	! out
    nDia,& 		!in
    diaSpec_SI,&	!in
    rho_air_SI,&	!in
    my_SI,&		!in
    mass_size_a_SI,&	!in
    mass_size_b,&	!in
    area_size_a_SI,&	!in
    area_size_b,&	!in
    velSpec)		!out

      !in
      !nDia: no of diameters
      !diaSpec_SI = diameter spectrum [m]
      !rho_air_SI = density of air [kg/m³]
      !my_SI = kinematic viscosity of air [m²/s]
      !mass_size_a_SI,mass_size_b parameters of mass size relation m = a*D_max^b [SI]
      !area_size_a_SI,area_size_b parameters of mass area relation A = a*D_max^b [SI]
      !out
      !velSpec: velocity spectrum [m/s]

      ! Khvorostyanov, V. I. & Curry, J. A. Terminal Velocities of Droplets and Crystals:
      ! Power Laws with Continuous Parameters over the Size Spectrum.
      ! Journal of the Atmospheric Sciences 59, 1872–1884 (2002).
      ! equation 3.3


      ! Modules used:
      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none

      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec_SI
      real(kind=dbl), intent(in) :: rho_air_SI, my_SI,mass_size_a_SI,mass_size_b,&
      area_size_a_SI,area_size_b
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
      real(kind=dbl) :: rho_air, g_cp, my
      real(kind=dbl) :: c1, delta0, mass_size_a, area_size_a
      real(kind=dbl), dimension(ndia):: X, bRe, aRe, Av, Bv
      real(kind=dbl), dimension(ndia)::diaSpec

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=33) :: nameOfRoutine = 'dia2vel_khvorostyanov01_particles'


      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

      ! no check for baundaries due to different possible particle types
      err = success

      !variables to cgs...
      mass_size_a = mass_size_a_SI * 10d0**(3.d0 - 2.d0*mass_size_b)
      area_size_a = area_size_a_SI * 10d0**(4.d0 - 2.d0*area_size_b)

      g_cp = g*100.d0 !cm/s
      rho_air = rho_air_SI/1000.d0 !g/cm³
      diaSpec = 100.d0*diaSpec_SI !cm
      my = my_SI*10000.d0 !cm² /s
      c1 = 0.0902d0
      delta0 = 9.06d0

      X = 2*mass_size_a*g_cp*diaSpec**(mass_size_b + 2.d0 - area_size_b)/&
      (area_size_a*rho_air*my**2) !eq. 17 of Mitchell et al 1996

      bRe = 0.5d0*c1*X**0.5d0*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**(-1) * &
      (1 + c1*X**0.5d0)**(-0.5d0) !2.12

      aRe = (delta0**2/4.d0)*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**2 / X**bRe !2.13


      Av = aRe * my**(1.d0-2.d0*bRe) * ((2.d0*mass_size_a*g_cp)/(rho_air*area_size_a))**bRe
      Bv = (bRe * (mass_size_b-area_size_b+2.d0)) - 1.d0

      velSpec = Av * diaSpec**Bv

      velSpec = velSpec/100.d0 !CGS to SI, now m/s

      if (verbose .gt. 1) print*, 'Exiting dia2vel_khvorostyanov01_particles in dia2vel.f90'

      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

      return

  end subroutine dia2vel_khvorostyanov01_particles


  subroutine dia2vel_khvorostyanov01_spheres &
    (errorstatus,&     	! out
    nDia,&     		! in
    diaSpec,&     	! in
    rho_air,&     	! in
    my,&     		! in
    rho_particle,&     	! in
    velSpec)		! out
    
      !in
      !nDia: no of diameters
      !diaSpec = diameter spectrum [m]
      !rho_air = density of air [kg/m³]
      !rho_particle = density of particle [kg/m³]
      !my = kinematic viscosity of air [m²/s]

      !
      !out
      !velSpec: velocity spectrum [m/s]

      ! Khvorostyanov, V. I. & Curry, J. A. Terminal Velocities of Droplets and Crystals:
      ! Power Laws with Continuous Parameters over the Size Spectrum.
      ! Journal of the Atmospheric Sciences 59, 1872–1884 (2002).
      ! equation 3.3


      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none

      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec, rho_particle
      real(kind=dbl), intent(in) :: rho_air, my
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
      real(kind=dbl) :: rho_air_cp, g_cp, my_cp
      real(kind=dbl) :: c1, delta0
      real(kind=dbl), dimension(ndia):: vB, A, X, bRe, aRe
      real(kind=dbl), dimension(ndia)::diaSpec_cp, rho_particle_cp

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=31) :: nameOfRoutine = 'dia2vel_khvorostyanov01_particles'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

      ! no check for baundaries due to different possible particle types
      err = success

      !variables to cgs...
      rho_particle_cp = rho_particle/1000.d0 !g/cm³
      g_cp = g*100.d0 !cm/s
      rho_air_cp = rho_air/1000.d0 !g/cm³
      diaSpec_cp = 100.d0*diaSpec !cm
      my_cp = my*10000.d0 !cm² /s
      c1 = 0.0902d0
      delta0 = 9.06d0

      vB = 4d0/3d0 * pi * (diaSpec_cp*0.5d0)**3
      A = pi * (diaSpec_cp*0.5d0)**2
      X = 2*vB*(rho_particle_cp-rho_air_cp)*g_cp*diaSpec_cp**2/&
      (A*rho_air_cp*my_cp**2) !eq 2.7

      bRe = 0.5d0*c1*X**0.5d0*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**(-1) * &
      (1 + c1*X**0.5d0)**(-0.5d0) !eq. 2.12

      aRe = (delta0**2/4.d0)*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**2 / X**bRe !2.13

      velSpec = aRe*my_cp**(1-2*bRe) * (2*vB*g_cp/A * ((rho_particle_cp-rho_air_cp)/rho_air_cp))**bRe *&
      diaSpec_cp**(2*bRe-1)             !2.20b

      velSpec = velSpec/100.d0 !CGS to SI, now m/s

      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

      return
  end subroutine dia2vel_khvorostyanov01_spheres


  subroutine dia2vel_khvorostyanov01_drops &
    (errorstatus,&     	! out
    nDia,&    	 	! in
    diaSpec,&     	! in
    rho_air,&     	! in
    my,&     		! in
    velSpec)		! out
    
      !in
      !nDia: no of diameters
      !diaSpec = diameter spectrum [m]
      !rho_air density of air [kg/m³]
      !my = kinematic viscosity of air [m²/s]

      !out
      !velSpec: velocity spectrum [m/s]

      ! Khvorostyanov, V. I. & Curry, J. A. Terminal Velocities of Droplets and Crystals:
      ! Power Laws with Continuous Parameters over the Size Spectrum.
      ! Journal of the Atmospheric Sciences 59, 1872–1884 (2002).
      ! equation 2.20c

      ! defined for diameters from 0 to 8.5 mm, non spherical shape is corrected

      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none

      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), intent(in) :: rho_air, my
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
      real(kind=dbl) :: rho_air_cp, g_cp, rho_water_cp, my_cp
      real(kind=dbl) :: c1, delta0, lam
      real(kind=dbl), dimension(ndia):: vB, A, X, bRe, aRe, xi
      real(kind=dbl), dimension(ndia)::diaSpec_cp

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=29) :: nameOfRoutine = 'dia2vel_khvorostyanov01_drops'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

      !check for boundaries (including 1% numerical tolerance)
      if (MAXVAL(diaSpec) > 8.5d-3*1.01d0) then
	  print*, "Largest diameter for dia2vel_khvorostyanov01_spheres is 8.5 mm, got",&
	  MAXVAL(diaSpec), "[m]"
	  errorstatus = fatal
	  msg = 'Diameter out of specs!'
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      else
	  err = success
      end if
      !variables to cgs...
      rho_water_cp = rho_water/1000.d0 !g/cm³
      g_cp = g*100.d0 !cm/s
      rho_air_cp = rho_air/1000.d0 !g/cm³
      diaSpec_cp = 100.d0*diaSpec !cm
      my_cp = my*10000.d0 !cm² /s
      c1 = 0.0902d0
      delta0 = 9.06d0


      lam = 0.47
      xi = exp(-diaSpec_cp/lam)+(1.d0-exp(-diaSpec_cp/lam))*(1.d0/(1.d0+(diaSpec_cp/lam))) !eq. 3.4

      vB = 4d0/3d0 * pi * (diaSpec_cp*0.5d0)**3 * xi !correct for non-spherity
      A = pi * (diaSpec_cp*0.5d0)**2 !no correction neccessary for A

      X = 2*vB*(rho_water_cp-rho_air_cp)*g_cp*diaSpec_cp**2/&
      (A*rho_air_cp*my_cp**2) !eq. 2.7

      bRe = 0.5d0*c1*X**0.5d0*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**(-1) * &
      (1 + c1*X**0.5d0)**(-0.5d0) !eq. 2.12

      aRe = (delta0**2/4.d0)*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**2 / X**bRe !eq. 2.13



      velSpec = aRe*my_cp**(1-2*bRe) * (4.d0/3.d0 * g_cp * xi * ((rho_water_cp-rho_air_cp)/rho_air_cp))**bRe *&
      diaSpec_cp**(3*bRe-1)  !eq. 2.20c

      velSpec = velSpec/100.d0 !CGS to SI, now m/s

      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

      return

  end subroutine dia2vel_khvorostyanov01_drops




  subroutine dia2vel_foote69_rain &
    (errorstatus,&     	! out
    nDia,&    	 	! in
    diaSpec,&     	! in
    rho_air,&     	! in
    temp,&    	 	! in
    velSpec)		! out
      

      !in
      !nDia: no of diameters
      !diaSpec = diameter spectrum [m]
      !rho_air density of air [kg/m³]
      !out
      !velSpec: velocity spectrum [m/s]

      !Foote, G. B. & Du Toit, P. S. Terminal Velocity of Raindrops Aloft. Journal of Applied Meteorology 8, 249–253 (1969).
      ! defined for 0.1-5.8 mm diameter

      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none

      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), intent(in) :: rho_air, temp
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
      real(kind=dbl) :: aj(10), Y, rho0
      integer :: jj

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=20) :: nameOfRoutine = 'dia2vel_foote69_rain'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

      !check for boundaries (including 1% numerical tolerance)
      if (MINVAL(diaSpec) < 1d-4/1.01d0) then
	  print*, "Smallest diameter for dia2vel_foote69_rain is 0.1mm, got",&
	  MINVAL(diaSpec), "[m]"
	  errorstatus = fatal
	  msg = 'Diameter out of specs!'
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      else if (MAXVAL(diaSpec) > 5.8d-3*1.01d0) then
	  print*, "Largest diameter for dia2vel_foote69_rain is 5.8mm, got",&
	  MAXVAL(diaSpec), "[m]"
	  errorstatus = fatal
	  msg = 'Diameter out of specs!'
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      else
	  err = success
      end if
      aj(:) = (/-8.5731540d-2, 3.3265862d0, 4.3843578d0, -6.8813414d0, 4.7570205d0,&
      -1.9046601d0, 4.6339978d-1, -6.7607898d-2, 5.4455480d-3, -1.8631087d-4/)

      velSpec = 0.d0
      do jj=1,10
	  velSpec = velSpec + (aj(jj)*(diaSpec*1d3)**(jj-1))
      end do


      rho0 = 1.2038631624242195d0 !p=1013 hPa, T=20°C
      if (rho0 < rho_air) rho0 = rho_air*1.00001 !it's numeric, stupid!

      !apply density correction
      Y = 0.43d0*log10(rho0/rho_air)-0.4d0*(log10(rho0/rho_air))**2.5d0
      velSpec = velSpec*10.d0**Y*(1d0+(0.0023d0*(1.1-(rho_air/rho0))*(293.15d0-temp)))

      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

      return

  end subroutine dia2vel_foote69_rain


  subroutine dia2vel_pavlos_cloud &
    (errorstatus,&     	! out
    nDia,&     		! in
    diaSpec,&     	! in
    velSpec)		! out
    
      !in
      !nDia: no of diameters
      !diaSpec = diameter spectrum [m]
      !out
      !velSpec: velocity spectrum [m/s]

      !from Pavlos Radar Simulator

      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none

      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=20) :: nameOfRoutine = 'dia2vel_pavlos_cloud'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
      err = success

      velSpec = 30.d0*(diaSpec*1.d3)**2 ! [m/s], diameter_spec_cp in m!

      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)


      return

  end subroutine dia2vel_pavlos_cloud


  subroutine dia2vel_metek_rain &
    (errorstatus,&     	! out
    nDia,&    	 	! in
    diaSpec,&     	! in
    rho_air,&     	! in
    temp,&     		! in
    velSpec)		! out
    
      !in
      !nDia: no of diameters
      !diaSpec = diameter spectrum [m]
      !out
      !velSpec: velocity spectrum [m/s]

      !v from Metek physical basics, density correction by foot et al
      ! for 0.109 ≤ D ≤ 6 mm


      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none



      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), intent(in) :: rho_air, temp
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
      real(kind=dbl) :: rho0, Y

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=20) :: nameOfRoutine = 'dia2vel_metek_rain'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

      !check for boundaries (including 1% numerical tolerance)
      if (MINVAL(diaSpec) < 1.09d-4/1.01d0) then
	  print*, "Smallest diameter for dia2vel_metek_rain is 0.109 mm, got", &
	  MINVAL(diaSpec), "[m]"
	  errorstatus = fatal
	  msg = 'Diameter out of specs!'
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      else if (MAXVAL(diaSpec) > 6d-3*1.01d0) then
	  print*, "Largest diameter for dia2vel_metek_rain is 6 mm, got", &
	  MAXVAL(diaSpec), "[m]"
	  errorstatus = fatal
	  msg = 'Diameter out of specs!'
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      else
	  err = success
      end if
      velSpec= ( 9.65d0 - 10.3d0 * exp(-0.6d0 *diaSpec*1.d3))

      rho0 = 1.2038631624242195d0 !p=1013 hPa, T=20°C
      if (rho0 < rho_air) rho0 = rho_air*1.00001 !it's numeric, stupid!

      !apply density correction
      Y = 0.43d0*log10(rho0/rho_air)-0.4d0*(log10(rho0/rho_air))**2.5d0
      velSpec = velSpec*10.d0**Y*(1d0+(0.0023d0*(1.1-(rho_air/rho0))*(293.15d0-temp)))

      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

      return

  end subroutine dia2vel_metek_rain


  subroutine dia2vel_rogers_drops &
    (errorstatus,&     	! out
    nDia,&   	 	! in
    diaSpec,&     	! in
    rho_air,&     	! in
    velSpec)		! out
    
      !in
      !nDia: no of diameters
      !diaSpec = diameter spectrum [m]
      !rho_air density of air [kg/m³]
      !out
      !velSpec: velocity spectrum [m/s]

      ! Rogers, R. R. & Yau, M. K. A short course in cloud physics. (Pergamon Press: 1989).
      ! D from 0 to 8.5 mm (breackup of drops)
      ! this function is NOT homegenous!

      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none

      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), intent(in) :: rho_air
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
      real(kind=dbl) :: rho0

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=18) :: nameOfRoutine = 'dia2vel_rogers_drops'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

      !check for boundaries (including 1% numerical tolerance)
      if (MAXVAL(diaSpec) > 8.5d-3*1.01d0) then
	  print*, "Largest diameter for dia2vel_rogers_drops is 8.5 mm, got", &
	  MAXVAL(diaSpec), "[m]"
	  errorstatus = fatal
	  msg = 'Diameter out of specs!'
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      else
	  err = success
      end if


      rho0 = 1.2038631624242195d0 !p=1013 hPa, T=20°C

      where (diaSpec < 80d-6) !small drops
	  velSpec = 1.19d6*(diaSpec*0.5*1d2)**2
      elsewhere (diaSpec > 1200d-6) !large drops
      velSpec = 2.2d3 * (rho0/rho_air)**0.5d0 * (diaSpec*0.5*1d2)**0.5d0 ! [m/s], diaSpec in m!
      elsewhere !intermediate
	  velSpec = 8d3*diaSpec*0.5*1d2
      end where

      velSpec = velSpec*1d-2 !CGS to Si

      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

    return
  end subroutine dia2vel_rogers_drops

  subroutine dia2vel_rogers_graupel &
    (errorstatus,&     	! out
    nDia,&     		! in
    diaSpec,&     	! in
    velSpec)		! out
    
      !in
      !nDia: no of diameters
      !diaSpec = diameter spectrum [m] in Rogers definition more dMin than dMax, which is diaSpec
      !out
      !velSpec: velocity spectrum [m/s]

      ! Rogers, R. R. & Yau, M. K. A short course in cloud physics. (Pergamon Press: 1989). eq. 9.5
      ! attention, this formula is without any density correction!
      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none

      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=20) :: nameOfRoutine = 'dia2vel_rogers_graupel'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
      err = success

      velSpec = 1.d-3 * 343.d0*(diaSpec*1.d3)**0.6d0 ! [m/s], diameter_spec_cp in m!
      
      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

      return
  end subroutine dia2vel_rogers_graupel
end module
