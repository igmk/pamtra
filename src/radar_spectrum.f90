subroutine radar_spectrum(nbins,diameter_spec, back,back_spec,&
  temp,press,hgt,frequency,particle_type,mass_size_a, mass_size_b,area_size_a,&
  area_size_b,particle_spec)
! this routine takes the backscattering spectrum depending on size and converts it
! into a spectrum depending on radar Doppler (=fall) velocity
! based on Spectra_simulator by P. Kollias
! converted from Matlab to Fortran by M. Maahn, IGMK (11/2012)
!in
!nbins: No of bins + 1
!diameter_spec: Diameter Spectrum (SI)
!back_spec: backscattering cross section per volume in m²/m⁴ (includes number density)
!temp: temperature in K
!press: air pressure in Pa
!frequency in GHz
!particle_type: cloud|rain|snow|ice|hail|graupel
!mass_size_a: a of mass size relation, needed for graupel, hail, snow, ice
!mass_size_b: b of mass size relation, needed for graupel, hail, snow, ice
!area_size_a: a of mass size relation, needed for graupel, hail, snow, ice
!area_size_b: b of mass size relation, needed for graupel, hail, snow, ice
!out
!particle_spec particle spectrum in dependence of radar Dopple velocity in m6m-3/ms-1

  use kinds
  use nml_params
  use constants
  implicit none

  integer,intent(in) ::  nbins !for Mie it is actually nbins+1
  real(kind=dbl), dimension(nbins),intent(in):: diameter_spec, back_spec
  character(5),intent(in) ::  particle_type
  real(kind=dbl), intent(in):: temp, frequency, press,hgt, back, mass_size_a, mass_size_b,&
    area_size_a, area_size_b

  real(kind=dbl), dimension(nbins):: vel_spec,dD_dU,back_vel_spec, back_spec_ref,&
				del_v_model, diameter_spec_cp
  real(kind=dbl), dimension(nbins) :: vel_spec_ext, back_vel_spec_ext
  real(kind=dbl), dimension(:,:), allocatable :: particle_spec_ext
  real(kind=dbl), intent(out), dimension(radar_nfft_aliased):: particle_spec
  real(kind=dbl), dimension(radar_nfft_aliased):: radar_velo_aliased, rho_particle
  real(kind=dbl):: del_v_radar, K2, dielec_water, wavelength, &
	delta_air, rho_air, rho, viscosity_air, my, del_d, Ze, K, &
	min_V_aliased, max_V_aliased
  integer :: ii, jj

  interface
    subroutine dia2vel_khvorostyanov01_particles(nDia,diaSpec_SI,rho_air_SI,my_SI,&
		  mass_size_a_SI,mass_size_b,area_size_a_SI,area_size_b,velSpec)
      use kinds
      implicit none
      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec_SI
      real(kind=dbl), intent(in) :: rho_air_SI, my_SI,mass_size_a_SI,mass_size_b,&
	      area_size_a_SI,area_size_b
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec      
    end subroutine dia2vel_khvorostyanov01_particles

    subroutine dia2vel_khvorostyanov01_spheres(nDia,diaSpec,rho_air,my,rho_particle,velSpec)
      use kinds
      implicit none
      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec, rho_particle
      real(kind=dbl), intent(in) :: rho_air, my
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec      
    end subroutine dia2vel_khvorostyanov01_spheres

    subroutine dia2vel_khvorostyanov01_drops(nDia,diaSpec,rho_air,my,velSpec)
      use kinds
      implicit none
      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), intent(in) :: rho_air, my
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec      
    end subroutine dia2vel_khvorostyanov01_drops

    subroutine dia2vel_foote69_rain(nDia,diaSpec,rho_air,temp,velSpec)
      use kinds
      implicit none
      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), intent(in) :: rho_air, temp
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec      
    end subroutine dia2vel_foote69_rain

    subroutine dia2vel_pavlos_cloud(nDia,diaSpec,velSpec)
      use kinds
      implicit none
      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec      
    end subroutine dia2vel_pavlos_cloud

    subroutine dia2vel_metek_rain(nDia,diaSpec,rho_air,temp,velSpec)
      use kinds
      implicit none
      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), intent(in) :: rho_air, temp
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec      
    end subroutine dia2vel_metek_rain

    subroutine dia2vel_rogers_drops(nDia,diaSpec,rho_air,velSpec)
      use kinds
      implicit none
      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), intent(in) :: rho_air
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
    end subroutine dia2vel_rogers_drops

    subroutine dia2vel_rogers_graupel(nDia,diaSpec,velSpec)
      use kinds
      implicit none
      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec      
    end subroutine dia2vel_rogers_graupel

    subroutine rescale_spectra(nx1,nx2,sort,x1,y1,x2,y2)
      use kinds
      implicit none
      integer, intent(in) :: nx1,nx2
      real(kind=dbl), intent(in), dimension(nx1) :: x1,y1
      real(kind=dbl), intent(in), dimension(nx2) :: x2
      real(kind=dbl), intent(out), dimension(nx2) :: y2
      logical, intent(in) :: sort
    end subroutine rescale_spectra
  end interface
  

  if (verbose .gt. 1) print*, 'Entering radar_spectrum.f90, particle type: ', particle_type

  ! get |K|**2 and lambda
  K2 = dielec_water(0.D0,temp-t_abs,frequency)
  wavelength = c / (frequency*1.d9)   ! m

  diameter_spec_cp(:) = diameter_spec(:)


  rho = rho_air(temp, press)
  my = viscosity_air(temp)/rho !kinematic viscosity

  !get the particle velocities depending on particle type!
  if (particle_type .eq. "cloud") then
    if (verbose .gt. 3) print*, 'using: ', radar_fallVel_cloud, 'to calculate fall velocity for cloud'
    if (radar_fallVel_cloud .eq. "khvorostyanov01_spheres") then
      rho_particle(:) = rho_water
      call dia2vel_khvorostyanov01_spheres(nbins,diameter_spec_cp,rho,my,rho_particle,vel_spec)
    else if (radar_fallVel_cloud .eq. "khvorostyanov01_drops") then
      call dia2vel_khvorostyanov01_drops(nbins,diameter_spec_cp,rho,my,vel_spec)
    else if (radar_fallVel_cloud .eq. "rogers_drops") then
      call dia2vel_rogers_drops(nbins,diameter_spec_cp,rho,vel_spec)
    else if (radar_fallVel_cloud .eq. "pavlos_cloud") then
      call dia2vel_pavlos_cloud(nbins,diameter_spec_cp,vel_spec)
    else
      print*, "did not understand radar_fallVel_cloud ", radar_fallVel_cloud
      stop
    end if

  else if (particle_type .eq. "rain") then 
    if (verbose .gt. 3) print*, 'using: ', radar_fallVel_rain, 'to calculate fall velocity for rain'
    if (radar_fallVel_rain .eq. "khvorostyanov01_drops") then
      call dia2vel_khvorostyanov01_drops(nbins,diameter_spec_cp,rho,my,vel_spec)
    else if (radar_fallVel_rain .eq. "foote69_rain") then
      call dia2vel_foote69_rain(nbins,diameter_spec_cp,rho,temp,vel_spec)
    else if (radar_fallVel_rain .eq. "metek_rain") then
      call dia2vel_metek_rain(nbins,diameter_spec_cp,rho,temp,vel_spec)
    else if (radar_fallVel_rain .eq. "rogers_drops") then
      call dia2vel_rogers_drops(nbins,diameter_spec_cp,rho,vel_spec)
    else
      print*, "did not understand radar_fallVel_rain ", radar_fallVel_rain
      stop
    end if

  else if (particle_type .eq. "ice") then 
    if (verbose .gt. 3) print*, 'using: ', radar_fallVel_ice, 'to calculate fall velocity for ice'
    if (radar_fallVel_ice .eq. "khvorostyanov01_particles") then
      call dia2vel_khvorostyanov01_particles(nbins,diameter_spec_cp,rho,my,&
	      mass_size_a,mass_size_b,area_size_a,area_size_b,vel_spec)
   else 
      print*, "did not understand radar_fallVel_ice ", radar_fallVel_ice
      stop      
    end if

  else if (particle_type .eq. "snow") then 
    if (verbose .gt. 3) print*, 'using: ', radar_fallVel_snow, 'to calculate fall velocity for snow'
    if (radar_fallVel_snow .eq. "khvorostyanov01_particles") then
      call dia2vel_khvorostyanov01_particles(nbins,diameter_spec_cp,rho,my,&
	      mass_size_a,mass_size_b,area_size_a,area_size_b,vel_spec)
   else 
      print*, "did not understand radar_fallVel_snow ", radar_fallVel_snow
      stop      
    end if

  else if (particle_type .eq. "graup") then 
    if (verbose .gt. 3) print*, 'using: ', radar_fallVel_graupel, 'to calculate fall velocity for cloud'
    if (radar_fallVel_graupel .eq. "khvorostyanov01_spheres") then
      !reverse the mass-size relation to get the density assuming spheric shape, nonsense for rain and clouds
      rho_particle = mass_size_a * diameter_spec_cp**( mass_size_b-3.d0) * 6/pi
      call dia2vel_khvorostyanov01_spheres(nbins,diameter_spec_cp,rho,my,rho_particle,vel_spec)
    else if (radar_fallVel_graupel .eq. "rogers_graupel") then
      call dia2vel_rogers_graupel(nbins,diameter_spec_cp,vel_spec)
    else 
      print*, "did not understand radar_fallVel_graupel ", radar_fallVel_graupel
      stop      
    end if

  else if (particle_type .eq. "hail") then 
    if (verbose .gt. 3) print*, 'using: ', radar_fallVel_hail, 'to calculate fall velocity for cloud'
    if (radar_fallVel_hail .eq. "khvorostyanov01_spheres") then
      !reverse the mass-size relation to get the density assuming spheric shape, nonsense for rain and clouds
      rho_particle = mass_size_a * diameter_spec_cp**( mass_size_b-3.d0) * 6/pi
      call dia2vel_khvorostyanov01_spheres(nbins,diameter_spec_cp,rho,my,rho_particle,vel_spec)
    else 
      print*, "did not understand radar_fallVel_hail ", radar_fallVel_hail
      stop  
    end if

  else
    print*, particle_type, ": did not understand particle_type in calc_radar_spectrum"
    stop
  end if

  back_spec_ref = (1d0/ (K2*pi**5) ) * back_spec * (wavelength)**4 ![m⁶/m⁴]
  back_spec_ref =  back_spec_ref * 1d18 !now non-SI: [mm⁶/m³/m]


  del_d = (diameter_spec_cp(2)-diameter_spec_cp(1)) * 1.d3 ![mm]
! print*,"SUM(back_spec_ref * 1d-3)*del_d",SUM(back_spec_ref * 1d-3)*del_d*0.5

! print*,particle_type," Ze log SUM(back_spec_ref)*del_d",10*log10(SUM(back_spec_ref * 1d-3)*del_d),&
!       "assumes equidistant d"
! 
! 
! 
! if (nbins > 2) then
!  call avint( back_spec_ref * 1d-3, diameter_spec_cp * 1.d3, SIZE(diameter_spec_cp), &
!   diameter_spec_cp(1) * 1.d3, diameter_spec_cp(SIZE(diameter_spec_cp)) * 1.d3, Ze )
!  !print*, "Ze", Ze
!  print*, "Ze log INT(back_spec_ref)", 10*log10(Ze)
! else
!  print*, "Ze log INT(back_spec_ref)", "too few bins to integrate"
! end if

Ze = 1d18* (1d0/ (K2*pi**5) ) * back * (wavelength)**4

! print*, "Ze_PAV", Ze*0.5
! print*, "Ze_PAV log", 10*log10(Ze*0.5)

  !move from dimension to velocity!
  do jj=1,nbins-1
    dD_dU(jj) = (diameter_spec_cp(jj+1)-diameter_spec_cp(jj))/(vel_spec(jj+1)-vel_spec(jj)) ![m/(m/s)]
    !is all particles fall with the same velocity, dD_dU gets infinitive!
    if (abs(dD_dU(jj)) .ge. huge(dD_dU(jj))) then
      print*, "Stop in calc_radar_spectrum: dD_dU is infinitive", jj,&
	(diameter_spec_cp(jj+1)-diameter_spec_cp(jj)), (vel_spec(jj+1)-vel_spec(jj))
      stop
    end if
    if (verbose .gt. 3) print*,"jj,dD_dU(jj)",jj,dD_dU(jj)
    del_v_model(jj) = ABS(vel_spec(jj+1)-vel_spec(jj))



  end do
  dD_dU(nbins) = dD_dU(nbins-1)
  del_v_model(nbins) = del_v_model(nbins-1)
  back_vel_spec = back_spec_ref * ABS(dD_dU)  !non-SI: [mm⁶/m³/m * m/(m/s)]

  !get delta velocity
  del_v_radar = (radar_max_V-radar_min_V)/radar_nfft ![m/s]

  min_V_aliased = radar_min_V - radar_aliasing_nyquist_interv*(radar_max_V-radar_min_V)
  max_V_aliased = radar_max_V + radar_aliasing_nyquist_interv*(radar_max_V-radar_min_V)

  !create array from min_v to max_v iwth del_v_radar spacing -> velocity spectrum of radar
  radar_velo_aliased = (/(((ii*del_v_radar)+min_V_aliased),ii=0,radar_nfft_aliased)/) ! [m/s]


  !add vertical air motion to the observations
  if (radar_airmotion) then
    if (verbose .gt. 2) print*, "Averaging spectrum and Adding vertical air motion: ", radar_airmotion_model
    !constant air motion
    if (radar_airmotion_model .eq. "constant") then
      vel_spec = vel_spec + radar_airmotion_vmin
      !interpolate OR average (depending who's bins size is greater) from N(D) bins to radar bins. 
      call rescale_spectra(nbins,radar_nfft_aliased,.false.,vel_spec,back_vel_spec,radar_velo_aliased,particle_spec) ! particle_spec in [mm⁶/m³/m * m/(m/s)]
    !step function
    else if (radar_airmotion_model .eq. "step") then

      allocate(particle_spec_ext(2,radar_nfft_aliased))
      !for vmin
      vel_spec_ext = vel_spec + radar_airmotion_vmin
      back_vel_spec_ext = back_vel_spec * radar_airmotion_step_vmin
      !interpolate OR average (depending who's bins size is greater) from N(D) bins to radar bins. 
      call rescale_spectra(nbins,radar_nfft_aliased,.false.,vel_spec_ext,back_vel_spec_ext,radar_velo_aliased,&
	   particle_spec_ext(1,:))! particle_spec in [mm⁶/m³/m * m/(m/s)]
      !for vmax
      vel_spec_ext = vel_spec + radar_airmotion_vmax
      back_vel_spec_ext = back_vel_spec *(1.d0-radar_airmotion_step_vmin)
      !interpolate OR average (depending who's bins size is greater) from N(D) bins to radar bins. 
      call rescale_spectra(nbins,radar_nfft_aliased,.false.,vel_spec_ext,back_vel_spec_ext,radar_velo_aliased,&
	   particle_spec_ext(2,:))
      !join results
      particle_spec = SUM(particle_spec_ext,1)
      deallocate(particle_spec_ext)
! 
    else if (radar_airmotion_model .eq. "linear") then
      allocate(particle_spec_ext(radar_airmotion_linear_steps,radar_nfft_aliased))
      delta_air = (radar_airmotion_vmax - radar_airmotion_vmin)/REAL(radar_airmotion_linear_steps -1)
!     loop for linear steps
      do jj=1, radar_airmotion_linear_steps
	vel_spec_ext = vel_spec + radar_airmotion_vmin + (jj-1)*delta_air
	back_vel_spec_ext = back_vel_spec / REAL(radar_airmotion_linear_steps)
	!interpolate OR average (depending whos bins size is greater) from N(D) bins to radar bins. 
	call rescale_spectra(nbins,radar_nfft_aliased,.false.,vel_spec_ext,back_vel_spec_ext,radar_velo_aliased,&
	     particle_spec_ext(jj,:))
      end do
      !join results
      particle_spec = SUM(particle_spec_ext,1)
      deallocate(particle_spec_ext)
    else
      print*, "unknown radar_airmotion_model: ", radar_airmotion_model
      stop
    end if
  else
      !no air motion, just rescale
      if (verbose .gt. 2) print*, "AVERAGING particle spectrum"
      call rescale_spectra(nbins,radar_nfft_aliased,.false.,vel_spec,back_vel_spec,radar_velo_aliased,particle_spec) ! particle_spec in [mm⁶/m³/m * m/(m/s)]
  end if



  K = (Ze/SUM(particle_spec*del_v_radar))
  particle_spec = K* particle_spec

  if (verbose .gt. 3) print*,particle_type,"K",K
  if (verbose .gt. 3) print*,particle_type," Ze SUM(back_vel_spec)*del_v_model",10*log10(SUM(back_vel_spec*del_v_model))
  if (verbose .gt. 3) print*,particle_type," Ze SUM(back_vel_spec_ext)*del_v_model",10*log10(SUM(back_vel_spec_ext*del_v_model))
  if (verbose .gt. 3) print*,particle_type," Ze SUM(particle_spec)*del_v_radar",10*log10(SUM(particle_spec)*del_v_radar)




  if (verbose .gt. 1) print*, 'Done radar_spectrum.f90'

  return

end subroutine radar_spectrum