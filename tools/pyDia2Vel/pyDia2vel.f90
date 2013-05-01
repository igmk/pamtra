subroutine pyDia2vel(model,nbins,diameter_spec,velSpec)
  use settings, only: verbose
  use kinds
  use constants
  use dia2vel
  implicit none

  character(100), intent(in) :: model
  integer, intent(in) :: nbins
  real, intent(in),dimension(nbins) :: diameter_spec
  real, intent(out),dimension(nbins) :: velSpec

  real(kind=dbl) :: temp, press, rho, my, rho_air, viscosity_air,&
    mass_size_a, mass_size_b,&
    area_size_a,area_size_b
  real(kind=dbl),dimension(nbins) :: diameter_spec_cp
  real(kind=dbl),dimension(nbins) :: vel_spec, rho_particle

  integer(kind=long) :: err = 0

  !settings
  !f2py intent(in) :: model,nbins,diameter_spec
  !f2py intent(out) :: velSpec


  temp = 255.37899780273438D0
  press = 63032.068049610127D0

!surface
  temp = 293.15
  press =101300

  rho = rho_air(temp,press)
  my = viscosity_air(temp)/rho

  diameter_spec_cp = DBLE(diameter_spec)

  if (model .eq. "khvorostyanov01_spheres") then
    rho_particle(:) = rho_water
    call dia2vel_khvorostyanov01_spheres(err,nbins,diameter_spec_cp,rho,my,rho_particle,vel_spec)

  else if (model .eq. "khvorostyanov01_graupel") then
    mass_size_b = 3.1d0
    mass_size_a = 169.6d0
    rho_particle = mass_size_a * diameter_spec_cp**( mass_size_b-3.d0) * 6/pi
    call dia2vel_khvorostyanov01_spheres(err,nbins,diameter_spec_cp,rho,my,rho_particle,vel_spec)

  else if (model .eq. "rogers_graupel") then
    call dia2vel_rogers_graupel(err,nbins,diameter_spec_cp,vel_spec)

  else if (model .eq. "khvorostyanov01_hail") then
    mass_size_b = 3.0000030000030002d0
    mass_size_a = 392.32907616514177d0
    rho_particle = mass_size_a * diameter_spec_cp**( mass_size_b-3.d0) * 6/pi
    call dia2vel_khvorostyanov01_spheres(err,nbins,diameter_spec_cp,rho,my,rho_particle,vel_spec)

  else if (model .eq. "khvorostyanov01_particles_waterspheres") then

    mass_size_b = 3.d0
    mass_size_a = 0.524d0
    area_size_b = 2.d0
    area_size_a = 0.785d0
    !CGS to SI
    mass_size_a =mass_size_a * 10.d0**(mass_size_b * 2.d0 - 3.d0)
    area_size_a = area_size_a  * 10.d0**( area_size_b * 2.d0 - 4.d0)

    call dia2vel_khvorostyanov01_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec)       


  else if (model .eq. "khvorostyanov01_particles_rosettes") then

    mass_size_b = 2.26d0
    mass_size_a = 0.00308d0
    area_size_b = 1.57d0
    area_size_a = 0.0869d0
    !CGS to SI
    mass_size_a =mass_size_a * 10.d0**(mass_size_b * 2.d0 - 3.d0)
    area_size_a = area_size_a  * 10.d0**( area_size_b * 2.d0 - 4.d0)

    call dia2vel_khvorostyanov01_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec)       

  else if (model .eq. "khvorostyanov01_particles_aggregates") then

    mass_size_b = 2.1d0
    mass_size_a = 0.0028d0
    area_size_b = 1.88d0
    area_size_a = 0.2285d0
    !CGS to SI
    mass_size_a =mass_size_a * 10.d0**( mass_size_b * 2.d0 - 3)
    area_size_a = area_size_a  * 10.d0**(area_size_b * 2.d0 - 4)

print*, "khvorostyanov01_particles_aggregates: area_size_a", area_size_a
print*, "khvorostyanov01_particles_aggregates: mass_size_a", mass_size_a


    call dia2vel_khvorostyanov01_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec)       

  else if (model .eq. "khvorostyanov01_particles_sector") then

    mass_size_b = 2.02d0
    mass_size_a = 0.00142d0
    area_size_b = 1.97d0
    area_size_a = 0.55d0
    !CGS to SI
    mass_size_a =mass_size_a * 10.d0**(2.d0 * mass_size_b -3.d0)
    area_size_a = area_size_a  * 10.d0**(2.d0 * area_size_b -4.d0)

    call dia2vel_khvorostyanov01_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec) 


  else if (model .eq. "khvorostyanov01_particles_hexPlates") then

    mass_size_b = 2.45d0
    mass_size_a = 0.00739d0
    area_size_b = 2.0d0
    area_size_a = 0.65d0
    !CGS to SI
    mass_size_a =mass_size_a * 10.d0**(2.d0 * mass_size_b -3.d0)
    area_size_a = area_size_a  * 10.d0**(2.d0 * area_size_b -4.d0)

print*, "khvorostyanov01_particles_hexPlates: area_size_a", area_size_a
print*, "khvorostyanov01_particles_hexPlates: mass_size_a", mass_size_a

    call dia2vel_khvorostyanov01_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec)

  else if (model .eq. "heymsfield10_particles_hexPlates") then

    mass_size_b = 2.45d0
    mass_size_a = 0.00739d0
    area_size_b = 1.85d0
    area_size_a = 0.65d0
    !CGS to SI
    mass_size_a =mass_size_a * 10.d0**(2.d0 * mass_size_b -3.d0)
    area_size_a = area_size_a  * 10.d0**(2.d0 * area_size_b -4.d0)

    call dia2vel_heymsfield10_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec) 


  else if (model .eq. "heymsfield10_particles_rosettes") then

    mass_size_b = 2.26d0
    mass_size_a = 0.00308d0
    area_size_b = 1.57d0
    area_size_a = 0.0869d0
    !CGS to SI
    mass_size_a =mass_size_a * 10.d0**(mass_size_b * 2.d0 - 3.d0)
    area_size_a = area_size_a  * 10.d0**( area_size_b * 2.d0 - 4.d0)

    call dia2vel_heymsfield10_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec)       

  else if (model .eq. "heymsfield10_particles_aggregates") then

    mass_size_b = 2.1d0
    mass_size_a = 0.0028d0
    area_size_b = 1.88d0
    area_size_a = 0.2285d0
    !CGS to SI
    mass_size_a =mass_size_a * 10.d0**(2.d0 * mass_size_b -3.d0)
    area_size_a = area_size_a  * 10.d0**(2.d0 * area_size_b -4.d0)

print*, "dia2vel_heymsfield10_particles_aggregates: area_size_a", area_size_a
print*, "dia2vel_heymsfield10_particles_aggregates: mass_size_a", mass_size_a


    call dia2vel_heymsfield10_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec)       

  else if (model .eq. "heymsfield10_particles_MPACE") then

    mass_size_b = 1.7d0
    mass_size_a = 1.07d-10
    area_size_b = 1.88d0
    area_size_a = 0.2285d0
    
    !CGS to SI
    mass_size_a =mass_size_a *   10.d0**(6.d0 * mass_size_b - 3.d0)
    area_size_a = area_size_a  * 10.d0**(2.d0 * area_size_b - 4.d0)

print*, "heymsfield10_particles_MPACE: area_size_a", area_size_a
print*, "heymsfield10_particles_MPACE: mass_size_a", mass_size_a

    call dia2vel_heymsfield10_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec)    

  else if (model .eq. "heymsfield10_particles_sector") then

    mass_size_b = 2.02d0
    mass_size_a = 0.00142d0
    area_size_b = 1.97d0
    area_size_a = 0.55d0
    !CGS to SI
    mass_size_a =mass_size_a * 10.d0**(2.d0 * mass_size_b -3.d0)
    area_size_a = area_size_a  * 10.d0**(2.d0 * area_size_b -4.d0)

    call dia2vel_heymsfield10_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec) 

 else if (model .eq. "heymsfield10_particles_cosmoIce") then

     mass_size_a = 130d0
     mass_size_b = 3.d0
     !area-size relation in SI
     area_size_a = 0.684 !0.684 also in CGS
     area_size_b = 2.d0 !from mitchell 1996 similar to a_msnow&b_snow

    !CGS to SI
area_size_a = area_size_a  * 10.d0**(2.d0 * area_size_b -4.d0)
print*, "heymsfield10_particles_cosmoIce: area_size_a", area_size_a
!      mass_size_a = 0.82d0
!      mass_size_b = 2.5d0
!      !area-size relation in SI
!      area_size_a = 0.4788629555925309d0 !0.24 in CGS
!      area_size_b = 1.85d0 !from mitchell 1996 similar to a_msnow&b_snow


    call dia2vel_heymsfield10_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec) 


  else if (model .eq. "khvorostyanov01_particles_cosmoIce") then

     mass_size_a = 130d0
     mass_size_b = 3.d0
     !area-size relation in SI
     area_size_a = 0.684 !0.684 also in CGS
     area_size_b = 2.d0 !from mitchell 1996 similar to a_msnow&b_snow

    !CGS to SI
area_size_a = area_size_a  * 10.d0**(2.d0 * area_size_b -4.d0)
print*, "khvorostyanov01_particles_cosmoIce: area_size_a", area_size_a
!      mass_size_a = 0.82d0
!      mass_size_b = 2.5d0
!      !area-size relation in SI
!      area_size_a = 0.4788629555925309d0 !0.24 in CGS
!      area_size_b = 1.85d0 !from mitchell 1996 similar to a_msnow&b_snow


    call dia2vel_khvorostyanov01_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec) 


  else if (model .eq. "khvorostyanov01_particles_cosmoSnow") then

        mass_size_b = 2.0d0
        mass_size_a = 0.038d0 ! Locatelli and Hobbs (1974)
	!area-size relation in SI
        area_size_a = 0.2285 
        area_size_b = 1.88d0 !from mitchell 1996 similar to a_msnow&b_snow

area_size_a = area_size_a  * 10.d0**(2.d0 * area_size_b -4.d0)
print*, "khvorostyanov01_particles_cosmoSnow: area_size_a", area_size_a

    call dia2vel_khvorostyanov01_particles(err,nbins,diameter_spec_cp,rho,my,mass_size_a,mass_size_b,&
         area_size_a,area_size_b,vel_spec) 



  else if (model .eq. "khvorostyanov01_drops") then
    call dia2vel_khvorostyanov01_drops(err,nbins,diameter_spec_cp,rho,my,vel_spec)

  else if (model .eq. "foote69_rain") then
    call dia2vel_foote69_rain(err,nbins,diameter_spec_cp,rho,temp,vel_spec)

  else if (model .eq. "metek_rain") then
    call dia2vel_metek_rain(err,nbins,diameter_spec_cp,rho,temp,vel_spec)

  else if (model .eq. "rogers_drops") then
    call dia2vel_rogers_drops(err,nbins,diameter_spec_cp,rho,vel_spec)

  else if (model .eq. "pavlos_cloud") then
    call dia2vel_pavlos_cloud(err,nbins,diameter_spec_cp,vel_spec)

  else
    print*, "did not understand model ", model
    return
  end if

  velSpec = SNGL(vel_spec)

end subroutine pyDia2vel