subroutine lambert_surface (nstokes, nummu, mode, mu_values,      &
     quad_weights, ground_albedo, reflect, trans, source)              
  !      lambert_surface makes the reflection matrix for a lambertian  
  !      surface with albedo (ground_albedo).  also makes the diagonal    
  !      transmission and zero source function.                           
  !     input:
  !         nstokes         number of stokes components 
  !         numu            number of zenith angles
  !         mode            order of fourier azimuth series:     
  !                          0 is azimuthially symmetric case. 
  !         mu_values       zenith angles
  !         quad_weights    weights of the selected quadrature
  !         ground_albedo   surface albedo
  !
  !     output:
  !         reflect         reflection matrix
  !         trans           diagonal transmission matrix
  !         source          zero source function

  use kinds
  integer nstokes, nummu, mode 
  real(kind=dbl) mu_values (nummu), quad_weights (nummu) 
  real(kind=dbl) ground_albedo 
  real(kind=dbl) reflect (nstokes, nummu, nstokes, nummu, 2) 
  real(kind=dbl) trans (nstokes, nummu, nstokes, nummu, 2) 
  real(kind=dbl) source (nstokes, nummu, 2) 
  integer j1, j2, n 


  n = nstokes * nummu 
  call mzero (2 * n, n, reflect) 
  call mzero (2 * n, 1, source) 
  call midentity (n, trans (1, 1, 1, 1, 1) ) 
  call midentity (n, trans (1, 1, 1, 1, 2) ) 
  !           the lambertian ground reflects the flux equally in all directions
  !             and completely unpolarizes the radiation                  
  if (mode.eq.0) then 
     do j1 = 1, nummu 
        do j2 = 1, nummu 
           reflect (1, j1, 1, j2, 2) = 2.0 * ground_albedo * mu_values (  &
                j2) * quad_weights (j2)                                        
        enddo
     enddo
  endif

  return 
end subroutine lambert_surface


subroutine lambert_radiance (nstokes, nummu, mode, &
     ground_albedo, ground_temp, wavelength,radiance)
  !        lambert_radiance calculates the ground radiance of a lambertian
  !      surface.  the radiance is due to thermal radiation.
  use kinds
  use rt_utilities, only: planck_function
  integer nstokes, nummu, mode
  real(kind=dbl) ground_albedo, ground_temp, wavelength 
  real(kind=dbl) radiance (nstokes, nummu) 
  integer j, n 
  real(kind=dbl) :: thermal, planck, pi
  parameter (pi = 3.1415926535897932384d0) 

  n = nstokes * nummu 
  call mzero (n, 1, radiance) 
  if (mode.eq.0) then 
     !           thermal radiation going up                                  
    call planck_function (ground_temp, 'r', wavelength, planck)
    thermal = (1.0 - ground_albedo) * planck
    do j = 1, nummu
       radiance (1, j) = thermal
    enddo
  endif

  return 
end subroutine lambert_radiance


subroutine fresnel_surface (nstokes, nummu, mu_values, index,reflect, trans, source)
  !         fresnel_reflect makes the reflection matrix for a             
  !      plane surface with index of refraction (index) using             
  !      the fresnel reflection formulae.  also makes the diagonal        
  !      transmission and zero source function.

  use kinds
  integer nstokes, nummu 
  real(kind=dbl) mu_values (nummu) 
  real(kind=dbl) reflect (nstokes, nummu, nstokes, nummu, 2) 
  real(kind=dbl) trans (nstokes, nummu, nstokes, nummu, 2) 
  real(kind=dbl) source (nstokes, nummu, 2) 
  complex(kind=dbl) index 
  integer j, n 
  real(kind=dbl) cosi, r1, r2, r3, r4 
  complex(kind=dbl) epsilon, d, rh, rv

  n = nstokes * nummu 
  call mzero (2 * n, n, reflect) 
  call mzero (2 * n, 1, source) 
  call midentity (n, trans (1, 1, 1, 1, 1) ) 
  call midentity (n, trans (1, 1, 1, 1, 2) ) 
  epsilon = index**2
  do j = 1, nummu 
     cosi = mu_values (j) 
     d = cdsqrt (epsilon - 1.0d0 + cosi**2) 
     rh = (cosi - d) / (cosi + d) 
     rv = (epsilon * cosi - d) / (epsilon * cosi + d)
     r1 = (cdabs(rv)**2+cdabs(rh)**2)/2.0d0
     r2 = (cdabs(rv)**2-cdabs(rh)**2)/2.0d0
     r3 = dreal (rv * conjg (rh) ) 
     r4 = dimag (rv * conjg (rh) )
     reflect (1, j, 1, j, 2) = r1 
     if (nstokes.gt.1) then 
        reflect (1, j, 2, j, 2) = r2 
        reflect (2, j, 1, j, 2) = r2 
        reflect (2, j, 2, j, 2) = r1 
     endif
     if (nstokes.gt.2) then 
        reflect (3, j, 3, j, 2) = r3 
     endif
     if (nstokes.gt.3) then 
        reflect (3, j, 4, j, 2) = - r4 
        reflect (4, j, 3, j, 2) = r4 
        reflect (4, j, 4, j, 2) = r3 
     endif
     !         print*, cosi,r1,r2
  enddo

  return 
end subroutine fresnel_surface


subroutine fresnel_radiance (nstokes, nummu, mode, mu_values,     &
     index, ground_temp, wavelength, radiance)
  !        fresnel_radiance calculates the ground radiance of a plane     
  !      surface using the fresnel formulae.  the radiance is due only to 
  !      thermal radiation (this subroutine cannot do specular reflection)
  use kinds
  use rt_utilities, only: planck_function
  integer nstokes, nummu, mode 
  real(kind=dbl) ground_temp, wavelength, mu_values (nummu) 
  real(kind=dbl) radiance (nstokes, nummu) 
  complex(kind=dbl) index 
  integer j, n 
  real(kind=dbl) zero, planck, cosi, r1, r2 
  complex(kind=dbl) epsilon, d, rh, rv 
  parameter (zero = 0.0d0) 

  ! thermal radiation going up                                  
  n = nstokes * nummu 
  call mzero (n, 1, radiance) 
  if (mode.eq.0) then 
     call planck_function (ground_temp, 'r', wavelength, planck) 
     epsilon = index**2 
     do j = 1, nummu 
        cosi = mu_values (j) 
        d = cdsqrt (epsilon - 1.0d0 + cosi**2) 
        rh = (cosi - d) / (cosi + d) 
        rv = (epsilon * cosi - d) / (epsilon * cosi + d)
        r1 = (cdabs(rv)**2+cdabs(rh)**2)/ 2.0d0
        r2 = (cdabs(rv)**2-cdabs(rh)**2)/ 2.0d0
        radiance (1, j) = (1.0 - r1) * planck 
        if (nstokes.gt.1) radiance (2, j) = - r2 * planck 
     enddo
  endif

  return 
end subroutine fresnel_radiance


subroutine thermal_radiance(nstokes, nummu, mode, temperature,   &
     albedo, wavelength, radiance)                                     
  !        thermal_radiance returns a polarized radiance vector for       
  !      thermal emission at wavelength (microns) for a body with         
  !      albedo and with a temperature (kelvins).  the radiance is in     
  !      the units:  watts / (meter^2 ster micron).  the emission is      
  !      isotropic and unpolarized: only the first term in azimuth fourier
  !      series is nonzero, and all mu terms are the same.                
  use kinds
  use rt_utilities, only: planck_function
  integer nstokes, nummu, mode 
  real(kind=dbl) temperature, wavelength, albedo 
  real(kind=dbl) radiance (nstokes, nummu, 2) 
  integer j, n 
  real(kind=dbl) planck, thermal 

  n = nstokes * nummu 
  call mzero (2 * n, 1, radiance) 
  if (mode.eq.0) then 
     call planck_function (temperature, 'r', wavelength, planck) 
     thermal = (1.0 - albedo) * planck 
     do j = 1, nummu 
        radiance (1, j, 1) = thermal 
        radiance (1, j, 2) = thermal 
     enddo
  endif
  return 
end subroutine thermal_radiance
