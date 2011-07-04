      SUBROUTINE LAMBERT_SURFACE (NSTOKES, NUMMU, MODE, MU_VALUES,      &
      QUAD_WEIGHTS, GROUND_ALBEDO, REFLECT, TRANS, SOURCE)              
!      LAMBERT_SURFACE makes the reflection matrix for a Lambertian  
!      surface with albedo (GROUND_ALBEDO).  Also makes the diagonal    
!      transmission and zero source function.                           
!     Input:
!         NSTOKES         number of Stokes components 
!         NUMU            number of zenith angles
!         MODE            Order of Fourier azimuth series:     
!                          0 is azimuthially symmetric case. 
!         MU_VALUES       zenith angles
!         QUAD_WEIGHTS    Weights of the selected quadrature
!         GROUND_ALBEDO   surface albedo
!
!     Output:
!         REFLECT         reflection matrix
!         TRANS           diagonal transmission matrix
!         SOURCE          zero source function

  use kinds
      INTEGER NSTOKES, NUMMU, MODE 
      REAL(kind=dbl) MU_VALUES (NUMMU), QUAD_WEIGHTS (NUMMU) 
      REAL(kind=dbl) GROUND_ALBEDO 
      REAL(kind=dbl) REFLECT (NSTOKES, NUMMU, NSTOKES, NUMMU, 2) 
      REAL(kind=dbl) TRANS (NSTOKES, NUMMU, NSTOKES, NUMMU, 2) 
      REAL(kind=dbl) SOURCE (NSTOKES, NUMMU, 2) 
      INTEGER J1, J2, N 
                                                                        
                                                                        
      N = NSTOKES * NUMMU 
      CALL MZERO (2 * N, N, REFLECT) 
      CALL MZERO (2 * N, 1, SOURCE) 
      CALL MIDENTITY (N, TRANS (1, 1, 1, 1, 1) ) 
      CALL MIDENTITY (N, TRANS (1, 1, 1, 1, 2) ) 
!           The Lambertian ground reflects the flux equally in all directions
!             and completely unpolarizes the radiation                  
      IF (MODE.EQ.0) THEN 
         DO J1 = 1, NUMMU 
         DO J2 = 1, NUMMU 
         REFLECT (1, J1, 1, J2, 2) = 2.0 * GROUND_ALBEDO * MU_VALUES (  &
         J2) * QUAD_WEIGHTS (J2)                                        
         ENDDO 
         ENDDO 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE LAMBERT_SURFACE                
                                                                        
                                                                        
      SUBROUTINE LAMBERT_RADIANCE (NSTOKES, NUMMU, MODE, SRC_CODE,      &
      GROUND_ALBEDO, GROUND_TEMP, WAVELENGTH, DIRECT_SFC_FLUX, RADIANCE)
!        LAMBERT_RADIANCE calculates the ground radiance of a Lambertian
!      surface.  The radiance is due to thermal radiation and reflected 
!      direct solar radiation if there is a solar source (SRC_CODE=1,3).
!      TOTAL_DEPTH is the total optical depth in the direction of the   
!      direct source (sun), with the path length effect already included
  use kinds
      INTEGER NSTOKES, NUMMU, MODE, SRC_CODE 
      REAL(kind=dbl) GROUND_ALBEDO, GROUND_TEMP, WAVELENGTH 
      REAL(kind=dbl) DIRECT_SFC_FLUX 
      REAL(kind=dbl) RADIANCE (NSTOKES, NUMMU) 
      INTEGER J, N 
      REAL(kind=dbl) TMP, THERMAL, PLANCK, PI 
      PARAMETER (PI = 3.1415926535897932384D0) 
                                                                        
      N = NSTOKES * NUMMU 
      CALL MZERO (N, 1, RADIANCE) 
      IF (MODE.EQ.0) THEN 
!           Thermal radiation going up                                  
         IF (SRC_CODE.EQ.2.OR.SRC_CODE.EQ.3) THEN 
            CALL PLANCK_FUNCTION (GROUND_TEMP, 'R', WAVELENGTH, PLANCK)
            THERMAL = (1.0 - GROUND_ALBEDO) * PLANCK 
            DO J = 1, NUMMU 
            RADIANCE (1, J) = THERMAL 
            ENDDO 
         ENDIF 
                                                                        
!           Direct solar reflection (unpolarized)                       
         IF (SRC_CODE.EQ.1.OR.SRC_CODE.EQ.3) THEN 
            TMP = DIRECT_SFC_FLUX * GROUND_ALBEDO / PI 
            DO J = 1, NUMMU 
            RADIANCE (1, J) = RADIANCE (1, J) + TMP 
            ENDDO 
         ENDIF 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE LAMBERT_RADIANCE               
                                                                        
                                                                        
      SUBROUTINE FRESNEL_SURFACE (NSTOKES, NUMMU, MU_VALUES, INDEX,     &
      WAVELENGTH, wind10, REFLECT, TRANS, SOURCE)                                           
!         FRESNEL_REFLECT makes the reflection matrix for a             
!      plane surface with index of refraction (INDEX) using             
!      the Fresnel reflection formulae.  Also makes the diagonal        
!      transmission and zero source function.
                      
  use kinds
      INTEGER NSTOKES, NUMMU 
      REAL(kind=dbl) MU_VALUES (NUMMU) 
      REAL(kind=dbl) REFLECT (NSTOKES, NUMMU, NSTOKES, NUMMU, 2) 
      REAL(kind=dbl) TRANS (NSTOKES, NUMMU, NSTOKES, NUMMU, 2) 
      REAL(kind=dbl) SOURCE (NSTOKES, NUMMU, 2) 
      COMPLEX(kind=dbl) INDEX 
      INTEGER J, N 
      REAL(kind=dbl) COSI, R1, R2, R3, R4 
      COMPLEX(kind=dbl) EPSILON, D, RH, RV
      real(kind=dbl) :: wavelength, freq, wind10, c1, ffoam
      real(kind=dbl) :: c2(2)
      
      logical :: emissivity_correction = .false.
      
      c1 = 1.
      c2 = 0.
      ffoam = 0.
                                                                  
      N = NSTOKES * NUMMU 
      CALL MZERO (2 * N, N, REFLECT) 
      CALL MZERO (2 * N, 1, SOURCE) 
      CALL MIDENTITY (N, TRANS (1, 1, 1, 1, 1) ) 
      CALL MIDENTITY (N, TRANS (1, 1, 1, 1, 2) ) 
      EPSILON = INDEX**2
      DO J = 1, NUMMU 
      COSI = MU_VALUES (J) 
      D = CDSQRT (EPSILON - 1.0D0 + COSI**2) 
      RH = (COSI - D) / (COSI + D) 
      RV = (EPSILON * COSI - D) / (EPSILON * COSI + D)
      ! apply correction due to small and large scale roughness and foam
      ! the corrections are applied to abs(Rv)**2 and abs(Rh)**2 respectively
      if (emissivity_correction) then
	freq = 299792.5/wavelength
	call small_scale(cosi, freq, wind10, c1)
	call large_scale(cosi, freq, wind10, c2)
	call foam(wind10, ffoam)
      end if
      R1 = ((cdABS(RV)**2)*c1-c2(1)+(cdABS(RH)**2)*c1-c2(2))*(1-ffoam) / 2.0D0 
      R2 = ((cdABS(RV)**2)*c1-c2(1)-(cdABS(RH)**2)*c1+c2(2))*(1-ffoam) / 2.0D0
      R3 = DREAL (RV * CONJG (RH) ) 
      R4 = DIMAG (RV * CONJG (RH) )
      REFLECT (1, J, 1, J, 2) = R1 
      IF (NSTOKES.GT.1) THEN 
         REFLECT (1, J, 2, J, 2) = R2 
         REFLECT (2, J, 1, J, 2) = R2 
         REFLECT (2, J, 2, J, 2) = R1 
      ENDIF 
      IF (NSTOKES.GT.2) THEN 
         REFLECT (3, J, 3, J, 2) = R3 
      ENDIF 
      IF (NSTOKES.GT.3) THEN 
         REFLECT (3, J, 4, J, 2) = - R4 
         REFLECT (4, J, 3, J, 2) = R4 
         REFLECT (4, J, 4, J, 2) = R3 
      ENDIF 
!         print*, cosi,R1,R2
      ENDDO 

      RETURN 
      END SUBROUTINE FRESNEL_SURFACE                
                                                                        
                                                                        
      SUBROUTINE FRESNEL_RADIANCE (NSTOKES, NUMMU, MODE, MU_VALUES,     &
      INDEX, GROUND_TEMP, WAVELENGTH, wind10, RADIANCE)                         
!        FRESNEL_RADIANCE calculates the ground radiance of a plane     
!      surface using the Fresnel formulae.  The radiance is due only to 
!      thermal radiation (this subroutine cannot do specular reflection)
      use kinds
      INTEGER NSTOKES, NUMMU, MODE 
      REAL(kind=dbl) GROUND_TEMP, WAVELENGTH, MU_VALUES (NUMMU) 
      REAL(kind=dbl) RADIANCE (NSTOKES, NUMMU) 
      COMPLEX(kind=dbl) INDEX 
      INTEGER J, N 
      REAL(kind=dbl) ZERO, PLANCK, COSI, R1, R2 
      COMPLEX(kind=dbl) EPSILON, D, RH, RV 
      PARAMETER (ZERO = 0.0D0) 
      real(kind=dbl) :: freq, wind10, c1, ffoam
      real(kind=dbl) :: c2(2)
      
      logical :: emissivity_correction = .false.

      c1 = 1.
      c2 = 0.
      ffoam = 0.

! Thermal radiation going up                                  
      N = NSTOKES * NUMMU 
      CALL MZERO (N, 1, RADIANCE) 
      IF (MODE.EQ.0) THEN 
         CALL PLANCK_FUNCTION (GROUND_TEMP, 'R', WAVELENGTH, PLANCK) 
         EPSILON = INDEX**2 
         DO J = 1, NUMMU 
         COSI = MU_VALUES (J) 
         D = CDSQRT (EPSILON - 1.0D0 + COSI**2) 
         RH = (COSI - D) / (COSI + D) 
         RV = (EPSILON * COSI - D) / (EPSILON * COSI + D)
	 ! apply correction due to small and large scale roughness and foam
	 ! the corrections are applied to abs(Rv)**2 and abs(Rh)**2 respectively
	 if (emissivity_correction) then
	   freq = 299792.5/wavelength
	   call small_scale(cosi, freq, wind10, c1)
	   call large_scale(cosi, freq, wind10, c2)
           call foam(wind10, ffoam)
	 end if
         R1 = ((cdABS(RV)**2)*c1-c2(1)+(cdABS(RH)**2)*c1-c2(2))*(1-ffoam) / 2.0D0 
         R2 = ((cdABS(RV)**2)*c1-c2(1)-(cdABS(RH)**2)*c1+c2(2))*(1-ffoam) / 2.0D0
         RADIANCE (1, J) = (1.0 - R1) * PLANCK 
         IF (NSTOKES.GT.1) RADIANCE (2, J) = - R2 * PLANCK 
         ENDDO 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE FRESNEL_RADIANCE
                                                                        
                                                                        
      SUBROUTINE THERMAL_RADIANCE (NSTOKES, NUMMU, MODE, TEMPERATURE,   &
      ALBEDO, WAVELENGTH, RADIANCE)                                     
!        THERMAL_RADIANCE returns a polarized radiance vector for       
!      thermal emission at WAVELENGTH (microns) for a body with         
!      ALBEDO and with a TEMPERATURE (Kelvins).  The radiance is in     
!      the units:  Watts / (meter^2 ster micron).  The emission is      
!      isotropic and unpolarized: only the first term in azimuth Fourier
!      series is nonzero, and all MU terms are the same.                
  use kinds
      INTEGER NSTOKES, NUMMU, MODE 
      REAL(kind=dbl) TEMPERATURE, WAVELENGTH, ALBEDO 
      REAL(kind=dbl) RADIANCE (NSTOKES, NUMMU, 2) 
      INTEGER J, N 
      REAL(kind=dbl) PLANCK, THERMAL 
                                                                        
      N = NSTOKES * NUMMU 
      CALL MZERO (2 * N, 1, RADIANCE) 
      IF (MODE.EQ.0) THEN 
         CALL PLANCK_FUNCTION (TEMPERATURE, 'R', WAVELENGTH, PLANCK) 
         THERMAL = (1.0 - ALBEDO) * PLANCK 
         DO J = 1, NUMMU 
         RADIANCE (1, J, 1) = THERMAL 
         RADIANCE (1, J, 2) = THERMAL 
         ENDDO 
      ENDIF 
      RETURN 
      END SUBROUTINE THERMAL_RADIANCE               
