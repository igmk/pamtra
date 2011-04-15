SUBROUTINE specular_surface(NSTOKES, NUMMU, GROUND_ALBEDO,     &
REFLECT, TRANS, SOURCE)                                           

use kinds
INTEGER NSTOKES, NUMMU 
REAL(kind=dbl) REFLECT(NSTOKES, NUMMU, NSTOKES, NUMMU, 2) 
REAL(kind=dbl) TRANS(NSTOKES, NUMMU, NSTOKES, NUMMU, 2) 
REAL(kind=dbl) SOURCE(NSTOKES, NUMMU, 2) 
REAL(kind=dbl) GROUND_ALBEDO
INTEGER J, N 
REAL(kind=dbl) R1, R2
							    
N = NSTOKES * NUMMU 
CALL MZERO (2 * N, N, REFLECT) 
CALL MZERO (2 * N, 1, SOURCE) 
CALL MIDENTITY (N, TRANS (1, 1, 1, 1, 1) ) 
CALL MIDENTITY (N, TRANS (1, 1, 1, 1, 2) ) 
								  
DO J = 1, NUMMU 
  R1 = GROUND_ALBEDO
  REFLECT (1, J, 1, J, 2) = R1 
  R2 = 0.
  IF (NSTOKES.GT.1) THEN 
      REFLECT (1, J, 2, J, 2) = R2 
      REFLECT (2, J, 1, J, 2) = R2 
      REFLECT (2, J, 2, J, 2) = R1 
  ENDIF
ENDDO 
RETURN 
END SUBROUTINE specular_surface                
								  
								  
SUBROUTINE specular_radiance (NSTOKES, NUMMU, MODE,    &
GROUND_ALBEDO, GROUND_TEMP, WAVELENGTH, RADIANCE)                         

use kinds
INTEGER NSTOKES, NUMMU, MODE 
REAL(kind=dbl) GROUND_TEMP, GROUND_ALBEDO, WAVELENGTH
REAL(kind=dbl) RADIANCE (NSTOKES, NUMMU) 
INTEGER J, N 
REAL(kind=dbl) PLANCK, thermal

! Thermal radiation going up
N = NSTOKES * NUMMU 
CALL MZERO (N, 1, RADIANCE) 
IF (MODE.EQ.0) THEN 
  CALL PLANCK_FUNCTION (GROUND_TEMP, 'R', WAVELENGTH, PLANCK) 
  thermal = (1.0 - GROUND_ALBEDO) * PLANCK
  DO J = 1, NUMMU 
    RADIANCE (1, J) = thermal 
  END DO 
END IF 
								  
RETURN 
END SUBROUTINE specular_radiance
