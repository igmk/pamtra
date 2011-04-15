      SUBROUTINE PLANCK_FUNCTION (TEMP, UNITS, WAVELENGTH, PLANCK) 
!        Calculates the Planck blackbody radiance in                    
!      [Watts /(meter^2 ster micron)] for a temperature in [Kelvins]    
!      at a wavelength in [microns].  If using temperature units then   
!      the Planck function is simple the temperature.                   
  use kinds
      REAL(kind=dbl) TEMP, WAVELENGTH, PLANCK 
      CHARACTER(1) UNITS 
                                                                        
      IF (UNITS.EQ.'T') THEN 
         PLANCK = TEMP 
      ELSE 
         IF (TEMP.GT.0.0) THEN 
            PLANCK = 1.1911D8 / WAVELENGTH**5 / (DEXP (1.4388D4 /       &
            (WAVELENGTH * TEMP) ) - 1)                                  
         ELSE 
            PLANCK = 0.0 
         ENDIF 
      ENDIF 
                                                                        
      RETURN 
      END SUBROUTINE PLANCK_FUNCTION                
                                                                        
                                                                        
      SUBROUTINE GAUSS_LEGENDRE_QUADRATURE (NUM, ABSCISSAS, WEIGHTS) 
!      Generates the abscissas and weights for an even 2*NUM point    
!      Gauss-Legendre quadrature.  Only the NUM positive points are returned
  use kinds
      INTEGER NUM 
      REAL(kind=dbl) ABSCISSAS (1), WEIGHTS (1) 
      INTEGER N, I, J, L 
      REAL(kind=dbl) X, XP, PL, PL1, PL2, DPL, TINY 
      PARAMETER (TINY = 3.0D-14) 
                                                                        
      N = 2 * NUM 
      DO J = 1, NUM 
      X = COS (3.141592654 * (J - .25) / (N + .5) ) 
      I = 0 
  100 CONTINUE 
      PL1 = 1 
      PL = X 
      DO L = 2, N 
      PL2 = PL1 
      PL1 = PL 
      PL = ( (2 * L - 1) * X * PL1 - (L - 1) * PL2) / L 
      ENDDO 
      DPL = N * (X * PL - PL1) / (X * X - 1) 
      XP = X 
      X = XP - PL / DPL 
      I = I + 1 
      IF (ABS (X - XP) .GT. TINY .AND. I .LT. 10) GOTO 100 
      ABSCISSAS (NUM + 1 - J) = X 
      WEIGHTS (NUM + 1 - J) = 2 / ( (1 - X * X) * DPL * DPL) 
      ENDDO 
                                                                        
      RETURN 
      END SUBROUTINE GAUSS_LEGENDRE_QUADRATURE      
                                                                        
                                                                        
      SUBROUTINE DOUBLE_GAUSS_QUADRATURE (NUM, ABSCISSAS, WEIGHTS) 
!        Generates the abscissas and weights for an even 2*NUM point    
!      Gauss-Legendre quadrature.  Only the NUM positive points are returned
  use kinds
      INTEGER NUM 
      REAL(kind=dbl) ABSCISSAS ( * ), WEIGHTS ( * ) 
      INTEGER N, K, I, J, L 
      REAL(kind=dbl) X, XP, PL, PL1, PL2, DPL, TINY 
      PARAMETER (TINY = 3.0D-14) 
                                                                        
      N = NUM 
      K = (N + 1) / 2 
      DO J = 1, K 
      X = COS (3.141592654 * (J - .25) / (N + .5) ) 
      I = 0 
  100 CONTINUE 
      PL1 = 1 
      PL = X 
      DO L = 2, N 
      PL2 = PL1 
      PL1 = PL 
      PL = ( (2 * L - 1) * X * PL1 - (L - 1) * PL2) / L 
      ENDDO 
      DPL = N * (X * PL - PL1) / (X * X - 1) 
      XP = X 
      X = XP - PL / DPL 
      I = I + 1 
      IF (ABS (X - XP) .GT.TINY.AND.I.LT.10) GOTO 100 
      ABSCISSAS (J) = (1D0 - X) / 2D0 
      ABSCISSAS (NUM + 1 - J) = (1D0 + X) / 2D0 
      WEIGHTS (NUM + 1 - J) = 1.0D0 / ( (1.0D0 - X * X) * DPL * DPL) 
      WEIGHTS (J) = 1.0D0 / ( (1.0D0 - X * X) * DPL * DPL) 
      ENDDO 
                                                                        
      RETURN 
      END SUBROUTINE DOUBLE_GAUSS_QUADRATURE        
                                                                        
                                                                        
      SUBROUTINE LOBATTO_QUADRATURE (NUM, ABSCISSAS, WEIGHTS) 
!        Generates the abscissas and weights for an even 2*NUM point    
!      Gauss-Legendre quadrature.  Only the NUM positive points are retu
  use kinds
      INTEGER NUM 
      REAL(kind=dbl) ABSCISSAS ( * ), WEIGHTS ( * ) 
      INTEGER N, N1, I, J, L 
      REAL(kind=dbl) X, XP, PL, PL1, PL2, DPL, D2PL, CI, TINY 
      PARAMETER (TINY = 3.0D-14) 
                                                                        
      N = 2 * NUM 
      N1 = N - 1 
      CI = 0.50 
      IF (MOD (N, 2) .EQ.1) CI = 1.00 
      DO J = 1, NUM - 1 
      X = SIN (3.141592654 * (J - CI) / (N - .5) ) 
      I = 0 
  100 CONTINUE 
      PL1 = 1 
      PL = X 
      DO L = 2, N1 
      PL2 = PL1 
      PL1 = PL 
      PL = ( (2 * L - 1) * X * PL1 - (L - 1) * PL2) / L 
      ENDDO 
      DPL = N1 * (X * PL - PL1) / (X * X - 1) 
      D2PL = (2.D0 * X * DPL - N1 * (N1 + 1) * PL) / (1D0 - X * X) 
      XP = X 
      X = XP - DPL / D2PL 
      I = I + 1 
      IF (ABS (X - XP) .GT.TINY.AND.I.LT.10) GOTO 100 
      ABSCISSAS (J) = X 
      WEIGHTS (J) = 2.0D0 / (N * N1 * PL * PL) 
      ENDDO 
      ABSCISSAS (NUM) = 1.D0 
      WEIGHTS (NUM) = 2.D0 / (N * N1) 
                                                                        
      RETURN 
      END SUBROUTINE LOBATTO_QUADRATURE             
                                                                        
!CC INTERPOLATE THE RESULTS OF RT3 TO the desired mu value              
      SUBROUTINE interp_RT3 (MU_RT3, umu, TB_RT3, TB, maxv) 
  use kinds
      IMPLICIT none 
      INTEGER i, maxv 
      REAL(kind=dbl) MU_RT3 (MAXV), TB_RT3 (MAXV), umu, TB, x1, x2, x3 
                                                                        
      IF (umu.eq.1.0d0) then 
         TB = TB_RT3 (1) 
      ELSE 
         CALL locate (MU_RT3, maxv, umu, i) 
         i = i + 1 
         x1 = MU_RT3 (I - 1) 
         x2 = umu 
         x3 = MU_RT3 (I) 
         TB = TB_RT3 (I - 1) + (TB_RT3 (I) - TB_RT3 (I - 1) ) * (x2 -   &
         x1) / (x3 - x1)                                                
      ENDIF 
      RETURN 
      END SUBROUTINE interp_RT3                     
