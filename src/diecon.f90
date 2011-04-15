!********************************************************************
      SUBROUTINE DIECON (S, T, FREQ, E1, E2) 
!***********************************************************************
! compute the dieletric constant (E1,E2) of a SURFACE given             
! INPUT                                                                 
!     S=salinity in PPM                                                 
!     T=temperature in Celsius                                          
!     freq=frequency in GHz                                             

  use kinds
                                                                        
      IMPLICIT none 
      REAL(kind=dbl) ST, S, T, FREQ, S2, T2, SST, STT, SSTT, ES, TAU, SIGMA,   &
      ZNU, OMEGA, DEN, TWOPI, EINF, SOLD, TOLD, E1, E2                  
                                                                        
      DATA TWOPI / 6.283185307 /, EINF / 4.9 /, SOLD / 0.0 /, TOLD /    &
      - 99. /                                                           
!      save es, tau,sigma                                               
!      IF (S .EQ. SOLD .AND. T .EQ. TOLD) GO TO 10                      
      ST = S * T 
      S2 = S * S 
      T2 = T * T 
      SST = S2 * T 
      STT = T2 * S 
      SSTT = S2 * T2 
      ES = 88. - 4.339E-01 * S + 1.71E-03 * S2 - 4.035E-01 * T +        &
      8.065E-04 * T2 + 6.170E-03 * ST - 8.910E-05 * SST - 6.934E-05 *   &
      STT + 1.439E-06 * SSTT                                            
      TAU = (18.70 - 7.924E-02 * S + 6.35E-04 * S2 - 5.489E-01 * T +    &
      5.758E-03 * T2 + 1.889E-03 * ST - 7.209E-06 * SST - 5.299E-07 *   &
      STT - 2.101E-07 * SSTT) * 1.0E-12                                 
      SIGMA = (7.788E-03 * S - 1.672E-06 * S2 - 8.570E-15 * T +         &
      2.996E-16 * T2 + 4.059E-04 * ST - 3.215E-06 * SST - 1.423E-06 *   &
      STT + 3.229E-08 * SSTT) * 1.0E11                                  
      ZNU = FREQ * 1.E09 
!      write(18,*)'tau,es,sigma',tau,es,sigma                           
      OMEGA = TWOPI * ZNU 
      DEN = 1. + (OMEGA * TAU) **2 
      E1 = (ES - EINF) / DEN + EINF 
      E2 = (ES - EINF) * OMEGA * TAU / DEN + 2. * SIGMA / ZNU 
      SOLD = S 
      TOLD = T 
      RETURN 
      END SUBROUTINE DIECON                         