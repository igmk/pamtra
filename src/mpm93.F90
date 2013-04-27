SUBROUTINE mpm93(errstat,freq, Pbkpa, Ekpa, Tc, W, abscof)
  !  MPM93 - subroutines adapted by Jeff Haferman (NASA/GSFC 5/97)
  !  from Liebe's MPM93 model.  His comments are included below.
  !  I've based this adaptation on Frank Evans' MPM92 extraction.
  !
  !---------------------------------------------------------------------
  !       ATMOSPHERIC ATTENUATION AND DELAY RATES UP TO 1000 GHz
  !       June 1993
  !       Hans J. Liebe     (303-497-3310)    
  !       George A. Hufford (       -3457)
  !       Michael G. Cotton (       -7346)
  !       Institute for Telecommunication Sciences
  !       NTIA/ITS.S3 
  !       325 BROADWAY
  !       Boulder, CO  80303,  USA
  !
  !       FAX   :  (303) 497-5993 (ITS), 497-3680 (ITS.S2)
  !       E-Mail:  HLIEBE@NTIA.ITS.BLDRDOC.GOV
  !
  ! COMMENTS:
  ! 
  !   The Millimeter-wave Propagation Model (MPM85) was reported in Ref.
  ! [1]. Molecular absorption by O2, H2O, and N2 is considered, as well as
  ! dielectric loss for haze and fog/cloud conditions (Rayleigh absorption
  ! approximation), and dielectric plus scatter losses (aR**b -
  ! approximation to Mie's theory) under rain conditions. The complex
  ! atmospheric refractivity N (or path-specific rates of attenuation A and
  ! delay B) were continued to be upgraded as discussed in [2] - [7].
  ! 
  !   Features of the current version, MPM93, are:
  ! 
  ! - Haze model to predict the water droplet density for 
  !       U = 80 to 99.95%RH , when a hygroscopic aerosol reference density
  !       wa(80%RH) and a climatic code ('A, B, C, or D') are provided [2],[3]   
  ! 
  ! - Improved model for the dielectric properties of liquid water to
  !       calculate RAYLEIGH absorption and delay by suspended water droplets
  !       for haze, fog, and cloud conditions [6],[7]
  ! 
  ! - Rain attenuation model for Laws & Parsons drop-sizes by Olsen et al. 
  !       [11], and associated dispersive delay, approximated from results 
  !       reported by Zuffery [12]
  ! 
  ! - New temperature-dependent linewidth data (b3 to b6) for the water
  !       vapor lines below 1 THz, and a 5 percent increase in the 
  !       strength b1 of the 22-GHz and 183-GHz lines [9]
  ! 
  ! - New set of line mixing coefficients (a5, a6) for dry air, and 
  !       their improved fit to the extensive 60-GHz lab. data [8],[9]
  ! 
  ! - Revised water vapor saturation pressure equation [10] 
  ! 
  ! - Approximation for Zeeman (O2) [4] and Doppler (H2O) line-broadening
  !       to cover heights up to 100 km.
  ! 
  ! - New pseudo-line water vapor continuum formulation [9]   
  ! 
  ! - Detailed treatment of the anisotropic, mesospheric Zeeman effect
  !   of O2 microwave lines [5]. The ZPM  code [9].
  ! 
  ! 
  !                                 REFERENCES
  ! 
  !  [1] H. Liebe, "An updated model for millimeter-wave propagation in
  !       moist air", Radio Science, vol. 20, no. 5, pp. 1069-1089, 1985.
  ! 
  !  [2] H. Liebe,"A contribution to modeling atmospheric mm-wave properties",
  !       FREQUENZ, vol.41, no. 1/2, pp. 31-36, 1987.
  ! 
  !  [3] H. Liebe and D. Layton, "MM-wave Properties of the Atmosphere:
  !       Laboratory Studies and Propagation Modeling",
  !       NTIA Report 87-224, 80p., Oct. 1987 (NTIS Order No. PB88-164215/AF).
  !       
  !  [4] H. Liebe,"MPM89 - An atmospheric mm-wave propagation model",
  !       Int. J. IR & MM Waves, vol.10, no.6, pp. 631-650, June 1989.
  ! 
  !  [5] G. Hufford and H. Liebe, "MM-Wave Propagation in the Mesosphere",
  !       NTIA Report 89-249, 67p., Sept. 1989 (NTIS Order No. PB90-119868/AS).
  ! 
  !  [6] H. Liebe, T. Manabe, and G. Hufford, "Mm-wave attenuation and delay
  !       rates due to fog/cloud conditions", IEEE Trans. Ant. Prop.,
  !       vol. 37, no. 12, pp. 1617-1623, Dec. 1989.
  ! 
  !  [7] H. Liebe, G. Hufford (ice), and T. Manabe, "A model for the complex
  !       refractivity of water (ice) at frequencies below 1 THz",
  !       Int. J. IR & MM Waves, vol. 12, no. 7, 659-682, 1991.
  !   
  !  [8] H. Liebe, P. Rosenkranz, and G. Hufford, "Atmospheric 60-GHz   
  !       oxygen spectrum: New laboratory measurements and line parameters", 
  !       J. Quant. Spectr. Rad. Transf., vol. 48, no. 5/6, pp. 629-643, 1992.
  ! 
  !  [9] H. Liebe, G. Hufford, and M. Cotton, "Propagation modeling of moist air 
  !       and suspended water/ice particles at frequencies below 1000 GHz", 
  !       Proc. AGARD Conf. Paper 3/1-10, Palma De Mallorca, Spain, May 1993.
  !  
  ! [10] W. Boegel, "Neue Naeherungsgleichungen fuer den Saettigungsdruck des
  !       Wasserdampfes, DFVLR Bericht DLR-FB 77-52, 1977.
  ! 
  ! [11] R.L. Olsen, D.V. Rogers, and D.B. Hodge, "The aRb relation in the
  !       calculation of rain attenuation",
  !       IEEE Trans. Ant. Prop., vol. AP-26, no. 2, pp. 318-329, 1978.
  ! 
  ! [12] C.H. Zuffery, "A study of rain effects on EM waves in the
  !       1 to 600 GHz range", MS-THesis, Dept. Electrical Eng.,
  !       University of Colorado, Boulder,  CO 80309, Feb., 1972.
  !-----------------------------------------------------------------------
  ! 
  !****************************************************************** 
  ! Computes volume absorption coefficient for an atmospheric
  ! layer given the meteorological properties. The allowed frequency
  ! range is from 1 to 1000 GHz.  This routine is hacked from Liebe's
  ! GAS1 subroutine in his MPM93 model, taking out rain and dispersion
  ! computations.  Included is dry air attenuation, oxygen and "psuedo"
  ! water vapor line-continuum absorption, and Rayleigh cloud droplet
  ! absorption. 
  !    Parameters:
  !      freq       frequency (GHz)
  !      Pbkpa   total pressure (kPa)
  !      Ekpa    water vapor pressure (kPa)
  !      Tc      temperature (C)
  !      W       cloud liquid water content (g/m^3)
  !      ABSCOF  absorption coefficient (km^-1)
  !****************************************************************** 
  use kinds, only: dbl, long
  use nml_params, only: verbose

  implicit none

#include "error_report.interface"

  REAL(kind=dbl), intent(in) :: freq, Pbkpa, Ekpa, Tc, W
  real(kind=dbl), intent(out) :: ABSCOF
  INTEGER IFIRST, I, ICE 

  REAL(kind=dbl) AT1, AT2, AT3, AT4 
  REAL(kind=dbl) GAMMA, S, DELTA, So, GAMMAo, Sn 
  REAL(kind=dbl) GAMH, GAMD2, DELH 
  REAL(kind=dbl) fD, fS, Eps, Epinf, Eopt 
  REAL(kind=dbl) Ai, Bi, fice 
  REAL(kind=dbl) V, P, Pb, E 

  COMPLEX ZN, ZNw, ZEp, ZF, ZFo, ZFn 

  ! Common block for oxygen and water vapor lines                         
  REAL(kind=dbl) F0O2 (44), A (6, 44) 
  REAL(kind=dbl) F0H2O (35), B (6, 35) 
  REAL(kind=dbl) A1 (44), A2 (44), A3 (44), A4 (44), A5 (44), A6 (44) 
  REAL(kind=dbl) B1 (35), B2 (35), B3 (35), B4 (35), B5 (35), B6 (35) 
  ! COMMON /MWLINES2/ F0O2,A, F0H2O,B                               

! Error handling

  integer(kind=long), intent(out) :: ErrStat
  character(len=80) :: ErrMsg
  character(len=14) :: NameOfRoutine = 'mpm93'

  DATA IFIRST / 0 / 
  ! hardcoded by JLH                                
  DATA ICE / 0 / 

  !                                                                       
  !     The following data was in mwlines.93.data                         
  DATA F0O2 / 50.474239, 50.987747, 51.503349, 52.021412, 52.542393,&
       53.066906, 53.595749, 54.130001, 54.671158, 55.221367, 55.783802, &
       56.264774, 56.363388, 56.968204, 57.612484, 58.323875, 58.446590, &
       59.164207, 59.590984, 60.306061, 60.434776, 61.150558, 61.800156, &
       62.411217, 62.486259, 62.997978, 63.568520, 64.127769, 64.678902, &
       65.224068, 65.764771, 66.302094, 66.836830, 67.369598, 67.900864, &
       68.431007, 68.960312, 118.750343, 368.498352, 424.763123,         &
       487.249359, 715.393127, 773.839661, 834.145325 /                  
  DATA A1 / 0.094, 0.246, 0.608, 1.414, 3.102, 6.410, 12.470,       &
       22.800, 39.180, 63.160, 95.350, 54.890, 134.400, 176.300, 214.100,&
       238.600, 145.700, 240.400, 211.200, 212.400, 246.100, 250.400,    &
       229.800, 193.300, 151.700, 150.300, 108.700, 73.350, 46.350,      &
       27.480, 15.300, 8.009, 3.946, 1.832, 0.801, 0.330, 0.128, 94.500, &
       6.790, 63.800, 23.500, 9.960, 67.100, 18.000 /                    
  DATA A2 / 9.694, 8.694, 7.744, 6.844, 6.004, 5.224, 4.484, 3.814, &
       3.194, 2.624, 2.119, 0.015, 1.660, 1.260, 0.915, 0.626, 0.084,    &
       0.391, 0.212, 0.212, 0.391, 0.626, 0.915, 1.260, 0.083, 1.665,    &
       2.115, 2.620, 3.195, 3.815, 4.485, 5.225, 6.005, 6.845, 7.745,    &
       8.695, 9.695, 0.009, 0.049, 0.044, 0.049, 0.145, 0.130, 0.147 /   
  DATA A3 / 0.890, 0.910, 0.940, 0.970, 0.990, 1.020, 1.050, 1.070, &
       1.100, 1.130, 1.170, 1.730, 1.200, 1.240, 1.280, 1.330, 1.520,    &
       1.390, 1.430, 1.450, 1.360, 1.310, 1.270, 1.230, 1.540, 1.200,    &
       1.170, 1.130, 1.100, 1.070, 1.050, 1.020, 0.990, 0.970, 0.940,    &
       0.920, 0.900, 1.630, 1.920, 1.930, 1.920, 1.810, 1.820, 1.810 /   
  DATA A4 / 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, &
       0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,    &
       0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,    &
       0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,    &
       0.000, 0.000, 0.000, 0.600, 0.600, 0.600, 0.600, 0.600, 0.600 /   
  DATA A5 / 0.240, 0.220, 0.197, 0.166, 0.136, 0.131, 0.230, 0.335, &
       0.374, 0.258, - 0.166, 0.390, - 0.297, - 0.416, - 0.613, - 0.205, &
       0.748, - 0.722, 0.765, - 0.705, 0.697, 0.104, 0.570, 0.360,       &
       - 0.498, 0.239, 0.108, - 0.311, - 0.421, - 0.375, - 0.267,        &
       - 0.168, - 0.169, - 0.200, - 0.228, - 0.240, - 0.250, - 0.036,    &
       0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /                        
  DATA A6 / 0.790, 0.780, 0.774, 0.764, 0.751, 0.714, 0.584, 0.431, &
       0.305, 0.339, 0.705, - 0.113, 0.753, 0.742, 0.697, 0.051, - 0.146,&
       0.266, - 0.090, 0.081, - 0.324, - 0.067, - 0.761, - 0.777, 0.097, &
       - 0.768, - 0.706, - 0.332, - 0.298, - 0.423, - 0.575, - 0.700,    &
       - 0.735, - 0.744, - 0.753, - 0.760, - 0.765, 0.009, 0.000, 0.000, &
       0.000, 0.000, 0.000, 0.000 /                                      
  DATA F0H2O / 22.235081, 67.803963, 119.995941, 183.310089,        &
       321.225647, 325.152924, 336.222595, 380.197357, 390.134521,       &
       437.346680, 439.150818, 443.018280, 448.001068, 470.888947,       &
       474.689117, 488.491119, 503.568542, 504.482697, 547.676453,       &
       552.020935, 556.935974, 620.700806, 645.866150, 658.005310,       &
       752.033203, 841.053955, 859.962341, 899.306702, 902.616150,       &
       906.207336, 916.171570, 923.118408, 970.315002, 987.926758,       &
       1780.000000 /                                                     
  DATA B1 / 0.01130, 0.00012, 0.00008, 0.24200, 0.00483, 0.14990,   &
       0.00011, 1.15200, 0.00046, 0.00650, 0.09218, 0.01976, 1.03200,    &
       0.03297, 0.12620, 0.02520, 0.00390, 0.00130, 0.97010, 1.47700,    &
       48.74000, 0.50120, 0.00713, 0.03022, 23.96000, 0.00140, 0.01472,  &
       0.00605, 0.00426, 0.01876, 0.83400, 0.00869, 0.89720, 13.21000,   &
       2230.00000 /                                                      
  DATA B2 / 2.143, 8.735, 8.356, 0.668, 6.181, 1.540, 9.829, 1.048, &
       7.350, 5.050, 3.596, 5.050, 1.405, 3.599, 2.381, 2.853, 6.733,    &
       6.733, 0.114, 0.114, 0.159, 2.200, 8.580, 7.820, 0.396, 8.180,    &
       7.989, 7.917, 8.432, 5.111, 1.442, 10.220, 1.920, 0.258, 0.952 /  
  DATA B3 / 2.811, 2.858, 2.948, 3.050, 2.303, 2.783, 2.693, 2.873, &
       2.152, 1.845, 2.100, 1.860, 2.632, 2.152, 2.355, 2.602, 1.612,    &
       1.612, 2.600, 2.600, 3.210, 2.438, 1.800, 3.210, 3.060, 1.590,    &
       3.060, 2.985, 2.865, 2.408, 2.670, 2.900, 2.550, 2.985, 17.620 /  
  DATA B4 / 4.80, 4.93, 4.78, 5.30, 4.69, 4.85, 4.74, 5.38, 4.81,   &
       4.23, 4.29, 4.23, 4.84, 4.57, 4.65, 5.04, 3.98, 4.01, 4.50, 4.50, &
       4.11, 4.68, 4.00, 4.14, 4.09, 5.76, 4.09, 4.53, 5.10, 4.70, 4.78, &
       5.00, 4.94, 4.55, 30.50 /                                         
  DATA B5 / 0.69, 0.69, 0.70, 0.64, 0.67, 0.68, 0.69, 0.54, 0.63,   &
       0.60, 0.63, 0.60, 0.66, 0.66, 0.65, 0.69, 0.61, 0.61, 0.70, 0.70, &
       0.69, 0.71, 0.60, 0.69, 0.68, 0.33, 0.68, 0.68, 0.70, 0.70, 0.70, &
       0.70, 0.64, 0.68, 2.00 /                                          
  DATA B6 / 1.00, 0.82, 0.79, 0.85, 0.54, 0.74, 0.61, 0.89, 0.55,   &
       0.48, 0.52, 0.50, 0.67, 0.65, 0.64, 0.72, 0.43, 0.45, 1.00, 1.00, &
       1.00, 0.68, 0.50, 1.00, 0.84, 0.45, 0.84, 0.90, 0.95, 0.53, 0.78, &
       0.80, 0.67, 0.90, 5.00 /                                          

  !---------------------------------------------------------------------  
  !                                                                       

  if (verbose > 1) print*, 'start of ' // NameOfRoutine

  DO i = 1, 44 
     A (1, i) = A1 (i) 
     A (2, i) = A2 (i) 
     A (3, i) = A3 (i) 
     A (4, i) = A4 (i) 
     A (5, i) = A5 (i) 
     A (6, i) = A6 (i) 
  enddo
  DO i = 1, 35 
     B (1, i) = B1 (i) 
     B (2, i) = B2 (i) 
     B (3, i) = B3 (i) 
     B (4, i) = B4 (i) 
     B (5, i) = B5 (i) 
     B (6, i) = B6 (i) 
  enddo

  ! Only read in line data the first time called                          
  !      IF (IFIRST.EQ.0) THEN                                            
  !        IFIRST = 1                                                     
  !	CALL READLINES2                                                       
  !      ENDIF                                                            

  ! Relative inverse temperature                                          
  V = 300. / (Tc + 273.15) 
  ! This version inputs E.                                                
  ! Note MPM93 has pressure in mb, whereas MPM92 uses kPA                 
  Pb = 10. * Pbkpa 
  E = 10. * Ekpa 
  P = Pb - E 
  IF (P.LT.0) THEN 
     P = 0. 
     Pb = E 
  ENDIF

  ! For OXYGEN                                                            
  ZN = CMPLX (0., 0.) 
  DO 10 I = 1, 44 
     GAMMA = 0. 
     S = A (1, I) * P * V**3 * EXP (A (2, I) * (1. - V) ) * 1.E-6 
     GAMMA = A (3, I) * (P * V** (0.8 - A (4, I) ) + 1.1 * E * V)   &
          * 1.E-3                                                        
     GAMMA = (GAMMA**2 + (25 * 0.6E-4) **2) **0.5 
     DELTA = (A (5, I) + A (6, I) * V) * (P + E) * (V**0.8) * 1.E-3 
     ZF = freq / F0O2 (I) * (CMPLX (1., - DELTA) / CMPLX (F0O2 (I)     &
          - freq, - GAMMA) - CMPLX (1., DELTA) / CMPLX (F0O2 (I) + freq, GAMMA)&
          )                                                              
     ZN = ZN + S * ZF 
10 END DO

  ! OXYGEN LINE ABSORPTION                                                
  ! Cannot be less than 0.                                                
  AT1 = .182 * freq * AIMAG (ZN)
  IF (AT1.LT.0.) AT1 = 0. 
  !                                                                       
  ! DRY AIR CONTINUUM                                                     
  ZN = CMPLX (0., 0.) 
  So = 6.14E-5 * P * V**2 
  GAMMAo = 0.56E-3 * (P + E) * V**.8 
  ZFo = - freq / CMPLX (freq, GAMMAo)
  Sn = 1.40E-12 * p**2 * V**3.5 
  ZFn = CMPLX (0., freq / (1.93E-5 * freq**1.5 + 1.) )
  ZN = So * ZFo + Sn * ZFn 

  ! NONRESONAT DRY AIR ABSORPTION                                         
  AT2 = .182 * freq * AIMAG (ZN)
  !                                                                       
  ! WATER VAPOR                                                           
  ZN = CMPLX (0., 0.) 
  DO 20 I = 1, 35 
     GAMH = 0. 
     S = B (1, I) * E * V**3.5 * EXP (B (2, I) * (1. - V) ) 
     ! Doppler approximation.                                                
     GAMH = B (3, I) * (P * V**B (5, I) + B (4, I) * E * V**B (6, I)&
          ) * 1.E-3                                                      
     GAMD2 = 1E-12 / V * (1.46 * F0H2O (I) ) **2 
     GAMH = 0.535 * GAMH + (0.217 * GAMH**2 + GAMD2) **0.5 
     DELH = 0. 
     ZF = freq / F0H2O (I) * (CMPLX (1., - DELH) / CMPLX (F0H2O (I)    &
          - freq, - GAMH) - CMPLX (1., DELH) / CMPLX (F0H2O (I) + freq, GAMH) )
     ZN = ZN + S * ZF 
20 END DO

  ! WATER VAPOR LINE ABSORPTION                                           
  ! SEE LIEBE'S COMMENT REGARDING "PSUEDO-LINE WATER VAPOR CONTINUUM" - JL
  AT3 = .182 * freq * AIMAG (ZN)

  !>>>>>>>>>>>>>>>> Not used since W is set equal 0 <<<<<<<<<<<<<<<<<<
  !                                                                       
  ! LIQUID WATER PERMITTIVITY [8]                                         
  ! Use exponential form for gamma for T<0 extrapolation (a la Frank Evans
  IF (ICE.EQ.0) THEN 
     !JLH    fD=20.20-146.4*(V-1)+316*(V-1)**2                               
     fD = 20.1 * exp (7.88 * (1 - V) ) 
     fS = 39.8 * fD 
     Eps = 103.3 * (V - 1) + 77.66 
     Epinf = 0.0671 * Eps 
     Eopt = 3.52 
     ! Complex Permittivity of water (double-Debye model)                    
     ZEp = Eps - freq * ( (Eps - Epinf) / CMPLX (freq, fD) + (Epinf -     &
          Eopt) / CMPLX (freq, fS) )
     !                                                                       
     ! ICE PERMITTIVITY [8]                                                  
  ELSE 
     Ai = (62. * V - 11.6) * 1.E-4 * EXP ( - 22.1 * (V - 1.) ) 
     Bi = .542E-6 * ( - 24.17 + 116.79 / V + (V / (V - .9927) ) **2) 
     Eps = 3.15 
     ! Complex Permittivity of Ice                                           
     fice = freq
     IF (freq.LT..001) fice = .001
     ZEp = CMPLX (3.15, Ai / fice+Bi * fice) 
  ENDIF
  ! SUSPENDED PARTICLE RAYLEIGH APPROXIMATION [6]                         
  ZNw = 1.5 * W * ( (ZEp - 1.) / (ZEp + 2.) - 1. + 3. / (Eps + 2.) ) 
  !                                                                       

  ! SUSPENDED WATER DROPLET EXTINCTION                                    
  AT4 = .182 * freq * AIMAG (ZNw)
  !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<                                                                        
  ABSCOF = 0.23026 * (AT1 + AT2 + AT3 + AT4) 

  if (verbose > 1) print*, 'finished in ' // NameOfRoutine

  RETURN 

END SUBROUTINE mpm93
