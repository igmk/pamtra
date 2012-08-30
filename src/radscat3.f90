SUBROUTINE SCATTERING (NUMMU, AZIORDER, NSTOKES, MU_VALUES,       &
     QUAD_WEIGHTS, NUMLEGENDRE, LEGENDRE_COEF, SCAT_NUM, SCATBUF)      
  !        SCATTERING calculates the polarization scattering matrix.      
  !      For each pair of quadrature angles (incoming and outgoing)       
  !      the routine evaluates the single scattering phase matrix         
  !      for many delta phi's.  This is done by calculating the           
  !      scattering angle and summing the Legendre series for             
  !      necessary matrix elements.  Then the polarization reference      
  !      is rotated from the scattering plane to the meridional planes.   
  !      A Fourier transform is done to transform the phase matrices from 
  !      phi space to azimuth mode space.  Because all of the azimuth mode
  !      are calculated at one time the scattering matrix is stored in a  
  !      temporary disk file to be accessed later one mode at a time.     
  !      The scattering matrix includes the quadrature weights for        
  !      integrating.                                                     
  use kinds
  INTEGER NSTOKES, NUMMU, AZIORDER, NUMLEGENDRE, SCAT_NUM 
  REAL(kind=dbl) MU_VALUES (NUMMU), QUAD_WEIGHTS (NUMMU) 
  REAL(kind=dbl) LEGENDRE_COEF (6, 1), SCATBUF ( * ) 
  INTEGER MAXLEG 
  PARAMETER (MAXLEG = 256) 
  INTEGER J1, J2, K, L, M, IX 
  INTEGER NUMPTS, IREC, SCATBASE, NUMELEM, DOSUM (6) 
  REAL(kind=dbl) TMP, MU1, MU2, DELPHI, COS_SCAT 
  REAL(kind=dbl) PHASE_MATRIX (4, 4), OUT_MATRIX (16) 
  REAL(kind=dbl) SCAT_MATRIX (4, 4, 4 * MAXLEG), BASIS_MATRIX (4, 4, 4 *   &
       MAXLEG)                                                           
  REAL(kind=dbl) ZERO, TWOPI 
  PARAMETER (ZERO = 0.0D0, TWOPI = 2.0D0 * 3.1415926535897932384D0) 


  NUMPTS = 2 * 2**INT (LOG (FLOAT (NUMLEGENDRE+4) ) / LOG (2.0)     &
       + 1.0)                                                            
  IF (AZIORDER.EQ.0) NUMPTS = 2 * INT ( (NUMLEGENDRE+1) / 2) + 4
  SCATBASE = (SCAT_NUM - 1) * (AZIORDER + 1) * 2 * (NUMMU**2) 

  !           Find how many Legendre series must be summed                
  CALL NUMBER_SUMS (NSTOKES, NUMLEGENDRE, LEGENDRE_COEF, DOSUM) 

  !       MU1 is the incoming direction, and MU2 is the outgoing direction
  DO J1 = 1, NUMMU 
     TMP = QUAD_WEIGHTS (J1) / 2.0 
     DO J2 = 1, NUMMU 
        DO L = 1, 2 
           MU1 = MU_VALUES (J1) 
           MU2 = MU_VALUES (J2) 
           IF (MOD (L, 2) .EQ.0) MU1 = - MU1 
           !                   Only need to calculate phase matrix for half of     
           !                     the delphi's, the rest come from symmetry.        
           DO K = 1, NUMPTS / 2 + 1 
              DELPHI = (TWOPI * (K - 1) ) / NUMPTS 
              COS_SCAT = MU1 * MU2 + DSQRT ( (1. - MU1**2) * (1. - MU2**2) )    &
                   * DCOS (DELPHI)                                                   
              CALL SUM_LEGENDRE (NUMLEGENDRE, LEGENDRE_COEF, COS_SCAT, DOSUM,   &
                   PHASE_MATRIX)                                                     
              CALL ROTATE_PHASE_MATRIX (PHASE_MATRIX, MU1, MU2, DELPHI,         &
                   COS_SCAT, SCAT_MATRIX (1, 1, K), NSTOKES)                         
              CALL MATRIX_SYMMETRY (NSTOKES, SCAT_MATRIX (1, 1, K), SCAT_MATRIX &
                   (1, 1, NUMPTS - K + 2) )                                          
           ENDDO
           CALL FOURIER_MATRIX (AZIORDER, NUMPTS, NSTOKES, SCAT_MATRIX,      &
                BASIS_MATRIX)                                                     
           NUMELEM = NSTOKES**2 
           DO M = 0, AZIORDER 
              CALL COMBINE_PHASE_MODES (NSTOKES, AZIORDER, M, TMP, BASIS_MATRIX,&
                   OUT_MATRIX)                                                       
              IREC = SCATBASE+M * 2 * (NUMMU**2) + (L - 1) * (NUMMU**2) +       &
                   (J1 - 1) * NUMMU + J2                                             
              IX = NUMELEM * (IREC - 1) 
              DO K = 1, NUMELEM 
                 SCATBUF (IX + K) = OUT_MATRIX (K) 
              ENDDO
           ENDDO

        ENDDO
     ENDDO
  ENDDO

  RETURN 
END SUBROUTINE SCATTERING


SUBROUTINE GET_SCATTERING (NSTOKES, NUMMU, MODE, AZIORDER,        &
     SCAT_NUM, SCATBUF, SCATTER_MATRIX)                                
  !        GET_SCATTERING retrieves the desired azimuth mode of the       
  !      scattering matrix from the temporary disk file.                  
  !      The scattering matrix is actually four matrices:                 
  !      1 is P++, 2 is P+-, 3 is P-+, 4 is P--  where the + and - refer  
  !      to the sign of the quadrature angles. The scattering matrix      
  !      includes the 1/(4*pi) factor, the quadrature weights, and the    
  !      fourier basis constants, as well as the phase matrix.            
  use kinds
  INTEGER NSTOKES, NUMMU, MODE, AZIORDER, SCAT_NUM 
  REAL(kind=dbl) SCATBUF ( * ) 
  REAL(kind=dbl) SCATTER_MATRIX (NSTOKES, NUMMU, NSTOKES, NUMMU, 4) 
  INTEGER J1, J2, I1, I2, L, K 
  INTEGER SCATBASE, IREC, NUMELEM, IX 
  REAL(kind=dbl) IN_MATRIX (16) 

  SCATBASE = (SCAT_NUM - 1) * (AZIORDER + 1) * 2 * (NUMMU**2) 
  NUMELEM = NSTOKES**2 

  DO J1 = 1, NUMMU 
     DO J2 = 1, NUMMU 
        DO L = 1, 2 

           !               Read in the phase matrix for the azimuth mode           
           IREC = SCATBASE+MODE * 2 * (NUMMU**2) + (L - 1) * (NUMMU**2)      &
                + (J1 - 1) * NUMMU + J2                                           
           IX = NUMELEM * (IREC - 1) 
           DO K = 1, NUMELEM 
              IN_MATRIX (K) = SCATBUF (IX + K) 
           ENDDO
           K = 1 
           DO I1 = 1, NSTOKES 
              DO I2 = 1, NSTOKES 
                 SCATTER_MATRIX (I2, J2, I1, J1, L) = IN_MATRIX (K) 
                 K = K + 1 
              ENDDO
           ENDDO

        ENDDO
     ENDDO
  ENDDO

  !           Use the symmetry of the scattering matrix to get            
  !             P-- from P++, and P-+ from P+-.                           
  CALL SCATTER_SYMMETRY (NSTOKES, NUMMU, SCATTER_MATRIX) 

  RETURN 
END SUBROUTINE GET_SCATTERING


SUBROUTINE SCATTER_SYMMETRY (NSTOKES, NUMMU, SCAT) 
  !        SCATTER_SYMMETRY generates the P-- matrix from the P++ matrix  
  !      and the P-+ from the P+- using the symmetry of the scattering mat
  !      For randomly oriented particles with a plane of symmetry the     
  !      diagonal 2 by 2 blocks remain the same while the off-diagonal    
  !      2 by 2 blocks change sign under negation of mu1 and mu2.         
  use kinds
  INTEGER NSTOKES, NUMMU 
  REAL(kind=dbl) SCAT (NSTOKES, NUMMU, NSTOKES, NUMMU, 4) 
  INTEGER I1, I2, J1, J2, SYMDAT (4, 4) 
  DATA SYMDAT / 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1 / 

  DO I1 = 1, NSTOKES 
     DO I2 = 1, NSTOKES 
        IF (SYMDAT (I2, I1) .EQ.1) THEN 
           DO J1 = 1, NUMMU 
              DO J2 = 1, NUMMU 
                 SCAT (I2, J2, I1, J1, 4) = SCAT (I2, J2, I1, J1, 1) 
                 SCAT (I2, J2, I1, J1, 3) = SCAT (I2, J2, I1, J1, 2) 
              ENDDO
           ENDDO
        ELSE 
           DO J1 = 1, NUMMU 
              DO J2 = 1, NUMMU 
                 SCAT (I2, J2, I1, J1, 4) = - SCAT (I2, J2, I1, J1, 1) 
                 SCAT (I2, J2, I1, J1, 3) = - SCAT (I2, J2, I1, J1, 2) 
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  RETURN 
END SUBROUTINE SCATTER_SYMMETRY


SUBROUTINE CHECK_NORM (NSTOKES, NUMMU, QUAD_WEIGHTS,              &
     SCATTER_MATRIX)                                                   
  !        CHECK_NORM checks the normalization of the scattering matrix   
  !      by integrating the I-I term in the polarization matrix over      
  !      the outgoing directions for each incoming direction.             
  !      If the phase function is not normalized it is due either to      
  !      the first term in the Legendre series not being 1, or due        
  !      to too many terms in the Legendre series for the number of       
  !      quadrature angles.                                               
  use kinds
  INTEGER NSTOKES, NUMMU 
  REAL(kind=dbl) QUAD_WEIGHTS (NUMMU) 
  REAL(kind=dbl) SCATTER_MATRIX (NSTOKES, NUMMU, 1, NSTOKES, NUMMU, 1, 4) 
  INTEGER J1, J2, L 
  REAL(kind=dbl) SUM, MAXSUM, TMP 

  MAXSUM = 0.0D0 
  DO J1 = 1, NUMMU 
     TMP = QUAD_WEIGHTS (J1) 
     IF (TMP.NE.0.0) THEN 
        DO L = 1, 2 
           SUM = - 1.0D0 
           DO J2 = 1, NUMMU 
              SUM = SUM + QUAD_WEIGHTS (J2) / TMP * (SCATTER_MATRIX (1, J2,  &
                   1, 1, J1, 1, L) + SCATTER_MATRIX (1, J2, 1, 1, J1, 1, L + 2) ) 
           ENDDO
           MAXSUM = DMAX1 (MAXSUM, DABS (SUM) ) 
        ENDDO
     ENDIF
  ENDDO
  IF (MAXSUM.GT.1.0D-7) THEN 
     WRITE ( * , * ) 'Phase function not normalized:', MAXSUM 
     STOP
  ENDIF
  RETURN 
END SUBROUTINE CHECK_NORM


SUBROUTINE DIRECT_SCATTERING (NUMMU, AZIORDER, NSTOKES, MU_VALUES,&
     NUMLEGENDRE, LEGENDRE_COEF, DIRECT_MU, SCAT_NUM, DIRECTBUF)       
  !        DIRECT_SCATTERING makes the direct (solar) pseudo-source       
  !      vector and stores the results in a temporary buffer.             
  !      The direct vector is the integral of the phase matrix times      
  !      the delta function in the direction of the sun, which is         
  !      the phase matrix, with the incoming direction set to the         
  !      sun's direction, times the {1,0,0,0} Stokes vector.              
  !      The phase matrix is made in the same way as in SCATTERING.       
  use kinds
  INTEGER NSTOKES, NUMMU, AZIORDER, NUMLEGENDRE, SCAT_NUM 
  REAL(kind=dbl) MU_VALUES (NUMMU) 
  REAL(kind=dbl) LEGENDRE_COEF (6, 1), DIRECTBUF ( * ) 
  REAL(kind=dbl) DIRECT_MU 
  INTEGER MAXLEG 
  PARAMETER (MAXLEG = 256) 
  INTEGER J, K, L, M, MC, MS, NUMPTS 
  INTEGER SCATBASE, IX, DOSUM (6) 
  REAL(kind=dbl) MU2, DELPHI, COS_SCAT 
  REAL(kind=dbl) PHASE_MATRIX (4, 4) 
  REAL(kind=dbl) SCAT_MATRIX (4, 4, 2 * MAXLEG), BASIS_MATRIX (4, 4, 2 *   &
       MAXLEG)                                                           
  REAL(kind=dbl) ZERO, TWOPI 
  PARAMETER (ZERO = 0.0D0, TWOPI = 2.0D0 * 3.1415926535897932384D0) 


  NUMPTS = 2 * 2**INT (LOG (FLOAT (NUMLEGENDRE+4) ) / LOG (2.0)     &
       + 1.0)                                                            
  IF (AZIORDER.EQ.0) NUMPTS = 2 * INT ( (NUMLEGENDRE+1) / 2)        &
       + 4                                                               
  SCATBASE = (SCAT_NUM - 1) * (AZIORDER + 1) * 2 * NUMMU 

  !           Find how many Legendre series must be summed                
  CALL NUMBER_SUMS (NSTOKES, NUMLEGENDRE, LEGENDRE_COEF, DOSUM) 


  DO J = 1, NUMMU 
     DO L = 1, 2 
        MU2 = MU_VALUES (J) 
        IF (L.EQ.2) MU2 = - MU2 
        DO K = 1, NUMPTS 
           DELPHI = 0.0 - (TWOPI * (K - 1) ) / NUMPTS 
           COS_SCAT = MU2 * DIRECT_MU + DSQRT ( (1. - MU2**2) * (1. -        &
                DIRECT_MU**2) ) * DCOS (DELPHI)                                   
           CALL SUM_LEGENDRE (NUMLEGENDRE, LEGENDRE_COEF, COS_SCAT, DOSUM,   &
                PHASE_MATRIX)                                                     
           CALL ROTATE_PHASE_MATRIX (PHASE_MATRIX, DIRECT_MU, MU2, DELPHI,   &
                COS_SCAT, SCAT_MATRIX (1, 1, K), NSTOKES)                         
        ENDDO
        CALL FOURIER_MATRIX (AZIORDER, NUMPTS, NSTOKES, SCAT_MATRIX,      &
             BASIS_MATRIX)                                                     

        !               Store away the first column of the combined mode phase m
        IX = NSTOKES * (SCATBASE+ (L - 1) * NUMMU + J - 1) 
        DIRECTBUF (IX + 1) = BASIS_MATRIX (1, 1, 1) 
        IF (NSTOKES.GE.2) DIRECTBUF (IX + 2) = BASIS_MATRIX (2, 1, 1) 
        IF (NSTOKES.GE.3) DIRECTBUF (IX + 3) = 0.0 
        IF (NSTOKES.GE.4) DIRECTBUF (IX + 4) = 0.0 
        DO M = 1, AZIORDER 
           MC = M + 1 
           MS = M + 1 + AZIORDER 
           IX = NSTOKES * (SCATBASE+M * 2 * NUMMU + (L - 1) * NUMMU + J - 1) 
           DIRECTBUF (IX + 1) = BASIS_MATRIX (1, 1, MC) 
           IF (NSTOKES.GE.2) DIRECTBUF (IX + 2) = BASIS_MATRIX (2, 1, MC) 
           IF (NSTOKES.GE.3) DIRECTBUF (IX + 3) = BASIS_MATRIX (3, 1, MS) 
           IF (NSTOKES.GE.4) DIRECTBUF (IX + 4) = BASIS_MATRIX (4, 1, MS) 
        ENDDO

     ENDDO
  ENDDO
  RETURN 
END SUBROUTINE DIRECT_SCATTERING


SUBROUTINE GET_DIRECT (NSTOKES, NUMMU, MODE, AZIORDER, SCAT_NUM,  &
     DIRECTBUF, DIRECT_VECTOR)                                         
  !        GET_DIRECT retrieves the direct (solar) pseudo-source          
  !      vector from the temporary buffer.                                
  use kinds
  INTEGER NSTOKES, NUMMU, MODE, AZIORDER, SCAT_NUM 
  REAL(kind=dbl) DIRECTBUF ( * ), DIRECT_VECTOR (NSTOKES, NUMMU, 2) 
  INTEGER I, J, L, SCATBASE, IX 

  SCATBASE = (SCAT_NUM - 1) * (AZIORDER + 1) * 2 * NUMMU 

  DO J = 1, NUMMU 
     DO L = 1, 2 
        IX = NSTOKES * (SCATBASE+MODE * 2 * NUMMU + (L - 1) * NUMMU + J - &
             1)                                                                
        DO I = 1, NSTOKES 
           DIRECT_VECTOR (I, J, L) = DIRECTBUF (IX + I) 
        ENDDO
     ENDDO
  ENDDO
  RETURN 
END SUBROUTINE GET_DIRECT


SUBROUTINE NUMBER_SUMS (NSTOKES, NLEGEN, COEF, DOSUM) 
  use kinds
  INTEGER NSTOKES, NLEGEN, DOSUM (6) 
  REAL(kind=dbl) COEF (6, * ) 
  INTEGER I, SCAT, CASE, SUMCASES (6, 5) 
  REAL(kind=dbl) ZERO 
  PARAMETER (ZERO = 0.0D0) 
  DATA SUMCASES / 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0,&
       0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1 /                           

  !           SCAT: 1 is Rayleigh, 2 is Mie, 3 is general                 
  SCAT = 1 
  DO I = 1, NLEGEN + 1 
     IF (COEF (4, I) .NE.ZERO) SCAT = 2 
  ENDDO
  DO I = 1, NLEGEN + 1 
     IF (COEF (1, I) .NE.COEF (5, I) .OR.COEF (3, I) .NE.COEF (6, I) ) &
          SCAT = 3                                                          
  ENDDO

  IF (NSTOKES.EQ.1) THEN 
     CASE = 1 
  ELSEIF (NSTOKES.LE.3) THEN 
     CASE = 2 
     IF (SCAT.EQ.3) CASE = 4 
  ELSE 
     CASE = 2 
     IF (SCAT.EQ.2) CASE = 3 
     IF (SCAT.EQ.3) CASE = 5 
  ENDIF
  DO I = 1, 6 
     DOSUM (I) = SUMCASES (I, CASE) 
  ENDDO

  RETURN 
END SUBROUTINE NUMBER_SUMS


SUBROUTINE SUM_LEGENDRE (NLEGEN, COEF, X, DOSUM, PHASE_MATRIX) 
  !       SUM_LEGENDRE sums the Legendre series for each element of the   
  !       phase matrix using X for the scattering angle.  There are       
  !       only six sets of Legendre coefficients in COEF since there are  
  !       six independant parameters for randomly oriented particles      
  !       with a plane of symmetry.  The constant array DOSUM controls which
  !       of the six series are summed.  The ROW and COL arrays give the  
  !       location of each of the series in the phase matrix.             
  use kinds
  INTEGER NLEGEN, DOSUM (6) 
  REAL(kind=dbl) COEF (6, 1), X, PHASE_MATRIX (4, 4) 
  INTEGER I, L, M 
  INTEGER ROW (6), COL (6) 
  REAL(kind=dbl) SUM, PL, PL1, PL2 
  DATA ROW / 1, 1, 3, 3, 2, 4 /, COL / 1, 2, 3, 4, 2, 4 / 

  !           Sum the Legendre series                                     
  DO I = 1, 6 
     SUM = 0.0 
     IF (DOSUM (I) .EQ.1) THEN 
        PL1 = 1.0 
        PL = 1.0 
        DO L = 0, NLEGEN 
           M = L + 1 
           IF (L.GT.0) PL = (2 * L - 1) * X * PL1 / L - (L - 1) * PL2 / L 
           SUM = SUM + COEF (I, M) * PL 
           PL2 = PL1 
           PL1 = PL 
        ENDDO
     ENDIF
     PHASE_MATRIX (ROW (I), COL (I) ) = SUM 
  ENDDO
  PHASE_MATRIX (2, 1) = PHASE_MATRIX (1, 2) 
  PHASE_MATRIX (4, 3) = - PHASE_MATRIX (3, 4) 
  IF (DOSUM (5) .EQ.0) PHASE_MATRIX (2, 2) = PHASE_MATRIX (1, 1) 
  IF (DOSUM (6) .EQ.0) PHASE_MATRIX (4, 4) = PHASE_MATRIX (3, 3) 

  RETURN 
END SUBROUTINE SUM_LEGENDRE


SUBROUTINE ROTATE_PHASE_MATRIX (PHASE_MATRIX1, MU1, MU2, DELPHI,  &
     COS_SCAT, PHASE_MATRIX2, NSTOKES)                                 
  !        ROTATE_PHASE_MATRIX applies the rotation of the polarization   
  !      basis from the incident plane into the scattering plane and      
  !      from the scattering plane to the outgoing plane.                 
  !      MU1 is the incoming direction, and MU2 is the outgoing direction.
  !      Currently, set up for a phase matrix for randomly oriented partic
  !      with a plane of symmetry - only 6 unique parameters.             
  use kinds
  INTEGER NSTOKES 
  REAL(kind=dbl) PHASE_MATRIX1 (4, 4), PHASE_MATRIX2 (4, 4) 
  REAL(kind=dbl) MU1, MU2, DELPHI, COS_SCAT 
  REAL(kind=dbl) SIN_SCAT, SIN_THETA1, SIN_THETA2 
  REAL(kind=dbl) SINPHI, COSPHI, SIN1, SIN2, COS1, COS2 
  REAL(kind=dbl) SIN21, COS21, SIN22, COS22, A1, A2, A3, A4, B1, B2 
  REAL(kind=dbl) ZERO 
  PARAMETER (ZERO = 0.0D0) 

  A1 = PHASE_MATRIX1 (1, 1) 
  PHASE_MATRIX2 (1, 1) = A1 
  IF (NSTOKES.EQ.1) RETURN 

  SIN_SCAT = DSQRT (MAX (0.0D0, 1. - COS_SCAT**2) ) 
  SIN_THETA1 = DSQRT (1. - MU1**2) 
  SIN_THETA2 = DSQRT (1. - MU2**2) 
  SINPHI = DSIN (DELPHI) 
  COSPHI = DCOS (DELPHI) 
  IF (SIN_SCAT.EQ.ZERO) THEN 
     SIN1 = 0.0 
     SIN2 = 0.0 
     COS1 = 1.0 
     COS2 = - 1.0 
  ELSE 
     SIN1 = SIN_THETA2 * SINPHI / SIN_SCAT 
     SIN2 = SIN_THETA1 * SINPHI / SIN_SCAT 
     COS1 = (SIN_THETA1 * MU2 - SIN_THETA2 * MU1 * COSPHI) /        &
          SIN_SCAT                                                       
     COS2 = (SIN_THETA2 * MU1 - SIN_THETA1 * MU2 * COSPHI) /        &
          SIN_SCAT                                                       
  ENDIF
  SIN21 = 2.0 * SIN1 * COS1 
  COS21 = 1.0 - 2.0 * SIN1**2 
  SIN22 = 2.0 * SIN2 * COS2 
  COS22 = 1.0 - 2.0 * SIN2**2 

  IF (NSTOKES.GT.1) THEN 
     A2 = PHASE_MATRIX1 (2, 2) 
     A3 = PHASE_MATRIX1 (3, 3) 
     B1 = PHASE_MATRIX1 (1, 2) 
     PHASE_MATRIX2 (1, 2) = B1 * COS21 
     PHASE_MATRIX2 (2, 1) = B1 * COS22 
     PHASE_MATRIX2 (2, 2) = A2 * COS21 * COS22 - A3 * SIN21 * SIN22 
  ENDIF
  IF (NSTOKES.GT.2) THEN 
     PHASE_MATRIX2 (1, 3) = - B1 * SIN21 
     PHASE_MATRIX2 (2, 3) = - A2 * SIN21 * COS22 - A3 * COS21 *     &
          SIN22                                                          
     PHASE_MATRIX2 (3, 1) = B1 * SIN22 
     PHASE_MATRIX2 (3, 2) = A2 * COS21 * SIN22 + A3 * SIN21 * COS22 
     PHASE_MATRIX2 (3, 3) = - A2 * SIN21 * SIN22 + A3 * COS21 *     &
          COS22                                                          
  ENDIF
  IF (NSTOKES.GT.3) THEN 
     A4 = PHASE_MATRIX1 (4, 4) 
     B2 = PHASE_MATRIX1 (3, 4) 
     PHASE_MATRIX2 (1, 4) = 0.0 
     PHASE_MATRIX2 (2, 4) = - B2 * SIN22 
     PHASE_MATRIX2 (3, 4) = B2 * COS22 
     PHASE_MATRIX2 (4, 1) = 0.0 
     PHASE_MATRIX2 (4, 2) = - B2 * SIN21 
     PHASE_MATRIX2 (4, 3) = - B2 * COS21 
     PHASE_MATRIX2 (4, 4) = A4 
  ENDIF

  RETURN 
END SUBROUTINE ROTATE_PHASE_MATRIX


SUBROUTINE MATRIX_SYMMETRY (NSTOKES, MATRIX1, MATRIX2) 
  !        MATRIX_SYMMETRY performs a symmetry operation on a             
  !      phase matrix.  The operation consists of negating the            
  !      off-diagonal 2 by 2 blocks.  This operation is equivalent        
  !      to negating (mu) and (mu') or negating (phi'-phi).               
  use kinds
  INTEGER NSTOKES 
  REAL(kind=dbl) MATRIX1 (4, 4), MATRIX2 (4, 4) 

  MATRIX2 (1, 1) = MATRIX1 (1, 1) 
  IF (NSTOKES.GT.1) THEN 
     MATRIX2 (1, 2) = MATRIX1 (1, 2) 
     MATRIX2 (2, 1) = MATRIX1 (2, 1) 
     MATRIX2 (2, 2) = MATRIX1 (2, 2) 
  ENDIF
  IF (NSTOKES.GT.2) THEN 
     MATRIX2 (1, 3) = - MATRIX1 (1, 3) 
     MATRIX2 (2, 3) = - MATRIX1 (2, 3) 
     MATRIX2 (3, 1) = - MATRIX1 (3, 1) 
     MATRIX2 (3, 2) = - MATRIX1 (3, 2) 
     MATRIX2 (3, 3) = MATRIX1 (3, 3) 
  ENDIF
  IF (NSTOKES.GT.3) THEN 
     MATRIX2 (1, 4) = - MATRIX1 (1, 4) 
     MATRIX2 (2, 4) = - MATRIX1 (2, 4) 
     MATRIX2 (3, 4) = MATRIX1 (3, 4) 
     MATRIX2 (4, 1) = - MATRIX1 (4, 1) 
     MATRIX2 (4, 2) = - MATRIX1 (4, 2) 
     MATRIX2 (4, 3) = MATRIX1 (4, 3) 
     MATRIX2 (4, 4) = MATRIX1 (4, 4) 
  ENDIF

  RETURN 
END SUBROUTINE MATRIX_SYMMETRY


SUBROUTINE COMBINE_PHASE_MODES (NSTOKES, AZIORDER, M, TMP,        &
     BASIS_MATRIX, OUT_MATRIX)                                         
  use kinds
  INTEGER NSTOKES, AZIORDER, M 
  REAL(kind=dbl) TMP, BASIS_MATRIX (4, 4, * ), OUT_MATRIX ( * ) 
  INTEGER I1, I2, K, MC, MS, SINFLAG (4, 4) 
  REAL(kind=dbl) C 
  DATA SINFLAG / 0, 0, - 1, - 1, 0, 0, - 1, - 1, 1, 1, 0, 0, 1, 1,  &
       0, 0 /                                                            

  IF (M.EQ.0) THEN 
     K = 1 
     DO I1 = 1, NSTOKES 
        DO I2 = 1, NSTOKES 
           IF (SINFLAG (I2, I1) .EQ.0) THEN 
              OUT_MATRIX (K) = TMP * BASIS_MATRIX (I2, I1, 1) 
           ELSE 
              OUT_MATRIX (K) = 0.0 
           ENDIF
           K = K + 1 
        ENDDO
     ENDDO
  ELSE 
     MC = M + 1 
     MS = M + 1 + AZIORDER 
     C = 0.5 * TMP 
     K = 1 
     DO I1 = 1, NSTOKES 
        DO I2 = 1, NSTOKES 
           IF (SINFLAG (I2, I1) .EQ.0) THEN 
              OUT_MATRIX (K) = C * BASIS_MATRIX (I2, I1, MC) 
           ELSEIF (SINFLAG (I2, I1) .EQ. - 1) THEN 
              OUT_MATRIX (K) = - C * BASIS_MATRIX (I2, I1, MS) 
           ELSE 
              OUT_MATRIX (K) = C * BASIS_MATRIX (I2, I1, MS) 
           ENDIF
           K = K + 1 
        ENDDO
     ENDDO
  ENDIF

  RETURN 
END SUBROUTINE COMBINE_PHASE_MODES


SUBROUTINE FOURIER_MATRIX (AZIORDER, NUMPTS, NSTOKES, REAL_MATRIX,&
     BASIS_MATRIX)                                                     
  use kinds
  INTEGER AZIORDER, NUMPTS, NSTOKES 
  REAL(kind=dbl) REAL_MATRIX (4, 4, 1), BASIS_MATRIX (4, 4, 1) 
  INTEGER MAXLEG 
  PARAMETER (MAXLEG = 256) 
  REAL(kind=dbl) BASIS_VECTOR (4 * MAXLEG), REAL_VECTOR (4 * MAXLEG) 
  INTEGER NUMAZI, I, J, K, M 

  NUMAZI = 2 * AZIORDER + 1 
  DO I = 1, NSTOKES 
     DO J = 1, NSTOKES 
        DO K = 1, NUMPTS 
           REAL_VECTOR (K) = REAL_MATRIX (I, J, K) 
        ENDDO
        CALL FOURIER_BASIS (NUMAZI, AZIORDER, NUMPTS, + 1, BASIS_VECTOR,  &
             REAL_VECTOR)                                                      
        DO M = 1, NUMAZI 
           BASIS_MATRIX (I, J, M) = BASIS_VECTOR (M) 
        ENDDO
     ENDDO
  ENDDO
  RETURN 
END SUBROUTINE FOURIER_MATRIX


SUBROUTINE FOURIER_BASIS (NUMBASIS, ORDER, NUMPTS, DIRECTION,     &
     BASIS_VECTOR, REAL_VECTOR)                                        
  !       FOURIER_BASIS converts a vector between the Fourier basis       
  !      and azimuth space.  The Fourier basis functions are:             
  !      { 1, cos(x), cos(2x), . ., cos(Mx), sin(x), sin(2x), . . sin(Mx) 
  !      where M is the order of the Fourier basis (ORDER).  The number   
  !      of elements in the basis in NUMBASIS.  Generally, NUMBASIS=2*ORDE
  !      but for even functions NUMBASIS=ORDER+1                          
  !      DIRECTION is negative for conversion to real space, and positive 
  !      for conversion into Fourier basis.                               
  !      NUMPTS is the number of points in the real space (must be power o
  use kinds
  INTEGER NUMBASIS, NUMPTS, ORDER, DIRECTION 
  REAL(kind=dbl) BASIS_VECTOR ( * ), REAL_VECTOR ( * ) 
  INTEGER I, BASISLEN 
  REAL(kind=dbl) SUM 

  BASISLEN = MIN (ORDER, NUMPTS / 2 - 1) 
  IF (DIRECTION.LT.0) THEN 
     IF (ORDER.EQ.0) THEN 
        DO I = 1, NUMPTS 
           REAL_VECTOR (I) = BASIS_VECTOR (1) 
        ENDDO
     ELSE 
        DO I = 1, NUMPTS 
           REAL_VECTOR (I) = 0.0 
        ENDDO
        REAL_VECTOR (1) = BASIS_VECTOR (1) 
        DO I = 1, BASISLEN 
           REAL_VECTOR (2 * I + 1) = BASIS_VECTOR (I + 1) / 2 
        ENDDO
        IF (NUMBASIS.GT.ORDER + 1) THEN 
           DO I = 1, BASISLEN 
              REAL_VECTOR (2 * I + 2) = BASIS_VECTOR (I + ORDER + 1)   &
                   / 2                                                      
           ENDDO
        ENDIF
        CALL FFT1DR (REAL_VECTOR, NUMPTS, - 1) 
     ENDIF

  ELSE 

     IF (ORDER.EQ.0) THEN 
        SUM = 0.0 
        DO I = 1, NUMPTS 
           SUM = SUM + REAL_VECTOR (I) 
        ENDDO
        BASIS_VECTOR (1) = SUM / NUMPTS 
     ELSE 
        CALL FFT1DR (REAL_VECTOR, NUMPTS, + 1) 
        BASIS_VECTOR (1) = REAL_VECTOR (1) / NUMPTS 
        DO I = 1, BASISLEN 
           BASIS_VECTOR (I + 1) = 2 * REAL_VECTOR (2 * I + 1) / NUMPTS 
        ENDDO
        DO I = BASISLEN + 1, ORDER 
           BASIS_VECTOR (I + 1) = 0.0D0 
        ENDDO
        IF (NUMBASIS.GT.ORDER + 1) THEN 
           DO I = 1, BASISLEN 
              BASIS_VECTOR (I + ORDER + 1) = 2 * REAL_VECTOR (2 * I +  &
                   2) / NUMPTS                                              
           ENDDO
           DO I = BASISLEN + 1, ORDER 
              BASIS_VECTOR (I + ORDER + 1) = 0.0D0 
           ENDDO
        ENDIF
     ENDIF
  ENDIF

  RETURN 
END SUBROUTINE FOURIER_BASIS


SUBROUTINE FFT1DR (DATA, N, ISIGN) 
  !        Real 1D FFT.  N must be a power of two.                        
  !      If ISIGN=+1 then real to complex conjugate FFT is done.          
  !      If ISIGN=-1 then complex conjugate to real FFT is done.          
  !      The Nyquist frequency component is returned in the first imaginar
  !      element.   No normalization (by N) is performed.                 
  use kinds
  INTEGER N, ISIGN 
  REAL(kind=dbl) DATA ( * ) 
  INTEGER MAXN 
  PARAMETER (MAXN = 512) 
  REAL(kind=dbl) PHASE (4 * MAXN) 
  INTEGER MN 
  REAL(kind=dbl) NYQUIST (2) 
  SAVE MN, PHASE 

  IF (MN.LT.N) THEN 
     MN = N 
     IF (MN.GT.MAXN) STOP 'Phase array too small' 
     CALL MAKEPHASE (PHASE, MN) 
  ENDIF

  IF (ISIGN.GT.0) THEN 
     !           Forward transform:  real to complex-conjugate               
     CALL FFTC (DATA, N / 2, PHASE (1) ) 
     CALL FIXREAL (DATA, NYQUIST, N / 2, + 1, PHASE (1) ) 
     DATA (2) = NYQUIST (1) 
  ELSE 
     !           Inverse transform:  complex-conjugate to real               
     NYQUIST (1) = DATA (2) 
     CALL FIXREAL (DATA, NYQUIST, N / 2, - 1, PHASE (2 * MN + 1) ) 
     CALL FFTC (DATA, N / 2, PHASE (2 * MN + 1) ) 
  ENDIF

  RETURN 
END SUBROUTINE FFT1DR


SUBROUTINE FFTC (DATA, N, PHASE) 
  use kinds
  INTEGER N 
  REAL(kind=dbl) DATA ( * ) 
  REAL(kind=dbl) PHASE ( * ) 
  INTEGER I, IREV, M, J, K, M0, M1 
  INTEGER JMAX, POWER, IPH 
  REAL(kind=dbl) TMPR, TMPI, PHR, PHI 

  IF (N.LE.1) RETURN 
  IREV = 0 
  DO I = 0, N - 1 
     IF (I.GT.IREV) THEN 
        M0 = 2 * I + 1 
        M1 = 2 * IREV + 1 
        TMPR = DATA (M0) 
        TMPI = DATA (M0 + 1) 
        DATA (M0) = DATA (M1) 
        DATA (M0 + 1) = DATA (M1 + 1) 
        DATA (M1) = TMPR 
        DATA (M1 + 1) = TMPI 
     ENDIF
     M = N 
50   CONTINUE 
     M = M / 2 
     IF (IREV.LT.M) GOTO 70 
     IREV = IREV - M 
     IF (M.GT.1) GOTO 50 
70   IREV = IREV + M 
  ENDDO

  JMAX = N 
  POWER = 1 
200 CONTINUE 
  JMAX = JMAX / 2 
  M0 = 1 
  M1 = POWER * 2 + 1 
  DO J = 1, JMAX 
     IPH = 2 * POWER 
     DO K = 1, POWER 
        PHR = PHASE (IPH - 1) 
        PHI = PHASE (IPH) 
        IPH = IPH + 2 
        TMPR = PHR * DATA (M1) - PHI * DATA (M1 + 1) 
        TMPI = PHI * DATA (M1) + PHR * DATA (M1 + 1) 
        DATA (M1) = DATA (M0) - TMPR 
        DATA (M1 + 1) = DATA (M0 + 1) - TMPI 
        DATA (M0) = DATA (M0) + TMPR 
        DATA (M0 + 1) = DATA (M0 + 1) + TMPI 
        M0 = M0 + 2 
        M1 = M1 + 2 
     ENDDO
     M0 = M0 + POWER * 2 
     M1 = M1 + POWER * 2 
  ENDDO
  POWER = 2 * POWER 
  IF (JMAX.GT.1) GOTO 200 

  RETURN 
END SUBROUTINE FFTC


SUBROUTINE MAKEPHASE (PHASE, NMAX) 
  use kinds
  INTEGER NMAX 
  REAL(kind=dbl) PHASE ( * ) 
  INTEGER I, J, N 
  DOUBLEPRECISION PI, F 
  PARAMETER (PI = 3.1415926535897932D0) 

  J = 1 
  N = 1 
100 CONTINUE 
  F = PI / N 
  DO I = 0, N - 1 
     PHASE (J) = DCOS (F * DFLOAT (I) ) 
     PHASE (J + 1) = DSIN (F * DFLOAT (I) ) 
     J = J + 2 
  ENDDO
  N = 2 * N 
  IF (N.LT.NMAX) GOTO 100 

  J = 2 * NMAX + 1 
  N = 1 
200 CONTINUE 
  F = - PI / N 
  DO I = 0, N - 1 
     PHASE (J) = DCOS (F * DFLOAT (I) ) 
     PHASE (J + 1) = DSIN (F * DFLOAT (I) ) 
     J = J + 2 
  ENDDO
  N = 2 * N 
  IF (N.LT.NMAX) GOTO 200 
  RETURN 
END SUBROUTINE MAKEPHASE


SUBROUTINE FIXREAL (DATA, NYQUIST, N, ISIGN, PHASE) 
  use kinds
  INTEGER N, ISIGN 
  REAL(kind=dbl) DATA ( * ), NYQUIST (2) 
  REAL(kind=dbl) PHASE ( * ) 
  INTEGER I, IPH, M, MC 
  REAL(kind=dbl) TMP0R, TMP0I, TMP1R, TMP1I, PHR, PHI 

  IPH = 2 * N + 1 
  M = 3 
  MC = 2 * N - 1 
  IF (ISIGN.GT.0) THEN 
     NYQUIST (1) = DATA (1) - DATA (2) 
     NYQUIST (2) = 0 
     DATA (1) = DATA (1) + DATA (2) 
     DATA (2) = 0 
     DO I = 2, N / 2 + 1 
        PHR = PHASE (IPH) 
        PHI = PHASE (IPH + 1) 
        IPH = IPH + 2 
        TMP0R = DATA (M) + DATA (MC) 
        TMP0I = DATA (M + 1) - DATA (MC + 1) 
        TMP1R = - PHI * (DATA (M) - DATA (MC) ) - PHR * (DATA (M + 1)  &
             + DATA (MC + 1) )                                              
        TMP1I = PHR * (DATA (M) - DATA (MC) ) - PHI * (DATA (M + 1)    &
             + DATA (MC + 1) )                                              
        DATA (M) = 0.5 * (TMP0R - TMP1R) 
        DATA (M + 1) = 0.5 * (TMP0I - TMP1I) 
        DATA (MC) = 0.5 * (TMP0R + TMP1R) 
        DATA (MC + 1) = - 0.5 * (TMP0I + TMP1I) 
        M = M + 2 
        MC = MC - 2 
     ENDDO
  ELSE 
     DATA (2) = DATA (1) - NYQUIST (1) 
     DATA (1) = DATA (1) + NYQUIST (1) 
     DO I = 2, N / 2 + 1 
        PHR = PHASE (IPH) 
        PHI = PHASE (IPH + 1) 
        IPH = IPH + 2 
        TMP0R = DATA (M) + DATA (MC) 
        TMP0I = DATA (M + 1) - DATA (MC + 1) 
        TMP1R = PHI * (DATA (M) - DATA (MC) ) + PHR * (DATA (M + 1)    &
             + DATA (MC + 1) )                                              
        TMP1I = - PHR * (DATA (M) - DATA (MC) ) + PHI * (DATA (M + 1)  &
             + DATA (MC + 1) )                                              
        DATA (M) = TMP0R - TMP1R 
        DATA (M + 1) = TMP0I - TMP1I 
        DATA (MC) = TMP0R + TMP1R 
        DATA (MC + 1) = - (TMP0I + TMP1I) 
        M = M + 2 
        MC = MC - 2 
     ENDDO
  ENDIF

  RETURN 
END SUBROUTINE FIXREAL
