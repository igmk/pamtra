!+ <A one line description of this module> 
! 
module rt_utilities 
 
! 
! Description: 
!   <Say what this module is for> 
! 
! Current Code Owner: <Name of person responsible for this code> 
! 
! History: 
!  
! Version   Date     Comment 
! -------   ----     ------- 
! <version> <date>   Original code. <Your name> 
! 
! Code Description: 
!   Language:		Fortran 90. 
!   Software Standards: "European Standards for Writing and  
!     Documenting Exchangeable Fortran 90 Code". 
! 
! Modules used: 
! 
use kinds, only : dbl, &
long
! Imported Type Definitions: 
 
! Imported Parameters: 
 
! Imported Scalar Variables with intent (in): 
 
! Imported Scalar Variables with intent (out): 
 
! Imported Array Variables with intent (in): 
 
! Imported Array Variables with intent (out): 
 
! Imported Routines: 
 
! <Repeat from Use for each module...> 
 
! Declarations must be of the form: 
! <type>   <VariableName>      ! Description/ purpose of variable 
 
 implicit none 
! Global (i.e. public) Declarations: 
! Global Type Definitions: 
 
! Global Parameters: 
 
! Global Scalars: 
 
! Global Arrays: 
 
! Local (i.e. private) Declarations: 
! Local Type Definitions: 
 
! Local Parameters: 
 
! Local Scalars: 
 
! Local Arrays: 
 
! Operator definitions: 
!   Define new operators or overload existing ones. 
 
contains 
! Define procedures contained in this module. 
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

SUBROUTINE GAUSS_LEGENDRE_QUADRATURE(NUM, ABSCISSAS, WEIGHTS)
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
100  CONTINUE
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
100  CONTINUE
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


SUBROUTINE LOBATTO_QUADRATURE(NUM, ABSCISSAS, WEIGHTS)
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
100  CONTINUE
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
  REAL(kind=dbl) COEF (6, *), X, PHASE_MATRIX (4, 4)
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

SUBROUTINE FOURIER_MATRIX (AZIORDER, NUMPTS, NSTOKES, REAL_MATRIX,&
     BASIS_MATRIX)
  use kinds
  INTEGER AZIORDER, NUMPTS, NSTOKES
  REAL(kind=dbl) REAL_MATRIX (4, 4, *), BASIS_MATRIX (4, 4, *)
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

SUBROUTINE CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUT, WAVELEN, &
     FLUXCODE, OUTPUT)
  !       Converts the output radiance or flux arrays to VH polarization
  !     and effective blackbody temperature if desired.  OUTPOL='VH'
  !     converts the polarization basis of the first two Stokes parameters
  !     to vertical/horizontal polarization.  If UNITS='T' the radiance is
  !     converted to effective blackbody brightness temperature, and if
  !     UNITS='R' the radiance is converted to Rayleigh-Jeans brightness
  !     temperature.  If the output is flux then FLUXCODE=1, and the flux
  !     is divided by pi before converting to brightness temperature.
  use kinds

  INTEGER NSTOKES, NOUT, FLUXCODE
  REAL(kind=dbl) WAVELEN, OUTPUT (NSTOKES, NOUT)
  CHARACTER UNITS * 1, OUTPOL * 2
  INTEGER I, J
  REAL(kind=dbl) IV, IH, RAD, TEMP

  DO J = 1, NOUT
     !           Convert to Vertical and Horizontal polarization if desired
     IF (OUTPOL.EQ.'VH') THEN
        IV = 0.5 * (OUTPUT (1, J) + OUTPUT (2, J) )
        IH = 0.5 * (OUTPUT (1, J) - OUTPUT (2, J) )
        OUTPUT (1, J) = IV
        OUTPUT (2, J) = IH
     ENDIF
     !           Convert to brightness temperature
     IF (UNITS.EQ.'T'.OR.UNITS.EQ.'R') THEN
        DO I = 1, NSTOKES
           RAD = OUTPUT (I, J)
           IF (OUTPOL.EQ.'VH'.AND.I.LE.2) RAD = 2.0 * RAD
           IF (FLUXCODE.EQ.1) RAD = RAD / ACOS ( - 1.0)
           IF (UNITS.EQ.'R') THEN
              TEMP = RAD * WAVELEN**4 * 1.4388D4 / 1.1911D8
           ELSE
              IF (RAD.GT.0.0) THEN
                 TEMP = 1.4388D4 / (WAVELEN * DLOG (1.0 + 1.1911D8 /      &
                      (RAD * WAVELEN**5) ) )
              ELSEIF (RAD.EQ.0.0) THEN
                 TEMP = 0.0D0
              ELSE
                 TEMP = - 1.4388D4 / (WAVELEN * DLOG (1.0 + 1.1911D8 /    &
                      ( - RAD * WAVELEN**5) ) )
              ENDIF
           ENDIF
           OUTPUT (I, J) = TEMP
        ENDDO
     ENDIF
  ENDDO
  RETURN
END SUBROUTINE CONVERT_OUTPUT
 
end module rt_utilities 
 
!- End of module header
