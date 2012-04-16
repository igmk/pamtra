program run_trmatrix
    implicit none

    call cal_scat_rt4('S','L',16,270.0d0,0.19986163866666665d04, 2,1.d0)

end program run_trmatrix

SUBROUTINE LOBATTO_QUADRATURE(NUM, ABSCISSAS, WEIGHTS)
  !        Generates the abscissas and weights for an even 2*NUM point
  !      Gauss-Legendre quadrature.  Only the NUM positive points are retu

  INTEGER NUM
  REAL*8 ABSCISSAS ( * ), WEIGHTS ( * )
  INTEGER N, N1, I, J, L
  REAL*8 X, XP, PL, PL1, PL2, DPL, D2PL, CI, TINY
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

SUBROUTINE GAUSS_LEGENDRE_QUADRATURE(NUM, ABSCISSAS, WEIGHTS)
  !      Generates the abscissas and weights for an even 2*NUM point
  !      Gauss-Legendre quadrature.  Only the NUM positive points are returned

  INTEGER NUM
  REAL*8 ABSCISSAS (1), WEIGHTS (1)
  INTEGER N, I, J, L
  REAL*8 X, XP, PL, PL1, PL2, DPL, TINY
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
