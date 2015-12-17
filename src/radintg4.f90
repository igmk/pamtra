SUBROUTINE INITIALIZE4(NSTOKES, NUMMU,&
     DELTA_Z, MU_VALUES, QUAD_WEIGHTS,&
     GAS_EXTINCT, EXTINCT_MATRIX,&
     SCATTER_MATRIX,  REFLECT, TRANS)
  INTEGER NSTOKES, NUMMU
  REAL*8  DELTA_Z, GAS_EXTINCT
  REAL*8  MU_VALUES(NUMMU), QUAD_WEIGHTS(NUMMU)
  REAL*8  EXTINCT_MATRIX(nstokes,nstokes,NUMMU,2)
  REAL*8  SCATTER_MATRIX(NSTOKES,NUMMU,NSTOKES,NUMMU,4)
  REAL*8  REFLECT(NSTOKES,NUMMU,NSTOKES,NUMMU,2)
  REAL*8  TRANS(NSTOKES,NUMMU,NSTOKES,NUMMU,2)
  INTEGER  I1, J1, I2, J2
  REAL*8   C, TMP, GEXT, EXT, DIAG

  C = 2.0D0 *3.1415926535897932384D0

  DO I2 = 1, NSTOKES
     DO J2 = 1, NUMMU
        TMP = DELTA_Z/MU_VALUES(J2)
        DO I1 = 1, NSTOKES
           GEXT = 0.0
           IF (I1 .EQ. I2)  GEXT = GAS_EXTINCT
           DO J1 = 1, NUMMU
              REFLECT(I2,J2,I1,J1,1) = C*TMP*QUAD_WEIGHTS(J1)&
                   *SCATTER_MATRIX(I2,J2,I1,J1,2)
              REFLECT(I2,J2,I1,J1,2) = C*TMP*QUAD_WEIGHTS(J1)&
                   *SCATTER_MATRIX(I2,J2,I1,J1,3)
              DIAG = 0.0D0
              IF (I1 .EQ. I2 .AND. J1 .EQ. J2)  DIAG = 1.0D0
              EXT = 0.0D0
              IF (J1 .EQ. J2) EXT = EXTINCT_MATRIX(I2,I1,J2,1)+GEXT
              TRANS(I2,J2,I1,J1,1) = DIAG&
                   - TMP*( EXT - C*QUAD_WEIGHTS(J1)&
                   *SCATTER_MATRIX(I2,J2,I1,J1,1) )
              IF (J1 .EQ. J2) EXT = EXTINCT_MATRIX(I2,I1,J2,2)+GEXT
              TRANS(I2,J2,I1,J1,2) = DIAG&
                   - TMP*( EXT - C*QUAD_WEIGHTS(J1)&
                   *SCATTER_MATRIX(I2,J2,I1,J1,4) )
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE INITIALIZE4





SUBROUTINE INITIAL_SOURCE4(NSTOKES, NUMMU,&
     DELTA_Z, MU_VALUES,&
     PLANCK, EMIS_VECTOR, GAS_EXTINCT,&
     SOURCE)
  INTEGER NSTOKES, NUMMU
  REAL*8  MU_VALUES(NUMMU), DELTA_Z, PLANCK, GAS_EXTINCT
  REAL*8  EMIS_VECTOR(nstokes,NUMMU,2)
  REAL*8  SOURCE(NSTOKES,NUMMU,2)
  INTEGER  I, J, N
  REAL*8   TMP, EXT

  N = NSTOKES*NUMMU
  CALL MZERO (2*N, 1, SOURCE)

  DO I = 1, NSTOKES
     DO J = 1, NUMMU
        EXT = 0.0D0
        IF (I .EQ. 1)  EXT = GAS_EXTINCT
        TMP = PLANCK*DELTA_Z/MU_VALUES(J)
        SOURCE(I,J,1) = TMP*(EMIS_VECTOR(I,J,1) + EXT)
        SOURCE(I,J,2) = TMP*(EMIS_VECTOR(I,J,2) + EXT)
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE INITIAL_SOURCE4





SUBROUTINE NONSCATTER_LAYER4(NSTOKES, NUMMU,&
     DELTATAU, MU_VALUES,&
     PLANCK0, PLANCK1,&
     REFLECT, TRANS, SOURCE)
  !        NONSCATTER_LAYER calculates the reflection and transmission
  !      matrices and the source vectors for a purely absorbing layer.
  !      The source function, which is the Planck function, is assumed
  !      to vary linearly with optical depth across the layer.
  INTEGER  NSTOKES, NUMMU
  REAL*8   DELTATAU, MU_VALUES(NUMMU), PLANCK0, PLANCK1
  REAL*8   REFLECT(NSTOKES,NUMMU,NSTOKES,NUMMU,2)
  REAL*8   TRANS(NSTOKES,NUMMU,NSTOKES,NUMMU,2)
  REAL*8   SOURCE(NSTOKES,NUMMU,2)
  INTEGER  N, I, J
  REAL*8   FACTOR, SLOPE, PATH

  N = NSTOKES*NUMMU
  CALL MZERO (2*N, N, REFLECT)

  CALL MZERO (2*N, N, TRANS)
  DO J = 1, NUMMU
     FACTOR = DEXP(-DELTATAU/MU_VALUES(J))
     DO I = 1, NSTOKES
        TRANS(I,J,I,J,1) = FACTOR
        TRANS(I,J,I,J,2) = FACTOR
     ENDDO
  ENDDO

  CALL MZERO (2*N, 1, SOURCE)
  IF (DELTATAU .GT. 0.0) THEN
     DO J = 1, NUMMU
        PATH = DELTATAU/MU_VALUES(J)
        SLOPE = (PLANCK1-PLANCK0)/PATH
        SOURCE(1,J,1) = PLANCK1 - SLOPE&
             - ( PLANCK1 - SLOPE*(1.0+PATH) )*DEXP(-PATH)
        SOURCE(1,J,2) = PLANCK0 + SLOPE&
             - ( PLANCK0 + SLOPE*(1.0+PATH) )*DEXP(-PATH)
     ENDDO
  ENDIF
  RETURN
END SUBROUTINE NONSCATTER_LAYER4

SUBROUTINE INTERNAL_RADIANCE4(N, UPREFLECT, UPTRANS, UPSOURCE,&
     DOWNREFLECT, DOWNTRANS, DOWNSOURCE, &
     INTOPRAD, INBOTTOMRAD,&
     UPRAD, DOWNRAD)
  !        INTERNAL_RADIANCE calculates the internal radiance at a level.
  !      The reflection and transmission matrices and source vector are
  !      given for the atmosphere above (UP) and below (DOWN) the level.
  !      The upwelling and downwelling radiance are computed from the
  !      two layer properties and the radiance incident on the top and bottom.
  INTEGER   N
  REAL*8    UPREFLECT(N,N,2), DOWNREFLECT(N,N,2)
  REAL*8    UPTRANS(N,N,2), DOWNTRANS(N,N,2)
  REAL*8    UPSOURCE(N,2), DOWNSOURCE(N,2)
  REAL*8    INTOPRAD(N), INBOTTOMRAD(N)
  REAL*8    UPRAD(N), DOWNRAD(N)
  INTEGER   MAXV, MAXM
  PARAMETER (MAXV=64, MAXM=4096)
  REAL*8    S(MAXV), V(MAXV)
  REAL*8    X(MAXM), Y(MAXM)
  COMMON /SCRATCH1/ X, Y


  !               Compute gamma plus
  CALL MMULT (N,N,N, UPREFLECT(1,1,1), DOWNREFLECT(1,1,2), X)
  CALL MIDENTITY (N, Y)
  CALL MSUB (N,N, Y, X, Y)
  CALL MINVERT (N, Y, X)
  !               Calculate the internal downwelling (plus) radiance vector
  CALL MMULT (N,N,1, DOWNTRANS(1,1,2), INBOTTOMRAD, V)
  CALL MMULT (N,N,1, UPREFLECT(1,1,1), V, S)
  CALL MMULT (N,N,1, UPTRANS(1,1,1), INTOPRAD, V)
  CALL MADD (N, 1, V, S, S)
  CALL MMULT (N,N,1, UPREFLECT(1,1,1), DOWNSOURCE(1,2), V)
  CALL MADD (N, 1, V, S, S)
  CALL MADD (N, 1, UPSOURCE(1,1), S, S)
  CALL MMULT (N,N,1, X, S, DOWNRAD)

  !               Compute gamma minus
  CALL MMULT (N,N,N, DOWNREFLECT(1,1,2), UPREFLECT(1,1,1), X)
  CALL MIDENTITY (N, Y)
  CALL MSUB (N,N, Y, X, Y)
  CALL MINVERT (N, Y, X)
  !               Calculate the internal upwelling (minus) radiance vector
  CALL MMULT (N,N,1, UPTRANS(1,1,1), INTOPRAD, V)
  CALL MMULT (N,N,1, DOWNREFLECT(1,1,2), V, S)
  CALL MMULT (N,N,1, DOWNTRANS(1,1,2), INBOTTOMRAD, V)
  CALL MADD (N, 1, V, S, S)
  CALL MMULT (N,N,1, DOWNREFLECT(1,1,2), UPSOURCE(1,1), V)
  CALL MADD (N, 1, V, S, S)
  CALL MADD (N, 1, DOWNSOURCE(1,2), S, S)
  CALL MMULT (N,N,1, X, S, UPRAD)

  RETURN
END SUBROUTINE INTERNAL_RADIANCE4

SUBROUTINE DOUBLING_INTEGRATION4 (N, NUM_DOUBLES,      &
     SYMMETRIC, REFLECT, TRANS, LIN_SOURCE,     &
     LINFACTOR, T_REFLECT, T_TRANS, T_SOURCE)
  !        DOUBLING_INTEGRATION integrates homogeneous thin layers using
  !      the doubling algorithm.  NUM_DOUBLES doubling steps are done.
  !      The initial reflection (REFLECT) and transmission (TRANS) matrices
  !      are input.  Depending on SRC_CODE linear (thermal) and exponential
  !      (solar) sources are doubled.  The EXP_SOURCE and LIN_SOURCE vectors
  !      are the source vectors at zero optical depth.  The EXPFACTOR is
  !      the single layer attenuation factor for the exponential source,
  !      while the LINFACTOR is the single layer slope for the linear source.
  !      The SYMMETRIC flag specifies whether the plus and minus parts
  !      of the reflection and transmission matrices are separately calculated
  !      or are assumed to be the same.  The integrated output is in
  !      T_REFLECT, T_TRANS, and T_SOURCE.

  use kinds
  INTEGER N, NUM_DOUBLES
  LOGICAL SYMMETRIC
  REAL(kind=dbl) REFLECT (N, N, 2), TRANS (N, N, 2)
  REAL(kind=dbl) LIN_SOURCE (N, 2), LINFACTOR
  REAL(kind=dbl) T_REFLECT (N, N, 2), T_TRANS (N, N, 2)
  REAL(kind=dbl) T_SOURCE (N, 2)
  REAL(kind=dbl) Y_N1(N,1), X_N1(N,1)
  INTEGER MAXV, MAXM
  PARAMETER (MAXV = 64, MAXM = 4096)
  INTEGER I, NM
  REAL(kind=dbl) T_LIN (2 * MAXV)
  REAL(kind=dbl) CONST (2 * MAXV), T_CONST (2 * MAXV)
  REAL(kind=dbl) LINFAC
  REAL(kind=dbl) Y_NN(N,N), X_NN(N,N), G_NN(N,N), ID_NN(N,N)

  LINFAC = LINFACTOR
  CALL MCOPY (N, 1, LIN_SOURCE (1, 1), CONST (1) )
  CALL MCOPY (N, 1, LIN_SOURCE (1, 2), CONST (1 + N) )
  CONST (1:N)=LIN_SOURCE (:, 1)
  CONST (N+1:2*N)=LIN_SOURCE (:, 2)

  DO I = 1, NUM_DOUBLES

     !           Make gamma plus matrix: GAMMA = inv[1 - Rp*Rm]
     X_NN=MATMUL(REFLECT (:, :, 1), REFLECT (:, :, 2))
     CALL MIDENTITY (N, ID_NN)
     Y_NN=ID_NN - X_NN
     CALL MINVERT (N, Y_NN, G_NN)

     !           Rp(2N) = Rp + Tp * GAMMA * Rp * Tm
     X_NN=MATMUL(REFLECT (:, :, 1), TRANS (:, :, 2))
     Y_NN=MATMUL(G_NN, X_NN)
     X_NN=MATMUL(TRANS (:, :, 1), Y_NN)
     T_REFLECT (:, :, 1)=REFLECT (:, :, 1) + X_NN

     !           Tp(2N) = Tp * GAMMA * Tp
     X_NN=MATMUL(G_NN, TRANS (:, :, 1))
     T_TRANS (:, :, 1)=MATMUL(TRANS (:, :, 1), X_NN)

     !             Sp(2N) = (Sp+f*Cp) + Tp * GAMMA * (Sp + Rp * (Sm+f*Cm))
     X_N1(:,1)=LINFAC * CONST (1+N:2*N)
     Y_N1(:,1)=LIN_SOURCE (:, 2) + X_N1(:,1)
     X_N1=MATMUL(REFLECT (:, :, 1), Y_N1)
     Y_N1(:,1)=LIN_SOURCE (:, 1) + X_N1(:,1)
     X_N1=MATMUL(G_NN, Y_N1)
     Y_N1=MATMUL(TRANS (:, :, 1), X_N1)
     X_N1(:,1)=LIN_SOURCE (:, 1) + Y_N1(:,1)
     Y_N1(:,1)=LINFAC * CONST(1:N)
     T_LIN (1:N)=X_N1(:,1) + Y_N1(:,1)

     !             Cp(2N) = Cp + Tp * GAMMA * (Cp + Rp * Cm)
     Y_N1(:,1)=CONST (1+N:2*N)
     X_N1=MATMUL(REFLECT (:, :, 1),Y_N1)
     Y_N1(:,1)=CONST (1:N) + X_N1(:,1)
     X_N1=MATMUL(G_NN, Y_N1)
     Y_N1=MATMUL(TRANS (:, :, 1), X_N1)
     T_CONST (1:N)=CONST (1:N) + Y_N1(:,1)


     IF (SYMMETRIC) THEN
        CALL MCOPY (N, N, T_REFLECT (1, 1, 1), T_REFLECT (1, 1, 2) )
        CALL MCOPY (N, N, T_TRANS (1, 1, 1), T_TRANS (1, 1, 2) )

     ELSE
        !             Make gamma minus matrix: GAMMA = inv[1 - Rm*Rp]
        X_NN=MATMUL( REFLECT (:, :, 2), REFLECT (:, :, 1))
        Y_NN=ID_NN - X_NN
        CALL MINVERT (N, Y_NN, G_NN)

        !             Rm(2N) = Rm + Tm * GAMMA * Rm * Tp
        X_NN=MATMUL(REFLECT (:, :, 2), TRANS (:, :, 1))
        Y_NN=MATMUL(G_NN, X_NN)
        X_NN=MATMUL(TRANS (:, :, 2), Y_NN)
        T_REFLECT (:, :, 2)=REFLECT (:, :, 2) + X_NN


        !             Tm(2N) = Tm * GAMMA * Tm
        X_NN=MATMUL(G_NN, TRANS (:, :, 2))
        T_TRANS (:, :, 2)=MATMUL(TRANS (:, :, 2), X_NN)
     ENDIF

     !           Linear source doubling
     !             Sm(2N) = Sm + Tm * GAMMA * (Sm+f*Cm + Rm * Sp)
     X_N1(:,1)=LINFAC * CONST (1+N:2*N)
     Y_N1(:,1)=LIN_SOURCE (:, 2) + X_N1(:,1)
     X_N1(:,1)=LIN_SOURCE (:, 1)
     X_N1=MATMUL(REFLECT (:, :, 2), X_N1)
     Y_N1(:,1)=Y_N1(:,1) + X_N1(:,1)
     X_N1=MATMUL(G_NN, Y_N1)
     Y_N1=MATMUL(TRANS (:, :, 2), X_N1)
     T_LIN (1+N:2*N)=LIN_SOURCE (:, 2) + Y_N1(:,1)

     !             Cm(2N) = Cm + Tm * GAMMA * (Cm + Rm * Cp)
     X_N1(:,1)=CONST (1:N)
     X_N1=MATMUL(REFLECT (:, :, 2), X_N1)
     Y_N1(:,1)=CONST (1+N:2*N) + X_N1(:,1)
     X_N1=MATMUL(G_NN, Y_N1)
     Y_N1=MATMUL(TRANS (:, :, 2), X_N1)
     T_CONST (1+N:2*N)=CONST (1+N:2*n) + Y_N1(:,1)

     LIN_SOURCE (:, 1) = T_LIN (1:N)
     LIN_SOURCE (:, 2) = T_LIN (N+1:2*N)
     CONST = T_CONST
     LINFAC = 2.0 * LINFAC

     REFLECT = T_REFLECT
     TRANS = T_TRANS
  ENDDO

  NM = 2 * N
  IF (NUM_DOUBLES.LE.0) THEN
     T_REFLECT = REFLECT
     T_TRANS = TRANS
  ENDIF

  T_SOURCE = LIN_SOURCE

  RETURN
END SUBROUTINE DOUBLING_INTEGRATION4

!      SUBROUTINE DOUBLING_INTEGRATION4(N, NUM_DOUBLES,&
!                            SYMMETRIC, REFLECT, TRANS,&
!                            LIN_SOURCE, LINFACTOR,&
!                            T_REFLECT, T_TRANS, T_SOURCE)
! !        DOUBLING_INTEGRATION integrates homogeneous thin layers using
! !      the doubling algorithm.  NUM_DOUBLES doubling steps are done.
! !      The initial reflection (REFLECT) and transmission (TRANS) matrices
! !      are input.  A linear (thermal) source is assumed.
! !      The LIN_SOURCE vector is the source vector at zero optical depth.
! !      The LINFACTOR is the single layer slope for the linear source.
! !      The SYMMETRIC flag specifies whether the plus and minus parts
! !      of the reflection and transmission matrices are separately calculated
! !      or are assumed to be the same.  The integrated output is in
! !      T_REFLECT, T_TRANS, and T_SOURCE.
!      INTEGER  N, NUM_DOUBLES
!      LOGICAL  SYMMETRIC
!      REAL*8   REFLECT(N,N,2), TRANS(N,N,2)
!      REAL*8   LIN_SOURCE(N,2), LINFACTOR
!      REAL*8   T_REFLECT(N,N,2), T_TRANS(N,N,2)
!      REAL*8   T_SOURCE(N,2)
!      INTEGER  MAXV, MAXM
!      PARAMETER (MAXV=64, MAXM=4096)
!      INTEGER  I, NM
!      REAL*8   T_LIN(2*MAXV)
!      REAL*8   CONST(2*MAXV), T_CONST(2*MAXV)
!      REAL*8   LINFAC, ZERO
!      REAL*8   X(MAXM), Y(MAXM)
!      REAL*8   GAMMA(MAXM)
!      COMMON /SCRATCH1/ X, Y
!      COMMON /SCRATCH2/ GAMMA
!      PARAMETER (ZERO=0.0D0)
!
!
!
!      LINFAC = LINFACTOR
!      CALL MCOPY (N, 1, LIN_SOURCE(1,1), CONST(1))
!      CALL MCOPY (N, 1, LIN_SOURCE(1,2), CONST(1+N))
!
!      DO I = 1, NUM_DOUBLES
!
! !           Make gamma plus matrix: GAMMA = inv[1 - Rp*Rm]
!          CALL MMULT (N, N, N, REFLECT(1,1,1), REFLECT(1,1,2), X)
!          CALL MIDENTITY (N, Y)
!          CALL MSUB (N, N, Y, X, Y)
!          CALL MINVERT (N, Y, GAMMA)
!
! !           Rp(2N) = Rp + Tp * GAMMA * Rp * Tm
!          CALL MMULT (N, N, N, REFLECT(1,1,1), TRANS(1,1,2), X)
!          CALL MMULT (N, N, N, GAMMA, X, Y)
!          CALL MMULT (N, N, N, TRANS(1,1,1), Y, X)
!          CALL MADD (N, N, REFLECT(1,1,1), X, T_REFLECT(1,1,1))
!
! !           Tp(2N) = Tp * GAMMA * Tp
!          CALL MMULT (N, N, N, GAMMA, TRANS(1,1,1), X)
!          CALL MMULT (N, N, N, TRANS(1,1,1), X, T_TRANS(1,1,1))
!
! !           Linear source doubling
! !             Sp(2N) = (Sp+f*Cp) + Tp * GAMMA * (Sp + Rp * (Sm+f*Cm))
!            CALL MSCALARMULT (N, 1, LINFAC, CONST(1+N), X)
!            CALL MADD (N, 1, LIN_SOURCE(1,2), X, Y)
!            CALL MMULT (N, N, 1, REFLECT(1,1,1), Y, X)
!            CALL MADD (N, 1, LIN_SOURCE(1,1), X, Y)
!            CALL MMULT (N, N, 1, GAMMA, Y, X)
!            CALL MMULT (N, N, 1, TRANS(1,1,1), X, Y)
!            CALL MADD (N, 1, LIN_SOURCE(1,1), Y, X)
!            CALL MSCALARMULT (N, 1, LINFAC, CONST(1), Y)
!            CALL MADD (N, 1, X, Y, T_LIN(1))
! !             Cp(2N) = Cp + Tp * GAMMA * (Cp + Rp * Cm)
!            CALL MMULT (N, N, 1, REFLECT(1,1,1), CONST(1+N), X)
!            CALL MADD (N, 1, CONST(1), X, Y)
!            CALL MMULT (N, N, 1, GAMMA, Y, X)
!            CALL MMULT (N, N, 1, TRANS(1,1,1), X, Y)
!            CALL MADD (N, 1, CONST(1), Y, T_CONST(1))
!
!
!          IF (SYMMETRIC) THEN
!            CALL MCOPY (N, N, T_REFLECT(1,1,1), T_REFLECT(1,1,2))
!            CALL MCOPY (N, N, T_TRANS(1,1,1), T_TRANS(1,1,2))
!
!          ELSE
! !             Make gamma minus matrix: GAMMA = inv[1 - Rm*Rp]
!            CALL MMULT (N,N,N, REFLECT(1,1,2), REFLECT(1,1,1), X)
!            CALL MIDENTITY (N, Y)
!            CALL MSUB (N, N, Y, X, Y)
!            CALL MINVERT (N, Y, GAMMA)
!
! !             Rm(2N) = Rm + Tm * GAMMA * Rm * Tp
!            CALL MMULT (N, N, N, REFLECT(1,1,2), TRANS(1,1,1), X)
!            CALL MMULT (N, N, N, GAMMA, X, Y)
!            CALL MMULT (N, N, N, TRANS(1,1,2), Y, X)
!            CALL MADD (N, N, REFLECT(1,1,2), X, T_REFLECT(1,1,2))
!
! !             Tm(2N) = Tm * GAMMA * Tm
!            CALL MMULT (N, N, N, GAMMA, TRANS(1,1,2), X)
!            CALL MMULT (N, N, N, TRANS(1,1,2), X, T_TRANS(1,1,2))
!          ENDIF
!
!
! !           Linear source doubling
! !             Sm(2N) = Sm + Tm * GAMMA * (Sm+f*Cm + Rm * Sp)
!            CALL MSCALARMULT (N, 1, LINFAC, CONST(1+N), X)
!            CALL MADD (N, 1, LIN_SOURCE(1,2), X, Y)
!            CALL MMULT (N,N,1, REFLECT(1,1,2), LIN_SOURCE(1,1), X)
!            CALL MADD (N, 1, Y, X, Y)
!            CALL MMULT (N, N, 1, GAMMA, Y, X)
!            CALL MMULT (N, N, 1, TRANS(1,1,2), X, Y)
!            CALL MADD (N, 1, LIN_SOURCE(1,2), Y, T_LIN(1+N))
! !             Cm(2N) = Cm + Tm * GAMMA * (Cm + Rm * Cp)
!            CALL MMULT (N, N, 1, REFLECT(1,1,2), CONST(1), X)
!            CALL MADD (N, 1, CONST(1+N), X, Y)
!            CALL MMULT (N, N, 1, GAMMA, Y, X)
!            CALL MMULT (N, N, 1, TRANS(1,1,2), X, Y)
!            CALL MADD (N, 1, CONST(1+N), Y, T_CONST(1+N))
!
!            CALL MCOPY (N,1, T_LIN(1), LIN_SOURCE(1,1))
!            CALL MCOPY (N,1, T_LIN(1+N), LIN_SOURCE(1,2))
!            CALL MCOPY (2*N,1, T_CONST, CONST)
!            LINFAC = 2.0*LINFAC
!
!          CALL MCOPY (N,N, T_REFLECT(1,1,1), REFLECT(1,1,1))
!          CALL MCOPY (N,N, T_REFLECT(1,1,2), REFLECT(1,1,2))
!          CALL MCOPY (N,N, T_TRANS(1,1,1), TRANS(1,1,1))
!          CALL MCOPY (N,N, T_TRANS(1,1,2), TRANS(1,1,2))
!      ENDDO
!
!      NM = 2*N
!      IF (NUM_DOUBLES .LE. 0) THEN
!          CALL MCOPY (NM,N, REFLECT, T_REFLECT)
!          CALL MCOPY (NM,N, TRANS, T_TRANS)
!      ENDIF
!
!      CALL MCOPY (NM,1, LIN_SOURCE, T_SOURCE)
!
!      RETURN
!      END

SUBROUTINE COMBINE_LAYERS4(N, REFLECT1, TRANS1, SOURCE1, REFLECT2,&
     TRANS2, SOURCE2, OUT_REFLECT, OUT_TRANS, OUT_SOURCE)
  !        COMBINE_LAYERS combines the reflection and transmission matrice
  !      and source vectors for two layers into a combined reflection,
  !      transmission and source.  The positive side (down) of the first
  !      layer is attached to the negative side (up) of the second layer;
  !      thus layer 1 is put on top of layer 2.
  use kinds
  IMPLICIT NONE
  INTEGER N
  REAL(kind=dbl) REFLECT1 (N, N, 2), TRANS1 (N, N, 2)
  REAL(kind=dbl) REFLECT2 (N, N, 2), TRANS2 (N, N, 2)
  REAL(kind=dbl) SOURCE1 (N, 2), SOURCE2 (N, 2)
  REAL(kind=dbl) OUT_REFLECT (N, N, 2), OUT_TRANS (N, N, 2)
  REAL(kind=dbl) OUT_SOURCE (N, 2)
  REAL(kind=dbl) Y_NN(N,N), X_NN(N,N), G_NN(N,N), ID_NN(N,N), WORK_N1(N)
  INTEGER MAXM
  PARAMETER (MAXM = 4096)

  !           GAMMAp = inv[1 - R1p * R2m]     (p for +,  m for -)
  X_NN=MATMUL(REFLECT1 (:, :, 1), REFLECT2 (:, :, 2))
  CALL MIDENTITY (N, ID_NN)
  Y_NN=ID_NN - X_NN

  CALL MINVERT (N, Y_NN, G_NN)

  !           RTp = R2p + T2p * GAMMAp * R1p * T2m
  X_NN=MATMUL(REFLECT1 (:, :, 1), TRANS2 (:, :, 2))
  Y_NN=MATMUL(G_NN, X_NN)
  X_NN=MATMUL(TRANS2 (:, :, 1), Y_NN)
  OUT_REFLECT(:, :, 1)=REFLECT2 (:, :, 1) + X_NN

  !           TTp = T2p * GAMMAp * T1p
  X_NN=MATMUL(G_NN, TRANS1(:, :, 1))
  OUT_TRANS (:, :, 1) = MATMUL(TRANS2 (:, :, 1), X_NN)

  !           STp = S2p + T2p * GAMMAp * (S1p + R1p * S2m)
  WORK_N1=MATMUL(REFLECT1 (:, :, 1), SOURCE2 (:, 2))
  WORK_N1=SOURCE1(:,1)+ WORK_N1
  WORK_N1=MATMUL(G_NN,WORK_N1)
  WORK_N1=MATMUL(TRANS2 (:, :, 1), WORK_N1)
  OUT_SOURCE (:, 1)=SOURCE2 (:, 1) + WORK_N1

  !           GAMMAm = inv[1 - R2m * R1p]
  X_NN=MATMUL(REFLECT2 (:, :, 2), REFLECT1 (:, :, 1))
  Y_NN= ID_NN - X_NN
  CALL MINVERT (N, Y_NN, G_NN)

  !           RTm = R1m + T1m * GAMMAm * R2m * T1p
  X_NN=MATMUL(REFLECT2 (:, :, 2), TRANS1 (:, :, 1))
  Y_NN=MATMUL(G_NN, X_NN)
  X_NN=MATMUL(TRANS1 (:, :, 2), Y_NN)
  OUT_REFLECT (:, :, 2)=REFLECT1 (:, :, 2) + X_NN

  !           TTm = T1m * GAMMAm * T2m
  X_NN=MATMUL(G_NN, TRANS2 (:, :, 2))
  OUT_TRANS (:, :, 2)=MATMUL(TRANS1 (:, :, 2), X_NN)

  !           STm = S1m + T1m * GAMMAm * (S2m + R2m * S1p)
  WORK_N1=MATMUL(REFLECT2 (:, :, 2), SOURCE1 (:, 1))
  WORK_N1=SOURCE2 (:, 2) + WORK_N1
  WORK_N1=MATMUL(G_NN, WORK_N1)
  WORK_N1=MATMUL(TRANS1 (:, :, 2), WORK_N1)
  OUT_SOURCE (:, 2)=SOURCE1 (:, 2) + WORK_N1

  RETURN
END SUBROUTINE COMBINE_LAYERS4




!      SUBROUTINE COMBINE_LAYERS4(N,&
!                             REFLECT1, TRANS1, SOURCE1,&
!                             REFLECT2, TRANS2, SOURCE2,&
!                             OUT_REFLECT, OUT_TRANS, OUT_SOURCE)
!!        COMBINE_LAYERS combines the reflection and transmission matrices
!!      and source vectors for two layers into a combined reflection,
!!      transmission and source.  The positive side (down) of the first
!!      layer is attached to the negative side (up) of the second layer;
!!      thus layer 1 is put on top of layer 2.
!      INTEGER  N
!      REAL*8   REFLECT1(N,N,2), TRANS1(N,N,2)
!      REAL*8   REFLECT2(N,N,2), TRANS2(N,N,2)
!      REAL*8   SOURCE1(N,2), SOURCE2(N,2)
!      REAL*8   OUT_REFLECT(N,N,2), OUT_TRANS(N,N,2)
!      REAL*8   OUT_SOURCE(N,2)
!      INTEGER  MAXM
!      PARAMETER (MAXM=4096)
!      REAL*8   X(MAXM), Y(MAXM)
!      REAL*8   GAMMA(MAXM)
!      COMMON /SCRATCH1/ X, Y
!      COMMON /SCRATCH2/ GAMMA
!
!!           GAMMAp = inv[1 - R1p * R2m]     (p for +,  m for -)
!      CALL MMULT (N, N, N, REFLECT1(1,1,1), REFLECT2(1,1,2), X)
!      CALL MIDENTITY (N, Y)
!      CALL MSUB (N, N, Y, X, Y)
!      CALL MINVERT (N, Y, GAMMA)
!
!!           RTp = R2p + T2p * GAMMAp * R1p * T2m
!      CALL MMULT (N, N, N, REFLECT1(1,1,1), TRANS2(1,1,2), X)
!      CALL MMULT (N, N, N, GAMMA, X, Y)
!      CALL MMULT (N, N, N, TRANS2(1,1,1), Y, X)
!      CALL MADD (N, N, REFLECT2(1,1,1), X, OUT_REFLECT(1,1,1))
!
!!           TTp = T2p * GAMMAp * T1p
!      CALL MMULT (N, N, N, GAMMA, TRANS1(1,1,1), X)
!      CALL MMULT (N, N, N, TRANS2(1,1,1), X, OUT_TRANS(1,1,1))
!
!!           STp = S2p + T2p * GAMMAp * (S1p + R1p * S2m)
!      CALL MMULT (N, N, 1, REFLECT1(1,1,1), SOURCE2(1,2), X)
!      CALL MADD (N, 1, SOURCE1(1,1), X, Y)
!      CALL MMULT (N, N, 1, GAMMA, Y, X)
!      CALL MMULT (N, N, 1, TRANS2(1,1,1), X, Y)
!      CALL MADD (N, 1, SOURCE2(1,1), Y, OUT_SOURCE(1,1))
!
!!           GAMMAm = inv[1 - R2m * R1p]
!      CALL MMULT (N,N,N, REFLECT2(1,1,2), REFLECT1(1,1,1), X)
!      CALL MIDENTITY (N, Y)
!      CALL MSUB (N, N, Y, X, Y)
!      CALL MINVERT (N, Y, GAMMA)
!
!!           RTm = R1m + T1m * GAMMAm * R2m * T1p
!      CALL MMULT (N, N, N, REFLECT2(1,1,2), TRANS1(1,1,1), X)
!      CALL MMULT (N, N, N, GAMMA, X, Y)
!      CALL MMULT (N, N, N, TRANS1(1,1,2), Y, X)
!      CALL MADD (N, N, REFLECT1(1,1,2), X, OUT_REFLECT(1,1,2))
!
!!           TTm = T1m * GAMMAm * T2m
!      CALL MMULT (N, N, N, GAMMA, TRANS2(1,1,2), X)
!      CALL MMULT (N, N, N, TRANS1(1,1,2), X, OUT_TRANS(1,1,2))
!
!!           STm = S1m + T1m * GAMMAm * (S2m + R2m * S1p)
!      CALL MMULT (N, N, 1, REFLECT2(1,1,2), SOURCE1(1,1), X)
!      CALL MADD (N, 1, SOURCE2(1,2), X, Y)
!      CALL MMULT (N, N, 1, GAMMA, Y, X)
!      CALL MMULT (N, N, 1, TRANS1(1,1,2), X, Y)
!      CALL MADD (N, 1, SOURCE1(1,2), Y, OUT_SOURCE(1,2))
!
!      RETURN
!      END


