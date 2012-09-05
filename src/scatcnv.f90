      subroutine scatcnv(out_file,nlegen,&
      coef,extinct,albedo,scatter_matrix,ext_matrix,emis_vec)
!       Converts a scattering file (e.g. from MIESCAT) to DDA radiative
!       transfer format file.  Scattering files have the six unique elements
!       of the Stokes scattering matrix for randomly oriented particles in
!       the form of a Legendre series in cosine of the scattering angle.
!       DDA radiative transfer files have the scattering matrix, extinction,
!       matrix, and emission vector listed and ready for use in a discrete
!       angle radiative transfer code (i.e. RT4).  The Stokes 4x4
!       scattering matrix is listed for the incident and outgoing quadrature
!       angles and for the Fourier modes in azimuth.  The extinction matrix
!       and emission vector are listed for the quadrature incident angles.
!       The guts of the code are taken from RADSCAT3.FOR.
!

      use nml_params, only: dump_to_file, nstokes,nummu, aziorder, quad_type

      implicit none

      INTEGER  NLEGEN
      REAL*8   MU_VALUES(32), QUAD_WEIGHTS(32)
      REAL*8   COEF(6,100),  EXTINCT, ALBEDO, CONST
      CHARACTER*64  SCAT_FILE, OUT_FILE
      CHARACTER*8   QUADTYPE
      real*8 scatter_matrix(nstokes,nummu,nstokes,nummu,4)
      real*8 ext_matrix(nstokes,nstokes,nummu,2)
      real*8 emis_vec(nstokes,nummu,2)

      IF (QUAD_TYPE(1:1) .EQ. 'L') THEN
        QUADTYPE = 'LOBATTO'
        CALL LOBATTO_QUADRATURE(NUMMU, MU_VALUES, QUAD_WEIGHTS)
      ELSE
        QUADTYPE = 'GAUSSIAN'
        CALL GAUSS_LEGENDRE_QUADRATURE (NUMMU, MU_VALUES, QUAD_WEIGHTS)
      ENDIF

! print *, NUMMU
! print *, "ff", (180.*acos(MU_VALUES(1:NUMMU))/3.14)
! print *, "fdd", 180-(180.*acos(MU_VALUES(NUMMU:1:-1))/3.14)
      CONST = EXTINCT*ALBEDO/ (4.0D0*3.1415926535897932384D0)

      if (dump_to_file) then
          OPEN (UNIT=2, FILE=OUT_FILE, STATUS='UNKNOWN')
          WRITE (2,11) NUMMU, AZIORDER, QUADTYPE
11        FORMAT (1X,I3,2X,I3,3X,1H',A8,1H')
          CALL SCATTERING_CNV(CONST, MU_VALUES, NLEGEN, COEF)

          CALL OUTPUT_EXTINCT_EMIS(NUMMU, MU_VALUES, EXTINCT, ALBEDO)

          CLOSE (2)
      end if

      call transform_scatter(const,mu_values, nlegen, coef, scatter_matrix)

      call transform_ext_emis(MU_VALUES, EXTINCT, ALBEDO,ext_matrix,emis_vec)

      return

      end subroutine scatcnv

      subroutine transform_ext_emis(mu_values,extinct,albedo,ext_matrix,emis_vec)

      use kinds
      use nml_params, only: nstokes, nummu
      implicit none

      integer :: j, l
      real(kind=dbl), dimension(nummu) :: mu_values
      real(kind=dbl) :: extinct, albedo
      real(kind=dbl) :: mu, absorb
      real(kind=dbl), dimension(nstokes,nstokes,nummu,2) :: ext_matrix
      real(kind=dbl), dimension(nstokes,nummu,2) :: emis_vec

      do l = 1, 2
        do j = 1, nummu
          ext_matrix(1,1:2,j,l) = (/extinct, 0.0d0/)
          ext_matrix(2,1:2,j,l) = (/0.0d0, extinct/)
!          ext_matrix(3,1:4,j,l) = (/0.0d0, 0.0d0, extinct, 0.0d0/)
!          ext_matrix(4,1:4,j,l) = (/0.0d0, 0.0d0, 0.0d0, extinct/)
        end do
      end do

      absorb = (1.0-albedo)*extinct
      do l = 1, 2
        do j = 1, nummu
          emis_vec(1:2,j,l) = (/absorb,0.0d0/)
        end do
      end do

      return

      end subroutine transform_ext_emis

      SUBROUTINE transform_scatter(CONST,MU_VALUES, NLEGEN, COEF,scatter_matrix)

      use kinds
      use nml_params, only: nummu, aziorder, nstokes

      implicit none

      INTEGER  NLEGEN, l
      REAL(kind=dbl) :: MU_VALUES(NUMMU), CONST, COEF(6,1)
      INTEGER  MAXLEG
      PARAMETER (MAXLEG=64)
      INTEGER  I1, I2, I, J1, J2, K, L1, L2, M
      INTEGER  NUMPTS, DOSUM(6)
      REAL(kind=dbl) :: C, MU1, MU2,  DELPHI, COS_SCAT
      REAL(kind=dbl) :: PHASE_MATRIX(4,4)
      REAL(kind=dbl) :: SCAT_MATRIX(4,4,4*MAXLEG), BASIS_MATRIX(4,4,4*MAXLEG)
      REAL(kind=dbl) :: ZERO, TWOPI
      PARAMETER (ZERO=0.0D0, TWOPI=2.0D0*3.1415926535897932384D0)

      real(kind=dbl) :: scatter_matrix(nstokes,nummu,nstokes,nummu,4)


      NUMPTS = 2* 2**INT(LOG(FLOAT(NLEGEN+4))/LOG(2.0d0)+1.0d0)
      IF (AZIORDER .EQ. 0)  NUMPTS = 2*INT((NLEGEN+1)/2) + 4

!           Find how many Legendre series must be summed
      CALL NUMBER_SUMS (NSTOKES, NLEGEN, COEF, DOSUM)

      DO I = 1, NLEGEN+1
        DO K = 1, 6
          COEF(K,I) = CONST*COEF(K,I)
        end do
      end do

!       MU1 is the incoming direction, and MU2 is the outgoing direction.
      DO L1 = 1, 2
        DO J1 = 1, NUMMU
          MU1 = MU_VALUES(J1)
          IF (L1 .EQ. 2)  MU1 = -MU1
          DO L2 = 1, 2
            L = 2*(L2-1)+L1
            DO J2 = 1, NUMMU
              MU2 = MU_VALUES(J2)
              IF (L2 .EQ. 2)  MU2 = -MU2
!                   Only need to calculate phase matrix for half of
!                     the delphi's, the rest come from symmetry.
              DO K = 1, NUMPTS/2 + 1
                DELPHI = (TWOPI*(K-1))/NUMPTS
                COS_SCAT = MU1*MU2 + DSQRT((1.-MU1**2)*(1.-MU2**2))*&
                                 DCOS(DELPHI)
                CALL SUM_LEGENDRE(NLEGEN, COEF, COS_SCAT,&
                                DOSUM, PHASE_MATRIX)
                CALL ROTATE_PHASE_MATRIX (PHASE_MATRIX, MU1, MU2,&
                      DELPHI, COS_SCAT, SCAT_MATRIX(1,1,K), NSTOKES)
                CALL MATRIX_SYMMETRY (NSTOKES, SCAT_MATRIX(1,1,K),&
                                   SCAT_MATRIX(1,1,NUMPTS-K+2) )
            end do
            CALL FOURIER_MATRIX (AZIORDER, NUMPTS, NSTOKES,&
                        SCAT_MATRIX, BASIS_MATRIX)
              do i2 = 1, nstokes
                do i1 = 1, nstokes
                  scatter_matrix(i2,j2,i1,j1,l) = BASIS_MATRIX(I2,I1,1)
                end do
              end do
            end do
          end do
        end do
      end do

      return

      end subroutine transform_scatter


      SUBROUTINE OUTPUT_EXTINCT_EMIS(NUMMU, MU_VALUES,EXTINCT, ALBEDO)

      implicit none

      INTEGER  NUMMU
      REAL*8   MU_VALUES(NUMMU), EXTINCT, ALBEDO
      INTEGER  J, L
      REAL*8   MU, ABSORB

      WRITE (2,*) 'C  EXTINCTION MATRIX'
      DO 100 L = 1, 2
        DO 100 J = 1, NUMMU
          MU = MU_VALUES(J)
          IF (L .EQ. 2) MU = -MU
          WRITE (2,'(1X,F11.8)') MU
          WRITE (2,'(4(2X,E15.8))') EXTINCT, 0.0, 0.0, 0.0
          WRITE (2,'(4(2X,E15.8))') 0.0, EXTINCT, 0.0, 0.0
          WRITE (2,'(4(2X,E15.8))') 0.0, 0.0, EXTINCT, 0.0
          WRITE (2,'(4(2X,E15.8))') 0.0, 0.0, 0.0, EXTINCT
100   CONTINUE

      WRITE (2,*) 'C  EMISSION VECTOR'
      ABSORB = (1.0-ALBEDO)*EXTINCT
      DO 150 L = 1, 2
        DO 150 J = 1, NUMMU
          MU = MU_VALUES(J)
          IF (L .EQ. 2) MU = -MU
          WRITE (2,'(1X,F11.8,2X,4(1X,E15.8))') MU, ABSORB,0.0,0.0,0.0
150   CONTINUE

      RETURN
      END

      SUBROUTINE SCATTERING_CNV(CONST, MU_VALUES, NLEGEN, COEF)

      use nml_params, only: nstokes, nummu, aziorder

      implicit none

      INTEGER  NLEGEN, l
      REAL*8   MU_VALUES(NUMMU), CONST, COEF(6,1)
      INTEGER  MAXLEG
      PARAMETER (MAXLEG=64)
      INTEGER  I1, I2, I, J1, J2, K, L1, L2, M
      INTEGER  NUMPTS, DOSUM(6)
      REAL*8   C, MU1, MU2,  DELPHI, COS_SCAT
      REAL*8   PHASE_MATRIX(4,4)
      REAL*8   SCAT_MATRIX(4,4,4*MAXLEG), BASIS_MATRIX(4,4,4*MAXLEG)
      REAL*8   ZERO, TWOPI
      PARAMETER (ZERO=0.0D0, TWOPI=2.0D0*3.1415926535897932384D0)


      NUMPTS = 2* 2**INT(LOG(FLOAT(NLEGEN+4))/LOG(2.0d0)+1.0d0)
      IF (AZIORDER .EQ. 0)  NUMPTS = 2*INT((NLEGEN+1)/2) + 4

!           Find how many Legendre series must be summed
      CALL NUMBER_SUMS (NSTOKES, NLEGEN, COEF, DOSUM)

      DO 50 I = 1, NLEGEN+1
        DO 50 K = 1, 6
          COEF(K,I) = CONST*COEF(K,I)
50    CONTINUE

      WRITE (2,*) 'C  SCATTERING MATRIX'

!       MU1 is the incoming direction, and MU2 is the outgoing direction.
      DO 150 L1 = 1, 2
        DO 150 J1 = 1, NUMMU
          MU1 = MU_VALUES(J1)
          IF (L1 .EQ. 2)  MU1 = -MU1
          DO 150 L2 = 1, 2
            L = 2*(L2-1)+L1
            DO 150 J2 = 1, NUMMU
              MU2 = MU_VALUES(J2)
              IF (L2 .EQ. 2)  MU2 = -MU2
!                   Only need to calculate phase matrix for half of
!                     the delphi's, the rest come from symmetry.
              DO 100 K = 1, NUMPTS/2 + 1
                DELPHI = (TWOPI*(K-1))/NUMPTS
                COS_SCAT = MU1*MU2 + DSQRT((1.-MU1**2)*(1.-MU2**2))*&
                                 DCOS(DELPHI)
                CALL SUM_LEGENDRE(NLEGEN, COEF, COS_SCAT,&
                                DOSUM, PHASE_MATRIX)
                CALL ROTATE_PHASE_MATRIX (PHASE_MATRIX, MU1, MU2,&
                      DELPHI, COS_SCAT, SCAT_MATRIX(1,1,K), NSTOKES)
                CALL MATRIX_SYMMETRY (NSTOKES, SCAT_MATRIX(1,1,K),&
                                   SCAT_MATRIX(1,1,NUMPTS-K+2) )
100           CONTINUE
            CALL FOURIER_MATRIX (AZIORDER, NUMPTS, NSTOKES,&
                        SCAT_MATRIX, BASIS_MATRIX)
            WRITE (2,111) MU1, MU2, 0, ' '
            DO 110 I2 = 1, NSTOKES
                WRITE (2,113) (BASIS_MATRIX(I2,I1,1), I1=1,NSTOKES)
110         CONTINUE

150   CONTINUE

111         FORMAT (1X,F11.8,2X,F11.8,3X,I3,1X,A)
113         FORMAT (1X,4(1X,E15.8))

      RETURN
      END
