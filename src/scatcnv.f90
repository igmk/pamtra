subroutine scatcnv(errorstatus,nlegen,&
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

  use settings, only: nstokes,nummu, aziorder, &
    quad_type, verbose
  use rt_utilities, only: lobatto_quadrature
  use report_module

  implicit none

  INTEGER  NLEGEN
  REAL*8   MU_VALUES(32), QUAD_WEIGHTS(32)
  REAL*8   COEF(6,100),  EXTINCT, ALBEDO, CONST
  CHARACTER*8   QUADTYPE
  real*8 scatter_matrix(nstokes,nummu,nstokes,nummu,4)
  real*8 ext_matrix(nstokes,nstokes,nummu,2)
  real*8 emis_vec(nstokes,nummu,2)

  integer*8 :: errorstatus
  integer*8 :: err = 0
  character(len=80) :: msg
  character(len=40) :: nameOfRoutine = 'scatcnv'

  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

  if (ANY(ISNAN(coef))) then
      msg = 'Error, got nan in legendre coef'
      call report(2, msg, nameOfRoutine)
      errorstatus = err
      return
  end if  
  
  if (verbose >= 8) print*, "nlegen", nlegen
  if (verbose >= 8) print*, "coef", coef


  IF (QUAD_TYPE(1:1) .EQ. 'L') THEN
    QUADTYPE = 'LOBATTO'
    CALL LOBATTO_QUADRATURE(NUMMU, MU_VALUES, QUAD_WEIGHTS)
  ELSE
    QUADTYPE = 'GAUSSIAN'
    CALL GAUSS_LEGENDRE_QUADRATURE (NUMMU, MU_VALUES, QUAD_WEIGHTS)
  ENDIF

  CONST = EXTINCT*ALBEDO/ (4.0D0*3.1415926535897932384D0)

  call transform_scatter(const,mu_values, nlegen, coef, scatter_matrix)

  call transform_ext_emis(MU_VALUES, EXTINCT, ALBEDO,ext_matrix,emis_vec)

  if (verbose >= 3) call report(info,'END of ', nameOfRoutine)

  return

end subroutine scatcnv

subroutine transform_ext_emis(mu_values,extinct,albedo,ext_matrix,emis_vec)

  use kinds
  use settings, only: nstokes, nummu
  
  implicit none

  integer :: j, l
  real(kind=dbl), dimension(nummu) :: mu_values
  real(kind=dbl) :: extinct, albedo
  real(kind=dbl) :: absorb
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

subroutine transform_scatter(CONST,MU_VALUES, NLEGEN, COEF,scatter_matrix)

  use kinds
  use settings, only: nummu, aziorder, nstokes
  use rt_utilities, only: number_sums,&
    sum_legendre,&
    rotate_phase_matrix,&
    matrix_symmetry,&
    fourier_matrix

  implicit none

  INTEGER  NLEGEN, l
  REAL(kind=dbl) :: MU_VALUES(NUMMU), CONST, COEF(6,1)
  INTEGER  MAXLEG
  PARAMETER (MAXLEG=64)
  INTEGER  I1, I2, I, J1, J2, K, L1, L2
  INTEGER  NUMPTS, DOSUM(6)
  REAL(kind=dbl) :: MU1, MU2,  DELPHI, COS_SCAT
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
