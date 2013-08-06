!        RADTRAN4 solves the plane-parallel polarized radiative transfer
!    equation for an inhomogenous atmosphere with particles oriented
!    in zenith angle but randomly oriented in azimuth (azimuthal symmetry).
!    Thermal emission is the only source of radiation.  The radiance is
!    computed in units of Watts/meter^2 micron ster.  The ground surface
!    has Lambertian or Fresnel reflection and emission.
!        Input is the relevant parameters for the atmospheric layers and
!    the boundary conditions.  Output is the upward and downward fluxes
!    and radiances at the discrete quadrature angles for the specified levels.
!        The heights and temperatures are specified at the layer interfaces,
!    i.e. N+1 values for N layers.  The gaseous extinction and scattering
!    files are specified for the N layers.  The layers are listed from the
!    top to the bottom of the atmosphere: HEIGHT(1) is top of the highest
!    layer and HEIGHT(N+1) is bottom of the lowest layer, SCAT_FILES(1) is
!    the scattering file of the top layer and SCAT_FILES(N) is the scattering
!    file of the bottom layer.  The format of the oriented scattering files
!    is given below.
!
!
!    Parameter         Type             Description
!
!  NSTOKES           INTEGER       Number of Stokes parameters: 1 for I
!                                    (no polarization), 2 for I,Q.
!  NUMMU             INTEGER       Number of quadrature angles
!                                    (per hemisphere).
!  MAX_DELTA_TAU     REAL          Initial layer thickness for doubling;
!                                    governs accuracy, 10E-5 should be
!                                    adequate.  Don't go beyond half the
!                                    real precision, i.e. 10e-8 for REAL*8.
!  QUAD_TYPE         CHAR*1        Type of quadrature used:
!                                    (G-gaussian, D-double gaussian,
!                                    L-Lobatto, E-extra-angle).
!                                    If extra-angles then end of
!                                    MU_VALUES(<=NUMMU) contains the extra
!                                    angles and rest is zero.
!  GROUND_TEMP       REAL          Ground surface temperature in Kelvin
!  GROUND_TYPE       CHAR*1        Type of ground surface:
!                                    L for Lambertian, F for Fresnel.
!  GROUND_ALBEDO     REAL          Albedo of Lambertian surface
!  GROUND_INDEX      COMPLEX       Index of refraction of Fresnel surface
!  SKY_TEMP          REAL          Temperature of blackbody radiation
!                                    incident on atmosphere from above
!  WAVELENGTH        REAL          Wavelength of radiation in microns.
!
!  NUM_LAYERS        INTEGER       Number of atmosphere layers input
!  HEIGHT            REAL array    Height of layer interfaces from top down
!                                    Units are inverse of units of extinction
!                                    and scattering, e.g. km.
!  TEMPERATURES      REAL array    Temperature (Kelvins) of layer interfaces
!  GAS_EXTINCT       REAL array    Gaseous (nonscattering) extinction of layers
!                                    For processes not in scattering file
!  SCAT_FILES        CHAR*64 array Names of oriented scattering files for layers
!                                    String format 'PLATE.DDA', for no
!                                    scattering use ' '.  See below for
!                                    format of scattering file.
!
!  NOUTLEVELS        INTEGER       Number of output levels
!  OUTLEVELS         INTEGER       The levels numbers to output at,
!                                    from 1 at top to NUM_LAYERS+1 at bottom.
!
!  MU_VALUES         REAL array    Output quadrature angle values
!                                    (also input for QUAD_TYPE='E')
!  UP_FLUX           REAL array    Upward flux for each Stokes parameter
!                                    at each output level
!                                    UP_FLUX(NSTOKES,NOUTLEVELS)
!  DOWN_FLUX         REAL array    Downward flux (NSTOKES,NOUTLEVELS)
!  UP_RAD            REAL array    Upward radiances
!                                    (NSTOKES,NUMMU,NOUTLEVELS)
!  DOWN_RAD          REAL array    Downward radiances
!                                    (NSTOKES,NUMMU,NOUTLEVELS)
!
!
!             Format of Scattering Files
!
!  Nmu  Maz  'quadtype'   } Nmu is number of quadrature angles in each
!                           hemisphere.  Maz is azimuth series order.
!                           'quadtype' is the quadrature type (e.g. 'GAUSSIAN').
!  C  SCATTERING MATRIX
!                        _   (2*Nmu)^2 (2*Maz+1) sets of scattering matrices
!  mu'  mu   m  c/s       |
!  P11  P12  P13  P14     |  m is azimuth mode: 0, 1, 2, ...
!  P21  P22  P23  P24     |    c/s is C for cosine term, S for sine term
!  P31  P32  P33  P34     |    azimuth ordering is: 0, 1 C, 1 S, 2 C, 2 S, ...
!  P41  P42  P43  P44    _|  mu is outgoing cosine of zenith angle
!                            mu' is incoming cosine of zenith angle
!                            Order of indexing is:  m (fastest), mu, mu'
!  C  EXTINCTION MATRIX
!                        _   2Nmu sets of extinction matrices
!  mu'                    |
!  K11  K12  K13  K14     |  mu' is incoming cosine of zenith angle
!  K21  K22  K23  K24     |
!  K31  K32  K33  K34     |
!  K41  K42  K43  K44    _|
!
!  C  EMISSION VECTOR
!
!  mu'  S1  S2  S3  S4    }  2Nmu sets of emission vectors
!
!  Further description:
!      All angles (mu's) are stored by hemisphere: first from 0 to 1, then
!  from 0 to -1.  There are no blank lines in the file.  There may be any
!  number of comment lines (beginning with "C") before the first data line,
!  but only the three comment lines shown (scattering, extinction, emission)
!  beyond that.  There are spaces between all elements in the format, so that
!  Fortran free format reading will work.
!      The (I,Q,U,V) Stokes basis is used for the polarization.  The frame
!  of reference for the polarization is meridional plane, i.e. the
!  vertical polarization vector is in the plane defined by the ray
!  direction and the vertical (z) axis.  Because of the azimuthal symmetry
!  of the medium, the scattering matrix depends only on the incident and
!  outgoing zenith angles and the difference between the incident and
!  outgoing azimuth angles (delphi).  The delphi dependence is represented
!  as a Fourier series, because the radiative transfer of azimuthal modes
!  is decoupled.  For this code only the m=0 (lowest, symmetric) mode is used.
!  The units of the scattering matrix, extinction matrix, and emission vector
!  are units of inverse length;  K11 = S1 + Integ{P11}, where Integ{} is
!  integration over azimuth and zenith angles.




      SUBROUTINE RADTRAN4(errorstatus, NSTOKES, NUMMU, MAX_DELTA_TAU,&
                    QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,&
                    GROUND_ALBEDO, GROUND_INDEX,&
                    SKY_TEMP, WAVELENGTH,&
                    NUM_LAYERS, HEIGHT, TEMPERATURES,&
                    GAS_EXTINCT,&
                    NOUTLEVELS, OUTLEVELS,&
                    MU_VALUES, UP_FLUX, DOWN_FLUX,&
                    UP_RAD, DOWN_RAD)

      use kinds
      use vars_atmosphere
     use settings, only: verbose
        use report_module
use rt_utilities, only: planck_function,&
gauss_legendre_quadrature,&
double_gauss_quadrature,&
lobatto_quadrature
      implicit none

      INTEGER   NSTOKES, NUMMU, NUM_LAYERS
      INTEGER   NOUTLEVELS, OUTLEVELS(*)
      REAL*8    GROUND_TEMP, GROUND_ALBEDO
      COMPLEX*16  GROUND_INDEX
      REAL*8    SKY_TEMP
      REAL*8    WAVELENGTH, MAX_DELTA_TAU
      REAL*8    HEIGHT(*), TEMPERATURES(*)
      REAL*8    GAS_EXTINCT(*)
      REAL*8    MU_VALUES(*)
      REAL*8    UP_FLUX(*), DOWN_FLUX(*)
      REAL*8    UP_RAD(*), DOWN_RAD(*)
      CHARACTER*1  QUAD_TYPE, GROUND_TYPE

      INTEGER   MAXV, MAXM, MAXLAY, MAXLM
      PARAMETER (MAXV=64, MAXM=4096, MAXLAY=200, maxlm=201 * (maxv)**2)!MAXLM=201*256)

      REAL*8    PI, TWOPI, ZERO
      PARAMETER (PI = 3.1415926535897932384D0, TWOPI=2.0D0*PI)
      PARAMETER (ZERO=0.0D0)

      INTEGER   LAYER, NUM_DOUBLES
      INTEGER   I, J, K, L, N, KRT, KS
      LOGICAL   SYMMETRIC
      REAL*8    LINFACTOR
      REAL*8    PLANCK0, PLANCK1
      REAL*8    ZDIFF, DELTA_Z, F, NUM_SUB_LAYERS, EXTINCT
      REAL*8    QUAD_WEIGHTS(MAXV)
      REAL*8    SCATTER_MATRIX(4*MAXM)
      REAL*8    LIN_SOURCE(2*MAXV)
      REAL*8    EXTINCT_MATRIX(16*2*MAXV), EMIS_VECTOR(4*2*MAXV)
      REAL*8    REFLECT1(2*MAXM),UPREFLECT(2*MAXM),DOWNREFLECT(2*MAXM)
      REAL*8    TRANS1(2*MAXM),  UPTRANS(2*MAXM),  DOWNTRANS(2*MAXM)
      REAL*8    SOURCE1(2*MAXV), UPSOURCE(2*MAXV), DOWNSOURCE(2*MAXV)
      REAL*8    REFLECT(2*MAXLM)
      REAL*8    TRANS(2*MAXLM)
      REAL*8    SOURCE(2*MAXV*(MAXLAY+1))
      REAL*8    GND_RADIANCE(MAXV), SKY_RADIANCE(2*MAXV)
      CHARACTER*64 SCAT_FILE

  real(kind=dbl) wind10,windratio,windangle
  integer :: iquadrant
  REAL(KIND=dbl), PARAMETER :: quadcof(4,2) =      &
       & Reshape((/0.0d0, 1.0d0, 1.0d0, 2.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0/), (/4, 2/))

  ! variables needed for fastem4

  real(kind=dbl) :: rel_azimuth, salinity
  real(kind=dbl), dimension(nummu) :: transmittance

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'radtran4'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
!  if (verbose .gt. 1) print*, "Entered radtran ...."
      SYMMETRIC = .TRUE.
      N = NSTOKES*NUMMU
      IF (N .GT. MAXV) THEN
          WRITE (*,'(1X,A,I3)')&
          'Vector size exceeded.  Maximum size :', MAXV
            msg = 'radtran check'
            call report(err,msg, nameOfRoutine)
            errorstatus = fatal
            return
      ENDIF
      IF (N*N .GT. MAXM) THEN
          WRITE (*,'(1X,A,I3)')&
          'Matrix size exceeded.  Maximum size :', MAXM
            msg = 'radtran check'
            call report(err,msg, nameOfRoutine)
            errorstatus = fatal
            return
      ENDIF
      IF (NUM_LAYERS .GT. MAXLAY) THEN
          WRITE (*,'(1X,A,A,I3)') 'Number of layers exceeded.',&
          '  Maximum number :', MAXLAY
            msg = 'radtran check'
            call report(err,msg, nameOfRoutine)
            errorstatus = fatal
            return
      ENDIF
      IF ((NUM_LAYERS+1)*N*N .GT. MAXLM) THEN
          WRITE (*,'(1X,A,A,I3)') 'Matrix layer size exceeded.',&
          '  Maximum number :', MAXLM
            msg = 'radtran check'
            call report(err,msg, nameOfRoutine)
            errorstatus = fatal
            return
      ENDIF

!           Make the desired quadrature abscissas and weights
      IF (QUAD_TYPE(1:1) .EQ. 'D') THEN
        CALL DOUBLE_GAUSS_QUADRATURE&
                            (NUMMU, MU_VALUES, QUAD_WEIGHTS)
      ELSE IF (QUAD_TYPE(1:1) .EQ. 'L') THEN
        CALL LOBATTO_QUADRATURE&
                            (NUMMU, MU_VALUES, QUAD_WEIGHTS)
      ELSE IF (QUAD_TYPE(1:1) .EQ. 'E') THEN
        J = NUMMU
        DO I = NUMMU, 1, -1
          IF (MU_VALUES(I) .NE. 0.0) THEN
            QUAD_WEIGHTS(I) = 0.0
            J = I - 1
          ENDIF
        ENDDO
        CALL GAUSS_LEGENDRE_QUADRATURE&
                            (J, MU_VALUES, QUAD_WEIGHTS)
      ELSE
        CALL GAUSS_LEGENDRE_QUADRATURE&
                            (NUMMU, MU_VALUES, QUAD_WEIGHTS)
      ENDIF

      SCAT_FILE = '&&&'
!     ------------------------------------------------------
!           Loop through the layers
!              Do doubling to make the reflection and transmission matrices
!              and soure vectors for each layer, which are stored.

      DO LAYER = 1, NUM_LAYERS
!                   Calculate the layer thickness
          ZDIFF = ABS(HEIGHT(LAYER) - HEIGHT(LAYER+1))
          GAS_EXTINCT(LAYER) = MAX(GAS_EXTINCT(LAYER),0.0D0)
!        if (rt4kexttot(layer) .gt. 0.0) then
        if (rt4hydros_present(layer)) then
            call get_scat_mat(layer,NSTOKES, NUMMU,SCATTER_MATRIX,EXTINCT_MATRIX, EMIS_VECTOR)

            call CHECK_NORM4(err,NSTOKES, NUMMU, QUAD_WEIGHTS,&
                                      SCATTER_MATRIX,&
                                      EXTINCT_MATRIX, EMIS_VECTOR)
        end if
        if (err /= 0) then
            msg = 'error in CHECK_NORM4'
            call report(err,msg, nameOfRoutine)
            errorstatus = err
            return
        end if
!          IF (SCAT_FILES(LAYER) .NE. SCAT_FILE .AND. SCAT_FILES(LAYER) .NE. ' ')  THEN
!              SCAT_FILE = SCAT_FILES(LAYER)
!       Read the scattering matrix from the file
!              CALL GET_SCAT_FILE(NSTOKES, NUMMU, QUAD_TYPE, SCAT_FILE,&
!                                 SCATTER_MATRIX,EXTINCT_MATRIX, EMIS_VECTOR)
!              CALL CHECK_NORM4(NSTOKES, NUMMU, QUAD_WEIGHTS,&
!                              SCATTER_MATRIX,&
!                              EXTINCT_MATRIX, EMIS_VECTOR)
!          ENDIF

!                   Do the stuff for thermal source in layer
!                   Calculate the thermal source for end of layer
          CALL PLANCK_FUNCTION (TEMPERATURES(LAYER+1), 'R',WAVELENGTH, PLANCK1)
!                   Calculate the thermal source for beginning of layer
          CALL PLANCK_FUNCTION (TEMPERATURES(LAYER), 'R',WAVELENGTH, PLANCK0)

          KRT = 1 + 2*N*N*(LAYER-1)
          KS = 1 + 2*N*(LAYER-1)

!          IF (SCAT_FILES(LAYER) .EQ. ' ') THEN
!          IF (rt4kexttot(layer) .EQ. 0.d0) THEN
          IF (.not. rt4hydros_present(layer)) THEN
!                   If the layer is purely absorbing then quickly
!                     make the reflection and transmission matrices
!                     and source vector instead of doubling.
              CALL NONSCATTER_LAYER4(NSTOKES, NUMMU, &
                               ZDIFF*GAS_EXTINCT(LAYER), MU_VALUES,&
                               PLANCK0, PLANCK1,&
                               REFLECT(KRT), TRANS(KRT), SOURCE(KS))
          ELSE

!                   Find initial thickness of sublayer and
!                     the number of times to double
              EXTINCT = EXTINCT_MATRIX(1)+GAS_EXTINCT(LAYER)
!     print*, gas_extinct(layer), extinct_matrix(1), extinct
              F =DLOG(MAX(EXTINCT*ZDIFF,1.0D-7)/MAX_DELTA_TAU)/LOG(2.)
              NUM_DOUBLES = 0
              IF (F .GT. 0.0)  NUM_DOUBLES = INT(F) + 1
              NUM_SUB_LAYERS = 2.0**NUM_DOUBLES
              DELTA_Z = ZDIFF / NUM_SUB_LAYERS

!                   Initialize the source vector
              CALL INITIAL_SOURCE4(NSTOKES, NUMMU, DELTA_Z, MU_VALUES,&
                         PLANCK0, EMIS_VECTOR, GAS_EXTINCT(LAYER),&
                         LIN_SOURCE)
              IF (PLANCK0 .EQ. 0.0) THEN
                  LINFACTOR = 0.0
              ELSE
                  LINFACTOR = (PLANCK1/PLANCK0-1.0D0) /NUM_SUB_LAYERS
              ENDIF

!                Generate the local reflection and transmission matrices
              CALL INITIALIZE4(NSTOKES, NUMMU, &
                              DELTA_Z, MU_VALUES, QUAD_WEIGHTS,&
                              GAS_EXTINCT(LAYER), EXTINCT_MATRIX,&
                              SCATTER_MATRIX,  REFLECT1, TRANS1)

!                   Double up to the thickness of the layer
              CALL DOUBLING_INTEGRATION4(N, NUM_DOUBLES, SYMMETRIC,&
                          REFLECT1, TRANS1, LIN_SOURCE, LINFACTOR,&
                          REFLECT(KRT), TRANS(KRT), SOURCE(KS))
          ENDIF

      ENDDO

!            End of layer loop


!           Get the surface reflection and transmission matrices
!             and the surface radiance
      KRT = 1 + 2*N*N*(NUM_LAYERS)
      KS = 1 + 2*N*(NUM_LAYERS)

     if (verbose .gt. 1) print*, "Calculating surface emissivity ...."

      IF (GROUND_TYPE .EQ. 'F') THEN
!               For a Fresnel surface
        CALL FRESNEL_SURFACE (NSTOKES, NUMMU, &
                             MU_VALUES, GROUND_INDEX,&
                             REFLECT(KRT), TRANS(KRT), SOURCE(KS))
!                The radiance from the ground is thermal
        CALL FRESNEL_RADIANCE (NSTOKES, NUMMU,0,&
                       MU_VALUES, GROUND_INDEX, GROUND_TEMP,&
                       WAVELENGTH, GND_RADIANCE)
     ELSEIF (GROUND_TYPE.EQ.'S') THEN
        ! For a specular surface
        CALL specular_surface(NSTOKES, NUMMU, GROUND_ALBEDO, &
             REFLECT (KRT), TRANS (KRT), SOURCE (KS) )
        ! The radiance from the ground is thermal
        CALL specular_radiance(NSTOKES, NUMMU, 0,       &
             GROUND_ALBEDO, GROUND_TEMP, WAVELENGTH, GND_RADIANCE)
     ELSEIF(GROUND_TYPE .EQ. 'O') THEN
    ! call fastem4 ocean emissivity model. the correction due to transmittance is not necessary in
    ! our multi-stream model (?!)
        wind10 = sqrt(wind10u**2+wind10v**2)
        IF (wind10u >= 0.0 .AND. wind10v >= 0.0 ) iquadrant = 1
        IF (wind10u >= 0.0 .AND. wind10v <  0.0 ) iquadrant = 2
        IF (wind10u <  0.0 .AND. wind10v >= 0.0 ) iquadrant = 4
        IF (wind10u <  0.0 .AND. wind10v <  0.0 ) iquadrant = 3
        IF (abs(wind10v) >= 0.0001) THEN
           windratio = wind10u / wind10v
        ELSE
           windratio = 0.0
           IF (abs(wind10u) > 0.0001) THEN
              windratio = 999999.0 * wind10u
           ENDIF
        ENDIF
        windangle        = atan(abs(windratio))
        rel_azimuth = (quadcof(iquadrant, 1) * pi + windangle * quadcof(iquadrant, 2))*180./pi

        ! azimuthal component
        !
        ! the azimuthal component is ignored (rel_azimuth > 360°) when doing simulations for COSMO runs
        ! since we do not know what direction does the satellite have in advance
        !
        rel_azimuth = 400.
        transmittance(:) = 1.
        salinity = 33.
        call fastem4(wavelength   , &  ! Input
             mu_values, &  ! Input
             nummu, &
             ground_temp , &  ! Input
             Salinity    , &  ! Input
             wind10  ,&
             transmittance,&  ! Input, may not be used
             Rel_Azimuth, &  ! Input
             ground_index,&
             GND_RADIANCE, &  ! Output
             REFLECT(KRT), &  ! Output
             trans(krt),&
             source(ks))
     ELSE
!               For a Lambertian surface
        CALL LAMBERT_SURFACE (NSTOKES, NUMMU, 0,&
                            MU_VALUES, QUAD_WEIGHTS, GROUND_ALBEDO,&
                            REFLECT(KRT), TRANS(KRT), SOURCE(KS))
!                The radiance from the ground is thermal and reflected direct
        CALL LAMBERT_RADIANCE (NSTOKES, NUMMU,0, &
              GROUND_ALBEDO, GROUND_TEMP, WAVELENGTH, GND_RADIANCE)
      ENDIF

     if (verbose .gt. 1) print*, ".... done!"

!           Assume the radiation coming from above is blackbody radiation
! 0 stands for mode = 0. this is required, since we use the routine from the former
! radutil3.f
      CALL THERMAL_RADIANCE (NSTOKES, NUMMU,0, SKY_TEMP, ZERO,  &
                            WAVELENGTH,  SKY_RADIANCE)

!         For each desired output level (1 thru NL+2) add layers
!           above and below level and compute internal radiance.
!           OUTLEVELS gives the desired output levels.
      DO I = 1, NOUTLEVELS
        LAYER = MIN( MAX( OUTLEVELS(I), 1), NUM_LAYERS+2)
        CALL MZERO (2*N, N, UPREFLECT)
        CALL MZERO (2*N, N, DOWNREFLECT)
        CALL MIDENTITY (N, UPTRANS(1))
        CALL MIDENTITY (N, UPTRANS(1+N*N))
        CALL MIDENTITY (N, DOWNTRANS(1))
        CALL MIDENTITY (N, DOWNTRANS(1+N*N))
        CALL MZERO (2*N, 1, UPSOURCE)
        CALL MZERO (2*N, 1, DOWNSOURCE)
        DO L = 1, LAYER-1
          KRT = 1 + 2*N*N*(L-1)
          KS = 1 + 2*N*(L-1)
          IF (L .EQ. 1) THEN
            CALL MCOPY (2*N,N, REFLECT(KRT), UPREFLECT)
            CALL MCOPY (2*N,N, TRANS(KRT), UPTRANS) 
            CALL MCOPY (2*N,1, SOURCE(KS), UPSOURCE)
          ELSE
            CALL MCOPY (2*N,N, UPREFLECT, REFLECT1)
            CALL MCOPY (2*N,N, UPTRANS, TRANS1)
            CALL MCOPY (2*N,1, UPSOURCE, SOURCE1)
            CALL COMBINE_LAYERS4(N, REFLECT1, TRANS1, SOURCE1,&
                             REFLECT(KRT), TRANS(KRT), SOURCE(KS),&
                             UPREFLECT, UPTRANS, UPSOURCE)
          ENDIF
        ENDDO
        DO L = LAYER, NUM_LAYERS+1
          KRT = 1 + 2*N*N*(L-1)
          KS = 1 + 2*N*(L-1)
          IF (L .EQ. LAYER) THEN
            CALL MCOPY (2*N,N, REFLECT(KRT), DOWNREFLECT)
            CALL MCOPY (2*N,N, TRANS(KRT), DOWNTRANS) 
            CALL MCOPY (2*N,1, SOURCE(KS), DOWNSOURCE)
          ELSE
            CALL MCOPY (2*N,N, DOWNREFLECT, REFLECT1)
            CALL MCOPY (2*N,N, DOWNTRANS, TRANS1)
            CALL MCOPY (2*N,1, DOWNSOURCE, SOURCE1)
            CALL COMBINE_LAYERS4(N, REFLECT1, TRANS1, SOURCE1,&
                             REFLECT(KRT), TRANS(KRT), SOURCE(KS),&
                             DOWNREFLECT, DOWNTRANS, DOWNSOURCE)
          ENDIF
        ENDDO
        CALL INTERNAL_RADIANCE4(N, UPREFLECT, UPTRANS, UPSOURCE,&
                               DOWNREFLECT, DOWNTRANS, DOWNSOURCE, &
                               SKY_RADIANCE, GND_RADIANCE,&
                               UP_RAD(1+(I-1)*N), DOWN_RAD(1+(I-1)*N))
      ENDDO



!           Integrate the mu times the radiance to find the fluxes
      DO L = 1, NOUTLEVELS
        DO I = 1, NSTOKES
          K = I+NSTOKES*(L-1)
          UP_FLUX(K) = 0.0
          DOWN_FLUX(K) = 0.0
          DO J = 1, NUMMU
            UP_FLUX(K) = UP_FLUX(K)&
                    + TWOPI*QUAD_WEIGHTS(J) * MU_VALUES(J)&
                    * UP_RAD(I+NSTOKES*(J-1)+N*(L-1))
            DOWN_FLUX(K) = DOWN_FLUX(K)&
                    + TWOPI*QUAD_WEIGHTS(J) * MU_VALUES(J)&
                    * DOWN_RAD(I+NSTOKES*(J-1)+N*(L-1))
          ENDDO
        ENDDO
      ENDDO

    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)


      RETURN
      END

