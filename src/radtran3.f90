!      SUBROUTINE RADTRAN (NSTOKES, NUMMU, AZIORDER, MAX_DELTA_TAU,     
!     .                    SRC_CODE, QUAD_TYPE, DELTAM,                 
!     .                    DIRECT_FLUX, DIRECT_MU,                      
!     .                    GROUND_TEMP, GROUND_TYPE,                    
!     .                    GROUND_ALBEDO, GROUND_INDEX,                 
!     .                    SKY_TEMP, WAVELENGTH,                        
!     .                    NUM_LAYERS, HEIGHT, TEMPERATURES,            
!     .                    GAS_EXTINCT, SCAT_FILES,                     
!     .                    NOUTLEVELS, OUTLEVELS,                       
!     .                    MU_VALUES, UP_FLUX, DOWN_FLUX,               
!     .                    UP_RAD, DOWN_RAD)                            
!                                                                       
!        RADTRAN solves the plane-parallel polarized radiative transfer
!    equation for an inhomogenous atmosphere with randomly oriented
!    particles.  The sources of radiation are solar direct beam and thermal
!    emission.  The ground surface has Lambertian or Fresnel reflection and
!    emission.
!        Input is the relevant parameters for the atmospheric layers and
!    the boundary conditions for a number of different cases.  Output is
!    a Fourier azimuth series of the upwelling and downwelling radiances 
!    at selected levels at the discrete angles.
!        The heights and temperatures are specified at the layer interfaces,
!    i.e. N+1 values for N layers.  The gaseous extinction and scattering
!    files are specified for the N layers.  The layers are listed from the
!    top to the bottom of the atmosphere: HEIGHT(1) is top of the highest
!    layer and HEIGHT(N+1) is bottom of the lowest layer, SCAT_FILES(1) is
!    the scattering file of the top layer and SCAT_FILES(N) is the scattering
!    file of the bottom layer.
!
!
!    Parameter         Type             Description
!
!  NSTOKES           INTEGER       Number of Stokes parameters: 1 for I
!                                    (no polarization), 2 for I,Q,
!                                    3 for I,Q,U,  4 for I,Q,U,V.
!  NUMMU             INTEGER       Number of quadrature angles
!                                    (per hemisphere).
!  AZIORDER          INTEGER       Order of Fourier azimuth series:
!                                    0 is azimuthially symmetric case.
!  MAX_DELTA_TAU     REAL          Initial layer thickness for doubling;
!                                    governs accuracy, 10E-5 should be
!                                    adequate.  Don't go beyond half the
!                                    real precision, i.e. 10e-8 for REAL*8.
!  SRC_CODE          INTEGER       Radiation sources included:
!                                    none=0, solar=1, thermal=2, both=3
!  QUAD_TYPE         CHAR*1        Type of quadrature used: 
!                                    (G-gaussian, D-double gaussian, 
!                                     L-Lobatto, E-extra-angle).
!                                    If extra-angles then end of
!                                    MU_VALUES(<=NUMMU) contains the extra
!                                    angles and rest is zero.
!  DELTAM            CHAR*1        Delta-M scaling flag (Y or N)
!                                    If DELTAM='Y' then the scattering
!                                    properties are delta-M scaled when read in.
!  DIRECT_FLUX       REAL          Flux on horizontal plane from direct
!                                    (solar) source;  units W/(m*m)/um or K.
!  DIRECT_MU         REAL          Cosine of solar zenith angle
!
!  GROUND_TEMP       REAL          Ground surface temperature in Kelvin
!  GROUND_TYPE       CHAR*1        Type of ground surface:
!                                    L for Lambertian, F for Fresnel.
!                                    Only Lambertian allowed for solar source.
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
!  SCAT_FILES        CHAR*64 array Names of scattering files for layers
!                                    String format 'RAIN.SCA', for no
!                                    scattering use ' '.  See example for
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
!                                    (NSTOKES,NUMMU,AZIORDER+1,NOUTLEVELS)
!  DOWN_RAD          REAL array    Downward radiances
!                                    (NSTOKES,NUMMU,AZIORDER+1,NOUTLEVELS)

SUBROUTINE RADTRAN(NSTOKES, NUMMU, AZIORDER, MAX_DELTA_TAU,      &
     SRC_CODE, QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU, GROUND_TEMP, &
     GROUND_TYPE, GROUND_ALBEDO, GROUND_INDEX, SKY_TEMP, WAVELENGTH,   &
     NUM_LAYERS, HEIGHT, TEMPERATURES, GAS_EXTINCT,  &
     NOUTLEVELS, OUTLEVELS, MU_VALUES, UP_FLUX, DOWN_FLUX, UP_RAD,     &
     DOWN_RAD,wind10u,wind10v)

  use kinds
  use vars_atmosphere
  use nml_params, only: verbose

  INTEGER NSTOKES, NUMMU, AZIORDER 
  INTEGER NUM_LAYERS, SRC_CODE 
  INTEGER NOUTLEVELS, OUTLEVELS ( * ) 
  REAL(kind=dbl) DIRECT_FLUX, DIRECT_MU 
  REAL(kind=dbl) GROUND_TEMP, GROUND_ALBEDO 
  COMPLEX(kind=dbl) GROUND_INDEX 
  REAL(kind=dbl) SKY_TEMP 
  REAL(kind=dbl) WAVELENGTH, MAX_DELTA_TAU 
  REAL(kind=dbl) HEIGHT ( * ), TEMPERATURES ( * ), GAS_EXTINCT ( * ) 
  REAL(kind=dbl) MU_VALUES ( * ) 
  REAL(kind=dbl) UP_FLUX ( * ), DOWN_FLUX ( * ) 
  REAL(kind=dbl) UP_RAD ( * ), DOWN_RAD ( * ) 
  CHARACTER(1) QUAD_TYPE, DELTAM, GROUND_TYPE 

  INTEGER MAXV, MAXM, MAXLM, MAXLEG, MAXLAY, MAXSBUF, MAXDBUF 
  !      PARAMETER (MAXV=64, MAXM=4096, MAXLM=201*256)                    
!   PARAMETER (MAXV = 64, MAXM = 4096, MAXLM = 201 * 512) 
!                                    maxlm = layer*(stokes*angles)**2
  PARAMETER (MAXV = 64, MAXM = 4096, MAXLM = 201 * (2*32)**2)
  PARAMETER (MAXLEG = 256, MAXLAY = 200) 
  PARAMETER (MAXSBUF = MAXLAY * 2 * MAXM, MAXDBUF = MAXLAY * 2 *    &
       MAXV)                                                             

  REAL(kind=dbl) PI, TWOPI, ZERO 
  PARAMETER (PI = 3.1415926535897932384D0, TWOPI = 2.0D0 * PI) 
  PARAMETER (ZERO = 0.0D0) 

  INTEGER NUMLEGEN, NLEGLIM, MODE, LAYER, NUM_DOUBLES 
  INTEGER SCAT_NUMS (MAXLAY), SCAT_NUM 
  INTEGER I, J, K, N, NA, KRT, KS, L, LI 
  LOGICAL SYMMETRIC 
  REAL(kind=dbl) EXTINCTIONS (MAXLAY), ALBEDOS (MAXLAY) 
  REAL(kind=dbl) EXPFACTOR, LINFACTOR 
  REAL(kind=dbl) PLANCK0, PLANCK1, TMP 
  REAL(kind=dbl) ALBEDO, EXTINCTION, EXTINCT, SCATTER 
  REAL(kind=dbl) ZDIFF, DELTA_Z, F, NUM_SUB_LAYERS, TAU 
  REAL(kind=dbl) QUAD_WEIGHTS (MAXV) 
  REAL(kind=dbl) LEGENDRE_COEF (6, MAXLEG) 
  REAL(kind=dbl) SCATBUF (MAXSBUF), DIRECTBUF (MAXDBUF) 
  REAL(kind=dbl) SCATTER_MATRIX (4 * MAXM) 
  REAL(kind=dbl) DIRECT_LEVEL_FLUX (MAXLAY + 1) 
  REAL(kind=dbl) DIRECT_VECTOR (2 * MAXV), EXP_SOURCE (2 * MAXV) 
  REAL(kind=dbl) THERMAL_VECTOR (2 * MAXV), LIN_SOURCE (2 * MAXV) 
  REAL(kind=dbl) REFLECT1 (2 * MAXM), UPREFLECT (2 * MAXM), DOWNREFLECT (2 &
       * MAXM)                                                           
  REAL(kind=dbl) TRANS1 (2 * MAXM), UPTRANS (2 * MAXM), DOWNTRANS (2 *     &
       MAXM)                                                             
  REAL(kind=dbl) SOURCE1 (2 * MAXV), UPSOURCE (2 * MAXV), DOWNSOURCE (2 *  &
       MAXV)                                                             
  REAL(kind=dbl) REFLECT (2 * MAXLM) 
  REAL(kind=dbl) TRANS (2 * MAXLM) 
  REAL(kind=dbl) SOURCE (2 * MAXV * (MAXLAY + 1) ) 
  REAL(kind=dbl) GND_RADIANCE (MAXV), SKY_RADIANCE (2 * MAXV) 
  real(kind=dbl) wind10,wind10u,wind10v,windratio,windangle
  integer :: iquadrant
  REAL(KIND=dbl), PARAMETER :: quadcof  (4, 2  ) =      &
    & Reshape((/0.0_dbl, 1.0_dbl, 1.0_dbl, 2.0_dbl, 1.0_dbl,  - 1.0_dbl, 1.0_dbl,  - 1.0_dbl/), (/4, 2/))

  ! variables needed for fastem4

  real(kind=dbl) :: rel_azimuth, salinity
  real(kind=dbl), dimension(nummu) :: transmittance

  if (verbose .gt. 1) print*, "Entered radtran ...."

  IF (SRC_CODE.EQ.1.OR.SRC_CODE.EQ.3) THEN 
     IF (GROUND_TYPE.NE.'L') THEN 
        WRITE (*, * ) 'Solar case requires Lambertian surface' 
        STOP 
     ENDIF
  ENDIF

  SYMMETRIC = .FALSE. 
  IF (NSTOKES.LE.2) SYMMETRIC = .TRUE. 

  N = NSTOKES * NUMMU 
  IF (N.GT.MAXV) THEN 
     WRITE (*, '(1X,A,I3)') 'Vector size exceeded.  Maximum size :', M&
          &AXV                                                               
     STOP 
  ELSEIF (N * N.GT.MAXM) THEN 
     WRITE (*, '(1X,A,I3)') 'Matrix size exceeded.  Maximum size :', M&
          &AXM                                                               
     STOP 
  ENDIF
  IF ( (SRC_CODE.EQ.1.OR.SRC_CODE.EQ.3) .AND.MAXDBUF.LT. (AZIORDER +&
       1) * 2 * N * NUM_LAYERS) THEN                                     
     WRITE (*, '(1X,A,I3)') 'Direct source buffer size exceeded.' 
     STOP 
  ENDIF
  IF (NUM_LAYERS.GT.MAXLAY) THEN 
     WRITE (*, '(1X,A,I3)') 'Too many layers.  Maximum number :', MAXL&
          &AY                                                                
     STOP 
  ENDIF
  IF ( (NUM_LAYERS + 1) * N * N.GT.MAXLM) THEN 
     WRITE (*, '(1X,A,A,I3)') 'Matrix layer size exceeded.', '  Maximu&
          &m number :', MAXLM   
     STOP 
  ENDIF

  !   Make the desired quadrature abscissas and weights           
  if (verbose .gt. 1) print*, "Make the desired quadrature abscissas and weights ...."

  IF (QUAD_TYPE (1:1) .EQ.'D') THEN 
     CALL DOUBLE_GAUSS_QUADRATURE (NUMMU, MU_VALUES, QUAD_WEIGHTS) 
     NLEGLIM = 2 * NUMMU - 3 
  ELSEIF (QUAD_TYPE (1:1) .EQ.'L') THEN 
     CALL LOBATTO_QUADRATURE (NUMMU, MU_VALUES, QUAD_WEIGHTS) 
     NLEGLIM = 4 * NUMMU - 5 
  ELSEIF (QUAD_TYPE (1:1) .EQ.'E') THEN 
     J = NUMMU
     DO I = NUMMU, 1, - 1 
        IF (MU_VALUES (I) .NE.0.0) THEN 
           QUAD_WEIGHTS (I) = 0.0 
           J = I - 1 
        ENDIF
     ENDDO
     CALL LOBATTO_QUADRATURE (J, MU_VALUES, QUAD_WEIGHTS) !$##
     !         CALL GAUSS_LEGENDRE_QUADRATURE (J, MU_VALUES, QUAD_WEIGHTS) 
     NLEGLIM = 4 * J - 3 
  ELSE 
     CALL GAUSS_LEGENDRE_QUADRATURE (NUMMU, MU_VALUES, QUAD_WEIGHTS) 
     NLEGLIM = 4 * NUMMU - 3 
  ENDIF
  NLEGLIM = MAX (NLEGLIM, 1) 

if (verbose .gt. 1) print*, ".... done!"
!       Make all of the scattering matrices ahead of time
!         and store them in memory

  SCAT_NUM = 0
!           Loop through the layers
  DO LAYER = 1, num_layers
     if (rt3kexttot(layer) .gt. 0.0) then
       SCAT_NUM = SCAT_NUM + 1
       IF (SCAT_NUM * (AZIORDER + 1) * 2 * (NUMMU * NSTOKES) **2.GT.MAXSBUF) THEN
          WRITE (*, '(1X,A,I3)') 'Scattering matrix buffer size exceeded.'
          STOP
       ENDIF
       ! get scattering coefficients
       CALL get_scat_coefs(layer, DELTAM, NUMMU, NUMLEGEN,LEGENDRE_COEF, EXTINCT, SCATTER)
           IF (NUMLEGEN.GT.MAXLEG) THEN
              WRITE (*, * ) 'Too many Legendre terms.'
              STOP
           ENDIF
           ! Truncate the Legendre series to enforce normalization
           IF (NUMLEGEN.GT.NLEGLIM) THEN
              WRITE ( * , * ) 'Truncating Legendre series for layer:',  &
                   layer, NUMLEGEN, NLEGLIM
              NUMLEGEN = NLEGLIM
           ENDIF
           ! Make the scattering matrix
           CALL SCATTERING(NUMMU, AZIORDER, NSTOKES, MU_VALUES,       &
                QUAD_WEIGHTS, NUMLEGEN, LEGENDRE_COEF, SCAT_NUM, SCATBUF)
     else
     ! Special case for a non-scattering layer
        EXTINCT = 0.0
        SCATTER = 0.0
     ENDIF

     SCAT_NUMS (LAYER) = SCAT_NUM
     EXTINCTIONS (LAYER) = EXTINCT + MAX (GAS_EXTINCT (LAYER), 0.0D0)
     IF (EXTINCTIONS (LAYER) .GT.0.0) THEN
        ALBEDOS (LAYER) = SCATTER / EXTINCTIONS (LAYER)
     ELSE
        ALBEDOS (LAYER) = 0.0
     ENDIF
  ENDDO

  !           Compute the direct beam flux at each level                  
  IF (SRC_CODE.EQ.1.OR.SRC_CODE.EQ.3) THEN 
     TAU = 0.0 
     DIRECT_LEVEL_FLUX (1) = DIRECT_FLUX 
     DO LAYER = 1, NUM_LAYERS 
        TAU = TAU + EXTINCTIONS (LAYER) / DIRECT_MU * ABS (HEIGHT (    &
             LAYER) - HEIGHT (LAYER + 1) )                                  
        DIRECT_LEVEL_FLUX (LAYER + 1) = DIRECT_FLUX * DEXP ( - TAU) 
     ENDDO
  ENDIF

  !       Loop through each azimuth mode                                  
  DO MODE = 0, AZIORDER 
     SCAT_NUM = 0 
     !     ------------------------------------------------------            
     !           Loop through the layers
     DO LAYER = 1, NUM_LAYERS
        !                   Calculate the layer thickness                       
        ZDIFF = ABS (HEIGHT (LAYER) - HEIGHT (LAYER + 1) ) 
        EXTINCTION = EXTINCTIONS (LAYER)
        ALBEDO = ALBEDOS (LAYER) 

        IF (SCAT_NUMS (LAYER) .NE.SCAT_NUM) THEN 
           SCAT_NUM = SCAT_NUMS (LAYER) 
           !                   Get the scattering matrix from the buffer           
           CALL GET_SCATTERING (NSTOKES, NUMMU, MODE, AZIORDER, SCAT_NUM, &
                SCATBUF, SCATTER_MATRIX)                                       
           !                   Check the normalization of the scattering matrix    
           IF (MODE.EQ.0) THEN 
              CALL CHECK_NORM (NSTOKES, NUMMU, QUAD_WEIGHTS, SCATTER_MATRIX)
           ENDIF
           !                   Get the direct (solar) vector from the buffer       
           IF (SRC_CODE.EQ.1.OR.SRC_CODE.EQ.3) THEN 
              CALL GET_DIRECT (NSTOKES, NUMMU, MODE, AZIORDER, SCAT_NUM,  &
                   DIRECTBUF, DIRECT_VECTOR)                                   
           ENDIF
        ENDIF


        !           Compute the thermal emission at top and bottom of layer
        IF (SRC_CODE.EQ.2.OR.SRC_CODE.EQ.3) THEN 
           !                   Calculate the thermal source for end of layer       
           CALL THERMAL_RADIANCE (NSTOKES, NUMMU, MODE, TEMPERATURES (    &
                LAYER + 1), ALBEDO, WAVELENGTH, THERMAL_VECTOR)                
           PLANCK1 = THERMAL_VECTOR (1) 
           !                   Calculate the thermal source for beginning of layer 
           CALL THERMAL_RADIANCE (NSTOKES, NUMMU, MODE, TEMPERATURES (    &
                LAYER), ALBEDO, WAVELENGTH, THERMAL_VECTOR)                    
           PLANCK0 = THERMAL_VECTOR (1) 
        ELSE 
           PLANCK0 = 0.0 
           PLANCK1 = 0.0 
        ENDIF


        KRT = 1 + 2 * N * N * (LAYER - 1) 
        KS = 1 + 2 * N * (LAYER - 1) 
        IF (ALBEDO.EQ.ZERO) THEN 
           !                   If the layer is purely absorbing then quickly       
           !                     make the reflection and transmission matrices     
           !                     and source vector instead of doubling.
           CALL NONSCATTER_LAYER (NSTOKES, NUMMU, MODE, ZDIFF *           &
                EXTINCTION, MU_VALUES, PLANCK0, PLANCK1, REFLECT (KRT),        &
                TRANS (KRT), SOURCE (KS) )
        ELSE 

           !                   Find initial thickness of sublayer and              
           !                     the number of times to double                     
           F = MAX (EXTINCTION * ZDIFF, 1.0D-7) 
           F = LOG (F / MAX_DELTA_TAU) / LOG (2.) 
           NUM_DOUBLES = 0 
           IF (F.GT.0.0) NUM_DOUBLES = INT (F) + 1 
           NUM_SUB_LAYERS = 2.0**NUM_DOUBLES 
           DELTA_Z = ZDIFF / NUM_SUB_LAYERS 

           !                   For a solar source make the pseudo source vector    
           !                     and initialize it                                 
           IF (SRC_CODE.EQ.1.OR.SRC_CODE.EQ.3) THEN 
              TMP = DIRECT_LEVEL_FLUX (LAYER) * ALBEDO / (4.0D0 * PI *    &
                   DIRECT_MU)                                                  
              CALL MSCALARMULT (2 * N, 1, TMP, DIRECT_VECTOR, SOURCE1) 
              CALL INITIAL_SOURCE (NSTOKES, NUMMU, N, DELTA_Z, MU_VALUES, &
                   EXTINCTION, SOURCE1, EXP_SOURCE)                            
              EXPFACTOR = DEXP ( - EXTINCTION * DELTA_Z / DIRECT_MU) 
           ENDIF

           !                   Initialize the thermal source vector                
           IF (SRC_CODE.EQ.2.OR.SRC_CODE.EQ.3) THEN 
              CALL INITIAL_SOURCE (NSTOKES, NUMMU, N, DELTA_Z, MU_VALUES, &
                   EXTINCTION, THERMAL_VECTOR, LIN_SOURCE)                     
              IF (PLANCK0.EQ.0.0) THEN 
                 LINFACTOR = 0.0 
              ELSE 
                 LINFACTOR = (PLANCK1 / PLANCK0 - 1.0) / NUM_SUB_LAYERS 
              ENDIF
           ENDIF

           !                Generate the local reflection and transmission matrices
           CALL INITIALIZE (NSTOKES, NUMMU, N, DELTA_Z, MU_VALUES,        &
                EXTINCTION, ALBEDO, SCATTER_MATRIX, REFLECT1, TRANS1)          

           !                   Double up to the thickness of the layer             
           CALL DOUBLING_INTEGRATION (N, NUM_DOUBLES, SRC_CODE, SYMMETRIC,&
                REFLECT1, TRANS1, EXP_SOURCE, EXPFACTOR, LIN_SOURCE, LINFACTOR,&
                REFLECT (KRT), TRANS (KRT), SOURCE (KS) )                      

        ENDIF
     ENDDO
     !            End of layer loop

	! transmittance and ground level for all zenith angles (needed for fastem)

	transmittance = trans(krt:krt + (nummu-1)*66:66)

     !           Get the surface reflection and transmission matrices        
     !             and the surface radiance                                  
     KRT = 1 + 2 * N * N * (NUM_LAYERS) 
     KS = 1 + 2 * N * (NUM_LAYERS) 

     if (verbose .gt. 1) print*, "Calculating surface emissivity ...."

      IF (GROUND_TYPE.EQ.'F') THEN
  ! For a Fresnel surface
    wind10 = sqrt(wind10u**2+wind10v**2)
	CALL FRESNEL_SURFACE (NSTOKES, NUMMU, MU_VALUES, GROUND_INDEX, &
	wavelength, wind10, REFLECT (KRT), TRANS (KRT), SOURCE (KS) )
  ! The radiance from the ground is thermal                
	CALL FRESNEL_RADIANCE (NSTOKES, NUMMU, MODE, MU_VALUES,        &
	GROUND_INDEX, GROUND_TEMP, WAVELENGTH, wind10, GND_RADIANCE)           
      ELSEIF (GROUND_TYPE.EQ.'S') THEN
        ! For a specular surface                                   
        CALL specular_surface(NSTOKES, NUMMU, GROUND_ALBEDO, &
             REFLECT (KRT), TRANS (KRT), SOURCE (KS) )
        ! The radiance from the ground is thermal                
        CALL specular_radiance(NSTOKES, NUMMU, MODE,       &
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
        ! For a Lambertian surface                                
        CALL LAMBERT_SURFACE (NSTOKES, NUMMU, MODE, MU_VALUES,         &
             QUAD_WEIGHTS, GROUND_ALBEDO, REFLECT (KRT), TRANS (KRT),       &
             SOURCE (KS) )                                                  
        ! The radiance from the ground is thermal and reflected direct
        CALL LAMBERT_RADIANCE (NSTOKES, NUMMU, MODE, SRC_CODE,         &
             GROUND_ALBEDO, GROUND_TEMP, WAVELENGTH, DIRECT_LEVEL_FLUX (    &
             NUM_LAYERS + 1), GND_RADIANCE)                                 
     ENDIF

     if (verbose .gt. 1) print*, ".... done!"

     ! Assume the radiation coming from above is blackbody radiation
     CALL THERMAL_RADIANCE (NSTOKES, NUMMU, MODE, SKY_TEMP, ZERO,      &
          WAVELENGTH, SKY_RADIANCE)                                         

     !               For each desired output level (1 thru NL+2) add layers  
     !               above and below level and compute internal radiance     
     DO I = 1, NOUTLEVELS 
        LAYER = MIN (MAX (OUTLEVELS (I), 1), NUM_LAYERS + 2) 
        CALL MZERO (2 * N, N, UPREFLECT) 
        CALL MZERO (2 * N, N, DOWNREFLECT) 
        CALL MIDENTITY (N, UPTRANS (1) ) 
        CALL MIDENTITY (N, UPTRANS (1 + N * N) ) 
        CALL MIDENTITY (N, DOWNTRANS (1) ) 
        CALL MIDENTITY (N, DOWNTRANS (1 + N * N) ) 
        CALL MZERO (2 * N, 1, UPSOURCE) 
        CALL MZERO (2 * N, 1, DOWNSOURCE) 
        DO L = 1, LAYER - 1 
           KRT = 1 + 2 * N * N * (L - 1) 
           KS = 1 + 2 * N * (L - 1)
           IF (L.EQ.1) THEN 
              CALL MCOPY (2 * N, N, REFLECT (KRT), UPREFLECT) 
              CALL MCOPY (2 * N, N, TRANS (KRT), UPTRANS) 
              CALL MCOPY (2 * N, 1, SOURCE (KS), UPSOURCE) 
           ELSE 
              CALL MCOPY (2 * N, N, UPREFLECT, REFLECT1) 
              CALL MCOPY (2 * N, N, UPTRANS, TRANS1) 
              CALL MCOPY (2 * N, 1, UPSOURCE, SOURCE1) 
              CALL COMBINE_LAYERS (N, REFLECT1, TRANS1, SOURCE1, REFLECT (   &
                   KRT), TRANS (KRT), SOURCE (KS), UPREFLECT, UPTRANS, UPSOURCE)  
           ENDIF
        ENDDO
        DO L = LAYER, NUM_LAYERS + 1 
           KRT = 1 + 2 * N * N * (L - 1) 
           KS = 1 + 2 * N * (L - 1) 
           IF (L.EQ.LAYER) THEN 
              CALL MCOPY (2 * N, N, REFLECT (KRT), DOWNREFLECT) 
              CALL MCOPY (2 * N, N, TRANS (KRT), DOWNTRANS) 
              CALL MCOPY (2 * N, 1, SOURCE (KS), DOWNSOURCE) 
           ELSE 
              CALL MCOPY (2 * N, N, DOWNREFLECT, REFLECT1) 
              CALL MCOPY (2 * N, N, DOWNTRANS, TRANS1) 
              CALL MCOPY (2 * N, 1, DOWNSOURCE, SOURCE1) 
              CALL COMBINE_LAYERS (N, REFLECT1, TRANS1, SOURCE1, REFLECT (   &
                   KRT), TRANS (KRT), SOURCE (KS), DOWNREFLECT, DOWNTRANS,        &
                   DOWNSOURCE)                                                    
           ENDIF
        ENDDO
        NA = N * (AZIORDER + 1) 
        CALL INTERNAL_RADIANCE (N, UPREFLECT, UPTRANS, UPSOURCE,          &
             DOWNREFLECT, DOWNTRANS, DOWNSOURCE, SKY_RADIANCE, GND_RADIANCE,   &
             UP_RAD (1 + MODE * N + (I - 1) * NA), DOWN_RAD (1 + MODE * N +    &
             (I - 1) * NA) )                                                   
     ENDDO

  ENDDO
  !         End of azimuth mode loop                                      

  !           Integrate mu times the radiance to find the fluxes          
  DO L = 1, NOUTLEVELS 
     NA = (AZIORDER + 1) * NSTOKES * NUMMU 
     DO I = 1, NSTOKES 
        K = I + NSTOKES * (L - 1) 
        UP_FLUX (K) = 0.0 
        DOWN_FLUX (K) = 0.0 
        DO J = 1, NUMMU 
           UP_FLUX (K) = UP_FLUX (K) + TWOPI * QUAD_WEIGHTS (J) * MU_VALUES (&
                J) * UP_RAD (I + NSTOKES * (J - 1) + NA * (L - 1) )               
           DOWN_FLUX (K) = DOWN_FLUX (K) + TWOPI * QUAD_WEIGHTS (J) *        &
                MU_VALUES (J) * DOWN_RAD (I + NSTOKES * (J - 1) + NA * (L - 1) )  
        ENDDO
     ENDDO
     !           Add in direct beam fluxes                                   
     IF (SRC_CODE.EQ.1.OR.SRC_CODE.EQ.3) THEN 
        K = NSTOKES * (L - 1) 
        LI = OUTLEVELS (L) 
        DOWN_FLUX (K + 1) = DOWN_FLUX (K + 1) + DIRECT_LEVEL_FLUX (LI) 
     ENDIF
  ENDDO

  RETURN 
END SUBROUTINE RADTRAN
