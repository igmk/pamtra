!       RT4 solves the plane parallel case of polarized monochromatic
!     radiative transfer for azimuthally symmetric scattering media.
!       This model is briefly described in Evans, K. F., and G. L. Stephens,
!     1995: Microwave radiative transfer through clouds composed of
!     realistically shaped ice crystals. Part II: Remote Sensing of Ice Clouds
!     J. Atmos. Sci., v. 52, 2058-2072.
!
!       This model now allow the output of radiances at any level in the
!     input layer file.  The output format has also been changed.
!
!        Frank Evans,  University of Colorado, Boulder,  May, 1995
!
!
!                          Method
!         This model is a modification of the model for randomly oriented
!     particles (RT3).  Because of the additional complexities of the
!     doubling method associated with the polarization of the direct beam,
!     the solar source parts have been removed, so thermal emission is
!     the only source of radiation.  In a plane-parallel geometry with
!     isotropic thermal emission and particles having random azimuthal
!     orientation, the radiation field is azimuthally symmetric.  This
!     symmetry also implies that the U and V Stokes parameters are zero.
!
!     There may be many arbitrary layers which are uniform and infinite in
!     horizontal extent and may be any thickness.  The properties of the
!     layers are read from a file.  The single scattering properties are
!     read in from other files containing the Stokes scattering (Meuller)
!     matrix, extinction matrix, and emission vector for appropriate sets
!     of discrete quadrature angles.  Linear thermal emission within each
!     layer is calculated.  Thermal emission and reflection from a
!     Lambertian or Fresnel ground surface is incorporated.
!         The doubling and adding technique is used to solve the plane-
!     parallel radiative transfer equation.  Each input layer is divided
!     into a number of homogeneous sublayers with each sublayer being
!     thin enough for the finite difference initialization to be accurate.
!     Infinitesimal generator initialization is used to relate the scattering
!     matrix to the reflection and transmission matrices. The sublayers
!     are integrated with the doubling algorithm.  For each desired output
!     level the transmission, reflection, and source of the layers above
!     and below the level are combined with the adding algorithm.  The
!     internal radiances are computed from the properties of the layers
!     above and below and the incident radiance from the boundaries.
!
!                          Operation
!         First a subroutine is called to get user input.  A bunch of
!     parameters are input with (hopefully) self-explanatory prompts.
!     See radtran4.f for explanation of input parameters.  Note that
!     letter inputs (except for filenames) must be in uppercase and
!     only the first letter is used.
!     The parameters for the layers are read in from the layer file.
!     All files, input and output, are Fortran formatted (text) files.
!     Each line of that file contains the height, temperature, gaseous
!     extinction, and scattering file name. The height and temperature
!     are specified for the interfaces between the layers, while the
!     gaseous extinction, and the scattering file are specified for each
!     layer.  The layers should start at the top and go down.  The
!     scattering file name is a Fortran string in single quotes.
!     The format of the scattering file is documented in radtran4.f and
!     can be found from the subroutine GET_SCAT_FILE or from an example.
!     The oriented scattering file must have the same quadrature angle setup
!     as used in the code.
!
!         The operation is considerably simplified over that of RT3.  There
!     is only the lowest azimuthal Fourier mode to consider.  The single
!     scattering properties are read directly in a form that can be used
!     so there is no conversion from the phase function of scattering angle
!     to the quadrature coordinates.  For each layer the single scattering
!     properties are read in and the normalization is check (to make sure
!     that scattering plus absorption equals extinction). For each layer the
!     infinitesimal generator initialization is used to calculate the
!     local reflection and transmission matrices and source vectors for a
!     very thin sublayer.  A doubling subroutine is then called to calculate
!     the matrices and vectors for the whole layer.  If the layer doesn't
!     scatter then a subroutine calculates the reflection and transmission
!     matrices and source vector rather than using initialization and
!     doubling.  The reflection and transmission matrices and source
!     vectors for each layer are stored in memory.
!
!         After the doubling has been done to find the properties of all
!     the layers, there is a loop over the output levels.  For each
!     output level an adding subroutine is called to combine the layers
!     above and below the output level.  Then an internal radiance subroutine
!     is called to compute the radiance at the output level from the
!     reflection and transmission matrices and source vectors for the
!     medium above and for below and the incident radiance.  There is
!     assumed to be thermal radiance from above and thermal radiance
!     from the lower surface.  The reflection from the lower surface
!     is simply treated as another layer in the medium (with unity
!     transmission and no source).  The radiances for each quadrature
!     direction for the m=0 azimuth mode are output, as are the integrated
!     fluxes, for each output level.  Note, that these output levels
!     have to be at input layer boundaries.
!
!         There are four types of numerical quadrature schemes available.
!     Gaussian, double Gaussian, and Lobatto are standard.  The 'extra-angle'
!     is the same as gaussian quadrature but with extra angles added in. The
!     weights for the extra angles are zero, so the radiances calculated at
!     the gaussian quadrature angles are uneffected.  The radiances at the
!     extra angles are effectively interpolated.
!
!         The output is a text file that contains the parameter values at
!     the beginning followed by the radiance and flux values.  The
!     polarized radiance may be output as I and Q Stokes parameters or
!     as V and H polarizations.  The radiance may be also converted to
!     brightness temperature, though it is always first computed using
!     the Planck function, in Watts/(meter^2 ster micron).  The brightness
!     temperature may be effective blackbody (UNITS=T) or Rayleigh-Jeans (R).
!     The output Stokes parameter are listed together for each angle
!     and height.  A mu of +2 or -2 indicates the hemispheric flux value.
!     Positive mu values are downwelling, and negative are upwelling angles.
!
!                        Program Structure
!         The radiative transfer program is contained in six files.
!     rt4.f has the main program, the user input routine, the layer reading
!     routine, and the output routine.  radtran4.f has the RADTRANO subroutine
!     which performs all of the radiative transfer calculation. It may be
!     called separately, and its parameter passing is documented in
!     radtran4.f.  radscat4.f reads in the oriented scattering file and
!     checks the normalization.  radintg4.f contains routines that perform
!     the initialization, doubling, and adding.  radutil4.f has the Lambertian
!     and Fresnel surface, Planck function, and quadrature routines.
!     radmat.f contains general purpose matrix routines.
!     The location of the individual subroutines is given below.
!
!         The Fortran used for the program is basically Fortran 77 with
!     a few major differences:  long descriptive variable and subroutine
!     names and the use of ENDDO.  In addition the floating point variables
!     are declared real(kind=dbl).  The program will compile and work under
!     VMS and most Unix Fortrans.
!
!                        Data Storage
!         The basis for the radiance vectors has NSTOKES*NUMMU elements.
!     The first dimension is the polarization vector, made up of the Stokes
!     parameters either just I or I and Q.  The second dimension is the
!     quadrature angles for the mu (cosine zenith angle) parameter with
!     NUMMU elements.
!         The scattering matrix variable is actually four matrices:
!     P++, P+-, P-+, and P-- in that order.  Note: positive angles are
!     in the direction of increasing optical depth (down).
!         All real numbers are real(kind=dbl).  The vectors and matrices are
!     one-dimensional arrays in the main program, and adjustable array
!     sizing is used in the subroutines.  The main storage requirement
!     is the scattering, reflection, and transmission matrices for all
!     the layers (6*MAXLAY*MAXM real numbers).  This is wasteful of
!     memory, but is done for speed at computing radiance at many levels.
!
!                        Limits and Caveats
!         The various array size limits are in Fortran parameter statements.
!     MAXV is the maximum size of the basis vector (NSTOKES*NUMMU), while
!     MAXM is the square of MAXV.  MAXLAY is the maximum number of layers.
!     MAXLM is the size of the largest arrays which hold the reflection
!     and transmission matrices for each layer.  The following table gives
!     the location of array size parameters:
!         rt4              MAXV, MAXLAY
!         radtran4.f       MAXV, MAXM, MAXLM, MAXLAY
!         radintg4.f       MAXV, MAXM
!
!     The fractional accuracy of the output radiances is about the size
!     of the MAX_DELTA_TAU parameter.  It is highly recommended that
!     double precision (15 decimals) rather than single precision
!     (7 decimals) be used.
!
!     It is assumed that the scattering and extinction matrices that are
!     read in are symmetric between up and down (P+-=P-+).  If this is
!     not true then the SYMMETRIC logical in RADTRANO should be set to FALSE.
!
!     Routine locations:
!          File         Routines
!        rt4.f          READ_LAYERS, USER_INPUT, OUTPUT_FILE
!        radtran4.f     RADTRANO
!        radutil4.f     LAMBERT_SURFACE, LAMBERT_RADIANCE,
!                       FRESNEL_SURFACE, FRESNEL_RADIANCE,
!                       THERMAL_RADIANCE, PLANCK_FUNCTION,
!                       GAUSS_LEGENDRE_QUADRATURE,
!                       DOUBLE_GAUSS_QUADRATURE, LOBATTO_QUADRATURE
!        radscat4.f     GET_SCAT_FILE, CHECK_NORM
!        radintg4.f     INITIALIZE, INITIAL_SOURCE,
!                       NONSCATTER_LAYER, INTERNAL_RADIANCE,
!                       DOUBLING_INTEGRATION, COMBINE_LAYERS
!        radmat.f       MCOPY, MADD, MSUB, MSCALARMULT, MZERO, MDIAG,
!                       MIDENTITY, MTRANSPOSE, MMULT, MINVERT
!
!
subroutine RT4(errorstatus, nstokes,nummu,mu_values,out_file,quad_type,ground_temp,&
ground_type,ground_albedo,ground_index,sky_temp,&
wavelength,units,outpol,noutlevels,outlevels,nx,ny,fi)

    use kinds
    use vars_atmosphere
    use settings, only: write_nc, in_python, numazimuths
    use report_module

    implicit none

    integer :: nx,ny,fi

    INTEGER   MAXV, MAXLAY
    PARAMETER (MAXV=64)
    PARAMETER (MAXLAY=200)

    INTEGER   NSTOKES, NUMMU
    INTEGER   NUM_LAYERS
    INTEGER   NOUTLEVELS, OUTLEVELS(MAXLAY)
    real(kind=dbl)    GROUND_TEMP, GROUND_ALBEDO
    COMPLEX*16 GROUND_INDEX
    real(kind=dbl)    SKY_TEMP, WAVELENGTH, MAX_DELTA_TAU
    real(kind=dbl)    MU_VALUES(MAXV)
    real(kind=dbl)    HEIGHT(MAXLAY), TEMPERATURES(MAXLAY)
    real(kind=dbl)    GAS_EXTINCT(MAXLAY)
    real(kind=dbl)    UP_RAD(MAXV*(MAXLAY+1)), DOWN_RAD(MAXV*(MAXLAY+1))
    real(kind=dbl)    UP_FLUX(4*(MAXLAY+1)), DOWN_FLUX(4*(MAXLAY+1))
    CHARACTER QUAD_TYPE*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1
    CHARACTER*64 LAYER_FILE, OUT_FILE

    ! Error handling

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'rt4'

    height = 0.
    temperatures = 0.
    gas_extinct = 0.
    MAX_DELTA_TAU = 1.0d-6

    LAYER_FILE=""

    if (verbose >= 1) call report(info, 'Start of ', nameOfRoutine)
    !scat_files = ''
    !scat_files(nlyr) = '1.txt'

    num_layers = nlyr
    height(1:nlyr+1) = hgt_lev(nlyr:0:-1)             ! [m]
    temperatures(1:nlyr+1) = temp_lev(nlyr:0:-1)      ! [K]
    gas_extinct(1:nlyr) = kextatmo(nlyr:1:-1)         ! [Np/m]

    rt4hydros_present(1:nlyr) = hydros_present(nlyr:1:-1)

    rt4scatter_matrix(1:nlyr,:,:,:,:,:) = scattermatrix(nlyr:1:-1,:,:,:,:,:)
    rt4ext_matrix(1:nlyr,:,:,:,:) = extmatrix(nlyr:1:-1,:,:,:,:)
    rt4emis_vec(1:nlyr,:,:,:) = emisvec(nlyr:1:-1,:,:,:)

    !  if (verbose .gt. 0) print*, ".... read_layers done!"

    CALL RADTRAN4(NSTOKES, NUMMU, MAX_DELTA_TAU,&
    QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,&
    GROUND_ALBEDO, GROUND_INDEX,&
    SKY_TEMP, WAVELENGTH,&
    NUM_LAYERS, HEIGHT, TEMPERATURES,&
    GAS_EXTINCT,&
    NOUTLEVELS, OUTLEVELS,&
    MU_VALUES, UP_FLUX, DOWN_FLUX,&
    UP_RAD, DOWN_RAD)

    !  if (verbose .gt. 0) print*, ".... radtran done!"

    if (write_nc .or. in_python) then
        call collect_output(NSTOKES, NUMMU, 0, &
        WAVELENGTH,   &
        UNITS, OUTPOL,NOUTLEVELS, OUTLEVELS,         &
        NUMAZIMUTHS,UP_RAD, DOWN_RAD,     &
        nx,ny,fi)
    else
        CALL OUTPUT_FILE4(NSTOKES, NUMMU,&
        LAYER_FILE, OUT_FILE,&
        QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,&
        GROUND_ALBEDO, GROUND_INDEX,&
        SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,&
        NUM_LAYERS, HEIGHT,&
        NOUTLEVELS, OUTLEVELS, &
        MU_VALUES, UP_FLUX, DOWN_FLUX,&
        UP_RAD, DOWN_RAD)
    end if

    errorstatus = err

    if (verbose >= 1) call report(info, 'End of ', nameOfRoutine)
    return

END subroutine rt4


SUBROUTINE READ_LAYERS (LAYER_FILE, MAXLAY, NUM_LAYERS,&
HEIGHT, TEMPERATURES,&
GAS_EXTINCT, SCAT_FILES)
    use kinds

    implicit none

    INTEGER  MAXLAY, NUM_LAYERS
    real(kind=dbl)   HEIGHT(*), TEMPERATURES(*)
    real(kind=dbl)   GAS_EXTINCT(*)
    CHARACTER*(*)  LAYER_FILE, SCAT_FILES(*)
    INTEGER   I

    !           Read in height, temperature, gaseous extinction, and
    !                 scattering file for the layers
    OPEN (UNIT=1, FILE=LAYER_FILE, STATUS='OLD')
    I = 1
100 CONTINUE
    READ (1,*,ERR=990,END=110) HEIGHT(I), TEMPERATURES(I),&
    GAS_EXTINCT(I), SCAT_FILES(I)
    I = I + 1
    IF (I .EQ. MAXLAY) THEN
        WRITE (*,*) 'Too many layers'
        STOP
    ENDIF
    GOTO 100
110 CONTINUE
    CLOSE(1)
    NUM_LAYERS = I - 2
    RETURN

990 CONTINUE
    WRITE (*,*) 'Error reading layers data file'

    RETURN
END

SUBROUTINE OUTPUT_FILE4(NSTOKES, NUMMU, &
LAYER_FILE, OUT_FILE,&
QUAD_TYPE, GROUND_TEMP, GROUND_TYPE,&
GROUND_ALBEDO, GROUND_INDEX,&
SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,&
NUM_LAYERS, HEIGHT,&
NOUTLEVELS, OUTLEVELS, &
MU_VALUES, UP_FLUX, DOWN_FLUX,&
UP_RAD, DOWN_RAD)

    use kinds

    implicit none

    INTEGER  NSTOKES, NUMMU, NUM_LAYERS
    INTEGER  NOUTLEVELS, OUTLEVELS(*)
    real(kind=dbl)   GROUND_TEMP, GROUND_ALBEDO
    real(kind=dbl)   SKY_TEMP, WAVELENGTH
    real(kind=dbl)   MU_VALUES(NUMMU)
    real(kind=dbl)   HEIGHT(NUM_LAYERS+1)
    real(kind=dbl)   UP_FLUX(NSTOKES,NOUTLEVELS)
    real(kind=dbl)   DOWN_FLUX(NSTOKES,NOUTLEVELS)
    real(kind=dbl)   UP_RAD(NSTOKES,NUMMU,NOUTLEVELS)
    real(kind=dbl)   DOWN_RAD(NSTOKES,NUMMU,NOUTLEVELS)
    COMPLEX*16  GROUND_INDEX
    CHARACTER*(*) LAYER_FILE, OUT_FILE
    CHARACTER  QUAD_TYPE*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1
    CHARACTER*32 QUAD_NAME, UNITS_NAME, GROUND_NAME
    CHARACTER*64 FORM1
    INTEGER  I, J, L, LI


    CALL CONVERT_OUTPUT4(UNITS, OUTPOL, NSTOKES, NUMMU*NOUTLEVELS, &
    WAVELENGTH, 0, UP_RAD)
    CALL CONVERT_OUTPUT4(UNITS, OUTPOL, NSTOKES, NUMMU*NOUTLEVELS, &
    WAVELENGTH, 0, DOWN_RAD)
    CALL CONVERT_OUTPUT4(UNITS, OUTPOL, NSTOKES, NOUTLEVELS, &
    WAVELENGTH, 1, UP_FLUX)
    CALL CONVERT_OUTPUT4(UNITS, OUTPOL, NSTOKES, NOUTLEVELS, &
    WAVELENGTH, 1, DOWN_FLUX)

    QUAD_NAME = 'GAUSSIAN'
    IF (QUAD_TYPE .EQ. 'D')  QUAD_NAME = 'DOUBLEGAUSS'
    IF (QUAD_TYPE .EQ. 'L')  QUAD_NAME = 'LOBATTO'
    IF (QUAD_TYPE .EQ. 'E')  QUAD_NAME = 'EXTRA-ANGLES'
    UNITS_NAME = 'WATTS/(M^2 MICRON STER)'
    IF (UNITS .EQ. 'T') UNITS_NAME = 'KELVINS - EBB'
    IF (UNITS .EQ. 'R') UNITS_NAME = 'KELVINS - RJ'
    GROUND_NAME = 'LAMBERTIAN'
    IF (GROUND_TYPE .EQ. 'F')  GROUND_NAME = 'FRESNEL'

    OPEN (UNIT=3, FILE=OUT_FILE, STATUS='UNKNOWN')

    !           Output the parameters
    WRITE (3,'(A,I3,A,I3,A,I3,A,I1)')&
    'C  NUMMU=', NUMMU,  '  NUMAZI=',1,&
    '  AZIORDER=',0, '  NSTOKES=',NSTOKES
    WRITE (3,'(A,A32)')&
    'C  LAYER_FILE=',    LAYER_FILE
    WRITE (3,'(A,I1,A,A16)')&
    'C  SRC_CODE=',      2,&
    '   QUAD_TYPE=',     QUAD_NAME
    WRITE (3,'(A,F8.2,A,A16)')&
    'C  GROUND_TEMP=',   GROUND_TEMP,&
    '   GROUND_TYPE=',   GROUND_NAME
    IF (GROUND_TYPE(1:1) .EQ. 'F') THEN
        WRITE (3,'(A,2F9.4,A,F8.2)')&
        'C  GROUND_INDEX=',  GROUND_INDEX,&
        '   SKY_TEMP=',      SKY_TEMP
    ELSE
        WRITE (3,'(A,F8.5,A,F8.2)')&
        'C  GROUND_ALBEDO=', GROUND_ALBEDO,&
        '   SKY_TEMP=',      SKY_TEMP
    ENDIF
    WRITE (3,'(A,E12.6)') 'C  WAVELENGTH=',    WAVELENGTH
    WRITE (3,'(A,A25,A,A2)') 'C  UNITS='     ,    UNITS_NAME,&
    '   OUTPUT_POLARIZATION=', OUTPOL


    IF (UNITS(1:1) .EQ. 'T') THEN
        FORM1 = '(F8.1,1X,F8.5,2(1X,F7.2),:)'
    ELSE
        FORM1 = '(F8.1,1X,F8.5,2(1X,E13.6),:)'
    ENDIF
 
    IF (OUTPOL .EQ. 'VH') THEN
        WRITE (3,'(A,A)') 'C    Z       MU    FLUX/RADIANCE (V,H)'
    ELSE
        WRITE (3,'(A,A)') 'C    Z       MU    FLUX/RADIANCE (I,Q)'
    ENDIF
 
    DO L = 1, NOUTLEVELS
        LI = OUTLEVELS(L)
        !               Output fluxes at this level
        WRITE (3,FORM1) HEIGHT(LI), -2.0,&
        (SNGL(UP_FLUX(I,L)),I=1,NSTOKES)
        WRITE (3,FORM1) HEIGHT(LI), +2.0,&
        (SNGL(DOWN_FLUX(I,L)),I=1,NSTOKES)
 
        !           For each zenith at this level output the Stokes parameters.
        !             Output upwelling radiance: -1 < mu < 0
        DO J = NUMMU, 1, -1
            WRITE (3,FORM1) HEIGHT(LI), -MU_VALUES(J),&
            (SNGL(UP_RAD(I,J,L)),I=1,NSTOKES)
        ENDDO
        !             Output downwelling radiance: 0 < mu < 1
        DO J = 1, NUMMU
            WRITE (3,FORM1) HEIGHT(LI), MU_VALUES(J),&
            (SNGL(DOWN_RAD(I,J,L)),I=1,NSTOKES)
        ENDDO
    ENDDO

    CLOSE (3)

    RETURN
END

SUBROUTINE CONVERT_OUTPUT4(UNITS, OUTPOL, NSTOKES, NOUT,&
WAVELEN, FLUXCODE, OUTPUT)
    !       Converts the output radiance or flux arrays to VH polarization
    !     and effective blackbody temperature if desired.  OUTPOL='VH'
    !     converts the polarization basis of the first two Stokes parameters
    !     to vertical/horizontal polarization.  If UNITS='T' the radiance is
    !     converted to effective blackbody brightness temperature, and if
    !     UNITS='R' the radiance is converted to Rayleigh-Jeans brightness
    !     temperature.  If the output is flux then FLUXCODE=1, and the flux
    !     is divided by pi before converting to brightness temperature.
    use kinds

    implicit none

    INTEGER NSTOKES, NOUT, FLUXCODE
    real(kind=dbl)  WAVELEN, OUTPUT(NSTOKES,NOUT)
    CHARACTER UNITS*1, OUTPOL*2
    INTEGER I, J
    real(kind=dbl)  IV, IH, RAD, TEMP

    DO J = 1, NOUT
        !           Convert to Vertical and Horizontal polarization if desired
        IF (OUTPOL .EQ. 'VH') THEN
            IV = 0.5*(OUTPUT(1,J) + OUTPUT(2,J))
            IH = 0.5*(OUTPUT(1,J) - OUTPUT(2,J))
            OUTPUT(1,J) = IV
            OUTPUT(2,J) = IH
        ENDIF
        !           Convert to brightness temperature
        IF (UNITS .EQ. 'T' .OR. UNITS .EQ. 'R') THEN
            DO I = 1, NSTOKES
                RAD = OUTPUT(I,J)
                IF (OUTPOL .EQ. 'VH' .AND. I .LE. 2)  RAD = 2.0*RAD
                IF (FLUXCODE .EQ. 1)  RAD = RAD/ACOS(-1.0)
                IF (UNITS .EQ. 'R') THEN
                    TEMP = RAD * WAVELEN**4 * 1.4388D4/1.1911D8
                ELSE
                    IF (RAD .GT. 0.0) THEN
                        TEMP = 1.4388D4 /&
                        (WAVELEN*DLOG(1.0+ 1.1911D8/(RAD*WAVELEN**5)))
                    ELSE IF (RAD .EQ. 0.0) THEN
                        TEMP = 0.0D0
                    ELSE
                        TEMP = -1.4388D4 /&
                        (WAVELEN*DLOG(1.0+ 1.1911D8/(-RAD*WAVELEN**5)))
                    ENDIF
                ENDIF
                OUTPUT(I,J) = TEMP
            ENDDO
        ENDIF
    ENDDO
    RETURN
END



