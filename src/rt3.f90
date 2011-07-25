!     RT3 solves the  plane parallel case of polarized monochromatic
!     radiative transfer for isotropic media.  This model is described
!     in K. F. Evans and G. L. Stephens, 1991: A New Polarized Atmospheric
!     Radiative Transfer Model, J. Quant. Spectrosc. Radiat. Transfer,
!     v. 46, no. 5, pp. 413-423, 1991.
!
!       This model now allows the output of radiances at any level in the
!     input layer file.  The output format has also been changed.  May 1995
!       Delta-M scaling added. July 1996
!
!        Frank Evans,  University of Colorado, Boulder,  July, 1996
!
!                          Method
!         The radiation field may have full angular dependence (polar and
!     azimuthal angles).  The scattering particles are assumed to be
!     randomly oriented (isotropic media) and to have a plane of symmetry.
!     The layers are uniform and infinite in horizontal extent and may be
!     any thickness.  The properties of the layers are read from a file.
!     Solar and/or thermal sources of radiation are treated.
!         For randomly oriented particles with a plane of symmetry there are
!     six unique elements (out of 16) in the scattering phase matrix.
!     The phase matrix information is read in from a file containing the
!     coefficients for Legendre polynomial expansions for the six elements.
!     Linear thermal emission within each layer is calculated.
!     Thermal emission and reflection from a Lambertian or Fresnel
!     ground surface is incorporated.
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
!     See radtran3.f for explanation of input parameters.  Note that
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
!     The format of the scattering file can be found from the subroutine
!     READ_SCAT_FILE or from an example.
!
!         The calculation is performed sequentially for each azimuth mode.
!     Before the loop over the azimuth modes, the Fourier modes of the
!     phase matrix for each scattering file is made. The phase matrices
!     are stored temporily either in memory.  Then there is a loop for
!     the azimuth modes, and inside that there is a loop for each layer.
!     For each layer the phase matrix and pseudo-source vector is retrieved,
!     and the normalization of the phase matrix is checked.  The thermal
!     source vector is made if needed.  The parameters for the doubling
!     of the sources are made.  The source is linear in optical depth
!     for the thermal case and exponential in optical depth for the
!     solar case.  For each layer the infinitesimal generator initialization
!     is used to calculate the local reflection and transmission matrices
!     and source vectors for a very thin sublayer.  A doubling subroutine
!     is then called to calculate the matrices and vectors for the whole
!     layer.  If the layer doesn't scatter then a subroutine calculates the
!     reflection and transmission matrices and source vector rather than
!     using initialization and doubling.  The reflection and transmission
!     matrices and source vectors for each layer are stored in memory.
!
!         After the doubling has been done to find the properties of all
!     the layers, there is a loop over the output levels.  For each
!     output level an adding subroutine is called to combine the layers
!     above and below the output level.  Then an internal radiance subroutine
!     is called to compute the radiance at the output level from the
!     reflection and transmission matrices and source vectors for the
!     medium above and for below and the incident radiance.  There is
!     assumed to be thermal radiance from above and thermal and/or
!     reflected direct solar radiance from the lower surface.  The
!     reflection from the lower surface is simply treated as another
!     layer in the medium (with unity transmission and no source).  The
!     RADTRAN subroutine computes the Fourier azimuthal modes for the
!     discrete quadrature zenith angles at each output level, which must
!     be at the input layer boundaries.
!
!         There are four types of numerical quadrature schemes available.
!     Gaussian, double Gaussian, and Lobatto are standard.  The 'extra-angle'
!     is the same as gaussian quadrature but with extra angles added in. The
!     weights for the extra angles are zero, so the radiances calculated at
!     the gaussian quadrature angles are uneffected.  The radiances at the
!     extra angles are effectively interpolated.  The user defined quadrature
!     was removed because of instabilities.
!
!         The output is a text file that contains the parameter values at
!     the beginning followed by the radiance and flux values.  The
!     polarized radiance may be output as I, Q, U, V Stokes parameters or
!     as V, H, U, V.  The radiance may be also converted to
!     brightness temperature, though it is always first computed using
!     the Planck function, in Watts/(meter^2 ster micron).  The brightness
!     temperature may be effective blackbody (UNITS=T) or Rayleigh-Jeans (R).
!     The output Stokes parameter are listed together for each angle
!     and height.  A mu of +2 or -2 indicates the hemispheric flux value.
!     Positive mu values are downwelling, and negative are upwelling angles.
!     For the output the Fourier azimuth series is summed and the radiances
!     output at a specified number of evenly spaced azimuthal angles between
!     0 and 180 degrees inclusive (NUMAZIMUTHS=3 means 0, 90, 180 degrees).
!
!                        Program Structure
!         The radiative transfer program is contained in six files.
!     rt3.f has the main program, the user input routine, the layer reading
!     routine, and the output routine. radtran3.f has the RADTRAN subroutine
!     which performs all of the radiative transfer calculation. It may be
!     called separately, and its parameter passing is documented
!     in radtran3.f.  radscat3.f has routines that convert the Legendre series
!     scattering file data into the phase matrices.  radintg3.f contains
!     routines that perform the initialization, doubling, and adding.
!     radutil3.f has the Lambertian and Fresnel surface, Planck function, and
!     quadrature routines.  radmat.f contains general purpose matrix routines.
!     The location of the individual subroutines is given below.
!
!         The Fortran used for the program is basically Fortran 77 with
!     a few major differences:  long descriptive variable and subroutine
!     names and the use of ENDDO.  In addition the floating point variables
!     are declared REAL*8.  The program will compile and work under
!     VMS and most Unix Fortrans.
!
!                        Data Storage
!         The basis for the radiance vectors has NSTOKES*NUMMU elements.
!     The first dimension is the polarization vector, made up of the Stokes
!     parameters.  As described in the paper the Stokes vector is the
!     cosine azimuth modes of I and Q, and the sine modes of U and V
!     (Ic,Qc,Us,Vs).  The second dimension is the quadrature angles for the
!     mu (cosine theta) parameter with NUMMU elements.
!         The scattering matrix variable is actually four matrices:
!     P++, P+-, P-+, and P-- in that order.  Note: positive angles are
!     in the direction of increasing optical depth (down).
!         All real numbers are REAL*8.  The vectors and matrices are
!     one-dimensional arrays in the main program, and adjustable array
!     sizing is used in the subroutines.  The main storage requirement
!     is the scattering, reflection, and transmission matrices for all
!     the layers (6*MAXLAY*MAXM real numbers).  This is wasteful of
!     memory, but is done for speed at computing radiance at many levels.
!
!                        Limits and Caveats
!         The various array size limits are in Fortran parameter statements.
!     MAXV is the maximum size of the basis vector (NSTOKES*NUMMU), while
!     MAXM is the square of MAXV.  MAXA is the maximum number of azimuth
!     modes.  MAXLEG is the maximum number of terms in the Legendre series.
!     MAXLAY is the maximum number of layers.  MAXLM is the size of the
!     largest arrays which hold the reflection and transmission matrices
!     for each layer.  The following table gives the location of array
!     size parameters:
!         rt3.f             MAXV, MAXA, MAXLAY
!         radtran3.f        MAXV, MAXM, MAXLM, MAXLEG, MAXLAY
!         radscat3.f        MAXLEG
!         radintg3.f        MAXV, MAXM
!     The RADTRAN subroutine also has array sizes for the scattering
!     matrices temporary buffer (MAXSBUF) and the direct source temporary
!     buffer (MAXDBUF).
!
!     The fractional accuracy of the output radiances is about the size
!     of the MAX_DELTA_TAU parameter.  It is highly recommended that
!     double precision (15 decimals) rather than single precision
!     (7 decimals) be used.
!
!     It is important that the phase matrix be normalized so energy is
!     conserved, by having enough quadrature angles for the number of
!     Legendre terms in the scattering file.  This limit on the number of
!     Legendre terms is enforced by truncating the series at the appropriate
!     degree (4*NUMMU-3 for Gaussian quadrature) and printing out a warning
!     message giving the scattering file name.  The Legendre series truncation
!     should prevent the normalization check from failing, but the
!     normalization check is not always sufficient for energy to be conserved.
!     If truncation occurs then the true shape of the phase function is
!     distorted, which may not be significant if the truncation is modest.
!     To avoid truncation the number of quadrature angles may be increased,
!     but note that the CPU time goes as the cube of number of angles.
!     Delta-M scaling may be used to avoid truncation for highly peaked
!     phase functions.  Note, that for collimated (solar) problems the
!     delta-M solution will be accurate for flux, but may still have
!     radiances which oscillate unphysically (though less than with
!     truncating the original phase function).
!     
!
!     Specular surface reflection is not treated so only the Lambertian
!     surface may be used with a solar source.
!
!
!     Routine locations:
!          File         Routines
!        rt3.f          READ_LAYERS, USER_INPUT, OUTPUT_FILE
!        radtran3.f     RADTRAN
!        radutil3.f     LAMBERT_SURFACE, LAMBERT_RADIANCE,
!                       FRESNEL_SURFACE, FRESNEL_RADIANCE,
!                       THERMAL_RADIANCE, PLANCK_FUNCTION,
!                       GAUSS_LEGENDRE_QUADRATURE,
!                       DOUBLE_GAUSS_QUADRATURE, LOBATTO_QUADRATURE
!        radscat3.f     READ_SCAT_FILE, SCATTERING, GET_SCATTER,
!                       CHECK_NORM, DIRECT_SCATTERING, GET_DIRECT
!                       (SUM_LEGENDRE, NUMBER_SUMS,
!                        ROTATE_PHASE_MATRIX, SUM_MATRIX)
!        radintg3.f     INITIALIZE, INITIAL_SOURCE,
!                       NONSCATTER_LAYER, INTERNAL_RADIANCE,
!                       DOUBLING_INTEGRATION, COMBINE_LAYERS
!        radmat.f       MCOPY, MADD, MSUB, MSCALARMULT, MZERO, MDIAG,
!                       MIDENTITY, MTRANSPOSE, MMULT, MINVERT
!
!
      SUBROUTINE RT3 (NSTOKES, NUMMU, AZIORDER, MU_VALUES, SRC_CODE,    &
      LAYER_FILE, OUT_FILE, QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,  &
      GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO, GROUND_INDEX, SKY_TEMP,  &
      WAVELENGTH, UNITS, OUTPOL, NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS,    &
      nx,ny,write_nc,verbose)

      use kinds
      use vars_atmosphere
  
      INTEGER MAXV, MAXA, MAXLAY 
      PARAMETER (MAXV = 64, MAXA = 32) 
      PARAMETER (MAXLAY = 200) 
                                                                        
      INTEGER NSTOKES, NUMMU, AZIORDER
      INTEGER NUM_LAYERS, SRC_CODE 
      INTEGER NOUTLEVELS, OUTLEVELS (MAXLAY), NUMAZIMUTHS 
      REAL(kind=dbl) GROUND_TEMP, GROUND_ALBEDO 
      COMPLEX(kind=dbl) GROUND_INDEX 
      REAL(kind=dbl) SKY_TEMP, WAVELENGTH, MAX_DELTA_TAU 
      REAL(kind=dbl) DIRECT_FLUX, DIRECT_MU 
      REAL(kind=dbl) MU_VALUES (MAXV) 
      REAL(kind=dbl) HEIGHT (MAXLAY), TEMPERATURES (MAXLAY) 
      REAL(kind=dbl) GAS_EXTINCT (MAXLAY) 
      REAL(kind=dbl) UP_FLUX (4 * MAXLAY), DOWN_FLUX (4 * MAXLAY) 
      REAL(kind=dbl) UP_RAD (MAXLAY * MAXA * MAXV), DOWN_RAD (MAXLAY * MAXA *  &
      MAXV)                                   
      CHARACTER QUAD_TYPE * 1, DELTAM * 1, UNITS * 1, OUTPOL * 2,       &
      GROUND_TYPE * 1                                                   
      CHARACTER(68) LAYER_FILE 
      CHARACTER(100) OUT_FILE
      CHARACTER(64) SCAT_FILES (MAXLAY) 
                                                                        
      integer :: model_i, model_j
      real lon,lat,lfrac,wind10u,wind10v,iwv,cwp,iwp,rwp,swp,gwp

      integer :: verbose
      logical :: write_nc

      ! initialization

      height = 0.
      temperatures = 0.
      gas_extinction = 0.

	  nlyr = profiles(nx,ny)%nlyr
      model_i = profiles(nx,ny)%isamp
      model_j = profiles(nx,ny)%jsamp
      lon = profiles(nx,ny)%longitude
      lat = profiles(nx,ny)%latitude
      lfrac = profiles(nx,ny)%land_fraction
      wind10u = profiles(nx,ny)%wind_10u
      wind10v = profiles(nx,ny)%wind_10v

      iwv = profiles(nx,ny)%iwv
      cwp = profiles(nx,ny)%cwp
      iwp = profiles(nx,ny)%iwp
      rwp = profiles(nx,ny)%rwp
      swp = profiles(nx,ny)%swp
      gwp = profiles(nx,ny)%gwp
!       write(27,*)'inputs rt3',NSTOKES, NUMMU,AZIORDER, MU_VALUES,     
!     .                 SRC_CODE, LAYER_FILE, OUT_FILE,                 
!     .                 QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,      
!     .                 GROUND_TEMP, GROUND_TYPE,                       
!     .                 GROUND_ALBEDO, GROUND_INDEX,                    
!     .                 SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,            
!     .                 NOUTLEVELS, OUTLEVELS,NUMAZIMUTHS               
                                                                        
      IF (1 + AZIORDER.GT.MAXA) THEN 
         WRITE ( * , * ) 'Maximum number of azimuth modes exceeded.' 
         STOP 
      ENDIF 
                                                                        
      if (verbose .gt. 1) print*, "reading layers"                            
      CALL READ_LAYERS (LAYER_FILE, MAXLAY, NUM_LAYERS, HEIGHT,         &
      TEMPERATURES, GAS_EXTINCT, SCAT_FILES)

      height(1:nlyr+1) = profiles(nx,ny)%hgt_lev(nlyr:0:-1)
      temperatures(1:nlyr+1) = profiles(nx,ny)%temp_lev(nlyr:0:-1)
      gas_extinct(1:nlyr+1) = kextatmo(nlyr:0:-1)

      if (verbose .gt. 1) print*, "read layers,now calling RT "                       
                                                                        
      MAX_DELTA_TAU = 1.0E-6 

      if (verbose .gt. 1) print*, "Entering radtran ...."

      CALL RADTRAN(NSTOKES, NUMMU, AZIORDER, MAX_DELTA_TAU, SRC_CODE,  &
      QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU, GROUND_TEMP,           &
      GROUND_TYPE, GROUND_ALBEDO, GROUND_INDEX, SKY_TEMP, WAVELENGTH,   &
      NUM_LAYERS, HEIGHT, TEMPERATURES, GAS_EXTINCT, SCAT_FILES,        &
      NOUTLEVELS, OUTLEVELS, MU_VALUES, UP_FLUX, DOWN_FLUX, UP_RAD,     &
      DOWN_RAD,dble(wind10u),dble(wind10v),verbose)

      if (verbose .gt. 1) print*, ".... radtran done!"

!            write(27,*)'RT done'                                       
!       N = NUMMU*(AZIORDER+1)*NOUTLEVELS                               
!      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, N,                  
!     .                     WAVELENGTH, 0, UP_RAD)                      
                                                                        

      if (verbose .gt. 1) print*, "Writing output ...."

      if (write_nc) then 
	call collect_output(NSTOKES, NUMMU, AZIORDER, &
	QUAD_TYPE, GROUND_TEMP, GROUND_TYPE, WAVELENGTH,   &
	UNITS, OUTPOL,NOUTLEVELS, OUTLEVELS,         &
	NUMAZIMUTHS, UP_FLUX, DOWN_FLUX, UP_RAD, DOWN_RAD,     &
	lon,lat,lfrac,wind10u,wind10v,iwv,cwp,iwp,rwp,swp,gwp,model_i,model_j,nx,ny)
      else
	CALL OUTPUT_FILE (NSTOKES, NUMMU, AZIORDER, SRC_CODE, LAYER_FILE, &
	OUT_FILE, QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU, GROUND_TEMP, &
	GROUND_TYPE, GROUND_ALBEDO, GROUND_INDEX, SKY_TEMP, WAVELENGTH,   &
	UNITS, OUTPOL, NUM_LAYERS, HEIGHT, NOUTLEVELS, OUTLEVELS,         &
	NUMAZIMUTHS, MU_VALUES, UP_FLUX, DOWN_FLUX, UP_RAD, DOWN_RAD,     &
	lon,lat,lfrac,wind10,iwv,cwp,iwp,rwp,swp,gwp)
      end if

      if (verbose .gt. 1) print*, ".... done!"
                                                                        
      END SUBROUTINE RT3                            
                                                                        
                                                                        
      SUBROUTINE RT3_oneangle (angle_view, NSTOKES, NUMMU, AZIORDER,    &
      MU_VALUES, SRC_CODE, LAYER_FILE, OUT_FILE, QUAD_TYPE, DELTAM,     &
      DIRECT_FLUX, DIRECT_MU, GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO,  &
      GROUND_INDEX, SKY_TEMP, WAVELENGTH, UNITS, OUTPOL, NOUTLEVELS,    &
      OUTLEVELS, NUMAZIMUTHS, I1, I2, Upar, Vpar)                       
                                                                        
      use kinds
  
      INTEGER MAXV, MAXA, MAXLAY 
      PARAMETER (MAXV = 64, MAXA = 32) 
      PARAMETER (MAXLAY = 200) 
                                                                        
      INTEGER NSTOKES, NUMMU, AZIORDER 
      INTEGER NUM_LAYERS, SRC_CODE 
      INTEGER NOUTLEVELS, OUTLEVELS (MAXLAY), NUMAZIMUTHS 
      REAL(kind=dbl) GROUND_TEMP, GROUND_ALBEDO, angle_view 
      COMPLEX(kind=dbl) GROUND_INDEX 
      REAL(kind=dbl) SKY_TEMP, WAVELENGTH, MAX_DELTA_TAU 
      REAL(kind=dbl) DIRECT_FLUX, DIRECT_MU 
      REAL(kind=dbl) MU_VALUES (MAXV) 
      REAL(kind=dbl) HEIGHT (MAXLAY), TEMPERATURES (MAXLAY) 
      REAL(kind=dbl) GAS_EXTINCT (MAXLAY) 
      REAL(kind=dbl) UP_FLUX (4 * MAXLAY), DOWN_FLUX (4 * MAXLAY) 
      REAL(kind=dbl) UP_RAD (MAXLAY * MAXA * MAXV), DOWN_RAD (MAXLAY * MAXA *  &
      MAXV)                        
      CHARACTER QUAD_TYPE * 1, DELTAM * 1, UNITS * 1, OUTPOL * 2,       &
      GROUND_TYPE * 1                                                   
      CHARACTER(68) LAYER_FILE 
      CHARACTER(70) OUT_FILE
      CHARACTER(64) SCAT_FILES (MAXLAY) 
      REAL(kind=dbl) Field_1, field_2, I1, I2, Upar, Vpar 
                                                                        
                                                                        
      IF (1 + AZIORDER.GT.MAXA) THEN 
         WRITE ( * , * ) 'Maximum number of azimuth modes exceeded.' 
         STOP 
      ENDIF 
                                                                        
                                                                        
      CALL READ_LAYERS (LAYER_FILE, MAXLAY, NUM_LAYERS, HEIGHT,         &
      TEMPERATURES, GAS_EXTINCT, SCAT_FILES)                            
                                                                        
                                                                        
      MAX_DELTA_TAU = 1.0E-6 
      CALL RADTRAN (NSTOKES, NUMMU, AZIORDER, MAX_DELTA_TAU, SRC_CODE,  &
      QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU, GROUND_TEMP,           &
      GROUND_TYPE, GROUND_ALBEDO, GROUND_INDEX, SKY_TEMP, WAVELENGTH,   &
      NUM_LAYERS, HEIGHT, TEMPERATURES, GAS_EXTINCT, SCAT_FILES,        &
      NOUTLEVELS, OUTLEVELS, MU_VALUES, UP_FLUX, DOWN_FLUX, UP_RAD,     &
      DOWN_RAD)                                                         
                                                                        
                                                                        
      CALL OUTPUT_FILE (NSTOKES, NUMMU, AZIORDER, SRC_CODE, LAYER_FILE, &
      OUT_FILE, QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU, GROUND_TEMP, &
      GROUND_TYPE, GROUND_ALBEDO, GROUND_INDEX, SKY_TEMP, WAVELENGTH,   &
      UNITS, OUTPOL, NUM_LAYERS, HEIGHT, NOUTLEVELS, OUTLEVELS,         &
      NUMAZIMUTHS, MU_VALUES, UP_FLUX, DOWN_FLUX, UP_RAD, DOWN_RAD,    &
      lon,lat,lfrac,wind10,iwv,cwp,iwp,rwp,swp,gwp)     
                                                                        
      CALL locate (mu_values, NUMMU, dcos (angle_view), i_loc) 
                   !x is bet xx(j) and xx(j+1)                          
      nshift = NSTOKES - 1 
      Field_1 = DOWN_RAD (NSTOKES * nummu * (Noutlevels - 1) + Nstokes *&
      i_loc - nshift)                                                   
      Field_2 = DOWN_RAD (NSTOKES * nummu * (Noutlevels - 1) + Nstokes *&
      (i_loc + 1) - nshift)                                             
!       write(23,*)'field1',Field_1,field_2                             
      nshift = NSTOKES - 2 
      I1 = Field_1 + (Field_2 - Field_1) * (dcos (angle_view) -         &
      mu_values (i_loc) ) / (mu_values (i_loc + 1) - mu_values (i_loc) )
      Field_1 = DOWN_RAD (NSTOKES * nummu * (Noutlevels - 1) + Nstokes *&
      i_loc - nshift)                                                   
      Field_2 = DOWN_RAD (NSTOKES * nummu * (Noutlevels - 1) + Nstokes *&
      (i_loc + 1) - nshift)                                             
      I2 = Field_1 + (Field_2 - Field_1) * (dcos (angle_view) -         &
      mu_values (i_loc) ) / (mu_values (i_loc + 1) - mu_values (i_loc) )
      IF (Nstokes.gt.2) then 
         nshift = NSTOKES - 3 
         Field_1 = DOWN_RAD (NSTOKES * nummu * (Noutlevels - 1) +       &
         Nstokes * i_loc - nshift)                                      
         Field_2 = DOWN_RAD (NSTOKES * nummu * (Noutlevels - 1) +       &
         Nstokes * (i_loc + 1) - nshift)                                
         Upar = Field_1 + (Field_2 - Field_1) * (dcos (angle_view)      &
         - mu_values (i_loc) ) / (mu_values (i_loc + 1) - mu_values (   &
         i_loc) )                                                       
         nshift = NSTOKES - 4 
         Field_1 = DOWN_RAD (NSTOKES * nummu * (Noutlevels - 1) +       &
         Nstokes * i_loc - nshift)                                      
         Field_2 = DOWN_RAD (NSTOKES * nummu * (Noutlevels - 1) +       &
         Nstokes * (i_loc + 1) - nshift)                                
         Vpar = Field_1 + (Field_2 - Field_1) * (dcos (angle_view)      &
         - mu_values (i_loc) ) / (mu_values (i_loc + 1) - mu_values (   &
         i_loc) )                                                       
      ENDIF 
                                                                        
      END SUBROUTINE RT3_oneangle                   
                                                                        
                                                                        
      SUBROUTINE READ_LAYERS (LAYER_FILE, MAXLAY, NUM_LAYERS, HEIGHT,   &
      TEMPERATURES, GAS_EXTINCT, SCAT_FILES)       
  use kinds
  
      INTEGER MAXLAY, NUM_LAYERS 
      REAL(kind=dbl) HEIGHT ( * ), TEMPERATURES ( * ) 
      REAL(kind=dbl) GAS_EXTINCT ( * ) 
      CHARACTER ( * ) LAYER_FILE, SCAT_FILES ( * ) 
      INTEGER I 
                                                                        
!           Read in height, temperature, gaseous extinction, and        
!                 scattering file for the layers                        
      OPEN (UNIT = 1, FILE = LAYER_FILE, STATUS = 'OLD') 
      I = 1 
  100 CONTINUE 
      READ (1, *, ERR = 990, END = 110) HEIGHT (I), TEMPERATURES (I),   &
      GAS_EXTINCT (I), SCAT_FILES (I)
      I = I + 1 
      IF (I.EQ.MAXLAY) THEN 
         WRITE ( * , * ) 'Too many layers' 
         STOP 
      ENDIF 
      GOTO 100 
  110 CONTINUE 
      CLOSE (1)!, status='delete')
      NUM_LAYERS = I - 2 
      RETURN 
                                                                        
  990 CONTINUE 
      WRITE ( * , * ) 'Error reading layers data file' 
      RETURN 
      END SUBROUTINE READ_LAYERS                    
                                                                        
                                                                        
      SUBROUTINE USER_INPUT (NSTOKES, NUMMU, AZIORDER, MU_VALUES,       &
      SRC_CODE, LAYER_FILE, OUT_FILE, QUAD_TYPE, DELTAM, DIRECT_FLUX,   &
      DIRECT_MU, GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO, GROUND_INDEX, &
      SKY_TEMP, WAVELENGTH, UNITS, OUTPOL, NOUTLEVELS, OUTLEVELS,       &
      NUMAZIMUTHS)
      
      use kinds
      
      INTEGER NSTOKES, NUMMU, AZIORDER, SRC_CODE 
      INTEGER NOUTLEVELS, OUTLEVELS ( * ), NUMAZIMUTHS
      REAL(kind=dbl) MU_VALUES ( * ), GROUND_TEMP, GROUND_ALBEDO 
      REAL(kind=dbl) SKY_TEMP, WAVELENGTH 
      REAL(kind=dbl) DIRECT_FLUX, DIRECT_MU 
      COMPLEX(kind=dbl) GROUND_INDEX 
      CHARACTER QUAD_TYPE * 1, DELTAM * 1, UNITS * 1, OUTPOL * 2,       &
      GROUND_TYPE * 1                                                   
      CHARACTER ( * ) LAYER_FILE, OUT_FILE 
      REAL(kind=dbl) THETA 
      INTEGER I 
                                                                        
      WRITE ( * , '(1X,A)') 'Number of Stokes parameters (1 - 4) : ' 
      READ ( *, * ) NSTOKES 
      WRITE ( * , '(1X,A)') 'Number of quadrature directions : ' 
      READ ( *, * ) NUMMU 
      WRITE ( * , '(1X,A)') 'Type of quadrature : ' 
      WRITE ( * , '(1X,A)') '(Gaussian, Double-Gauss, Lobatto, Extra-ang&
     &les) : '                                                          
      READ ( * , '(A)') QUAD_TYPE 
      IF (QUAD_TYPE (1:1) .EQ.'E') THEN 
      WRITE ( * ,  * ) 'Enter extra quadrature mu values (end with 0):' 
         I = NUMMU 
   50    CONTINUE 
         WRITE ( * , '(1X,A)') 'Mu value : ' 
         READ ( *, * ) MU_VALUES (I) 
         I = I - 1 
         IF (MU_VALUES (I + 1) .NE.0.0) GOTO 50 
      ENDIF 
                                                                        
      WRITE ( * , '(1X,A)') 'Order of azimuth expansion (0,1,...) : ' 
      READ ( *, * ) AZIORDER 
                                                                        
      WRITE ( * , '(1X,A)') 'Layers data file name : ' 
      READ ( * , '(A)') LAYER_FILE 
                                                                        
      WRITE ( * , '(1X,A)') 'Delta-M scaling (Y or N) : ' 
      READ ( * , '(A)') DELTAM 
                                                                        
      WRITE ( * , '(1X,A)') 'Source code (none=0, solar=1, thermal=2, bo&
     &th=3) : '                                                         
      READ ( *, * ) SRC_CODE 
      SRC_CODE = MIN0 (MAX0 (SRC_CODE, 0), 3) 
                                                                        
      DIRECT_MU = 1.0 
      IF (SRC_CODE.EQ.1.OR.SRC_CODE.EQ.3) THEN 
         WRITE ( * , '(1X,A)') 'Direct flux (W/(m*m)/um or K) : ' 
         READ ( *, * ) DIRECT_FLUX 
      WRITE ( * , '(1X,A)') 'Direct flux direction (zenith angle) (deg) &
     &: '                                                               
         READ ( *, * ) THETA 
         DIRECT_MU = DABS (DCOS (0.017453292D0 * (THETA) ) ) 
      ENDIF 
                                                                        
      WRITE ( * , '(1X,A)') 'Ground temperature : ' 
      READ ( *, * ) GROUND_TEMP 
      WRITE ( * , '(1X,A)') 'Ground type (Lambertian or Fresnel) : ' 
      READ ( * , '(A)') GROUND_TYPE 
      IF (GROUND_TYPE (1:1) .EQ.'F') THEN 
      WRITE ( * , '(1X,A)') 'Complex index of refraction of ground : ' 
         READ ( *, * ) GROUND_INDEX 
      ELSE 
         WRITE ( * , '(1X,A)') 'Ground albedo : ' 
         READ ( *, * ) GROUND_ALBEDO 
      ENDIF 
      WRITE ( * , '(1X,A)') 'Sky temperature : ' 
      READ ( *, * ) SKY_TEMP 
                                                                        
      WRITE ( * , '(1X,A)') 'Wavelength (microns) : ' 
      READ ( *, * ) WAVELENGTH 
      WRITE ( * , '(1X,A)') 'Output radiance units :' 
      WRITE ( * , '(1X,A,A)') '(W-W/m^2 um sr, ', 'T-EBB brightness temp&
     &erature, R-Rayleigh-Jeans Tb) : '                                 
      READ ( * , '(A)') UNITS 
      WRITE ( * , '(1X,A)') 'Output polarization (IQ or VH) : ' 
      READ ( * , '(A)') OUTPOL 
                                                                        
      WRITE ( * , '(1X,A)') 'Number of output levels : ' 
      READ ( *, * ) NOUTLEVELS 
      WRITE ( * , '(1X,A)') 'Output level numbers : ' 
      READ ( *, * ) (OUTLEVELS (I), I = 1, NOUTLEVELS) 
      WRITE ( * , '(1X,A)') 'Number of output azimuths : ' 
      READ ( *, * ) NUMAZIMUTHS 
                                                                        
      WRITE ( * , '(1X,A)') 'Output data file name : ' 
      READ ( * , '(A)') OUT_FILE 
                                                                        
      RETURN 
      END SUBROUTINE USER_INPUT                     
                                                                        
                                                                        
      SUBROUTINE OUTPUT_FILE (NSTOKES, NUMMU, AZIORDER, SRC_CODE,       &
      LAYER_FILE, OUT_FILE, QUAD_TYPE, DELTAM, DIRECT_FLUX, DIRECT_MU,  &
      GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO, GROUND_INDEX, SKY_TEMP,  &
      WAVELENGTH, UNITS, OUTPOL, NUM_LAYERS, HEIGHT, NOUTLEVELS,        &
      OUTLEVELS, NUMAZIMUTHS, MU_VALUES, UP_FLUX, DOWN_FLUX, UP_RAD,    &
      DOWN_RAD,lon,lat,lfrac,wind10,iwv,cwp,iwp,rwp,swp,gwp)

  use kinds
  
      INTEGER NSTOKES, NUMMU, NUMAZI, AZIORDER, SRC_CODE, NUM_LAYERS 
      INTEGER NOUTLEVELS, OUTLEVELS ( * ), NUMAZIMUTHS 
      REAL(kind=dbl) GROUND_TEMP, GROUND_ALBEDO 
      REAL(kind=dbl) SKY_TEMP, WAVELENGTH 
      REAL(kind=dbl) DIRECT_FLUX, DIRECT_MU 
      REAL(kind=dbl) HEIGHT (NUM_LAYERS + 1) 
      REAL(kind=dbl) MU_VALUES (NUMMU) 
      REAL(kind=dbl) UP_FLUX (NSTOKES, NOUTLEVELS) 
      REAL(kind=dbl) DOWN_FLUX (NSTOKES, NOUTLEVELS) 
      REAL(kind=dbl) UP_RAD (NSTOKES, NUMMU, AZIORDER + 1, NOUTLEVELS) 
      REAL(kind=dbl) DOWN_RAD (NSTOKES, NUMMU, AZIORDER + 1, NOUTLEVELS)
      COMPLEX(kind=dbl) GROUND_INDEX
      CHARACTER ( * ) LAYER_FILE, OUT_FILE 
      CHARACTER QUAD_TYPE * 1, DELTAM * 1, UNITS * 1, OUTPOL * 2,       &
      GROUND_TYPE * 1                                                   
      CHARACTER(32) QUAD_NAME, UNITS_NAME, GROUND_NAME 
      CHARACTER(64) FORM1 
      INTEGER I, J, K, L, LI, M, N 
      REAL OUT (4), PHI, PHID, PI 
      PARAMETER (PI = 3.1415926535897932384D0) 
                                                                        
      real lon,lat,lfrac,wind10,iwv,cwp,iwp,rwp,swp,gwp
                                                                        
      N = NUMMU * (AZIORDER + 1) * NOUTLEVELS 
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, N, WAVELENGTH, 0,    &
      UP_RAD)                                                           
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, N, WAVELENGTH, 0,    &
      DOWN_RAD)                                                         
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUTLEVELS,          &
      WAVELENGTH, 1, UP_FLUX)                                           
      CALL CONVERT_OUTPUT (UNITS, OUTPOL, NSTOKES, NOUTLEVELS,          &
      WAVELENGTH, 1, DOWN_FLUX)                                         
                                                                        
      NUMAZI = 2 * AZIORDER + 1 
      IF (NSTOKES.LE.2) NUMAZI = AZIORDER + 1 
      QUAD_NAME = 'GAUSSIAN' 
      IF (QUAD_TYPE.EQ.'D') QUAD_NAME = 'DOUBLEGAUSS' 
      IF (QUAD_TYPE.EQ.'L') QUAD_NAME = 'LOBATTO' 
      IF (QUAD_TYPE.EQ.'E') QUAD_NAME = 'EXTRA-ANGLES' 
      UNITS_NAME = 'WATTS/(M^2 MICRON STER)' 
      IF (UNITS.EQ.'T') UNITS_NAME = 'KELVINS - EBB' 
      IF (UNITS.EQ.'R') UNITS_NAME = 'KELVINS - RJ' 
      GROUND_NAME = 'LAMBERTIAN' 
      IF (GROUND_TYPE.EQ.'F') GROUND_NAME = 'FRESNEL' 
      IF (GROUND_TYPE.EQ.'O') GROUND_NAME = 'OCEAN'
      IF (GROUND_TYPE.EQ.'S') GROUND_NAME = 'SPECULAR' 
                                                                       
      OPEN (UNIT = 3, FILE = OUT_FILE, STATUS = 'UNKNOWN') 
      if (verbose .gt. 1) print *,"writing... ", OUT_FILE                                                     
!           Output the parameters                                       
      WRITE (3, '(A,I3,A,I3,A,I3,A,I1)') 'C  NUMMU=', NUMMU, '  NUMAZI='&
     &, NUMAZI, '  AZIORDER=', AZIORDER, '  NSTOKES=', NSTOKES          
      WRITE (3, '(A,A32,A,A1)') 'C  LAYER_FILE=', LAYER_FILE, '   DELTA-&
     &M=', DELTAM                                                       
      WRITE (3, '(A,I1,A,A16)') 'C  SRC_CODE=', SRC_CODE, '   QUAD_TYPE=&
     &', QUAD_NAME                                                      
      IF (SRC_CODE.EQ.1.OR.SRC_CODE.EQ.3) THEN 
      WRITE (3, '(A,E11.5,A,F8.6)') 'C  DIRECT_FLUX=', DIRECT_FLUX, '   &
     &DIRECT_MU=', DIRECT_MU                                            
      ENDIF 
      write(3,'(A,F8.3,A,F8.3,A,F8.3)') 'C  LON=',lon,' LAT=',lat,&
      ' LFRAC=',lfrac
      WRITE (3, '(A,F8.2,A,A16,A,F8.3)') 'C  GROUND_TEMP=', GROUND_TEMP,&
     & '   GROUND_TYPE=', GROUND_NAME,' W10=',wind10                                    
      IF (GROUND_TYPE (1:1) .EQ.'F' .OR. GROUND_TYPE (1:1) .EQ. 'O' ) THEN
      WRITE (3, '(A,2F9.4,A,F8.2)') 'C  GROUND_INDEX=', GROUND_INDEX, ' &
     &  SKY_TEMP=', SKY_TEMP                                            
      ELSE 
      WRITE (3, '(A,F8.5,A,F8.2)') 'C  GROUND_ALBEDO=', GROUND_ALBEDO, '&
     &   SKY_TEMP=', SKY_TEMP                                           
      ENDIF 
      WRITE (3, '(A,E12.6)') 'C  WAVELENGTH=', WAVELENGTH 
      WRITE (3, '(A,A25,A,A2)') 'C  UNITS=', UNITS_NAME, '   OUTPUT_POLA&
     &RIZATION=', OUTPOL                                                
!      Output integrated quantities
      write(3,'(A,ES12.6)') 'C  IWV= ', iwv
      write(3,'(A,ES12.6)') 'C  CWP= ', cwp
      write(3,'(A,ES12.6)') 'C  IWP= ', iwp
      write(3,'(A,ES12.6)') 'C  RWP= ', rwp
      write(3,'(A,ES12.6)') 'C  SWP= ', swp
      write(3,'(A,ES12.6)') 'C  GWP= ', gwp

      IF (UNITS (1:1) .EQ.'T') THEN 
!         FORM1 = '(F8.3,1X,F5.1,1X,F8.5,4(1X,F7.2),:)' 
         FORM1 = '(F8.1,1X,F5.1,1X,F8.5,4(1X,F7.2),:)' 
      ELSE 
!         FORM1 = '(F8.3,1X,F5.1,1X,F8.5,4(1X,E13.6),:)' 
         FORM1 = '(F8.1,1X,F5.1,1X,F8.5,4(1X,E13.6),:)' 
      ENDIF 
                                                                        
      IF (OUTPOL.EQ.'VH') THEN 
      WRITE (3, '(A)') 'C    Z      PHI     MU    FLUX/RADIANCE (V,H,U,V&
     &)'                                                                
      ELSE 
      WRITE (3, '(A)') 'C    Z      PHI     MU    FLUX/RADIANCE (I,Q,U,V&
     &)'                                                                
      ENDIF 
                                                                        
      DO L = 1, NOUTLEVELS 
      LI = OUTLEVELS (L) 
!               Output fluxes at this level                             
      WRITE (3, FORM1) HEIGHT(LI), 0., - 2.0, (SNGL (UP_FLUX (I, L) ), &
      I = 1, NSTOKES)                                                   
      WRITE (3, FORM1) HEIGHT(LI), 0., + 2.0, (SNGL (DOWN_FLUX (I, L) )&
      , I = 1, NSTOKES)                                                 
                                                                        
!    For each azimuth and zenith at this level sum the Fourier
!    azimuth series appropriate for the particular Stokes parameter
!    and output the radiance.                                
      DO K = 1, NUMAZIMUTHS 
      IF (NUMAZIMUTHS.EQ.1) THEN 
         PHID = 0.0 
      ELSE 
         PHID = 180.0 * FLOAT (K - 1) / (NUMAZIMUTHS - 1) 
      ENDIF 
      PHI = PI * PHID / 180.0 
!               Output upwelling radiance: -1 < mu < 0                  
      DO J = NUMMU, 1, - 1 
      DO I = 1, NSTOKES 
      OUT (I) = 0.0 
      DO M = 0, AZIORDER 
      IF (I.LE.2) THEN 
         OUT (I) = OUT (I) + COS (M * PHI) * UP_RAD (I, J, M + 1, L) 
      ELSE 
         OUT (I) = OUT (I) + SIN (M * PHI) * UP_RAD (I, J, M + 1, L) 
      ENDIF 
      ENDDO 
      ENDDO 
      WRITE (3, FORM1) HEIGHT (LI), PHID, - MU_VALUES (J), (OUT (I),    &
      I = 1, NSTOKES)                                                   
      ENDDO 
!               Output downwelling radiance: 0 < mu < 1                 
      DO J = 1, NUMMU 
      DO I = 1, NSTOKES 
      OUT (I) = 0.0 
      DO M = 0, AZIORDER 
      IF (I.LE.2) THEN 
         OUT (I) = OUT (I) + COS (M * PHI) * DOWN_RAD (I, J, M + 1, L) 
      ELSE 
         OUT (I) = OUT (I) + SIN (M * PHI) * DOWN_RAD (I, J, M + 1, L) 
      ENDIF 
      ENDDO 
      ENDDO 
      WRITE (3, FORM1) HEIGHT (LI), PHID, MU_VALUES (J), (OUT (I),      &
      I = 1, NSTOKES)                                                   
      ENDDO 
      ENDDO 
      ENDDO 
                                                                        
      CLOSE (3) 
                                                                        
      RETURN 
      END SUBROUTINE OUTPUT_FILE                    
                                                                        
                                                                        
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
