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
subroutine rt4(errorstatus,out_file,&
     ground_type,ground_albedo,ground_index,sky_temp,&
     wavelength,outlevels)

  use kinds
  use vars_atmosphere
  use vars_rt, only : &
       rt_kextatmo,&
       rt_hydros_present, &
       rt_scattermatrix, &
       rt_extmatrix, &
       rt_emisvec, &
       rt_hydros_present_reverse, &
       rt_scattermatrix_reverse, &
       rt_extmatrix_reverse, &
       rt_emisvec_reverse
  use vars_index, only: i_x,i_y,i_f

  use settings, only: write_nc,&
       in_python,&
       numazimuths,&
       verbose,&
       nstokes,&
       nummu,&
       mu_values,&
       quad_type,&
       noutlevels,&
       outpol,&
       units

  use report_module

  implicit none

  integer, parameter :: maxv=64, maxlay=300

  !    real(kind=dbl), intent(in) ::  mu_values(maxv)
  character*64, intent(in) :: out_file
  character, intent(in) ::  ground_type*1
  real(kind=dbl), intent(in) ::  ground_albedo
  complex*16, intent(in) ::   ground_index
  real(kind=dbl), intent(in) ::  sky_temp
  real(kind=dbl), intent(in) ::   wavelength
  integer, intent(in), dimension(maxlay) ::  outlevels

  integer :: num_layers

  real(kind=dbl) :: ground_temp
  real(kind=dbl) :: max_delta_tau
  real(kind=dbl), dimension(maxlay) :: height, temperatures
  real(kind=dbl), dimension(maxlay) :: gas_extinct
  real(kind=dbl), dimension(maxv*(maxlay+1)) :: up_rad, down_rad
  real(kind=dbl), dimension(4*(maxlay+1)) :: up_flux, down_flux

  character(len=64) :: layer_file


  ! error handling

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'rt4'

  err = 0

  if (verbose >= 1) call report(info, 'start of ', nameOfRoutine)

  height = 0.
  temperatures = 0.
  gas_extinct = 0.
  max_delta_tau = 1.0d-6

  layer_file=""

  !scat_files = ''
  !scat_files(atmo_nlyrs(i_x,i_y)) = '1.txt'


  ground_temp = atmo_groundtemp(i_x,i_y)
  num_layers = atmo_nlyrs(i_x,i_y)
  height(1:atmo_nlyrs(i_x,i_y)+1) = atmo_hgt_lev(i_x,i_y,atmo_nlyrs(i_x,i_y)+1:1:-1)             ! [m]
  temperatures(1:atmo_nlyrs(i_x,i_y)+1) = atmo_temp_lev(i_x,i_y,atmo_nlyrs(i_x,i_y)+1:1:-1)      ! [k]
  gas_extinct(1:atmo_nlyrs(i_x,i_y)) = rt_kextatmo(atmo_nlyrs(i_x,i_y):1:-1)         ! [np/m]

  !do some tests
  call assert_true(err,(maxlay>=num_layers),&
       "maxlay>=num_layers")  
  call assert_true(err,(ground_temp>1),&
       "ground_temp must be greater 1")  
  call assert_true(err,(num_layers>0),&
       "num_layers must be greater 0")   
  call assert_true(err,all(height(1:atmo_nlyrs(i_x,i_y)+1)>=-370),&
       "height must be greater -370 (depth of tagebau hambach :-))")  
  call assert_true(err,all(temperatures(1:atmo_nlyrs(i_x,i_y)+1)>1),&
       "temperatures must be greater 1")   
  call assert_false(err,all(isnan(gas_extinct(1:atmo_nlyrs(i_x,i_y)+1))),&
       "gas_extinct must be non nan zero")       
  if (err > 0) then
     errorstatus = fatal
     msg = "assertation error"
     call report(errorstatus, msg, nameOfRoutine)
     return
  end if


  rt_hydros_present_reverse(1:atmo_nlyrs(i_x,i_y)) = rt_hydros_present(atmo_nlyrs(i_x,i_y):1:-1)

  rt_scattermatrix_reverse(1:atmo_nlyrs(i_x,i_y),:,:,:,:,:) = rt_scattermatrix(atmo_nlyrs(i_x,i_y):1:-1,:,:,:,:,:)
  rt_extmatrix_reverse(1:atmo_nlyrs(i_x,i_y),:,:,:,:) = rt_extmatrix(atmo_nlyrs(i_x,i_y):1:-1,:,:,:,:)
  rt_emisvec_reverse(1:atmo_nlyrs(i_x,i_y),:,:,:) = rt_emisvec(atmo_nlyrs(i_x,i_y):1:-1,:,:,:)

  !  if (verbose .gt. 0) print*, ".... read_layers done!"

  call radtran4(err,max_delta_tau,&
       ground_temp, ground_type,&
       ground_albedo, ground_index,&
       sky_temp, wavelength,&
       num_layers, height, temperatures,&
       gas_extinct,&
       outlevels,&
       up_flux, down_flux,&
       up_rad, down_rad)
       
  if (err /= 0) then
     msg = 'error in radtran4'
     call report(err,msg, nameOfRoutine)
     errorstatus = err
     return
  end if

  !  if (verbose .gt. 0) print*, ".... radtran done!"

  if (write_nc .or. in_python) then
     call collect_output(nstokes, nummu, 0, &
          wavelength,   &
          units, outpol,noutlevels, outlevels,         &
          numazimuths,up_rad, down_rad,     &
          i_x,i_y,i_f)
  else
     call output_file4(nstokes, nummu,&
          layer_file, out_file,&
          quad_type, ground_temp, ground_type,&
          ground_albedo, ground_index,&
          sky_temp, wavelength, units, outpol,&
          num_layers, height,&
          noutlevels, outlevels, &
          mu_values, up_flux, down_flux,&
          up_rad, down_rad)
  end if

  errorstatus = err

  if (verbose >= 1) call report(info, 'end of ', nameOfRoutine)
  return

end subroutine rt4


subroutine read_layers (layer_file, maxlay, num_layers,&
     height, temperatures,&
     gas_extinct, scat_files)
  use kinds

  implicit none

  integer  maxlay, num_layers
  real(kind=dbl)   height(*), temperatures(*)
  real(kind=dbl)   gas_extinct(*)
  character*(*)  layer_file, scat_files(*)
  integer   i

  !           read in height, temperature, gaseous extinction, and
  !                 scattering file for the layers
  open (unit=1, file=layer_file, status='old')
  i = 1
100 continue
  read (1,*,err=990,end=110) height(i), temperatures(i),&
       gas_extinct(i), scat_files(i)
  i = i + 1
  if (i .eq. maxlay) then
     write (*,*) 'too many layers'
     stop
  endif
  goto 100
110 continue
  close(1)
  num_layers = i - 2
  return

990 continue
  write (*,*) 'error reading layers data file'

  return
end subroutine read_layers

subroutine output_file4(nstokes, nummu, &
     layer_file, out_file,&
     quad_type, ground_temp, ground_type,&
     ground_albedo, ground_index,&
     sky_temp, wavelength, units, outpol,&
     num_layers, height,&
     noutlevels, outlevels, &
     mu_values, up_flux, down_flux,&
     up_rad, down_rad)

  use kinds

  implicit none

  integer  nstokes, nummu, num_layers
  integer  noutlevels, outlevels(*)
  real(kind=dbl)   ground_temp, ground_albedo
  real(kind=dbl)   sky_temp, wavelength
  real(kind=dbl)   mu_values(nummu)
  real(kind=dbl)   height(num_layers+1)
  real(kind=dbl)   up_flux(nstokes,noutlevels)
  real(kind=dbl)   down_flux(nstokes,noutlevels)
  real(kind=dbl)   up_rad(nstokes,nummu,noutlevels)
  real(kind=dbl)   down_rad(nstokes,nummu,noutlevels)
  complex*16  ground_index
  character*(*) layer_file, out_file
  character  quad_type*1, units*1, outpol*2, ground_type*1
  character*32 quad_name, units_name, ground_name
  character*64 form1
  integer  i, j, l, li


  call convert_output4(units, outpol, nstokes, nummu*noutlevels, &
       wavelength, 0, up_rad)
  call convert_output4(units, outpol, nstokes, nummu*noutlevels, &
       wavelength, 0, down_rad)
  call convert_output4(units, outpol, nstokes, noutlevels, &
       wavelength, 1, up_flux)
  call convert_output4(units, outpol, nstokes, noutlevels, &
       wavelength, 1, down_flux)

  quad_name = 'gaussian'
  if (quad_type .eq. 'd')  quad_name = 'doublegauss'
  if (quad_type .eq. 'l')  quad_name = 'lobatto'
  if (quad_type .eq. 'e')  quad_name = 'extra-angles'
  units_name = 'watts/(m^2 micron ster)'
  if (units .eq. 't') units_name = 'kelvins - ebb'
  if (units .eq. 'r') units_name = 'kelvins - rj'
  ground_name = 'lambertian'
  if (ground_type .eq. 'f')  ground_name = 'fresnel'

  open (unit=3, file=out_file, status='unknown')

  !           output the parameters
  write (3,'(a,i3,a,i3,a,i3,a,i1)')&
       'c  nummu=', nummu,  '  numazi=',1,&
       '  aziorder=',0, '  nstokes=',nstokes
  write (3,'(a,a32)')&
       'c  layer_file=',    layer_file
  write (3,'(a,i1,a,a16)')&
       'c  src_code=',      2,&
       '   quad_type=',     quad_name
  write (3,'(a,f8.2,a,a16)')&
       'c  ground_temp=',   ground_temp,&
       '   ground_type=',   ground_name
  if (ground_type(1:1) .eq. 'f') then
     write (3,'(a,2f9.4,a,f8.2)')&
          'c  ground_index=',  ground_index,&
          '   sky_temp=',      sky_temp
  else
     write (3,'(a,f8.5,a,f8.2)')&
          'c  ground_albedo=', ground_albedo,&
          '   sky_temp=',      sky_temp
  endif
  write (3,'(a,e12.6)') 'c  wavelength=',    wavelength
  write (3,'(a,a25,a,a2)') 'c  units='     ,    units_name,&
       '   output_polarization=', outpol


  if (units(1:1) .eq. 't') then
     form1 = '(f8.1,1x,f8.5,2(1x,f7.2),:)'
  else
     form1 = '(f8.1,1x,f8.5,2(1x,e13.6),:)'
  endif

  if (outpol .eq. 'vh') then
     write (3,'(a,a)') 'c    z       mu    flux/radiance (v,h)'
  else
     write (3,'(a,a)') 'c    z       mu    flux/radiance (i,q)'
  endif

  do l = 1, noutlevels
     li = outlevels(l)
     !               output fluxes at this level
     write (3,form1) height(li), -2.0,&
          (sngl(up_flux(i,l)),i=1,nstokes)
     write (3,form1) height(li), +2.0,&
          (sngl(down_flux(i,l)),i=1,nstokes)

     !           for each zenith at this level output the stokes parameters.
     !             output upwelling radiance: -1 < mu < 0
     do j = nummu, 1, -1
        write (3,form1) height(li), -mu_values(j),&
             (sngl(up_rad(i,j,l)),i=1,nstokes)
     enddo
     !             output downwelling radiance: 0 < mu < 1
     do j = 1, nummu
        write (3,form1) height(li), mu_values(j),&
             (sngl(down_rad(i,j,l)),i=1,nstokes)
     enddo
  enddo

  close (3)

  return
end subroutine output_file4

subroutine convert_output4(units, outpol, nstokes, nout,&
     wavelen, fluxcode, output)
  !       converts the output radiance or flux arrays to vh polarization
  !     and effective blackbody temperature if desired.  outpol='vh'
  !     converts the polarization basis of the first two stokes parameters
  !     to vertical/horizontal polarization.  if units='t' the radiance is
  !     converted to effective blackbody brightness temperature, and if
  !     units='r' the radiance is converted to rayleigh-jeans brightness
  !     temperature.  if the output is flux then fluxcode=1, and the flux
  !     is divided by pi before converting to brightness temperature.
  use kinds

  implicit none

  integer nstokes, nout, fluxcode
  real(kind=dbl)  wavelen, output(nstokes,nout)
  character units*1, outpol*2
  integer i, j
  real(kind=dbl)  iv, ih, rad, temp

  do j = 1, nout
     !           convert to vertical and horizontal polarization if desired
     if (outpol .eq. 'vh') then
        iv = 0.5*(output(1,j) + output(2,j))
        ih = 0.5*(output(1,j) - output(2,j))
        output(1,j) = iv
        output(2,j) = ih
     endif
     !           convert to brightness temperature
     if (units .eq. 't' .or. units .eq. 'r') then
        do i = 1, nstokes
           rad = output(i,j)
           if (outpol .eq. 'vh' .and. i .le. 2)  rad = 2.0*rad
           if (fluxcode .eq. 1)  rad = rad/acos(-1.0)
           if (units .eq. 'r') then
              temp = rad * wavelen**4 * 1.4388d4/1.1911d8
           else
              if (rad .gt. 0.0) then
                 temp = 1.4388d4 /&
                      (wavelen*dlog(1.0+ 1.1911d8/(rad*wavelen**5)))
              else if (rad .eq. 0.0) then
                 temp = 0.0d0
              else
                 temp = -1.4388d4 /&
                      (wavelen*dlog(1.0+ 1.1911d8/(-rad*wavelen**5)))
              endif
           endif
           output(i,j) = temp
        enddo
     endif
  enddo
  return
end subroutine convert_output4



