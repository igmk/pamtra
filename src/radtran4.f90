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

subroutine radtran4(errorstatus, max_delta_tau,&
     ground_temp, ground_type,&
     ground_albedo,&
     sky_temp, wavelength,&
     num_layers, height, temperatures,&
     gas_extinct,&
     outlevels,&
     up_flux, down_flux,&
     up_rad, down_rad)

  use kinds
  use vars_rt, only : &
       rt_hydros_present_reverse
  use settings, only: verbose, maxlay, &
    nstokes, nummu, mu_values, quad_weights, noutlevels
  use report_module
  use rt_utilities, only: planck_function!,&
  !gauss_legendre_quadrature,&
  !double_gauss_quadrature,&
  !lobatto_quadrature
  use sfc_matrices, only: get_sfc_matrices

  implicit none

  integer   num_layers
  integer   outlevels(*)
  real*8    ground_temp, ground_albedo
  real*8    sky_temp
  real*8    wavelength, max_delta_tau
  real*8    height(*), temperatures(*)
  real*8    gas_extinct(*)
  !      real*8    mu_values(*)
  real*8    up_flux(*), down_flux(*)
  real*8    up_rad(*), down_rad(*)
  character*1  ground_type

  integer   maxv, maxm, maxlm
  parameter (maxv=64, maxm=4096, maxlm=2001 * (maxv)**2)!maxlm=201*256) !maxlm=201 * (maxv)**2

  real*8    pi, twopi, zero
  parameter (pi = 3.1415926535897932384d0, twopi=2.0d0*pi)
  parameter (zero=0.0d0)

  integer   layer, num_doubles
  integer   i, j, k, l, n, krt, ks
  logical   symmetric
  real*8    linfactor
  real*8    planck0, planck1
  real*8    zdiff, delta_z, f, num_sub_layers, extinct
  !      real*8    quad_weights(maxv)
  real*8    scatter_matrix(4*maxm)
  real*8    lin_source(2*maxv)
  real*8    extinct_matrix(16*2*maxv), emis_vector(4*2*maxv)
  real*8    reflect1(2*maxm),upreflect(2*maxm),downreflect(2*maxm)
  real*8    trans1(2*maxm),  uptrans(2*maxm),  downtrans(2*maxm)
  real*8    source1(2*maxv), upsource(2*maxv), downsource(2*maxv)
  real*8    reflect(2*maxlm)
  real*8    trans(2*maxlm)
  real*8    source(2*maxv*(maxlay+1))
  real*8    gnd_radiance(maxv), sky_radiance(2*maxv)
  character*64 scat_file

  !real(kind=dbl), parameter :: quadcof(4,2) =      & ! (unused)
  !     & reshape((/0.0d0, 1.0d0, 1.0d0, 2.0d0, 1.0d0, -1.0d0, 1.0d0, -1.0d0/), (/4, 2/))

  ! variables needed for fastem4

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'radtran4'

  !Without setting these values to zero, PyPamtra results might when called several times!
  source(:) = 0.d0
  reflect(:) = 0.d0
  trans(:) = 0.d0

  err = 0
  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)


  if (verbose >= 99) print*, "RADTRAN4"
  if (verbose >= 99) print*, max_delta_tau,&
     ground_temp, ground_type,&
     ground_albedo,&
     sky_temp, wavelength,&
     num_layers, sum(height(1:num_layers)), sum(temperatures(1:num_layers)),&
     SUM(gas_extinct(1:num_layers)),&
     outlevels(1:noutlevels)

  symmetric = .true.
  n = nstokes*nummu
  if (n .gt. maxv) then
     write (*,'(1x,a,i3)')&
          'vector size exceeded.  maximum size :', maxv
     msg = 'radtran check'
     call report(err,msg, nameOfRoutine)
     errorstatus = fatal
     return
  end if
  if (n*n .gt. maxm) then
     write (*,'(1x,a,i3)')&
          'matrix size exceeded.  maximum size :', maxm
     msg = 'radtran check'
     call report(err,msg, nameOfRoutine)
     errorstatus = fatal
     return
  end if
  if (num_layers .gt. maxlay) then
     write (*,'(1x,a,a,i3)') 'number of layers exceeded.',&
          '  maximum number :', maxlay
     msg = 'radtran check'
     call report(err,msg, nameOfRoutine)
     errorstatus = fatal
     return
  end if
  if ((num_layers+1)*n*n .gt. maxlm) then
     write (*,'(1x,a,a,i3)') 'matrix layer size exceeded.',&
          '  maximum number :', maxlm
     msg = 'radtran check'
     call report(err,msg, nameOfRoutine)
     errorstatus = fatal
     return
  end if

  ! !           make the desired quadrature abscissas and weights
  !       if (quad_type(1:1) .eq. 'd') then
  !         call double_gauss_quadrature&
  !                             (nummu, mu_values, quad_weights)
  !       else if (quad_type(1:1) .eq. 'l') then
  !         call lobatto_quadrature&
  !                             (nummu, mu_values, quad_weights)
  !       else if (quad_type(1:1) .eq. 'e') then
  !         j = nummu
  !         do i = nummu, 1, -1
  !           if (mu_values(i) .ne. 0.0) then
  !             quad_weights(i) = 0.0
  !             j = i - 1
  !           endif
  !         enddo
  !         call gauss_legendre_quadrature&
  !                             (j, mu_values, quad_weights)
  !       else
  !         call gauss_legendre_quadrature&
  !                             (nummu, mu_values, quad_weights)
  !       endif

  scat_file = '&&&'
  !     ------------------------------------------------------
  !           loop through the layers
  !              do doubling to make the reflection and transmission matrices
  !              and soure vectors for each layer, which are stored.

  do layer = 1, num_layers
     !                   calculate the layer thickness
     zdiff = abs(height(layer) - height(layer+1))
     gas_extinct(layer) = max(gas_extinct(layer),0.0d0)
     !        if (rt4kexttot(layer) .gt. 0.0) then

     if (rt_hydros_present_reverse(layer)) then
        call get_scat_mat(layer,nstokes, nummu,scatter_matrix,extinct_matrix, emis_vector)

        call check_norm4(err,nstokes, nummu, quad_weights,&
             scatter_matrix,&
             extinct_matrix, emis_vector)
     end if
     if (err /= 0) then
        msg = 'error in check_norm4'
        call report(err,msg, nameOfRoutine)
        errorstatus = err
        return
     end if
     !          if (scat_files(layer) .ne. scat_file .and. scat_files(layer) .ne. ' ')  then
     !              scat_file = scat_files(layer)
     !       read the scattering matrix from the file
     !              call get_scat_file(nstokes, nummu, quad_type, scat_file,&
     !                                 scatter_matrix,extinct_matrix, emis_vector)
     !              call check_norm4(nstokes, nummu, quad_weights,&
     !                              scatter_matrix,&
     !                              extinct_matrix, emis_vector)
     !          endif

     !                   do the stuff for thermal source in layer
     !                   calculate the thermal source for end of layer
     call planck_function (temperatures(layer+1), 'r',wavelength, planck1)
     !                   calculate the thermal source for beginning of layer
     call planck_function (temperatures(layer), 'r',wavelength, planck0)

     krt = 1 + 2*n*n*(layer-1)
     ks = 1 + 2*n*(layer-1)

     !          if (scat_files(layer) .eq. ' ') then
     !          if (rt4kexttot(layer) .eq. 0.d0) then
     if (.not. rt_hydros_present_reverse(layer)) then
        !                   if the layer is purely absorbing then quickly
        !                     make the reflection and transmission matrices
        !                     and source vector instead of doubling.
        call nonscatter_layer4(nstokes, nummu, &
             zdiff*gas_extinct(layer), mu_values,&
             planck0, planck1,&
             reflect(krt), trans(krt), source(ks))
     else

        !                   find initial thickness of sublayer and
        !                     the number of times to double
        extinct = extinct_matrix(1)+gas_extinct(layer)
        !     print*, gas_extinct(layer), extinct_matrix(1), extinct
        f = dlog(max(extinct*zdiff,1.0d-7)/max_delta_tau)/log(2.)
        num_doubles = 0
        if (f .gt. 0.0)  num_doubles = int(f) + 1
        num_sub_layers = 2.0**num_doubles
        delta_z = zdiff / num_sub_layers

        !                   initialize the source vector
        call initial_source4(nstokes, nummu, delta_z, mu_values,&
             planck0, emis_vector, gas_extinct(layer),&
             lin_source)
        if (planck0 .eq. 0.0) then
           linfactor = 0.0
        else
           linfactor = (planck1/planck0-1.0d0) /num_sub_layers
        endif

        !                generate the local reflection and transmission matrices
        call initialize4(nstokes, nummu, &
             delta_z, mu_values, quad_weights,&
             gas_extinct(layer), extinct_matrix,&
             scatter_matrix,  reflect1, trans1)

        !                   double up to the thickness of the layer
        call doubling_integration4(n, num_doubles, symmetric,&
             reflect1, trans1, lin_source, linfactor,&
             reflect(krt), trans(krt), source(ks))
     endif
  enddo
  !            end of layer loop

  !           get the surface reflection and transmission matrices
  !             and the surface radiance
  krt = 1 + 2*n*n*(num_layers)
  ks = 1 + 2*n*(num_layers)

  if (verbose > 4) print*, "calculating surface emissivity ...."

  call get_sfc_matrices(err,ground_type, ground_temp, reflect(krt), trans(krt), gnd_radiance)

  if (verbose > 4) print*, ".... done!"

  !           assume the radiation coming from above is blackbody radiation
  ! 0 stands for mode = 0. this is required, since we use the routine from the former
  ! radutil3.f
  call thermal_radiance (nstokes, nummu,0, sky_temp, zero,  &
       wavelength,  sky_radiance)

  !         for each desired output level (1 thru nl+2) add layers
  !           above and below level and compute internal radiance.
  !           outlevels gives the desired output levels.
  do i = 1, noutlevels
     layer = min( max( outlevels(i), 1), num_layers+2)
     upreflect(:) = 0.D0
     downreflect(:) = 0.D0
     call midentity (n, uptrans(1))
     call midentity (n, uptrans(1+n*n))
     call midentity (n, downtrans(1))
     call midentity (n, downtrans(1+n*n))
     upsource(:) = 0.D0
     downsource(:) = 0.D0

     do l = 1, layer-1
        krt = 1 + 2*n*n*(l-1)
        ks = 1 + 2*n*(l-1)
        if (l .eq. 1) then
           call mcopy (2*n,n, reflect(krt), upreflect)
           call mcopy (2*n,n, trans(krt), uptrans)
           call mcopy (2*n,1, source(ks), upsource)
        else
           call mcopy (2*n,n, upreflect, reflect1)
           call mcopy (2*n,n, uptrans, trans1)
           call mcopy (2*n,1, upsource, source1)
           call combine_layers4(n, reflect1, trans1, source1,&
                reflect(krt), trans(krt), source(ks),&
                upreflect, uptrans, upsource)
        endif
     enddo
     do l = layer, num_layers+1
        krt = 1 + 2*n*n*(l-1)
        ks = 1 + 2*n*(l-1)
        if (l .eq. layer) then
           call mcopy (2*n,n, reflect(krt), downreflect)
           call mcopy (2*n,n, trans(krt), downtrans)
           call mcopy (2*n,1, source(ks), downsource)
        else
           call mcopy (2*n,n, downreflect, reflect1)
           call mcopy (2*n,n, downtrans, trans1)
           call mcopy (2*n,1, downsource, source1)
           call combine_layers4(n, reflect1, trans1, source1,&
                reflect(krt), trans(krt), source(ks),&
                downreflect, downtrans, downsource)
        endif
     enddo
     call internal_radiance4(n, upreflect, uptrans, upsource,&
          downreflect, downtrans, downsource, &
          sky_radiance, gnd_radiance,&
          up_rad(1+(i-1)*n), down_rad(1+(i-1)*n))
  enddo



  !           integrate the mu times the radiance to find the fluxes
  do l = 1, noutlevels
     do i = 1, nstokes
        k = i+nstokes*(l-1)
        up_flux(k) = 0.0
        down_flux(k) = 0.0
        do j = 1, nummu
           up_flux(k) = up_flux(k)&
                + twopi*quad_weights(j) * mu_values(j)&
                * up_rad(i+nstokes*(j-1)+n*(l-1))
           down_flux(k) = down_flux(k)&
                + twopi*quad_weights(j) * mu_values(j)&
                * down_rad(i+nstokes*(j-1)+n*(l-1))
        enddo
     enddo
  enddo

  errorstatus = err
  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

  return
end subroutine radtran4
