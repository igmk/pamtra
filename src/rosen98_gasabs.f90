!+ gaseous absorption after Rosenkranz 98 model
subroutine rosen98_gasabs &
 (errorstatus,&     ! out
  freq,  &      ! in
  tempK, &      ! in
  rhoWv, &      ! in
  pres,  &      ! in
  absAir,&      ! out
  absWv  )      ! out

! Description:
!  Based on frequency, temperature, water vapor density, and pressure, this routine
!  calculates the absorption due to air (N2 and O2) and water vapor in the frequency
!  range from 0.1 to 800.0 GHz, pressure range from 10 to 1.2e5 Pa, and absolute
!  temperatures larger than 100 K.
!
! Method:
!  ROSENKRANZ (1998) model  -- reference "Water vapor microwave
!  continuum absorption: a comparison of measurements and results"
!  To appear in Radio Science
!
! Owner: IGMK
!
! History:
!
! Version Date         Comment
! ------- ----         -------
! 0.01    02/08/2001   Fixed division by zero when rhowv = 0 - G. Petty
! 0.1     21/09/2009   Code adaption from G. Petty - M. Mech
! 0.2     27/10/2012   Corrections to WV bands and continuum - E. Orlandi
! 0.3     13/11/2012   Application of European Standards for Writing and
!                      Documenting Exchangeable Fortran 90 Code - M. Mech
! 0.4     27/02/2013   functions in module gasabs_module
!
! Code Description:
!  Language: Fortran 90.
!  Software Standards: "European Standards for Writing and
!   Documenting Exchangeable Fortran 90 Code".
!
! Parent Module: get_gasabs
!
! Declarations:

! Modules used:

  use kinds, only: dbl, &     ! integer parameter specifying double precision
                    long       ! integer parameter specifying long integer

  use gasabs_module            ! functions for calculating absorption by gases

!  use settings, only: verbose

  use report_module
! Imported Scalar Variables with intent (in):

  implicit none

!- End of header ---------------------------------------------------------------

! Subroutine arguments
! Scalar arguments with intent(in):

  real(kind=dbl), intent(in) :: freq     ! frequency [GHz]
  real(kind=dbl), intent(in) :: tempK    ! temperature [K]
  real(kind=dbl), intent(in) :: pres     ! pressure [Pa]
  real(kind=dbl), intent(in) :: rhoWv    ! water vapor density [kg/m**3]

! Scalar arguments with intent(out):

  real(kind=dbl), intent(out) :: absAir  ! extinction by dry air [Np/km]
  real(kind=dbl), intent(out) :: absWv   ! extinction by water vapor [Np/km]

! End of Subroutine arguments

! Local scalars:

  real(kind=dbl) :: pmb         ! pressure [mb]
  real(kind=dbl) :: vapden      ! water vpor density [g/m**3]
  real(kind=dbl) :: e           ! water vapor pressure [Pa]
  real(kind=dbl) :: q           ! specific humidity
  real(kind=dbl) :: Tv          ! virtuel temperature [K]
  real(kind=dbl) :: rhoair      ! moist air density [kg/m**3]

! Used Functions

!  real(kind=dbl) :: absn2       ! function to calculate extinction by n_2
!  real(kind=dbl) :: o2abs       ! function to calculate extinction by o2
!  real(kind=dbl) :: abh2o       ! function to calculate extinction by h_2o


! Error handling

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'rosen98_gasabs'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  ! check for "reasonable" input values

!  print*, freq
  if ((freq <= 0.0_dbl) .or. (freq > 800.0_dbl)) then
     errorstatus = fatal
     msg = 'Frequency not between 0 and 800 GHz in rosen98_gasabs!'
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif (tempK <= 100.0_dbl) then
     errorstatus = fatal
     msg = 'Temperature lower than 100 K in rosen98_gasabs!'
     call report(errorstatus, msg, nameOfRoutine)
     return
  elseif ((pres < 10.0_dbl) .or. (pres > 1.2d5)) then
     print*, pres
     errorstatus = fatal
     msg = 'Pressure not between 10 and 1.2d5 Pa in rosen98_gasabs!'
     call report(errorstatus, msg, nameOfRoutine)
     return
  else
     err = success
  end if

  ! convert pressure from Pa to Mb
  pmb = pres / 100.0_dbl

  ! convert vapor density from kg/m**3 to g/m**3

  vapden = rhoWv * 1000.0_dbl

  ! get volume extinction coefficients
  absair = absn2(tempK,pmb,freq) + o2abs(tempK,pmb,vapden,freq)
  abswv = abh2o(tempK,pmb,vapden,freq)

  ! convert vapor density to vapor pressure
  e = rhoWv * (tempK * 461.5_dbl)
  ! calculate specific humidity
  q = 0.622_dbl * e / pres
  ! calculate virtual temperature
  Tv = (1._dbl + 0.61_dbl * q) * tempK
  ! moist air density
  rhoair = pres / (Tv * 287.06_dbl)

  errorstatus = err

  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

  return

end subroutine rosen98_gasabs
