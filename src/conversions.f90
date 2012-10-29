module conversions

contains

real function q2abs(spec_var,t,p,qv,qc,qi,qr,qs,qg,qh)

  use kinds

  implicit none

  real(kind=dbl), parameter :: r_d = 287.05d0,& ! gas constant of dry air
							   r_v = 461.5d0    ! gas constant of water vapor

  real(kind=dbl), intent(in) :: t,&
  								p

  real(kind=dbl), intent(in) :: spec_var,& ! specific variable to convert
  								qv,&
  								qc,&
  								qi,&
  								qr,&
  								qs,&
  								qg

  real(kind=dbl), intent(in), optional :: qh

  if (present(qh)) then
	q2abs = spec_var*p/(r_d*(1.d0+(r_v/r_d-1.d0)*qv-qs-qc-qi-qr-qg-qh)*t)
  else
	q2abs = spec_var*p/(r_d*(1.d0+(r_v/r_d-1.d0)*qv-qs-qc-qi-qr-qg)*t)
  end if


end function q2abs

real function vapor2rh(temp_p,pres_p,hum_massmix)

  implicit none
  !     XPABSM air pressure in Pa
  !     ZTEMP air temperature in K
  !     XRM water vapor mass mixing ratio kg/kg

  real, intent(in) :: temp_p,pres_p,hum_massmix

  REAL :: XMD,XMV            ! Molar mass of dry air and molar mass of vapor
  REAL :: XBOLTZ             ! Boltzman constant
  REAL :: XAVOGADRO          ! Avogadro number
  REAL :: XRD,XRV            ! Gaz constant for dry air, gaz const for vapor
  REAL :: XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
  REAL :: XCL,XCI            ! Cl (liquid), Ci (ice)
  REAL :: XTT                ! Triple point temperature
  REAL :: XLVTT              ! Vaporization heat constant
  REAL :: XLSTT              ! Sublimation heat constant
  REAL :: XLMTT              ! Melting heat constant
  REAL :: XESTT              ! Saturation vapor pressure triple point temp.
  REAL :: XALPW,XBETAW,XGAMW ! Const saturation vapor pressure  (liquid)
  REAL :: XALPI,XBETAI,XGAMI ! Consts saturation vapor pressure  (solid ice)
  real :: ztemp,xrm,xpabsm,zwork31,zwork32

!  print*, temp_p, pres_p, hum_massmix
  XBOLTZ      = 1.380658E-23
  XAVOGADRO   = 6.0221367E+23
  XMD    = 28.9644E-3
  XMV    = 18.0153E-3
  XRD    = XAVOGADRO * XBOLTZ / XMD
  XRV    = XAVOGADRO * XBOLTZ / XMV
  XCPV   = 4.* XRV
  XCL    = 4.218E+3
  XCI    = 2.106E+3
  XTT    = 273.16
  XLVTT  = 2.5008E+6
  XLSTT  = 2.8345E+6
  XLMTT  = XLSTT - XLVTT
  XESTT  = 611.14
  XGAMW  = (XCL - XCPV) / XRV
  XBETAW = (XLVTT/XRV) + (XGAMW * XTT)
  XALPW  = LOG(XESTT) + (XBETAW /XTT) + (XGAMW *LOG(XTT))
  XGAMI  = (XCI - XCPV) / XRV
  XBETAI = (XLSTT/XRV) + (XGAMI * XTT)
  XALPI  = LOG(XESTT) + (XBETAI /XTT) + (XGAMI *LOG(XTT))

  ZTEMP=temp_p
  XRM=hum_massmix
  XPABSM=pres_p

  !     * humidit351 relative par rapport 340 l'eau liquide

  if (ZTEMP.ge.XTT) then

     ZWORK31=EXP(XALPW - XBETAW/ZTEMP - XGAMW*ALOG(ZTEMP))
     ZWORK31=(XMV/XMD)*ZWORK31/(XPABSM-ZWORK31)
     ZWORK32=100.*XRM/ZWORK31

  end if

  !     * humidit351 relative par rapport 340 la glace

  if (ZTEMP.lt.XTT) then
     ZWORK31= EXP( XALPI - XBETAI/ZTEMP - XGAMI*ALOG(ZTEMP) )
     ZWORK31=(XMV/XMD)*ZWORK31/(XPABSM-ZWORK31)
     ZWORK32=100.*XRM/ZWORK31
  end if

  vapor2rh=ZWORK32

end function vapor2rh
end module conversions
