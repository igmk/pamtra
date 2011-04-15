subroutine small_scale(mu, freq, wind10, xcorr1)
!**** *SMALL_SCALE* -  MODIFIES FRESNEL
!       REFLECTION COEFFICIENTS FOR SMALL SCALE RIPPLES.
! PURPOSE.
! --------
!       TO CALCULATE A CORRECTION TO THE FRESNEL REFLECTION COEFFICIENTS
!       ALLOWING FOR THE PRESENCE OF SMALL RIPPLES (Bragg Scattering).
!
!** INTERFACE.
!   ----------
!	mu	     (Input)     cosine of theta
!       freq         (Input)     frequency in GHz
!       wind10       (INPUT) R4  10M NEUTRAL WIND SPEED (M/S)
!       XCORR1     (OUTPUT)  R4  REFLECTIVITY CORRECTION   
!
! METHOD.
! -------
!       SEE REFERENCES.
! EXTERNALS.
! ----------
!       NONE
! REFERENCE.
! ----------
!       [1] ENGLISH S.J., 1997: A FAST EMISSIVITY MODEL FOR
!              ATOVS. NWP TECHNICAL REPORT (BEING WRITTEN).
!
!       [2] WU, S. T., AND FUNG A.K., 1972:  A NONCOHERENT MODEL
!	       FOR MICROWAVE EMISSION AND BACKSCATTERING FROM THE SEA
!              SURFACE. J.G.R., 77, 5917-5929.
!
!       [3] GUILLOU C., 1996: ETUDE DU TRANSFERT RADIATIF DANS
!              L'ATMOSPHERE EN MICRO-ONDES, A PARTIR D'OBSERVATIONS
!              RADIOMETRIQUES AEROPORTEES. THESE, UNIVERSITE DE PARIS.
! AUTHOR.
! -------
!       S.J.ENGLISH        *UKMO*     6/11/97
! MODIFICATIONS.
! --------------
!
!        

use data, only: emc
use kinds
  implicit none

  real(kind=dbl) :: mu, freq, freq_sq
  real(kind=dbl) :: mu_sq, wind10
  real(kind=dbl) :: xcorr1

  mu_sq = mu*mu
  freq_sq = freq*freq
  if (freq.gt.0.1) then
    xcorr1=exp(emc(21)*wind10*mu_sq/freq_sq)
  else
    xcorr1=1.0
  endif

  return

end subroutine small_scale

subroutine large_scale(mu, freq, wind10, xcorr2)
!**** *LARGE_SCALE_CORRECTION* -  MODIFIES FRESNEL
!       REFLECTION COEFFICIENTS FOR GEOMETRIC REFLECTIONS
! PURPOSE.
! --------
!       ROUTINE: TO CALCULATE A CORRECTION TO THE FRESNEL REFLECTION COEFFICIENTS
!       ALLOWING FOR THE PRESENCE OF LARGE SCALE ROUGHNESS      
!
!** INTERFACE.
!   ----------
!       *CALL* *LARGE_SCALE_CORRECTION.F
!       U10MPS       (INPUT) R4  10M NEUTRAL WIND SPEED (M/S)
!       XCORR2     (OUTPUT)  R4  REFLECTIVITY CORRECTION   
!
!       ALSO TAKES INPUT FROM *COMMON/(VARIABLES)/ *
!       AND *COMMON/CONSTANTS/ *.
!       AND PUTS OUTPUT IN *COMMON/(OUTPUT/ *.
! METHOD.
! -------
!       SEE REFERENCES.
! EXTERNALS.
! ----------
!         NONE
! REFERENCE.
! ----------
!        [1] ENGLISH S.J., 1997: A FAST EMISSIVITY MODEL FOR
!              ATOVS. NWP TECHNICAL REPORT (BEING WRITTEN).
! AUTHOR.
! -------
!       S.J.ENGLISH        *UKMO*     3/11/97
! MODIFICATIONS.
! --------------
!
!           SEE *COMMON* DECK.
!#include "cparam.h"
!#include "emismw.h"

  use data, only: emc
  use kinds

  implicit none

  integer :: jp, jc, iemc
  real(kind=dbl), dimension(6) :: zc
  real(kind=dbl), dimension(2) :: xcorr2
  real(kind=dbl) :: freq, wind10, mu, usec
  real(kind=dbl) :: wind10_sq, freq_sq
  real(kind=dbl) :: sec, sec_sq
! EFFECTIVE SPECULAR EMISSIVITY CALCULATED
!
! CALCULATE FREQUENCY, WINDSPEED AND VIEW ANGLE SQUARED FOR
! POLYNOMIALS
  freq_sq = freq*freq
  wind10_sq = wind10*wind10

  sec=1.0/mu
  sec_sq=sec*sec
  usec=wind10*sec
! JP=1 => V POL, JP=2 => H POL 
! JP: TWO POLARISTIONS
  do jp=1,2
! JC: SIX COEFFICIENTS (CONSTANT, U, U^2, SEC, SEC^2, U*SEC)	
    do jc=1,6
! SELECT ELEMENT FROM DATA ARRAY EMC ELEMENTS 24-59 FOR THIS MODEL
      iemc=24+(jc-1)*3+(jp-1)*18
! COEFFICIENTS 1-5 FOR THIS POLARISATION STORED IN ZC
      zc(jc)=emc(iemc)+emc(iemc+1)*freq   &
           +emc(iemc+2)*freq_sq
    end do
! CALCULATE CORRECTION FOR THIS POLARISATION
    xcorr2(jp)=zc(1)+zc(2)*sec+zc(3)*sec_sq  &
             +zc(4)*wind10+zc(5)*wind10_sq+zc(6)*usec
    xcorr2(jp)=xcorr2(jp)/100.0
  end do

  return

end subroutine large_scale

subroutine foam(wind10, ffoam)
!**** *FOAM* -  MODEL OF FOAM COVERAGE FOR *RT* CALCULATION
!
! PURPOSE.
! --------
!       ROUTINE
!       TO CALCULATE FOAM COVERAGE AT A GIVEN NEUTRAL
!       STABILITY WINDSPEED. IN PRACTISE THE 10M WINDSPEED
!       IS USUALLY USED WITHOUT A STABILITY CORRECTION.
!
!** INTERFACE.
!   ----------
!       *CALL* *FOAM_CORRECTION
!       U10MPS       (INPUT)   R4  10M NEUTRAL WIND SPEED (M/S)
!       FFOAM        (OUTPUT)  R4  FOAM COVERAGE FRACTIONC
!       ALSO TAKES INPUT FROM *COMMON/(VARIABLES)/ *
!       AND *COMMON/CONSTANTS/ *.
!       AND PUTS OUTPUT IN *COMMON/(OUTPUT/ *.
! METHOD.
! -------
!       SEE REFERENCES.
! EXTERNALS.
! ----------
!       NONE.
! REFERENCE.
! ----------
!       NWP TECH MEMO #?
!
!       [1] MONAHAN E.C. AND O'MUIRCHEARTAIGH I.G., 1986: WHITECAPS
!             AND THE PASSIVE REMOTE SENSING OF THE OCEAN'S SURFACE., 
!             INT. J. OF REMOTE SENSING, 7 (NO 5), 627-642.
!
!       [2] WU S.T., 1979: OCEANIC WHITECAPS AND SEA STATE., J. PHYS.
!              OCEANOGRAPHY, 9, 1064-1078.
!
! AUTHOR.
! -------
!       S.J.ENGLISH        *UKMO*     6/11/97
! MODIFICATIONS.
! --------------
!
! 
!*     *COMMON*C           SEE *COMMON* DECK.
!#include "cparam.h"
!#include "emismw.h"
!        
  use data, only: emc
  use kinds

  implicit none

  real(kind=dbl) :: ffoam, wind10

  ffoam=emc(22)*(wind10**emc(23))

  return

end subroutine foam
