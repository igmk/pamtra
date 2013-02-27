subroutine get_gasabs &
    (freq,errstat) ! in

! Description:
!  Based on on the namelist parameter gas_mod, the routine for calculating the
!  gas absorption due to air (N2 and O2) and water vapor for the input frequency
!  is selected.
!
!  The absorption coefficient is converted from Np/km to Np/m in the end.
!
! Input parameters:
!    frequency                    GHz
!    temperature                   K
!    water vapor density         kg/m**3
!    total air pressure            Pa
! Output parameters:
!    absorption coefficient       Np/m
!
! Owner: IGMK
!
! History:
!
! Version Date         Comment
! ------- ----         -------
! 0.01    21/09/2008   Code adaption from A. Battaglia - M. Mech
! 0.1     21/09/2009   Include Rosenkranz gas absorption model - M. Mech
! 0.2     26/02/2013   Application of European Standards for Writing and
!                      Documenting Exchangeable Fortran 90 Code - M. Mech
!
! Code Description:
!  Language: Fortran 90.
!  Software Standards: "European Standards for Writing and
!   Documenting Exchangeable Fortran 90 Code".
!
! Parent Module: run_rt
!
! Declarations:

! Modules used:

  use kinds, only: dbl, & ! integer parameter specifying double precision
                    long   ! integer parameter specifying long integer
  use vars_atmosphere, only: nlyr, press, temp, vapor_pressure, rho_vap, kextatmo
  use constants, only: t_abs, &
                        errorstatus_fatal
  use nml_params, only: gas_mod, verbose

  implicit none

#include "error_report.interface"
!- End of header ---------------------------------------------------------------

! Subroutine arguments
! Scalar arguments with intent(in):
  real(kind=dbl), intent(in) :: freq
! End of Subroutine arguments

! Local scalars:
  integer(kind=long) :: nz

  real(kind=dbl) :: absair, abswv
  real(kind=dbl) :: tc

! Error handling

  integer(kind=long), intent(out) :: ErrStat
  character(len=80) :: ErrMsg
  character(len=14) :: NameOfRoutine = 'get_gasabs'

  do nz = 1, nlyr          
     tc = temp(nz) - t_abs
     if (gas_mod .eq. 'L93') then
        call mpm93(freq, press(nz)/1.d3, vapor_pressure(nz)/1.d3,tc, 0.d0, kextatmo(nz),errstat)
        if (errstat == errorstatus_fatal) Then
           errmsg = 'error in mpm93'
           call error_report(errstat, errmsg, NameOfRoutine)
           return
        end if
        kextatmo(nz) = kextatmo(nz)/1.d3
     else if (gas_mod .eq. 'R98') then
        call rosen98_gasabs(freq,temp(nz),rho_vap(nz),press(nz),absair,abswv,errstat)
        if (errstat == errorstatus_fatal) Then
           errmsg = 'error in rosen98_gasabs'
           call error_report(errstat, errmsg, NameOfRoutine)
           return
        end if
        kextatmo(nz) = (absair + abswv)/1.d3    ! conversion to Np/m
     else
        ! here is the question whether we want the possibility to switch off the gas absorption
        kextatmo(nz) = 0
        errmsg = 'No gas absorption model specified!'
        errstat = errorstatus_fatal
        call error_report(errstat, errmsg, NameOfRoutine)
        return
     end if
  end do
  if (verbose > 1) print*, 'variables filled up!'
  return 

end subroutine get_gasabs
