!
Subroutine error_report(ErrStatus, ErrMessage, NameOfRoutine)
  ! Description:
  !   Write out fatal and warning error messages to unit 6.
  !   Execution is stopped in the event of a fatal error.
  !
  ! Copyright:
  !
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  !
  !
  ! Current Code Owner: IGMK
  !
  ! History:
  !
  ! Version  Date       Comment
  ! -------  ----       -------
  ! 0.1      13/02/2013 Adaption of the idea from rttov_ErrorReport within RTTOV - M. Mech
  !
  ! Code Description:
  !  Language: Fortran 90.
  !  Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations
  !
  ! Global variables:
  ! Modules used:
  !
  Use kinds, Only : long

  Use constants, Only : &
      errorstatus_fatal,&
      ErrorStatus_text,&
      NerrorStatus
  use nml_params, only: verbose
!  Use rttov_global, Only : &
!      verbose_message,&
!      err_init,&
!      error_unit
!INTF_ON


  Implicit None
  !
  ! Subroutine arguments
  !   Scalar arguments with intent(in):
  Integer(Kind=long) , Intent (in) :: ErrStatus     ! +ve => fatal error, -ve => warning
  Character (len=*) , Intent (in) :: ErrMessage    ! ..to output
  Character (len=*) , Intent (in) :: NameOfRoutine ! ..calling this one

  ! local
  integer(kind=long), parameter :: error_unit = 6
  Character (len=8) :: date
  Character (len=10):: time
  !- End of header --------------------------------------------------------

  Call DATE_AND_Time(date, time)

!  ! If globlal variables not defined then use default values
!  if( .not. err_init ) then
!     call rttov_errorhandling ( -1_jpim , -1_jpim)
!  endif


  If ( ErrStatus >= 0 .And. ErrStatus <=nerrorstatus ) Then
     ! Display message only if allowed by verbose flag
     if( verbose > 0 .or.  errstatus == errorstatus_fatal) then
        Write(Error_Unit,"(1X,a4,'/',a2,'/',a2,2x,2(a2,':'),a2,2x,a,a,a)") &
                & date(1:4), date(5:6), date(7:8), &
                & time(1:2), time(3:4), time(5:6), &
                & ErrorStatus_text(errstatus),&
                & " in module ",Trim(NameOfRoutine)
        Write(Error_Unit,"(5X,A)") Trim(ErrMessage)
     Endif
  Else
     ! This error level is different from predefined
     ! Output it anyway
     Write(Error_Unit,"(1X,a4,'/',a2,'/',a2,2x,2(a2,':'),a2,2x,i6,a,a)") &
             & date(1:4), date(5:6), date(7:8), &
             & time(1:2), time(3:4), time(5:6), &
             & errstatus                      ,&
             & " in module ",Trim(NameOfRoutine)
     Write(Error_Unit,"(5X,A)") Trim(ErrMessage)
  Endif

End Subroutine error_report
