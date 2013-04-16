!
Module report_module

    use kinds, only: long
    use settings, only: verbose

    implicit none
    save
    ! status values
    Integer(Kind=long), Parameter :: nstatus = 3
    Integer(Kind=long), Parameter :: success = 0
    Integer(Kind=long), Parameter :: warning = 1
    Integer(Kind=long), Parameter :: fatal   = 2
    Integer(Kind=long), Parameter :: info    = 3
    Character(len=*), Parameter :: status_text(0:nstatus) = &
    (/ 'success', &
    'warning', &
    'fatal  ', &
    'info   '  /)

contains
    subroutine report(status, message, nameOfRoutine)
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
        !
        ! Subroutine arguments
        !   Scalar arguments with intent(in):
        Integer(Kind=long) , Intent (in) :: status     !
        Character (len=*) , Intent (in) :: message    ! ..to output
        Character (len=*) , Intent (in) :: nameOfRoutine ! ..calling this one

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


        If ( status > 0 .And. status < nstatus ) Then
            if (verbose >= 0) then
                ! Display message for warning and fatal
                Write(Error_Unit,"(1X,a4,'/',a2,'/',a2,2x,2(a2,':'),a2,2x,a,a,a)") &
                date(1:4), date(5:6), date(7:8), &
                time(1:2), time(3:4), time(5:6), &
                status_text(status),&
                " in module ",Trim(nameOfRoutine)
                Write(Error_Unit,"(5X,A)") Trim(message)
            end if
        Else
            ! This is the verbose and success message
            Write(Error_Unit,"(1X,a4,'/',a2,'/',a2,2x,2(a2,':'),a2,2x,a,a,x,a)") &
            date(1:4), date(5:6), date(7:8), &
            time(1:2), time(3:4), time(5:6), &
            status_text(status),&
            Trim(message),Trim(NameOfRoutine)
        Endif

    End Subroutine report
end module report_module
