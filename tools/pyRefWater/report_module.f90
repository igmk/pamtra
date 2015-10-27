!
Module report_module

    use kinds
    use vars_index, only: i_x, i_y, i_z, i_f, i_h
    implicit none
    save
    ! status values
    Integer(kind=long) :: verbose=0
    Integer(kind=long), Parameter :: nstatus = 3
    Integer(kind=long), Parameter :: success = 0
    Integer(kind=long), Parameter :: warning = 1
    Integer(kind=long), Parameter :: fatal   = 2
    Integer(kind=long), Parameter :: info    = 3
    character(len=7), Parameter :: status_text(0:nstatus) = &
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


!to do: add i_x etc to error messages

        Integer(Kind=long) , Intent (in) :: status     !
        Character (len=*) , Intent (in) :: message    ! ..to output
        Character (len=*) , Intent (in) :: nameOfRoutine ! ..calling this one

        ! local
        integer(kind=long), parameter :: error_unit = 6
        Character (len=8) :: date
        Character (len=10):: time
        !- End of header --------------------------------------------------------

        Call DATE_AND_Time(date, time)

        If ( (status > 0) .and. (status < nstatus) ) Then
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
    
    Subroutine assert_true(error,logic,message)
       ! Description:
        ! Subroutine to simplify unit checks  
        !        
        ! Current Code Owner: IGMK
        !
        ! History:
        !
        ! Version  Date       Comment
        ! -------  ----       -------
        ! 0.1      11/07/2013 Initial Implemetation M.Maahn
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
        implicit none
        ! Subroutine arguments
        !   Scalar arguments with intent(inout):
        integer, intent(inout) :: error
        !   Logic arguments with intent(in):
	logical, intent(in) :: logic
        !   String arguments with intent(in):
        Character (len=*) , Intent (in) :: message    ! ..to output
  
	if (.not. logic) then
	  call report(fatal, message, "report_module: assert_true")
	  error = fatal
        else if (verbose >= 5) then
          call report(info, "PASSED: "//message, "report_module: assert_true")
	end if
	return
    End Subroutine assert_true
    
    Subroutine assert_false(error,logic,message)
       ! Description see assert_true
        implicit none
        integer, intent(inout) :: error
        !   Logic arguments with intent(in):
	logical, intent(in) :: logic
        !   String arguments with intent(in):
        Character (len=*) , Intent (in) :: message    ! ..to output
	if (logic) then
	  call report(fatal, message, "report_module: assert_false")
	  error = fatal
        else
          if (verbose >= 5) call report(info, "PASSED: "//message, "report_module: assert_false")
	end if
       

	return
    End Subroutine assert_false

end module report_module
