module sfc_optics

  use kinds, only: dbl
  use report_module
  use settings, only: ground_type, emissivity, nstokes, nummu
  use vars_index, only: i_x, i_y
  use vars_atmosphere, only : atmo_lfrac
  use ocean_sfc_optics, only: ocean_sfc_optics_fastem5
  use land_sfc_optics, only: get_land_sfc_optics
  use vars_rt, only: rt_sfc_emissivity, rt_sfc_reflectivity
  
  implicit none
  
contains
  subroutine set_sfc_optics(errorstatus,freq)
  
    ! Description:
    !   <Say what this routine does>
    !
    ! Method:
    !   <Say how it does it: include references to external documentation>
    !   <If this routine is divided into sections, be brief here,
    !        and put Method comments at the start of each section>
    !
    ! Current Code Owner: <Name of person responsible for this code>
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    ! <version> <date>   Original code. <Your name>
    !
    ! Code Description:
    !   Language:		Fortran 90.
    !   Software Standards: "European Standards for Writing and  
    !     Documenting Exchangeable Fortran 90 Code". 
    !
    ! Declarations:
    ! Modules used:
    
    real(dbl), intent(in) :: freq

    ! Error handling

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'set_sfc_optics'

    if (verbose >= 1) call report(info,'Start of ', nameOfRoutine)

    if ((atmo_lfrac(i_x,i_y) >= 0._dbl) .and. (atmo_lfrac(i_x,i_y) < 0.5_dbl)) then
      call ocean_sfc_optics_fastem5(err,freq)
      ground_type = 'O'
    elseif ((atmo_lfrac(i_x,i_y) >= 0.5_dbl) .and. (atmo_lfrac(i_x,i_y) <= 1.0_dbl)) then
      call get_land_sfc_optics(err,freq)
      ground_type = 'S'
    else
      rt_sfc_emissivity(:,:) = emissivity
      rt_sfc_reflectivity(:,:) = 1._dbl - emissivity
      ground_type = 'L'
    end if

    if (err > 0) then
      errorstatus = fatal
      msg = "assertation error"
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if 

    errorstatus = err
    if (verbose >= 1) call report(info,'End of ', nameOfRoutine)
    
  end subroutine set_sfc_optics

end module sfc_optics

