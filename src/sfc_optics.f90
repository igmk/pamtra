module sfc_optics

  ! variables
  use kinds, only: dbl
  use report_module
  use settings, only: ground_type, emissivity!, nstokes, nummu
  use vars_index, only: i_x, i_y
  use vars_atmosphere, only : atmo_lfrac
  use vars_rt, only: rt_sfc_emissivity, rt_sfc_reflectivity
  ! routines
  use ocean_sfc_optics, only: ocean_sfc_optics_fastemx
  use land_sfc_optics, only: get_land_sfc_optics
  
  implicit none
  
contains
  subroutine set_sfc_optics(errorstatus,freq)
    ! Description:
    !   Within this routine the type of surface reflection is determined.
    !
    ! Method:
    !   This detemination is done based on the land-sea fraction parameter.
    !   Future developments should do this based on a surface type to include
    !   different surface, e.g., ice, forest, sand.
    !
    ! Current Code Owner: Mario Mech
    !
    ! History:
    ! Version   Date     Comment
    ! -------   ----     -------
    ! 0.1     24/02/15   Original code. Mario Mech
    !
    ! Code Description:
    !   Language:            Fortran 90.
    !   Software Standards: "European Standards for Writing and  
    !     Documenting Exchangeable Fortran 90 Code". 
    !
    ! Declarations:
    ! Modules used: ocean_sfc_optics_fastem5, get_land_sfc_optics
    
    real(dbl), intent(in) :: freq

    ! Error handling

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'set_sfc_optics'

    if (verbose >= 1) call report(info,'Start of ', nameOfRoutine)
    err = 0
    if ((atmo_lfrac(i_x,i_y) >= 0._dbl) .and. (atmo_lfrac(i_x,i_y) < 0.5_dbl)) then
      call ocean_sfc_optics_fastemx(err,freq)
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
      msg = "error in get_land_sfc_optics or ocean_sfc_optics_fastemx"
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if 

    errorstatus = err
    if (verbose >= 1) call report(info,'End of ', nameOfRoutine)
    
  end subroutine set_sfc_optics

end module sfc_optics

