module sfc_optics

  ! variables
  use kinds
  use report_module
  use settings, only: emissivity!, nstokes, nummu
  use vars_index, only: i_x, i_y, i_f
  use vars_atmosphere, only : sfc_type, sfc_model
  use vars_rt, only: rt_sfc_emissivity, rt_sfc_reflectivity
  use vars_output, only: out_emissivity
  ! routines
  use ocean_sfc_optics
  use land_sfc_optics
  
  implicit none
  
contains
  subroutine set_sfc_optics(errorstatus,freq)
    ! Description:
    !   Within this routine the surface and the reflection type are determined.
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
    ! 0.9     27/10/16   Complete redesign
    !
    ! Code Description:
    !   Language:            Fortran 90.
    !   Software Standards: "European Standards for Writing and  
    !     Documenting Exchangeable Fortran 90 Code". 
    !
    ! Declarations:
    ! Modules used: ocean_sfc_optics, land_sfc_optics
    
    real(dbl), intent(in) :: freq

    ! Error handling

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'set_sfc_optics'

    if (verbose >= 1) call report(info,'Start of ', nameOfRoutine)
    
    if (sfc_type(i_x,i_y) == 0) then ! water
        if ((sfc_model(i_x,i_y) == 0) .or. (sfc_model(i_x,i_y) == -9999)) then ! TESSEM2
            call ocean_sfc_optics_tessem2(err,freq)
        elseif (sfc_model(i_x,i_y) == 1) then ! FASTEM
            call ocean_sfc_optics_fastemx(err,freq)
        end if
    elseif (sfc_type(i_x,i_y) == 1) then ! land
        if ((sfc_model(i_x,i_y) == 0) .or. (sfc_model(i_x,i_y) == -9999)) then ! TELSEM2
            call land_sfc_optics_telsem2(err,freq)
        elseif (sfc_model(i_x,i_y) == 1) then ! SSMI
            call land_sfc_optics_ssmi(err,freq)
        end if
    else ! default sfc_type == -9999, sfc_model == -9999 and sfc_refl == 'L'
        rt_sfc_emissivity(:,:) = emissivity
        rt_sfc_reflectivity(:,:) = 1._dbl - emissivity
    end if

    if (err /= 0) then
      errorstatus = fatal
      msg = "error in land_sfc_optics_xxx or ocean_sfc_optics_xxx"
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if
    
    out_emissivity(i_x,i_y,:,i_f,:) = rt_sfc_emissivity(:,:)

    errorstatus = err
   
    if (verbose >= 1) call report(info,'End of ', nameOfRoutine)
    
  end subroutine set_sfc_optics

end module sfc_optics

