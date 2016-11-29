module land_sfc_optics

use kinds
use report_module
use settings, only: nummu, nstokes
use vars_index, only: i_x, i_y
use vars_atmosphere, only: atmo_month, atmo_lat, atmo_lon, atmo_model_i, atmo_model_j
use vars_rt, only: rt_sfc_emissivity, rt_sfc_reflectivity
  
implicit none

contains
    subroutine land_sfc_optics_telsem2(errorstatus, freq)

        real(dbl), intent(in) :: freq
        integer(kind=long), intent(out) :: errorstatus
        integer(kind=long) :: err = 0
        character(len=80) :: msg    
        character(23) :: nameOfRoutine = 'land_sfc_optics_telsem2'

        real(dbl), dimension(2,nummu) :: land_emissivity

        if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

        call telsem2(err,atmo_month(i_x,i_y),atmo_lon(i_x,i_y),atmo_lat(i_x,i_y),freq,land_emissivity)
!        print*, land_emissivity(:,0)

        if (err /= 0) then
            print*, i_x, i_y,atmo_lon(i_x,i_y),atmo_lat(i_x,i_y),&
            atmo_model_i(i_x,i_y), atmo_model_j(i_x,i_y)
            msg = 'error in land surface emissivity'
            call report(err, msg, nameOfRoutine)
            errorstatus = err
            return
        end if

        rt_sfc_emissivity(:,:) = land_emissivity(:,:)
        rt_sfc_reflectivity(:,:) = 1._dbl - land_emissivity(:,:)

        errorstatus = err

        if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

    end subroutine land_sfc_optics_telsem2

    subroutine land_sfc_optics_ssmi(errorstatus, freq)
    
        real(dbl), intent(in) :: freq
            
        integer(kind=long), intent(out) :: errorstatus
        integer(kind=long) :: err = 0
        character(len=80) :: msg    
        character(19) :: nameOfRoutine = 'land_sfc_optics_ssmi'

        real(dbl) :: land_emissivity

        if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

        call land_emis_ssmi(err,atmo_month(i_x,i_y),atmo_lon(i_x,i_y),atmo_lat(i_x,i_y),freq,land_emissivity)
        print*, land_emissivity

        if (land_emissivity < 0.01_dbl) then
            land_emissivity = 0.94_dbl
        end if
        rt_sfc_emissivity(:,:) = land_emissivity
        rt_sfc_reflectivity(:,:) = 1._dbl - land_emissivity

        if (err /= 0) then
            print*, i_x, i_y
            msg = 'error in land surface emissivity'
            call report(err, msg, nameOfRoutine)
            errorstatus = err
            return
        end if
        errorstatus = err

        if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    
    end subroutine land_sfc_optics_ssmi
end module land_sfc_optics