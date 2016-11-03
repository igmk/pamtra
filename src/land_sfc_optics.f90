module land_sfc_optics
use kinds
use report_module
use settings, only: data_path, verbose
use mod_io_strings, only: nxstr, nystr
use vars_index, only: i_x, i_y
use vars_atmosphere, only:atmo_month,atmo_lat, atmo_lon
use vars_rt, only: rt_sfc_emissivity, rt_sfc_reflectivity
  
implicit none

contains
subroutine land_sfc_optics_telsem2(errorstatus, freq)

  real(dbl),intent(in) :: freq
  
  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg    
  character(14) :: nameOfRoutine = 'land_sfc_optics_telsem2'

  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
  err= 0
  
  errorstatus = err

  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

end subroutine land_sfc_optics_telsem2

subroutine land_sfc_optics_ssmi(errorstatus, freq)
  
  integer(long) :: imonth
  real(dbl) :: land_emissivity, freq
  character(80) :: femis ! filename for the emissivity databases
  
  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg    
  character(14) :: nameOfRoutine = 'land_sfc_optics_ssmi'

  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
  err= 0

  read(atmo_month(i_x,i_y),'(i2)') imonth
  if ((imonth .ge. 7) .and. (imonth .le. 12)) then
      femis = data_path(:len_trim(data_path))//'/emissivity/ssmi_mean_emis_92'//atmo_month(i_x,i_y)//'_direct'
  else if ((imonth .ge. 1 ).and. (imonth .lt. 7)) then
      femis = data_path(:len_trim(data_path))//'/emissivity/ssmi_mean_emis_93'//atmo_month(i_x,i_y)//'_direct'
  else
      msg = "Warning: No emissivity data found for "//nxstr//" and "//nystr
      errorstatus = fatal
      call report(errorstatus,msg,nameOfRoutine)
      return
  end if

  ! land_emis could give polarized reflectivities
  call land_emis_ssmi(err,femis,atmo_lon(i_x,i_y),atmo_lat(i_x,i_y),freq,land_emissivity)

  if (land_emissivity < 0.01_dbl) then
      land_emissivity = 0.94_dbl
  end if
  rt_sfc_emissivity(:,:) = land_emissivity
  rt_sfc_reflectivity(:,:) = 1._dbl - land_emissivity

  errorstatus = err

  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
  
end subroutine land_sfc_optics_ssmi
end module land_sfc_optics