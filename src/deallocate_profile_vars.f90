subroutine deallocate_profile_vars

  use vars_atmosphere
  use vars_output
  use settings
  use mod_io_strings
  use report_module
  
  implicit none

!   integer(kind=long), intent(out) :: errorstatus
!   integer(kind=long) :: err = 0
  character(len=30) :: nameOfRoutine = 'deallocate_profile_vars'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  if (allocated(hgt_lev)) deallocate(hgt_lev)
  if (allocated(press_lev)) deallocate(press_lev)
  if (allocated(press)) deallocate(press)
  if (allocated(temp_lev)) deallocate(temp_lev)
  if (allocated(temp)) deallocate(temp)
  if (allocated(relhum_lev)) deallocate(relhum_lev)
  if (allocated(relhum)) deallocate(relhum)
  if (allocated(vapor_pressure)) deallocate(vapor_pressure)
  if (allocated(rho_vap)) deallocate(rho_vap)
  if (allocated(q_hum)) deallocate(q_hum)
  if (allocated(cwc_q)) deallocate(cwc_q)
  if (allocated(iwc_q)) deallocate(iwc_q)
  if (allocated(rwc_q)) deallocate(rwc_q)
  if (allocated(swc_q)) deallocate(swc_q)
  if (allocated(gwc_q)) deallocate(gwc_q)
  if (allocated(hgt)) deallocate(hgt)
  if (allocated(delta_hgt_lev)) deallocate(delta_hgt_lev)
  if (allocated(hwc_q)) deallocate(hwc_q)
  if (allocated(cwc_n)) deallocate(cwc_n)
  if (allocated(iwc_n)) deallocate(iwc_n)
  if (allocated(rwc_n)) deallocate(rwc_n)
  if (allocated(swc_n)) deallocate(swc_n)
  if (allocated(gwc_n)) deallocate(gwc_n)
  if (allocated(hwc_n)) deallocate(hwc_n)
  if (allocated(q_hydro)) deallocate(q_hydro)
  if (allocated(nlegen)) deallocate(nlegen)
  if (allocated(kextatmo)) deallocate(kextatmo)
  if (allocated(kexttot)) deallocate(kexttot)
  if (allocated(salbtot)) deallocate(salbtot)
  if (allocated(g_coeff)) deallocate(g_coeff)
  if (allocated(back)) deallocate(back)
  if (allocated(legen)) deallocate(legen)
  if (allocated(legen2)) deallocate(legen2)
  if (allocated(legen3)) deallocate(legen3)
  if (allocated(legen4)) deallocate(legen4)
  if (allocated(rt4hydros_present)) deallocate(rt4hydros_present)
  if (allocated(hydros_present)) deallocate(hydros_present)
  if (allocated(file_ph)) deallocate(file_ph)

  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

end subroutine deallocate_profile_vars
