subroutine deallocate_jacobian_vars
  use vars_atmosphere
  use report_module

  implicit none
  !   integer(kind=long), intent(out) :: errorstatus
  !   integer(kind=long) :: err = 0
  character(len=30) :: nameOfRoutine = 'deallocate_jacobian_vars'

  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

  if (allocated(jac_scattermatrix)) deallocate(jac_scattermatrix)
  if (allocated(jac_extmatrix)) deallocate(jac_extmatrix)
  if (allocated(jac_emisvec)) deallocate(jac_emisvec)
  if (allocated(jac_hydros_present)) deallocate(jac_hydros_present)
  if (allocated(jac_kextsn)) deallocate(jac_kextsn)
  if (allocated(jac_backsn)) deallocate(jac_backsn)
  if (allocated(jac_kextcw)) deallocate(jac_kextcw)
  if (allocated(jac_backcw)) deallocate(jac_backcw)
  if (allocated(jac_kextrr)) deallocate(jac_kextrr)
  if (allocated(jac_backrr)) deallocate(jac_backrr)
  if (allocated(jac_kextgr)) deallocate(jac_kextgr)
  if (allocated(jac_backgr)) deallocate(jac_backgr)
  if (allocated(jac_kextci)) deallocate(jac_kextci)
  if (allocated(jac_backci)) deallocate(jac_backci)
  if (allocated(jac_kextha)) deallocate(jac_kextha)
  if (allocated(jac_backha)) deallocate(jac_backha)
  if (allocated(jac_relhum_lev)) deallocate(jac_relhum_lev)
  if (allocated(jac_temp_lev)) deallocate(jac_temp_lev)
  if (allocated(jac_cwc_q)) deallocate(jac_cwc_q)
  if (allocated(jac_iwc_q)) deallocate(jac_iwc_q)
  if (allocated(jac_rwc_q)) deallocate(jac_rwc_q)
  if (allocated(jac_swc_q)) deallocate(jac_swc_q)
  if (allocated(jac_gwc_q)) deallocate(jac_gwc_q)
  if (allocated(jac_hwc_q)) deallocate(jac_hwc_q)
  if (allocated(jac_cwc_n)) deallocate(jac_cwc_n)
  if (allocated(jac_iwc_n)) deallocate(jac_iwc_n)
  if (allocated(jac_rwc_n)) deallocate(jac_rwc_n)
  if (allocated(jac_swc_n)) deallocate(jac_swc_n)
  if (allocated(jac_gwc_n)) deallocate(jac_gwc_n)
  if (allocated(jac_hwc_n)) deallocate(jac_hwc_n)

  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

end subroutine deallocate_jacobian_vars
