subroutine deallocate_jacobian_vars
  use kinds
  use vars_atmosphere
  use nml_params

  implicit none



  deallocate(jac_scattermatrix,jac_extmatrix,&
      jac_emisvec)
  deallocate(jac_kextsn, jac_backsn, jac_kextcw, &
      jac_backcw, jac_kextrr, jac_backrr, jac_kextgr,&
      jac_backgr, jac_kextci, jac_backci, jac_kextha,&
      jac_backha)
  deallocate(jac_relhum_lev,jac_temp_lev,&
      jac_cwc_q, jac_iwc_q, jac_rwc_q, jac_swc_q,&
      jac_gwc_q)
  if (n_moments .eq. 2) then
      deallocate(jac_hwc_q, jac_cwc_n, jac_iwc_n,&
	jac_rwc_n, jac_swc_n, jac_gwc_n, jac_hwc_n)
  end if
end subroutine deallocate_jacobian_vars