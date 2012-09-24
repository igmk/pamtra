subroutine allocate_jacobian_vars
  use kinds
  use vars_atmosphere
  use nml_params

  implicit none


  allocate(jac_scattermatrix(nlyr,nstokes,nummu,nstokes,nummu,4),stat=alloc_status)
  allocate(jac_extmatrix(nlyr,nstokes,nstokes,nummu,4),stat=alloc_status)
  allocate(jac_emisvec(nlyr,nstokes,nummu,4),stat=alloc_status)
  allocate(jac_kextsn(nlyr),stat=alloc_status)
  allocate(jac_backsn(nlyr),stat=alloc_status)
  allocate(jac_kextcw(nlyr),stat=alloc_status)
  allocate(jac_backcw(nlyr),stat=alloc_status)
  allocate(jac_kextrr(nlyr),stat=alloc_status)
  allocate(jac_backrr(nlyr),stat=alloc_status)
  allocate(jac_kextgr(nlyr),stat=alloc_status)
  allocate(jac_backgr(nlyr),stat=alloc_status)
  allocate(jac_kextci(nlyr),stat=alloc_status)
  allocate(jac_backci(nlyr),stat=alloc_status)
  allocate(jac_kextha(nlyr),stat=alloc_status)
  allocate(jac_backha(nlyr),stat=alloc_status)
  allocate(jac_salbtot(nlyr),stat=alloc_status)
  allocate(jac_hydros_present(nlyr),stat=alloc_status)
  
  allocate(jac_temp_lev(0:nlyr),stat=alloc_status)
  allocate(jac_relhum_lev(0:nlyr),stat=alloc_status)
  allocate(jac_cwc_q(nlyr),stat=alloc_status)
  allocate(jac_iwc_q(nlyr),stat=alloc_status)
  allocate(jac_rwc_q(nlyr),stat=alloc_status)
  allocate(jac_swc_q(nlyr),stat=alloc_status)
  allocate(jac_gwc_q(nlyr),stat=alloc_status)
  if (n_moments .eq. 2) then
    allocate(jac_hwc_q(nlyr),stat=alloc_status)
    allocate(jac_cwc_n(nlyr),stat=alloc_status)
    allocate(jac_iwc_n(nlyr),stat=alloc_status)
    allocate(jac_rwc_n(nlyr),stat=alloc_status)
    allocate(jac_swc_n(nlyr),stat=alloc_status)
    allocate(jac_gwc_n(nlyr),stat=alloc_status)
    allocate(jac_hwc_n(nlyr),stat=alloc_status)
  end if


end subroutine allocate_jacobian_vars