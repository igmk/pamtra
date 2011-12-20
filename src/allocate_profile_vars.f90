subroutine allocate_profile_vars

  use vars_atmosphere
  use vars_output
  use nml_params
  use mod_io_strings

  implicit none


  allocate(hgt_lev(0:nlyr),stat=alloc_status)
  allocate(press_lev(0:nlyr),stat=alloc_status)
  allocate(press(nlyr),stat=alloc_status)
  allocate(temp_lev(0:nlyr),stat=alloc_status)
  allocate(temp(nlyr),stat=alloc_status)
  allocate(relhum_lev(0:nlyr),stat=alloc_status)
  allocate(relhum(nlyr),stat=alloc_status)

  allocate(vapor_pressure(nlyr),stat=alloc_status)
  allocate(rho_vap(nlyr),stat=alloc_status)
  allocate(q_hum(nlyr),stat=alloc_status)

  allocate(cwc_q(nlyr))
  allocate(iwc_q(nlyr))
  allocate(rwc_q(nlyr))
  allocate(swc_q(nlyr))
  allocate(gwc_q(nlyr))
  if (n_moments .eq. 2) then
     allocate(hwc_q(nlyr))
     allocate(cwc_n(nlyr))
     allocate(iwc_n(nlyr))
     allocate(rwc_n(nlyr))
     allocate(swc_n(nlyr))
     allocate(gwc_n(nlyr))
     allocate(hwc_n(nlyr))
  end if

  allocate(nlegen(nlyr),stat=alloc_status)
  allocate(rt3nlegen(nlyr),stat=alloc_status)
  
  allocate(kextatmo(nlyr), stat=alloc_status)
  allocate(kexttot(nlyr), stat=alloc_status)
  allocate(kextsn(nlyr), stat=alloc_status)
  allocate(kextcw(nlyr), stat=alloc_status)
  allocate(kextrr(nlyr), stat=alloc_status)
  allocate(kextgr(nlyr), stat=alloc_status)
  allocate(kextci(nlyr), stat=alloc_status)
  allocate(kextha(nlyr), stat=alloc_status)
  
  allocate(salbtot(nlyr), stat=alloc_status)
  allocate(rt3kexttot(nlyr), stat=alloc_status)
  allocate(rt3salbtot(nlyr), stat=alloc_status)
  allocate(g_coeff(nlyr), stat=alloc_status)
  
  allocate(back(nlyr), stat=alloc_status)
  allocate(backcw(nlyr), stat=alloc_status)
  allocate(backrr(nlyr), stat=alloc_status)
  allocate(backci(nlyr), stat=alloc_status)
  allocate(backsn(nlyr), stat=alloc_status)
  allocate(backgr(nlyr), stat=alloc_status)
  allocate(backha(nlyr), stat=alloc_status)

  allocate(legen(nlyr,200), stat=alloc_status)
  allocate(legen2(nlyr,200), stat=alloc_status)
  allocate(legen3(nlyr,200), stat=alloc_status)
  allocate(legen4(nlyr,200), stat=alloc_status)
  allocate(rt3legen(nlyr,200), stat=alloc_status)
  allocate(rt3legen2(nlyr,200), stat=alloc_status)
  allocate(rt3legen3(nlyr,200), stat=alloc_status)
  allocate(rt3legen4(nlyr,200), stat=alloc_status)

!   allocate(ics(ngridx, ngridy))
  if (dump_to_file) then
	allocate(file_ph(nlyr))
end if

!   allocate(angles_deg(2*NUMMU))

  ! set them to zero, just in case they are not calculated but used for Ze/PIA calculation
  kexttot(:) = 0d0
  kextatmo(:) = 0d0
  back(:) = 0d0


end subroutine allocate_profile_vars
