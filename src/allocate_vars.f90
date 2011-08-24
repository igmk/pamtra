subroutine allocate_vars

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
  allocate(salbtot(nlyr), stat=alloc_status)
  allocate(rt3kexttot(nlyr), stat=alloc_status)
  allocate(rt3salbtot(nlyr), stat=alloc_status)
  allocate(g_coeff(nlyr), stat=alloc_status)
  allocate(back(nlyr), stat=alloc_status)
  allocate(legen(nlyr,200), stat=alloc_status)
  allocate(legen2(nlyr,200), stat=alloc_status)
  allocate(legen3(nlyr,200), stat=alloc_status)
  allocate(legen4(nlyr,200), stat=alloc_status)
  allocate(rt3legen(nlyr,200), stat=alloc_status)
  allocate(rt3legen2(nlyr,200), stat=alloc_status)
  allocate(rt3legen3(nlyr,200), stat=alloc_status)
  allocate(rt3legen4(nlyr,200), stat=alloc_status)

  allocate(ics(ngridx, ngridy))
  allocate(file_ph(nlyr))

  if (write_nc) then
     allocate(is(ngridy,ngridx),js(ngridy,ngridx))
     allocate(lons(ngridy,ngridx),lats(ngridy,ngridx),lfracs(ngridy,ngridx))
     allocate(iwvs(ngridy,ngridx))
     allocate(cwps(ngridy,ngridx),iwps(ngridy,ngridx),rwps(ngridy,ngridx),&
          swps(ngridy,ngridx),gwps(ngridy,ngridx),hwps(ngridy,ngridx))
     allocate(tb(nstokes,nfrq,2*nummu,noutlevels,ngridy,ngridx))
     lons = 0.; lats = 0.; lfracs = 0.;
     iwvs = 0.; cwps = 0.; iwps = 0.; rwps = 0.; swps = 0.; gwps = 0.; hwps = 0.;
     tb = 0.

  end if

  if (active) then
     allocate(Ze(ngridx,ngridy,nlyr,nfrq))
     allocate(Attenuation_hydro(ngridx,ngridy,nlyr,nfrq))
     allocate(Attenuation_atmo(ngridx,ngridy,nlyr,nfrq))
     allocate(hgt(ngridx,ngridy,nlyr))
  end if

!   allocate(angles_deg(2*NUMMU))

  ! set them to zero, just in case they are not calculated but used for Ze/PIA calculation
  kexttot(:) = 0d0
  kextatmo(:) = 0d0
  back(:) = 0d0


end subroutine allocate_vars
