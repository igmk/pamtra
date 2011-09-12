subroutine allocate_output_vars(no_allocated_lyrs)


  use vars_atmosphere
  use vars_output
  use nml_params
  use mod_io_strings
  implicit none
  integer, intent(in) :: no_allocated_lyrs



  if (write_nc) then
     allocate(is(ngridy,ngridx),js(ngridy,ngridx))
     allocate(lons(ngridy,ngridx),lats(ngridy,ngridx),lfracs(ngridy,ngridx))
     allocate(iwvs(ngridy,ngridx))
     allocate(cwps(ngridy,ngridx),iwps(ngridy,ngridx),rwps(ngridy,ngridx),&
          swps(ngridy,ngridx),gwps(ngridy,ngridx),hwps(ngridy,ngridx))
     
      lons = 0.; lats = 0.; lfracs = 0.;
      iwvs = 0.; cwps = 0.; iwps = 0.; rwps = 0.; swps = 0.; gwps = 0.; hwps = 0.;
      
  end if


  if (write_nc .or. in_python) then
     allocate(tb(nstokes,nfrq,2*nummu,noutlevels,ngridy,ngridx))
     tb = 0.
  end if

  if (active) then
     allocate(Ze(ngridx,ngridy,no_allocated_lyrs,nfrq))
     allocate(Attenuation_hydro(ngridx,ngridy,no_allocated_lyrs,nfrq))
     allocate(Attenuation_atmo(ngridx,ngridy,no_allocated_lyrs,nfrq))
     allocate(hgt(ngridx,ngridy,no_allocated_lyrs))
  end if


end subroutine allocate_output_vars
