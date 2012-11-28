subroutine allocate_output_vars(no_allocated_lyrs)


  use vars_atmosphere
  use vars_output
  use nml_params
  use mod_io_strings
  implicit none
  integer, intent(in) :: no_allocated_lyrs

  if (verbose .gt. 1) print*, 'Entering allocate_output_vars'

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

  if ((active) .and. ((radar_mode .eq. "simple") .or. (radar_mode .eq. "splitted"))) then
     allocate(Ze(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Ze_cw(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Ze_rr(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Ze_ci(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Ze_sn(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Ze_gr(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Ze_ha(ngridx,ngridy,no_allocated_lyrs,nfrq))
     allocate(Att_hydro(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Att_cw(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Att_rr(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Att_ci(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Att_sn(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Att_gr(ngridx,ngridy,no_allocated_lyrs,nfrq),&
             Att_ha(ngridx,ngridy,no_allocated_lyrs,nfrq))
     allocate(Att_atmo(ngridx,ngridy,no_allocated_lyrs,nfrq))
     allocate(radar_hgt(ngridx,ngridy,no_allocated_lyrs))
  end if
  if((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments")))  then
    allocate(&
	  radar_spectra(ngridx,ngridy,no_allocated_lyrs,nfrq,radar_nfft),&
	  radar_snr(ngridx,ngridy,no_allocated_lyrs,nfrq),&
	  radar_moments(ngridx,ngridy,no_allocated_lyrs,nfrq,4),&
	  radar_slope(ngridx,ngridy,no_allocated_lyrs,nfrq,2),&
	  radar_quality(ngridx,ngridy,no_allocated_lyrs,nfrq),&
	  radar_vel(radar_nfft))
	  radar_spectra = -9999.d0
	  radar_snr = -9999.d0
	  radar_vel = -9999.d0
	  radar_moments = -9999.d0
	  radar_slope = -9999.d0
	  radar_quality = -9999
  end if
  if (verbose .gt. 1) print*, 'Done allocate_output_vars'

end subroutine allocate_output_vars
