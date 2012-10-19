subroutine deallocate_output_vars

  use vars_atmosphere
  use vars_output
  use nml_params
  use mod_io_strings

  implicit none

  if (write_nc)  deallocate(is,js,lons,lats,lfracs,iwvs,cwps,iwps,rwps,swps,gwps,hwps)
  if (write_nc .or. in_python) deallocate(tb)
  if (active) deallocate(Ze,Ze_cw,Ze_rr,Ze_ci,Ze_sn,Ze_gr,Ze_ha, &
               Att_hydro,Att_atmo,Att_cw,Att_rr,Att_ci,Att_sn,Att_gr,Att_ha, &
               hgt)
  if(radar_spectrum) deallocate(radar_spectra, radar_snr, radar_vel)

end subroutine deallocate_output_vars
