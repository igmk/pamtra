subroutine deallocate_output_vars

  use vars_atmosphere
  use vars_output
  use nml_params
  use mod_io_strings

  implicit none

  if (write_nc)  deallocate(is,js,lons,lats,lfracs,iwvs,cwps,iwps,rwps,swps,gwps,hwps)
  if (write_nc .or. in_python) deallocate(tb)
  if (active) deallocate(Ze,Attenuation_hydro,Attenuation_atmo,hgt)


end subroutine deallocate_output_vars
