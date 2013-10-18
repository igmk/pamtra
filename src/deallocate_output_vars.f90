subroutine deallocate_output_vars

  use vars_atmosphere
  use vars_output
  use report_module

  implicit none
  !   integer(kind=long), intent(out) :: errorstatus
  !   integer(kind=long) :: err = 0
  character(len=30) :: nameOfRoutine = 'deallocate_output_vars'

  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

  if (allocated(is)) deallocate(is)
  if (allocated(js)) deallocate(js)
  if (allocated(lons)) deallocate(lons)
  if (allocated(lats)) deallocate(lats)
  if (allocated(lfracs)) deallocate(lfracs)
  if (allocated(iwvs)) deallocate(iwvs)
  if (allocated(cwps)) deallocate(cwps)
  if (allocated(iwps)) deallocate(iwps)
  if (allocated(rwps)) deallocate(rwps)
  if (allocated(swps)) deallocate(swps)
  if (allocated(gwps)) deallocate(gwps)
  if (allocated(hwps)) deallocate(hwps)
  if (allocated(tb)) deallocate(tb)
  if (allocated(Ze)) deallocate(Ze)
  if (allocated(radar_hgt)) deallocate(radar_hgt)
  if (allocated(Att_hydro)) deallocate(Att_hydro)
  if (allocated(Att_atmo)) deallocate(Att_atmo)
  if (allocated(Ze_cw)) deallocate(Ze_cw)
  if (allocated(Ze_rr)) deallocate(Ze_rr)
  if (allocated(Ze_ci)) deallocate(Ze_ci)
  if (allocated(Ze_sn)) deallocate(Ze_sn)
  if (allocated(Ze_gr)) deallocate(Ze_gr)
  if (allocated(Ze_ha)) deallocate(Ze_ha)
  if (allocated(Att_cw)) deallocate(Att_cw)
  if (allocated(Att_rr)) deallocate(Att_rr)
  if (allocated(Att_ci)) deallocate(Att_ci)
  if (allocated(Att_sn)) deallocate(Att_sn)
  if (allocated(Att_gr)) deallocate(Att_gr)
  if (allocated(Att_ha)) deallocate(Att_ha)
  if (allocated(radar_spectra)) deallocate(radar_spectra)
  if (allocated(radar_snr)) deallocate(radar_snr)
  if (allocated(radar_vel)) deallocate(radar_vel)
  if (allocated(radar_moments)) deallocate(radar_moments)
  if (allocated(radar_slope)) deallocate(radar_slope)
  if (allocated(radar_slope)) deallocate(radar_edge)
  if (allocated(radar_quality)) deallocate(radar_quality)

  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)


end subroutine deallocate_output_vars
