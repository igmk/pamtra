! subroutine deallocate_output_vars
! 
!   use vars_atmosphere
!   use vars_output
!   use report_module
! 
!   implicit none
!   !   integer(kind=long), intent(out) :: errorstatus
!   !   integer(kind=long) :: err = 0
!   character(len=30) :: nameOfRoutine = 'deallocate_output_vars'
! ! 
!   if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
! 
!   if (allocated(is)) deallocate(is)
!   if (allocated(js)) deallocate(js)
!   if (allocated(lons)) deallocate(lons)
!   if (allocated(lats)) deallocate(lats)
!   if (allocated(lfracs)) deallocate(lfracs)
!   if (allocated(iwvs)) deallocate(iwvs)
!   if (allocated(cwps)) deallocate(cwps)
!   if (allocated(iwps)) deallocate(iwps)
!   if (allocated(rwps)) deallocate(rwps)
!   if (allocated(swps)) deallocate(swps)
!   if (allocated(gwps)) deallocate(gwps)
!   if (allocated(hwps)) deallocate(hwps)
!   if (allocated(tb)) deallocate(tb)
!   if (allocated(Ze)) deallocate(Ze)
!   if (allocated(radar_hgt)) deallocate(radar_hgt)
!   if (allocated(Att_hydro)) deallocate(Att_hydro)
!   if (allocated(Att_atmo)) deallocate(Att_atmo)
!   if (allocated(radar_spectra)) deallocate(radar_spectra)
!   if (allocated(radar_snr)) deallocate(radar_snr)
!   if (allocated(radar_vel)) deallocate(radar_vel)
!   if (allocated(radar_moments)) deallocate(radar_moments)
!   if (allocated(radar_slopes)) deallocate(radar_slopes)
!   if (allocated(radar_slopes)) deallocate(radar_edge)
!   if (allocated(radar_quality)) deallocate(radar_quality)
!   if (allocated(psd_d_bound)) deallocate(psd_d_bound)
!   if (allocated(psd_f)) deallocate(psd_f)
!   if (allocated(psd_mass)) deallocate(psd_mass)
!   if (allocated(psd_area)) deallocate(psd_area)
! 
!   if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
! 
! 
! end subroutine deallocate_output_vars
