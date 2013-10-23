! subroutine allocate_output_vars(no_allocated_lyrs)
! 
!     use vars_atmosphere
!     use vars_output
!     use settings, only: write_nc,&
!       active,&
!       passive,&
!       radar_mode,&
!       save_psd,&
!       in_python,&
!       nstokes,&
!       nfrq,&
!       nummu,&
!       noutlevels, &
!       radar_nfft
! 
!     use mod_io_strings
!     use descriptor_file, only: n_hydro, nbin_arr
!     use report_module
!     implicit none
!     integer, intent(in) :: no_allocated_lyrs
!     integer(kind=long)  :: max_nbin1                                           ! maximum number of nbins+1
! 
!     if (verbose .gt. 1) print*, 'Entering allocate_output_vars'
! 
!     if (write_nc) then
!         allocate(is(ngridy,ngridx),js(ngridy,ngridx))
!         allocate(lons(ngridy,ngridx),lats(ngridy,ngridx),lfracs(ngridy,ngridx))
!         allocate(iwvs(ngridy,ngridx))
!         allocate(cwps(ngridy,ngridx),iwps(ngridy,ngridx),rwps(ngridy,ngridx),&
!         swps(ngridy,ngridx),gwps(ngridy,ngridx),hwps(ngridy,ngridx))
!      
!         lons = 0.; lats = 0.; lfracs = 0.;
!         iwvs = 0.; cwps = 0.; iwps = 0.; rwps = 0.; swps = 0.; gwps = 0.; hwps = 0.;
!       
!     end if
! 
! 
! !     if (write_nc .or. in_python) then
!     if (passive) then
!         allocate(tb(nstokes,nfrq,2*nummu,noutlevels,ngridy,ngridx))
!         tb = 0._dbl
!     end if
! 
!     if (active) then
!         allocate(Ze(ngridx,ngridy,no_allocated_lyrs,nfrq))
!         allocate(Att_hydro(ngridx,ngridy,no_allocated_lyrs,nfrq))
!         allocate(Att_atmo(ngridx,ngridy,no_allocated_lyrs,nfrq))
!         allocate(radar_hgt(ngridx,ngridy,no_allocated_lyrs))
!         !set to -9999, because height of profiles can vary!
!         Ze = -9999._dbl
!         Att_hydro= 0._dbl
!         Att_atmo= 0._dbl
!         radar_hgt= -9999._dbl
!     end if
! 
!     if((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments")))  then
!         allocate(&
!         radar_spectra(ngridx,ngridy,no_allocated_lyrs,nfrq,radar_nfft),&
!         radar_snr(ngridx,ngridy,no_allocated_lyrs,nfrq),&
!         radar_moments(ngridx,ngridy,no_allocated_lyrs,nfrq,4),&
!         radar_slopes(ngridx,ngridy,no_allocated_lyrs,nfrq,2),&
!         radar_edge(ngridx,ngridy,no_allocated_lyrs,nfrq,2),&
!         radar_quality(ngridx,ngridy,no_allocated_lyrs,nfrq),&
!         radar_vel(radar_nfft))
!         !set to -9999, because height of profiles can vary!
!         radar_spectra = -9999.d0
!         radar_snr = -9999.d0
!         radar_vel = -9999.d0
!         radar_moments = -9999.d0
!         radar_slopes = -9999.d0
!         radar_edge = -9999.d0
!         radar_quality = -9999
!     end if
! 
!     if (save_psd) then
!         !how much space do we need?
!         max_nbin1 = MAXVAL(nbin_arr) +1
!     
!         allocate(&
!           psd_d_bound(ngridx,ngridy,no_allocated_lyrs,n_hydro,max_nbin1),&
!           psd_f(ngridx,ngridy,no_allocated_lyrs,n_hydro,max_nbin1),&
!           psd_mass(ngridx,ngridy,no_allocated_lyrs,n_hydro,max_nbin1),&
!           psd_area(ngridx,ngridy,no_allocated_lyrs,n_hydro,max_nbin1))
! 
!           psd_d_bound = -9999.d0
!           psd_f = -9999.d0
!           psd_mass = -9999.d0
!           psd_area = -9999.d0
! 
!     end if
! 
! 
!     if (verbose .gt. 1) print*, 'Done allocate_output_vars'
! 
! end subroutine allocate_output_vars
