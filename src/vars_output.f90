module vars_output

! TODO here: clean up modules, change variable names to out_xxx


  use kinds
  use report_module
  implicit none

  save

  integer, allocatable, dimension(:,:) :: is, js
  real(kind=sgl), allocatable, dimension(:,:) :: lons, lats, lfracs, iwvs, cwps, &
       iwps, rwps, swps, gwps,hwps

  !for passive
  real(kind=dbl), allocatable, dimension(:,:,:,:,:,:) :: tb
  real(kind=dbl), dimension(32) :: angles_deg !2*NUMMU=32

  !for active 
  real(kind=dbl), allocatable, dimension(:,:,:) :: radar_hgt
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Ze
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Att_atmo
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Att_hydro
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: radar_spectra
  real(kind=dbl), allocatable, dimension(:,:,:,:) ::    radar_snr
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) ::    radar_moments
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) ::    radar_slopes
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) ::    radar_edge
  integer, allocatable, dimension(:,:,:,:) ::    radar_quality
  real(kind=dbl), allocatable, dimension(:) :: radar_vel

  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: psd_d_bound
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: psd_f
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: psd_mass
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: psd_area

  contains

  subroutine allocate_output_vars(no_allocated_lyrs)

    use vars_atmosphere
    use settings, only: write_nc,&
      active,&
      passive,&
      radar_mode,&
      save_psd,&
      in_python,&
      nstokes,&
      nfrq,&
      nummu,&
      noutlevels, &
      radar_nfft

    use mod_io_strings
    use descriptor_file, only: n_hydro, nbin_arr

    integer, intent(in) :: no_allocated_lyrs
    integer(kind=long)  :: max_nbin1                                           ! maximum number of nbins+1

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


!     if (write_nc .or. in_python) then
    if (passive) then
        allocate(tb(nstokes,nfrq,2*nummu,noutlevels,ngridy,ngridx))
        tb = 0._dbl
    end if

    if (active) then
        allocate(Ze(ngridx,ngridy,no_allocated_lyrs,nfrq))
        allocate(Att_hydro(ngridx,ngridy,no_allocated_lyrs,nfrq))
        allocate(Att_atmo(ngridx,ngridy,no_allocated_lyrs,nfrq))
        allocate(radar_hgt(ngridx,ngridy,no_allocated_lyrs))
        !set to -9999, because height of profiles can vary!
        Ze = -9999._dbl
        Att_hydro= 0._dbl
        Att_atmo= 0._dbl
        radar_hgt= -9999._dbl
    end if

    if((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments")))  then
        allocate(&
        radar_spectra(ngridx,ngridy,no_allocated_lyrs,nfrq,radar_nfft),&
        radar_snr(ngridx,ngridy,no_allocated_lyrs,nfrq),&
        radar_moments(ngridx,ngridy,no_allocated_lyrs,nfrq,4),&
        radar_slopes(ngridx,ngridy,no_allocated_lyrs,nfrq,2),&
        radar_edge(ngridx,ngridy,no_allocated_lyrs,nfrq,2),&
        radar_quality(ngridx,ngridy,no_allocated_lyrs,nfrq),&
        radar_vel(radar_nfft))
        !set to -9999, because height of profiles can vary!
        radar_spectra = -9999.d0
        radar_snr = -9999.d0
        radar_vel = -9999.d0
        radar_moments = -9999.d0
        radar_slopes = -9999.d0
        radar_edge = -9999.d0
        radar_quality = -9999
    end if

    if (save_psd) then
        !how much space do we need?
        max_nbin1 = MAXVAL(nbin_arr) +1
    
        allocate(&
          psd_d_bound(ngridx,ngridy,no_allocated_lyrs,n_hydro,max_nbin1),&
          psd_f(ngridx,ngridy,no_allocated_lyrs,n_hydro,max_nbin1),&
          psd_mass(ngridx,ngridy,no_allocated_lyrs,n_hydro,max_nbin1),&
          psd_area(ngridx,ngridy,no_allocated_lyrs,n_hydro,max_nbin1))

          psd_d_bound = -9999.d0
          psd_f = -9999.d0
          psd_mass = -9999.d0
          psd_area = -9999.d0

    end if

    if (verbose .gt. 1) print*, 'Done allocate_output_vars'

  end subroutine allocate_output_vars

  subroutine deallocate_output_vars

    use vars_atmosphere

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
    if (allocated(radar_spectra)) deallocate(radar_spectra)
    if (allocated(radar_snr)) deallocate(radar_snr)
    if (allocated(radar_vel)) deallocate(radar_vel)
    if (allocated(radar_moments)) deallocate(radar_moments)
    if (allocated(radar_slopes)) deallocate(radar_slopes)
    if (allocated(radar_slopes)) deallocate(radar_edge)
    if (allocated(radar_quality)) deallocate(radar_quality)
    if (allocated(psd_d_bound)) deallocate(psd_d_bound)
    if (allocated(psd_f)) deallocate(psd_f)
    if (allocated(psd_mass)) deallocate(psd_mass)
    if (allocated(psd_area)) deallocate(psd_area)

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)


  end subroutine deallocate_output_vars


end module vars_output

