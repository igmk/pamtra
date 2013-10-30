module vars_output

! TODO here: clean up modules, change variable names to out_xxx


  use kinds
  use report_module
  implicit none

  save

!   integer, allocatable, dimension(:,:) :: is, js
!   real(kind=sgl), allocatable, dimension(:,:) :: lons, lats, lfracs, iwvs  
!   !find smarter solution here: 
  real(kind=sgl), allocatable, dimension(:,:) :: cwps, iwps, rwps, swps, gwps,hwps


  !for passive
  real(kind=dbl), allocatable, dimension(:,:,:,:,:,:) :: out_tb
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

  subroutine allocate_output_vars(errorstatus, no_allocated_lyrs)

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
    real(kind=dbl) :: nan

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'allocate_output_vars' 

    call assert_true(err,(atmo_ngridx>0),&
        "atmo_ngridx must be greater zero")   
    call assert_true(err,(atmo_ngridy>0),&
        "atmo_ngridy must be greater zero")   
    call assert_true(err,(no_allocated_lyrs>0),&
        "no_allocated_lyrs must be greater zero")   
    call assert_true(err,(radar_nfft>0),&
        "radar_nfft must be greater zero")   
    call assert_true(err,(nummu>0),&
        "nummu must be greater zero")   
    call assert_true(err,(nstokes>0),&
        "nstokes must be greater zero")   
    call assert_true(err,(nfrq>0),&
        "nfrq must be greater zero")  
    if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
    end if  

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    if (write_nc) then
!         allocate(is(atmo_ngridy,atmo_ngridx),js(atmo_ngridy,atmo_ngridx))
!         allocate(lons(atmo_ngridy,atmo_ngridx),lats(atmo_ngridy,atmo_ngridx),lfracs(atmo_ngridy,atmo_ngridx))
!         allocate(iwvs(atmo_ngridy,atmo_ngridx))
        allocate(cwps(atmo_ngridy,atmo_ngridx),iwps(atmo_ngridy,atmo_ngridx),rwps(atmo_ngridy,atmo_ngridx),&
        swps(atmo_ngridy,atmo_ngridx),gwps(atmo_ngridy,atmo_ngridx),hwps(atmo_ngridy,atmo_ngridx))
      
!         lons = 0.; lats = 0.; lfracs = 0.;
        
        cwps(:,:) = nan()
        iwps(:,:)= nan()
        rwps(:,:)= nan()
        swps(:,:)= nan()
        gwps(:,:)= nan()
        hwps(:,:)= nan()

      
    end if


!     if (write_nc .or. in_python) then
    if (passive) then
        allocate(out_tb(atmo_ngridx,atmo_ngridy,noutlevels,2*nummu,nfrq,nstokes))
        out_tb = 0._dbl
    end if

    if (active) then
        allocate(Ze(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq))
        allocate(Att_hydro(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq))
        allocate(Att_atmo(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq))
        allocate(radar_hgt(atmo_ngridx,atmo_ngridy,no_allocated_lyrs))
        !set to -9999, because height of profiles can vary!
        Ze = -9999._dbl
        Att_hydro= 0._dbl
        Att_atmo= 0._dbl
        radar_hgt= -9999._dbl
    end if

    if((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments")))  then
        allocate(&
        radar_spectra(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq,radar_nfft),&
        radar_snr(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq),&
        radar_moments(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq,4),&
        radar_slopes(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq,2),&
        radar_edge(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq,2),&
        radar_quality(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq),&
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
          psd_d_bound(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,n_hydro,max_nbin1),&
          psd_f(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,n_hydro,max_nbin1),&
          psd_mass(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,n_hydro,max_nbin1),&
          psd_area(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,n_hydro,max_nbin1))

          psd_d_bound = -9999.d0
          psd_f = -9999.d0
          psd_mass = -9999.d0
          psd_area = -9999.d0

    end if

    errorstatus = err
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    return
  end subroutine allocate_output_vars

  subroutine deallocate_output_vars

    use vars_atmosphere

    implicit none
    !   integer(kind=long), intent(out) :: errorstatus
    !   integer(kind=long) :: err = 0
    character(len=30) :: nameOfRoutine = 'deallocate_output_vars'

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

!     if (allocated(is)) deallocate(is)
!     if (allocated(js)) deallocate(js)
!     if (allocated(lons)) deallocate(lons)
!     if (allocated(lats)) deallocate(lats)
!     if (allocated(lfracs)) deallocate(lfracs)
!     if (allocated(iwvs)) deallocate(iwvs)
    if (allocated(cwps)) deallocate(cwps)
    if (allocated(iwps)) deallocate(iwps)
    if (allocated(rwps)) deallocate(rwps)
    if (allocated(swps)) deallocate(swps)
    if (allocated(gwps)) deallocate(gwps)
    if (allocated(hwps)) deallocate(hwps)

    if (allocated(out_tb)) deallocate(out_tb)
    if (allocated(Ze)) deallocate(Ze)
    if (allocated(radar_hgt)) deallocate(radar_hgt)
    if (allocated(Att_hydro)) deallocate(Att_hydro)
    if (allocated(Att_atmo)) deallocate(Att_atmo)
    if (allocated(radar_spectra)) deallocate(radar_spectra)
    if (allocated(radar_snr)) deallocate(radar_snr)
    if (allocated(radar_vel)) deallocate(radar_vel)
    if (allocated(radar_moments)) deallocate(radar_moments)
    if (allocated(radar_slopes)) deallocate(radar_slopes)
    if (allocated(radar_edge)) deallocate(radar_edge)
    if (allocated(radar_quality)) deallocate(radar_quality)
    if (allocated(psd_d_bound)) deallocate(psd_d_bound)
    if (allocated(psd_f)) deallocate(psd_f)
    if (allocated(psd_mass)) deallocate(psd_mass)
    if (allocated(psd_area)) deallocate(psd_area)

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)


  end subroutine deallocate_output_vars


end module vars_output

