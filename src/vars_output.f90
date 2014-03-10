module vars_output

! TODO here: clean up modules, change variable names to out_xxx


  use kinds
  use report_module
  use settings, only: radar_nfft
  implicit none

  save

!   integer, allocatable, dimension(:,:) :: is, js
!   real(kind=sgl), allocatable, dimension(:,:) :: lons, lats, lfracs, iwvs  
!   !find smarter solution here: 
!   real(kind=sgl), allocatable, dimension(:,:) :: cwps, iwps, rwps, swps, gwps,hwps


  !for passive
  real(kind=dbl), allocatable, dimension(:,:,:,:,:,:) :: out_tb
  real(kind=dbl), allocatable, dimension(:) :: out_angles_deg 

  !for active 
  real(kind=dbl), allocatable, dimension(:,:,:) :: out_radar_hgt
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: out_Ze
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: out_att_atmo
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: out_att_hydro
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: out_radar_spectra
  real(kind=dbl), allocatable, dimension(:,:,:,:) ::    out_radar_snr
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) ::    out_radar_moments
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) ::    out_radar_slopes
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) ::    out_radar_edges
  integer, allocatable, dimension(:,:,:,:) ::    out_radar_quality
  real(kind=dbl), allocatable, dimension(:) :: out_radar_vel

  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: out_psd_d
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: out_psd_n
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: out_psd_mass
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: out_psd_area

  real(kind=dbl), dimension(300) :: out_debug_diameter
  real(kind=dbl), dimension(300) :: out_debug_back_of_d
  real(kind=dbl), allocatable, dimension(:) :: out_debug_radarvel
  real(kind=dbl), allocatable, dimension(:) :: out_debug_radarback
  real(kind=dbl), allocatable, dimension(:) :: out_debug_radarback_wturb
  real(kind=dbl), allocatable, dimension(:) :: out_debug_radarback_wturb_wnoise

  
  
  contains

  subroutine allocate_output_vars(errorstatus, no_allocated_lyrs)

    use vars_atmosphere
    use settings, only: write_nc,&
      active,&
      passive,&
      radar_mode,&
      save_psd,&
!       in_python,&
      nstokes,&
      nfrq,&
      nummu,&
      noutlevels, &
      radar_nfft

    use mod_io_strings
    use descriptor_file, only: n_hydro, nbin_arr

    integer, intent(in) :: no_allocated_lyrs
    integer(kind=long)  :: max_nbin                                           ! maximum number of nbins+1
    real(kind=dbl) :: nan

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'allocate_output_vars' 

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

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

    if (verbose >= 5) print*, atmo_ngridx, atmo_ngridy, no_allocated_lyrs, &
      radar_nfft, nummu, nstokes, nfrq, MAXVAL(nbin_arr) +1

    if (write_nc) then
!         allocate(is(atmo_ngridy,atmo_ngridx),js(atmo_ngridy,atmo_ngridx))
!         allocate(lons(atmo_ngridy,atmo_ngridx),lats(atmo_ngridy,atmo_ngridx),lfracs(atmo_ngridy,atmo_ngridx))
!         allocate(iwvs(atmo_ngridy,atmo_ngridx))
!         allocate(cwps(atmo_ngridy,atmo_ngridx),iwps(atmo_ngridy,atmo_ngridx),rwps(atmo_ngridy,atmo_ngridx),&
!         swps(atmo_ngridy,atmo_ngridx),gwps(atmo_ngridy,atmo_ngridx),hwps(atmo_ngridy,atmo_ngridx))
      
!         lons = 0.; lats = 0.; lfracs = 0.;
        
!         cwps(:,:) = nan()
!         iwps(:,:)= nan()
!         rwps(:,:)= nan()
!         swps(:,:)= nan()
!         gwps(:,:)= nan()
!         hwps(:,:)= nan()

      
    end if

    allocate(out_angles_deg(2*nummu))
    out_angles_deg = nan()
!     if (write_nc .or. in_python) then
    if (passive) then
        allocate(out_tb(atmo_ngridx,atmo_ngridy,noutlevels,2*nummu,nfrq,nstokes))
        out_tb = 0._dbl
    end if

    if (active) then
        allocate(out_Ze(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq))
        allocate(out_att_hydro(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq))
        allocate(out_att_atmo(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq))
        allocate(out_radar_hgt(atmo_ngridx,atmo_ngridy,no_allocated_lyrs))
        !set to -9999, because height of profiles can vary!
        out_Ze = -9999._dbl
        out_att_hydro= 0._dbl
        out_att_atmo= 0._dbl
        out_radar_hgt= -9999._dbl
    end if

    if((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments")))  then
        allocate(&
        out_radar_spectra(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq,radar_nfft),&
        out_radar_snr(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq),&
        out_radar_moments(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq,4),&
        out_radar_slopes(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq,2),&
        out_radar_edges(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq,2),&
        out_radar_quality(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,nfrq),&
        out_radar_vel(radar_nfft))
        !set to -9999, because height of profiles can vary!
        out_radar_spectra = -9999.d0
        out_radar_snr = -9999.d0
        out_radar_vel = -9999.d0
        out_radar_moments = -9999.d0
        out_radar_slopes = -9999.d0
        out_radar_edges = -9999.d0
        out_radar_quality = -9999
    end if

    if (save_psd) then
        !how much space do we need?
        max_nbin = MAXVAL(nbin_arr) 
    
        allocate(&
          out_psd_d(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,n_hydro,max_nbin),&
          out_psd_n(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,n_hydro,max_nbin),&
          out_psd_mass(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,n_hydro,max_nbin),&
          out_psd_area(atmo_ngridx,atmo_ngridy,no_allocated_lyrs,n_hydro,max_nbin))

          out_psd_d = -9999.d0
          out_psd_n = -9999.d0
          out_psd_mass = -9999.d0
          out_psd_area = -9999.d0

    end if

    !debuging stuff, only deallocated if necessary:
    if (allocated(out_debug_radarvel)) deallocate(out_debug_radarvel)
    if (allocated(out_debug_radarback)) deallocate(out_debug_radarback)
    if (allocated(out_debug_radarback_wturb)) deallocate(out_debug_radarback_wturb)
    if (allocated(out_debug_radarback_wturb_wnoise)) deallocate(out_debug_radarback_wturb_wnoise)
    
    allocate(out_debug_radarvel(radar_nfft))
    allocate(out_debug_radarback(radar_nfft))
    allocate(out_debug_radarback_wturb(radar_nfft))
    allocate(out_debug_radarback_wturb_wnoise(radar_nfft))
    out_debug_diameter(:) = 0.d0
    out_debug_back_of_d(:) = 0.d0
    
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
!     if (allocated(cwps)) deallocate(cwps)
!     if (allocated(iwps)) deallocate(iwps)
!     if (allocated(rwps)) deallocate(rwps)
!     if (allocated(swps)) deallocate(swps)
!     if (allocated(gwps)) deallocate(gwps)
!     if (allocated(hwps)) deallocate(hwps)

    if (allocated(out_angles_deg)) deallocate(out_angles_deg)
    if (allocated(out_tb)) deallocate(out_tb)
    if (allocated(out_Ze)) deallocate(out_Ze)
    if (allocated(out_radar_hgt)) deallocate(out_radar_hgt)
    if (allocated(out_att_hydro)) deallocate(out_att_hydro)
    if (allocated(out_att_atmo)) deallocate(out_att_atmo)
    if (allocated(out_radar_spectra)) deallocate(out_radar_spectra)
    if (allocated(out_radar_snr)) deallocate(out_radar_snr)
    if (allocated(out_radar_vel)) deallocate(out_radar_vel)
    if (allocated(out_radar_moments)) deallocate(out_radar_moments)
    if (allocated(out_radar_slopes)) deallocate(out_radar_slopes)
    if (allocated(out_radar_edges)) deallocate(out_radar_edges)
    if (allocated(out_radar_quality)) deallocate(out_radar_quality)
    if (allocated(out_psd_d)) deallocate(out_psd_d)
    if (allocated(out_psd_n)) deallocate(out_psd_n)
    if (allocated(out_psd_mass)) deallocate(out_psd_mass)
    if (allocated(out_psd_area)) deallocate(out_psd_area)

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)


  end subroutine deallocate_output_vars


end module vars_output

