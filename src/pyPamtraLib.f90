subroutine pyPamtraLib(&!in
set_namelist_file,&
in_ngridx,&
 in_ngridy,&
 max_in_nlyrs,&
 in_nlyrs,&
 in_nfreq,&
 in_freqs,&
in_nfft,&
in_timestamp,&
in_deltax,&
in_deltay,&
 in_lat,&
in_lon,&
in_model_i,&
in_model_j,&
in_wind10u,&
in_wind10v,&
in_lfrac,&
in_relhum_lev,&
in_press_lev,&
in_temp_lev,&
in_hgt_lev,&
in_cwc_q,&
in_iwc_q,&
in_rwc_q,&
in_swc_q,&
in_gwc_q,&
in_hwc_q,&
 in_cwc_n,&
 in_iwc_n,&
 in_rwc_n,&
 in_swc_n,&
 in_gwc_n,&
 in_hwc_n&
,& !meta out
out_gitVersion,&
out_gitHash &
,&!data_out
out_Ze,&
out_Ze_cw,&
out_Ze_rr,&
out_Ze_ci,&
out_Ze_sn,&
out_Ze_gr,&
out_Ze_ha,&
out_Att_hydro,&
out_Att_cw,&
out_Att_rr,&
out_Att_ci,&
out_Att_sn,&
out_Att_gr,&
out_Att_ha,&
out_Att_atmo,&
out_radar_hgt,&
out_tb,&
out_radar_spectra,&
out_radar_snr,&
out_radar_moments,&
out_radar_slope,&
out_radar_quality,&
out_radar_vel,&
out_angles&
)

! python Version of the Pamtra model.

! differences between Pamtra and PyPamtraLib:

! 1) PyPamtra expects relative humidity in Pa/Pa, Pamtra wants relative humidity in % for backwards compatibility
! 2) Only PyPamtra can deal with variing height numbers (Pamtra file format has to be changed for that, otherwise implementation is easy!)
! 3) PyPamtra wants the time as unix timestamp (seconds since 1970) and can handle a different time for each gridpoint

  use kinds
  use constants !physical constants live here
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
  use vars_profile
  use double_moments_module !double moments variables are stored here
  use mod_io_strings !some strings for nice filenames
  use file_mod, only: namelist_file
        use report_module

 !                                                                      
  !     By convention, the quantities followed  by "_lev"
  !     are given at the layer heights while the quantitites w/o
  !     "_lev" are layer average quantities                           
  !                                                                       

  implicit none

!NOTE: f2py cannot deal with kind=dbl or kind=sgl declarations!
!see http://cens.ioc.ee/projects/f2py2e/FAQ.html#q-what-if-fortran-90-code-uses-type-spec-kind-kinds

 !!Settings
!   integer,intent(in) :: set_verbose, set_n_moments, set_isnow_n0, set_liu_type, 
! 
!   real,intent(in) :: set_obs_height     ! upper level output height [m] (> 100000. for satellite)
!   real,intent(in) :: set_emissivity
!   real,intent(in) :: set_N_0snowDsnow, set_N_0grauDgrau, set_N_0rainD, set_N_0hailDhail
!   real,intent(in) :: set_SP, set_as_ratio
!   real,intent(in) :: set_salinity         ! sea surface salinity
!   real,intent(in) :: set_snow_density, set_graupel_density, set_hail_density   
! 
!   logical,intent(in) :: set_dump_to_file   ! flag for profile and ssp dump
!   logical,intent(in) :: set_lphase_flag, &        ! flag for phase function calculation
!        set_lgas_extinction, &    ! gas extinction desired
!        set_lhyd_extinction, &    ! hydrometeor extinction desired
!        set_active, &       ! calculate active stuff
!        set_passive, &         ! calculate passive stuff (with RT3/4)
!        set_use_rain_db, &
!        set_use_snow_db
! 
!   character(5),intent(in) :: set_EM_snow, set_EM_grau, set_EM_ice, set_EM_hail
!   character(1),intent(in) :: set_SD_snow, set_SD_grau, set_SD_rain, set_SD_cloud, set_SD_ice, set_SD_hail
!   character(3),intent(in) :: set_gas_mod, set_rt_mode
!   character(20),intent(in) :: set_moments_file
!   character(100),intent(in) :: set_tmp_path,set_creator, set_data_path
!   character(2),intent(in) :: set_OUTPOL
!   character(1),intent(in) :: set_GROUND_TYPE, set_UNITS

character(300),intent(in) :: set_namelist_file 

!Input


  integer, intent(in) :: in_nfreq, max_in_nlyrs, in_ngridx, in_ngridy, in_nfft
  real, dimension(in_nfreq), intent(in) :: in_freqs


  integer, dimension(in_ngridx,in_ngridy), intent(in) :: in_timestamp
  real, intent(in) :: in_deltax, in_deltay

  integer, dimension(in_ngridx,in_ngridy), intent(in) :: in_nlyrs

  real, dimension(in_ngridx,in_ngridy), intent(in) :: in_lat,in_lon
  integer, dimension(in_ngridx,in_ngridy), intent(in) :: in_model_i,in_model_j
  real, dimension(in_ngridx,in_ngridy), intent(in) :: in_wind10u,in_wind10v,in_lfrac
  real, dimension(in_ngridx,in_ngridy,1+max_in_nlyrs), intent(in) :: in_relhum_lev,in_press_lev,in_temp_lev,in_hgt_lev
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs), intent(in) :: in_cwc_q,in_iwc_q,in_rwc_q,in_swc_q,in_gwc_q
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs), intent(in) :: in_hwc_q,in_cwc_n,in_iwc_n,in_rwc_n
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs), intent(in) :: in_swc_n,in_gwc_n,in_hwc_n
    character(40),intent(out) :: out_gitHash, out_gitVersion

  !Output
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq),intent(out) :: out_Ze,out_Att_hydro,out_Att_atmo
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq),intent(out)::out_Ze_cw,out_Ze_rr,out_Ze_ci,&
                  out_Ze_sn,out_Ze_gr,out_Ze_ha
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq),intent(out)::out_Att_cw,out_Att_rr,out_Att_ci,&
                  out_Att_sn,out_Att_gr,out_Att_ha

  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs),intent(out) :: out_radar_hgt
  real, dimension(32),intent(out) :: out_angles !2*NUMMU instead of 32 does not work, because f2py does not know dimensions!
  real, dimension(in_ngridx,in_ngridy,2,32,in_nfreq,2),intent(out) :: out_tb !same here: noutlevels=2, 2*NUMMU = 32, NSTOKES = 2
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq,in_nfft),intent(out):: out_radar_spectra
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq),intent(out):: out_radar_snr
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq,4),intent(out):: out_radar_moments
  real, dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq,2),intent(out):: out_radar_slope
  integer, dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq),intent(out):: out_radar_quality
  real, dimension(in_nfft),intent(out):: out_radar_vel

  !settings
  !f2py intent(in) :: set_namelist_file
  !input
  !f2py intent(in) :: max_in_nlyrs, in_nlyrs, in_ngridx, in_ngridy,in_nfreq, in_freqs
  !f2py intent(in) :: in_timestamp, in_nfft
  !f2py intent(in) :: in_deltax,in_deltay, in_lat,in_lon,in_model_i,in_model_j
  !f2py intent(in) :: in_wind10u,in_wind10v,in_lfrac
  !f2py intent(in) :: in_relhum_lev,in_press_lev,in_temp_lev,in_hgt_lev
  !f2py intent(in) :: in_cwc_q,in_iwc_q,in_rwc_q,in_swc_q,in_gwc_q
  !f2py intent(in) :: in_hwc_q, in_cwc_n, in_iwc_n, in_rwc_n, in_swc_n, in_gwc_n, in_hwc_n
  !meta out
  !f2py intent(out) :: out_gitVersion,out_gitHash
  !data out
  !f2py intent(out) :: out_Ze,out_Att_hydro,out_Att_atmo,out_radar_hgt,out_tb
  !f2py intent(out) :: out_Ze_cw,out_Ze_rr,out_Ze_ci,out_Ze_sn,out_Ze_gr,out_Ze_ha
  !f2py intent(out) :: out_Att_cw,out_Att_rr,out_Att_ci,out_Att_sn,out_Att_gr,out_Att_ha
  !f2py intent(out) :: out_radar_spectra,out_radar_snr,out_radar_vel,out_angles
  !f2py intent(out) :: out_radar_moments,out_radar_slope,out_radar_quality


  !!!loop variables
  integer ::  fi,nx, ny

  integer,dimension(9) :: timestamp

! 		set_sd_cloud, set_sd_ice, set_em_ice, set_sd_rain, set_n_0raind, set_sd_snow, set_n_0snowdsnow, set_em_snow, set_snow_density, set_sp, set_isnow_n0, set_liu_type, set_sd_grau, set_n_0graudgrau, set_em_grau, set_graupel_density, set_sd_hail, set_n_0haildhail, set_em_hail, set_hail_density,


  interface
    subroutine versionNumber(gitVersion,gitHash)
      implicit none
      character(40), intent(out) ::gitVersion,gitHash
    end subroutine versionNumber

    subroutine allocate_output_vars(no_allocated_lyrs)
      implicit none
      integer, intent(in) :: no_allocated_lyrs
    end  subroutine allocate_output_vars

    subroutine run_rt(nx,ny,fi,freq,frq_str)
      use kinds
      implicit none
      integer, intent(in) :: nx,ny,fi 
      real(kind=dbl), intent(in) :: freq ! frequency [GHz]
      character(8), intent(in) :: frq_str !from commandline
    end subroutine run_rt
  end interface



  !get git data
  call versionNumber(out_gitVersion,out_gitHash)

  namelist_file = set_namelist_file !temporary solution!

  !!! read variables from namelist file
  call nml_params_read !from nml_params.f90

  if ((radar_mode .eq. "spectrum") .and. (radar_nfft .ne. in_nfft)) stop "nfft in python input and nml file must be equal!"

  if (verbose .gt. 1) print*,in_freqs, in_nlyrs, max_in_nlyrs

  in_python = .true.


  ngridx = in_ngridx
  ngridy = in_ngridy
  nfrq = in_nfreq
  freqs = in_freqs(1:nfrq)


  deltax = in_deltax
  deltay = in_deltay



!! read n-moments file
  if (n_moments .eq. 2) call double_moments_module_read(moments_file) !from double_moments_module.f90

  ! now allocate output variables
   call allocate_output_vars(max_in_nlyrs)



   out_Ze = -9999.
   out_Att_hydro = -9999.
   out_Att_atmo = -9999.

   out_Ze_cw = -9999.
   out_Ze_rr = -9999.
   out_Ze_ci = -9999.
   out_Ze_sn = -9999.
   out_Ze_gr = -9999.
   out_Ze_ha = -9999.
   out_Att_cw = -9999.
   out_Att_rr = -9999.
   out_Att_ci = -9999.
   out_Att_sn = -9999.
   out_Att_gr = -9999.
   out_Att_ha = -9999.
   out_radar_spectra = -9999.
   out_radar_snr = -9999.
   out_radar_vel = -9999.
   out_radar_hgt = -9999.
   out_radar_moments = -9999.
   out_radar_slope = -9999.
   out_radar_quality = -9999
   out_angles = -9999.
   out_tb = -9999


!   if (write_nc .eqv. .false.) call mod_io_strings_get_filename()


  if (verbose .gt. 1) print*, 'Start loop over frequencies & profiles!'


  grid_f: do fi =1, nfrq
      if (jacobian_mode) then
      !for jacobian mode. non disturbed profile is expected in grid 1,1!
        call allocate_jacobian_vars
      end if
      grid_y: do ny = 1, ngridy !nx_in, nx_fin   
          grid_x: do nx = 1, ngridx !ny_in, ny_fin  

          call GMTIME(in_timestamp(nx,ny),timestamp)

          write(year,"(i4.4)") timestamp(6)+1900
          write(month,"(i2.2)") timestamp(5)+1
          write(day,"(i2.2)") timestamp(4)
          write(time(1:2),"(i2.2)") timestamp(3)
          write(time(3:4),"(i2.2)") timestamp(2)

          nlyr = in_nlyrs(nx,ny)  
          
          call allocate_profile_vars

          !   ground_temp = profiles(nx,ny)%temp_lev(0)       ! K
          lat = in_lat(nx,ny)                  ! °
          lon = in_lon(nx,ny)                 ! °
          lfrac = in_lfrac(nx,ny)
          model_i = in_model_i (nx,ny)
          model_j = in_model_j(nx,ny)
          wind10u = in_wind10u(nx,ny)
          wind10v = in_wind10v(nx,ny)


          relhum_lev = 100.*in_relhum_lev(nx,ny,1:nlyr+1)         ! %
          press_lev = in_press_lev(nx,ny,1:nlyr+1)           ! Pa
          temp_lev = in_temp_lev(nx,ny,1:nlyr+1)             ! K
          hgt_lev = in_hgt_lev(nx,ny,1:nlyr+1)               ! m

          ! to avoid possible strange effects which are not understood integrated values are set but not used...
          iwv = -1.
          cwp = -1.
          iwp = -1.
          rwp = -1.
          swp = -1.
          gwp = -1.
          hwp = -1.

          cwc_q = in_cwc_q(nx,ny,1:nlyr)           ! kg/kg
          iwc_q = in_iwc_q(nx,ny,1:nlyr)             ! kg/kg
          rwc_q = in_rwc_q(nx,ny,1:nlyr)                  ! kg/kg
          swc_q = in_swc_q(nx,ny,1:nlyr)                  ! kg/kg
          gwc_q = in_gwc_q(nx,ny,1:nlyr)               ! kg/kg

         if (n_moments .eq. 2) then
            hwc_q = in_hwc_q(nx,ny,1:nlyr)              ! kg/kg
            cwc_n = in_cwc_n(nx,ny,1:nlyr)              ! #/kg
            iwc_n = in_iwc_n(nx,ny,1:nlyr)              ! #/kg
            rwc_n = in_rwc_n(nx,ny,1:nlyr)              ! #/kg
            swc_n = in_swc_n(nx,ny,1:nlyr)              ! #/kg
            gwc_n = in_gwc_n(nx,ny,1:nlyr)              ! #/kg
            hwc_n = in_hwc_n(nx,ny,1:nlyr)              ! #/kg
         end if
           !run the model
           call run_rt(nx,ny,fi,freqs(fi),freq_str)
!           if ((active) .and. ((radar_mode .eq. "simple") .or. (radar_mode .eq. "splitted"))) then
!             out_Ze_cw(nx,ny,1:nlyr,:) = REAL(Ze_cw(nx,ny,1:nlyr,:))
!             out_Ze_rr(nx,ny,1:nlyr,:) = REAL(Ze_rr(nx,ny,1:nlyr,:))
!             out_Ze_ci(nx,ny,1:nlyr,:) = REAL(Ze_ci(nx,ny,1:nlyr,:))
!             out_Ze_sn(nx,ny,1:nlyr,:) = REAL(Ze_sn(nx,ny,1:nlyr,:))
!             out_Ze_gr(nx,ny,1:nlyr,:) = REAL(Ze_gr(nx,ny,1:nlyr,:))
!             out_Ze_ha(nx,ny,1:nlyr,:) = REAL(Ze_ha(nx,ny,1:nlyr,:))
!             out_Att_cw(nx,ny,1:nlyr,:) = REAL(Att_cw(nx,ny,1:nlyr,:))
!             out_Att_rr(nx,ny,1:nlyr,:) = REAL(Att_rr(nx,ny,1:nlyr,:))
!             out_Att_ci(nx,ny,1:nlyr,:) = REAL(Att_ci(nx,ny,1:nlyr,:))
!             out_Att_sn(nx,ny,1:nlyr,:) = REAL(Att_sn(nx,ny,1:nlyr,:))
!             out_Att_gr(nx,ny,1:nlyr,:) = REAL(Att_gr(nx,ny,1:nlyr,:))
!             out_Att_ha(nx,ny,1:nlyr,:) = REAL(Att_ha(nx,ny,1:nlyr,:))
!           end if
          call deallocate_profile_vars()
         end do grid_x
    end do grid_y
  if (jacobian_mode) then
  !for jacobian mode
      call deallocate_jacobian_vars
  end if
  end do grid_f

  if (active) then
    out_Ze(:,:,:,:) = REAL(Ze(:,:,:,:))
    out_Att_hydro(:,:,:,:) = REAL(Att_hydro(:,:,:,:))
    out_Att_atmo(:,:,:,:) = REAL(Att_atmo(:,:,:,:))
    out_radar_hgt(:,:,:) = REAL(radar_hgt(:,:,:))
  end if


  if ((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments"))) then
    out_radar_spectra(:,:,:,:,:) = REAL(radar_spectra(:,:,:,:,:))
    out_radar_snr(:,:,:,:) = REAL(radar_snr(:,:,:,:))
    out_radar_vel(:) = REAL(radar_vel(:))
    out_radar_moments(:,:,:,:,:) = REAL(radar_moments(:,:,:,:,:))
    out_radar_slope(:,:,:,:,:) = REAL(radar_slope(:,:,:,:,:))
    out_radar_quality(:,:,:,:) = radar_quality(:,:,:,:)
  end if


  if (passive) then
    out_angles = REAL(angles_deg(:))
    out_tb = RESHAPE( REAL(tb), (/ngridx, ngridy, noutlevels, 2*nummu, nfrq,nstokes /),&
          ORDER = (/6,5,4,3,2,1/))
  end if


  call deallocate_output_vars()

end subroutine pyPamtraLib
