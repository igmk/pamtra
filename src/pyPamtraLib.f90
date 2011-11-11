subroutine pyPamtraLib(&!settings
set_verbose,set_dump_to_file,set_tmp_path,&
set_data_path,set_obs_height,set_units,set_outpol,set_creator,&
set_active,set_passive,set_ground_type,set_salinity,set_emissivity,set_lgas_extinction,&
set_gas_mod,set_lhyd_extinction,set_lphase_flag,set_SD_snow,set_N_0snowDsnow,set_EM_snow,&
set_SP,set_isnow_n0,set_liu_type,set_SD_grau,set_N_0grauDgrau,set_EM_grau,set_EM_ice,set_SD_rain,&
set_N_0rainD,set_n_moments,set_moments_file&
,& !in
in_ngridx, in_ngridy, max_in_nlyrs, in_nlyrs, in_nfreq, in_freqs,&
in_timestamp,&
in_deltax,in_deltay, in_lat,in_lon,in_model_i,in_model_j,&
in_wind10u,in_wind10v,in_lfrac,&
in_relhum_lev,in_press_lev,in_temp_lev,in_hgt_lev,&
in_iwv,in_cwp,in_iwp,in_rwp,in_swp,in_gwp,in_hwp,&
in_cwc_q,in_iwc_q,in_rwc_q,in_swc_q,in_gwc_q,&
in_hwc_q, in_cwc_n, in_iwc_n, in_rwc_n, in_swc_n, in_gwc_n, in_hwc_n&
,& !meta out
out_gitVersion,out_gitHash &
,& !data_out
out_Ze,&
out_Ze_cw,out_Ze_rr,out_Ze_ci,out_Ze_sn,out_Ze_gr,out_Ze_ha,&
out_Att_hydro,&
out_Att_cw,out_Att_rr,out_Att_ci,out_Att_sn,out_Att_gr,out_Att_ha,&
out_Att_atmo,out_hgt,out_tb,&
out_angles&
)




! python Version of the Pamtra model.

! differences between Pamtra and PyPamtraLib:

! 1) PyPamtra expects relative humidity in Pa/Pa, Pamtra wants relative humidity in % for backwards compatibility
! 2) Only PyPamtra can deal with variing height numbers (Pamtra file format has to be changed for that, otherwise implementation is easy!)
! 3) PyPamtra wants the time as unix timestamp (seconds since 1070) and can handle different times per gridpoint

  use kinds
  use constants !physical constants live here
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
  use vars_profile
  use double_moments_module !double moments variables are stored here
  use mod_io_strings !some strings for nice filenames

 !                                                                      
  !     By convention, the quantities followed  by "_lev"
  !     are given at the layer heights while the quantitites w/o
  !     "_lev" are layer average quantities                           
  !                                                                       

  implicit none


 !!Settings
  integer,intent(in) :: set_verbose, set_n_moments, set_isnow_n0, set_liu_type

  real(kind=sgl),intent(in) :: set_obs_height     ! upper level output height [m] (> 100000. for satellite)
  real(kind=sgl),intent(in) :: set_emissivity
  real(kind=sgl),intent(in) :: set_N_0snowDsnow, set_N_0grauDgrau, set_N_0rainD, set_SP
  real(kind=sgl),intent(in) :: set_salinity         ! sea surface salinity

  logical,intent(in) :: set_dump_to_file   ! flag for profile and ssp dump
  logical,intent(in) :: set_lphase_flag, &        ! flag for phase function calculation
       set_lgas_extinction, &    ! gas extinction desired
       set_lhyd_extinction, &    ! hydrometeor extinction desired
       set_active, &       ! calculate active stuff
       set_passive         ! calculate passive stuff (with RT3)

  character(5),intent(in) :: set_EM_snow, set_EM_grau, set_EM_ice
  character(1),intent(in) :: set_SD_snow, set_SD_grau, set_SD_rain
  character(3),intent(in) :: set_gas_mod
  character(20),intent(in) :: set_moments_file
  character(100),intent(in) :: set_tmp_path,set_creator, set_data_path
  character(2),intent(in) :: set_OUTPOL
  character(1),intent(in) :: set_GROUND_TYPE, set_UNITS
!Input

  integer, intent(in) :: in_nfreq, max_in_nlyrs, in_ngridx, in_ngridy
  real(kind=sgl), dimension(in_nfreq), intent(in) :: in_freqs


  integer, dimension(in_ngridx,in_ngridy), intent(in) :: in_timestamp
  real(kind=sgl), intent(in) :: in_deltax, in_deltay

  integer, dimension(in_ngridx,in_ngridy), intent(in) :: in_nlyrs

  real(kind=sgl), dimension(in_ngridx,in_ngridy), intent(in) :: in_lat,in_lon
  integer, dimension(in_ngridx,in_ngridy), intent(in) :: in_model_i,in_model_j
  real(kind=sgl), dimension(in_ngridx,in_ngridy), intent(in) :: in_wind10u,in_wind10v,in_lfrac
  real(kind=sgl), dimension(in_ngridx,in_ngridy,1+max_in_nlyrs), intent(in) :: in_relhum_lev,in_press_lev,in_temp_lev,in_hgt_lev
  real(kind=sgl), dimension(in_ngridx,in_ngridy), intent(in) :: in_iwv,in_cwp,in_iwp,in_rwp,in_swp,in_gwp,in_hwp
  real(kind=sgl), dimension(in_ngridx,in_ngridy,max_in_nlyrs), intent(in) :: in_cwc_q,in_iwc_q,in_rwc_q,in_swc_q,in_gwc_q
  real(kind=sgl), dimension(in_ngridx,in_ngridy,max_in_nlyrs), intent(in) :: in_hwc_q,in_cwc_n,in_iwc_n,in_rwc_n
  real(kind=sgl), dimension(in_ngridx,in_ngridy,max_in_nlyrs), intent(in) :: in_swc_n,in_gwc_n,in_hwc_n
    character(40),intent(out) :: out_gitHash, out_gitVersion

  !Output
  real(kind=sgl), dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq),intent(out) :: out_Ze,out_Att_hydro,out_Att_atmo
  real(kind=sgl), dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq),intent(out)::out_Ze_cw,out_Ze_rr,out_Ze_ci,&
                  out_Ze_sn,out_Ze_gr,out_Ze_ha
  real(kind=sgl), dimension(in_ngridx,in_ngridy,max_in_nlyrs,in_nfreq),intent(out)::out_Att_cw,out_Att_rr,out_Att_ci,&
                  out_Att_sn,out_Att_gr,out_Att_ha
  real(kind=sgl), dimension(in_ngridx,in_ngridy,max_in_nlyrs),intent(out) :: out_hgt
  real(kind=sgl), dimension(32),intent(out) :: out_angles !2*NUMMU instead of 32 does not work, because f2py does not know dimensions!
  real(kind=sgl), dimension(in_ngridx,in_ngridy,2,32,in_nfreq,2),intent(out) :: out_tb !same here: noutlevels=2, 2*NUMMU = 32, NSTOKES = 2

  !settings
  !f2py intent(in) :: set_verbose,set_dump_to_file,set_tmp_path
  !f2py intent(in) :: set_data_path,set_obs_height,set_units,set_outpol,set_creator
  !f2py intent(in) :: set_active,set_passive,set_ground_type,set_salinity,set_emissivity,set_lgas_extinction
  !f2py intent(in) :: set_gas_mod,set_lhyd_extinction,set_lphase_flag,set_SD_snow,set_N_0snowDsnow,set_EM_snow
  !f2py intent(in) :: set_SP,set_isnow_n0,set_liu_type,set_SD_grau,set_N_0grauDgrau,set_EM_grau,set_EM_ice,set_SD_rain
  !f2py intent(in) :: set_N_0rainD,set_n_moments,set_moments_file
  !input
  !f2py intent(in) :: max_in_nlyrs, in_nlyrs, in_ngridx, in_ngridy,in_nfreq, in_freqs
  !f2py intent(in) :: in_timestamp
  !f2py intent(in) :: in_deltax,in_deltay, in_lat,in_lon,in_model_i,in_model_j
  !f2py intent(in) :: in_wind10u,in_wind10v,in_lfrac
  !f2py intent(in) :: in_relhum_lev,in_press_lev,in_temp_lev,in_hgt_lev
  !f2py intent(in) :: in_iwv,in_cwp,in_iwp,in_rwp,in_swp,in_gwp,in_hwp
  !f2py intent(in) :: in_cwc_q,in_iwc_q,in_rwc_q,in_swc_q,in_gwc_q
  !f2py intent(in) :: in_hwc_q, in_cwc_n, in_iwc_n, in_rwc_n, in_swc_n, in_gwc_n, in_hwc_n
  !meta out
  !f2py intent(out) :: out_gitVersion,out_gitHash
  !data out
  !f2py intent(out) :: out_Ze,out_Att_hydro,out_Att_atmo,out_hgt,out_tb
  !f2py intent(out) :: out_Ze_cw,out_Ze_rr,out_Ze_ci,out_Ze_sn,out_Ze_gr,out_Ze_ha
  !f2py intent(out) :: out_Att_cw,out_Att_rr,out_Att_ci,out_Att_sn,out_Att_gr,out_Att_ha
  !f2py intent(out) :: out_angles



  !!!loop variables
  integer ::  fi,nx, ny

  integer,dimension(9) :: timestamp




if (set_verbose .gt. 1) print*,in_freqs, in_nlyrs, max_in_nlyrs


  !get git data
  call versionNumber(out_gitVersion,out_gitHash)


  in_python = .true.
  !write_nc must be true to collect the results
  write_nc = .false.

  !these are not(?) needed any more
  input_path = ""
  output_path = "/tmp"
  freq_str = "pythonFrequen"
  file_desc = ""



  !load settings, uggly but neccessary!
  verbose = set_verbose
  dump_to_file = set_dump_to_file
  tmp_path = set_tmp_path
  data_path = set_data_path
  obs_height = set_obs_height
  units = set_units
  outpol = set_outpol
  creator = set_creator
  active = set_active
  passive = set_passive
  ground_type = set_ground_type
  salinity = set_salinity
  emissivity = set_emissivity
  lgas_extinction = set_lgas_extinction
  gas_mod = set_gas_mod
  lhyd_extinction = set_lhyd_extinction
  lphase_flag = set_lphase_flag
  SD_snow = set_SD_snow
  N_0snowDsnow = set_N_0snowDsnow
  EM_snow = set_EM_snow
  SP = set_SP
  isnow_n0 = set_isnow_n0
  liu_type = set_liu_type
  SD_grau = set_SD_grau
  N_0grauDgrau = set_N_0grauDgrau
  EM_grau = set_EM_grau
  EM_ice = set_EM_ice
  SD_rain = set_SD_rain
  N_0rainD = set_N_0rainD
  n_moments = set_n_moments
  moments_file = set_moments_file

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


   out_hgt = -9999.
   out_angles = -9999.
   out_tb = -9999


!   if (write_nc .eqv. .false.) call mod_io_strings_get_filename()


  if (verbose .gt. 1) print*, 'Start loop over frequencies & profiles!'


  grid_f: do fi =1, nfrq
     grid_y: do ny = 1, ngridy !ny_in, ny_fin  
        grid_x: do nx = 1, ngridx !nx_in, nx_fin   

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

          iwv = in_iwv(nx,ny)
          cwp = in_cwp(nx,ny)
          iwp = in_iwp(nx,ny)
          rwp = in_rwp(nx,ny)
          swp = in_swp(nx,ny)
          gwp = in_gwp(nx,ny)
          hwp = in_hwp(nx,ny)

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
           call run_rt3(nx,ny,fi,freqs(fi),freq_str)

          if (active) then
            out_Ze(nx,ny,1:nlyr,:) = REAL(Ze(nx,ny,1:nlyr,:))
            out_Att_hydro(nx,ny,1:nlyr,:) = REAL(Att_hydro(nx,ny,1:nlyr,:))
            out_Att_atmo(nx,ny,1:nlyr,:) = REAL(Att_atmo(nx,ny,1:nlyr,:))

            out_Ze_cw(nx,ny,1:nlyr,:) = REAL(Ze_cw(nx,ny,1:nlyr,:))
            out_Ze_rr(nx,ny,1:nlyr,:) = REAL(Ze_rr(nx,ny,1:nlyr,:))
            out_Ze_ci(nx,ny,1:nlyr,:) = REAL(Ze_ci(nx,ny,1:nlyr,:))
            out_Ze_sn(nx,ny,1:nlyr,:) = REAL(Ze_sn(nx,ny,1:nlyr,:))
            out_Ze_gr(nx,ny,1:nlyr,:) = REAL(Ze_gr(nx,ny,1:nlyr,:))
            out_Ze_ha(nx,ny,1:nlyr,:) = REAL(Ze_ha(nx,ny,1:nlyr,:))
            out_Att_cw(nx,ny,1:nlyr,:) = REAL(Att_cw(nx,ny,1:nlyr,:))
            out_Att_rr(nx,ny,1:nlyr,:) = REAL(Att_rr(nx,ny,1:nlyr,:))
            out_Att_ci(nx,ny,1:nlyr,:) = REAL(Att_ci(nx,ny,1:nlyr,:))
            out_Att_sn(nx,ny,1:nlyr,:) = REAL(Att_sn(nx,ny,1:nlyr,:))
            out_Att_gr(nx,ny,1:nlyr,:) = REAL(Att_gr(nx,ny,1:nlyr,:))
            out_Att_ha(nx,ny,1:nlyr,:) = REAL(Att_ha(nx,ny,1:nlyr,:))


            out_hgt(nx,ny,1:nlyr) = REAL(hgt(nx,ny,1:nlyr))
          end if

          call deallocate_profile_vars()

        end do grid_x
     end do grid_y
  end do grid_f



  if (passive) then
    out_angles = REAL(angles_deg(:))
    out_tb = RESHAPE( REAL(tb), (/ngridx, ngridy, noutlevels, 2*nummu, nfrq,nstokes /),&
          ORDER = (/6,5,4,3,2,1/))
  end if


  call deallocate_output_vars()

end subroutine pyPamtraLib