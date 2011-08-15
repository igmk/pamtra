program pamtra

  use kinds
  use constants
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
  use double_moments_module 
  use mod_io_strings

  !     Radiative transfer code to process COSMO-model derived profiles   
  !     The code reads a full COSMO grid and computes for each profile the  
  !     radiative transfer for the given frequency                        
  !                                                                      
  !     This code is completely self contained and needs no databases     
  !     or lookup tables to run.  By convention, the quantities followed  
  !     by "_lev" are given at the layer heights while the quantitites
  !     w/o "_lev" are layer average quantities                           
  !                                                                       
  !                    **  RADTRAN I/O SPECIFICATIONS  **

  implicit none





!!! NOT YET SORTED STUFF

  integer :: jj, nf, nx, ny, nz, nlev, fi,&
             offset1,offset2, length1, length2,&
             ise, imonth ! filehandle for the emissivity data

  integer :: NUMRAD, NLEGENcw, NLEGENci, NLEGENgr, &
       NLEGENsn, NLEGENrr


  integer :: NP_LUT, N_temp

  integer :: j_temp, ind_temp

  integer :: i_bot, i_top, nnz, N_layer_new, ss, tt, length


  integer, dimension(maxlay) :: OUTLEVELS



   real(kind=dbl) :: freq,   & ! frequency [GHz]
              gammln




  real(kind=dbl) :: E1, E2

  real(kind=dbl) :: RAD1, RAD2, refre, refim,&
       N_0sr, Coeff_corr, AD, BD,  &
       n0S, lambda_D, tmp, N_0snowD, N_0grauD,              &
       Coeff_snow, a_mgraup, b_g, b_snow,        &
       a_msnow, Coeff_grau, den_liq, drop_mass, del_r, den_ice,&
       densnow, fvol_ice

  real(kind=dbl) ::  ABSIND, ABSCOF

  real(kind=dbl) :: atm_ext, kextcw, salbcw, asymcw, kextrr, salbrr, asymrr,  &
       kextci, salbci, asymci, kextsn, salbsn, asymsn, kextgr, salbgr,   &
       asymgr, salbhl, asymhl, backcw, backrr, backci, backsn, backgr
                                                      

  real(kind=dbl) :: mu, D_0, D_max, D_min, A1, A2

  real(kind=dbl) :: GROUND_TEMP, ground_albedo


  real(kind=dbl) :: wavelength       ! microns

!  real(kind=dbl) :: emissivity     ! land surface reflectivity

  real(kind=dbl) :: Angle_view, Angle_zenith, I1, I2,     &
       Angle_viewdeg, Upar, Vpar, Ds_cloud_slant, Ds_obs_point,&
       DH_bot_intersect, H_bot, H_top, DH_top_intersect, K_extbot,&
       K_exttop, T_bot, T_top                                            

  real(kind=dbl), dimension(maxleg) :: LEGENcw, LEGENrr, LEGENci, LEGENgr, LEGENsn,       &
       LEGEN2cw, LEGEN2rr, LEGEN2ci, LEGEN2gr, LEGEN2sn,  &
       LEGEN3cw, LEGEN3rr, LEGEN3ci, LEGEN3gr, LEGEN3sn,  &
       LEGEN4cw, LEGEN4rr, LEGEN4ci, LEGEN4gr, LEGEN4sn


  real(kind=dbl), dimension(maxv) :: MU_VALUES

  complex(kind=dbl) :: GROUND_INDEX
  complex(kind=dbl) :: eps_water, & ! function to calculate the dielectic properties of (salt)water
		epsi         ! result of function eps_water
  complex(kind=dbl) :: MINDEX, m_air, m_MG, m_ice

  character :: rLWC_str*4                                                      

  character(300) :: OUT_FILE_PAS, OUT_FILE_ACT, tmp_file1, nc_out_file

  character :: ssstr*1, ttstr*1, Anglestr*4, FILEOUT3D*65




  character(80) :: femis ! filename for the emissivity databases




! temporary variables

  real(kind=sgl) :: lat, lon, lfrac
  real(kind=dbl) ::  wind10u, wind10v

  real(kind=dbl) :: lwc, iwc, rwc, gwc, swc

  real(kind=dbl) :: spec2abs

!!! INTERNAL "HANDLE COMMAND LINE PARAMETERS" !!! 

  integer :: inarg, ff
  character(40) :: gitHash, gitVersion
  character(6) :: formatted_frqstr !function call

!!! SET BY "HANDLE COMMAND LINE PARAMETERS" !!! 

  character(99)  :: input_file !name of profile
  character(300) :: namelist_file
  character(6), dimension(maxfreq) :: frqs_str !from commandline
  character(5*7) :: frq_str_list ! for the filename only!


!!! SET BY mod_io_strings_get_file_name !!!





!!! HANDLE COMMAND LINE PARAMETERS !!!

!get git data
call versionNumber(gitVersion,gitHash)

!get command line parameters
inarg = iargc()

if (inarg .lt. 3) then
   print *,'Usage: pamtra profile_file namelist_file (list of frequencies)'
   print *,'Example: ./pamtra rt_comp_single.dat run_params.nml 35 94'
   print *,'See namelist file for further pamtra options'
   print *,''
   print *,'Version:  '//gitVersion
   print *,'Git Hash: '//gitHash
   stop
else if (inarg .gt. maxfreq + 2) then
   print *,'Too many frequencies! Increase maxfreq!'
   stop
end if
   call getarg(1,input_file)
   call getarg(2,namelist_file)

nfrq = inarg - 2
allocate(freqs(nfrq))

frq_str_list = "" 
!get integer and character frequencies
do ff = 1, inarg-2
    call getarg(ff+2,frqs_str(ff))
    read(frqs_str(ff),*) freqs(ff)
    frqs_str(ff) = formatted_frqstr(frqs_str(ff))
    frq_str_list = frq_str_list(:len_trim(frq_str_list)) // "_" //  frqs_str(ff)
end do


if (verbose .gt. 1) print *,"input_file: ",input_file(:len_trim(input_file)),&
                           " namelist file: ",namelist_file," freq: ",frqs_str

!!! READ NAMELIST FILE !!!
call nml_params_read(namelist_file) !from nml_params.f90

!!! READ MOMENTS FILE !!!
if (n_moments .eq. 2) call double_moments_module_read(moments_file) !from double_moments_module.f90

!!! READ PROFILES !!!
call vars_atmosphere_read_profile(input_file)

! now allocate variables
call allocate_vars


  !                                                                       
  !     This GCE model format does not have all the fields expected by    
  !     the radiative transfer code (i.e. total pressure, and water vapor 
  !     pressure for this model).  Assign/compute the missing fields first
  !make layer averages
  call get_atmosG0

  if (write_nc .eqv. .false.) call mod_io_strings_get_filename()

  if (verbose .gt. 1) print*, 'Start loop over frequencies!'

grid_f: do fi =1, nfrq

  freq = freqs(fi)
  frq_str = frqs_str(fi)

  wavelength = c / (freq*1.d3)   ! microns

  if (verbose .gt. 1) print*, 'Start loop over profiles!'

  grid_y: do ny = 1, ngridy !ny_in, ny_fin  
    grid_x: do nx = 1, ngridx !nx_in, nx_fin   

    write(xstr, '(i3.3)') profiles(nx,ny)%isamp
    write(ystr, '(i3.3)') profiles(nx,ny)%jsamp

    if (verbose .gt. 0) print*, "calculating: ", frq_str, " Y:",ny, " of ", ngridy, "X:", nx, " of ", ngridx

    ground_temp = profiles(nx,ny)%temp_lev(0)       ! K
    lat = profiles(nx,ny)%latitude                  ! °
    lon = profiles(nx,ny)%longitude                 ! °
    lfrac = profiles(nx,ny)%land_fraction
    relhum_lev = profiles(nx,ny)%relhum_lev         ! %
    press_lev = profiles(nx,ny)%press_lev           ! Pa
    temp_lev = profiles(nx,ny)%temp_lev             ! K
    hgt_lev = profiles(nx,ny)%hgt_lev               ! m

    cwc_q = profiles(nx,ny)%cloud_water_q           ! kg/kg
    iwc_q = profiles(nx,ny)%cloud_ice_q             ! kg/kg
    rwc_q = profiles(nx,ny)%rain_q                  ! kg/kg
    swc_q = profiles(nx,ny)%snow_q                  ! kg/kg
    gwc_q = profiles(nx,ny)%graupel_q               ! kg/kg

    if (n_moments .eq. 2) then
    hwc_q = profiles(nx,ny)%hail_q              ! kg/kg
    cwc_n = profiles(nx,ny)%cloud_water_n       ! #/kg
    iwc_n = profiles(nx,ny)%cloud_ice_n         ! #/kg
    rwc_n = profiles(nx,ny)%rain_n              ! #/kg
    swc_n = profiles(nx,ny)%snow_n              ! #/kg
    gwc_n = profiles(nx,ny)%graupel_n           ! #/kg
    hwc_n = profiles(nx,ny)%hail_n              ! #/kg
    end if

    press = profiles(nx,ny)%press                   ! Pa
    temp = profiles(nx,ny)%temp                     ! K
    relhum = profiles(nx,ny)%relhum                 ! %
    vapor_pressure = profiles(nx,ny)%vapor_pressure ! Pa
    rho_vap = profiles(nx,ny)%rho_vap               ! kg/m^3
    q_hum = profiles(nx,ny)%q_hum                   ! kg/kg


    if (verbose .gt. 1) print*, nx,ny, 'type to local variables done' 

    ! Determine surface properties

    if (lfrac .ge. 0.5 .and. lfrac .le. 1.0) then
    ground_type = 'S' ! changed to specular after advice of cathrine prigent
    ise=13
    read(month,'(i2)') imonth
    if (imonth .ge. 7 .and. imonth .le. 12) then
        femis = data_path(:len_trim(data_path))//'/emissivity/ssmi_mean_emis_92'//month//'_direct'
    else if (imonth .ge. 1 .and. imonth .lt. 7) then
        femis = data_path(:len_trim(data_path))//'emissivity/ssmi_mean_emis_93'//month//'_direct'
    else
        print*, nx,ny, "Warning: No emissivity data found!"
        stop
    end if
    open(ise,file=trim(femis),status='old',form='unformatted',&
                access='direct',recl=28)
	  ! land_emis could give polarized reflectivities
      call land_emis(ise,lon,lat,real(freq),emissivity)
	  close(ise)
	  ground_albedo = 1.d0 - emissivity
	else if (lfrac .ge. 0.0 .and. lfrac .lt. 0.5) then
      ! computing the refractive index of the sea (Fresnel) surface
	  ground_type = 'O'
	  ground_albedo = 1.0d0
	  epsi = eps_water(salinity, ground_temp - 273.15d0, freq)
	  ground_index = dconjg(sqrt(epsi))
	else
	! this is for ground_type specified in run_params.nml
	  ground_albedo = 1.d0 - emissivity
	end if

	if (verbose .gt. 1) print*, nx,ny, 'Surface emissivity calculated!'

    ! gaseous absorption
    ! 
    ! kextatmo   extinction by moist air [Np/m]
    !
    if (lgas_extinction) then
    !returns kextatmo!
    call get_atmosg(freq)
    else
    kextatmo = 0.0D0 ! for the whole column
    end if

    if (verbose .gt. 1) print*, nx,ny, 'Gas absorption calculated'


    ! hydrometeor extinction desired


    if (lhyd_extinction) call hydrometeor_extinction(freq)

!
    if (dump_to_file) call dump_profile

        !&&&&&&&&   I/O FILE NAMES   &&&&&&&&&&&&&&&&&&

        OUT_FILE_PAS = output_path(:len_trim(output_path))//"/"//&
        micro_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_passive"

        OUT_FILE_ACT = output_path(:len_trim(output_path))//"/"//&
        micro_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_active"


    if (active) then
        call calculate_active(OUT_FILE_ACT,freq,hgt(nx,ny,:),Ze(nx,ny,:,fi),Attenuation_atmo(nx,ny,:,fi),&
            Attenuation_hydro(nx,ny,:,fi))
        if (verbose .gt. 1) print*, nx,ny, 'calculate_active done'
    end if

	if (write_nc) then
		!      Output integrated quantities
		call collect_boundary_output(lon,lat,lfrac,&
			profiles(nx,ny)%iwv, profiles(nx,ny)%cwp,profiles(nx,ny)%iwp,profiles(nx,ny)%rwp,profiles(nx,ny)%swp, &
			profiles(nx,ny)%gwp,profiles(nx,ny)%hwp,profiles(nx,ny)%isamp,profiles(nx,ny)%jsamp,nx,ny)
        if (verbose .gt. 1) print*, nx,ny, 'collect_boundary_output done'
    end if

! find the output level
! in rt3 layers are reversed

    if (obs_height .gt. 99999. .or. obs_height .gt. hgt_lev(nlyr)) then
        outlevels(1) = 1
    else if (obs_height .lt. 0.1 .or. obs_height .lt. hgt_lev(1)) then
        outlevels(1) = nlyr + 1
    else
        out_search: do nz = 1, nlyr
        if (hgt_lev(nz) .ge. obs_height) then
            if (abs(hgt_lev(nz) - obs_height) .lt. abs(hgt_lev(nz-1) - obs_height)) then
            outlevels(1) = nlyr-nz+1
            else
            outlevels(1) = nlyr-nz+2
            end if
            exit out_search
        end if
        end do out_search
    end if

    OUTLEVELS(2) = nlyr+1    ! this is the bottom

    if (passive .eqv. .true.) then

    if (verbose .gt. 1) print*, nx,ny, "Entering rt3 ...."
    
    call RT3(NSTOKES, NUMMU, AZIORDER, MU_VALUES, src_code,     &
        out_file_pas, QUAD_TYPE, deltam, DIRECT_FLUX,     &
        DIRECT_MU, GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO,         &
        GROUND_INDEX, SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,          &
        NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS,&
        nx,ny,fi,write_nc,verbose)

    !calculate human readable angles!
    angles_deg(1:NUMMU) = 180-(180.*acos(MU_VALUES(NUMMU:1:-1))/pi)
    angles_deg(1+NUMMU:2*NUMMU) = (180.*acos(MU_VALUES(1:NUMMU))/pi)

    if (verbose .gt. 1) print*, nx,ny, "....rt3 finished"

    end if

    end do grid_x
  end do grid_y
end do grid_f

  if (write_nc) then
    nc_out_file = output_path(1:len_trim(output_path))//"/"//trim(input_file(1:len_trim(input_file)-4))//&
                    frq_str_list(1:len_trim(frq_str_list))//'_res.nc'
    if (verbose .gt. 0) print*,"writing: ", nc_out_file
    call write_nc_results(nc_out_file)
  end if

end program pamtra
