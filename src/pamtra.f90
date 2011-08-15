program pamtra

  use kinds
  use constants
  use nml_params
  use vars_atmosphere
  use vars_output
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

  integer, parameter :: nummu = 16,         & ! no. of observation angles
       ntheta_i = 2*nummu, &
       ntheta_s = 2*nummu

  integer, parameter :: MAXV = 64,   &
       MAXA = 32,   &
       MAXLAY = 200, &
       maxleg = 200, &
       maxfreq = 40

  integer, parameter :: NSTOKES = 2

  integer, parameter :: NOUTLEVELS = 2

  integer :: i, j, k, jj, nf, nx, ny, nz, nlev,&
             offset1,offset2, length1, length2,&
             inarg, fi, &
             ise, imonth ! filehandle for the emissivity data

  integer :: NUMRAD, NLEGENcw, NLEGENci, NLEGENgr, &
       NLEGENsn, NLEGENrr, aziorder, NUMAZIMUTHS

  integer :: SRC_CODE     ! describes the type of radiation 
  ! 1: 
  ! 2: 

  integer :: NP_LUT, N_temp

  integer :: j_temp, ind_temp

  integer :: i_bot, i_top, nnz, N_layer_new, ss, tt, length


  integer, dimension(maxlay) :: OUTLEVELS

  real(kind=dbl) :: DIRECT_FLUX, DIRECT_MU

   real(kind=dbl) :: freq,   & ! frequency [GHz]
 	             gammln
  character(6), dimension(maxfreq) :: frqs_str



  real(kind=dbl) :: E1, E2

  real(kind=dbl) :: RAD1, RAD2, refre, refim,&
       N_0sr, Coeff_corr, AD, BD,  &
       n0S, lambda_D, tmp, N_0snowD, N_0grauD,              &
       Coeff_snow, a_mgraup, b_g, b_snow,        &
       a_msnow, Coeff_grau, den_liq, drop_mass, del_r, den_ice,&
       densnow, fvol_ice

  real(kind=dbl) :: deltax, deltay, ABSIND, ABSCOF

  real(kind=dbl) :: atm_ext, kextcw, salbcw, asymcw, kextrr, salbrr, asymrr,  &
       kextci, salbci, asymci, kextsn, salbsn, asymsn, kextgr, salbgr,   &
       asymgr, salbhl, asymhl, backcw, backrr, backci, backsn, backgr
                                                      

  real(kind=dbl) :: mu, D_0, D_max, D_min, A1, A2

  real(kind=dbl) :: GROUND_TEMP, ground_albedo

  real(kind=dbl) :: SKY_TEMP         ! cosmic background

  real(kind=dbl) :: wavelength       ! microns

  real(kind=dbl) :: salinity         ! sea surface salinity

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

  real(kind=dbl), dimension(ntheta_i) :: TTheta_i, TTheta_s


  complex(kind=dbl) :: GROUND_INDEX
  complex(kind=dbl) :: eps_water, & ! function to calculate the dielectic properties of (salt)water
		epsi         ! result of function eps_water
  complex(kind=dbl) :: MINDEX, m_air, m_MG, m_ice

  character :: QUAD_TYPE*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1,  &
       rLWC_str*4                                                      

  character(300) :: OUT_FILE_PAS, OUT_FILE_ACT, tmp_file1, nc_out_file, namelist_file

  character :: ssstr*1, ttstr*1, Anglestr*4, FILEOUT3D*65

  character(6) :: formatted_frqstr
  character(5*7) :: frq_str_list
  character(99) :: input_file

  character :: micro_str*31, SP_str*3,  &
       DELTAM*1, file_profile2*78, &
        N0snowstr*3, N0graustr*3, N0rainstr*3

  character(80) :: femis ! filename for the emissivity databases

  character(20) :: dummy

  integer :: istat


! temporary variables

  real(kind=sgl) :: lat, lon, lfrac
  real(kind=dbl) ::  wind10u, wind10v

  real(kind=dbl) :: lwc, iwc, rwc, gwc, swc

  real(kind=dbl) :: spec2abs

  real(kind=dbl), dimension(2) :: P11, ang

!  integer, dimension(ngridx,ngridy) :: isamp, jsamp ! temporary for naming output

  integer, allocatable, dimension(:,:) :: ics

 ! name list declarations
 
  namelist / verbose_mode / verbose
  namelist / inoutput_mode / write_nc, input_path, output_path, tmp_path, dump_to_file
  namelist / output / obs_height,units,outpol,creator
  namelist / run_mode / active, passive
  namelist / surface_params / ground_type,salinity, emissivity
  namelist / gas_abs_mod / lgas_extinction, gas_mod
  namelist / hyd_opts / lhyd_extinction, lphase_flag
  namelist / snow_params / SD_snow, N_0snowDsnow, EM_snow, SP, isnow_n0, liu_type
  namelist / graupel_params / SD_grau, N_0grauDgrau, EM_grau
  namelist / ice_params / EM_ice
  namelist / rain_params / SD_rain, N_0rainD
  namelist / moments / n_moments, moments_file

!parse command line parameters
inarg = iargc()

if (inarg .lt. 3) then
   print *,'Usage: pamtra profile_file namelist_file (list of frequencies)'
   print *,'Example: ./pamtra rt_comp_single.dat run_params.nml 35 94'
   print *,'See namelist file for further pamtra options'
   stop
else if (inarg .gt. maxfreq + 2) then
   print *,'Too many frequencies! Increase maxfreq!'
   stop
end if
   call getarg(1,input_file)
   call getarg(2,namelist_file)

nfrq = inarg - 2
allocate(freqs(nfrq))

do fi = 1, inarg-2
    call getarg(fi+2,frqs_str(fi))
    read(frqs_str(fi),*) freqs(fi)
    frqs_str(fi) = formatted_frqstr(frqs_str(fi))
end do


  ! read name list parameter file

  open(7, file=namelist_file,delim='APOSTROPHE')
  read(7,nml=verbose_mode)
  read(7,nml=inoutput_mode)
  read(7,nml=output)
  read(7,nml=run_mode)
  read(7,nml=surface_params)
  read(7,nml=gas_abs_mod)
  read(7,nml=hyd_opts)
  read(7,nml=snow_params)
  read(7,nml=graupel_params)
  read(7,nml=ice_params)
  read(7,nml=rain_params)
  read(7,nml=moments)
  close(7)
  
  if (n_moments .ne. 1 .and. n_moments .ne. 2) stop "n_moments is not 1 or 2"



  if (n_moments .eq. 2) then
    open(118,file=moments_file)
    !read NU & Mu parameter of the drop size distr. as a function of MASS
    !and Alpha & Beta parameter of Diameter-Mass function
    !gamma_xxx(1)=nu	gamma_xxx(2)=mu		gamma_xxx(3)=alpha		gamma_xxx(4)=beta
    read(118,'(a20,4(x,d13.6))') dummy,gamma_cloud
    read(118,'(a20,4(x,d13.6))') dummy,gamma_rain
    read(118,'(a20,4(x,d13.6))') dummy,gamma_ice
    read(118,'(a20,4(x,d13.6))') dummy,gamma_snow
    read(118,'(a20,4(x,d13.6))') dummy,gamma_graupel
    read(118,'(a20,4(x,d13.6))') dummy,gamma_hail
	close(118)
  end if

  !                                                                       
  !     read atmospheric profiles                 
  !  
  !     quantities followed  
  !     by "_lev" are given at the layer heights while the quantitites    
  !     w/o "_lev" are layer average quantities
  !
  !  Variables:
  !     name           unit        description
  !     hgt_lev         m
  !     press_lev       Pa
  !     temp_lev        K
  !     relhum_lev      %
  !
  !     cloud_water_q kg/kg
  !     cloud_ice_q   kg/kg 
  !     rain_q        kg/kg
  !     snow_q        kg/kg
  !     graupel_q     kg/kg

  ! Perform calculation through dispatch or not

   if (verbose .gt. 1) print *,"input_file: ",input_file(:len_trim(input_file))," namelist file: ",namelist_file," freq: ",frqs_str

   if (verbose .gt. 1) print *,"PASSIVE: ", passive, "ACTIVE: ", active

   if (verbose .gt. 0) print *,"opening: ",input_path(:len_trim(input_path))//"/"//input_file

  open(UNIT=14, FILE=input_path(:len_trim(input_path))//"/"//input_file, STATUS='OLD', form='formatted',iostat=istat)
 
 if (istat .ne. 0) then
    call error_msg(input_file,0,0)
    stop "Read error 1: Cannot open file"
    end if

  read(14,*,iostat=istat) year, month, day, time, ngridx, ngridy, nlyr, deltax, deltay

  if (istat .ne. 0) then
    call error_msg(input_file,0,0)
    stop "Read error 2: Cannot read first line of file"
    end if
	! now allocate variables



	allocate(ics(ngridx, ngridy))
    allocate(file_ph(nlyr))

  !    some inputs variable                                               
  QUAD_TYPE = 'L'                                             
  Aziorder = 0 
  NUMAZIMUTHS = 1 
  DELTAM = 'N' 
  SRC_CODE = 2 
  DIRECT_FLUX = 0.d0 
  DIRECT_MU = 0.0d0 
  SKY_TEMP = 2.73d0

  ! initialize values

!   tau = 0.0d0
!   tau_hydro = 0.0d0 


  call allocate_vars_atmosphere

  ! $## think about order of reading
  do i = 1, ngridx
     do j = 1, ngridy 
        read(14,*,iostat=istat) profiles(i,j)%isamp, profiles(i,j)%jsamp !
        if (istat .ne. 0) then
            call error_msg(input_file,i,j)
            stop "Read error 3: Cannot read profile index i,j"
        end if

        read(14,*,iostat=istat) &
           profiles(i,j)%latitude, &	     ! degree
		   profiles(i,j)%longitude,&         ! degree
		   profiles(i,j)%land_fraction,&     !
		   profiles(i,j)%wind_10u,&          ! m/s
		   profiles(i,j)%wind_10v            ! m/s
        if (istat .ne. 0) then
            call error_msg(input_file,i,j)
            stop "Read error 4: Cannot read profile lat/lon/lfrac/wind"
        end if
        ! integrated quantities
        if (n_moments .eq. 1) then
           read(14,*,iostat=istat) &
           profiles(i,j)%iwv,&               ! kg/m^2
		   profiles(i,j)%cwp,&               ! kg/m^2
		   profiles(i,j)%iwp,&               ! kg/m^2
		   profiles(i,j)%rwp,&               ! kg/m^2
		   profiles(i,j)%swp,&               ! kg/m^2
		   profiles(i,j)%gwp                 ! kg/m^2
		   profiles(i,j)%hwp = 0.
		end if
		if (n_moments .eq. 2) then
           read(14,*,iostat=istat) &
           profiles(i,j)%iwv,&               ! kg/m^2
		   profiles(i,j)%cwp,&               ! kg/m^2
		   profiles(i,j)%iwp,&               ! kg/m^2
		   profiles(i,j)%rwp,&               ! kg/m^2
		   profiles(i,j)%swp,&               ! kg/m^2
		   profiles(i,j)%gwp,&               ! kg/m^2
		   profiles(i,j)%hwp                 ! kg/m^2
		end if

        if (istat .ne. 0) then
            call error_msg(input_file,i,j)
            stop "Read error 5: Cannot read profile integrated quantities"
        end if

        ! surface values
        read(14,*,iostat=istat) &
           profiles(i,j)%hgt_lev(0),&
		   profiles(i,j)%press_lev(0),&
		   profiles(i,j)%temp_lev(0),&
		   profiles(i,j)%relhum_lev(0)
	    if (istat .ne. 0) call error_msg(input_file, i, j)
        do k = 1, nlyr 
		   if (n_moments .eq. 1) then
              read(14,*,iostat=istat) &

              profiles(i,j)%hgt_lev(k), &             ! m
		      profiles(i,j)%press_lev(k), &           ! Pa
		      profiles(i,j)%temp_lev(k), &            ! K
		      profiles(i,j)%relhum_lev(k), &          ! %
		      profiles(i,j)%cloud_water_q(k), &       ! kg/kg
		      profiles(i,j)%cloud_ice_q(k), &         ! kg/kg
		      profiles(i,j)%rain_q(k), &              ! kg/kg
		      profiles(i,j)%snow_q(k), &              ! kg/kg
		      profiles(i,j)%graupel_q(k)              ! kg/kg
		    end if
			if (n_moments .eq. 2) then
		      read(14,*,iostat=istat) &
              profiles(i,j)%hgt_lev(k), &             ! m
		      profiles(i,j)%press_lev(k), &           ! Pa
		      profiles(i,j)%temp_lev(k), &            ! K
		      profiles(i,j)%relhum_lev(k), &          ! %
		      profiles(i,j)%cloud_water_q(k), &       ! kg/kg
		      profiles(i,j)%cloud_ice_q(k), &         ! kg/kg
		      profiles(i,j)%rain_q(k), &              ! kg/kg
		      profiles(i,j)%snow_q(k), &              ! kg/kg
		      profiles(i,j)%graupel_q(k), &           ! kg/kg
		      profiles(i,j)%hail_q(k), &       		  ! kg/kg
		      profiles(i,j)%cloud_water_n(k), &       ! #/kg
		      profiles(i,j)%cloud_ice_n(k), &         ! #/kg
		      profiles(i,j)%rain_n(k), &              ! #/kg
		      profiles(i,j)%snow_n(k), &              ! #/kg
		      profiles(i,j)%graupel_n(k), &           ! #/kg
		      profiles(i,j)%hail_n(k)                 ! #/kg
	        end if
        if (istat .ne. 0) then
            call error_msg(input_file,i,j,k)
            print *,"Error at layer", k
            stop "Read error 6: Cannot read profile values"
        end if
		end do
     end do
  end do

  close(14)



  if (verbose .gt. 0) print*, 'profile reading done!'

  if (write_nc) then
    allocate(is(ngridy,ngridx),js(ngridy,ngridx))
    allocate(lons(ngridy,ngridx),lats(ngridy,ngridx),lfracs(ngridy,ngridx))
    allocate(iwvs(ngridy,ngridx))
    allocate(cwps(ngridy,ngridx),iwps(ngridy,ngridx),rwps(ngridy,ngridx),&
    swps(ngridy,ngridx),gwps(ngridy,ngridx),hwps(ngridy,ngridx))
    allocate(tb(nstokes,nfrq,2*nummu,noutlevels,ngridy,ngridx))
    lons = 0.; lats = 0.; lfracs = 0.;
    iwvs = 0.; cwps = 0.; iwps = 0.; rwps = 0.; swps = 0.; gwps = 0.; hwps = 0.;
    tb = 0.

  end if

if (active) then
    allocate(Ze(ngridx,ngridy,nlyr,nfrq))
    allocate(PIA_hydro_bottomup(ngridx,ngridy,nlyr,nfrq))
    allocate(PIA_atmo_bottomup(ngridx,ngridy,nlyr,nfrq))
    allocate(PIA_hydro_topdown(ngridx,ngridy,nlyr,nfrq))
    allocate(PIA_atmo_topdown(ngridx,ngridy,nlyr,nfrq))
    allocate(hgt(ngridx,ngridy,nlyr))
end if
  !                                                                       
  !     This GCE model format does not have all the fields expected by    
  !     the radiative transfer code (i.e. total pressure, and water vapor 
  !     pressure for this model).  Assign/compute the missing fields first
  !

  !make layer averages
  call get_atmosG0
  if (verbose .gt. 1) print*, 'variables filled up!'

  write (SP_str (1:3) , '(f3.1)') SP

  date_str = year//month//day//time

  if (N_0snowDsnow .le. 9.95) then 
     write (N0snowstr, '(f3.1)') N_0snowDsnow 
  else 
     write (N0snowstr, '(f3.0)') N_0snowDsnow 
  end if

  if (N_0rainD .le. 9.95) then 
     write (N0rainstr, '(f3.1)') N_0rainD 
  else 
     write (N0rainstr, '(f3.0)') N_0rainD 
  end if
  if (N_0grauDgrau .le. 9.95) then 
     write (N0graustr, '(f3.1)') N_0grauDgrau 
  else 
     write (N0graustr, '(f3.0)') N_0grauDgrau 
  end if

  micro_str = SD_snow//N0snowstr//EM_snow//SP_str//SD_grau//        &
       N0graustr//EM_grau//SD_rain//N0rainstr                            

  frq_str_list = ""
  if (verbose .gt. 1) print*, 'Start loop over frequencies!'

grid_f: do fi =1, nfrq

  freq = freqs(fi)
  frq_str = frqs_str(fi)
  frq_str_list = frq_str_list(:len_trim(frq_str_list)) // "_" // frq_str
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
        femis = 'data/emissivity/ssmi_mean_emis_92'//month//'_direct'
    else if (imonth .ge. 1 .and. imonth .lt. 7) then
        femis = 'data/emissivity/ssmi_mean_emis_93'//month//'_direct'
    else
        print*, nx,ny, "Warning: No emissivity data found!"
        stop
    end if
    open(ise,file=trim(femis),status='old',form='unformatted',&
                access='direct',recl=28) ! on cliot recl=7
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
        date_str//micro_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_passive"

        OUT_FILE_ACT = output_path(:len_trim(output_path))//"/"//&
        date_str//micro_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_active"


    if (active) then
        call calculate_active(OUT_FILE_ACT,freq,hgt(nx,ny,:),Ze(nx,ny,:,fi),PIA_atmo_bottomup(nx,ny,:,fi),&
            PIA_hydro_bottomup(nx,ny,:,fi), PIA_atmo_topdown(nx,ny,:,fi),PIA_hydro_topdown(nx,ny,:,fi))
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
