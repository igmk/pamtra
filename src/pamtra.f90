program pamtra

  use kinds
  use constants
  use nml_params
  use vars_atmosphere
  use vars_output
  use double_moments_module

  !     Radiative transfer code to process COSMO-model derived profiles   
  !     The code read a full COSMO grid and compute for each profile the  
  !     radiative transfer for the given frequency                        
  !                                                                      
  !     This code is completely self contained and needs no databases     
  !     or lookup tables to run.  By convention, the quantities followed  
  !     by "_lev" are given at the layer heights while the quantitites
  !     w/o "_lev" are layer average quantities                           
  !                                                                       
  !                    **  RADTRAN I/O SPECIFICATIONS  **

  implicit none

  integer, parameter :: mxgridx = 880, &  ! max grid dimension in x (421)
                        mxgridy = 880, &  ! max grid dimension in y (461)
                         mxlyr = 50       ! max grid dimension in z

  integer, parameter :: nummu = 16,         & ! no. of observation angles
       ntheta_i = 2*nummu, &
       ntheta_s = 2*nummu

  integer, parameter :: MAXV = 64,   &
       MAXA = 32,   &
       MAXLAY = 200

  integer, parameter :: NSTOKES = 2

  integer, parameter :: NOUTLEVELS = 2

  integer :: i, j, k, jj, nf, nx, ny, nz, nlev,&
             nx_in, nx_fin, ny_in, ny_fin,     &
             offset1,offset2, length1, length2,&
             ise, imonth ! filehandle for the emissivity data

  integer :: MAXLEG, NLEGEN, NUMRAD, NLEGENcw, NLEGENci, NLEGENgr, &
       NLEGENsn, NLEGENrr, aziorder, NUMAZIMUTHS

  integer :: SRC_CODE     ! describes the type of radiation 
  ! 1: 
  ! 2: 

  integer :: NP_LUT, N_temp

  integer :: j_temp, ind_temp, N_lay_cut

  integer :: i_bot, i_top, nnz, N_layer_new, ss, tt, length

!  integer, dimension(mxgridx,mxgridy) :: isamp, jsamp ! temporary for naming output

  integer, dimension(mxgridx,mxgridy) :: ics

  integer, dimension(maxlay) :: OUTLEVELS

  real(kind=dbl) :: DIRECT_FLUX, DIRECT_MU

  real(kind=dbl) :: freq,   & ! frequency [GHz]
             lam,    & ! wavelength [mm] 
             gammln    

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
       asymgr, salbhl, asymhl, backcw, backrr, backci, backsn, backgr,   &
       tau_min, tau_max                                                  

  real(kind=dbl) :: mu, D_0, D_max, D_min, A1, A2

  real(kind=dbl) :: GROUND_TEMP, ground_albedo

  real(kind=dbl) :: SKY_TEMP         ! cosmic background

  real(kind=dbl) :: WAVELENGTH

  real(kind=dbl) :: salinity         ! sea surface salinity

!  real(kind=dbl) :: emissivity     ! land surface reflectivity

  real(kind=dbl) :: Angle_view, Angle_zenith, I1, I2,     &
       Angle_viewdeg, Upar, Vpar, Ds_cloud_slant, Ds_obs_point,&
       DH_bot_intersect, H_bot, H_top, DH_top_intersect, K_extbot,&
       K_exttop, T_bot, T_top                                            

  real(kind=dbl), dimension(200) :: LEGEN, LEGEN2, LEGEN3, LEGEN4,&
       LEGENcw, LEGENrr, LEGENci, LEGENgr, LEGENsn,       &
       LEGEN2cw, LEGEN2rr, LEGEN2ci, LEGEN2gr, LEGEN2sn,  &
       LEGEN3cw, LEGEN3rr, LEGEN3ci, LEGEN3gr, LEGEN3sn,  &
       LEGEN4cw, LEGEN4rr, LEGEN4ci, LEGEN4gr, LEGEN4sn

  real(kind=dbl), dimension(mxgridx,mxgridy) :: tskin ! surface temperature

  real(kind=dbl), dimension(mxlyr) :: KEXTATMO   ! atmospheric extinction due to gases
!  real(kind=dbl), dimension(mxlyr) :: abscoef_o2,abscoef_h2o,abscoef_n2
  real(kind=dbl), dimension(2) :: P11, ang

  real(kind=dbl), dimension(mxgridx,mxgridy,mxlyr) :: &
       g_coeff,    &
       kexttot,    &
       kextcloud,  &
       kextrain,   &
       kextice,    &
       kextgraupel,&
       kextsnow,   &
       salbtot,    &
       absorp,     & ! might be unnecessary
       asymtot,    &
       back        

  real(kind=dbl), dimension(mxgridx, mxgridy) :: &
       tau,       &
       tau_hydro

  real(kind=dbl), dimension(maxv) :: MU_VALUES

  real(kind=dbl), dimension(ntheta_i) :: TTheta_i, TTheta_s

  real(kind=dbl), dimension(mxlyr+1) :: H_levs

  complex(kind=dbl) :: GROUND_INDEX
  complex(kind=dbl) :: eps_water, & ! function to calculate the dielectic properties of (salt)water
		epsi         ! result of function eps_water
  complex(kind=dbl) :: MINDEX, m_air, m_MG, m_ice

  character(2) :: month, day
  character(4) :: year, time
  character(12) :: date_str

  character :: QUAD_TYPE*1, UNITS*1, OUTPOL*2, GROUND_TYPE*1,  &
       rLWC_str*4                                                      

  character(100) :: OUT_FILE, tmp_file1, nc_out_file

  character(64) :: file_PH(mxlyr), file_PH2(mxlyr)

  character(74) :: file_profile

  character :: ssstr*1, ttstr*1, Anglestr*4, FILEOUT3D*65

  character :: Nzstr*2, xstr*3, ystr*3,&
       frq_str*6, theta_str*3, H_str*3,&
       surf_type*10

  character(99) :: input_file

  character :: micro_str*31, SP_str*3, str1*1,  &
       DELTAM*1, file_profile2*78, &
        N0snowstr*3, N0graustr*3, N0rainstr*3

  character(100) :: FILEOUT1

  character(80) :: femis ! filename for the emissivity databases

  character(20) :: moments_filein, dummy

  ! namelist parameters

!   real(kind=dbl) :: obs_height     ! upper level output height [m] (> 100000. for satellite)
!  
!   logical :: lphase_flag, &        ! flag for phase function calculation
! 	     lgas_extinction, &    ! gas extinction desired
! 	     lhyd_extinction, &       ! no hydrometeors desired
! 	     grid_calc, &	   ! calculate with dispatch or not
! 	     write_nc		   ! write netcdf output
! 
!   character(3) :: gas_mod

  !     The following parameters are set for the TMI sensor

!  real(kind=dbl), parameter :: pi = 3.141592653589793

!  complex(kind=dbl), parameter :: Im = (0.0, 1.0)

! temporary variables

  real :: lat, lon, lfrac, wind10

  real(kind=dbl) :: lwc, iwc, rwc, gwc, swc

  real(kind=dbl) :: spec2abs

  ! name list declarations
 
  namelist / verbose_mode / verbose
  namelist / process_mode / grid_calc
  namelist / output_mode / write_nc
  namelist / output / obs_height,units,outpol
  namelist / surface_params / ground_type,salinity, emissivity
  namelist / gas_abs_mod / lgas_extinction, gas_mod
  namelist / hyd_opts / lhyd_extinction, lphase_flag
  namelist / snow_params / SD_snow, N_0snowDsnow, EM_snow, SP
  namelist / graupel_params / SD_grau, N_0grauDgrau, EM_grau
  namelist / rain_params / SD_rain, N_0rainD
  namelist / moments / n_moments

  !    some inputs variable                                               
  N_lay_cut = 50
  QUAD_TYPE = 'L'                                             
  Aziorder = 0 
  NUMAZIMUTHS = 1 
  DELTAM = 'N' 
  SRC_CODE = 2 
  DIRECT_FLUX = 0.d0 
  DIRECT_MU = 0.0d0 
  maxleg = 200 
  SKY_TEMP = 2.7 
  ! WAVELENGTH = 1000000.0 ! check this

  ! initialize values

  tau = 0.0d0
  tau_hydro = 0.0d0 
  file_ph2(:) = ''

  ! read name list parameter file

  open(7, file='run_params.nml',delim='APOSTROPHE')
  read(7,nml=verbose_mode)
  read(7,nml=process_mode)
  read(7,nml=output_mode)
  read(7,nml=output)
  read(7,nml=surface_params)
  read(7,nml=gas_abs_mod)
  read(7,nml=hyd_opts)
  read(7,nml=snow_params)
  read(7,nml=graupel_params)
  read(7,nml=rain_params)
  read(7,nml=moments)
  close(7)

  !
  !     Get input/output file names from run file.
  !
  read ( *, * ) input_file
  read ( *, * ) nx_in, nx_fin, ny_in, ny_fin, tau_min, tau_max
  read ( *, * ) freq               ! frequency [Ghz]
  if (n_moments .eq. 2) read ( *, * ) moments_filein
  lam = 299.7925 / freq   ! mm

  if (n_moments .eq. 2) then
    open(118,file=moments_filein)
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

  if (grid_calc) then
    open(UNIT=14, FILE='/net/roumet/mech/pamtra/profiles/'//input_file, STATUS='OLD', form='formatted')
  else
    open(UNIT=14, FILE='profiles/'//input_file, STATUS='OLD', form='formatted')
!    open(UNIT=14, FILE='/work/mech/pamtra/profiles/'//input_file, STATUS='OLD', form='formatted')
  end if

  read(14,*) year, month, day, time, ngridx, ngridy, nlyr, deltax, deltay

  n_lay_cut = nlyr  ! in future n_lay_cut can be delete  $##

  call allocate_vars_atmosphere

  ! $## think about order of reading
  do i = 1, ngridx
     do j = 1, ngridy 
        read(14,*) profiles(i,j)%isamp, profiles(i,j)%jsamp ! 
        read(14,*) &
           profiles(i,j)%latitude, &	     ! �
		   profiles(i,j)%longitude,&         ! �
		   profiles(i,j)%land_fraction,&     !
		   profiles(i,j)%wind_10m            ! m/s
        ! integrated quantities
        if (n_moments .eq. 1) then
           read(14,*) &
           profiles(i,j)%iwv,&               ! kg/m^2
		   profiles(i,j)%cwp,&               ! kg/m^2
		   profiles(i,j)%iwp,&               ! kg/m^2
		   profiles(i,j)%rwp,&               ! kg/m^2
		   profiles(i,j)%swp,&               ! kg/m^2
		   profiles(i,j)%gwp                 ! kg/m^2
		end if
		if (n_moments .eq. 2) then
           read(14,*) &
           profiles(i,j)%iwv,&               ! kg/m^2
		   profiles(i,j)%cwp,&               ! kg/m^2
		   profiles(i,j)%iwp,&               ! kg/m^2
		   profiles(i,j)%rwp,&               ! kg/m^2
		   profiles(i,j)%swp,&               ! kg/m^2
		   profiles(i,j)%gwp,&               ! kg/m^2
		   profiles(i,j)%hwp                 ! kg/m^2
		end if

        ! surface values
        read(14,*) &
           profiles(i,j)%hgt_lev(0),&
		   profiles(i,j)%press_lev(0),&
		   profiles(i,j)%temp_lev(0),&
		   profiles(i,j)%relhum_lev(0)
        do k = 1, nlyr 
		   if (n_moments .eq. 1) then
              read(14,*) &
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
		      read(14,*) &
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
		end do
     end do
  end do

  close(14)

  if (verbose .gt. 0) print*, 'profile reading done!'

  if (write_nc) then
    allocate(is(ngridy,ngridx),js(ngridy,ngridx))
    allocate(lons(ngridy,ngridx),lats(ngridy,ngridx),lfracs(ngridy,ngridx))
    allocate(t_g(ngridy,ngridx),w10s(ngridy,ngridx),iwvs(ngridy,ngridx))
    allocate(cwps(ngridy,ngridx),iwps(ngridy,ngridx),rwps(ngridy,ngridx),swps(ngridy,ngridx),gwps(ngridy,ngridx))
    allocate(flux_up(nstokes,noutlevels,ngridy,ngridx),flux_down(nstokes,noutlevels,ngridy,ngridx))
    allocate(tb(nstokes,2*nummu,noutlevels,ngridy,ngridx))
  end if

  !                                                                       
  !     This GCE model format does not have all the fields expected by    
  !     the radiative transfer code (i.e. total pressure, and water vapor 
  !     pressure for this model).  Assign/compute the missing fields first
  !

  call get_atmosG0

  if (verbose .gt. 0) print*, 'variables filled up!'

  write (frq_str, '(f6.2)') freq
  frq_str = repeat('0', 6-len_trim(adjustl(frq_str)))//adjustl(frq_str)

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

  if (verbose .gt. 0) print*, 'Start loop over profiles!'

  grid_y: do ny = 1, ngridy !ny_in, ny_fin  
    grid_x: do nx = 1, ngridx !nx_in, nx_fin   

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
	  hwc_q = profiles(nx,ny)%hail_q				! kg/kg
	  cwc_n = profiles(nx,ny)%cloud_water_n			! #/kg
	  iwc_n = profiles(nx,ny)%cloud_ice_n			! #/kg
	  rwc_n = profiles(nx,ny)%rain_n				! #/kg
	  swc_n = profiles(nx,ny)%snow_n				! #/kg
	  gwc_n = profiles(nx,ny)%graupel_n				! #/kg
	  hwc_n = profiles(nx,ny)%hail_n				! #/kg
	end if

	press = profiles(nx,ny)%press                   ! Pa
	temp = profiles(nx,ny)%temp                     ! K
	relhum = profiles(nx,ny)%relhum                 ! %
	vapor_pressure = profiles(nx,ny)%vapor_pressure ! Pa
	rho_vap = profiles(nx,ny)%rho_vap               ! kg/m^3
	q_hum = profiles(nx,ny)%q_hum                   ! kg/kg

	if (verbose .gt. 0) print*, 'type to local variables done' 

	! Determine surface properties

	if (lfrac .ge. 0.5 .and. lfrac .le. 1.0) then
	  ground_type = 'L'
	  ise=13
	  read(month,'(i2)') imonth
	  if (imonth .ge. 7 .and. imonth .le. 12) then
	    femis = 'data/emissivity/ssmi_mean_emis_92'//month//'_direct'
	  else if (imonth .ge. 1 .and. imonth .lt. 7) then
	    femis = 'data/emissivity/ssmi_mean_emis_93'//month//'_direct'
	  else
	    print*, "Warning: No emissivity data found!"
	    stop
	  end if
      open(ise,file=trim(femis),status='old',form='unformatted',&
                access='direct',recl=28)
	  ! land_emis could give polarized reflectivities
      call land_emis(ise,lon,lat,real(freq),emissivity)
	  close(ise)
	  ground_albedo = 1. - emissivity
	else if (lfrac .ge. 0.0 .and. lfrac .lt. 0.5) then
      ! computing the refractive index of the sea (Fresnel) surface
	  ground_type = 'F'
	  ground_albedo = 1.0
	  epsi = eps_water(salinity, ground_temp - 273.15, freq)
	  ground_index = dconjg(sqrt(epsi))
	else
	! this is for specular reflection (only testing)
	  ground_type = 'S'
	  ground_albedo = 1 - emissivity
	end if

	if (verbose .gt. 0) print*, 'Surface emissivity calculated!'

    ! gaseous absorption
	! 
	! kextatmo   extinction by moist air [Np/m]
    !

    if (lgas_extinction) then
      call get_atmosg(press, temp, vapor_pressure, rho_vap, nlyr, freq, gas_mod, kextatmo)!,abscoef_o2,abscoef_h2o,abscoef_n2)
    else
      kextatmo = 0.0D0 ! for the whole column
    end if

	if (verbose .gt. 0) print*, 'Gas absorption calculated'

    write(xstr, '(i3.3)') profiles(nx,ny)%isamp
    write(ystr, '(i3.3)') profiles(nx,ny)%jsamp
    file_profile = '/tmp/Profilex'//xstr//'y'//ystr//'f'//frq_str

	! hydrometeor extinction desired

	if (lhyd_extinction) call hydrometeor_extinction(freq, n_lay_cut,xstr,ystr,frq_str,file_ph)

    !      Preparation of the PROFILE file  (needed by RT3)
    open(21, file = file_profile, form = 'FORMATTED', status =  &
         'unknown')
 
    do nz = N_lay_cut, 1, - 1 !nlyr,1,-1
       str1 = ''''
       ! position of the first blank space
       offset1 = index(FILE_PH(nz) , ' ')
       tmp_file1 = FILE_PH(nz)
       FILE_PH2(nz) = str1//tmp_file1(1:offset1 - 1)//str1
       write(21,1013) hgt_lev(nz), temp_lev(nz), kextatmo(nz), FILE_PH2(nz)
1013   format(f7.1,1x,f6.2,1x,E9.4,1x,a38)
    end do !end of cycle over the vertical layers
    write(21,1012) hgt_lev(0) , temp_lev(0), KEXTATMO (1) , ''' '''
1012 format(f7.1,1x,f6.2,1x,E9.4,1x,a3)
        close(21)


        file_profile2 = micro_str//'date'//date_str//'x'//xstr//    &
             'y'//ystr//'f'//frq_str                                     


1110    format   (i3,1x,i3,1x,4(1x,f7.3))
1111    format   (i3,1x,i3,1x,i2,1x,f6.3,10(1x,e9.3),1x,e9.4,2(1x,f7.4),  &
             & 1x,e9.3,1x,e9.4,1x,f5.1,1x,e9.3)                                 


        !&&&&&&&&   I/O FILE NAMES for the MC&&&&&&&&&&&&&&&&&&                 

!        FILEOUT1 = '/work/mech/pamtra/output/RT3TB'//date_str//micro_str//'x'//xstr//'y'//&
!             ystr//'f'//frq_str
       FILEOUT1 = 'output/RT3TB'//date_str//micro_str//'x'//xstr//'y'//&
            ystr//'f'//frq_str

        OUT_FILE = FILEOUT1


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

      wavelength = lam*1.e3    ! in microns

      if (verbose .gt. 0) print*, "Entering rt3 ...."
  
      call RT3(NSTOKES, NUMMU, AZIORDER, MU_VALUES, src_code,     &
	    FILE_profile, out_file, QUAD_TYPE, deltam, DIRECT_FLUX,     &
	    DIRECT_MU, GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO,         &
	    GROUND_INDEX, SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,          &
	    NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS,&
	    nx,ny,write_nc,verbose)

      if (verbose .gt. 0) print*, "....rt3 finished"

     end do grid_x
  end do grid_y

  if (write_nc) then
    nc_out_file = trim(input_file(1:len_trim(input_file)-4))//'_'//frq_str//'_res.nc'
    call write_nc_results(nc_out_file)
  end if

end program pamtra
