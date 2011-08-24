program pamtra

  use kinds
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
  use double_moments_module !double moments variables are stored here
  use mod_io_strings !some strings for nice filenames


  !     The code reads a full (e.g. COSMO) grid and computes for each 
  !     profile the radiative transfer for the given frequencies                       
  !                                                                      
  !     By convention, the quantities followed  by "_lev"
  !     are given at the layer heights while the quantitites w/o
  !     "_lev" are layer average quantities                           
  !                                                                       

  implicit none


!!! internal "handle command line parameters" !!! 

  integer :: inarg, ff
  character(40) :: gitHash, gitVersion
  character(6) :: formatted_frqstr !function call

!!! set by "handle command line parameters" !!! 

  character(99)  :: input_file !name of profile
  character(300) :: namelist_file
  character(6), dimension(maxfreq) :: frqs_str !from commandline
  character(7) :: frq_str_s,frq_str_e

!!!loop variables
  integer ::  fi,nx, ny

!!!output variables
  character(300) ::nc_out_file

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

  do ff = 1, nfrq
     call getarg(ff+2,frqs_str(ff))
    read(frqs_str(ff),*) freqs(ff)
     frqs_str(ff) = formatted_frqstr(frqs_str(ff))
  end do

!!! read variables from namelist file
  call nml_params_read(namelist_file) !from nml_params.f90

  ! create frequency string of not set in pamtra
  if (freq_str .eq. "") then
     ! get integer and character frequencies
    frq_str_s = "_"//frqs_str(1)
     if (nfrq .eq. 1) then
       frq_str_e = ""
     else
       frq_str_e = "-"//frqs_str(nfrq)
     end if
     freq_str = frq_str_s//frq_str_e
  end if
!      frq_str_list = frq_str_list(:len_trim(frq_str_list)) // "_" //  frqs_str(ff)

  if (verbose .gt. 1) print *,"input_file: ",input_file(:len_trim(input_file)),&
       " namelist file: ",namelist_file," freq: ",freq_str

!!! read n-moments file
  if (n_moments .eq. 2) call double_moments_module_read(moments_file) !from double_moments_module.f90

!!! read the data
  call vars_atmosphere_read_profile(input_file) !from vars_atmosphere.f90

  ! now allocate variables
  call allocate_vars


  ! This GCE model format does not have all the fields expected by    
  ! the radiative transfer code (i.e. total pressure, and water vapor 
  ! pressure for this model).  Assign/compute the missing fields first
  ! make layer averages
  call get_atmosG0

  if (write_nc .eqv. .false.) call mod_io_strings_get_filename()


  if (verbose .gt. 1) print*, 'Start loop over frequencies & profiles!'

  grid_f: do fi =1, nfrq
     grid_y: do ny = 1, ngridy !ny_in, ny_fin  
        grid_x: do nx = 1, ngridx !nx_in, nx_fin   
            
            model_i = profiles(nx,ny)%isamp
            model_j = profiles(nx,ny)%jsamp

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



           !run the model
           call run_rt3(nx,ny,fi,freqs(fi),frqs_str(fi),&
            model_i,model_j,ground_temp,lat,lon,lfrac,relhum_lev,press_lev,&
            temp_lev,hgt_lev,cwc_q,iwc_q,rwc_q,swc_q,gwc_q,hwc_q,cwc_n,&
            iwc_n,rwc_n,swc_n,gwc_n,hwc_n,press,temp,relhum,vapor_pressure,&
            rho_vap,q_hum)

        end do grid_x
     end do grid_y
  end do grid_f

  if (write_nc) then
     nc_out_file = trim(output_path)//"/"//trim(input_file(1:len_trim(input_file)-4))//&
          trim(freq_str)//trim(file_desc)//'.nc'
     call write_nc_results(nc_out_file)
  end if

end program pamtra
