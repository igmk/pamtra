program pamtra

    use kinds
    use constants !physical constants live here
    use nml_params !all settings go here
    use file_mod
    use vars_atmosphere !input variables and reading routine
    use vars_output !output variables
    use vars_profile
    use double_moments_module !double moments variables are stored here


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
    character(8) :: formatted_frqstr !function call

    !!! set by "handle command line parameters" !!!

    character(8), dimension(maxfreq) :: frqs_str !from commandline
    character(9) :: frq_str_s,frq_str_e

    !!!loop variables
    integer ::  fi,nx, ny,i

    !get git data
    call versionNumber(gitVersion,gitHash)

    !get and process command line parameters
    call parse_options(gitVersion,gitHash,frqs_str,nfrq)

    !  inarg = iargc()
    !
    !  if (inarg .lt. 3) then
    !     print *,'Usage: pamtra profile_file namelist_file (list of frequencies)'
    !     print *,'Example: ./pamtra rt_comp_single.dat run_params.nml 35 94'
    !     print *,'See namelist file for further pamtra options'
    !     print *,''
    !     print *,'Version:  '//gitVersion
    !     print *,'Git Hash: '//gitHash
    !     stop
    !  else if (inarg .gt. maxfreq + 2) then
    !     print *,'Too many frequencies! Increase maxfreq!'
    !     stop
    !  end if
    !  call getarg(1,input_file)
    !  call getarg(2,namelist_file)
    !
    !  nfrq = inarg - 2
    !
    do ff = 1, nfrq
        !     call getarg(ff+2,frqs_str(ff))
        read(frqs_str(ff),*) freqs(ff)
        frqs_str(ff) = formatted_frqstr(frqs_str(ff))
    end do

    !!! read variables from namelist file
    call nml_params_read !from nml_params.f90

    ! create frequency string if not set in pamtra
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
    call get_atmosphere

    year = profiles_year
    month = profiles_month
    day = profiles_day
    time = profiles_time
    ngridx = profiles_ngridx
    ngridy = profiles_ngridy
    nlyr = profiles_nlyr
    deltax = profiles_deltax
    deltay = profiles_deltay
    date_str = year//month//day//time

    ! now allocate variables
    call allocate_output_vars(nlyr)

    if (verbose .gt. 1) print*, 'Start loop over frequencies & profiles!'

    grid_f: do fi =1, nfrq
        grid_y: do ny = 1, ngridy !ny_in, ny_fin
            grid_x: do nx = 1, ngridx !nx_in, nx_fin
         
                call allocate_profile_vars
         
                !   ground_temp = profiles(nx,ny)%temp_lev(0)       ! K
                lat = profiles(nx,ny)%latitude                  ! °
                lon = profiles(nx,ny)%longitude                 ! °
                lfrac = profiles(nx,ny)%land_fraction
                relhum_lev = profiles(nx,ny)%relhum_lev         ! %
                press_lev = profiles(nx,ny)%press_lev           ! Pa
                temp_lev = profiles(nx,ny)%temp_lev             ! K
                hgt_lev = profiles(nx,ny)%hgt_lev               ! m

                model_i = profiles(nx,ny)%isamp
                model_j = profiles(nx,ny)%jsamp
                wind10u = profiles(nx,ny)%wind_10u
                wind10v = profiles(nx,ny)%wind_10v

                iwv = profiles(nx,ny)%iwv
                cwp = profiles(nx,ny)%cwp
                iwp = profiles(nx,ny)%iwp
                rwp = profiles(nx,ny)%rwp
                swp = profiles(nx,ny)%swp
                gwp = profiles(nx,ny)%gwp
                hwp = profiles(nx,ny)%hwp


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

                !run the model
                call run_rt(nx,ny,fi,freqs(fi),frqs_str(fi))

                call deallocate_profile_vars()

            end do grid_x
        end do grid_y
    end do grid_f

    if (write_nc) then
        call write_nc_results
    end if

    call deallocate_output_vars()

end program pamtra
