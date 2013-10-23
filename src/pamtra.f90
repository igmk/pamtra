program pamtra

    use kinds, only: long
    !    use constants !physical constants live here
    use settings !all settings go here
    use vars_atmosphere !input variables and reading routine
    use vars_output !output variables
    use vars_profile
    use double_moments_module !double moments variables are stored here
    use report_module
    use descriptor_file
    use deallocate_everything, only : do_deallocate_everything

    !     The code reads a full (e.g. COSMO) grid and computes for each
    !     profile the radiative transfer for the given frequencies
    !
    !     By convention, the quantities followed  by "_lev"
    !     are given at the layer heights while the quantitites w/o
    !     "_lev" are layer average quantities
    !

    implicit none

    !!! internal "handle command line parameters" !!!

    integer(kind=long) :: inarg
    character(40) :: gitHash, gitVersion

    !!! set by "handle command line parameters" !!!


    !!!loop variables
    integer(kind=long) ::  fi,nx, ny,i

    ! Error handling

    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=200) :: msg
    character(len=14) :: nameOfRoutine = 'pamtra'

    !get git data
    call versionNumber(gitVersion,gitHash)

    !get and process command line parameters
    call parse_options(gitVersion,gitHash)

    !!! read variables from namelist file
    call settings_read  !from settings.f90

    in_python = .false.! we are _not_ in python

    if (verbose >= 1) then
        msg = "input_file: "//input_file(:len_trim(input_file))//&
        " namelist file: "//trim(namelist_file)//&
        " freqs: "//trim(frqs_str(1))//" to "//trim(frqs_str(nfrq))
        call report(info, msg, nameOfRoutine)
    end if

    !!! read n-moments file
    if (n_moments == 2) call double_moments_module_read(moments_file) !from double_moments_module.f90

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


    call read_descriptor_file(err)
    if (err /= 0) then
        msg = 'Error in read_descriptor_file!'
        call report(fatal, msg, nameOfRoutine)
        errorstatus = err
        go to 666
    end if
    ! now allocate variables
    call allocate_output_vars(nlyr)

    msg = 'Start loop over frequencies & profiles!'
    if (verbose >= 2)  call report(info, msg, nameOfRoutine)



    grid_f: do fi =1, nfrq
        if (jacobian_mode) then
            !for jacobian mode. non disturbed profile is expected in grid 1,1!
            call allocate_jacobian_vars
        end if
        grid_y: do ny = 1, ngridy !nx_in, nx_fin
            grid_x: do nx = 1, ngridx !ny_in, ny_fin
         
                call allocate_profile_vars(err)
                if (err /= 0) then
                    msg = 'Error in allocate_profile_vars!'
                    call report(fatal, msg, nameOfRoutine)
                  errorstatus = err
                  go to 666
                end if
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
                call run_rt(err,nx,ny,fi)
                if (err /= 0) then
                    msg = 'Error in run_rt!'
                    call report(fatal, msg, nameOfRoutine)
                    errorstatus = err
                    go to 666
                end if
                !DEALLOCATE profile variables
                call deallocate_profile_vars()
            end do grid_x
        end do grid_y
        if (jacobian_mode) then
                  !for jacobian mode
            call deallocate_jacobian_vars
        end if
    end do grid_f

    
    if (write_nc) then
      call write_nc_results
    end if

    !now clean up and deallocate ALL variables

666 call do_deallocate_everything(err)
    if (err /= 0) then
        msg = 'Error in do_deallocate_everything!'
        call report(fatal, msg, nameOfRoutine)
        errorstatus = err
    end if


    if (verbose >= 1 .and. errorstatus == 0) then
        msg = 'Progam finished successfully'
        call report(info, msg, nameOfRoutine)
    end if
    if (errorstatus /= 0) then
        msg = 'Something went wrong!'
        call report(fatal, msg, nameOfRoutine)
    end if

end program pamtra
