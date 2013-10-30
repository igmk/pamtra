program pamtra

    use kinds, only: long
    !    use constants !physical constants live here
    use settings !all settings go here
    use vars_atmosphere !input variables and reading routine
    use vars_output !output variables
    use vars_profile
    use vars_jacobian, only: allocate_jacobian_vars, deallocate_jacobian_vars
    use double_moments_module !double moments variables are stored here
    use report_module
    use descriptor_file
    use deallocate_everything, only : do_deallocate_everything
    use vars_index, only: i_x, i_y, i_f

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

    real(kind=dbl) :: nan


    !!!loop variables
    integer(kind=long) :: i

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
    call settings_read(err)  !from settings.f90
    if (err /= 0) then
        msg = 'error in settings_read!'
        call report(err, msg, nameOfRoutine)
        errorstatus = err
        return
    end if

    in_python = .false.! we are _not_ in python

    if (verbose >= 1) then
        msg = "input_file: "//input_file(:len_trim(input_file))//&
        " namelist file: "//trim(namelist_file)//&
        " freqs: "//trim(frqs_str(1))//" to "//trim(frqs_str(nfrq))
        call report(info, msg, nameOfRoutine)
    end if

    !!! read n-moments file
!     if (n_moments == 2) call double_moments_module_read(moments_file) !from double_moments_module.f90


    call read_descriptor_file(err)
    if (err /= 0) then
        msg = 'Error in read_descriptor_file!'
        call report(fatal, msg, nameOfRoutine)
        errorstatus = err
        go to 666
    end if

!!! read the data
call get_atmosphere
! 1tmporary: this should go into the call get_atmosphere routine!


atmo_max_nlyrs =profiles_nlyr
atmo_ngridx = profiles_ngridx
atmo_ngridy = profiles_ngridy

call allocate_atmosphere_vars(err)
if (err /= 0) then
    msg = 'Error in allocate_atmosphere_vars!'
    call report(fatal, msg, nameOfRoutine)
  errorstatus = err
  return
end if

atmo_nlyrs(:,:) = profiles_nlyr

!temporary loop to fill atmosphere array:
do i_y = 1, profiles_ngridy !i_x_in, i_x_fin
  do i_x = 1, profiles_ngridx

      atmo_relhum_lev(i_x,i_y,:) = profiles(i_x,i_y)%relhum_lev
      atmo_press_lev(i_x,i_y,:) = profiles(i_x,i_y)%press_lev
      atmo_temp_lev(i_x,i_y,:) = profiles(i_x,i_y)%temp_lev
      atmo_hgt_lev(i_x,i_y,:) = profiles(i_x,i_y)%hgt_lev

      atmo_hydro_reff(i_x,i_y,:,:) = nan()
      atmo_hydro_n(i_x,i_y,:,:) = nan()


      atmo_hydro_q(i_x,i_y,:,1) = profiles(i_x,i_y)%cloud_water_q
!       atmo_hydro_n(i_x,i_y,:,1) = profiles(i_x,i_y)%cloud_water_n

      atmo_hydro_q(i_x,i_y,:,2) = profiles(i_x,i_y)%cloud_ice_q
!       atmo_hydro_n(i_x,i_y,:,2) = profiles(i_x,i_y)%cloud_ice_n

      atmo_hydro_q(i_x,i_y,:,3) = profiles(i_x,i_y)%rain_q
!       atmo_hydro_n(i_x,i_y,:,3) = profiles(i_x,i_y)%rain_n

      atmo_hydro_q(i_x,i_y,:,4) = profiles(i_x,i_y)%snow_q
!       atmo_hydro_n(i_x,i_y,:,4) = profiles(i_x,i_y)%snow_n

      atmo_hydro_q(i_x,i_y,:,5) = profiles(i_x,i_y)%graupel_q
!       atmo_hydro_n(i_x,i_y,:,5) = profiles(i_x,i_y)%graupel_n
    atmo_month(i_x,i_y) = profiles_month
    atmo_day(i_x,i_y) = profiles_day
    atmo_year(i_x,i_y) = profiles_year
    atmo_time(i_x,i_y) =profiles_time
    atmo_deltax(i_x,i_y) = profiles_deltax
    atmo_deltay(i_x,i_y) = profiles_deltay
    atmo_model_i(i_x,i_y) = profiles(i_x,i_y)%isamp
    atmo_model_j(i_x,i_y) = profiles(i_x,i_y)%jsamp
    atmo_lon(i_x,i_y) = profiles(i_x,i_y)%longitude       
    atmo_lat(i_x,i_y) = profiles(i_x,i_y)%latitude       
    atmo_lfrac(i_x,i_y) = profiles(i_x,i_y)%land_fraction
    atmo_wind10u(i_x,i_y) = profiles(i_x,i_y)%wind_10u
    atmo_wind10v(i_x,i_y) = profiles(i_x,i_y)%wind_10v

    atmo_iwv(i_x,i_y) = profiles(i_x,i_y)%iwv

    end do
end do


    ! make sure that all the levels and layer variables are present
    call fillMissing_atmosphere_vars(err)
    if (err /= 0) then
        msg = 'Error in fillMissing_atmosphere_vars!'
        call report(fatal, msg, nameOfRoutine)
      errorstatus = err
      return
    end if


    ! now allocate variables
    call allocate_output_vars(err, atmo_max_nlyrs)
    if (err /= 0) then
        msg = 'Error in allocate_output_vars!'
        call report(fatal, msg, nameOfRoutine)
        errorstatus = err
        go to 666
    end if
    msg = 'Start loop over frequencies & profiles!'
    if (verbose >= 2)  call report(info, msg, nameOfRoutine)



    grid_f: do i_f =1, nfrq
        if (jacobian_mode) then
            !for jacobian mode. non disturbed profile is expected in grid 1,1!
            call allocate_jacobian_vars(atmo_nlyrs(i_x,i_y))
        end if
          grid_y: do i_y = 1, atmo_ngridy !i_x_in, i_x_fin
              grid_x: do i_x = 1, atmo_ngridx !i_y_in, i_y_fin

                !run the model
                call run_rt(err)
                if (err /= 0) then
                    msg = 'Error in run_rt!'
                    call report(fatal, msg, nameOfRoutine)
                    errorstatus = err
                    go to 666
                end if

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
