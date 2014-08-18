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

    errorstatus = success
    
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

    call screen_input(err)
    if (err /= 0) then
        msg = 'Error in screen_input!'
        call report(fatal, msg, nameOfRoutine)
      errorstatus = err
      go to 666
    end if

    call allocate_atmosphere_vars(err)
    if (err /= 0) then
        msg = 'Error in allocate_atmosphere_vars!'
        call report(fatal, msg, nameOfRoutine)
      errorstatus = err
      go to 666
    end if

    if ((atmo_input_type == 'lev') .or. (atmo_input_type == 'lay')) &
        call read_new_fill_variables(err)
    if (atmo_input_type == 'cla') &
        call read_classic_fill_variables(err)

    ! make sure that all the levels and layer variables are present
    call fillMissing_atmosphere_vars(err)
    if (err /= 0) then
        msg = 'Error in fillMissing_atmosphere_vars!'
        call report(fatal, msg, nameOfRoutine)
      errorstatus = err
      go to 666
    end if

    if (add_obs_height_to_layer) then
      call add_obs_height(err)
      if (err /= 0) then
          msg = 'Error in add_obs_height!'
          call report(fatal, msg, nameOfRoutine)
        errorstatus = err
        go to 666
      end if
    end if

      ! now allocate output variables
      call allocate_output_vars(err,atmo_max_nlyrs)
      if (err /= 0) then
          msg = 'Error in allocate_output_vars!'
          call report(fatal, msg, nameOfRoutine)
        errorstatus = err
      go to 666
      end if


! do i_x = 1, atmo_ngridx
!   do i_y = 1, atmo_ngridy
! print*,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    LAYER"
!     call print_out_layer(i_x,i_y)
! print*,'-----------  ',atmo_nlyrs(i_x,i_y)
! print*,"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@    LEVEL"
!     call print_out_level(i_x,i_y)
!   enddo
! enddo

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
!       call collect_boundary_output
      call write_nc_results
    end if

    !now clean up and deallocate ALL variables

666 call do_deallocate_everything(err)
    if (err /= 0) then
        msg = 'Error in do_deallocate_everything!'
        call report(fatal, msg, nameOfRoutine)
        errorstatus = err
    end if

    if ((verbose >= 1) .and. (errorstatus == 0)) then
        msg = 'Progam finished successfully'
        call report(info, msg, nameOfRoutine)
    end if
    if (errorstatus /= 0) then
        msg = 'Something went wrong!'
        call report(fatal, msg, nameOfRoutine)
    end if

end program pamtra
