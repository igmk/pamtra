  module pyPamtraLib
  
    use kinds, only: long
    !    use constants !physical constants live here
    use settings !all settings go here
    use vars_atmosphere !input variables 
    use vars_output !output variables
    use report_module
    use descriptor_file
    use vars_index, only: i_x, i_y, i_z, i_f, i_h

    implicit none
    character(40) :: gitHash, gitVersion

    contains


  subroutine run_pamtra(errorstatus)

    use mod_io_strings, only: formatted_frqstr


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
      real(kind=dbl) :: nan
  

      !!!loop variables
      integer(kind=long) ::  fi,i

      ! Error handling

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=200) :: msg
      character(len=14) :: nameOfRoutine = 'run_pamtra'
      character(len=8) :: frq_str

     if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
     
      errorstatus = 0
      
      !get git data
      call versionNumber(gitVersion,gitHash)

      !get and process command line parameters
!       call parse_options(gitVersion,gitHash)


!fill frqs_str array!
    do fi = 1, nfrq
        write (frq_str, '(f8.3)') freqs(fi)
    frqs_str(fi) = repeat('0', 8-len_trim(adjustl(frq_str)))//adjustl(frq_str)
    end do


      in_python = .true.! we are _in_ python


    !!! add missing settings and test them
    call add_settings(err)  !from settings.f90
    if (err /= 0) then
        msg = 'error in add_settings!'
        call report(err, msg, nameOfRoutine)
        errorstatus = err
        return
    end if
    call test_settings(err)  !from settings.f90
    if (err /= 0) then
        msg = 'error in test_settings!'
        call report(err, msg, nameOfRoutine)
        errorstatus = err
        return
    end if

      if (verbose >= 1) then
          msg = "freqs: "//trim(frqs_str(1))//" to "//trim(frqs_str(nfrq))
          call report(info, msg, nameOfRoutine)
      end if


      if (add_obs_height_to_layer) then
          msg = 'add_obs_height_to_layer not implemented yet in python'
          call report(fatal, msg, nameOfRoutine)
          errorstatus = err
          return
        end if
!         call add_obs_height(errorstatus)
!         if (err /= 0) then
!             msg = 'Error in add_obs_height!'
!             call report(fatal, msg, nameOfRoutine)
!           errorstatus = err
!           return
!         end if


      ! now allocate variables
      call allocate_output_vars(err,atmo_max_nlyrs)
      if (err /= 0) then
          msg = 'Error in allocate_output_vars!'
          call report(fatal, msg, nameOfRoutine)
        errorstatus = err
        return
      end if


      msg = 'Start loop over frequencies & profiles!'
      if (verbose >= 2)  call report(info, msg, nameOfRoutine)



      grid_f: do i_f =1, nfrq
          grid_y: do i_y = 1, atmo_ngridy !i_x_in, i_x_fin
              grid_x: do i_x = 1, atmo_ngridx !i_y_in, i_y_fin
                  !run the model
                  call run_rt(err)
                  if (err /= 0) then
                      msg = 'Error in run_rt!'
                      call report(fatal, msg, nameOfRoutine)
                      errorstatus = err
                      return
                  end if

              end do grid_x
          end do grid_y
      end do grid_f

      if ((verbose >= 1) .and. (errorstatus == 0)) then
          msg = 'Progam finished successfully'
          call report(info, msg, nameOfRoutine)
      end if
      if (errorstatus /= 0) then
          msg = 'Something went wrong!'
          call report(fatal, msg, nameOfRoutine)
      end if

     if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
      
      
  end subroutine run_pamtra

end module pyPamtraLib