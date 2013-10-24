  module pyPamtraLib
  
    use kinds, only: long
    !    use constants !physical constants live here
    use settings !all settings go here
    use vars_atmosphere !input variables and reading routine
    use vars_rt, only: allocate_rt_vars, deallocate_rt_vars
    use vars_output !output variables
    use vars_profile
    use vars_jacobian, only: allocate_jacobian_vars, deallocate_jacobian_vars
    use double_moments_module !double moments variables are stored here
    use report_module
    use descriptor_file
    use vars_index, only: i_x, i_y, i_f

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

     if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
     
      errorstatus = 0
      
      !get git data
      call versionNumber(gitVersion,gitHash)

      !get and process command line parameters
!       call parse_options(gitVersion,gitHash)


!fill frqs_str array!
    do fi = 1, nfrq
        write(frqs_str(fi),"(F8.2)") freqs(fi)
    end do


      !!! read variables from namelist file
!       call settings_read  !from settings.f90

      in_python = .false.! we are _not_ in python

      if (verbose >= 1) then
          msg = "input_file: "//input_file(:len_trim(input_file))//&
          " namelist file: "//trim(namelist_file)//&
          " freqs: "//trim(frqs_str(1))//" to "//trim(frqs_str(nfrq))
          call report(info, msg, nameOfRoutine)
      end if





!!! read the data
call get_atmosphere
! 1tmporary: this should go into the call get_atmosphere routine!

year = profiles_year
month = profiles_month
day = profiles_day
time = profiles_time
ngridx = profiles_ngridx
ngridy = profiles_ngridy
deltax = profiles_deltax
deltay = profiles_deltay
date_str = year//month//day//time


!python input can have different number of nlayer
atmo_max_nlyr = profiles_nlyr

call allocate_atmosphere_vars()

atmo_nlyrs(:,:) = profiles_nlyr

!temporary loop to fill atmosphere array:
do i_y = 1, ngridy !i_x_in, i_x_fin
  do i_x = 1, ngridx
print*, i_y, i_x
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
    end do
end do
      ! now allocate variables
      call allocate_output_vars(atmo_max_nlyr)



      msg = 'Start loop over frequencies & profiles!'
      if (verbose >= 2)  call report(info, msg, nameOfRoutine)



      grid_f: do i_f =1, nfrq
          if (jacobian_mode) then
              !for jacobian mode. non disturbed profile is expected in grid 1,1!
            call allocate_jacobian_vars(atmo_nlyrs(i_x,i_y))
          end if
          grid_y: do i_y = 1, ngridy !i_x_in, i_x_fin
              grid_x: do i_x = 1, ngridx !i_y_in, i_y_fin

                  call allocate_rt_vars(err)
                  if (err /= 0) then
                      msg = 'Error in allocate_rt_vars!'
                      call report(fatal, msg, nameOfRoutine)
                    errorstatus = err
                    return
                  end if
                  !   ground_temp = profiles(i_x,i_y)%temp_lev(0)       ! K
                  lat = profiles(i_x,i_y)%latitude                  ! °
                  lon = profiles(i_x,i_y)%longitude                 ! °
                  lfrac = profiles(i_x,i_y)%land_fraction
!                   relhum_lev = profiles(i_x,i_y)%relhum_lev         ! %
!                   press_lev = profiles(i_x,i_y)%press_lev           ! Pa
!                   temp_lev = profiles(i_x,i_y)%temp_lev             ! K
!                   hgt_lev = profiles(i_x,i_y)%hgt_lev               ! m

                  model_i = profiles(i_x,i_y)%isamp
                  model_j = profiles(i_x,i_y)%jsamp
                  wind10u = profiles(i_x,i_y)%wind_10u
                  wind10v = profiles(i_x,i_y)%wind_10v

                  iwv = profiles(i_x,i_y)%iwv
                  cwp = profiles(i_x,i_y)%cwp
                  iwp = profiles(i_x,i_y)%iwp
                  rwp = profiles(i_x,i_y)%rwp
                  swp = profiles(i_x,i_y)%swp
                  gwp = profiles(i_x,i_y)%gwp
                  hwp = profiles(i_x,i_y)%hwp

!                   cwc_q = profiles(i_x,i_y)%cloud_water_q           ! kg/kg
!                   iwc_q = profiles(i_x,i_y)%cloud_ice_q             ! kg/kg
!                   rwc_q = profiles(i_x,i_y)%rain_q                  ! kg/kg
!                   swc_q = profiles(i_x,i_y)%snow_q                  ! kg/kg
!                   gwc_q = profiles(i_x,i_y)%graupel_q               ! kg/kg
!                   if (n_moments .eq. 2) then
!                       hwc_q = profiles(i_x,i_y)%hail_q              ! kg/kg
!                       cwc_n = profiles(i_x,i_y)%cloud_water_n       ! #/kg
!                       iwc_n = profiles(i_x,i_y)%cloud_ice_n         ! #/kg
!                       rwc_n = profiles(i_x,i_y)%rain_n              ! #/kg
!                       swc_n = profiles(i_x,i_y)%snow_n              ! #/kg
!                       gwc_n = profiles(i_x,i_y)%graupel_n           ! #/kg
!                       hwc_n = profiles(i_x,i_y)%hail_n              ! #/kg
!                   end if

                  !run the model
                  call run_rt(err)
                  if (err /= 0) then
                      msg = 'Error in run_rt!'
                      call report(fatal, msg, nameOfRoutine)
                      errorstatus = err
                      return
                  end if
                  !DEALLOCATE rt variables
                  call deallocate_rt_vars()
              end do grid_x
          end do grid_y
          if (jacobian_mode) then
                    !for jacobian mode
              call deallocate_jacobian_vars
          end if
      end do grid_f

      if (verbose >= 1 .and. errorstatus == 0) then
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