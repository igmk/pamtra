subroutine run_rt(errorstatus)

    use kinds, only: long, dbl
    use constants, only: c,&
    pi,&
    sky_temp

    use settings !all settings go here
    use vars_atmosphere !input variables and reading routine
    use vars_output !output variables
    use vars_rt, only: rt_kextatmo, allocate_rt_vars, deallocate_rt_vars
    use double_moments_module
    use mod_io_strings, only: xstr, nxstr, ystr, nystr, frq_str
    use report_module
    use vars_index, only: i_x, i_y, i_f

    implicit none

    real(kind=dbl) :: freq ! frequency [GHz]

    integer(kind=long), dimension(maxlay) :: OUTLEVELS
    integer(kind=long) :: nz

    real(kind=dbl), dimension(maxv) :: MU_VALUES
    real(kind=dbl) :: wavelength       ! microns
    real(kind=dbl) :: ground_albedo

    complex(kind=dbl) :: ground_index

    character(300) :: OUT_FILE_PAS, OUT_FILE_ACT !file names if no nc

    ! Error handling

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'run_rt'

    interface
      subroutine RT4(errorstatus,mu_values,out_file,&
      ground_type,ground_albedo,ground_index,sky_temp,&
      wavelength,outlevels)
        use kinds
        implicit none
    integer   maxv, maxlay
    parameter (maxv=64)
    parameter (maxlay=200)
        real(kind=dbl), intent(in) ::  mu_values(maxv)
        character*64, intent(in) :: out_file
        character, intent(in) ::  ground_type*1
        real(kind=dbl), intent(in) ::  ground_albedo
        complex*16, intent(in) ::   ground_index
        real(kind=dbl), intent(in) ::  sky_temp
        real(kind=dbl), intent(in) ::   wavelength
        integer, intent(in) ::  outlevels(maxlay)
        integer(kind=long), intent(out) :: errorstatus
      end subroutine RT4

    end interface


    if (verbose >= 1) call report(info,'Start of ', nameOfRoutine)

    call allocate_rt_vars(err)
    if (err /= 0) then
        msg = 'Error in allocate_rt_vars!'
        call report(fatal, msg, nameOfRoutine)
      errorstatus = err
      return
    end if

    freq = freqs(i_f)
    frq_str = frqs_str(i_f)
    wavelength = c / (freq*1.d3)   ! microns

    !  if (verbose .gt. 0) print*, "calculating: ", frq_str, " Y:",i_y, " of ", ngridy, "X:", i_x, " of ", ngridx

    write(xstr, '(i3.3)') atmo_model_i(i_x,i_y)
    write(ystr, '(i3.3)') atmo_model_j(i_x,i_y)
    write(nxstr, '(i4)') atmo_ngridx
    write(nystr, '(i4)') atmo_ngridy

    msg = "calculating: "// frq_str// " Y: "//ystr//" of "//nystr//" X: "//xstr//" of "//nxstr

    if (verbose >= 2) call report(info,msg, nameOfRoutine)
    ! This GCE model format does not have all the fields expected by
    ! the radiative transfer code (i.e. total pressure, and water vapor
    ! pressure for this model).  Assign/compute the missing fields first
    ! make layer averages

    if (verbose >= 2) call report(info,nxstr//' '//nystr//'type to local variables done',nameOfRoutine)

    call get_surface(err,freq, atmo_groundtemp(i_x,i_y), salinity, ground_albedo,ground_index,ground_type)
    if (err /= 0) then
        msg = 'error in get_surface'
        call report(err,msg, nameOfRoutine)
        errorstatus = err
        return
    end if

    ! Determine surface properties
    if (verbose >= 2) call report(info,'surface emissivity calculated:', nameOfRoutine)

    ! gaseous absorption
    !
    ! rt_kextatmo   extinction by moist air [Np/m]
    !
    if (lgas_extinction) then
        !returns rt_kextatmo!
        call get_gasabs(err,freq)
        if (err /= 0) then
            msg = 'error in get_gasabs'
            call report(err,msg, nameOfRoutine)
            errorstatus = err
            return
        end if
    else
        rt_kextatmo = 0._dbl ! for the whole column
    end if
    !save atmospheric attenuation and height for radar
    if (active) then
        Att_atmo(i_x,i_y,:,i_f)  = 10._dbl*log10(exp(rt_kextatmo*atmo_delta_hgt_lev(i_x,i_y,:)))
        radar_hgt(i_x,i_y,:) = atmo_hgt(i_x,i_y,:)
    end if


    if (verbose >= 2) print*, i_x,i_y, 'Gas absorption calculated'

    ! hydrometeor extinction desired
    if (lhyd_extinction) then
    
!       call hydrometeor_extinction_rt4(err,freq,i_x,i_y,i_f)!hier i_x, i_y
            call hydrometeor_extinction(err)!hier i_x, i_y

    end if    
    
      if (err == 2) then
	msg = 'Error in run_drop_size_dist'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	return
      end if

    !
    if (dump_to_file) call dump_profile()

    !&&&&&&&&   I/O FILE NAMES   &&&&&&&&&&&&&&&&&&

    OUT_FILE_PAS = output_path(:len_trim(output_path))//"/"//&
      atmo_year(i_x, i_y)//atmo_month(i_x, i_y)//atmo_day(i_x, i_y)//atmo_time(i_x, i_y)//&
      'x'//xstr//'y'//ystr//'f'//trim(frq_str)//"_passive"

    OUT_FILE_ACT = output_path(:len_trim(output_path))//"/"//&
      atmo_year(i_x, i_y)//atmo_month(i_x, i_y)//atmo_day(i_x, i_y)//atmo_time(i_x, i_y)//&
      'x'//xstr//'y'//ystr//'f'//trim(frq_str)//"_active"


    !save active to ASCII
    if (active .and. (write_nc .eqv. .false.) .and. (in_python .eqv. .false.)) then
        call save_active(OUT_FILE_ACT,i_x,i_y,i_f)
    end if

    ! find the output level
print*, "EMILIANO, double check here height index of atom_hgt_lev, please (run_rt.f90)"
    if (atmo_obs_height(i_x,i_y) > 99999._dbl .or. atmo_obs_height(i_x,i_y) > atmo_hgt_lev(i_x,i_y,atmo_nlyrs(i_x,i_y)+1)) then
        outlevels(1) = 1
    else if (atmo_obs_height(i_x,i_y) .lt. 0.1 .or. atmo_obs_height(i_x,i_y) .lt. atmo_hgt_lev(i_x,i_y,2)) then
        outlevels(1) = atmo_nlyrs(i_x,i_y) + 1
    else
        out_search: do nz = 1, atmo_nlyrs(i_x,i_y)
            if (atmo_hgt_lev(i_x,i_y,nz+1) .ge. atmo_obs_height(i_x,i_y)) then
                if (abs(atmo_hgt_lev(i_x,i_y,nz+1) - atmo_obs_height(i_x,i_y)) .lt. &
                    abs(atmo_hgt_lev(i_x,i_y,nz) - atmo_obs_height(i_x,i_y))) then
                    outlevels(1) = atmo_nlyrs(i_x,i_y)-nz+1
                else
                    outlevels(1) = atmo_nlyrs(i_x,i_y)-nz+2
                end if
                exit out_search
            end if
        end do out_search
    end if

    OUTLEVELS(2) = atmo_nlyrs(i_x,i_y)+1    ! this is the bottom

    if (passive .eqv. .true.) then

        if (verbose >= 2) print*, i_x,i_y, "Entering rt4 ...."

        call rt4(err, mu_values,out_file_pas,&
        ground_type,ground_albedo,ground_index,sky_temp,&
        wavelength,outlevels)

        if (verbose >= 2) print*, i_x,i_y, "....rt4 finished"
        !calculate human readable angles!
        angles_deg(1:NUMMU) = 180-(180.*acos(MU_VALUES(NUMMU:1:-1))/pi)
        angles_deg(1+NUMMU:2*NUMMU) = (180.*acos(MU_VALUES(1:NUMMU))/pi)

    end if

    !DEALLOCATE rt variables
    call deallocate_rt_vars()

    if (verbose >= 1) call report(info,'End of ', nameOfRoutine)
    errorstatus = err
    return

end subroutine run_rt
