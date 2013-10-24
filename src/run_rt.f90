subroutine run_rt(errorstatus)

    use kinds, only: long, dbl
    use constants, only: c,&
    pi,&
    sky_temp

    use settings !all settings go here
    use vars_atmosphere !input variables and reading routine
    use vars_output !output variables
    use vars_rt, only: rt_kextatmo
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
    real(kind=dbl) :: GROUND_TEMP, ground_albedo

    complex(kind=dbl) :: ground_index

    character(300) :: OUT_FILE_PAS, OUT_FILE_ACT !file names if no nc

    ! Error handling

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'run_rt'

    if (verbose >= 1) call report(info,'Start of ', nameOfRoutine)

    freq = freqs(i_f)
    frq_str = frqs_str(i_f)
    wavelength = c / (freq*1.d3)   ! microns
    GROUND_TEMP = atmo_temp_lev(i_x,i_y,1)
print*, "replace by atmo_groundtemp?"
    !  if (verbose .gt. 0) print*, "calculating: ", frq_str, " Y:",i_y, " of ", ngridy, "X:", i_x, " of ", ngridx

    write(xstr, '(i3.3)') model_i
    write(ystr, '(i3.3)') model_j
    write(nxstr, '(i4)') ngridx
    write(nystr, '(i4)') ngridy

    msg = "calculating: "// frq_str// " Y: "//ystr//" of "//nystr//" X: "//xstr//" of "//nxstr

    if (verbose >= 2) call report(info,msg, nameOfRoutine)
    ! This GCE model format does not have all the fields expected by
    ! the radiative transfer code (i.e. total pressure, and water vapor
    ! pressure for this model).  Assign/compute the missing fields first
    ! make layer averages

    call get_atmosG0()

    if (verbose >= 2) call report(info,nxstr//' '//nystr//'type to local variables done',nameOfRoutine)

    call get_surface(err,freq, ground_temp, salinity, ground_albedo,ground_index,ground_type)
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
    date_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_passive"

    OUT_FILE_ACT = output_path(:len_trim(output_path))//"/"//&
    date_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_active"

    !save active to ASCII
    if (active .and. (write_nc .eqv. .false.) .and. (in_python .eqv. .false.)) then
        call save_active(OUT_FILE_ACT,i_x,i_y,i_f)
    end if


    if (write_nc) then
        !      Output integrated quantities
        call collect_boundary_output(lon,lat,lfrac,&
        iwv, cwp,iwp,rwp,swp, &
        gwp,hwp,model_i,model_j,i_x,i_y)
        if (verbose >= 2) print*, i_x,i_y, 'collect_boundary_output done'
    end if

    ! find the output level
print*, "EMILIANO, double check here height index of atom_hgt_lev, please (run_rt.f90)"
    if (obs_height > 99999._dbl .or. obs_height > atmo_hgt_lev(i_x,i_y,atmo_nlyrs(i_x,i_y)+1)) then
        outlevels(1) = 1
    else if (obs_height .lt. 0.1 .or. obs_height .lt. atmo_hgt_lev(i_x,i_y,2)) then
        outlevels(1) = atmo_nlyrs(i_x,i_y) + 1
    else
        out_search: do nz = 1, atmo_nlyrs(i_x,i_y)
            if (atmo_hgt_lev(i_x,i_y,nz+1) .ge. obs_height) then
                if (abs(atmo_hgt_lev(i_x,i_y,nz+1) - obs_height) .lt. abs(atmo_hgt_lev(i_x,i_y,nz) - obs_height)) then
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

        call rt4(err, nstokes,nummu,mu_values,out_file_pas,quad_type,ground_temp,&
        ground_type,ground_albedo,ground_index,sky_temp,&
        wavelength,units,outpol,noutlevels,outlevels)

        if (verbose >= 2) print*, i_x,i_y, "....rt4 finished"
        !calculate human readable angles!
        angles_deg(1:NUMMU) = 180-(180.*acos(MU_VALUES(NUMMU:1:-1))/pi)
        angles_deg(1+NUMMU:2*NUMMU) = (180.*acos(MU_VALUES(1:NUMMU))/pi)

    end if

    if (verbose >= 1) call report(info,'End of ', nameOfRoutine)
    errorstatus = err
    return

end subroutine run_rt
