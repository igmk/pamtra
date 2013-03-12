subroutine run_rt(errorstatus, nx,ny,fi)

    use kinds, only: long, dbl
    use constants, only: c,&
    pi,&
    sky_temp

    use settings !all settings go here
    use vars_atmosphere !input variables and reading routine
    use vars_output !output variables
    use double_moments_module
    use mod_io_strings, only: xstr, nxstr, ystr, nystr, frq_str
    use report_module

    implicit none

    integer(kind=long), intent(in) :: nx,ny,fi
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

    freq = freqs(fi)
    frq_str = frqs_str(fi)
    wavelength = c / (freq*1.d3)   ! microns
    GROUND_TEMP = temp_lev(0)

    !  if (verbose .gt. 0) print*, "calculating: ", frq_str, " Y:",ny, " of ", ngridy, "X:", nx, " of ", ngridx

    write(xstr, '(i3.3)') model_i
    write(ystr, '(i3.3)') model_j
    write(nxstr, '(i4)') ngridx
    write(nystr, '(i4)') ngridy

    msg = "calculating: "// frq_str// " Y: "//ystr//" of "//nystr//" X: "//xstr//" of "//nxstr

    if (verbose >= 1) call report(info,msg, nameOfRoutine)
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
    ! kextatmo   extinction by moist air [Np/m]
    !
    if (lgas_extinction) then
        !returns kextatmo!
        call get_gasabs(err,freq)
        if (err /= 0) then
            msg = 'error in get_gasabs'
            call report(err,msg, nameOfRoutine)
            errorstatus = err
            return
        end if
    else
        kextatmo = 0._dbl ! for the whole column
    end if
    !save atmospheric attenuation and height for radar
    if (active) then
        Att_atmo(nx,ny,:,fi)  = 10._dbl*log10(exp(kextatmo*delta_hgt_lev))
        radar_hgt(nx,ny,:) = hgt(:)
    end if


    if (verbose >= 2) print*, nx,ny, 'Gas absorption calculated'


    ! hydrometeor extinction desired
    if (lhyd_extinction) then
        if (rt_mode .eq. 'rt3') then
            call hydrometeor_extinction_rt3(freq)
        elseif (rt_mode .eq. 'rt4') then
            call hydrometeor_extinction_rt4(freq,nx,ny,fi)!hier nx, ny
        end if
    end if



    !
    if (dump_to_file) call dump_profile()

    !&&&&&&&&   I/O FILE NAMES   &&&&&&&&&&&&&&&&&&

    OUT_FILE_PAS = output_path(:len_trim(output_path))//"/"//&
    date_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_passive"

    OUT_FILE_ACT = output_path(:len_trim(output_path))//"/"//&
    date_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_active"
    !
    !   if ((active) .and. ((radar_mode == "simple") .or. (radar_mode == "splitted")))  then
    !      call calculate_active(OUT_FILE_ACT,freq,&
    !           Ze(nx,ny,:,fi),Ze_cw(nx,ny,:,fi),Ze_rr(nx,ny,:,fi),Ze_ci(nx,ny,:,fi),&
    !           Ze_sn(nx,ny,:,fi),Ze_gr(nx,ny,:,fi),Ze_ha(nx,ny,:,fi),&
    !           Att_atmo(nx,ny,:,fi),Att_hydro(nx,ny,:,fi),Att_cw(nx,ny,:,fi),Att_rr(nx,ny,:,fi),&
    !           Att_ci(nx,ny,:,fi),Att_sn(nx,ny,:,fi),Att_gr(nx,ny,:,fi),Att_ha(nx,ny,:,fi))
    !      if (verbose .gt. 1) print*, nx,ny, 'calculate_active done'
    !
    !   end if

    !save active to ASCII
    if (active .and. (write_nc .eqv. .false.) .and. (in_python .eqv. .false.)) then
        call save_active(OUT_FILE_ACT,nx,ny,fi)
    end if


    if (write_nc) then
        !      Output integrated quantities
        call collect_boundary_output(lon,lat,lfrac,&
        iwv, cwp,iwp,rwp,swp, &
        gwp,hwp,model_i,model_j,nx,ny)
        if (verbose >= 2) print*, nx,ny, 'collect_boundary_output done'
    end if

    ! find the output level
    ! in rt3 and rt4 layers are reversed

    if (obs_height > 99999._dbl .or. obs_height > hgt_lev(nlyr)) then
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

    if (passive .eqv. .true.) then

        if (rt_mode .eq. 'rt3') then
            if (verbose >= 2) print*, nx,ny, "Entering rt3 ...."

            call RT3(NSTOKES, NUMMU, AZIORDER, MU_VALUES, src_code, &
            out_file_pas, QUAD_TYPE, deltam, DIRECT_FLUX,     &
            DIRECT_MU, GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO,  &
            GROUND_INDEX, SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,  &
            NOUTLEVELS, OUTLEVELS, nx,ny,fi)

            if (verbose >= 2) print*, nx,ny, "....rt3 finished"
        elseif (rt_mode .eq. 'rt4') then
            if (verbose >= 2) print*, nx,ny, "Entering rt4 ...."

            call rt4(nstokes,nummu,mu_values,out_file_pas,quad_type,ground_temp,&
            ground_type,ground_albedo,ground_index,sky_temp,&
            wavelength,units,outpol,noutlevels,outlevels,nx,ny,fi)

            if (verbose >= 2) print*, nx,ny, "....rt4 finished"

        else
            msg = 'no rt_mode selected'
            errorstatus = fatal
            call report(errorstatus,msg,nameOfRoutine)
            return
        end if
        !calculate human readable angles!
        angles_deg(1:NUMMU) = 180-(180.*acos(MU_VALUES(NUMMU:1:-1))/pi)
        angles_deg(1+NUMMU:2*NUMMU) = (180.*acos(MU_VALUES(1:NUMMU))/pi)

    end if

    if (verbose >= 1) call report(info,'End of ', nameOfRoutine)
    errorstatus = err
    return

end subroutine run_rt
