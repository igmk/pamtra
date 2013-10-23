subroutine allocate_profile_vars(errorstatus)

    use vars_atmosphere
    use vars_output
    use settings
    use mod_io_strings
    use kinds, only: long
    use report_module

    implicit none

    ! Error handling

    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=200) :: msg
    character(len=30) :: nameOfRoutine = 'allocate_profile_vars'

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    allocate(hgt_lev(0:nlyr),stat=alloc_status)
    allocate(press_lev(0:nlyr),stat=alloc_status)
    allocate(press(nlyr),stat=alloc_status)
    allocate(temp_lev(0:nlyr),stat=alloc_status)
    allocate(temp(nlyr),stat=alloc_status)
    allocate(hgt(nlyr),stat=alloc_status)
    allocate(delta_hgt_lev(nlyr),stat=alloc_status)
    allocate(relhum_lev(0:nlyr),stat=alloc_status)
    allocate(relhum(nlyr),stat=alloc_status)

    allocate(vapor_pressure(nlyr),stat=alloc_status)
    allocate(rho_vap(nlyr),stat=alloc_status)
    allocate(q_hum(nlyr),stat=alloc_status)

    allocate(q_hydro(5,nlyr))


    allocate(cwc_q(nlyr))
    allocate(iwc_q(nlyr))
    allocate(rwc_q(nlyr))
    allocate(swc_q(nlyr))
    allocate(gwc_q(nlyr))
    if (n_moments .eq. 2) then
        allocate(hwc_q(nlyr))
        allocate(cwc_n(nlyr))
        allocate(iwc_n(nlyr))
        allocate(rwc_n(nlyr))
        allocate(swc_n(nlyr))
        allocate(gwc_n(nlyr))
        allocate(hwc_n(nlyr))
    end if


    if (dump_to_file) then
        allocate(file_ph(nlyr))
    end if



    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

    errorstatus = err

    return

end subroutine allocate_profile_vars
