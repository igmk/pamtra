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

    allocate(nlegen(nlyr),stat=alloc_status)
  
    allocate(kextatmo(nlyr), stat=alloc_status)
    allocate(kexttot(nlyr), stat=alloc_status)
    allocate(kextsn(nlyr), stat=alloc_status)
    allocate(kextcw(nlyr), stat=alloc_status)
    allocate(kextrr(nlyr), stat=alloc_status)
    allocate(kextgr(nlyr), stat=alloc_status)
    allocate(kextci(nlyr), stat=alloc_status)
    allocate(kextha(nlyr), stat=alloc_status)
  
    allocate(salbtot(nlyr), stat=alloc_status)
    allocate(g_coeff(nlyr), stat=alloc_status)
  
    allocate(back(nlyr), stat=alloc_status)
    allocate(backcw(nlyr), stat=alloc_status)
    allocate(backrr(nlyr), stat=alloc_status)
    allocate(backci(nlyr), stat=alloc_status)
    allocate(backsn(nlyr), stat=alloc_status)
    allocate(backgr(nlyr), stat=alloc_status)
    allocate(backha(nlyr), stat=alloc_status)

    allocate(legen(nlyr,200), stat=alloc_status)
    allocate(legen2(nlyr,200), stat=alloc_status)
    allocate(legen3(nlyr,200), stat=alloc_status)
    allocate(legen4(nlyr,200), stat=alloc_status)

    allocate(rt4salbtot(nlyr),stat=alloc_status)
    allocate(rt4scatter_matrix(nlyr,nstokes,nummu,nstokes,nummu,4),stat=alloc_status)
    allocate(scattermatrix(nlyr,nstokes,nummu,nstokes,nummu,4),stat=alloc_status)
    allocate(rt4ext_matrix(nlyr,nstokes,nstokes,nummu,4),stat=alloc_status)
    allocate(extmatrix(nlyr,nstokes,nstokes,nummu,4),stat=alloc_status)
    allocate(rt4emis_vec(nlyr,nstokes,nummu,4),stat=alloc_status)
    allocate(emisvec(nlyr,nstokes,nummu,4),stat=alloc_status)




    allocate(hydros_present(nlyr),stat=alloc_status)
    allocate(rt4hydros_present(nlyr),stat=alloc_status)

    if (dump_to_file) then
        allocate(file_ph(nlyr))
    end if

    ! set them to zero, just in case they are not calculated but used for Ze/PIA calculation
    kexttot(:) = 0d0
    kextatmo(:) = 0d0
    back(:) = 0d0

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

    errorstatus = err

    return

end subroutine allocate_profile_vars
