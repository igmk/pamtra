subroutine get_gasabs &
(errorstatus,&      ! out
freq)              ! in
    ! Description:
    !  Based on on the namelist parameter gas_mod, the routine for calculating the
    !  gas absorption due to air (N2 and O2) and water vapor for the input frequency
    !  is selected.
    !
    !  The absorption coefficient is converted from Np/km to Np/m in the end.
    !
    !  The result is stored in kextatmo from vars_atmosphere module.
    !
    ! Input parameters:
    !    frequency                    GHz
    ! Output parameters:
    !    errorstatus
    !
    ! Owner: IGMK
    !
    ! History:
    !
    ! Version Date         Comment
    ! ------- ----         -------
    ! 0.01    21/09/2008   Code adaption from A. Battaglia - M. Mech
    ! 0.1     21/09/2009   Include Rosenkranz gas absorption model - M. Mech
    ! 0.2     26/02/2013   Application of European Standards for Writing and
    !                      Documenting Exchangeable Fortran 90 Code - M. Mech
    !
    ! Code Description:
    !  Language: Fortran 90.
    !  Software Standards: "European Standards for Writing and
    !   Documenting Exchangeable Fortran 90 Code".
    !
    ! Parent Module: run_rt
    !
    ! Child modules: mpm93
    !                rosen_gasabs98
    !
    ! Declarations:

    ! Imported Parameters:
    ! Modules used:
    use kinds, only: dbl, & ! integer parameter specifying double precision
    long   ! integer parameter specifying long integer
    use vars_atmosphere, only: nlyr, press, temp, vapor_pressure, rho_vap, kextatmo
    use constants, only: t_abs
    use settings, only: gas_mod

    use report_module

    implicit none

    !- End of header ---------------------------------------------------------------

    ! Subroutine arguments
    ! Scalar arguments with intent(in):
    real(kind=dbl), intent(in) :: freq
    ! End of Subroutine arguments

    ! Local scalars:
    integer(kind=long) :: nz

    real(kind=dbl) :: absair, abswv
    real(kind=dbl) :: tc

    ! Error handling

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'get_gasabs'

    if (verbose >= 1) call report(info,'Start of ', nameOfRoutine)

    do nz = 1, nlyr
        tc = temp(nz) - t_abs
        if (gas_mod .eq. 'L93') then
            call mpm93(err,freq, press(nz)/1.d3, vapor_pressure(nz)/1.d3,tc, 0.d0, kextatmo(nz))
            if (err /= 0) then
                msg = 'error in mpm93!'
                call report(err, msg, nameOfRoutine)
                errorstatus = err
                return
            end if
            kextatmo(nz) = kextatmo(nz)/1.d3
        else if (gas_mod .eq. 'R98') then
            call rosen98_gasabs(err,freq,temp(nz),rho_vap(nz),press(nz),absair,abswv)
            if (err /= 0) then
                msg = 'Error in rosen98_gasabs!'
                call report(err, msg, nameOfRoutine)
                errorstatus = err
                return
            end if
            kextatmo(nz) = (absair + abswv)/1.d3    ! conversion to Np/m
        else
            kextatmo(nz) = 0
            msg = 'No gas absorption model specified!'
            err = fatal
            call report(err, msg, nameOfRoutine)
            errorstatus = err
            return
        end if
    end do

    errorstatus = err

    if (verbose >= 1) call report(info,'End of ', nameOfRoutine)

    return

end subroutine get_gasabs
