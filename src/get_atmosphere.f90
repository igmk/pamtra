subroutine get_atmosphere(errorstatus)

    use nml_params, only: input_type
    use vars_profile
    use kinds, only: long

    implicit none

    integer(kind=long), intent(out) :: errorstatus

    if (trim(input_type) .eq. 'profile') then
        call vars_profile_read_profile(errorstatus)
    else if (trim(input_type) .eq. 'cosmo') then
        call vars_profile_read_cosmo
    else
        stop 'input type not recognized!'
    end if

    return

end subroutine get_atmosphere
