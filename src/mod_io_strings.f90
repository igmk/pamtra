module mod_io_strings

  implicit none
  save

  character(2) :: Nzstr
  character(3) :: xstr, ystr,&
       theta_str, H_str
  character(10) :: surf_type
  character(78) :: file_profile

  character(64), allocatable, dimension(:) :: file_PH

  character(43) :: micro_str


contains
  subroutine mod_io_strings_get_filename()
    use nml_params
    use vars_atmosphere, only: year,month,day,time

    character(12) :: date_str
    !          character(43), intent(out) :: micro_str
    character :: SP_str*3,&
         N0snowstr*3, N0graustr*3, N0rainstr*3


    write (SP_str (1:3) , '(f3.1)') SP

    date_str = year//month//day//time

    if (N_0snowDsnow .le. 9.95) then 
       write (N0snowstr, '(f3.1)') N_0snowDsnow 
    else 
       write (N0snowstr, '(f3.0)') N_0snowDsnow 
    end if

    if (N_0rainD .le. 9.95) then 
       write (N0rainstr, '(f3.1)') N_0rainD 
    else 
       write (N0rainstr, '(f3.0)') N_0rainD 
    end if
    if (N_0grauDgrau .le. 9.95) then 
       write (N0graustr, '(f3.1)') N_0grauDgrau 
    else 
       write (N0graustr, '(f3.0)') N_0grauDgrau 
    end if

    micro_str = date_str//SD_snow//N0snowstr//EM_snow//SP_str//SD_grau//        &
         N0graustr//EM_grau//SD_rain//N0rainstr                            


  end subroutine mod_io_strings_get_filename

end module mod_io_strings
