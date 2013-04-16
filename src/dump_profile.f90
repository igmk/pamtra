subroutine dump_profile()

  use vars_atmosphere
  use mod_io_strings
  use settings, only: tmp_path

  implicit none

  integer :: nz
  character(1) :: str1

  !      Preparation of the PROFILE file
  str1 = ''''
  file_profile = trim(tmp_path)//'/Profilex'//xstr//'y'//ystr//'f'//frq_str
  open(21, file = file_profile, form = 'FORMATTED', status =  'unknown')
  do nz = nlyr, 1, - 1 !nlyr,1,-1
     write(21,1013) hgt_lev(nz), temp_lev(nz), kextatmo(nz), str1//trim(file_ph(nz))//str1
  end do !end of cycle over the vertical layers
  write(21,1012) hgt_lev(0), temp_lev(0), KEXTATMO (1) , ''' '''
  close(21)

1012 format(f7.1,1x,f6.2,1x,E9.4,1x,a3)
1013 format(f7.1,1x,f6.2,1x,E9.4,1x,a38)

  return

end subroutine dump_profile
