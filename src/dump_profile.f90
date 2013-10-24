subroutine dump_profile()

  use vars_atmosphere
  use vars_rt, only: rt_kextatmo, rt_file_ph
  use mod_io_strings
  use settings, only: tmp_path
  use vars_index

  implicit none

  integer :: nz
  character(1) :: str1

  !      Preparation of the PROFILE file
  str1 = ''''
  file_profile = trim(tmp_path)//'/Profilex'//xstr//'y'//ystr//'f'//frq_str
  open(21, file = file_profile, form = 'FORMATTED', status =  'unknown')
  do nz = atmo_nlyrs(i_x,i_y), 1, - 1 !nlyr,1,-1
     write(21,1013) atmo_hgt_lev(i_x,i_y,nz+1), atmo_temp_lev(i_x,i_y,nz+1), rt_kextatmo(nz), str1//trim(rt_file_ph(nz))//str1
  end do !end of cycle over the vertical layers
  write(21,1012) atmo_hgt_lev(i_x,i_y,1), atmo_temp_lev(i_x,i_y,1), rt_kextATMO (1) , ''' '''
  close(21)

1012 format(f7.1,1x,f6.2,1x,E9.4,1x,a3)
1013 format(f7.1,1x,f6.2,1x,E9.4,1x,a38)

  return

end subroutine dump_profile
