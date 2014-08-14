subroutine dump_scatProfile(errorstatus,mod)

  use kinds, only: long
  use vars_atmosphere
  use vars_rt, only: rt_kextatmo!, rt_file_ph
  use mod_io_strings
  use settings, only: tmp_path
  use vars_index
  use report_module

  implicit none

  integer :: nz, mod
  character(1) :: str1
  
   integer(kind=long) :: errorstatus
   integer(kind=long) :: err = 0
   character(len=80) :: msg
  character(len=40) :: nameOfRoutine = 'dump_scatProfile'

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)
    
  !      Preparation of the PROFILE file
  select case (mod)
    case (0)
      str1 = ''''
      file_profile = trim(tmp_path)//'Profilex'//xstr//'y'//ystr//'f'//frq_str
      open(21, file = file_profile, form = 'FORMATTED', status =  'unknown')
!  do nz = atmo_nlyrs(i_x,i_y), 1, - 1 !nlyr,1,-1
    case (1)
      write(nzstr,'(i2.2)') i_z
      write(21,1013) atmo_hgt_lev(i_x,i_y,i_z+1), atmo_temp_lev(i_x,i_y,i_z+1), rt_kextatmo(i_z),&
       &trim(file_profile)//'_'//nzstr//str1
!  end do !end of cycle over the vertical layers
    case (2)
      write(21,1012) atmo_hgt_lev(i_x,i_y,1), atmo_temp_lev(i_x,i_y,1), rt_kextATMO (1) , ''' '''
      close(21)
    end select
    
  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

1012 format(f7.1,1x,f6.2,1x,E9.4,1x,a3)
1013 format(f7.1,1x,f6.2,1x,E9.4,1x,a50)
  return

end subroutine dump_scatProfile
