module mod_io_strings

  implicit none
  save

  character(2) :: Nzstr
  character(3) :: xstr, ystr,&
       theta_str, H_str
  character(10) :: surf_type
  character(78) :: file_profile

  character(64), allocatable, dimension(:) :: file_PH

end module mod_io_strings
