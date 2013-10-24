module mod_io_strings

  implicit none
  save

  character(2) :: Nzstr
  character(3) :: xstr, ystr,&
       theta_str, H_str
  character(4) :: nxstr, nystr
  character(8) :: frq_str
  character(10) :: surf_type
  character(78) :: file_profile

contains

function formatted_frqstr(input_string)

  implicit none

  integer :: pos

  character(4) :: tail
  character(3) :: head
  character(8) :: input_string, formatted_frqstr

  pos = index(input_string,'.')

  if (pos == 0) then
     input_string = input_string//'.'
     pos = len_trim(input_string)+1
  end if

  head = repeat('0', 3-len_trim(adjustl(input_string(:pos-1))))//input_string(:pos-1)
  tail = trim(input_string(pos+1:))//repeat('0', 4-len_trim(input_string(pos+1:)))

  formatted_frqstr = head//'.'//tail

  return

end function formatted_frqstr

end module mod_io_strings
