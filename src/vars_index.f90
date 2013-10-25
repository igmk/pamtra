module vars_index
  !small simple module that every function has access to the current indices
  use kinds
  implicit none

  integer(kind=long) :: i_x= -99
  integer(kind=long) :: i_y= -99
  integer(kind=long) :: i_z= -99
  integer(kind=long) :: i_f= -99
  integer(kind=long) :: i_h= -99

end module vars_index