module vars_index
  !small simple module that every function has access to the current indices
  use kinds
  implicit none

  integer(kind=long) :: i_x
  integer(kind=long) :: i_y
  integer(kind=long) :: i_z
  integer(kind=long) :: i_f
  integer(kind=long) :: i_h

end module vars_index