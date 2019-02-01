
module vars_index
  !small simple module that every function has access to the current indices
  use kinds
  implicit none

  integer(kind=long) :: i_x= -99 !x direction
  integer(kind=long) :: i_y= -99 !y direction
  integer(kind=long) :: i_z= -99 !height
  integer(kind=long) :: i_f= -99 !frequency
  integer(kind=long) :: i_h= -99 !hydrometeor
  integer(kind=long) :: i_p= -99 !radar polarisation
  integer(kind=long) :: i_n= -99 !radar peak number

end module vars_index