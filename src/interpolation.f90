subroutine interpolation(nx1,nx2,x1,y1,x2,y2)

  use kinds
  use report_module

  implicit none

  integer :: i
  integer :: nx1,nx2

  integer:: ix2

  real(kind=dbl), intent(in), dimension(nx1) :: x1,y1
  real(kind=dbl), intent(in), dimension(nx2) :: x2
  real(kind=dbl), intent(out), dimension(nx2) :: y2
  interface
    subroutine find_index(xx, n, x, j) 
      use kinds                                                                 
      integer :: j, n 
      real(kind=dbl) :: x
      real(kind=dbl), dimension(n) :: xx 
    end subroutine find_index 
  end interface


  ix2 = 0

  do i = 1, nx2
     call find_index(x1,nx1,x2(i),ix2)
     y2(i) = (x2(i)-x1(ix2))*(y1(ix2+1)-y1(ix2))/(x1(ix2+1)-x1(ix2))+y1(ix2)
  end do

  return

end subroutine interpolation

subroutine find_index(xx,n,x,j)

  !    Given an array xx(1:n) and given a value x, returns a value j such 
  !    that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, either
  !    increasing or decreasing. j=0 or j=n is returned to indicate         
  !    that x is out of range                                             

   use kinds

   implicit none
   
   integer :: j, n, st, si 
   real(kind=dbl) :: x
   real(kind=dbl), dimension(n) :: xx 

   st = minloc(abs(xx - x), 1)
      
   if ((xx(1) < xx(n)) .eqv. (xx(st) < x)) then
     si = st
   else
     si = st - 1 
   end if
   if (x == xx(1)) then
     j = 1
   elseif (x == xx(n)) then
     j = n - 1
   else
     j = si
   end if

   return

end subroutine find_index