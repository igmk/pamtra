subroutine interpolation(nx1,nx2,x1,y1,x2,y2)

!  use settings, only: verbose
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
    SUBROUTINE locate (xx, n, x, j) 
      use kinds                                                                 
      INTEGER j, n 
      REAL(kind=dbl) x, xx (n) 
    end SUBROUTINE locate 
  end interface


  if (verbose .gt. 1) print*, 'entering interpolation'

  ix2 = 0

  do i = 1, nx2
     call locate(x1,nx1,x2(i),ix2)
     y2(i) = (x2(i)-x1(ix2))*(y1(ix2+1)-y1(ix2))/(x1(ix2+1)-x1(ix2))+y1(ix2)
  end do

  return

end subroutine interpolation
