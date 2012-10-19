
subroutine interpolate_spectra(nx1,nx2,x1,y1,x2,y2)

  !interpolation which gives back ZERO if x2 is more than 1bin out of range of x1!
  !values MUST be sorted
  use nml_params, only: verbose
  use kinds
  implicit none

  integer :: i
  integer :: nx1,nx2

  integer:: ix2

  real(kind=dbl), intent(in), dimension(nx1) :: x1,y1
  real(kind=dbl), dimension(nx1+2) :: x1_ext,y1_ext
  real(kind=dbl), intent(in), dimension(nx2) :: x2
  real(kind=dbl), intent(out), dimension(nx2) :: y2

  real(kind=dbl) :: nan

  if (verbose .gt. 1) print*, 'entering interpolation'

  nan = 0.d0
  nan = nan/nan

  !extend x1 and y1 to make sure that interpolated bins distance is not larger than whole original dataset
  x1_ext(1)     = x1(1) - (x1(2)-x1(1))
  x1_ext(nx1+2) = x1(nx1) +  (x1(2)-x1(1))
  y1_ext(:) = 0.d0 

  x1_ext(2:nx1+1) = x1
  y1_ext(2:nx1+1) = y1

  ix2 = 0

  do i = 1, nx2
     !points our of range?
     if ((x2(i) .gt. MAXVAL(x1_ext)) .or. (x2(i) .lt. MINVAL(x1_ext))) then
       y2(i) = 0.d0
     else
      call locate(x1_ext,nx1+2,x2(i),ix2)
      y2(i) = (x2(i)-x1_ext(ix2))*(y1_ext(ix2+1)-y1_ext(ix2))/(x1_ext(ix2+1)-x1_ext(ix2))+y1_ext(ix2)
     end if
  end do

  return

end subroutine interpolate_spectra