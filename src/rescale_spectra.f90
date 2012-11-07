
subroutine rescale_spectra(nx1,nx2,x1,y1,x2,y2)
  !(c) M.Maahn, IGMK, 11/2012

  !nx1,in: length of x1, y1
  !nx2,in: length of x2, y2
  !x1, in: x-values original data
  !y1, in: y-values original values
  !x2, in: center of new bins
  !y2, out: result

  !this routine averages and interpolates for every interval depending what is needed
  !in praxis, at the edges of the interval, values are interpolated  
  ! (thus we can be sure, that there are at least 2 values in every interval)
  !interpolated and original values are combined and sorted
  !then we average each interval.

  use nml_params, only: verbose
  use kinds
  implicit none


  integer, intent(in) :: nx1,nx2
  real(kind=dbl), intent(in), dimension(nx1) :: x1,y1
  real(kind=dbl), intent(in), dimension(nx2) :: x2
  real(kind=dbl), intent(out), dimension(nx2) :: y2

  real(kind=dbl), dimension(nx2+1) :: x2_shift

  integer :: i
  real(kind=dbl), dimension(nx1+nx2+1) :: x12,x12_sorted
  real(kind=dbl), dimension(nx1+nx2+1) :: y12,y12_sorted
  real(kind=dbl), dimension(nx2+1) :: y2_interp

  x2_shift(2:nx2) = x2(2:) - 0.5d0*(x2(2:) - x2(1:nx2-1))
  x2_shift(1) = x2(1) -  0.5d0*(x2(2) - x2(1))
  x2_shift(nx2+1) = x2(nx2) +  0.5d0*(x2(nx2) - x2(nx2-1))  


  call interpolate_spectra(nx1,nx2,x1,y1,x2_shift,y2_interp)

  !join interpolated and original array
  x12(1:nx1)=x1
  x12(nx1+1:nx1+nx2+1)=x2_shift
  y12(1:nx1)=y1
  y12(nx1+1:nx1+nx2+1)=y2_interp

  !make order right
  x12_sorted = x12
  y12_sorted=y12
  call dsort(x12_sorted, y12_sorted, nx1+nx2+1, 2)


  call average_spectra(nx1+nx2+1,nx2+1,x12_sorted,y12_sorted,x2_shift,y2)

end subroutine rescale_spectra


subroutine average_spectra(nx12,nx2,x12_sorted,y12_sorted,x2,y_result)
  !averages the spectrum
  !borders of the averaged intervalls must be already present in in x12_sorted!
  !works only in combination with average_spectra

  !(c) M.Maahn, IGMK, 11/2012

  use nml_params, only: verbose
  use kinds
  implicit none

  integer, intent(in) :: nx12,nx2
  real(kind=dbl), intent(in), dimension(nx12) :: x12_sorted,y12_sorted
  real(kind=dbl), intent(in), dimension(nx2) :: x2
  real(kind=dbl), intent(out), dimension(nx2-1) :: y_result

  integer :: ii,jj1,jj2

  if (verbose .gt. 1) print*, 'entering rescale_spectra averaging'

  !step zero
  call locate (x12_sorted, nx12, x2(1), jj1)
  do ii=1,nx2-1
    !find indices
    call locate (x12_sorted, nx12, x2(ii+1), jj2) 
    !locate does not work properly if an entry is searched EQUAL to the last
    if (jj1 .eq. jj2) then
      jj2 = jj1 +1 
    end if

    if ((jj2 .gt. nx12) .or. (jj1 .gt. nx12)) then
      print*, "WARNING: Exiting averaging loop"
      exit
    end if
!     end if 
    !make the averaging, first width half of the weights on the left side
    y_result(ii) = SUM(y12_sorted(jj1:jj2-1) * 0.5d0*(x12_sorted(jj1+1:jj2)-x12_sorted(jj1:jj2-1)  ) ) 
    !now right-side weights
    y_result(ii) = y_result(ii) + SUM(y12_sorted(jj1+1:jj2) * 0.5d0*(x12_sorted(jj1+1:jj2)-x12_sorted(jj1:jj2-1)  ) ) 
    !devide by weights
    y_result(ii) = y_result(ii) / SUM((x12_sorted(jj1+1:jj2)-x12_sorted(jj1:jj2-1)))
    !save idnex for next iteration
    jj1=jj2
  end do

  if (verbose .gt. 1) print*, 'exiting rescale_spectra averaging'

  return

end subroutine average_spectra

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


  if (verbose .gt. 1) print*, 'entering interpolation'



  !extend x1 and y1 to give a "0" reference point
  x1_ext(1)     = x1(1) - (x1(2)-x1(1))
  x1_ext(nx1+2) = x1(nx1) +  (x1(2)-x1(1))
  y1_ext(:) = 0.d0 

  x1_ext(2:nx1+1) = x1
  y1_ext(2:nx1+1) = y1

  ix2 = 0

  do i = 1, nx2
     call locate(x1_ext,nx1+2,x2(i),ix2)
     !points our of range?
     if ((ix2 .eq. 0) .or. (ix2 .eq. nx1+2)) then
      y2(i) = 0.d0
     else
      y2(i) = (x2(i)-x1_ext(ix2))*(y1_ext(ix2+1)-y1_ext(ix2))/(x1_ext(ix2+1)-x1_ext(ix2))+y1_ext(ix2)
     end if
  if (verbose .gt. 5) print*, "y2(i),i",y2(i),i
  end do

  if (verbose .gt. 1) print*, 'exiting interpolation'

  return

end subroutine interpolate_spectra