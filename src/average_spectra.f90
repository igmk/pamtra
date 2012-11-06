subroutine average_spectra(nx1,nx2,x1,y1,x2,y2)

  !interpolation which gives back ZERO if x2 is more than 1bin out of range of x1!
  !values MUST be sorted

  !assumes evenly spaced x2 and delta(x1) < delta(x2)

  !(c) M.Maahn, IGMK, 11/2012

  use nml_params, only: verbose
  use kinds
  implicit none

  integer :: ii,jj,kk,nn, counter
  integer, intent(in) :: nx1,nx2

  integer:: ix2
  real(kind=dbl):: del_x1_lower, del_x1_upper, weight

  real(kind=dbl), intent(in), dimension(nx1) :: x1,y1
  real(kind=dbl), dimension(:), allocatable :: x1_ext,y1_ext,weights
  real(kind=dbl), intent(in), dimension(nx2) :: x2
  real(kind=dbl), dimension(nx2+1) :: x2_shift
  real(kind=dbl), intent(out), dimension(nx2) :: y2
  logical :: first

  if (verbose .gt. 1) print*, 'entering average_spectra'



  !extend x1 and y1 to give a "0" reference point

  if (MAXVAL(x1(2:) - x1(:nx1-1)) .ge. MINVAL(x2(2:) - x2(:nx2-1))) then
    print*, "Error in average_spectra.f90: function average_spectra must not be used for interpolation"
    stop
  end if

  !we have too make this complicated to exclude numerical noise...
  if (ANY(ABS((x2(2:) - x2(:nx2-1)) - (x2(2)-x2(1))) .gt. ABS(x2(2)-x2(1))*0.00001 )) then
    print*, "Error in average_spectra.f90: function average_spectra assumes evenly spaced x2!"
    print*, (x2(2:) - x2(:nx2-1)) .ne. (x2(2)-x2(1))
    print*, (x2(2:) - x2(:nx2-1)) 
    stop
  end if  
  
  del_x1_lower = x1(2)-x1(1)
  del_x1_upper = x1(nx1)-x1(nx1-1)

  ! find out how much we have to extend x1 to cover the bins, which are going to be averaged, completely.
  nn = 1
  do while ((x1(nx1) + (nn* del_x1_upper) - (x1(1) - (nn*del_x1_lower))) .le. 3 * (x2(2)-x2(1)))
    nn = nn + 2 !save time by using 2 instead of 1
  end do 
  !now allocate the found length
  allocate(x1_ext(nx1+2*nn))
  allocate(y1_ext(nx1+2*nn))
  allocate(weights(nx1+2*nn))
  !fill the extended x1
  do kk = 1, nn
    if (verbose .gt. 4) print*,"kk,nn+1-kk,nx1+kk,x1(1),(kk*del_x1_lower),x1(nx1),",&
        "(kk* del_x1_upper),x1(1) - (kk*del_x1_lower),x1(nx1) + (kk* del_x1_upper)"
    if (verbose .gt. 4) print*,kk,nn+1-kk,nx1+kk,x1(1),(kk*del_x1_lower),x1(nx1),&
	(kk* del_x1_upper),x1(1) - (kk*del_x1_lower),x1(nx1) + (kk* del_x1_upper)
    x1_ext(nn+1-kk)     = x1(1) - (kk*del_x1_lower)
    x1_ext(nn+nx1+kk) = x1(nx1) + (kk* del_x1_upper)
  end do 
  y1_ext(:) = 0.d0 

  !put in the original values
  x1_ext(nn+1:nx1+nn) = x1
  y1_ext(nn+1:nx1+nn) = y1

  !their weights are the width of the spectrum they are covering, thus delta x1
  weights(1:nx1+2*nn-1) = x1_ext(2:nx1+2*nn) - x1_ext(1:nx1+2*nn-1)
  weights(nx1+2*nn) = weights(nx1+2*nn-1)

  if (ANY(weights .eq. 0.d0)) then
    print*,"Error in average_spectra.f90: weights must not be 0!"
    stop
  end if

  !shift the bins of x2 from their center to the left edge
  x2_shift(1:nx2) = x2 - (x2(2)-x2(1))/2.d0  !now shifted from the center to the edge!
  x2_shift(1+nx2) = x2(nx2) + (x2(2)-x2(1))/2.d0  !last entry


  y2(:) = 0.d0
  weight = 0.d0
  jj = 1
  first = .true.
  !find a place for every y1:
  do ii = 1, nx1+2*nn
   if (verbose .gt. 4) print*,"iterate",x1_ext(ii), y1_ext(ii)
   !skip x2 where no x1 is present
   do while (ALL(x2_shift(jj+1) .lt. x1_ext) .and. (jj .le. nx2))
     if (verbose .gt. 4) print*,"skip",ii,jj
     jj = jj+1
   end do
   !if x1 was found to be between xs(jj) and x2(jj+1)
   if ((x1_ext(ii) .ge. x2_shift(jj)) .and. (x1_ext(ii) .lt. x2_shift(jj+1)) ) then
    !add y1 to y2
    if (verbose .gt. 4) print*,"add",ii,jj,y2(jj),y1_ext(ii),weights(ii),y1_ext(ii)*weights(ii)
    y2(jj) = y2(jj) + y1_ext(ii)*weights(ii)
    weight = weight +weights(ii)
    first = .false.
   !skip x1 which are smaller than range of x2
   else if (first) then
     print*,"Cycle",ii,jj
     CYCLE
   else
    !now calculate the mean
    if (verbose .gt. 4) print*,"CALC",y2(jj)/weight
    y2(jj) = y2(jj)/weight
    !for the next x2 position, iterate:
    jj = jj+1
    if (verbose .gt. 4) print*,"new",ii,jj
    !exit if jj is greater than avaiable x2 (this also removes x1 greater than range of x2)
    if (jj .gt. nx2) then
    if (verbose .gt. 4) print*,"EXIT",ii,jj
      EXIT
    end if
    !start to collect x1 for the next x2 position
    if (verbose .gt. 4) print*,"add",ii,jj,y2(jj),y1_ext(ii),weights(ii),y1_ext(ii)*weights(ii)
    y2(jj) = y2(jj) + y1_ext(ii)*weights(ii)
    weight = weights(ii)

    end if

  end do
  !for the very last bin:
  if (jj .le. nx2) y2(jj) = y2(jj)/weight



  deallocate(x1_ext,y1_ext)

  if (verbose .gt. 1) print*, 'exiting average_spectra'



  return

end subroutine average_spectra