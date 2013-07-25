function trapez_int(nbins,xs,ys,delta_xs)

    ! integrate using the trapez rule
    !
    ! Input:
    !  xs: x-values at the borders (nbins+1)
    !  ys: y-values at the borders (nbins+1)
    !  delta_xs: width of x values (nbins)
    !
    ! Result:
    !   trapez_int: result
    !
    ! References:

    use kinds
    implicit none

    integer, intent(in) :: nbins
    real(kind=dbl), dimension(nbins+1), intent(in) :: xs
    real(kind=dbl), dimension(nbins+1), intent(in) :: ys
    real(kind=dbl), dimension(nbins), intent(in) :: delta_xs

    integer :: ii
    real(kind=dbl) :: trapez_int
  
    trapez_int = 0.d0
    do ii = 1,nbins
      trapez_int = trapez_int + (ys(ii) + ys(ii+1)) / 2.d0 * delta_xs(ii) 
    end do
 


end function trapez_int