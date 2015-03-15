

subroutine gausquad (n, xa, wt) 
    !      generates the abscissas (x) and weights (w) for an n point
    !      gauss-legendre quadrature.
    use kinds

    implicit none
    integer :: n
    real(kind=dbl) :: xa(*), wt(*)
    integer :: k, i, j, l
    real(kind=dbl) :: x, xp, pl, pl1, pl2, dpl
    real(kind=dbl), parameter :: tiny = 3.0d-13

    k = (n + 1) / 2
    do j = 1, k
        x = cos(3.141592654 * (j - .25) / (n + .5))
        i = 0
100 continue
    pl1 = 1
    pl = x
    do l = 2, n
        pl2 = pl1
        pl1 = pl
        pl = ( (2 * l - 1) * x * pl1 - (l - 1) * pl2) / l
    end do
    dpl = n * (x * pl - pl1) / (x * x - 1)
    xp = x
    x = xp - pl / dpl
    i = i + 1
    if (abs (x - xp) .gt. tiny .and. i .lt. 10) goto 100
    xa(j) = - x
    xa(n - j + 1) = x
    wt(j) = 2.0d0 / ( (1.0d0 - x * x) * dpl * dpl)
    wt(n - j + 1) = wt(j)
end do

return
end subroutine gausquad
