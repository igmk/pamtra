! ----------------------------------------------------------------------------
! subroutine AVINT
! ----------------------------------------------------------------------------
  subroutine avint ( ftab, xtab, ntab, a_in, b_in, result )
  implicit none
!
! Purpose:
!   estimate the integral of unevenly spaced data
!
! Inputs:
!   [ftab]     functional values
!   [xtab]     abscissa values
!   [ntab]     number of elements of [ftab] and [xtab]
!   [a]        lower limit of integration
!   [b]        upper limit of integration
!
! Outputs:
!   [result]   approximate value of integral
!
! Reference:
!   From SLATEC libraries, in public domain
!   Via Quickbeam
!
!***********************************************************************
!
!  AVINT estimates the integral of unevenly spaced data.
!
!  Discussion:
!
!    The method uses overlapping parabolas and smoothing.
!
!  Modified:
!
!    30 October 2000
!
!  Reference:
!
!    Philip Davis and Philip Rabinowitz,
!    Methods of Numerical Integration,
!    Blaisdell Publishing, 1967.
!
!    P E Hennion,
!    Algorithm 77,
!    Interpolation, Differentiation and Integration,
!    Communications of the Association for Computing Machinery,
!    Volume 5, page 96, 1962.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) FTAB(NTAB), the function values,
!    FTAB(I) = F(XTAB(I)).
!
!    Input, real ( kind = 8 ) XTAB(NTAB), the abscissas at which the
!    function values are given.  The XTAB's must be distinct
!    and in ascending order.
!
!    Input, integer NTAB, the number of entries in FTAB and
!    XTAB.  NTAB must be at least 3.
!
!    Input, real ( kind = 8 ) A, the lower limit of integration.  A should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Input, real ( kind = 8 ) B, the upper limit of integration.  B should
!    be, but need not be, near one endpoint of the interval
!    (X(1), X(NTAB)).
!
!    Output, real ( kind = 8 ) RESULT, the approximate value of the integral.

  integer, intent(in) :: ntab

  real ( kind = 8 ), intent(in) :: a_in
  real ( kind = 8 ) a
  real ( kind = 8 ) atemp
  real ( kind = 8 ), intent(in) :: b_in
  real ( kind = 8 ) b
  real ( kind = 8 ) btemp
  real ( kind = 8 ) ca
  real ( kind = 8 ) cb
  real ( kind = 8 ) cc
  real ( kind = 8 ) ctemp
  real ( kind = 8 ), intent(in) :: ftab(ntab)
  integer i
  integer ihi
  integer ilo
  integer ind
  real ( kind = 8 ), intent(out) :: result
  real ( kind = 8 ) sum1
  real ( kind = 8 ) syl
  real ( kind = 8 ) term1
  real ( kind = 8 ) term2
  real ( kind = 8 ) term3
  real ( kind = 8 ) x1
  real ( kind = 8 ) x2
  real ( kind = 8 ) x3
  real ( kind = 8 ), intent(in) :: xtab(ntab)

  a = a_in
  b = b_in  
  
  if ( ntab < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Fatal error!'
    write ( *, '(a,i6)' ) '  NTAB is less than 3.  NTAB = ', ntab
    stop
  end if
 
  do i = 2, ntab
 
    if ( xtab(i) <= xtab(i-1) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'AVINT - Fatal error!'
      write ( *, '(a)' ) '  XTAB(I) is not greater than XTAB(I-1).'
      write ( *, '(a,i6)' ) '  Here, I = ', I
      write ( *, '(a,g14.6)' ) '  XTAB(I-1) = ', xtab(i-1)
      write ( *, '(a,g14.6)' ) '  XTAB(I) =   ', xtab(i)
      stop
    end if
 
  end do
 
  result = 0.0D+00
 
  if ( a == b ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'AVINT - Warning!'
    write ( *, '(a)' ) '  A = B, integral=0.'
    return
  end if
!
!  If B < A, temporarily switch A and B, and store sign.
!
  if ( b < a ) then
    syl = b
    b = a
    a = syl
    ind = -1
  else
    syl = a
    ind = 1
  end if
!
!  Bracket A and B between XTAB(ILO) and XTAB(IHI).
!
  ilo = 1
  ihi = ntab

  do i = 1, ntab
    if ( a <= xtab(i) ) then
      exit
    end if
    ilo = ilo + 1
  end do

  ilo = max ( 2, ilo )
  ilo = min ( ilo, ntab - 1 )

  do i = 1, ntab
    if ( xtab(i) <= b ) then
      exit
    end if
    ihi = ihi - 1
  end do
  
  ihi = min ( ihi, ntab - 1 )
  ihi = max ( ilo, ihi - 1 )
!
!  Carry out approximate integration from XTAB(ILO) to XTAB(IHI).
!
  sum1 = 0.0D+00
 
  do i = ilo, ihi
 
    x1 = xtab(i-1)
    x2 = xtab(i)
    x3 = xtab(i+1)
    
    term1 = ftab(i-1) / ( ( x1 - x2 ) * ( x1 - x3 ) )
    term2 = ftab(i)   / ( ( x2 - x1 ) * ( x2 - x3 ) )
    term3 = ftab(i+1) / ( ( x3 - x1 ) * ( x3 - x2 ) )
 
    atemp = term1 + term2 + term3

    btemp = - ( x2 + x3 ) * term1 &
            - ( x1 + x3 ) * term2 &
            - ( x1 + x2 ) * term3

    ctemp = x2 * x3 * term1 + x1 * x3 * term2 + x1 * x2 * term3
 
    if ( i <= ilo ) then
      ca = atemp
      cb = btemp
      cc = ctemp
    else
      ca = 0.5D+00 * ( atemp + ca )
      cb = 0.5D+00 * ( btemp + cb )
      cc = 0.5D+00 * ( ctemp + cc )
    end if
 
    sum1 = sum1 &
          + ca * ( x2**3 - syl**3 ) / 3.0D+00 &
          + cb * 0.5D+00 * ( x2**2 - syl**2 ) &
          + cc * ( x2 - syl )
 
    ca = atemp
    cb = btemp
    cc = ctemp
 
    syl = x2
 
  end do
 
  result = sum1 &
        + ca * ( b**3 - syl**3 ) / 3.0D+00 &
        + cb * 0.5D+00 * ( b**2 - syl**2 ) &
        + cc * ( b - syl )
!
!  Restore original values of A and B, reverse sign of integral
!  because of earlier switch.
!
  if ( ind /= 1 ) then
    ind = 1
    syl = b
    b = a
    a = syl
    result = -result
  end if
 
  return
  end subroutine avint
  