
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

  use kinds
  implicit none
  integer,parameter :: verbose = 1

  integer, intent(in) :: nx1,nx2
  real(kind=dbl), intent(in), dimension(nx1) :: x1,y1
  real(kind=dbl), intent(in), dimension(nx2) :: x2
  real(kind=dbl), intent(out), dimension(nx2) :: y2

  real(kind=dbl), dimension(nx2+1) :: x2_shift

  integer :: i
  integer, dimension(nx1+nx2+1) :: indarr
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
  indarr = (/ (i,i=1,nx1+nx2+1) /)
  x12_sorted = x12
  call Quicksort(nx1+nx2+1,x12_sorted, 1, nx1+nx2+1, indarr)
  y12_sorted=y12(indarr)

  call average_spectra(nx1+nx2+1,nx2+1,x12_sorted,y12_sorted,x2_shift,y2)
  print*, "x2_shift",x2_shift
  print*, "y_result",y2
  print*,"x2",x2
end subroutine rescale_spectra


subroutine average_spectra(nx12,nx2,x12_sorted,y12_sorted,x2,y_result)
  !averages the spectrum
  !borders of the averaged intervalls must be already present in in x12_sorted!
  !works only in combination with average_spectra

  !(c) M.Maahn, IGMK, 11/2012

!   use nml_params, only: verbose
  use kinds
  implicit none

  integer,parameter :: verbose = 1

  integer, intent(in) :: nx12,nx2
  real(kind=dbl), intent(in), dimension(nx12) :: x12_sorted,y12_sorted
  real(kind=dbl), intent(in), dimension(nx2) :: x2
  real(kind=dbl), intent(out), dimension(nx2-1) :: y_result

  integer :: ii,jj1,jj2
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
!       y_result(ii) = 0.d0
!       print*, "EXIT",y_result(ii)
!       CYCLE
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

end subroutine average_spectra

subroutine interpolate_spectra(nx1,nx2,x1,y1,x2,y2)

  !interpolation which gives back ZERO if x2 is more than 1bin out of range of x1!
  !values MUST be sorted
!   use nml_params, only: verbose
  use kinds
  implicit none
  integer,parameter :: verbose = 1

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
     print*,x2(i),ix2
     !points our of range?
     if ((ix2 .eq. 0) .or. (ix2 .eq. nx1+2)) then
      y2(i) = 0.d0
     else
      y2(i) = (x2(i)-x1_ext(ix2))*(y1_ext(ix2+1)-y1_ext(ix2))/(x1_ext(ix2+1)-x1_ext(ix2))+y1_ext(ix2)
     end if
  end do

  return

end subroutine interpolate_spectra


!                                                                       
      SUBROUTINE locate (xx, n, x, j) 
!    Given an array xx(1:n) and given a value x, returns a value j such 
!    that x is between xx(j) and xx(j+1) xx(1:n) must be monotonic, either
!    increasing or decreasing. j=0 or j=n is returned to indicate         
!    that x is out of range                                             
      use kinds                                                                 
      INTEGER j, n 
      REAL(kind=dbl) x, xx (n) 
      INTEGER jl, jm, ju 
                !initialise lower and                                   
      jl = 0 
                ! upper boundaries                                      
      ju = n + 1 
                           !if we are nmot yet done                     
   10 IF (ju - jl.gt.1) then 
                           !compute a midpoint                          
         jm = (ju + jl) / 2 
         IF ( (xx (n) .ge.xx (1) ) .eqv. (x.ge.xx (jm) ) ) then 
                   !and replace either the lower                        
            jl = jm 
         ELSE 
                    !or the upper limit                                 
            ju = jm 
         ENDIF 
         GOTO 10 
      ENDIF 
      IF (x.eq.xx (1) ) then 
         j = 1 
      ELSEIF (x.eq.xx (n) ) then 
         j = n - 1 
      ELSE 
         j = jl 
      ENDIF 
      RETURN 
      END SUBROUTINE locate                         


   !----------------------------------------------------------------------------
   !
   ! This file is based on the the routine in "Fortran 90 for Engineers & 
   ! Scientists" by Nyhoff and Leestma
   !
   ! Note: In the following subroutines, Item is an assumed-shape array
   !       so a program unit that calls these subroutines must:
   !	   1. contain this subroutine as an internal subprogram,
   !	   2. import this subroutine from a module, or
   !	   3. contain an interface block for this subroutine.
   !
   !----------------------------------------------------------------------------

   !-Quicksort------------------------------------------------------------------
   !
   ! Subroutine to sort a list using the quicksort method. Call it with 
   ! First = the lower bound on the subscripts of the array and 
   ! Last  = the upper bound. 
   !
   ! Accepts : Array "Item", array "Indices"
   ! Returns : Array "Item"    (modified) with elements in ascending order
   !           array "Indices" (modified) with elements 
   !----------------------------------------------------------------------------
   
   RECURSIVE SUBROUTINE Quicksort(Length,Item, First, Last, Indices)
   !----------------------------------------------------------------------------
   ! This routine is based on a similar routine in "Fortran 90 for Engineers & 
   ! Scientists" by Nyhoff and Leestma.  I modified it to return an integer 
   ! array sorted based on the relationship of the real data in "Item".
   !
   ! Example:
   ! real,    dimension(100) :: randvals,randcopy
   ! integer, dimension(100) :: indarr = (/ (i, i=1,100) /)
   !  ...
   ! call random_number(randvals)		! F90 intrinsic subroutine
   ! randcopy = randvals			! save for comparison
   ! call Quicksort(randvals,1,size(randvals),indarr)
   ! print *,'sorted - indexed original is ',SUM(randvals - randcopy(indarr))
   !
   ! TJH 21 Oct 1998
   !----------------------------------------------------------------------------
  use kinds
      INTEGER                     :: Length
      REAL(kind=dbl),    DIMENSION(Length), INTENT(INOUT) :: Item	! array of values
      INTEGER,               INTENT(IN)    :: First,Last
      INTEGER, DIMENSION(Length), INTENT(INOUT) :: Indices
   
      !--------------------------------------------------------------------
      ! Interface block(s) & Local Variables
      !--------------------------------------------------------------------

      INTERFACE
       SUBROUTINE Split(Item, Low, High, Mid, Indices)
use kinds
          REAL(kind=dbl),    DIMENSION(:), INTENT(INOUT) :: Item
          INTEGER,               INTENT(IN)    :: Low, High
          INTEGER,               INTENT(OUT)   :: Mid
          INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices
       END SUBROUTINE Split
      END INTERFACE

      INTEGER	:: Mid
      
      !--------------------------------------------------------------------

      IF (First < Last) THEN				! IF list size >= 2
         CALL Split(Item, First, Last, Mid, Indices)	! Split it
         CALL Quicksort(Length,Item, First, Mid-1, Indices)	! Sort left  half
         CALL Quicksort(Length,Item, Mid+1, Last,  Indices)	! Sort right half
      END IF

   END SUBROUTINE Quicksort

   !-Split----------------------------------------------------------------------
   !
   ! Subroutine to split a list into two sublists, using the first element 
   ! as a pivot, and return the position of the element about which the 
   ! list was divided. Local variables used are:
   ! Left	: position of the first element
   ! Right	: position of the last element
   ! Pivot	: pivot element
   ! Swap	: used to swap elements
   !
   ! Accepts:	Array Item and positions Low and High of the first and 
   !            last elements
   ! Returns:	Array Item (modified) with elements in ascending order
   !
   ! Note:	Item is an assumed-shape array so a program unit that calls
   !		this subroutine must:
   !		1. contain this subroutine as an internal subprogram,
   !		2. import this subroutine from a module
   !		3. contain an interface block for this subroutine.
   !----------------------------------------------------------------------------

   SUBROUTINE Split(Item, Low, High, Mid, Indices)
use kinds
      REAL(kind=dbl),    DIMENSION(:), INTENT(INOUT) :: Item
      INTEGER,               INTENT(IN)    :: Low, High
      INTEGER,               INTENT(OUT)   :: Mid
      INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices


      INTEGER ::   Left, Right
      REAL(kind=dbl)    ::  Pivot,  Swap
      INTEGER :: iPivot, iSwap

      Left   = Low
      Right  = High
      Pivot  = Item(Low)
      iPivot = Indices(Low)
   
      ! Repeat the following while Left and Right haven't met
   
      DO
         IF ( Left >= Right ) Exit
   
         ! Scan right to left to find element < Pivot
   
         DO
   	    IF ( Left >= Right .OR. Item(Right) < Pivot ) EXIT
   	    Right = Right - 1
         END DO

         ! Scan left to right to find element > Pivot

         DO
   	    IF (Item(Left) > Pivot) EXIT
   	    Left = Left + 1
         END DO
   
         ! If Left and Right haven't met, exchange the items
   
         IF (Left < Right) THEN
   	    Swap        = Item(Left)		! EXCHANGE THE ARRAY ITEMS
   	    Item(Left)  = Item(Right)
   	    Item(Right) = Swap

   	    iSwap          = Indices(Left)	! EXCHANGE THE INDICES ITEMS
   	    Indices(Left)  = Indices(Right)
   	    Indices(Right) = iSwap
         END IF
   
      END DO
   
      ! Switch element in split position with pivot
   
      Item(Low)   = Item(Right)			! SWITCH ARRAY ELEMS
      Item(Right) = Pivot
      Mid         = Right

      Indices(Low)   = Indices(Right)		! SWITCH ARRAY ELEMS
      Indices(Right) = iPivot

   END SUBROUTINE Split