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