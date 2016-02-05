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


!  if (verbose .gt. 1) print*, 'entering interpolation'

  ix2 = 0

  do i = 1, nx2
     call locate(x1,nx1,x2(i),ix2)
     y2(i) = (x2(i)-x1(ix2))*(y1(ix2+1)-y1(ix2))/(x1(ix2+1)-x1(ix2))+y1(ix2)
  end do

  return

end subroutine interpolation

                                                                       
!********************************************************************** 
!    CALCULATION location for interpolation                           * 
!    Copr. 1986-92 Numerical Recipes Software +>k-5V1`                * 
!                                                                     * 
!********************************************************************** 
!                                                                       
      SUBROUTINE locate (xx, n, x, j) 
!    Given an array xx(1:n) and given a value x, returns a value j such 
!    that x is between xx(j) and xx(j+1) xx(1:n) must be monotonic, eith
!  increasing or decreasing. j=0 or j=n is returned to indicate         
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

