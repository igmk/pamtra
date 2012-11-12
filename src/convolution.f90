! Convolution demonstration using a bidirectional full-systolic convolution array.
!
! Author : Erik Urbach (0906204)
!
! Course : Systolic Parallel Processing
!
! Date   : Tue Nov  2 11:08:46 MET 1999
!
! Note   : See figure 6.4 from "Systolic Parallel Processing" by N.Petkov.
!

subroutine convolution(X,M,A,N,Y)

  use kinds
  implicit none
  
   INTEGER, intent(in) :: M  ! Size of input vector X
   INTEGER, intent(in) :: N   ! Size of convolution filter A
   INTEGER :: I

   REAL(kind=dbl), intent(in), DIMENSION(M) :: X
   REAL(kind=dbl), intent(in), DIMENSION(N) :: A
   REAL(kind=dbl), DIMENSION(N) :: DX, DY
   REAL(kind=dbl), intent(out), DIMENSION(M+N-1) :: Y

! !  In the next two lines you can specify the contents of X and A.
! !  Don't forget to change M and N when you change the size of X and A.
!    X = (/3, 5, 2, 7/)  ! Contents of input vector X
!    A = (/9, 7, 1, 8/)  ! Contents of convolution filter A

!  Initialize dummy arrays DX and DY.
   DX = 0.d0
   DY = 0.d0

!  Do the convolution
   DO I=0, 2*M+N
      DX = EOSHIFT(DX, SHIFT=-1, DIM=1)
      DY = EOSHIFT(DY, SHIFT=1, DIM=1)
      IF ((MOD(I,2)==0) .AND. (I<2*M)) THEN
         DX(1) = X(I/2 + 1)
      ENDIF
      DY = DY + A*DX
      IF (MOD(I,2)==0) THEN
         Y(I/2 + 1) = DY(1)
      ENDIF
   ENDDO


  return
end subroutine convolution