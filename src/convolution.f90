subroutine convolution(X,M,A,N,Y)
    ! in
    ! X input vector
    ! M size of X
    ! A convolution filter
    ! N size of convolution filter
    ! Y result, length M+N-1

    !uses either standard approach or fft method

    use kinds
    use settings, only: radar_convolution_fft
    use report_module
    implicit none

    INTEGER, intent(in) :: M  ! Size of input vector X
    INTEGER, intent(in) :: N   ! Size of convolution filter A

    REAL(kind=dbl), intent(in), DIMENSION(M) :: X
    REAL(kind=dbl), intent(in), DIMENSION(N) :: A
    REAL(kind=dbl), intent(out), DIMENSION(M+N-1) :: Y

    if (radar_convolution_fft) then
        if (verbose > 2) print*, "Entering FFT-Comvolution"
        call convolutionFFT(X,M,A,N,Y)
        if (verbose > 2) print*, "Done FFT-Comvolution"
    else
        if (verbose > 2) print*, "Entering Non-FFT-Comvolution"
        call convolution_slow(X,M,A,N,Y)
        if (verbose > 2) print*, "Done Non-FFT-Comvolution"
    end if
end subroutine convolution

subroutine convolution_slow(X,M,A,N,Y)


    ! convolution demonstration using a bidirectional full-systolic convolution array.
    !
    ! Author : Erik Urbach (0906204)
    !
    ! Course : Systolic Parallel Processing
    !
    ! Date   : Tue Nov  2 11:08:46 MET 1999
    !
    ! Note   : See figure 6.4 from "Systolic Parallel Processing" by N.Petkov.
    !

    use kinds
    implicit none
  
    INTEGER, intent(in) :: M  ! Size of input vector X
    INTEGER, intent(in) :: N   ! Size of convolution filter A
    INTEGER :: I

    REAL(kind=dbl), intent(in), DIMENSION(M) :: X
    REAL(kind=dbl), intent(in), DIMENSION(N) :: A
    REAL(kind=dbl), DIMENSION(N) :: DX, DY
    REAL(kind=dbl), intent(out), DIMENSION(M+N-1) :: Y

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
end subroutine convolution_slow

! convolutionFFT 
! requires FFTPACK

! subroutine convolution_fast(Xin,M,Ain,N,Yout)
! ! in
! ! X input vector
! ! M size of X
! ! A convolutionFFT filter
! ! N size of convolutionFFT filter
! ! Y result, length M+N-1
! 
! !FFT works best for arrays of length 2**x , this function extends the arrays to the needed length
! 
!   use kinds
!   implicit none
!    INTEGER, intent(in) :: M  ! Size of input vector X
!    INTEGER, intent(in) :: N   ! Size of convolutionFFT filter A
! 
!    REAL(kind=dbl), intent(in), DIMENSION(M) :: Xin
!    REAL(kind=dbl), intent(in), DIMENSION(N) :: Ain
!    REAL(kind=dbl), intent(out), DIMENSION(M+N-1) :: Yout
! 
!    INTEGER:: MN, MNext, Mext, Next
!    REAL(kind=dbl), dimension(:), allocatable :: Xext, Aext, Yext
! 
!   MN = M+N-1
!   MNext  = 2**CEILING(log(DBLE(MN))/log(2.d0))
! 
! 
!   allocate(Xext(MNext),Aext(MNext),Yext(MNext))
! 
!   Xext = 0.d0
!   Aext = 0.d0
!   Xext(1:M) = Xin
!   Aext(1:N) = Ain
! 
!   call convolutionFFT(Xext,Mext,Aext,Next,Yext)
! 
!   Yout(1:MN) = Yext
! 
!   deallocate(Xext,Aext,Yext)
! 
! end subroutine convolution_fast

subroutine convolutionFFT(Xin,M,Ain,N,Yout)
    ! in
    ! X input vector
    ! M size of X, must be a power of 2
    ! A convolution filter
    ! N size of convolution filter, must be a power of 2
    ! Y result, length M+N-1

    ! based on scipy/scipy/signal/signaltools.py
    ! and ffttest of pda:
    ! https://starlink.jach.hawaii.edu/svn/trunk/libraries/pda/Ffttest.f

    use kinds
    implicit none
  
    INTEGER, intent(in) :: M  ! Size of input vector X
    INTEGER, intent(in) :: N   ! Size of convolution filter A

    REAL(kind=dbl), intent(in), DIMENSION(M) :: Xin
    REAL(kind=dbl), intent(in), DIMENSION(N) :: Ain
    REAL(kind=dbl), intent(out), DIMENSION(M+N-1) :: Yout
    INTEGER :: MN, MNext
    INTEGER :: I
    REAL(kind=dbl) :: A, B, C, D
    REAL(kind=dbl),allocatable :: R1(:),R2(:),RF(:), &
    WSAVE(:)

    !increase input to same length
    MN = M+N-1
    !fft works best for power of 2 length
    MNext  = 2**CEILING(log(DBLE(MN))/log(2.d0))

    ! print*, M, N, MN, MNext

    allocate(R1(MNext),R2(MNext),RF(MNext), WSAVE(4*(MNext)+15))

    R1 = 0.d0
    R2 = 0.d0

    R1(1:M) = Xin(:)
    R2(1:N) = Ain(:)


    CALL DFFTI( MNext, WSAVE )
    CALL DFFTF( MNext, R1, WSAVE )
    CALL DFFTF( MNext, R2, WSAVE )


    !  Multiply the 2 transforms together. First multiply the zeroth term
    !  for which all imaginary parts are zero.
    RF( 1 ) = R1( 1 ) * R2( 1 )

    !  Now do the remaining terms. Real and imaginary terms are stored in
    !  adjacent elements of the arrays.
    DO I = 2, MNext - 1, 2

        A = R1( I )
        B = R1( I + 1 )
        C = R2( I )
        D = R2( I + 1 )

        RF( I ) = A*C - B*D
        RF( I + 1 ) = B*C + A*D

    END DO

    !  If there are an even number of elements, do the last term, for which
    !  the imaginary parts are again zero.
    IF( MOD( MNext, 2 ) .EQ. 0 ) RF( MNext ) = R1( MNext ) * R2( MNext )

    !  Now take the inverse FFT.
    CALL DFFTB( MNext, RF, WSAVE )


    !  Divide the results by MN to take account of the different
    !  normalisation of the FFTPACK results.
    RF = RF/( DBLE( MNext ))


    Yout(1:MN) = RF(:)
    deallocate(R1,R2,RF, WSAVE)
    return
end subroutine convolutionFFT
