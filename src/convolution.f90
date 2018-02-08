subroutine convolution(errorstatus,X,M,A,N,Y)
    ! convolve X with filter A
    ! uses either standard approach or fft method

    use kinds
    use settings, only: radar_convolution_fft
    use report_module
    implicit none

    integer, intent(in) :: M  ! Size of input vector X
    integer, intent(in) :: N ! Size of convolution filter A
    real(kind=dbl), intent(in), dimension(M) :: X ! X input vector
    real(kind=dbl), intent(in), dimension(N) :: A ! A convolution filter
    real(kind=dbl), intent(out), dimension(M+N-1) :: Y ! Y result, length M+N-1

    integer(kind=long), intent(out) :: errorstatus ! error reported to report module
    integer(kind=long) :: err
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'convolution' 
    
    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    err = 0

    call assert_false(err,(ALL(X == 0)),&
      "all x values zero")
    call assert_false(err,(ALL(A == 0)),&
      "all A values zero")
    call assert_false(err,(ANY(ISNAN(X))),&
      "found nan in x")
    call assert_false(err,(ANY(ISNAN(A))),&
      "found nan in a")
    if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
    end if


    if (radar_convolution_fft) then
        if (verbose > 2) print*, "Entering FFT-Convolution"
        call convolutionFFT(X,M,A,N,Y)
        if (verbose > 2) print*, "Done FFT-Convolution"
    else
        if (verbose > 2) print*, "Entering Non-FFT-Convolution"
        call convolution_slow(X,M,A,N,Y)
        if (verbose > 2) print*, "Done Non-FFT-Convolution"
    end if
    
    call assert_false(err,(ALL(X == Y)),&
      "all Y values zero")
    call assert_false(err,(ANY(ISNAN(Y))),&
      "found nan in Y")
    if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
    end if


    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    include 'fftw3.f' ! To be included in the same directory of this file. 
                      ! Fortran compiler does not look into /usr/include/
                      ! Copy it from /usr/include/ or provide -I $FFTW3_INC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    INTEGER, intent(in) :: M  ! Size of input vector X
    INTEGER, intent(in) :: N   ! Size of convolution filter A

    REAL(kind=dbl), intent(in), DIMENSION(M) :: Xin
    REAL(kind=dbl), intent(in), DIMENSION(N) :: Ain
    REAL(kind=dbl), intent(out), DIMENSION(M+N-1) :: Yout
    INTEGER :: MN, MNext
!    INTEGER :: I
    REAL(kind=dbl),allocatable :: R1(:),R2(:),RF(:)!, &

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    double complex, allocatable :: R1F(:), R2F(:), RFF(:) ! intermidiate stage
    integer*8 plan                                    ! We always need a plan!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !increase input to same length
    MN = M+N-1
    !fft works best for power of 2 length
    MNext  = 2**CEILING(log(DBLE(MN))/log(2.d0))

    allocate(R1(MNext),R2(MNext),RF(MNext))
    allocate(R1F(MNext/2+1),R2F(MNext/2+1),RFF(MNext/2+1))

    R1 = 0.d0
    R2 = 0.d0

    R1(1:M) = Xin(:)
    R2(1:N) = Ain(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call dfftw_plan_dft_r2c_1d(plan, MNext, R1, R1F, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan,R1,R1F)

    call dfftw_plan_dft_r2c_1d(plan, MNext, R2, R2F, FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan, R2, R2F)

    RFF = R1F*R2F ! complex vector arithmetics is cool and super efficient  !!
    call dfftw_plan_dft_c2r_1d(plan,MNext,RFF,RF,FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan, RFF, RF)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! I keep this version as a comment since it involves halfcomplex formatted
! vectors which imply roughly half of memory occupancy since only three vectors
! are needed in this implementation and of course three less allocations
! but unfortunately, due to the complex arithmetics of the inner loop to compute
! the transform convolution the efficiency of the algorithm is lower than the r2c
! current implementation (which internally uses halfcomplex values... ).
! The vector-vector multiplication is what makes the current implementation so
! efficient, maybe, future compiler optimizations or a better halfcomplex
! multiplication algorithm will bring this back to production

!    call dfftw_plan_r2r_1d(plan, MNext, R1, R1, FFTW_R2HC, FFTW_ESTIMATE)
!    call dfftw_execute_r2r(plan,R1,R1)
!    call dfftw_plan_r2r_1d(plan, MNext, R2, R2, FFTW_R2HC, FFTW_ESTIMATE)
!    call dfftw_execute_r2r(plan,R2,R2)

!    RF(1) = R1(1)*R2(1)
!    RF(MNext/2+1) = R1(MNext/2+1)*R2(MNext/2+1)

!    do I = 2,MNext/2
!        RF(I) = R1(I)*R2(I) - R1(MNext+2-I)*R2(MNext+2-I)
!        RF(MNext+2-I) = R1(I)*R2(MNext+2-I) + R1(MNext+2-I)*R2(I)
!    enddo
!    call dfftw_plan_r2r_1d(plan, MNext, RF, RF, FFTW_HC2R, FFTW_ESTIMATE)
!    call dfftw_execute_r2r(plan,RF,RF)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Divide the results by MN to normalize
    RF = RF/( DBLE( MNext ))

    Yout(:) = RF(1:MN)
    deallocate(R1,R2,RF)
    deallocate(R1F, R2F, RFF)
    
    return
end subroutine convolutionFFT
