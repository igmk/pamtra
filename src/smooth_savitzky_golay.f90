SUBROUTINE SMOOTH_SAVITZKY_GOLAY(errorstatus,dataIn,length,dataOut)

!     Smooth data with a Savitzky-Golay filter.
!     The Savitzky-Golay filter removes high frequency noise from data.
!     It has the advantage of preserving the original shape and
!     features of the signal better than other types of filtering
!     approaches, such as moving averages techniques.
!
!     In the implementation here, the Savitzky-Golay filter is used 
!     with a 7 point window and second order!
!
!     Parameters
!     ----------
!     dataIn : array_like, shape (length)
!         the values of the time history of the signal.
!     length : int
!         the size of the data
!     Returns
!     -------
!     dataOut : ndarray, shape (N)
!         the smoothed signal 
!     Notes
!     -----
!     The Savitzky-Golay is a type of low-pass filter, particularly
!     suited for smoothing noisy data. The main idea behind this
!     approach is to make for each point a least-square fit with a
!     polynomial of high order over a odd-sized window centered at
!     the point.
!     ----------
!     .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
!        Data by Simplified Least Squares Procedures. Analytical
!        Chemistry, 1964, 36 (8), pp 1627-1639.
!     .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
!        W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
!        Cambridge University Press ISBN-13: 9780521880688
! 

  use kinds
  use report_module
  
  implicit none

  integer, parameter :: window = 7

  integer, intent(in) :: length
  real(kind=dbl), intent(in), dimension(length) :: dataIn
  real(kind=dbl), intent(out), dimension(length) :: dataOut

  real(kind=dbl), dimension(window) :: m
  real(kind=dbl), dimension(2*window+length-2) :: dataTmp
  real(kind=dbl), dimension(window+length-1) :: dataExt
  integer :: half_window

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'SMOOTH_SAVITZKY_GOLAY' 
  
  interface
    subroutine convolution(errorstatus,X,M,A,N,Y)
      use kinds
      implicit none
      integer(kind=long), intent(out) :: errorstatus
      INTEGER, intent(in) :: M  ! Size of input vector X
      INTEGER, intent(in) :: N   ! Size of convolution filter A
      REAL(kind=dbl), intent(in), DIMENSION(M) :: X
      REAL(kind=dbl), intent(in), DIMENSION(N) :: A
      REAL(kind=dbl), intent(out), DIMENSION(M+N-1) :: Y
    end subroutine convolution
  end interface

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
  
! coefficients, gained from http://www.scipy.org/Cookbook/SavitzkyGolay
  m = (/ -0.0952381d0 ,  0.14285714d0,  0.28571429d0,  0.33333333d0, &
        0.28571429d0, 0.14285714d0, -0.0952381d0 /)
  half_window = (window-1)/2 

! to deal with the borders, interpolate additional values
  dataExt((window+1)/2:length+(window-1)/2) = dataIn
  dataExt(1:half_window) = dataIn(1) - ABS(dataIn(half_window+1:2:-1) - dataIn(1))
  dataExt(length+half_window+1:window+length-1) = dataIn(length) &
      +ABS(dataIn(length-1:length - half_window:-1) - dataIn(length))

! convolve data with precalculated coefficients
  call convolution(err,dataExt,window+length-1,m,window,dataTmp)
  if (err /= 0) then
      msg = 'error in convolution!'
      call report(err, msg, nameOfRoutine)
      errorstatus = err
      return
  end if   
  
! trim result to original length
  dataOut(:) = dataTmp(window:window+length-1)

  errorstatus = err
  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return
  
END SUBROUTINE SMOOTH_SAVITZKY_GOLAY
