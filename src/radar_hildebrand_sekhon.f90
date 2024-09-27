! 
!  Description:  
!  After sorting the spectral points in order of ascending 
!  power, the algorithm goes through all values of n from 0 to numPts-1 
!  and for each n computes the values of both sides of the above equation.
!  For early values of n, the left side is generally less than the 
!  right and then becomes greater than the right side as higher-valued 
!  points that do not belong to the original distribution are included
!  in the n points. The last crossing from less than to greater than
!  is considered the division between the desired set of points and
!  "outliers" (e.g. noise vs. signal noise or signal vs. interference,
!  etc.).
! 
!  Hildebrand, P. H., and R. S. Sekhon, Objective determination of
!  the noise level in Doppler spectra, J. Appl. Meteorol., 13, 808, 1974.
!  
!  sum(S^2)/N - (sum(S)/N)^2 = (mean(S))^2 / navg
!  navg*[sum(S^2)/N - (sum(S))^2 / N^2] = (sum(S))^2 / N^2
!  navg*[N*sum(S^2) - sum(S) * sum(S)] = sum(S) * sum(S)
! 
! 

subroutine radar_hildebrand_sekhon(&
  errorstatus,&
  spectrum,&
  n_ave,&
  n_ffts,&
  noise_mean,&
  noise_max)

! written by P. Kollias, tranlated to Fortran by M. Maahn (12.2012)
! 
! in
! spectrum: radar sectrum [mm6/m3]
! n_ave: no of averaged spectra [-]
! n_ffts: no fft points [-]
!
! out
! noise_mean: neam noise level according to Hildebrand [mm6/m3]


  use kinds
!  use settings, only: verbose
        use report_module
  implicit none


  integer, intent(in) :: n_ave, n_ffts
  real(kind=dbl), dimension(n_ffts), intent(in) :: spectrum
  real(kind=dbl), intent(out) :: noise_mean
  real(kind=dbl), intent(out) :: noise_max

  real(kind=dbl), dimension(n_ffts) :: dummy, a1, a3,spectrum_sorted
  real(kind=dbl) :: sumLi, sumSq, sumNs, maxNs
  integer :: n, i, numNs

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'radar_hildebrand_sekhon' 
      
  interface
    SUBROUTINE DSORT (errorstatus,DX, DY, N, KFLAG)
      use kinds
      implicit none
      integer(kind=long), intent(out) :: errorstatus
      real(kind=dbl), dimension(N), intent(inout) :: DX, DY
      integer, intent(in) :: N, KFLAG
    END SUBROUTINE DSORT
  end interface

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
  
  spectrum_sorted = spectrum
  dummy=0.d0
  call dsort(err,spectrum_sorted, dummy, n_ffts, 1)
  if (err /= 0) then
      msg = 'error in dsort!'
      call report(err, msg, nameOfRoutine)
      errorstatus = err
      stop !return
  end if 
  sumLi = 0.d0
  sumSq = 0.d0
  n     = 0  
  a1 = 0.d0
  a3 = 0.d0

  do i= 1, n_ffts
    sumLi = sumLi + spectrum_sorted(i)
    sumSq = sumSq + spectrum_sorted(i)**2
    n = n+1
    a3(i) = sumLi*sumLi
    a1(i) = DBLE(n_ave)*(n*sumSq-a3(i))
      if (a1(i) <= a3(i)) then
        sumNs = sumLi
        numNs = n
        maxNs = spectrum_sorted(i)
      else
        !partial spectrum no longer has characteristics of white noise
        EXIT
      end if

  end do

  noise_mean   = sumNs/numNs
  noise_max    = maxNs
!   N_points = numNs

  if (verbose >= 5) print*, "hildebrand found",  noise_mean, noise_max


  errorstatus = err
  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return

end subroutine radar_hildebrand_sekhon
