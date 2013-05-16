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

subroutine radar_hildebrand_sekhon(spectrum,n_ave,n_ffts,noise_mean)

! written by P. Kollias, tranlated to Fortran by M. Maahn (12.2012)
! 
! in
! spectrum: radar sectrum [mm⁶/m³]
! n_ave: no of averaged spectra [-]
! n_ffts: no fft points [-]
!
! out
! noise_mean: neam noise leveö according to Hildebrand [mm⁶/m³]


  use kinds
!  use settings, only: verbose
        use report_module
  implicit none


  integer, intent(in) :: n_ave, n_ffts
  real(kind=dbl), dimension(n_ffts), intent(in) :: spectrum
  real(kind=dbl), intent(out) :: noise_mean

  real(kind=dbl), dimension(n_ffts) :: dummy, a1, a3,spectrum_sorted
  real(kind=dbl) :: sumLi, sumSq, sumNs, maxNs
  integer :: n, i, numNs


  interface
    SUBROUTINE DSORT (DX, DY, N, KFLAG)
      use kinds
      implicit none
      real(kind=dbl), dimension(N), intent(inout) :: DX, DY
      integer, intent(in) :: N, KFLAG
    END SUBROUTINE DSORT
  end interface

  spectrum_sorted = spectrum
  dummy=0.d0
  call dsort(spectrum_sorted, dummy, n_ffts, 1)

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
    a1(i) = SQRT(DBLE(n_ave))*(n*sumSq-a3(i))
    if (n > n_ffts/4) then
      if (a1(i) <= a3(i)) then
        sumNs = sumLi
        numNs = n
        maxNs = spectrum_sorted(i)
      end if
    else
      sumNs = sumLi
      numNs = n
      maxNs = spectrum_sorted(i)
    end if
  end do

  noise_mean   = sumNs/numNs
!   N_max    = maxNs
!   N_points = numNs


end subroutine radar_hildebrand_sekhon
