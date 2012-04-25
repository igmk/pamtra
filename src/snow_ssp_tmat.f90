! Subroutine for the setup of the parameters of the snow particle size distribution.
!
!
subroutine snow_ssp_tmat(f,swc,t,salb, back,scatter_matrix,extinct_matrix, emis_vector, nc)

    use kinds
    use nml_params, only: verbose, lphase_flag, n_0snowDsnow, EM_snow, &
    n_moments, isnow_n0, SD_snow, snow_density
    use constants, only: pi, im
    use double_moments_module
    use conversions

    implicit none

    integer :: nbins

    integer, parameter :: nstokes = 2, nquad = 16

    real(kind=dbl), intent(in) :: &
       swc,&
       t,&
       f

    real(kind=dbl), optional, intent(in) :: nc


    real(kind=dbl) :: dia1, dia2, ad, bd, alpha, gamma, b_snow, a_msnow

    real(kind=dbl), intent(out) :: salb, back

    real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
    real(kind=dbl), dimension(4,4,nquad,2), intent(out) :: extinct_matrix
    real(kind=dbl), dimension(4,nquad,2), intent(out) :: emis_vector

    if (verbose .gt. 1) print*, 'Entering snow_ssp'

    if (.not. present(nc)) stop 'STOP in routine snow_ssp_tmat'
    call double_moments(swc,nc,gamma_snow(1),gamma_snow(2),gamma_snow(3),gamma_snow(4), &
      ad,bd,alpha,gamma,a_msnow,b_snow)
!    dia1 = 1.d-4
    dia1 = 1.d-4
    dia2 = 2.d-2
    nbins = 30
!    dia1 = 150.0d-5
!    dia2 = 10000.0d-6
!    nbins = (dia2-dia1)/250.d-6+2
    call tmatrix(f, swc, t, nc, dia1, dia2, nbins, &
          ad, bd, alpha, gamma, a_msnow,b_snow,SD_snow, scatter_matrix,extinct_matrix, emis_vector)

    if (verbose .gt. 1) print*, 'Exiting snow_ssp_tmat'

    return

end subroutine snow_ssp_tmat
