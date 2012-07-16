! Subroutine for the setup of the parameters of the snow particle size distribution.
!
!
subroutine rain_ssp_tmat(f,rwc,cwc,t,salb, back,scatter_matrix,extinct_matrix, emis_vector, nc)

    use kinds
    use nml_params, only: verbose, lphase_flag, SD_rain
    use constants, only: pi, im
    use double_moments_module
    use conversions

    implicit none

    integer, parameter :: nstokes = 2, nquad = 16

    real(kind=dbl), intent(in) :: &
       rwc,&
       cwc,&
       t,&
       f

    real(kind=dbl), optional, intent(in) :: nc


    real(kind=dbl) :: ad, bd, alpha, gamma, b, a_m

    real(kind=dbl), intent(out) :: salb, back

    real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
    real(kind=dbl), dimension(4,4,nquad,2), intent(out) :: extinct_matrix
    real(kind=dbl), dimension(4,nquad,2), intent(out) :: emis_vector

    if (verbose .gt. 1) print*, 'Entering rain_ssp_tmat'

    if (.not. present(nc)) stop 'STOP in routine rain_ssp_tmat'
    call double_moments(rwc,nc,gamma_rain(1),gamma_rain(2),gamma_rain(3),gamma_rain(4), &
      ad,bd,alpha,gamma,a_m,b,cwc)
    call tmatrix_rain(f, rwc, t, nc, &
          ad, bd, alpha, gamma, a_m, b, SD_rain, scatter_matrix,extinct_matrix, emis_vector)
    salb = 0.d0
    back = 0.d0
    if (verbose .gt. 1) print*, 'Exiting rain_ssp_tmat'

    return

end subroutine rain_ssp_tmat
