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

    real(kind=dbl) :: refre, refim

    real(kind=dbl) :: dia1, dia2, ad, bd, alpha, gamma, b_snow, a_msnow

    real(kind=dbl), intent(out) :: salb, back

    complex(kind=dbl) :: mindex

  real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
  real(kind=dbl), dimension(4,4,nquad,2), intent(out) :: extinct_matrix
  real(kind=dbl), dimension(4,nquad,2), intent(out) :: emis_vector

    if (verbose .gt. 1) print*, 'Entering snow_ssp'

    call ref_ice(t, f, refre, refim)

    mindex = refre+Im*refim
    mindex = (  1.3404916923064127     , 1.06728917227096910E-003)
    mindex = (  1.0436438962756640     , 1.50318928653283049E-004)
    mindex = (  1.0658472188833250     , 2.27883690173356359E-004)
    if (.not. present(nc)) stop 'STOP in routine snow_ssp_tmat'
    call double_moments(swc,nc,gamma_snow(1),gamma_snow(2),gamma_snow(3),gamma_snow(4), &
      ad,bd,alpha,gamma,a_msnow,b_snow)
    nbins = 1
    dia1 = 2.5d-4
    dia2 = 1.d-2
    call tmatrix(f, mindex, dia1, dia2, nbins, &
          ad, bd, alpha, gamma, SD_snow,snow_density,swc,scatter_matrix,extinct_matrix, emis_vector)

    if (verbose .gt. 1) print*, 'Exiting snow_ssp_tmat'

    return

end subroutine snow_ssp_tmat
