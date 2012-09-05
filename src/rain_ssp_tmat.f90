! Subroutine for the setup of the parameters of the snow particle size distribution.
!
!
subroutine rain_ssp_tmat(f,rwc,cwc,t,kext, back,scatter_matrix,extinct_matrix, emis_vector, nc)

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

    real(kind=dbl), intent(out) :: kext, back

    real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
    real(kind=dbl), dimension(nstokes,nstokes,nquad,2), intent(out) :: extinct_matrix
    real(kind=dbl), dimension(nstokes,nquad,2), intent(out) :: emis_vector

    if (verbose .gt. 1) print*, 'Entering rain_ssp_tmat'

    if (.not. present(nc)) stop 'STOP in routine rain_ssp_tmat'
    call double_moments(rwc,nc,gamma_rain(1),gamma_rain(2),gamma_rain(3),gamma_rain(4), &
      ad,bd,alpha,gamma,a_m,b,cwc)
    call tmatrix_rain(f, rwc, t, nc, &
          ad, bd, alpha, gamma, a_m, b, SD_rain, scatter_matrix,extinct_matrix, emis_vector)

    !back is for NOT polarized radiation only, if you want to simulate a polarized Radar, use the full scattering matrix!
    back = scatter_matrix(1,16,1,16,2) !scatter_matrix(A,B;C;D;E) backscattering is M11 of Mueller or Scattering Matrix (A;C=1), in quadrature 2 (E) first 16 (B) is 180deg (upwelling), 2nd 16 (D) 0deg (downwelling). this definition is lokkiing from BELOW, scatter_matrix(1,16,1,16,3) would be from above!
    back = 4*pi*back!/k**2 !eq 4.82 Bohren&Huffman without k**2 (because of different definition of Mueller matrix according to Mishenko AO 2000). note that scatter_matrix contains already squard entries!

    kext = extinct_matrix(1,1,16,1) !11 of extinction matrix (=not polarized), at 0°, first quadrature. equal to extinct_matrix(1,1,16,2)

    if (verbose .gt. 1) print*, 'Exiting rain_ssp_tmat'

    return

end subroutine rain_ssp_tmat
