! Subroutine for the setup of the parameters of the snow particle size distribution.
!
!
subroutine snow_ssp_tmat(f,swc,t,kext, back,scatter_matrix,extinct_matrix, emis_vector, nc)

    use kinds
    use nml_params, only: verbose, lphase_flag, n_0snowDsnow, EM_snow, &
    n_moments, isnow_n0, SD_snow
    use constants, only: pi, im, c
    use double_moments_module
    use conversions

    implicit none

    integer, parameter :: nstokes = 2, nquad = 16

    real(kind=dbl), intent(in) :: &
       swc,&
       t,&
       f

    real(kind=dbl), optional, intent(in) :: nc


    real(kind=dbl) :: ad, bd, alpha, gamma, b, a_m, k

    real(kind=dbl), intent(out) :: kext, back

    real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
    real(kind=dbl), dimension(nstokes,nstokes,nquad,2), intent(out) :: extinct_matrix
    real(kind=dbl), dimension(nstokes,nquad,2), intent(out) :: emis_vector


    !wavenumber
    k = 2* pi*f*1d9/c
    if (verbose .gt. 1) print*, 'Entering snow_ssp_tmat'

    if (.not. present(nc)) stop 'STOP in routine snow_ssp_tmat'
    call double_moments(swc,nc,gamma_snow(1),gamma_snow(2),gamma_snow(3),gamma_snow(4), &
      ad,bd,alpha,gamma,a_m,b)
    call tmatrix_snow(f, swc, t, nc, &
          ad, bd, alpha, gamma, a_m, b, SD_snow, scatter_matrix,extinct_matrix, emis_vector)



    !back is for NOT polarized radiation only, if you want to simulate a polarized Radar, use the full scattering matrix!
    back = scatter_matrix(1,16,1,16,2) !scatter_matrix(A,B;C;D;E) backscattering is M11 of Mueller or Scattering Matrix (A;C=1), in quadrature 2 (E) first 16 (B) is 180deg (upwelling), 2nd 16 (D) 0deg (downwelling). this definition is lokkiing from BELOW, scatter_matrix(1,16,1,16,3) would be from above!
    back = 4*pi*back!/k**2 !eq 4.82 Bohren&Huffman without k**2 (because of different definition of Mueller matrix according to Mishenko AO 2000). note that scatter_matrix contains already squard entries!

    kext = extinct_matrix(1,1,16,1) !11 of extinction matrix (=not polarized), at 0°, first quadrature. equal to extinct_matrix(1,1,16,2)
     
    if (verbose .gt. 1) print*, 'Exiting snow_ssp_tmat'

    return

end subroutine snow_ssp_tmat
