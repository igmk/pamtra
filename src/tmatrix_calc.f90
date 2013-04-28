

subroutine tmatrix_calc(quad,qua_num,frequency,wave_num,ref_index,axi, nstokes,&
    as_ratio, alpha, beta, azimuth_num, azimuth0_num,&
    scatter_matrix,extinct_matrix,emis_vector)

!DB VERSION

    !  calculate the matrix and vectors, for a single particle with a single orientation
    !
    !   input:
    !       quad            char    name of quadrature
    !       quad_num        int     number of quadrature angles
    !       frequency       double  frequency [Hz]
    !       wave_num        double  wave number [1/m]
    !       ref_index       complex refractive index
    !       axi             double  equivalent sphere radius [m]
    !       nstokes         int     number of stokes components
    !       as_ratio        double  aspect ratio
    !       alpha           double  orientation of the particle [°]
    !       beta            double  orientation of the particle [°]
    !       azimuth_num     int     number of azimuth angles
    !       azimuth0_num    int     number of azimuth angles
    !
    !   output:
    !       scatter_matrix  double  scattering matrix []
    !       extinct_matrix  double  extinction matrix []
    !       emis_vector     double  emission vector []

    use kinds
    use rt_utilities, only: lobatto_quadrature

    implicit none


    character(1), intent(in) :: quad
    integer, intent(in) :: qua_num
    real(kind=dbl), intent(in) :: frequency
    real(kind=dbl), intent(in) :: wave_num
    complex(kind=ext) :: ref_index
    real(kind=dbl), intent(in) :: axi
    integer, intent(in) :: nstokes
    real(kind=dbl), intent(in) :: as_ratio
    real(kind=dbl), intent(in) :: alpha
    real(kind=dbl), intent(in) :: beta
    integer, intent(in) :: azimuth_num
    integer, intent(in) :: azimuth0_num
    real(kind=dbl), intent(out), dimension(nstokes,qua_num,nstokes,qua_num,2) :: scatter_matrix
    real(kind=dbl), intent(out), dimension(nstokes,nstokes,qua_num) :: extinct_matrix
    real(kind=dbl), intent(out), dimension(nstokes,qua_num) :: emis_vector


    real(kind=dbl) :: maximum_size, minimum_size, &
        num_0, lambda_0, temperature, &
        particle_size, particle_mass, pi
    real(kind=dbl) ::thet0, thet, phi, phi0,&
        rat, sum_tmp, phi_weights, phi0_weights
    integer :: i, j, k, l, m, n, &
       particle_num,&
       ii,jj,kk,ll,kkk,kkk1,kkk2

    !       parameter(qua_num=8, nstokes = 4)


    complex(kind=dbl) :: s11, s12, s21, s22
    real(kind=dbl) :: qua_angle(qua_num), qua_weights(qua_num)
    real(kind=dbl) :: &
        scatt_matrix_tmp1_11, scatt_matrix_tmp1_12,&
        scatt_matrix_tmp1_21, scatt_matrix_tmp1_22
    real(kind=dbl) ::&
        emis_vector_tmp1_11(2*qua_num), emis_vector_tmp1_12(2*qua_num),&
        emis_vector_tmp2_11,emis_vector_tmp2_12
    real(kind=dbl) :: thet0_weights, thet_weights


    integer :: nmax, np

    real(kind=ext) :: mrr, mri, lam
    ! some factors that stay constant during calculations
    real(kind=dbl) :: fact_sca
    complex(kind=dbl) :: fact_ext


    np = -1
    rat = 1.e0

    LAM = 3.0e8_ext/(frequency)!*1e6
!     axi = axi!*1.e6
    mrr = REAL(ref_index)
    mri = abs(IMAG(ref_index))


!      RAT=.963711
!      LAM=6.283185   
!      MRR=1.5D0
!      MRI=0.02D0
!      NP=-1

! !ICHOICE=1  NCHECK=1
! RAT= .963711
! !PROLATE SPHEROIDS, A/B=   .5000000
! 
! LAM=  6.283185   
! MRR= .1500D+01   
! MRI= .2000D-01


!ACCURACY OF COMPUTATIONS DDELT =  .10D-02
!EQUAL-SURFACE-AREA-SPHERE RADIUS= 10.0000
!thet0= 56.00  thet= 65.00  phi0=114.00  phi=128.00  alpha=145.00  beta= 52.00
!AMPLITUDE MATRIX
!S11=-.50941D+01 + i* .24402D+02
!S12=-.19425D+01 + i* .19971D+01
!S21=-.11521D+01 + i*-.30977D+01
!S22=-.69323D+01 + i* .24748D+02
!PHASE MATRIX
!  650.3172  -17.9846   10.0498  -12.7580
!  -21.1462  631.6323 -127.3059   87.2144
!    6.8321  132.6132  635.2768  -34.7730
!   -9.6629  -78.1229   51.4094  643.1738
! time =     .35 min
!


! call the tmatrix routine amplq -> fills common block /TMAT/
    call tmatrix_amplq(lam, mrr,mri, AXI, AS_RATIO, RAT, NP,nmax)

    extinct_matrix = 0.e0
    scatter_matrix = 0.e0
    emis_vector = 0.e0
    emis_vector_tmp2_11 = 0.e0
    emis_vector_tmp2_12 = 0.e0
    ! if the particle is rotationally-symmetric, reduce calculation time for orientation-averaging
    ! if not, do orientation averaging for incident and scatterred directions

    !      write(*,*)wave_num
    pi = acos(-1.d0)

    fact_sca = 0.5e0/(wave_num**2)
    fact_ext = 2.0e0*pi*cmplx(0.,1.)/wave_num**2.e0
    ! calculate the quadrature angle, number and weight according to quadrature method
    ! subroutine lobatto_quadrature and gauss_legendre_quadrature in the file named 'refractive_index.f'
    if (quad(1:1).eq.'l'.or.quad(1:1).eq.'L') then
        call lobatto_quadrature(qua_num,qua_angle,qua_weights)
    endif

    ! for each quadrature angle
    ii = 1
    do 1241 jj = 1, qua_num
        thet0=acos(qua_angle(jj)*(-1.)**(real(ii)-1))*180.d0/pi
        thet0_weights = qua_weights(jj)
        if(thet0.gt.179.9999)thet0=180.0d0
        ! initializing the emis vector summation
        emis_vector_tmp1_11 = 0.d0
        emis_vector_tmp1_12 = 0.d0

        do 1242 kk = 1, 2
            kkk1 = kk
            do 1243 ll = 1, qua_num
                thet=acos(qua_angle(ll)*(-1.)**(real(kk)-1))*180.d0/pi
                thet_weights=qua_weights(ll)
                if(thet.gt.179.9999)thet=180.d0

                do 1244 m = 1, azimuth0_num ! 1
                    phi0 = 360.0d0/(real(azimuth0_num))*(real(m)-1.d0)
                    phi0_weights = 1.d0/360.d0*(360.d0/azimuth0_num)
                    !		     if(azimuth0_num.eq.1)phi0 = 0.0

                    scatt_matrix_tmp1_11 = 0.d0
                    scatt_matrix_tmp1_12 = 0.d0

                    scatt_matrix_tmp1_21 = 0.d0
                    scatt_matrix_tmp1_22 = 0.d0

                    do 1245 n = 1, azimuth_num ! 30
                        phi = 360.d0/real(azimuth_num)*(real(n)-1.d0)
                        phi_weights = 1.d0/360.d0*(360.d0/azimuth_num)



		CALL tmatrix_AMPL(NMAX,dble(LAM),THET0,THET,PHI0,PHI,ALPHA,BETA,&
		  S11,S12,S21,S22)


                        s11 = s11*wave_num
                        s12 = s12*wave_num
                        s21 = s21*wave_num
                        s22 = s22*wave_num

                        scatt_matrix_tmp1_11 = scatt_matrix_tmp1_11 + (fact_sca*&
                            (s11*dconjg(s11)+s12*dconjg(s12)+s21*dconjg(s21)+s22*dconjg(s22)))*phi_weights

                        scatt_matrix_tmp1_12 = scatt_matrix_tmp1_12 + (fact_sca*&
                            (s11*dconjg(s11)-s12*dconjg(s12)+s21*dconjg(s21)-s22*dconjg(s22)))*phi_weights

                        scatt_matrix_tmp1_21 = scatt_matrix_tmp1_21 + (fact_sca*&
                            (s11*dconjg(s11)+s12*dconjg(s12)-s21*dconjg(s21)-s22*dconjg(s22)))*phi_weights

                        scatt_matrix_tmp1_22 = scatt_matrix_tmp1_22 + (fact_sca*&
                            (s11*dconjg(s11)-s12*dconjg(s12)-s21*dconjg(s21)+s22*dconjg(s22)))*phi_weights

                        if (phi0 .eq. phi .and. thet0 .eq. thet) then
                            extinct_matrix(1,1,jj) = extinct_matrix(1,1,jj)+phi0_weights*(-real((s11 + s22)*fact_ext))
                            extinct_matrix(1,2,jj) = extinct_matrix(1,2,jj)+phi0_weights*(-real((s11 - s22)*fact_ext))
                            extinct_matrix(2,1,jj) = extinct_matrix(2,1,jj)+phi0_weights*(-real((s11 - s22)*fact_ext))
                            extinct_matrix(2,2,jj) = extinct_matrix(2,2,jj)+phi0_weights*(-real((s11 + s22)*fact_ext))
                        end if
1245                    continue   ! phi

                    scatter_matrix(1,ll,1,jj,kkk1) = scatter_matrix(1,ll,1,jj,kkk1) + scatt_matrix_tmp1_11*phi0_weights
                    scatter_matrix(1,ll,2,jj,kkk1) = scatter_matrix(1,ll,2,jj,kkk1) + scatt_matrix_tmp1_12*phi0_weights
                    scatter_matrix(2,ll,1,jj,kkk1) = scatter_matrix(2,ll,1,jj,kkk1) + scatt_matrix_tmp1_21*phi0_weights
                    scatter_matrix(2,ll,2,jj,kkk1) = scatter_matrix(2,ll,2,jj,kkk1) + scatt_matrix_tmp1_22*phi0_weights

1244                continue  ! phi0
                ! calculate the summation of the scattering matrix in the whole sphere
                emis_vector_tmp1_11(ll+(kk-1)*qua_num) = scatter_matrix(1,ll,1,jj,kkk1)*thet_weights*2.*pi
                emis_vector_tmp1_12(ll+(kk-1)*qua_num) = scatter_matrix(1,ll,2,jj,kkk1)*thet_weights*2.*pi

1243            continue ! thet ll
1242        continue

        emis_vector(1,jj) = extinct_matrix(1,1,jj) - sum(emis_vector_tmp1_11)
        emis_vector(2,jj) = extinct_matrix(1,2,jj) - sum(emis_vector_tmp1_12)

1241    continue ! thet0 jj



     return
 end subroutine tmatrix_calc

