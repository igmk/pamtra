subroutine matrix_cal(quad,qua_num,frequency,wave_num,snow_ref,axi, nstokes,&
as_ratio, alpha, beta, azimuth_num, azimuth0_num,&
scatter_matrix,extinct_matrix,emis_vector)
    !  calculate the matrix and vectors, for a single particle with a single orientation
    !       program matrix_cal
    !
    !   input:
    !       quad            char    name of quadrature
    !       quad_num        int     number of quadrature angles
    !       frequency       double  frequency [Hz]
    !       wave_num        double  wave number [1/m]
    !       snow_ref        complex refractive index
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

    use kinds, only: dbl
    use rt_utilities, only: gauss_legendre_quadrature,&
    lobatto_quadrature
    implicit none

    real*8 frequency, &
    as_ratio, &
    pi
    real*8  wave_num, axi, thet0, thet, phi, phi0,&
    alpha, beta, rat, phi_weights, phi0_weights
    integer m, n, azimuth_num, &
    azimuth0_num, qua_num,&
    ii,jj,kk,ll,nstokes,kkk1,kkk2

    !       parameter(qua_num=8, nstokes = 4)

    character*1 quad
    complex*16 snow_ref, s11, s12, s21, s22
    real*8 qua_angle(qua_num), qua_weights(qua_num)
    real*8 &
    scatt_matrix_tmp1_11, scatt_matrix_tmp1_12,&! scatt_matrix_tmp1_13, scatt_matrix_tmp1_14,&
    scatt_matrix_tmp1_21, scatt_matrix_tmp1_22 ! scatt_matrix_tmp1_23, scatt_matrix_tmp1_24,&
    !        scatt_matrix_tmp1_31, scatt_matrix_tmp1_32, scatt_matrix_tmp1_33, scatt_matrix_tmp1_34,&
    !        scatt_matrix_tmp1_41, scatt_matrix_tmp1_42, scatt_matrix_tmp1_43, scatt_matrix_tmp1_44
    real*8&
    emis_vector_tmp1_11(2*qua_num), emis_vector_tmp1_12(2*qua_num),&
    emis_vector_tmp2_11,emis_vector_tmp2_12 !,&
    !        emis_vector_tmp1_13(2*qua_num), emis_vector_tmp1_14(2*qua_num),&
    !        emis_vector_tmp2_13,emis_vector_tmp2_14
    real*8 thet0_weights, thet_weights
    real*8  &
    scatter_matrix(nstokes,qua_num,nstokes,qua_num,4),&
    extinct_matrix(nstokes,nstokes,qua_num,2),&
    emis_vector(nstokes,qua_num,2)

    ! some factors that stay constant during calculations
    real*8 :: fact_sca
    complex*8 :: fact_ext

print*, "IN"
print*, quad,qua_num,frequency,wave_num,snow_ref,axi, nstokes,&
as_ratio, alpha, beta, azimuth_num, azimuth0_num

    extinct_matrix = 0.d0
    scatter_matrix = 0.d0
    emis_vector = 0.d0
    emis_vector_tmp2_11 = 0.d0
    emis_vector_tmp2_12 = 0.d0
    !    emis_vector_tmp2_13 = 0.d0
    !    emis_vector_tmp2_14 = 0.d0
    ! if the particle is rotationally-symmetric, reduce calculation time for orientation-averaging
    ! if not, do orientation averaging for incident and scatterred directions

    !      write(*,*)wave_num
    pi = acos(-1.d0)

    fact_sca = 0.5d0/(wave_num**2)
    fact_ext = 2.0d0*pi*dcmplx(0,1)/wave_num**2.d0

    ! calculate the quadrature angle, number and weight according to quadrature method
    ! subroutine lobatto_quadrature and gauss_legendre_quadrature in the file named 'refractive_index.f'
    if ((quad(1:1).eq.'l').or.(quad(1:1).eq.'L')) then
        call lobatto_quadrature(qua_num,qua_angle(1:qua_num),&
        qua_weights(1:qua_num))
    endif
    if ((quad(1:1).eq.'g').or.(quad(1:1).eq.'g')) then
        call gauss_legendre_quadrature(qua_num,qua_angle(1:qua_num),&
        qua_weights(1:qua_num))
    endif
 
    !	 write(*,*)quad,qua_num,qua_angle

    ! for each quadrature angle
    ii = 1
    do 1241 jj = 1, qua_num
        thet0=acos(qua_angle(jj)*(-1.)**(real(ii)-1))*180.d0/pi
        thet0_weights = qua_weights(jj)
        if(thet0.gt.179.9999)thet0=180.0d0
        !	       write(*,*)thet0, qua_angle(jj),qua_weights(jj)
        ! initializing the emis vector summation
        emis_vector_tmp1_11 = 0.d0
        emis_vector_tmp1_12 = 0.d0
        !        emis_vector_tmp1_13 = 0.d0
        !        emis_vector_tmp1_14 = 0.d0

        do 1242 kk = 1, 2
            kkk1 = (kk-1)*2 + 1
            kkk2 = 4 - (kk-1)*2
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
                    !                    scatt_matrix_tmp1_13 = 0.d0
                    !                    scatt_matrix_tmp1_14 = 0.d0

                    scatt_matrix_tmp1_21 = 0.d0
                    scatt_matrix_tmp1_22 = 0.d0
                    !                    scatt_matrix_tmp1_23 = 0.d0
                    !                    scatt_matrix_tmp1_24 = 0.d0
                    !
                    !                    scatt_matrix_tmp1_31 = 0.d0
                    !                    scatt_matrix_tmp1_32 = 0.d0
                    !                    scatt_matrix_tmp1_33 = 0.d0
                    !                    scatt_matrix_tmp1_34 = 0.d0
                    !
                    !                    scatt_matrix_tmp1_41 = 0.d0
                    !                    scatt_matrix_tmp1_42 = 0.d0
                    !                    scatt_matrix_tmp1_43 = 0.d0
                    !                    scatt_matrix_tmp1_44 = 0.d0

                    do 1245 n = 1, azimuth_num ! 30
                        phi = 360.d0/real(azimuth_num)*(real(n)-1.d0)
                        phi_weights = 1.d0/360.d0*(360.d0/azimuth_num)

                        call tmatrix_cal(frequency, snow_ref, axi, as_ratio, alpha, beta, phi0, &
                        thet0, phi, thet, rat, s11, s12, s21, s22)

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

                        !                        scatt_matrix_tmp1_13 = scatt_matrix_tmp1_13 + (fact_sca*&
                        !                            real(s11*dconjg(s12)+s21*dconjg(s22)))*phi_weights
                        !
                        !                        scatt_matrix_tmp1_14 = scatt_matrix_tmp1_14 + (fact_sca*&
                        !                            (-imag(s11*dconjg(s12)+s21*dconjg(s22))))*phi_weights

                        !                            scatt_matrix_tmp1_23 = scatt_matrix_tmp1_23 + (fact_sca*&
                        !                                real(s11*dconjg(s12)-s21*dconjg(s22)))*phi_weights
                        !
                        !                            scatt_matrix_tmp1_24 = scatt_matrix_tmp1_24 + (fact_sca*&
                        !                                (-imag(s11*dconjg(s12)-s21*dconjg(s22))))*phi_weights
                        !
                        !                            scatt_matrix_tmp1_31 = scatt_matrix_tmp1_31 + (fact_sca*&
                        !                                real(s11*dconjg(s21)+s12*dconjg(s22)))*phi_weights
                        !
                        !                            scatt_matrix_tmp1_32 = scatt_matrix_tmp1_32 + (fact_sca*&
                        !                                real(s11*dconjg(s21)-s12*dconjg(s22)))*phi_weights
                        !
                        !                            scatt_matrix_tmp1_33 = scatt_matrix_tmp1_33 + (fact_sca*&
                        !                                real(s11*dconjg(s22)+s12*dconjg(s21)))*phi_weights
                        !
                        !                            scatt_matrix_tmp1_34 = scatt_matrix_tmp1_34 + (fact_sca*&
                        !                                (-imag(s11*dconjg(s22)-s12*dconjg(s21))))*phi_weights
                        !
                        !                            scatt_matrix_tmp1_41 = scatt_matrix_tmp1_41 + (fact_sca*&
                        !                                imag(s11*dconjg(s21)+s12*dconjg(s22)))*phi_weights
                        !
                        !                            scatt_matrix_tmp1_42 = scatt_matrix_tmp1_42 + (fact_sca*&
                        !                                imag(s11*dconjg(s21)-s12*dconjg(s22)))*phi_weights
                        !
                        !                            scatt_matrix_tmp1_43 = scatt_matrix_tmp1_43 + (fact_sca*&
                        !                                imag(s11*dconjg(s22)+s12*dconjg(s21)))*phi_weights
                        !
                        !                            scatt_matrix_tmp1_44 = scatt_matrix_tmp1_44 + (fact_sca*&
                        !                                real(s11*dconjg(s22)-s12*dconjg(s21)))*phi_weights

                        if ((phi0 .eq. phi) .and. (thet0 .eq. thet)) then
                            !			  write(*,*)phi0,phi,thet0,thet, wave_num
                            !			  write(*,*)'matrix_cal.f',s11,s22
                            extinct_matrix(1,1,jj,:) = extinct_matrix(1,1,jj,:)+phi0_weights*(-real((s11 + s22)*fact_ext))
                            extinct_matrix(1,2,jj,:) = extinct_matrix(1,2,jj,:)+phi0_weights*(-real((s11 - s22)*fact_ext))
                            extinct_matrix(2,1,jj,:) = extinct_matrix(2,1,jj,:)+phi0_weights*(-real((s11 - s22)*fact_ext))
                            extinct_matrix(2,2,jj,:) = extinct_matrix(2,2,jj,:)+phi0_weights*(-real((s11 + s22)*fact_ext))
                        !                            extinct_matrix(1,3,jj,:) = extinct_matrix(1,3,jj,:)+phi0_weights*(-real((s12 + s21)*fact_ext))
                        !                            extinct_matrix(1,4,jj,:) = extinct_matrix(1,4,jj,:)+phi0_weights*(-imag((s12 - s21)*fact_ext))
                        !                            extinct_matrix(2,3,jj,:) = extinct_matrix(2,3,jj,:)+phi0_weights*(-real((s12 - s21)*fact_ext))
                        !                            extinct_matrix(2,4,jj,:) = extinct_matrix(2,4,jj,:)+phi0_weights*(-imag((s12 + s21)*fact_ext))
                        !                            extinct_matrix(3,1,jj,:) = extinct_matrix(3,1,jj,:)+phi0_weights*(-real((s21 + s12)*fact_ext))
                        !                            extinct_matrix(3,2,jj,:) = extinct_matrix(3,2,jj,:)+phi0_weights*(-real((s21 - s12)*fact_ext))
                        !                            extinct_matrix(3,3,jj,:) = extinct_matrix(3,3,jj,:)+phi0_weights*(-real((s11 + s22)*fact_ext))
                        !                            extinct_matrix(3,4,jj,:) = extinct_matrix(3,4,jj,:)+phi0_weights*(-imag((s22 - s11)*fact_ext))
                        !                            extinct_matrix(4,1,jj,:) = extinct_matrix(4,1,jj,:)+phi0_weights*(-imag((s21 - s12)*fact_ext))
                        !                            extinct_matrix(4,2,jj,:) = extinct_matrix(4,2,jj,:)+phi0_weights*(imag((s12 - s21)*fact_ext))
                        !                            extinct_matrix(4,3,jj,:) = extinct_matrix(4,3,jj,:)+phi0_weights*(-imag((s11 - s22)*fact_ext))
                        !                            extinct_matrix(4,4,jj,:) = extinct_matrix(4,4,jj,:)+phi0_weights*(-real((s11 + s22)*fact_ext))
                        end if
                    !                    write(1234,*)thet0, phi0, thet,phi
                    !                    write(1234,*)thet0, phi0, thet,phi
                    !                    write(1234,*)s11, s12
                    !                    write(1234,*)s21, s22
1245                continue   ! phi

                    scatter_matrix(1,ll,1,jj,kkk1) = scatter_matrix(1,ll,1,jj,kkk1) + scatt_matrix_tmp1_11*phi0_weights
                    scatter_matrix(1,ll,2,jj,kkk1) = scatter_matrix(1,ll,2,jj,kkk1) + scatt_matrix_tmp1_12*phi0_weights
                    scatter_matrix(2,ll,1,jj,kkk1) = scatter_matrix(2,ll,1,jj,kkk1) + scatt_matrix_tmp1_21*phi0_weights
                    scatter_matrix(2,ll,2,jj,kkk1) = scatter_matrix(2,ll,2,jj,kkk1) + scatt_matrix_tmp1_22*phi0_weights

                    scatter_matrix(1,ll,1,jj,kkk2) = scatter_matrix(1,ll,1,jj,kkk1)
                    scatter_matrix(1,ll,2,jj,kkk2) = scatter_matrix(1,ll,2,jj,kkk1)
                    scatter_matrix(2,ll,1,jj,kkk2) = scatter_matrix(2,ll,1,jj,kkk1)
                    scatter_matrix(2,ll,2,jj,kkk2) = scatter_matrix(2,ll,2,jj,kkk1)

                !                        scatter_matrix(1,ll,3,jj,kkk) = scatter_matrix(1,ll,3,jj,kkk) + scatt_matrix_tmp1_13/azimuth0_num
                !                        scatter_matrix(1,ll,4,jj,kkk) = scatter_matrix(1,ll,4,jj,kkk) + scatt_matrix_tmp1_14/azimuth0_num
                !                        scatter_matrix(2,ll,3,jj,kkk) = scatter_matrix(2,ll,3,jj,kkk) + scatt_matrix_tmp1_23/azimuth0_num
                !                        scatter_matrix(2,ll,4,jj,kkk) = scatter_matrix(2,ll,4,jj,kkk) + scatt_matrix_tmp1_24/azimuth0_num
                !
                !                        scatter_matrix(3,ll,1,jj,kkk) = scatter_matrix(3,ll,1,jj,kkk) + scatt_matrix_tmp1_31/azimuth0_num
                !                        scatter_matrix(3,ll,2,jj,kkk) = scatter_matrix(3,ll,2,jj,kkk) + scatt_matrix_tmp1_32/azimuth0_num
                !                        scatter_matrix(3,ll,3,jj,kkk) = scatter_matrix(3,ll,3,jj,kkk) + scatt_matrix_tmp1_33/azimuth0_num
                !                        scatter_matrix(3,ll,4,jj,kkk) = scatter_matrix(3,ll,4,jj,kkk) + scatt_matrix_tmp1_34/azimuth0_num
                !
                !                        scatter_matrix(4,ll,1,jj,kkk) = scatter_matrix(4,ll,1,jj,kkk) + scatt_matrix_tmp1_41/azimuth0_num
                !                        scatter_matrix(4,ll,2,jj,kkk) = scatter_matrix(4,ll,2,jj,kkk) + scatt_matrix_tmp1_42/azimuth0_num
                !                        scatter_matrix(4,ll,3,jj,kkk) = scatter_matrix(4,ll,3,jj,kkk) + scatt_matrix_tmp1_43/azimuth0_num
                !                        scatter_matrix(4,ll,4,jj,kkk) = scatter_matrix(4,ll,4,jj,kkk) + scatt_matrix_tmp1_44/azimuth0_num

1244            continue  ! phi0
                ! calculate the summation of the scattering matrix in the whole sphere
                emis_vector_tmp1_11(ll+(kk-1)*qua_num) = scatter_matrix(1,ll,1,jj,kkk1)*thet_weights*2.*pi
                emis_vector_tmp1_12(ll+(kk-1)*qua_num) = scatter_matrix(1,ll,2,jj,kkk1)*thet_weights*2.*pi
            !	 emis_vector_tmp1_13(ll+(kk-1)*qua_num) =scatter_matrix(1,ll,3,jj,kkk)*thet_weights*2.*pi
            ! 	 emis_vector_tmp1_14(ll+(kk-1)*qua_num) = scatter_matrix(1,ll,4,jj,kkk)*thet_weights*2.*pi
            !c		write(*,*)'tmp1:',emis_vector_tmp1_11(ll+(kk-1)*qua_num)

1243        continue ! thet
1242    continue

        emis_vector(1,jj,:) = extinct_matrix(1,1,jj,:) - sum(emis_vector_tmp1_11)
        emis_vector(2,jj,:) = extinct_matrix(1,2,jj,:) - sum(emis_vector_tmp1_12)
    !            emis_vector(3,jj,ii) = extinct_matrix(1,3,jj,ii) - sum(emis_vector_tmp1_13)
    !            emis_vector(4,jj,ii) = extinct_matrix(1,4,jj,ii) - sum(emis_vector_tmp1_14)

1241 continue ! thet0



     return
 end subroutine matrix_cal

