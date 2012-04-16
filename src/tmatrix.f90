subroutine tmatrix(f, mindex, dia1, dia2, nbins, &
     ad, bd, alpha, gamma, aerodist,density,wc,scatter_matrix,extinct_matrix, emis_vector)

  ! note that mindex has the convention with negative imaginary part
  !     
  ! computes the mie scattering properties for a gamma or lognormal 
  ! distribution of spheres.                                          

  use kinds
  use constants, only: pi, c

  implicit none

  integer :: nbins

  real(kind=dbl), intent(in) :: f  ! frequency [ghz]

  real(kind=dbl) :: wavelength, dia1, dia2, wave_num, freq
  real(kind=dbl) :: ad, bd, alpha, gamma
  complex(kind=dbl) :: mindex 
  real(kind=dbl) :: extinction, albedo, back_scatt
  integer, parameter :: nquad = 16
  integer, parameter :: nstokes = 2
  integer :: i, l, m, ir, azimuth_num, azimuth0_num
  real(kind=dbl) :: x, del_d, diameter, ndens, tmp, tot_mass, wc, density
  real(kind=dbl) :: qext, qscat, qback, scatter
  real(kind=dbl) :: distribution
  real(kind=dbl) :: eu_alpha, eu_beta, bin_wgt, as_ratio, equiv_radius

  character :: aerodist*1 

  real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
  real(kind=dbl), dimension(4,4,nquad,2), intent(out) :: extinct_matrix
  real(kind=dbl), dimension(4,nquad,2), intent(out) :: emis_vector
  real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4) :: scat_mat_sgl
  real(kind=dbl), dimension(4,4,nquad,2) :: ext_mat_sgl
  real(kind=dbl), dimension(4,nquad,2) :: emis_vec_sgl
  real(kind=dbl), dimension(nquad) :: qua_angle, qua_weights
  integer :: l1, j1, l2, j2, j, i1, i2

  freq = f*1.e9_dbl
  wavelength = c/freq !
  wave_num = 2.0_dbl*pi/wavelength

  eu_alpha = 0.0_dbl
  eu_beta = 0.0_dbl
  azimuth_num = 30
  azimuth0_num = 1
  as_ratio = 0.4_dbl

  scatter_matrix = 0.0_dbl
  extinct_matrix = 0.0_dbl
  emis_vector = 0.0_dbl
  !   integration loop over diameter of spheres
!  if (nbins .gt. 0) del_d = (dia2 - dia1) / nbins
!  tot_mass = 0.0_dbl
!  do ir = 1, nbins+1
!     diameter = dia1 + (ir - 1) * del_d
!     ndens = distribution(ad, bd, alpha, gamma, diameter, aerodist)  ! number density
!
!     if ((ir .eq. 1 .or. ir .eq. nbins+1) .and. nbins .gt. 0) then
!        ndens = 0.5_dbl * ndens
!     end if
!     tot_mass = tot_mass + ndens*del_d*pi/6.0_dbl*density*diameter**3.0_dbl
!     if ((ir .eq. nbins+1) .and. (tot_mass/wc*100.0_dbl .lt. 99.9_dbl)) then
!      ndens = ndens + (wc-tot_mass)/(del_d*pi/6.0_dbl*density*diameter**3.0_dbl)
!      tot_mass = wc
!     end if
     del_d = 250.0d-6
     diameter = 150.0d-5
     ndens = 4776.1196698715538_dbl
     ndens = 14814.814557263885_dbl
     bin_wgt = ndens!*del_d
     equiv_radius = 0.5_dbl*diameter*as_ratio**(1.0_dbl/3.0_dbl)
     print*, 'L',nquad,freq,wave_num,mindex,equiv_radius,ndens
     call matrix_cal('L',nquad,freq,wave_num,mindex,equiv_radius,nstokes,&
             as_ratio, eu_alpha, eu_beta, azimuth_num, azimuth0_num, &
           scat_mat_sgl,ext_mat_sgl,emis_vec_sgl)
     scatter_matrix = scatter_matrix + scat_mat_sgl*bin_wgt
     extinct_matrix = extinct_matrix + ext_mat_sgl*bin_wgt
     emis_vector = emis_vector + emis_vec_sgl*bin_wgt

!     ! WRITE OUT THE MATRIX FILE IN THE FILE OPENED AT THE BEGINNING OF THIS SUBROUTINE
!             call lobatto_quadrature(nquad,qua_angle(1:nquad),qua_weights(1:nquad))
!
!      WRITE(1232,*) 'scatter'
!      DO L1 = 1, 2
!        DO J1 = 1, nquad
!          DO L2 = 1, 2
!            L = 2*(L2-1)+L1
!            DO J2 = 1, nquad
!              WRITE(1232,*) (-1.0)**(L1+1.0)*QUA_ANGLE(J1),&
!                (-1.0)**(L2+1.0)*QUA_ANGLE(J2), 0E0
!              DO I2 = 1, NSTOKES
!                WRITE(1232,*)&
!            (SCAtter_MATRIX(I2,J2,I1,J1,L), I1=1,NSTOKES), 0E0, 0
!              ENDDO
!              DO I2 = NSTOKES+1,4
!                WRITE(1232,*) 0E0, 0E0, 0E0, 0E0
!              ENDDO
!            ENDDO
!          ENDDO
!        ENDDO
!      ENDDO

!      WRITE(1232,*) 'extinct'
!      DO L = 1, 2
!        DO J = 1, nquad
!          WRITE(1232,*) (-1.0)**(L+1.0)*QUA_ANGLE(J)
!          DO I2 = 1, NSTOKES
!           WRITE(1232,*)(EXTinct_MATRIX(I2,I1,J,L), I1=1,NSTOKES), 0E0, 0E0
!          ENDDO
!          DO I2 = NSTOKES+1,4
!             WRITE(1232,*) 0E0, 0E0, 0E0, 0E0
!          ENDDO
!        ENDDO
!      ENDDO

!      WRITE(1232,*) 'emis'
!      DO L = 1, 2
!        DO J = 1, nquad
!          WRITE(1232,*) (-1.0)**(L+1.0)*QUA_ANGLE(J),&
!            (EMIs_VECTOR(I,J,L), I=1,NSTOKES), 0E0, 0E0
!        ENDDO
!      ENDDO


!  end do

  return 
end subroutine tmatrix
