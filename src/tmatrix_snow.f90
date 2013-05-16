subroutine tmatrix_snow(f, wc, t, nc, &
     ad, bd, alpha, gamma, a_m,b,aerodist,nbins,scatter_matrix,extinct_matrix, emis_vector,&
     diameter, back_spec)

  ! note that mindex has the convention with negative imaginary part
  !     
  ! computes the mie scattering properties for a gamma or lognormal 
  ! distribution of spheres.
  !
  ! input:
  !     f       frequency [GHz]
  !     t       temperature [K]

  !out
  !...
  !diameter: diameter spectrum (maximum extent) [m]
  !back_spec: backscattering cross section per volume per del_d [m²/m⁴]



  use kinds
  use constants, only: pi, c
  use nml_params, only: use_snow_db, as_ratio

  use tmat_snow_db

  implicit none

  integer, intent(in) :: nbins

  real(kind=dbl), intent(in) :: f, &  ! frequency [ghz]
                                t
  real(kind=dbl) :: wavelength, wave_num, freq
  real(kind=dbl) :: ad, bd, alpha, gamma, a_m,b
  complex(kind=dbl) :: mindex
  real(kind=dbl) :: extinction, albedo, back_scatt
  integer, parameter :: nquad = 16
  integer, parameter :: nstokes = 2
  integer :: i, l, m, ir, azimuth_num, azimuth0_num
  real(kind=dbl) :: del_d, ndens, tmp, tot_mass, wc, density
  real(kind=dbl) :: qext, qscat, qback, scatter
  real(kind=dbl) :: distribution
  real(kind=dbl) :: eu_alpha, eu_beta, bin_wgt, equiv_radius, particle_mass

  character :: aerodist*1 

  real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
  real(kind=dbl), dimension(nstokes,nstokes,nquad,2), intent(out) :: extinct_matrix
  real(kind=dbl), dimension(nstokes,nquad,2), intent(out) :: emis_vector
  real(kind=dbl), intent(out) :: diameter(nbins)
  real(kind=dbl), intent(out) :: back_spec(nbins)
  real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4) :: scat_mat_sgl
  real(kind=dbl), dimension(nstokes,nstokes,nquad,2) :: ext_mat_sgl
  real(kind=dbl), dimension(nstokes,nquad,2) :: emis_vec_sgl
  real(kind=dbl), dimension(nquad) :: qua_angle, qua_weights
  integer :: l1, j1, l2, j2, j, i1, i2
  real(kind=dbl) :: gammln, ntot,nc

  freq = f*1.d9
  wavelength = c/freq !
  wave_num = 2.0_dbl*pi/wavelength

  eu_alpha = 0.0_dbl
  eu_beta = 0.0_dbl
  azimuth_num = 30
  azimuth0_num = 1

  scatter_matrix = 0.0_dbl
  extinct_matrix = 0.0_dbl
  emis_vector = 0.0_dbl
  !   integration loop over diameter of spheres
  tot_mass = 0.0_dbl
  ntot = 0.d0
  del_d = 2.d-4
  do ir = 1, nbins
     diameter(ir) = ir * del_d
     ndens = distribution(ad, bd, alpha, gamma, diameter(ir), aerodist)  ! number density [1/m^4]
     particle_mass = a_m*diameter(ir)**b ! particle mass [kg]
!         IF (diameter(ir)*1.e6.GE.0.005E4.AND.diameter(ir)*1.e6.LE.0.2E4)THEN ! 50um~200um
!            PARTICLE_MASS = 0.003D0*(diameter(ir)*1.e6/1.0E4)**2.0D0*1.e-3
!         ELSEIF (diameter(ir)*1.e6.GT.0.2E4.AND.diameter(ir)*1.e6.LE.2.0E4)THEN ! 200um~2cm
!            PARTICLE_MASS = 0.0067D0*(diameter(ir)*1.e6/1.0E4)**2.5D0*1.e-3
!         ELSEIF (diameter(ir)*1.e6.GT.2.0E4)THEN ! >2cm
!            PARTICLE_MASS = 0.0047D0*(diameter(ir)*1.e6/1.0E4)**3.0D0*1.e-3
!         ENDIF
!
!     density = 6.d0*particle_mass/(diameter(ir)**3.d0*pi*as_ratio)
     if (ir .eq. 1 .or. ir .eq. nbins) then
        ndens = 0.5_dbl * ndens
     end if
     tot_mass = tot_mass + ndens*del_d*particle_mass   ! [kg/m^3]
     if ((ir .eq. nbins) .and. (tot_mass/wc*100.0_dbl .lt. 99.9_dbl)) then
        ndens = ndens + (wc-tot_mass)/(particle_mass*del_d)
     end if
     bin_wgt = ndens*del_d
     ntot = ntot + ndens*del_d
     if (use_snow_db) then
        call get_snow_data(real(f),real(t),ir,scat_mat_sgl,ext_mat_sgl,emis_vec_sgl)
!        print*, ir,diameter(ir),ndens*del_d, tot_mass/wc*100.,ntot/nc*100.
     else
        equiv_radius = 0.5_dbl*diameter(ir)*as_ratio**(1.0_dbl/3.0_dbl)
        CALL CAL_REFRACTIVE_INDEX('S',t,freq, diameter(ir), as_ratio, particle_mass*1.d3, equiv_radius, mindex)
        print*, ir,diameter(ir),ndens*del_d, tot_mass/wc*100.,ntot/nc*100., mindex
    !     write(25,*) ir,equiv_radius,diameter(ir),density,ntot,ntot/nc*100.0_dbl,tot_mass,tot_mass/wc*100.0_dbl
        call matrix_cal('S',nquad,freq,wave_num,mindex,equiv_radius,nstokes,&
            as_ratio, eu_alpha, eu_beta, azimuth_num, azimuth0_num, &
            scat_mat_sgl,ext_mat_sgl,emis_vec_sgl)
     end if
     scatter_matrix = scatter_matrix + scat_mat_sgl*bin_wgt
     extinct_matrix = extinct_matrix + ext_mat_sgl*bin_wgt
     emis_vector = emis_vector + emis_vec_sgl*bin_wgt
     back_spec(ir) = 4*pi*ndens*scat_mat_sgl(1,16,1,16,2) !scatter_matrix(A,B;C;D;E) backscattering is M11 of Mueller or Scattering Matrix (A;C=1), in quadrature 2 (E) first 16 (B) is 180deg (upwelling), 2nd 16 (D) 0deg (downwelling). this definition is lokkiing from BELOW, scatter_matrix(1,16,1,16,3) would be from above!





!     ! WRITE OUT THE MATRIX FILE IN THE FILE OPENED AT THE BEGINNING OF THIS SUBROUTINE
!             call lobatto_quadrature(nquad,qua_angle(1:nquad),qua_weights(1:nquad))
!
!      WRITE(1232,*) ir, ' scatter'
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
!
!      WRITE(1232,*) ir, ' extinct'
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
!
!      WRITE(1232,*) ir, ' emis'
!      DO L = 1, 2
!        DO J = 1, nquad
!          WRITE(1232,*) (-1.0)**(L+1.0)*QUA_ANGLE(J),&
!            (EMIs_VECTOR(I,J,L), I=1,NSTOKES), 0E0, 0E0
!        ENDDO
!      ENDDO
end do

  return 
end subroutine tmatrix_snow
