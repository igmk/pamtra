subroutine tmatrix_ice(f, wc, t, nc, &
     ad, bd, alpha, gamma, a_m,b,aerodist,dia1,dia2,nbins,scatter_matrix,extinct_matrix, emis_vector,&
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
  use constants, only: pi, c, im
  use report_module
  use settings, only:as_ratio_ice!,  use_ice_db, 

!   use tmat_ice_db

  implicit none

  integer, intent(in) :: nbins

  real(kind=dbl), intent(in) :: f, &  ! frequency [ghz]
                                t, dia1, dia2
  real(kind=dbl) :: wavelength, wave_num, freq
  real(kind=dbl) :: ad, bd, alpha, gamma, a_m,b
  complex(kind=ext) :: mindex
  complex(kind=dbl) :: msphere, eps_mix,m_ice

  integer, parameter :: nquad = 16
  integer, parameter :: nstokes = 2
  integer :: ir, azimuth_num, azimuth0_num
  real(kind=dbl) :: del_d, ndens, tot_mass, wc

  real(kind=dbl) :: distribution
  real(kind=dbl) :: eu_alpha, eu_beta, bin_wgt, equiv_radius, particle_mass

  character :: aerodist*1 

  real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
  real(kind=dbl), dimension(nstokes,nstokes,nquad,2), intent(out) :: extinct_matrix
  real(kind=dbl), dimension(nstokes,nquad,2), intent(out) :: emis_vector
  real(kind=dbl), intent(out) :: diameter(nbins)
  real(kind=dbl), intent(out) :: back_spec(nbins)
  real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,2) :: scat_mat_sgl
  real(kind=dbl), dimension(nstokes,nstokes,nquad) :: ext_mat_sgl
  real(kind=dbl), dimension(nstokes,nquad) :: emis_vec_sgl


  real(kind=dbl) :: ntot,nc, density_eff, refre, refim

    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'tmatrix_ice'


  


  interface

    subroutine tmatrix_calc(quad,qua_num,frequency,wave_num,ref_index,axi, nstokes,&
    as_ratio_ice, alpha, beta, azimuth_num, azimuth0_num,&
    scatter_matrix,extinct_matrix,emis_vector)
      use kinds
      implicit none
      character(1), intent(in) :: quad
      integer, intent(in) :: qua_num
      real(kind=dbl), intent(in) :: frequency
      real(kind=dbl), intent(in) :: wave_num
      complex(kind=ext) :: ref_index
      real(kind=dbl), intent(in) :: axi
      integer, intent(in) :: nstokes
      real(kind=dbl), intent(in) :: as_ratio_ice
      real(kind=dbl), intent(in) :: alpha
      real(kind=dbl), intent(in) :: beta
      integer, intent(in) :: azimuth_num
      integer, intent(in) :: azimuth0_num
      real(kind=dbl), intent(out), dimension(nstokes,qua_num,nstokes,qua_num,2) :: scatter_matrix
      real(kind=dbl), intent(out), dimension(nstokes,nstokes,qua_num) :: extinct_matrix
      real(kind=dbl), intent(out), dimension(nstokes,qua_num) :: emis_vector
    end subroutine tmatrix_calc
  end interface


    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)


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

  if (nbins .gt. 0) then
    del_d = (dia2 - dia1) / nbins
  else
    del_d = 1.d0
  end if


  do ir = 1, nbins
     diameter(ir) = dia1 + (ir - 1) * del_d
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
!     density = 6.d0*particle_mass/(diameter(ir)**3.d0*pi*as_ratio_ice)
     if (ir .eq. 1 .or. ir .eq. nbins) then
        ndens = 0.5_dbl * ndens
     end if
     tot_mass = tot_mass + ndens*del_d*particle_mass   ! [kg/m^3]
!      if ((ir .eq. nbins) .and. (tot_mass/wc*100.0_dbl .lt. 99.9_dbl)) then
!         ndens = ndens + (wc-tot_mass)/(particle_mass*del_d)
!      end if
     bin_wgt = ndens*del_d
     ntot = ntot + ndens*del_d
!      if (use_ice_db) then
!         call get_ice_data(real(f),real(t),ir,scat_mat_sgl,ext_mat_sgl,emis_vec_sgl)
! !        print*, ir,diameter(ir),ndens*del_d, tot_mass/wc*100.,ntot/nc*100.
!      else

	  !XINXINs code
!         CALL CAL_REFRACTIVE_INDEX('S',t,freq, diameter(ir), as_ratio_ice, particle_mass*1.d3, equiv_radius, mindex)
    !     write(25,*) ir,equiv_radius,diameter(ir),density,ntot,ntot/nc*100.0_dbl,tot_mass,tot_mass/wc*100.0_dbl

	  !diameter of sphere with same mass
	  equiv_radius = 0.5_dbl*diameter(ir)*as_ratio_ice**(1.0_dbl/3.0_dbl)
	  density_eff = (3.d0/4.d0 * particle_mass) / (pi * equiv_radius**3)
	  call ref_ice(t, f, refre, refim)
	  m_ice = refre-Im*refim  ! mimicking a

	if (density_eff > 917.d0) then
	  if (verbose >= 0) print*, "WARNING changed density from ", density_eff, "kg/m3 to 917 kg/m3 for d=", diameter(ir)
	  density_eff = 917.d0
	end if

	  msphere = eps_mix((1.d0,0.d0),m_ice,density_eff)
	  mindex =conjg(msphere) !different convention

    if (verbose >= 4) print*, "density_eff, equiv_radius, diameter(ir),mindex"
    if (verbose >= 4) print*, density_eff, equiv_radius, diameter(ir),mindex

call tmatrix_calc('l',nquad,freq,wave_num,mindex,equiv_radius,nstokes,&
            as_ratio_ice, eu_alpha, eu_beta, azimuth_num, azimuth0_num, &
            scat_mat_sgl,ext_mat_sgl,emis_vec_sgl)
!      end if
     scatter_matrix(:,:,:,:,1:2) = scatter_matrix(:,:,:,:,1:2) + scat_mat_sgl*bin_wgt
     extinct_matrix(:,:,:,1) = extinct_matrix(:,:,:,1) + ext_mat_sgl*bin_wgt
     emis_vector(:,:,1) = emis_vector(:,:,1) + emis_vec_sgl*bin_wgt

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

!fill up the matrices
scatter_matrix(:,:,:,:,4) = scatter_matrix(:,:,:,:,1) 
scatter_matrix(:,:,:,:,3) = scatter_matrix(:,:,:,:,2)
extinct_matrix(:,:,:,2) = extinct_matrix(:,:,:,1)
emis_vector(:,:,2) = emis_vector(:,:,1)

    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)


  return 
end subroutine tmatrix_ice
