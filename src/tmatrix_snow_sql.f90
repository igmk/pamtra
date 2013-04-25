subroutine tmatrix_snow_sql(f, wc, t, nc, &
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
  use settings, only: use_snow_db, as_ratio, softsphere_adjust, nummu, nstokes
!   use rt_utilities, only: lobatto_quadrature
  use report_module

  use tmat_snow_db

  implicit none

  integer, intent(in) :: nbins

  real(kind=dbl), intent(in) :: f, &  ! frequency [ghz]
                                t, dia1, dia2
  real(kind=dbl) :: wavelength, wave_num, freq
  real(kind=dbl) :: ad, bd, alpha, gamma, a_m,b

  integer :: ir
  real(kind=dbl) :: del_d, tot_mass, wc, bin_wgt

  real(kind=dbl) :: distribution
  real(kind=dbl), dimension(nbins) :: particle_mass, ndens
  character :: aerodist*1 

  real(kind=dbl), dimension(nstokes,nummu,nstokes,nummu,4), intent(out) :: scatter_matrix
  real(kind=dbl), dimension(nstokes,nstokes,nummu,2), intent(out) :: extinct_matrix
  real(kind=dbl), dimension(nstokes,nummu,2), intent(out) :: emis_vector
  real(kind=dbl), intent(out) :: diameter(nbins)
  real(kind=dbl), intent(out) :: back_spec(nbins)
  real(kind=dbl), dimension(nstokes,nummu,nstokes,nummu,2) :: scat_mat_sgl
  real(kind=dbl), dimension(nstokes,nstokes,nummu) :: ext_mat_sgl
  real(kind=dbl), dimension(nstokes,nummu) :: emis_vec_sgl

!   integer :: j1, j2, l, j, i1, i2, l1, l2, i


  real(kind=dbl) :: ntot,nc

    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'tmatrix_ice_sql'

interface
subroutine tmatrix_sql(freq,t,as_ratio,diameter,particle_mass,scatter_matrix,extinct_matrix,emis_vector)

  use kinds
  use settings, only : nummu, nstokes
  implicit none

  real(kind=dbl), intent(in) :: freq, t
  real(kind=dbl), intent(in) :: as_ratio
  real(kind=dbl), intent(in) :: particle_mass
  real(kind=dbl), intent(in) :: diameter

  real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector
end subroutine tmatrix_sql
end interface

  freq = f*1.d9
  wavelength = c/freq !
  wave_num = 2.0_dbl*pi/wavelength

  scatter_matrix = 0.0_dbl
  extinct_matrix = 0.0_dbl
  emis_vector = 0.0_dbl
  !   integration loop over diameter of spheres
  tot_mass = 0.0_dbl
!   del_d = 2.d-4

  if (nbins .gt. 0) then
    del_d = (dia2 - dia1) / nbins
  else
    del_d = 1.d0
  end if

  do ir = 1, nbins
!       diameter(ir) = ir * del_d
     diameter(ir) = dia1 + (ir - 1) * del_d
     ndens(ir) = distribution(ad, bd, alpha, gamma, diameter(ir), aerodist)  ! number density [1/m^4]
     particle_mass(ir) = a_m*diameter(ir)**b ! particle mass [kg]
     tot_mass = tot_mass + ndens(ir)*del_d*particle_mass(ir)   ! [kg/m^3]
     if ((ir .eq. nbins) .and. (tot_mass/wc*100.0_dbl .lt. 99.9_dbl)) then
        ndens(ir) = ndens(ir) + (wc-tot_mass)/(particle_mass(ir)*del_d)
     end if


end do


  do ir = 1, nbins


    call tmatrix_sql(freq,t,as_ratio,diameter(ir),particle_mass(ir),&
      scat_mat_sgl,ext_mat_sgl,emis_vec_sgl)

     bin_wgt = ndens(ir)*del_d

     scatter_matrix(:,:,:,:,1:2) = scatter_matrix(:,:,:,:,1:2) + scat_mat_sgl*bin_wgt
     extinct_matrix(:,:,:,1) = extinct_matrix(:,:,:,1) + ext_mat_sgl*bin_wgt
     emis_vector(:,:,1) = emis_vector(:,:,1) + emis_vec_sgl*bin_wgt

    back_spec(ir) = 4*pi*ndens(ir)*scat_mat_sgl(1,16,1,16,2) !scatter_matrix(A,B;C;D;E) backscattering is M11 of Mueller or Scattering Matrix (A;C=1), in quadrature 2 (E) first 16 (B) is 180deg (upwelling), 2nd 16 (D) 0deg (downwelling). this definition is lokkiing from BELOW, scatter_matrix(1,16,1,16,3) would be from above!


  end do



!fill up the matrices
scatter_matrix(:,:,:,:,4) = scatter_matrix(:,:,:,:,1) 
scatter_matrix(:,:,:,:,3) = scatter_matrix(:,:,:,:,2)
extinct_matrix(:,:,:,2) = extinct_matrix(:,:,:,1)
emis_vector(:,:,2) = emis_vector(:,:,1)

    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)


  return 
end subroutine tmatrix_snow_sql	
