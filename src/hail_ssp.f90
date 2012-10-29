subroutine hail_ssp(f,hwc,t,maxleg,nc,kext, salb, back,  &
     nlegen, legen, legen2, legen3, legen4,&
     scatter_matrix,extinct_matrix, emis_vector)

  use kinds
  use nml_params, only: verbose, lphase_flag, EM_hail, n_moments, SD_hail, hail_density,&
	nstokes
  use constants, only: pi, im
  use double_moments_module
  use conversions

  implicit none

  integer :: nbins, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
       hwc,&
       t,&
       f

  real(kind=dbl), intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, ad, bd, alpha, gamma, b_hail, a_mhail

  real(kind=dbl), intent(out) :: &
       kext,&
       salb,&
       back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4
    integer, parameter ::  nquad = 16
    real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: scatter_matrix
    real(kind=dbl), dimension(nstokes,nstokes,nquad,2), intent(out) :: extinct_matrix
    real(kind=dbl), dimension(nstokes,nquad,2), intent(out) :: emis_vector
  complex(kind=dbl) :: mindex, m_air

  if (verbose .gt. 1) print*, 'Entering hail_ssp'

  call ref_ice(t, f, refre, refim)
  mindex = refre-Im*refim
  m_air = 1.0d0 - 0.0d0 * Im

  if (nc .eq. 0) stop 'STOP in routine hail_ssp'
  call double_moments(hwc,nc,gamma_hail(1),gamma_hail(2),gamma_hail(3),gamma_hail(4), &
       ad,bd,alpha,gamma,a_mhail, b_hail)
  nbins = 100
  dia1 = 1.d-8	! minimum diameter [m]
  dia2 = 2.d-2	! maximum diameter [m]

  if (EM_hail .eq. 'densi' .or. EM_hail .eq. 'surus') then
  	if (EM_hail .eq. 'surus') hail_density = 0.815*f+11.2d0
     call mie_densitydep_spheremasseq(f, t, mindex,      &
          a_mhail, b_hail, dia1, dia2, nbins, maxleg,   &
          ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
          back, nlegen, legen, legen2, legen3,        &
          legen4, SD_hail,hail_density,hwc)
      scatter_matrix= 0.d0
      extinct_matrix= 0.d0
      emis_vector= 0.d0
  else
     write (*, *) 'no em mod for hail'
     stop
  end if
  if (verbose .gt. 1) print*, 'Exiting hail_ssp'

  return

end subroutine hail_ssp
