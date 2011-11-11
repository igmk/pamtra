subroutine grau_ssp(f,gwc,t,maxleg,kext, salb, back,  &
     nlegen, legen, legen2, legen3, legen4, nc)

  use kinds
  use nml_params, only: verbose, lphase_flag, n_0grauDgrau, EM_grau, n_moments, SD_grau, graupel_density
  use constants, only: pi, im
  use double_moments_module
  use conversions

  implicit none

  integer :: nbins, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
       gwc,&
       t,&
       f

  real(kind=dbl), optional, intent(in) :: nc

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, ad, bd, alpha, gamma, b_grau, a_mgrau

  real(kind=dbl), intent(out) :: &
       kext,&
       salb,&
       back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex, m_air

  real(kind=dbl) :: gammln

  if (verbose .gt. 1) print*, 'Entering grau_ssp'

  call ref_ice(t, f, refre, refim)
  mindex = refre-Im*refim
  m_air = 1.0d0 - 0.0d0 * Im

  if (n_moments .eq. 1) then
	if (SD_grau .eq. 'C') then
     dia1 = 1.d-5 	! minimum diameter [m]
     dia2 = 1.d-2	! minimum diameter [m]

     b_grau = 3.1d0
     a_mgrau = 169.6d0
     ad = n_0grauDgrau*1.d6
     bd = (exp(gammln(b_grau + 1)) * a_mgrau * ad/gwc)**(1.0d0 /(1.0d0 + b_grau)) !  [m**-1]
     nbins = 100
     alpha = 0.d0
     gamma = 1.d0
	else if (SD_grau .eq. 'M') then
     dia1 = 1.d-10 	! minimum diameter [m]
     dia2 = 1.d-2	! minimum diameter [m]
     b_grau = 2.8d0
     a_mgrau = 19.6d0
     bd = (gwc/(a_mgrau*5.d5*exp(gammln(1.d0+b_grau))))**(1.d0/(-0.5d0-b_grau)) !  [m**-1]
     ad = 5.d5*bd**0.5
	 dia2 = log(ad)/bd
     if (dia2 .gt. 1.d-2) dia2 = 1.d-2
	 nbins = 100
     alpha = 0.d0
     gamma = 1.d0
	end if
  else if (n_moments .eq. 2) then
     if (.not. present(nc)) stop 'STOP in routine grau_ssp'
     call double_moments(gwc,nc,gamma_graupel(1),gamma_graupel(2),gamma_graupel(3),gamma_graupel(4), &
          ad,bd,alpha,gamma,a_mgrau, b_grau)
     nbins = 100
     dia1 = 1.d-5	! minimum diameter [m]
     dia2 = 1.d-2	! maximum diameter [m]
  else
     stop 'Number of moments is not specified'
  end if

  if (EM_grau .eq. 'densi' .or. EM_grau .eq. 'surus') then
     if (EM_grau .eq. 'surus') graupel_density =  0.815*f+11.2d0
     call mie_densitydep_spheremasseq(f, t,mindex,      &
          a_mgrau, b_grau, dia1, dia2, nbins, maxleg,   &
          ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
          back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
          LEGEN4, SD_grau,graupel_density,gwc)
  else 
     write (*, *) 'no em mod for grau'
     stop
  end if
  if (verbose .gt. 1) print*, 'Exiting grau_ssp'

  return

end subroutine grau_ssp
