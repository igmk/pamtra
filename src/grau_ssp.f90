subroutine grau_ssp(f,qg,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4)
  
  use kinds
  use nml_params, only: verbose, lphase_flag, n_0grauDgrau, EM_grau
  use constants, only: pi, im

  implicit none

  integer :: numrad, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qg,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: rad1, rad2, gwc, ad, bd, alpha, gamma, b_grau, a_mgrau

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex, m_air

  real(kind=dbl) :: spec2abs, gammln

  if (verbose .gt. 1) print*, 'Entering grau_ssp'

	b_grau = 3.1d0
	a_mgrau = 169.6d0 

	call ref_ice(t, f, refre, refim)
	mindex = refre-Im*refim
	m_air = 1.0d0 - 0.0d0 * Im

	rad1 = 1.d-5 
	rad2 = 1.d-2
	numrad = 100
	alpha = 0.
	gamma = 1.

	gwc =  spec2abs(qg,t,p,q) ! [kg/m^3]

	ad = n_0grauDgrau*1.d6
	bd = (exp(gammln(b_grau + 1)) * a_mgrau * ad/gwc)**(1.0d0 /(1.0d0 + b_grau)) !  [m**-1]

	if (EM_grau .eq. 'icesf') then
	  call mie_densitysizedep_spheremasseq(f, mindex,      &
		a_mgrau, b_grau, rad1/2., rad2/2., numrad, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, nlegen, legen, legen2, legen3,        &
		legen4, 'C')
	elseif (EM_grau .eq. 'surus') then 
	  call mie_icefactor(f, t,mindex,      &
		a_mgrau, b_grau, rad1/2., rad2/2., numrad, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
		LEGEN4, 'C',0.815*1.e-3*f+0.0112,44)
	else 
	    write (*, *) 'no em mod for grau'
	    stop
	end if
  if (verbose .gt. 1) print*, 'Exiting grau_ssp'

  return

end subroutine grau_ssp