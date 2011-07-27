subroutine grau_ssp(f,qg,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4, ng)
  
  use kinds
  use nml_params, only: verbose, lphase_flag, n_0grauDgrau, EM_grau, n_moments
  use constants, only: pi, im
  use double_moments_module

  implicit none

  integer :: nbins, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qg,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl), optional, intent(in) :: ng

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, gwc, ad, bd, alpha, gamma, b_grau, a_mgrau, ng_abs

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex, m_air

  real(kind=dbl) :: spec2abs, gammln

  character(1) :: dist_name

  if (verbose .gt. 1) print*, 'Entering grau_ssp'

	call ref_ice(t, f, refre, refim)
	mindex = refre-Im*refim
	m_air = 1.0d0 - 0.0d0 * Im

	gwc =  spec2abs(qg,t,p,q) ! [kg/m^3]

  if (n_moments .eq. 1) then
	dia1 = 1.d-5 	! minimum diameter [m]
	dia2 = 1.d-2	! minimum diameter [m]

	b_grau = 3.1d0
	a_mgrau = 169.6d0
	ad = n_0grauDgrau*1.d6
	bd = (exp(gammln(b_grau + 1)) * a_mgrau * ad/gwc)**(1.0d0 /(1.0d0 + b_grau)) !  [m**-1]
	nbins = 100
	alpha = 0.d0
	gamma = 1.d0
	dist_name='C'
  else if (n_moments .eq. 2) then
    if (.not. present(ng)) stop 'STOP in routine grau_ssp'
    ng_abs = spec2abs(ng,t,p,q) 							! [#/m^3]
    call double_moments(gwc,ng_abs,gamma_graupel(1),gamma_graupel(2),gamma_graupel(3),gamma_graupel(4), &
    	ad,bd,alpha,gamma,a_mgrau, b_grau)
    nbins = 100
    dia1 = 1.d-5	! minimum diameter [m]
    dia2 = 1.d-2	! maximum diameter [m]
    dist_name='G'
  else
    stop'Number of moments is not specified'
  end if

	if (EM_grau .eq. 'icesf') then
	  call mie_densitysizedep_spheremasseq(f, mindex,      &
		a_mgrau, b_grau, dia1, dia2, nbins, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, nlegen, legen, legen2, legen3,        &
		legen4, dist_name)
	elseif (EM_grau .eq. 'surus') then 
	  call mie_icefactor(f, t,mindex,      &
		a_mgrau, b_grau, dia1, dia2, nbins, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
		LEGEN4, dist_name,0.815*1.e-3*f+0.0112,44)
	else 
	    write (*, *) 'no em mod for grau'
	    stop
	end if
  if (verbose .gt. 1) print*, 'Exiting grau_ssp'

  return

end subroutine grau_ssp
