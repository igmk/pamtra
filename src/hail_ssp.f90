subroutine hail_ssp(f,qh,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4, nh)

  use kinds
  use nml_params, only: verbose, lphase_flag, EM_grau, n_moments
  use constants, only: pi, im
  use double_moments_module

  implicit none

  integer :: nbins, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qh,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl), optional, intent(in) :: nh

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, hwc, ad, bd, alpha, gamma, b_hail, a_mhail, nh_abs

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex, m_air

  real(kind=dbl) :: spec2abs

  character(1) :: dist_name

  if (verbose .gt. 1) print*, 'Entering hail_ssp'

	call ref_ice(t, f, refre, refim)
	mindex = refre-Im*refim
	m_air = 1.0d0 - 0.0d0 * Im
	hwc =  spec2abs(qh,t,p,q) ! [kg/m^3]


    if (.not. present(nh)) stop 'STOP in routine hail_ssp'
    nh_abs = spec2abs(nh,t,p,q) 							! [#/m^3]
    call double_moments(hwc,nh_abs,gamma_hail(1),gamma_hail(2),gamma_hail(3),gamma_hail(4), &
    	ad,bd,alpha,gamma,a_mhail, b_hail)
    nbins = 100
    dia1 = 1.d-5	! minimum diameter [m]
    dia2 = 2.d-2	! maximum diameter [m]
    dist_name='G'

	if (EM_grau .eq. 'icesf') then
	  call mie_densitysizedep_spheremasseq(f, mindex,      &
		a_mhail, b_hail, dia1, dia2, nbins, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, nlegen, legen, legen2, legen3,        &
		legen4, dist_name)
	elseif (EM_grau .eq. 'surus') then
	  call mie_icefactor(f, t,mindex,      &
		a_mhail, b_hail, dia1, dia2, nbins, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
		LEGEN4, dist_name,0.815*1.e-3*f+0.0112,44)
	else
	    write (*, *) 'no em mod for hail'
	    stop
	end if
  if (verbose .gt. 1) print*, 'Exiting hail_ssp'

  return

end subroutine hail_ssp
