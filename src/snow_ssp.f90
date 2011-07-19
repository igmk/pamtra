! Subroutine for the setup of the parameters of the snow particle size distribution.
!
!
subroutine snow_ssp(f,qs,t,p,q,maxleg,kext, salb, back,  &
	nlegen, legen, legen2, legen3, legen4)
  
  use kinds
  use nml_params, only: verbose, lphase_flag, n_0snowDsnow, EM_snow
  use constants, only: pi, im

  implicit none

  integer :: nbins, nlegen
  integer, intent(in) :: maxleg

  real(kind=dbl), intent(in) :: &
    qs,&
    t,&
    p,&
    q,&
    f

  real(kind=dbl) :: refre, refim

  real(kind=dbl) :: dia1, dia2, swc, ad, bd, alpha, gamma, b_snow, a_msnow

  real(kind=dbl), intent(out) :: &
    kext,&
    salb,&
    back

  real(kind=dbl), dimension(200), intent(out) :: legen, legen2, legen3, legen4

  complex(kind=dbl) :: mindex, m_air

  real(kind=dbl) :: spec2abs, gammln
  

  if (verbose .gt. 1) print*, 'Entering snow_ssp'

!	b_snow = 2.0d0
!	a_msnow = 0.038d0
	b_snow = 2.2850d0
	a_msnow = 0.2124d0


	call ref_ice(t, f, refre, refim)
	mindex = refre-Im*refim  ! mimicking a

	m_air = 1.0d0 - 0.0d0 * Im 

	dia1 = 0.51d-5 ! minimum maximum diameter [m] after kneifel
	dia2 = 1.d-2 ! maximum maximum diameter [m] after kneifel

	swc =  spec2abs(qs,t,p,q) ! [kg/m^3]
	! Field param. ! multiplied by 10^6 is 1/m^4
	ad = n_0snowDsnow * 1.d6 * exp(-0.107d0 * (t - 273.15))

	ad = 1.d6

	bd = (exp(gammln(b_snow + 1)) * a_msnow * ad/swc)**(1.0d0/(1.0d0 + b_snow))  ! [m**-1]   

	nbins = 100

	if (EM_snow .eq. 'icesf') then 
	  call mie_densitysizedep_spheremasseq(f, mindex,      &
		a_msnow, b_snow, dia1, dia2, nbins, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
		LEGEN4, 'C')
	elseif (EM_snow .eq. 'surus') then 
	  call mie_icefactor(f, t,mindex,      &
		a_msnow, b_snow, dia1, dia2, nbins, maxleg,   &
		ad, bd, alpha, gamma, lphase_flag, kext, salb,      &
		back, NLEGEN, LEGEN, LEGEN2, LEGEN3,        &
		LEGEN4, 'C',0.863*1.e-3*f+0.115,42)
	elseif (EM_snow(1:3) .eq. 'liu') then
	    call dda_db_liu(f,t,9,mindex, &
		dia1,dia2,nbins,maxleg,ad,&
		bd, alpha, gamma, lphase_flag,kext, salb,&
		back, nlegen, legen, legen2, legen3,&
		legen4, 'C')
	else 
	    write (*, *) 'no em mod', EM_snow
	    stop
	endif

  if (verbose .gt. 1) print*, 'Exiting snow_ssp'

  return

end subroutine snow_ssp
