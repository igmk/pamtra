subroutine allocate_vars_atmosphere

use vars_atmosphere
use nml_params

implicit none

integer :: i,j


      allocate(hgt_lev(0:nlyr),stat=alloc_status)
      allocate(press_lev(0:nlyr),stat=alloc_status)
      allocate(press(nlyr),stat=alloc_status)
      allocate(temp_lev(0:nlyr),stat=alloc_status)
      allocate(temp(nlyr),stat=alloc_status)
      allocate(relhum_lev(0:nlyr),stat=alloc_status)
      allocate(relhum(nlyr),stat=alloc_status)

      allocate(vapor_pressure(nlyr),stat=alloc_status)
      allocate(rho_vap(nlyr),stat=alloc_status)
      allocate(q_hum(nlyr),stat=alloc_status)

      allocate(cwc_q(nlyr))
      allocate(iwc_q(nlyr))
      allocate(rwc_q(nlyr))
      allocate(swc_q(nlyr))
      allocate(gwc_q(nlyr))
      if (n_moments .eq. 2) then
        allocate(hwc_q(nlyr))
        allocate(cwc_n(nlyr))
        allocate(iwc_n(nlyr))
        allocate(rwc_n(nlyr))
        allocate(swc_n(nlyr))
        allocate(gwc_n(nlyr))
        allocate(hwc_n(nlyr))
	  end if

  allocate(nlegen(nlyr),stat=alloc_status)
  allocate(rt3nlegen(nlyr),stat=alloc_status)
  allocate(kextatmo(nlyr), stat=alloc_status)
  allocate(kexttot(nlyr), stat=alloc_status)
  allocate(salbtot(nlyr), stat=alloc_status)
  allocate(rt3kexttot(nlyr), stat=alloc_status)
  allocate(rt3salbtot(nlyr), stat=alloc_status)
  allocate(g_coeff(nlyr), stat=alloc_status)
  allocate(back(nlyr), stat=alloc_status)
  allocate(legen(nlyr,200), stat=alloc_status)
  allocate(legen2(nlyr,200), stat=alloc_status)
  allocate(legen3(nlyr,200), stat=alloc_status)
  allocate(legen4(nlyr,200), stat=alloc_status)
  allocate(rt3legen(nlyr,200), stat=alloc_status)
  allocate(rt3legen2(nlyr,200), stat=alloc_status)
  allocate(rt3legen3(nlyr,200), stat=alloc_status)
  allocate(rt3legen4(nlyr,200), stat=alloc_status)

  allocate(profiles(ngridx,ngridy),stat=alloc_status)
  do i = 1, ngridx
     do j = 1, ngridy
	allocate(profiles(i,j)%hgt_lev(0:nlyr), stat=alloc_status)
	allocate(profiles(i,j)%press_lev(0:nlyr), stat=alloc_status)
	allocate(profiles(i,j)%press(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%temp_lev(0:nlyr), stat=alloc_status)
	allocate(profiles(i,j)%temp(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%relhum_lev(0:nlyr), stat=alloc_status)
	allocate(profiles(i,j)%relhum(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%cloud_water_q(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%cloud_ice_q(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%rain_q(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%snow_q(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%graupel_q(nlyr), stat=alloc_status)
	if (n_moments .eq. 2) then
	  allocate(profiles(i,j)%hail_q(nlyr), stat=alloc_status)
	  allocate(profiles(i,j)%cloud_water_n(nlyr), stat=alloc_status)
	  allocate(profiles(i,j)%cloud_ice_n(nlyr), stat=alloc_status)
	  allocate(profiles(i,j)%rain_n(nlyr), stat=alloc_status)
	  allocate(profiles(i,j)%snow_n(nlyr), stat=alloc_status)
	  allocate(profiles(i,j)%graupel_n(nlyr), stat=alloc_status)
	  allocate(profiles(i,j)%hail_n(nlyr), stat=alloc_status)
	end if
	allocate(profiles(i,j)%vapor_pressure(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%rho_vap(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%q_hum(nlyr), stat=alloc_status)
	end do
  end do

end subroutine allocate_vars_atmosphere
