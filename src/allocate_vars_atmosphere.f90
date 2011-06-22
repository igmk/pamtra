subroutine allocate_vars_atmosphere

use vars_atmosphere
use nml_params

implicit none

integer :: i,j

      allocate(hgt_lev(0:nlyr))
      allocate(press_lev(0:nlyr))
      allocate(press(nlyr))
      allocate(temp_lev(0:nlyr))
      allocate(temp(nlyr))
      allocate(relhum_lev(0:nlyr))
      allocate(relhum(nlyr))

      allocate(vapor_pressure(nlyr))
      allocate(rho_vap(nlyr))
      allocate(q_hum(nlyr))

!      allocate(lyr_temp(0:nlyr))
!      allocate(lyr_pres(0:nlyr))
!      allocate(rel_hum(nlyr))
!      allocate(avgpressure(nlyr))
!      allocate(vaporpressure(nlyr))

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

!       allocate(iwv(ngridx, ngridy))
!       allocate(cwp(ngridx, ngridy))
!       allocate(iwp(ngridx, ngridy))
!       allocate(rwp(ngridx, ngridy))
!       allocate(swp(ngridx, ngridy))
!       allocate(gwp(ngridx, ngridy))

!      allocate(lon(ngridx, ngridy),lat(ngridx, ngridy))
!      allocate(lfrac(ngridx, ngridy),wind10(ngridx, ngridy))

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
	allocate(profiles(i,j)%kextatmo(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%kexttot(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%salbtot(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%g_coeff(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%back(nlyr), stat=alloc_status)
	allocate(profiles(i,j)%legen(nlyr,200), stat=alloc_status)
	allocate(profiles(i,j)%legen2(nlyr,200), stat=alloc_status)
	allocate(profiles(i,j)%legen3(nlyr,200), stat=alloc_status)
	allocate(profiles(i,j)%legen4(nlyr,200), stat=alloc_status)
     end do
  end do

!   do i = 1, ngridx
!      do j = 1, ngridy
! 	nullify(profiles(i,j)%hgt_lev)
! 	nullify(profiles(i,j)%press_lev)
! 	nullify(profiles(i,j)%temp_lev)
! 	nullify(profiles(i,j)%relhum_lev)
! 	nullify(profiles(i,j)%cloud_water_q)
! 	nullify(profiles(i,j)%cloud_ice_q)
! 	nullify(profiles(i,j)%rain_q)
! 	nullify(profiles(i,j)%snow_q)
! 	nullify(profiles(i,j)%graupel_q)
! 	nullify(profiles(i,j)%kextatmo)
! 	nullify(profiles(i,j)%kexttot)
! 	nullify(profiles(i,j)%salbtot)
! 	nullify(profiles(i,j)%g_coeff)
! 	nullify(profiles(i,j)%back)
! 	nullify(profiles(i,j)%legen)
! 	nullify(profiles(i,j)%legen2)
! 	nullify(profiles(i,j)%legen3)
! 	nullify(profiles(i,j)%legen4)
!      end do
!   end do

end subroutine allocate_vars_atmosphere
