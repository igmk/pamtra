module vars_atmosphere

use kinds
use profile_type

implicit none
save

integer :: ngridx, ngridy, nlyr
			  
real(kind=dbl), allocatable, dimension(:) :: relhum_lev,&
                                          press_lev, &
                                          temp_lev, &
					  hgt_lev

real(kind=dbl), allocatable, dimension(:) :: press, &
					     temp,&
					     relhum,&
                                             vapor_pressure, &
                                             rho_vap, &
					     q_hum

!real(kind=dbl), dimension(:), allocatable :: lyr_temp, &
!                                      lyr_pres

!real(kind=dbl), dimension(:), allocatable :: rel_hum, &
!                                      avgpressure, &
!                                      vaporpressure

real(kind=dbl), allocatable, dimension(:) :: cwc_q, &
                                          iwc_q, &
                                          rwc_q, &
                                          swc_q, &
                                          gwc_q

! real, allocatable, dimension(:,:) :: lon, &
!                                      lat, &
!                                      lfrac, &
!                                      wind10

! real, allocatable, dimension(:,:) :: iwv,&
!                                      cwp,&
!                                      iwp,&
!                                      rwp,&
!                                      swp,&
!                                      gwp

  integer :: alloc_status

  type(profile), allocatable :: profiles(:,:)



end module vars_atmosphere