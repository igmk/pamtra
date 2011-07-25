module vars_atmosphere

use kinds
use profile_type

implicit none
save

integer :: nlyr

integer, allocatable, dimension(:) :: nlegen
			  
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

real(kind=dbl), allocatable, dimension(:) :: cwc_q, &
                              	             iwc_q, &
                              	             rwc_q, &
                               	             swc_q, &
                                             gwc_q

 real(kind=dbl), allocatable, dimension(:) :: kextatmo, &
                                  		      kexttot, &
                                      		  salbtot, &
                                      		  back, &
                                      		  g_coeff

 real(kind=dbl), allocatable, dimension(:,:) :: legen, &
												legen2, &
												legen3, &
												legen4

  integer :: alloc_status

  type(profile), allocatable :: profiles(:,:)

end module vars_atmosphere
