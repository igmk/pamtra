module vars_output

use kinds

implicit none
save

integer, allocatable, dimension(:,:) :: is, js
real(kind=sgl), allocatable, dimension(:,:) :: lons, lats, lfracs, iwvs, cwps, &
						iwps, rwps, swps, gwps,hwps

real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: tb

!for active 
real(kind=dbl), allocatable, dimension(:,:,:) :: hgt, Ze, PIA_atmo_bottomup, PIA_hydro_bottomup, &
						PIA_atmo_topdown, PIA_hydro_topdown

end module vars_output
