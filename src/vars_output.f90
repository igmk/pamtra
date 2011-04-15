module vars_output

use kinds

implicit none
save

integer, allocatable, dimension(:,:) :: is, js
real(kind=dbl), allocatable, dimension(:,:) :: lons, lats, lfracs, w10s, iwvs, cwps, iwps, rwps, swps, gwps

real(kind=dbl), allocatable, dimension(:,:,:,:) :: flux_up, flux_down
real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: tb_up, tb_down

end module vars_output