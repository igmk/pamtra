module vars_output

  use kinds
  use nml_params, only: NUMMU
  implicit none
  save

  integer, allocatable, dimension(:,:) :: is, js
  real(kind=sgl), allocatable, dimension(:,:) :: lons, lats, lfracs, iwvs, cwps, &
       iwps, rwps, swps, gwps,hwps

  !for passive
  real(kind=dbl), allocatable, dimension(:,:,:,:,:,:) :: tb
  real(kind=dbl), dimension(2*NUMMU) :: angles_deg

  !for active 
  real(kind=dbl), allocatable, dimension(:,:,:) :: hgt
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Ze, Attenuation_atmo, Attenuation_hydro!, &

end module vars_output
