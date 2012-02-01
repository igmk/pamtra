module vars_output

  use kinds
  implicit none
  save

  integer, allocatable, dimension(:,:) :: is, js
  real(kind=sgl), allocatable, dimension(:,:) :: lons, lats, lfracs, iwvs, cwps, &
       iwps, rwps, swps, gwps,hwps

  !for passive
  real(kind=dbl), allocatable, dimension(:,:,:,:,:,:) :: tb
  real(kind=dbl), dimension(32) :: angles_deg !2*NUMMU=32

  !for active 
  real(kind=dbl), allocatable, dimension(:,:,:) :: hgt
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Ze,&
                               Ze_cw,Ze_rr,Ze_ci,Ze_sn,Ze_gr,Ze_ha
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Att_atmo
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Att_hydro, &
                               Att_cw,Att_rr,Att_ci,Att_sn,Att_gr,Att_ha
  
end module vars_output

