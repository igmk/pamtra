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
  real(kind=dbl), allocatable, dimension(:,:,:) :: radar_hgt
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Ze
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Att_atmo
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: Att_hydro
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: radar_spectra
  real(kind=dbl), allocatable, dimension(:,:,:,:) ::    radar_snr
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) ::    radar_moments
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) ::    radar_slopes
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) ::    radar_edge
  integer, allocatable, dimension(:,:,:,:) ::    radar_quality
  real(kind=dbl), allocatable, dimension(:) :: radar_vel

  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: psd_d_bound
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: psd_f
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: psd_mass
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: psd_area


end module vars_output

