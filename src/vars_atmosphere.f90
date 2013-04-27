module vars_atmosphere

  use kinds
  use settings, only: maxfreq

  implicit none
  save

  integer(kind=long) :: nlyr
  integer(kind=long) :: ngridx, ngridy
  character(2) :: month, day
  character(4) :: year, time
  character(12) :: date_str
  real(kind=sgl) :: deltax, deltay

  integer(kind=long), allocatable, dimension(:) :: nlegen

  !is allocated in pamtra.f90!

  real(kind=dbl), allocatable, dimension(:) :: relhum_lev,&
       press_lev, &
       temp_lev, &
       hgt_lev

  real(kind=dbl), allocatable, dimension(:) :: press, &
       temp,&
       relhum,&
       vapor_pressure, &
       rho_vap, &
       q_hum,&
       hgt

  real(kind=dbl), allocatable, dimension(:) :: cwc_q, &
       iwc_q, &
       rwc_q, &
       swc_q, &
       gwc_q, &
       hwc_q

  real(kind=dbl), allocatable, dimension(:) :: cwc_n, &
       iwc_n, &
       rwc_n, &
       swc_n, &
       gwc_n, &
       hwc_n

  real(kind=dbl), allocatable, dimension(:) :: delta_hgt_lev


  real(kind=dbl), allocatable, dimension(:) :: kextatmo, &
       kexttot, kextsn, kextcw, kextrr, kextgr, kextci, kextha, &
       salbtot, &
       back, backcw, backrr, backci, backsn, backgr, backha, &
       g_coeff,&
       rt4salbtot

  real(kind=dbl), allocatable, dimension(:,:) :: legen, &
       legen2, &
       legen3, &
       legen4

  real(kind=dbl), allocatable, dimension(:,:,:,:,:,:) :: rt4scatter_matrix,scattermatrix
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: rt4ext_matrix,extmatrix
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: rt4emis_vec,emisvec

  !for jacobain mode
  real(kind=dbl), allocatable, dimension(:,:,:,:,:,:) :: jac_scattermatrix
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: jac_extmatrix
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: jac_emisvec
  real(kind=dbl), allocatable, dimension(:) :: jac_kextsn
  real(kind=dbl), allocatable, dimension(:) :: jac_backsn
  real(kind=dbl), allocatable, dimension(:) :: jac_kextcw
  real(kind=dbl), allocatable, dimension(:) :: jac_backcw
  real(kind=dbl), allocatable, dimension(:) :: jac_kextrr
  real(kind=dbl), allocatable, dimension(:) :: jac_backrr
  real(kind=dbl), allocatable, dimension(:) :: jac_kextgr
  real(kind=dbl), allocatable, dimension(:) :: jac_backgr
  real(kind=dbl), allocatable, dimension(:) :: jac_kextci
  real(kind=dbl), allocatable, dimension(:) :: jac_backci
  real(kind=dbl), allocatable, dimension(:) :: jac_kextha
  real(kind=dbl), allocatable, dimension(:) :: jac_backha
  logical, allocatable, dimension(:) :: jac_hydros_present

  !jacobian mode
  real(kind=dbl), allocatable, dimension(:) :: jac_cwc_q, &
       jac_iwc_q, &
       jac_rwc_q, &
       jac_swc_q, &
       jac_gwc_q, &
       jac_hwc_q

  real(kind=dbl), allocatable, dimension(:) :: jac_cwc_n, &
       jac_iwc_n, &
       jac_rwc_n, &
       jac_swc_n, &
       jac_gwc_n, &
       jac_hwc_n

  real(kind=dbl), allocatable, dimension(:) :: jac_relhum_lev,&
       jac_temp_lev

  logical, allocatable, dimension(:) :: hydros_present, rt4hydros_present

  integer(kind=long) :: alloc_status

  integer(kind=long) :: model_i, model_j
  real(kind=sgl) :: lon,lat,lfrac,wind10u,wind10v,iwv,cwp,iwp,rwp,swp,gwp,hwp

end module vars_atmosphere
