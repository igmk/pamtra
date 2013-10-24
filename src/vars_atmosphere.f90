module vars_atmosphere

  use kinds
  use settings, only: maxfreq

  implicit none
  save

!   integer(kind=long) :: nlyr
  integer(kind=long) :: ngridx, ngridy
  character(2) :: month, day
  character(4) :: year, time
  character(12) :: date_str
  real(kind=sgl) :: deltax, deltay

  integer(kind=long), allocatable, dimension(:) :: rt_nlegen

  !is allocated in pamtra.f90!

!   real(kind=dbl), allocatable, dimension(:) :: relhum_lev,&
!        press_lev, &
!        temp_lev, &
!        hgt_lev

!   real(kind=dbl), allocatable, dimension(:) :: hgt
! press, &
!        temp,&
!        relhum,&
!        vapor_pressure, &
!        rho_vap, &
!        q_hum,&
       

!   real(kind=dbl), allocatable, dimension(:) :: cwc_q, &
!        iwc_q, &
!        rwc_q, &
!        swc_q, &
!        gwc_q, &
!        hwc_q
! 
!   real(kind=dbl), allocatable, dimension(:) :: cwc_n, &
!        iwc_n, &
!        rwc_n, &
!        swc_n, &
!        gwc_n, &
!        hwc_n

!   real(kind=dbl), allocatable, dimension(:) :: delta_hgt_lev


  integer(kind=long) :: alloc_status

  integer(kind=long) :: model_i, model_j
  real(kind=sgl) :: lon,lat,lfrac,wind10u,wind10v,iwv,cwp,iwp,rwp,swp,gwp,hwp

  integer(kind=long), allocatable, dimension(:,:) :: atmo_nlyrs
  integer(kind=long) :: atmo_max_nlyr

  real(kind=dbl), allocatable, dimension(:,:) :: atmo_groundtemp

  real(kind=dbl), allocatable, dimension(:,:,:) :: atmo_relhum_lev,&
       atmo_press_lev, &
       atmo_temp_lev, &
       atmo_hgt_lev

  real(kind=dbl), allocatable, dimension(:,:,:) :: atmo_relhum,&
       atmo_press, &
       atmo_temp, &
       atmo_hgt, &
       atmo_delta_hgt_lev, &
       atmo_vapor_pressure, &
       atmo_rho_vap, &
       atmo_q_hum 

  real(kind=dbl), allocatable, dimension(:,:,:,:) :: atmo_hydro_q,&
       atmo_hydro_reff, &
       atmo_hydro_n

  contains

  subroutine allocate_atmosphere_vars()


    use descriptor_file, only: n_hydro
    real(kind=dbl) :: nan

    !todo: add assert ngrid_X > 0 etc etc
    allocate(atmo_nlyrs(ngridx,ngridy))

    allocate(atmo_relhum_lev(ngridx,ngridy,atmo_max_nlyr+1))
    allocate(atmo_press_lev(ngridx,ngridy,atmo_max_nlyr+1))
    allocate(atmo_temp_lev(ngridx,ngridy,atmo_max_nlyr+1))
    allocate(atmo_hgt_lev(ngridx,ngridy,atmo_max_nlyr+1))

    allocate(atmo_relhum(ngridx,ngridy,atmo_max_nlyr))
    allocate(atmo_press(ngridx,ngridy,atmo_max_nlyr))
    allocate(atmo_temp(ngridx,ngridy,atmo_max_nlyr))

    allocate(atmo_vapor_pressure(ngridx,ngridy,atmo_max_nlyr))
    allocate(atmo_rho_vap(ngridx,ngridy,atmo_max_nlyr))
    allocate(atmo_q_hum(ngridx,ngridy,atmo_max_nlyr))

    allocate(atmo_hgt(ngridx,ngridy,atmo_max_nlyr))
    allocate(atmo_delta_hgt_lev(ngridx,ngridy,atmo_max_nlyr))


    allocate(atmo_hydro_q(ngridx,ngridy,atmo_max_nlyr,n_hydro))
    allocate(atmo_hydro_reff(ngridx,ngridy,atmo_max_nlyr,n_hydro))
    allocate(atmo_hydro_n(ngridx,ngridy,atmo_max_nlyr,n_hydro))


    atmo_relhum_lev(:,:,:) = nan()
    atmo_press_lev(:,:,:) = nan()
    atmo_temp_lev(:,:,:) = nan()
    atmo_hgt_lev(:,:,:) = nan()

    atmo_relhum(:,:,:) = nan()
    atmo_press(:,:,:) = nan()
    atmo_temp(:,:,:) = nan()
    atmo_vapor_pressure(:,:,:) = nan()
    atmo_rho_vap(:,:,:) = nan()
    atmo_q_hum(:,:,:) = nan()
    atmo_hgt(:,:,:) = nan()
    atmo_delta_hgt_lev(:,:,:) = nan()

    atmo_hydro_q(:,:,:,:) = nan()
    atmo_hydro_reff(:,:,:,:) = nan()
    atmo_hydro_n(:,:,:,:) = nan()

  end subroutine allocate_atmosphere_vars

  subroutine deallocate_atmosphere_vars()

    if (allocated(atmo_nlyrs)) deallocate(atmo_nlyrs)

    if (allocated(atmo_relhum_lev)) deallocate(atmo_relhum_lev)
    if (allocated(atmo_press_lev)) deallocate(atmo_press_lev)
    if (allocated(atmo_temp_lev)) deallocate(atmo_temp_lev)
    if (allocated(atmo_hgt_lev)) deallocate(atmo_hgt_lev)

    if (allocated(atmo_relhum)) deallocate(atmo_relhum)
    if (allocated(atmo_press)) deallocate(atmo_press)
    if (allocated(atmo_vapor_pressure)) deallocate(atmo_vapor_pressure)
    if (allocated(atmo_rho_vap)) deallocate(atmo_rho_vap)
    if (allocated(atmo_q_hum)) deallocate(atmo_q_hum)
    if (allocated(atmo_temp)) deallocate(atmo_temp)
    if (allocated(atmo_hgt)) deallocate(atmo_hgt)
    if (allocated(atmo_delta_hgt_lev)) deallocate(atmo_delta_hgt_lev)

    if (allocated(atmo_hydro_q)) deallocate(atmo_hydro_q)
    if (allocated(atmo_hydro_reff)) deallocate(atmo_hydro_reff)
    if (allocated(atmo_hydro_n)) deallocate(atmo_hydro_n)

  end subroutine deallocate_atmosphere_vars

end module vars_atmosphere
