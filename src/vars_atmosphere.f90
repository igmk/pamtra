module vars_atmosphere

  use kinds
  use settings, only: maxfreq
  use report_module

  implicit none
  save

!   integer(kind=long) :: nlyr


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



!   integer(kind=long) :: ngridx, ngridy
!   character(2) :: month, day
!   character(4) :: year, time
!   real(kind=sgl) :: deltax, deltay
!   integer(kind=long) :: model_i, model_j
!   real(kind=sgl) :: lon,lat,lfrac,wind10u,wind10v


  real(kind=sgl) :: iwv,cwp,iwp,rwp,swp,gwp,hwp


  integer(kind=long) :: atmo_ngridx, atmo_ngridy
  character(len=2), allocatable, dimension(:,:) :: atmo_month, atmo_day
  character(len=4), allocatable, dimension(:,:) :: atmo_year, atmo_time
  character(len=12), allocatable, dimension(:,:) :: atmo_date_str
  real(kind=sgl), allocatable, dimension(:,:) :: atmo_deltax, atmo_deltay
  integer(kind=long), allocatable, dimension(:,:) :: atmo_model_i, atmo_model_j
  real(kind=sgl), allocatable, dimension(:,:) :: atmo_lon,atmo_lat,atmo_lfrac,atmo_wind10u,atmo_wind10v
  !number of layer can vary (from python)
  integer(kind=long), allocatable, dimension(:,:) :: atmo_nlyrs
  integer(kind=long) :: atmo_max_nlyr

  real(kind=dbl), allocatable, dimension(:,:) :: atmo_groundtemp, atmo_iwv

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

  subroutine allocate_atmosphere_vars(errorstatus)


    use descriptor_file, only: n_hydro
    real(kind=dbl) :: nan

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'allocate_atmosphere_vars' 

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    call assert_true(err,(atmo_ngridx>0),&
        "atmo_ngridx must be greater zero")   
    call assert_true(err,(atmo_ngridy>0),&
        "atmo_ngridy must be greater zero")   
    call assert_true(err,(atmo_max_nlyr>0),&
        "atmo_max_nlyr must be greater zero")   
    if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
    end if  

    allocate(atmo_month(atmo_ngridx,atmo_ngridy))
    allocate(atmo_day(atmo_ngridx,atmo_ngridy))
    allocate(atmo_year(atmo_ngridx,atmo_ngridy))
    allocate(atmo_time(atmo_ngridx,atmo_ngridy))
    allocate(atmo_date_str(atmo_ngridx,atmo_ngridy))
    allocate(atmo_deltax(atmo_ngridx,atmo_ngridy))
    allocate(atmo_deltay(atmo_ngridx,atmo_ngridy))
    allocate(atmo_model_i(atmo_ngridx,atmo_ngridy))
    allocate(atmo_model_j(atmo_ngridx,atmo_ngridy))
    allocate(atmo_lon(atmo_ngridx,atmo_ngridy))
    allocate(atmo_lat(atmo_ngridx,atmo_ngridy))
    allocate(atmo_lfrac(atmo_ngridx,atmo_ngridy))
    allocate(atmo_wind10u(atmo_ngridx,atmo_ngridy))
    allocate(atmo_wind10v(atmo_ngridx,atmo_ngridy))

    !todo: add assert ngrid_X > 0 etc etc
    allocate(atmo_nlyrs(atmo_ngridx,atmo_ngridy))
    allocate(atmo_groundtemp(atmo_ngridx,atmo_ngridy))
    allocate(atmo_iwv(atmo_ngridx,atmo_ngridy))

    allocate(atmo_relhum_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyr+1))
    allocate(atmo_press_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyr+1))
    allocate(atmo_temp_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyr+1))
    allocate(atmo_hgt_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyr+1))

    allocate(atmo_relhum(atmo_ngridx,atmo_ngridy,atmo_max_nlyr))
    allocate(atmo_press(atmo_ngridx,atmo_ngridy,atmo_max_nlyr))
    allocate(atmo_temp(atmo_ngridx,atmo_ngridy,atmo_max_nlyr))

    allocate(atmo_vapor_pressure(atmo_ngridx,atmo_ngridy,atmo_max_nlyr))
    allocate(atmo_rho_vap(atmo_ngridx,atmo_ngridy,atmo_max_nlyr))
    allocate(atmo_q_hum(atmo_ngridx,atmo_ngridy,atmo_max_nlyr))

    allocate(atmo_hgt(atmo_ngridx,atmo_ngridy,atmo_max_nlyr))
    allocate(atmo_delta_hgt_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyr))


    allocate(atmo_hydro_q(atmo_ngridx,atmo_ngridy,atmo_max_nlyr,n_hydro))
    allocate(atmo_hydro_reff(atmo_ngridx,atmo_ngridy,atmo_max_nlyr,n_hydro))
    allocate(atmo_hydro_n(atmo_ngridx,atmo_ngridy,atmo_max_nlyr,n_hydro))


    atmo_month(:,:) = "na"
    atmo_day(:,:) = "na"
    atmo_year(:,:) = "nan"
    atmo_time(:,:) = "nan"
    atmo_date_str(:,:) = "nan"
    atmo_deltax(:,:) = nan()
    atmo_deltay(:,:) = nan()
    atmo_model_i(:,:) = -9999
    atmo_model_j(:,:) = -9999
    atmo_lon(:,:) = nan()
    atmo_lat(:,:) = nan()
    atmo_lfrac(:,:) = nan()
    atmo_wind10u(:,:) = nan()
    atmo_wind10v(:,:) = nan()

    atmo_relhum_lev(:,:,:) = nan()
    atmo_press_lev(:,:,:) = nan()
    atmo_temp_lev(:,:,:) = nan()
    atmo_hgt_lev(:,:,:) = nan()

    atmo_groundtemp(:,:) = nan()
    atmo_iwv(:,:) = nan()


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

    atmo_nlyrs(:,:) = -9999

    errorstatus = err
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    return

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
 
    if (allocated(atmo_groundtemp)) deallocate(atmo_groundtemp)
    if (allocated(atmo_iwv)) deallocate(atmo_iwv)

    if (allocated(atmo_month)) deallocate(atmo_month)
    if (allocated(atmo_day)) deallocate(atmo_day)
    if (allocated(atmo_year)) deallocate(atmo_year)
    if (allocated(atmo_time)) deallocate(atmo_time)
    if (allocated(atmo_date_str)) deallocate(atmo_date_str)
    if (allocated(atmo_deltax)) deallocate(atmo_deltax)
    if (allocated(atmo_deltay)) deallocate(atmo_deltay)
    if (allocated(atmo_model_i)) deallocate(atmo_model_i)
    if (allocated(atmo_model_j)) deallocate(atmo_model_j)
    if (allocated(atmo_lon)) deallocate(atmo_lon)
    if (allocated(atmo_lat)) deallocate(atmo_lat)
    if (allocated(atmo_lfrac)) deallocate(atmo_lfrac)
    if (allocated(atmo_wind10u)) deallocate(atmo_wind10u)
    if (allocated(atmo_wind10v)) deallocate(atmo_wind10v)



  end subroutine deallocate_atmosphere_vars

end module vars_atmosphere
