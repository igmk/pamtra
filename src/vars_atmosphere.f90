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


  real(kind=sgl) :: cwp,iwp,rwp,swp,gwp,hwp


  integer(kind=long) :: atmo_ngridx, atmo_ngridy
  character(len=2), allocatable, dimension(:,:) :: atmo_month, atmo_day
  character(len=4), allocatable, dimension(:,:) :: atmo_year, atmo_time
  real(kind=sgl), allocatable, dimension(:,:) :: atmo_deltax, atmo_deltay
  integer(kind=long), allocatable, dimension(:,:) :: atmo_model_i, atmo_model_j
  real(kind=sgl), allocatable, dimension(:,:) :: atmo_lon,atmo_lat,atmo_lfrac,atmo_wind10u,atmo_wind10v
  !number of layer can vary (from python)
  integer(kind=long), allocatable, dimension(:,:) :: atmo_nlyrs, atmo_unixtime
  integer(kind=long) :: atmo_max_nlyrs

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
       atmo_airturb, &
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
    call assert_true(err,(atmo_max_nlyrs>0),&
        "atmo_max_nlyrs must be greater zero")   
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
    allocate(atmo_unixtime(atmo_ngridx,atmo_ngridy))
    allocate(atmo_groundtemp(atmo_ngridx,atmo_ngridy))
    allocate(atmo_iwv(atmo_ngridx,atmo_ngridy))

    allocate(atmo_relhum_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs+1))
    allocate(atmo_press_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs+1))
    allocate(atmo_temp_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs+1))
    allocate(atmo_hgt_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs+1))

    allocate(atmo_relhum(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs))
    allocate(atmo_press(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs))
    allocate(atmo_temp(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs))

    allocate(atmo_vapor_pressure(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs))
    allocate(atmo_rho_vap(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs))
    allocate(atmo_q_hum(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs))

    allocate(atmo_hgt(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs))
    allocate(atmo_delta_hgt_lev(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs))
    allocate(atmo_airturb(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs))


    allocate(atmo_hydro_q(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs,n_hydro))
    allocate(atmo_hydro_reff(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs,n_hydro))
    allocate(atmo_hydro_n(atmo_ngridx,atmo_ngridy,atmo_max_nlyrs,n_hydro))


    atmo_month(:,:) = "na"
    atmo_day(:,:) = "na"
    atmo_year(:,:) = "nan"
    atmo_time(:,:) = "nan"
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
    atmo_airturb(:,:,:) = 0.d0

    atmo_hydro_q(:,:,:,:) = nan()
    atmo_hydro_reff(:,:,:,:) = nan()
    atmo_hydro_n(:,:,:,:) = nan()

    atmo_nlyrs(:,:) = -9999
    atmo_unixtime(:,:) = -9999

    errorstatus = err
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    return

  end subroutine allocate_atmosphere_vars

  subroutine deallocate_atmosphere_vars()

    if (allocated(atmo_nlyrs)) deallocate(atmo_nlyrs)
    if (allocated(atmo_unixtime)) deallocate(atmo_unixtime)

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
    if (allocated(atmo_airturb)) deallocate(atmo_airturb)

    if (allocated(atmo_hydro_q)) deallocate(atmo_hydro_q)
    if (allocated(atmo_hydro_reff)) deallocate(atmo_hydro_reff)
    if (allocated(atmo_hydro_n)) deallocate(atmo_hydro_n)
 
    if (allocated(atmo_groundtemp)) deallocate(atmo_groundtemp)
    if (allocated(atmo_iwv)) deallocate(atmo_iwv)

    if (allocated(atmo_month)) deallocate(atmo_month)
    if (allocated(atmo_day)) deallocate(atmo_day)
    if (allocated(atmo_year)) deallocate(atmo_year)
    if (allocated(atmo_time)) deallocate(atmo_time)
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

  subroutine fillMissing_atmosphere_vars(errorstatus)

    integer :: nz, nx, ny
    real(kind=dbl) :: e_sat_gg_water
    integer,dimension(9) :: timestamp

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'fillMissing_atmosphere_vars' 

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    call assert_true(err,(atmo_ngridx>0),&
        "atmo_ngridx must be greater zero")   
    call assert_true(err,(atmo_ngridy>0),&
        "atmo_ngridy must be greater zero")   
    call assert_true(err,(all(atmo_nlyrs>0)),&
        "atmo_nlyrs must be greater zero")   
    if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
    end if  

    do nx = 1, atmo_ngridx
      do ny = 1, atmo_ngridy
          if (verbose >= 5) print*,  "processing", nx, ny, nz

          !these one depend not on z:
          if (atmo_unixtime(nx,ny) /= -9999) then
                call GMTIME(atmo_unixtime(nx,ny),timestamp)
                write(atmo_year(nx,ny),"(i4.4)") timestamp(6)+1900
                write(atmo_month(nx,ny),"(i2.2)") timestamp(5)+1
                write(atmo_day(nx,ny),"(i2.2)") timestamp(4)
                write(atmo_time(nx,ny)(1:2),"(i2.2)") timestamp(3)
                write(atmo_time(nx,ny)(3:4),"(i2.2)") timestamp(2)
          end if

          ! make the lev vairables which are needed
          ! not needed anywhere as of today are press_lev and relhum_lev:

          nz = 1
          if (isnan(atmo_temp_lev(nx,ny,1))) &
            atmo_temp_lev(nx,ny,1) = atmo_temp(nx,ny,1) - &
            (atmo_temp(nx,ny,2) - atmo_temp(nx,ny,1))*0.25
          if (isnan(atmo_hgt_lev(nx,ny,1))) &
            atmo_hgt_lev(nx,ny,1) = atmo_hgt(nx,ny,1) - &
            (atmo_hgt(nx,ny,2) - atmo_hgt(nx,ny,1))*0.25

          do nz = 2, atmo_nlyrs(nx,ny) 
           if (isnan(atmo_temp_lev(nx,ny,nz))) &
              atmo_temp_lev(nx,ny,nz) = atmo_temp(nx,ny,nz-1) + &
              (atmo_temp(nx,ny,nz) - atmo_temp(nx,ny,nz-1))/2
            if (isnan(atmo_hgt_lev(nx,ny,nz))) &
              atmo_hgt_lev(nx,ny,nz) = atmo_hgt(nx,ny,nz-1) + &
              (atmo_hgt(nx,ny,nz) - atmo_hgt(nx,ny,nz-1))/2
          end do 
          nz = atmo_nlyrs(nx,ny) +1
          if (isnan(atmo_temp_lev(nx,ny,nz))) &
            atmo_temp_lev(nx,ny,nz) = atmo_temp(nx,ny,nz-1) + &
            (atmo_temp(nx,ny,nz-1) - atmo_temp(nx,ny,nz-2))*0.25

          if (isnan(atmo_hgt_lev(nx,ny,nz))) &
            atmo_hgt_lev(nx,ny,nz) = atmo_hgt(nx,ny,nz-1) + &
            (atmo_hgt(nx,ny,nz-1) - atmo_hgt(nx,ny,nz-2))*0.25
          
          !now teh non-lev variables
          do nz = 1, atmo_nlyrs(nx,ny)

          ! first test the non lev variables for nan and calculate them otherwise
          if (isnan(atmo_temp(nx,ny,nz))) &
                atmo_temp(nx,ny,nz) = 0.5 * (atmo_temp_lev(nx,ny,nz) + atmo_temp_lev(nx,ny,nz+1)) 
          if (isnan(atmo_relhum(nx,ny,nz))) &
                atmo_relhum(nx,ny,nz) = 0.5 * (atmo_relhum_lev(nx,ny,nz) + atmo_relhum_lev(nx,ny,nz+1))
          if (isnan(atmo_press(nx,ny,nz))) &
                atmo_press(nx,ny,nz) = (atmo_press_lev(nx,ny,nz+1) - atmo_press_lev(nx,ny,nz))/ &
                log(atmo_press_lev(nx,ny,nz+1) / atmo_press_lev(nx,ny,nz))  ! [Pa]
          if (isnan(atmo_hgt(nx,ny,nz))) &
                atmo_hgt(nx,ny,nz) = (atmo_hgt_lev(nx,ny,nz)+atmo_hgt_lev(nx,ny,nz+1))*0.5d0
          if (isnan(atmo_delta_hgt_lev(nx,ny,nz))) &
                atmo_delta_hgt_lev(nx,ny,nz) = atmo_hgt_lev(nx,ny,nz+1) - atmo_hgt_lev(nx,ny,nz)
         ! now the ones which depend on several
          if (isnan(atmo_vapor_pressure(nx,ny,nz))) &
                atmo_vapor_pressure(nx,ny,nz) = atmo_relhum(nx,ny,nz) * &
                e_sat_gg_water(atmo_temp(nx,ny,nz)) ! Pa
          if (isnan(atmo_rho_vap(nx,ny,nz))) &
                atmo_rho_vap(nx,ny,nz) = atmo_vapor_pressure(nx,ny,nz)/(atmo_temp(nx,ny,nz) * 461.5)  ! [kg/m3]
          if (isnan(atmo_q_hum(nx,ny,nz))) &
                atmo_q_hum(nx,ny,nz) = 0.622*atmo_vapor_pressure(nx,ny,nz)/&
                (atmo_press(nx,ny,nz) - (1.- 0.622) * atmo_vapor_pressure(nx,ny,nz))  ! [kg/kg]

        end do
    !for ground temp, simply use the lowes avalable one if not provided
    if (isnan(atmo_groundtemp(nx,ny))) &
        atmo_groundtemp(nx,ny) = atmo_temp_lev(nx,ny,nz)

    !test whether we still have nans in our data!
    call assert_false(err,ANY(ISNAN(atmo_q_hum(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_q_hum")   
    call assert_false(err,ANY(ISNAN(atmo_rho_vap(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_rho_vap")   
    call assert_false(err,ANY(ISNAN(atmo_vapor_pressure(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_vapor_pressure")   
    call assert_false(err,ANY(ISNAN(atmo_delta_hgt_lev(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_delta_hgt_lev")   
    call assert_false(err,ANY(ISNAN(atmo_airturb(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_airturb")   
    call assert_false(err,ANY(ISNAN(atmo_hgt(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_hgt")   
    call assert_false(err,ANY(ISNAN(atmo_press(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_press")   
    call assert_false(err,ANY(ISNAN(atmo_relhum(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_relhum")   
    call assert_false(err,ANY(ISNAN(atmo_temp_lev(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_temp_lev")   
    call assert_false(err,ANY(ISNAN(atmo_hgt_lev(nx,ny,1:atmo_nlyrs(nx,ny)))),&
        "found nan in atmo_hgt_lev")   
!     call assert_false(err,ANY(ISNAN(atmo_relhum_lev(nx,ny,1:atmo_nlyrs(nx,ny)))),&
!         "found nan in atmo_relhum_lev")   
!     call assert_false(err,ANY(ISNAN(atmo_press_lev(nx,ny,1:atmo_nlyrs(nx,ny)))),&
!         "found nan in atmo_press_lev")   
      end do
    end do

call assert_false(err,ANY(ISNAN(atmo_groundtemp(:,:))),&
        "found nan in atmo_groundtemp")

    if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
    end if  


    errorstatus = err
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    return

  end subroutine fillMissing_atmosphere_vars

  subroutine print_vars_atmosphere()
    !for debuging prurposes


    print*, "atmo_nlyrs", atmo_nlyrs
    print*, "atmo_unixtime", atmo_unixtime

    print*, "atmo_relhum_lev", atmo_relhum_lev
    print*, "atmo_press_lev", atmo_press_lev
    print*, "atmo_temp_lev", atmo_temp_lev
    print*, "atmo_hgt_lev", atmo_hgt_lev

    print*, "atmo_relhum", atmo_relhum
    print*, "atmo_press", atmo_press
    print*, "atmo_vapor_pressure", atmo_vapor_pressure
    print*, "atmo_rho_vap", atmo_rho_vap
    print*, "atmo_q_hum", atmo_q_hum
    print*, "atmo_temp", atmo_temp
    print*, "atmo_hgt", atmo_hgt
    print*, "atmo_delta_hgt_lev", atmo_delta_hgt_lev
    print*, "atmo_airturb", atmo_airturb

    print*, "atmo_hydro_q", atmo_hydro_q
    print*, "atmo_hydro_reff", atmo_hydro_reff
    print*, "atmo_hydro_n", atmo_hydro_n
 
    print*, "atmo_groundtemp", atmo_groundtemp
    print*, "atmo_iwv", atmo_iwv

    print*, "atmo_month", atmo_month
    print*, "atmo_day", atmo_day
    print*, "atmo_year", atmo_year
    print*, "atmo_time", atmo_time
    print*, "atmo_deltax", atmo_deltax
    print*, "atmo_deltay", atmo_deltay
    print*, "atmo_model_i", atmo_model_i
    print*, "atmo_model_j", atmo_model_j
    print*, "atmo_lon", atmo_lon
    print*, "atmo_lat", atmo_lat
    print*, "atmo_lfrac", atmo_lfrac
    print*, "atmo_wind10u", atmo_wind10u
    print*, "atmo_wind10v", atmo_wind10v

  end subroutine print_vars_atmosphere

end module vars_atmosphere
