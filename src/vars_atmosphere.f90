module vars_atmosphere

  use kinds
  implicit none
  save

  type profile
     integer :: isamp
     integer :: jsamp
     real(kind=sgl) :: latitude
     real(kind=sgl) :: longitude
     real(kind=sgl) :: land_fraction
     real(kind=sgl) :: wind_10u, wind_10v
     real(kind=sgl) :: iwv
     real(kind=sgl) :: cwp
     real(kind=sgl) :: iwp
     real(kind=sgl) :: rwp
     real(kind=sgl) :: swp
     real(kind=sgl) :: gwp
     real(kind=sgl) :: hwp
     !  real(kind=dbl) :: ground_albedo
     !  real(kind=dbl) :: ground_index
     !  real(kind=dbl) :: ground_temperature
     real(kind=dbl), pointer :: hgt_lev(:)
     real(kind=dbl), pointer :: press_lev(:)
     real(kind=dbl), pointer :: temp_lev(:)
     real(kind=dbl), pointer :: relhum_lev(:)

     real(kind=dbl), pointer :: cloud_water_q(:)
     real(kind=dbl), pointer :: cloud_ice_q(:)
     real(kind=dbl), pointer :: rain_q(:)
     real(kind=dbl), pointer :: snow_q(:)
     real(kind=dbl), pointer :: graupel_q(:)
     real(kind=dbl), pointer :: hail_q(:)

     real(kind=dbl), pointer :: cloud_water_n(:)
     real(kind=dbl), pointer :: cloud_ice_n(:)
     real(kind=dbl), pointer :: rain_n(:)
     real(kind=dbl), pointer :: snow_n(:)
     real(kind=dbl), pointer :: graupel_n(:)
     real(kind=dbl), pointer :: hail_n(:)

     real(kind=dbl), pointer :: press(:)
     real(kind=dbl), pointer :: temp(:)
     real(kind=dbl), pointer :: relhum(:)
     real(kind=dbl), pointer :: vapor_pressure(:)
     real(kind=dbl), pointer :: rho_vap(:)
     real(kind=dbl), pointer :: q_hum(:)
  end type profile



  integer :: nlyr, nfrq
  integer :: ngridx, ngridy

  integer, allocatable, dimension(:) :: nlegen
  integer, allocatable, dimension(:) :: rt3nlegen

  !is allocated in pamtra.f90!
  real(kind=dbl), allocatable, dimension(:) :: freqs

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
       gwc_q, &
       hwc_q

  real(kind=dbl), allocatable, dimension(:) ::	cwc_n, &
       iwc_n, &
       rwc_n, &
       swc_n, &
       gwc_n, &
       hwc_n

  real(kind=dbl), allocatable, dimension(:) :: kextatmo, &
       kexttot, &
       rt3kexttot,&
       salbtot, &
       rt3salbtot,&
       back, &
       g_coeff

  real(kind=dbl), allocatable, dimension(:,:) :: legen, &
       legen2, &
       legen3, &
       legen4

  real(kind=dbl), allocatable, dimension(:,:) :: rt3legen, &
       rt3legen2, &
       rt3legen3, &
       rt3legen4

  integer :: alloc_status

  type(profile), allocatable :: profiles(:,:)
  integer, allocatable, dimension(:,:) :: ics

  character(2) :: month, day
  character(4) :: year, time
  character(12) :: date_str
  real(kind=sgl) :: deltax, deltay
contains

  subroutine vars_atmosphere_read_profile(input_file)
    use nml_params
    integer :: i,j,k,istat
    character(99),intent(in) :: input_file

    if (verbose .gt. 0) print *,"opening: ",input_file


    !                                                                       
    !     read atmospheric profiles                 
    !  
    !     quantities followed  
    !     by "_lev" are given at the layer heights while the quantitites    
    !     w/o "_lev" are layer average quantities
    !
    !  Variables:
    !     name           unit        description
    !     hgt_lev         m
    !     press_lev       Pa
    !     temp_lev        K
    !     relhum_lev      %
    !
    !     cloud_water_q kg/kg
    !     cloud_ice_q   kg/kg 
    !     rain_q        kg/kg
    !     snow_q        kg/kg
    !     graupel_q     kg/kg


    open(UNIT=14, FILE=input_path(:len_trim(input_path))//"/"//input_file(:len_trim(input_file)),&
         STATUS='OLD', form='formatted',iostat=istat)

    if (istat .ne. 0) then
       call error_msg(input_file,0,0)
       stop "Read error 1: Cannot open file"
    end if

    read(14,*,iostat=istat) year, month, day, time, ngridx, ngridy, nlyr, deltax, deltay

    date_str = year//month//day//time

    if (istat .ne. 0) then
       call error_msg(input_file,0,0)
       stop "Read error 2: Cannot read first line of file"
    end if



    allocate(profiles(ngridx,ngridy),stat=alloc_status)
    do i = 1, ngridx
       do j = 1, ngridy
          allocate(profiles(i,j)%hgt_lev(0:nlyr), stat=alloc_status)
          allocate(profiles(i,j)%press_lev(0:nlyr), stat=alloc_status)
          allocate(profiles(i,j)%press(nlyr), stat=alloc_status)
          allocate(profiles(i,j)%temp_lev(0:nlyr), stat=alloc_status)
          allocate(profiles(i,j)%temp(nlyr), stat=alloc_status)
          allocate(profiles(i,j)%relhum_lev(0:nlyr), stat=alloc_status)
          allocate(profiles(i,j)%relhum(nlyr), stat=alloc_status)
          allocate(profiles(i,j)%cloud_water_q(nlyr), stat=alloc_status)
          allocate(profiles(i,j)%cloud_ice_q(nlyr), stat=alloc_status)
          allocate(profiles(i,j)%rain_q(nlyr), stat=alloc_status)
          allocate(profiles(i,j)%snow_q(nlyr), stat=alloc_status)
          allocate(profiles(i,j)%graupel_q(nlyr), stat=alloc_status)
          if (n_moments .eq. 2) then
             allocate(profiles(i,j)%hail_q(nlyr), stat=alloc_status)
             allocate(profiles(i,j)%cloud_water_n(nlyr), stat=alloc_status)
             allocate(profiles(i,j)%cloud_ice_n(nlyr), stat=alloc_status)
             allocate(profiles(i,j)%rain_n(nlyr), stat=alloc_status)
             allocate(profiles(i,j)%snow_n(nlyr), stat=alloc_status)
             allocate(profiles(i,j)%graupel_n(nlyr), stat=alloc_status)
             allocate(profiles(i,j)%hail_n(nlyr), stat=alloc_status)
          end if
          allocate(profiles(i,j)%vapor_pressure(nlyr), stat=alloc_status)
          allocate(profiles(i,j)%rho_vap(nlyr), stat=alloc_status)
          allocate(profiles(i,j)%q_hum(nlyr), stat=alloc_status)
       end do
    end do

    ! $## think about order of reading
    do i = 1, ngridx
       do j = 1, ngridy 
          read(14,*,iostat=istat) profiles(i,j)%isamp, profiles(i,j)%jsamp !
          if (istat .ne. 0) then
             call error_msg(input_file,i,j)
             stop "Read error 3: Cannot read profile index i,j"
          end if

          read(14,*,iostat=istat) &
               profiles(i,j)%latitude, &         ! degree
               profiles(i,j)%longitude,&         ! degree
               profiles(i,j)%land_fraction,&     !
               profiles(i,j)%wind_10u,&          ! m/s
               profiles(i,j)%wind_10v            ! m/s
          if (istat .ne. 0) then
             call error_msg(input_file,i,j)
             stop "Read error 4: Cannot read profile lat/lon/lfrac/wind"
          end if
          ! integrated quantities
          if (n_moments .eq. 1) then
             read(14,*,iostat=istat) &
                  profiles(i,j)%iwv,&               ! kg/m^2
                  profiles(i,j)%cwp,&               ! kg/m^2
                  profiles(i,j)%iwp,&               ! kg/m^2
                  profiles(i,j)%rwp,&               ! kg/m^2
                  profiles(i,j)%swp,&               ! kg/m^2
                  profiles(i,j)%gwp                 ! kg/m^2
             profiles(i,j)%hwp = 0.
          end if
          if (n_moments .eq. 2) then
             read(14,*,iostat=istat) &
                  profiles(i,j)%iwv,&               ! kg/m^2
                  profiles(i,j)%cwp,&               ! kg/m^2
                  profiles(i,j)%iwp,&               ! kg/m^2
                  profiles(i,j)%rwp,&               ! kg/m^2
                  profiles(i,j)%swp,&               ! kg/m^2
                  profiles(i,j)%gwp,&               ! kg/m^2
                  profiles(i,j)%hwp                 ! kg/m^2
          end if

          if (istat .ne. 0) then
             call error_msg(input_file,i,j)
             stop "Read error 5: Cannot read profile integrated quantities"
          end if

          ! surface values
          read(14,*,iostat=istat) &
               profiles(i,j)%hgt_lev(0),&
               profiles(i,j)%press_lev(0),&
               profiles(i,j)%temp_lev(0),&
               profiles(i,j)%relhum_lev(0)
          if (istat .ne. 0) call error_msg(input_file, i, j)
          do k = 1, nlyr 
             if (n_moments .eq. 1) then
                read(14,*,iostat=istat) &
                     
                     profiles(i,j)%hgt_lev(k), &             ! m
                     profiles(i,j)%press_lev(k), &           ! Pa
                     profiles(i,j)%temp_lev(k), &            ! K
                     profiles(i,j)%relhum_lev(k), &          ! %
                     profiles(i,j)%cloud_water_q(k), &       ! kg/kg
                     profiles(i,j)%cloud_ice_q(k), &         ! kg/kg
                     profiles(i,j)%rain_q(k), &              ! kg/kg
                     profiles(i,j)%snow_q(k), &              ! kg/kg
                     profiles(i,j)%graupel_q(k)              ! kg/kg
             end if
             if (n_moments .eq. 2) then
                read(14,*,iostat=istat) &
                     profiles(i,j)%hgt_lev(k), &             ! m
                     profiles(i,j)%press_lev(k), &           ! Pa
                     profiles(i,j)%temp_lev(k), &            ! K
                     profiles(i,j)%relhum_lev(k), &          ! %
                     profiles(i,j)%cloud_water_q(k), &       ! kg/kg
                     profiles(i,j)%cloud_ice_q(k), &         ! kg/kg
                     profiles(i,j)%rain_q(k), &              ! kg/kg
                     profiles(i,j)%snow_q(k), &              ! kg/kg
                     profiles(i,j)%graupel_q(k), &           ! kg/kg
                     profiles(i,j)%hail_q(k), &              ! kg/kg
                     profiles(i,j)%cloud_water_n(k), &       ! #/kg
                     profiles(i,j)%cloud_ice_n(k), &         ! #/kg
                     profiles(i,j)%rain_n(k), &              ! #/kg
                     profiles(i,j)%snow_n(k), &              ! #/kg
                     profiles(i,j)%graupel_n(k), &           ! #/kg
                     profiles(i,j)%hail_n(k)                 ! #/kg
             end if
             if (istat .ne. 0) then
                call error_msg(input_file,i,j,k)
                print *,"Error at layer", k
                stop "Read error 6: Cannot read profile values"
             end if
          end do
       end do
    end do

!	do i = 1, ngridx
!       do j = 1, ngridy
!                  profiles(i,j)%cwp  =0.             ! kg/m^2
!                  profiles(i,j)%iwp  =0.               ! kg/m^2
!                  profiles(i,j)%rwp  =0.              ! kg/m^2
!                  profiles(i,j)%swp  =0.              ! kg/m^2
!                  profiles(i,j)%gwp  =0.
!     if (n_moments .eq. 2)              profiles(i,j)%hwp  =0.
!       do k = 1, nlyr
!	profiles(i,j)%cloud_water_q(k) = 0.
!	profiles(i,j)%cloud_ice_q(k) = 0.
!	profiles(i,j)%rain_q(k) = 0.
!	profiles(i,j)%snow_q(k) = 0.
!	profiles(i,j)%graupel_q(k) = 0.
!if (n_moments .eq. 2) 	profiles(i,j)%hail_q(k) = 0.
!	end do
!	end do
!	end do

    close(14)



    if (verbose .gt. 0) print*, 'profile reading done!'
    ! 
    return
  end subroutine vars_atmosphere_read_profile
end module vars_atmosphere
