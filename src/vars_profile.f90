module vars_profile

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

    type(profile), allocatable :: profiles(:,:)


    integer :: profiles_nlyr
    integer :: profiles_ngridx, profiles_ngridy
    integer, dimension(4) :: coords
    character(2) :: profiles_month, profiles_day
    character(4) :: profiles_year, profiles_time
    real(kind=sgl) :: profiles_deltax, profiles_deltay

contains

    subroutine allocate_profiles

        use settings, only: n_moments
        use report_module

        implicit none

        integer :: i,j,alloc_status

        allocate(profiles(profiles_ngridx,profiles_ngridy),stat=alloc_status)
        do i = 1, profiles_ngridx
            do j = 1, profiles_ngridy
                allocate(profiles(i,j)%hgt_lev(0:profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%press_lev(0:profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%press(profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%temp_lev(0:profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%temp(profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%relhum_lev(0:profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%relhum(profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%cloud_water_q(profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%cloud_ice_q(profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%rain_q(profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%snow_q(profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%graupel_q(profiles_nlyr), stat=alloc_status)
                if (n_moments .eq. 2) then
                    allocate(profiles(i,j)%hail_q(profiles_nlyr), stat=alloc_status)
                    allocate(profiles(i,j)%cloud_water_n(profiles_nlyr), stat=alloc_status)
                    allocate(profiles(i,j)%cloud_ice_n(profiles_nlyr), stat=alloc_status)
                    allocate(profiles(i,j)%rain_n(profiles_nlyr), stat=alloc_status)
                    allocate(profiles(i,j)%snow_n(profiles_nlyr), stat=alloc_status)
                    allocate(profiles(i,j)%graupel_n(profiles_nlyr), stat=alloc_status)
                    allocate(profiles(i,j)%hail_n(profiles_nlyr), stat=alloc_status)
                end if
                allocate(profiles(i,j)%vapor_pressure(profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%rho_vap(profiles_nlyr), stat=alloc_status)
                allocate(profiles(i,j)%q_hum(profiles_nlyr), stat=alloc_status)
            end do
        end do

    end subroutine allocate_profiles

    subroutine vars_profile_read_profile(errorstatus)

        use settings, only: verbose, input_path, input_file, &
        output_path, nc_out_file, file_desc, n_moments, freq_str

        use report_module

        integer :: i,j,k

        ! Error handling

        integer(kind=long),intent(out) :: errorstatus
        integer(kind=long) :: err = 0
        character(len=200) :: msg
        character(len=25) :: nameOfRoutine = 'vars_profile_read_profile'


        if (verbose >= 1) print *,"opening: ",input_file

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
        STATUS='OLD', form='formatted',iostat=err)

        if (err /= 0) then
            msg = "Read error 1: Cannot open file "//input_path(:len_trim(input_path))//"/"//input_file(:len_trim(input_file))
            call report(err,msg,nameOfRoutine)
            errorstatus = err
            return
        end if

        read(14,*,iostat=err) profiles_year, profiles_month, profiles_day, profiles_time, &
        profiles_ngridx, profiles_ngridy, profiles_nlyr, profiles_deltax, profiles_deltay

        if (err /= 0) then
            msg = "Read error 2: Cannot read first line of file "//trim(input_file)
            call report(err,msg,nameOfRoutine)
            errorstatus = err
            return
        end if

        call allocate_profiles!(profiles_ngridx,profiles_ngridy,profiles_nlyr)

        ! $## think about order of reading
        do i = 1, profiles_ngridx
            do j = 1, profiles_ngridy
                read(14,*,iostat=err) profiles(i,j)%isamp, profiles(i,j)%jsamp !
                if (err /= 0) then
                    msg = "Read error 3: Cannot read profile index i,j in"//trim(input_file)
                    call report(err,msg,nameOfRoutine)
                    errorstatus = err
                    return
                end if
                read(14,*,iostat=err) &
                profiles(i,j)%latitude, &         ! degree
                profiles(i,j)%longitude,&         ! degree
                profiles(i,j)%land_fraction,&     !
                profiles(i,j)%wind_10u,&          ! m/s
                profiles(i,j)%wind_10v            ! m/s
                if (err /= 0) then
                    msg = "Read error 4: Cannot read profile lat/lon/lfrac/wind in"//trim(input_file)
                    call report(err,msg,nameOfRoutine)
                    errorstatus = err
                    return
                end if

                ! integrated quantities
                if (n_moments .eq. 1) then
                    read(14,*,iostat=err) &
                    profiles(i,j)%iwv,&               ! kg/m^2
                    profiles(i,j)%cwp,&               ! kg/m^2
                    profiles(i,j)%iwp,&               ! kg/m^2
                    profiles(i,j)%rwp,&               ! kg/m^2
                    profiles(i,j)%swp,&               ! kg/m^2
                    profiles(i,j)%gwp                 ! kg/m^2
                    profiles(i,j)%hwp = 0.
                end if
                if (n_moments .eq. 2) then
                    read(14,*,iostat=err) &
                    profiles(i,j)%iwv,&               ! kg/m^2
                    profiles(i,j)%cwp,&               ! kg/m^2
                    profiles(i,j)%iwp,&               ! kg/m^2
                    profiles(i,j)%rwp,&               ! kg/m^2
                    profiles(i,j)%swp,&               ! kg/m^2
                    profiles(i,j)%gwp,&               ! kg/m^2
                    profiles(i,j)%hwp                 ! kg/m^2
                end if
                if (err /= 0) then
                    msg = "Read error 5: Cannot read profile integrated quantities in"//trim(input_file)
                    call report(err,msg,nameOfRoutine)
                    errorstatus = err
                    return
                end if
                ! surface values
                read(14,*,iostat=err) &
                profiles(i,j)%hgt_lev(0),&
                profiles(i,j)%press_lev(0),&
                profiles(i,j)%temp_lev(0),&
                profiles(i,j)%relhum_lev(0)
                if (err /= 0) then
                    write(msg,'(a,i3,x,i3)') "Error in reading profile ",i,j
                    call report(err,msg,nameOfRoutine)
                    errorstatus = err
                    return
                end if
                do k = 1, profiles_nlyr
                    if (n_moments .eq. 1) then
                        read(14,*,iostat=err) &
                     
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
                        read(14,*,iostat=err) &
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
                    if (err /= 0) then
                        write(msg, '(a,i3)') "Read error 6: Cannot read profile values in layer ", k
                        call report(err,msg,nameOfRoutine)
                        errorstatus = err
                        return
                    end if
                end do
            end do
        end do
        close(14)

        nc_out_file = trim(output_path)//"/"//trim(input_file(1:len_trim(input_file)-4))//&
        trim(freq_str)//trim(file_desc)//'.nc'

        if (verbose >= 1) print*, 'profile reading done!'

        return

    end subroutine vars_profile_read_profile


end module vars_profile
