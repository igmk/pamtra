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

    subroutine vars_profile_read_cosmo

        use settings, only: crm_case, n_moments, freq_str, output_path, file_desc, nc_out_file
        use conversions
        use cosmo_netcdf
        use double_moments_module
        use report_module

        implicit none

        integer :: i,j,k, k2
        character(12) :: date_str
        character(16) :: grid_str

        if (verbose .gt. 0) print *,"reading from cosmo output"

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

        if (crm_case .eq. 'acp') then
            call read_cosmo_netcdf_acp(coords,profiles_nlyr)
        elseif (crm_case .eq. 'narval') then
            call read_cosmo_netcdf_narval(coords,profiles_nlyr)
        else
            stop 'no case selected'
        end if

        profiles_ngridx = coords(2)-coords(1)+1
        profiles_ngridy = coords(4)-coords(3)+1

        call allocate_profiles

        profiles_year = yyyy
        profiles_month = mm
        profiles_day = dd
        profiles_time = hhmm

        do i = 1, profiles_ngridx
            do j = 1, profiles_ngridy
                profiles(i,j)%isamp = coords(1) + i - 1
                profiles(i,j)%jsamp = coords(3) + j - 1
                profiles(i,j)%longitude = lon(i,j)        ! degree
                profiles(i,j)%latitude = lat(i,j)         ! degree
                profiles(i,j)%land_fraction = lfrac(i,j)  !
                profiles(i,j)%wind_10u = u(i,j,1)          ! m/s
                profiles(i,j)%wind_10v = v(i,j,1)          ! m/s
                ! integrated quantities
                call integrate_profiles(hhl(i,j,:),pres(i,j,:),temp(i,j,:),cloud(i,j,:),rain(i,j,:),&
                ice(i,j,:),snow(i,j,:),graupel(i,j,:),hail(i,j,:),vapor(i,j,:),&
                profiles(i,j)%cwp,profiles(i,j)%iwp,profiles(i,j)%rwp,profiles(i,j)%swp,&
                profiles(i,j)%gwp,profiles(i,j)%hwp,profiles(i,j)%iwv,profiles_nlyr)
                ! surface values
                profiles(i,j)%hgt_lev(0)=hhl(i,j,profiles_nlyr+1)
                profiles(i,j)%press_lev(0)=psurf(i,j)
                profiles(i,j)%temp_lev(0)=tsurf(i,j)
                profiles(i,j)%relhum_lev(0)=vapor2rh(temphl(i,j,profiles_nlyr+1),psurf(i,j),vaporhl(i,j,profiles_nlyr+1))
                do k = 1, profiles_nlyr
                    k2 = profiles_nlyr+1-k
                    if (n_moments .eq. 1) then
                        profiles(i,j)%hgt_lev(k) = hhl(i,j,k2)             ! m
                        profiles(i,j)%press_lev(k) = preshl(i,j,k2)           ! Pa
                        profiles(i,j)%temp_lev(k) = temphl(i,j,k2)            ! K
                        profiles(i,j)%relhum_lev(k) = vapor2rh(temphl(i,j,k2),preshl(i,j,k2),vaporhl(i,j,k2))          ! %
                        profiles(i,j)%cloud_water_q(k) = cloud(i,j,k2)       ! kg/kg
                        profiles(i,j)%cloud_ice_q(k) = ice(i,j,k2)         ! kg/kg
                        profiles(i,j)%rain_q(k) = rain(i,j,k2)              ! kg/kg
                        profiles(i,j)%snow_q(k) = snow(i,j,k2)              ! kg/kg
                        profiles(i,j)%graupel_q(k) = graupel(i,j,k2)         ! kg/kg
                    end if
                    if (n_moments .eq. 2) then
                        profiles(i,j)%hgt_lev(k) = hhl(i,j,k2)             ! m
                        profiles(i,j)%press_lev(k) = preshl(i,j,k2)           ! Pa
                        profiles(i,j)%temp_lev(k) = temphl(i,j,k2)            ! K
                        profiles(i,j)%relhum_lev(k) = vapor2rh(temphl(i,j,k2),preshl(i,j,k2),vaporhl(i,j,k2)) ! %
                        profiles(i,j)%cloud_water_q(k) = cloud(i,j,k2)       ! kg/kg
                        profiles(i,j)%cloud_ice_q(k) = ice(i,j,k2)         ! kg/kg
                        profiles(i,j)%rain_q(k) = rain(i,j,k2)              ! kg/kg
                        profiles(i,j)%snow_q(k) = snow(i,j,k2)              ! kg/kg
                        profiles(i,j)%graupel_q(k) = graupel(i,j,k2)           ! kg/kg
                        profiles(i,j)%hail_q(k) = hail(i,j,k2)              ! kg/kg
                        profiles(i,j)%cloud_water_n(k) = max(cloud_n(i,j,k2),cloud(i,j,k2)/gamma_cloud(5))       ! #/kg
                        profiles(i,j)%cloud_ice_n(k) = max(ice_n(i,j,k2),ice(i,j,k2)/gamma_ice(5))         ! #/kg
                        profiles(i,j)%rain_n(k) = max(rain_n(i,j,k2),rain(i,j,k2)/gamma_rain(5))              ! #/kg
                        profiles(i,j)%snow_n(k) = max(snow_n(i,j,k2),snow(i,j,k2)/gamma_snow(5))              ! #/kg
                        profiles(i,j)%graupel_n(k) = max(graupel_n(i,j,k2),graupel(i,j,k2)/gamma_graupel(5))           ! #/kg
                        profiles(i,j)%hail_n(k) = max(hail_n(i,j,k2),hail(i,j,k2)/gamma_hail(5))              ! #/kg
                    end if
                end do
            end do
        end do

        !call write_profile
        date_str = yyyy//mm//dd//hhmm
        write(grid_str,'(I4.4,I4.4,I4.4,I4.4)') coords
        nc_out_file = trim(output_path)//"/"//'cosmo_'//date_str//'_'//grid_str//&
        trim(freq_str)//trim(file_desc)//'.nc'

        if (verbose .gt. 0) print*, 'profile reading from cosmo done!'

        return
    end subroutine vars_profile_read_cosmo

    subroutine write_profile
        use cosmo_netcdf
        implicit none
        integer :: i,j,k
        write(33,'(A4,1X,A2,1X,A2,1X,A4,1X,I3,1X,I3,1X,i2,1X,F3.1,1X,F3.1)')  &
        yyyy, mm, dd, hhmm, profiles_ngridx, profiles_ngridy, 50,2.8,2.8
        do i = 1, profiles_ngridx
            do j = 1, profiles_ngridy
                write(33,*)&
                profiles(i,j)%isamp,&
                profiles(i,j)%jsamp
                write(33,'(2(F8.3,1X),F5.3,1X,F7.3,1X,F7.3)')&
                profiles(i,j)%latitude,&
                profiles(i,j)%longitude,&
                profiles(i,j)%land_fraction,&
                profiles(i,j)%wind_10u,&
                profiles(i,j)%wind_10v
                ! integrated quantities
                write(33,'(29X,F7.3,1X,6(ES12.6,1X))')&
                profiles(i,j)%iwv,&
                profiles(i,j)%cwp,&
                profiles(i,j)%iwp,&
                profiles(i,j)%rwp,&
                profiles(i,j)%swp,&
                profiles(i,j)%gwp,&
                profiles(i,j)%hwp
                ! surface values
                write(33,'(F9.2,1X,F9.2,2X,F7.3,1X,F7.3)')&
                profiles(i,j)%hgt_lev(0),&
                profiles(i,j)%press_lev(0),&
                profiles(i,j)%temp_lev(0),&
                profiles(i,j)%relhum_lev(0)
                do k = profiles_nlyr,1,-1
                    write(33,'(F9.2,1X,F9.2,2X,F7.3,1X,F7.3,1X,12(ES12.6,1X))') &
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
                    profiles(i,j)%hail_n(k)
                end do
            end do
        end do
    end subroutine write_profile
end module vars_profile
