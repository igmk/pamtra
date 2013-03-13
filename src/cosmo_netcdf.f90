module cosmo_netcdf
    use netcdf
    implicit none
    save

    character(4) :: yyyy, hhmm
    character(2) :: mm, dd

    ! dimensions(rlon,rlat)
    !  real, dimension(:,:),allocatable :: lats, lons

    ! The variables for the 4D fields but we have only 1 time step
    ! dimensions(ntime,nlvl,rlat,nlon)
    ! on full-levels nlvl (50)
    real, dimension(:,:,:),allocatable :: pres      ! [Pa]
    real, dimension(:,:,:),allocatable :: temp      ! [K]
    real, dimension(:,:,:),allocatable :: u, v    ! [m/s]
    ! mixing ratios of
    real, dimension(:,:,:),allocatable :: vapor       ! [kg/kg]
    real, dimension(:,:,:),allocatable :: cloud       ! [kg/kg]
    real, dimension(:,:,:),allocatable :: ice         ! [kg/kg]
    real, dimension(:,:,:),allocatable :: rain        ! [kg/kg]
    real, dimension(:,:,:),allocatable :: snow        ! [kg/kg]
    real, dimension(:,:,:),allocatable :: graupel     ! [kg/kg]
    real, dimension(:,:,:),allocatable :: hail        ! [kg/kg]

    real, dimension(:,:,:),allocatable :: cloud_n     ! [#/kg]
    real, dimension(:,:,:),allocatable :: ice_n       ! [#/kg]
    real, dimension(:,:,:),allocatable :: rain_n      ! [#/kg]
    real, dimension(:,:,:),allocatable :: snow_n      ! [#/kg]
    real, dimension(:,:,:),allocatable :: graupel_n   ! [#/kg]
    real, dimension(:,:,:),allocatable :: hail_n      ! [#/kg]



    real, dimension(:,:,:),allocatable :: h           ! [m]
    ! on half-levels nlvl1 (51)
    real, dimension(:,:,:),allocatable :: preshl,temphl,vaporhl,hhl
    real, dimension(:,:,:),allocatable :: cloudhl,icehl,snowhl,graupelhl,rainhl,hailhl
    real, dimension(:,:,:),allocatable :: cloudhl_n,icehl_n,snowhl_n,graupelhl_n,rainhl_n,hailhl_n

    ! The variables for the 3D fields with 1 time step
    ! dimensions (ntime,rlat,rlon)
    real, dimension(:,:),allocatable :: tsurf       ! K
    real, dimension(:,:),allocatable :: psurf       ! Pa
    real, dimension(:,:),allocatable :: vaporsurf

    real, dimension(:,:),allocatable :: u10, v10    ! [m/s]

    real, dimension(:,:),allocatable :: lon, lat    ! °
    real, dimension(:,:),allocatable :: lfrac

contains
    subroutine read_cosmo_netcdf_acp(coords, nlvl)

        use settings, only: crm_data, crm_constants, n_moments
        use report_module

        implicit none

        ! This will be the netCDF ID for the file and data variable.
        integer :: nx,ny,nlvl,nlvl1
        integer :: ncid,ncidc

        integer, dimension(2) :: start_2d, count_2d
        integer, dimension(3) :: start_3d, count_fl, count_hl
        integer, dimension(4) :: coords

        if (verbose .gt. 0) print*, 'start reading cosmo_data from acp'

        call check(nf90_open(trim(crm_data), nf90_nowrite, ncid))

        call get_dimensions(ncid,"lv_HYBY2",nlvl)
        nlvl1 = nlvl+1

        nx = coords(2)-coords(1)+1
        ny = coords(4)-coords(3)+1
        start_2d = (/coords(1),coords(3)/)
        start_3d = (/coords(1),coords(3),1/)
        count_2d = (/nx,ny/)
        count_fl = (/nx,ny,nlvl/)
        count_hl = (/nx,ny,nlvl1/)


        !  offset = 1105228800 ! this is the time offset with respect to
        !                      ! 01-01-1970 00UTC for the seconds stored in the time variable

        ! we start with the time


        !  call check(nf90_inq_varid(ncid, 'time', timeID))

        !  call check(nf90_get_var(ncid, timeID, time))
  
        !  call gmtime(time+offset,gmtvalues)

        !  write(yyyy,'(i4)') gmtvalues(6)+1900
        !  write(mm,'(i2.2)') gmtvalues(5)+1
        !  write(dd,'(i2.2)') gmtvalues(4)+2
        !  write(hhmm,'(i2.2,i2.2)') gmtvalues(3),gmtvalues(2)

        !  write(yyyy,'(a4)') nc_file(len-12:len-9)
        !  write(mm,'(a2)') nc_file(len-8:len-7)
        !  write(dd,'(a2)') nc_file(len-6:len-5)
        !  write(hhmm,'(a4)') nc_file(len-4:len-3)//'00'

        yyyy = '2010'
        mm = '07'
        dd = '16'
        hhmm = '1800'

        call get_2d_vars(ncid, 'T_GDS10_SFC', tsurf, nx,ny, start_2d, count_2d)
        call get_2d_vars(ncid, 'PS_GDS10_SFC', psurf, nx,ny, start_2d, count_2d)
        call get_3d_vars(ncid, 'T_GDS10_HYBY', temp, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'PS_GDS10_HYBY', pres, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'U_GDS10_HYBY', u, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'V_GDS10_HYBY', v, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QV_GDS10_HYBY', vapor, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QC_GDS10_HYBY', cloud, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QI_GDS10_HYBY', ice, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QR_GDS10_HYBY', rain, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QG_GDS10_HYBY', graupel, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QS_GDS10_HYBY', snow, nx,ny, nlvl, start_3d, count_fl)
        if (n_moments == 2) then
            call get_3d_vars(ncid, 'QH_GDS10_HYBY', hail, nx,ny, nlvl, start_3d, count_fl)
              ! nc file with total number concentration
            !   call check(nf90_open(trim(nc_filen), nf90_nowrite, ncidn))
            call get_3d_vars(ncid, 'QNC_GDS10_HYBY', cloud_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncid, 'QNI_GDS10_HYBY', ice_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncid, 'QNR_GDS10_HYBY', rain_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncid, 'QNG_GDS10_HYBY', graupel_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncid, 'QNS_GDS10_HYBY', snow_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncid, 'QNH_GDS10_HYBY', hail_n, nx,ny, nlvl, start_3d, count_fl)

        !   call check(nf90_close(ncidn))
        else
            allocate(hail(nx,ny,nlvl))
            hail = 0.
        endif

        call check(nf90_close(ncid))

        call check(nf90_open(trim(crm_constants), nf90_nowrite, ncidc))

        call get_2d_vars(ncid, 'lon', lon, nx,ny, start_2d, count_2d)
        call get_2d_vars(ncid, 'lat', lat, nx,ny, start_2d, count_2d)
        call get_3d_vars(ncidc, 'HHL', hhl, nx,ny, nlvl1, start_3d, count_hl)
        call get_2d_vars(ncidc, 'FR_LAND', lfrac, nx,ny, start_2d, count_2d)

        call check(nf90_close(ncidc))

        allocate(h(nx,ny,nlvl))

        allocate(cloudhl(nx,ny,nlvl1),vaporhl(nx,ny,nlvl1))
        allocate(icehl(nx,ny,nlvl1),snowhl(nx,ny,nlvl1))
        allocate(graupelhl(nx,ny,nlvl1),rainhl(nx,ny,nlvl1))
        if (n_moments == 2) then
            allocate(hailhl(nx,ny,nlvl1))
            allocate(icehl_n(nx,ny,nlvl1),snowhl_n(nx,ny,nlvl1))
            allocate(graupelhl_n(nx,ny,nlvl1),rainhl_n(nx,ny,nlvl1))
            allocate(cloudhl_n(nx,ny,nlvl1),hailhl_n(nx,ny,nlvl1))
        endif

        allocate(preshl(nx,ny,nlvl1),temphl(nx,ny,nlvl1))

        call half_levels(nx,ny,nlvl)

        ! If we got this far, everything worked as expected. Yipee!
        !  print *,"*** SUCCESS reading example file ", trim(nc_file)//".nc", "! "

        if (verbose .gt. 0) print*, 'finished reading cosmo_data from acp'

    end subroutine read_cosmo_netcdf_acp

    subroutine read_cosmo_netcdf_narval(coords, nlvl)

        use settings, only: crm_data, crm_data2, crm_constants, n_moments
        use report_module

        implicit none

        ! This will be the netCDF ID for the file and data variable.
        integer :: nx,ny,nlvl,nlvl1
        integer :: ncid,ncidc, ncidn

        integer, dimension(2) :: start_2d, count_2d
        integer, dimension(3) :: start_3d, count_fl, count_hl
        integer, dimension(4) :: coords

        real, allocatable, dimension(:,:) :: u_dum, v_dum

        if (verbose .gt. 0) print*, 'start reading cosmo_data from narval'

        call check(nf90_open(trim(crm_data), nf90_nowrite, ncid))

        call get_dimensions(ncid,"level",nlvl)
        nlvl1 = nlvl+1

        nx = coords(2)-coords(1)+1
        ny = coords(4)-coords(3)+1
        start_2d = (/coords(1),coords(3)/)
        start_3d = (/coords(1),coords(3),1/)
        count_2d = (/nx,ny/)
        count_fl = (/nx,ny,nlvl/)
        count_hl = (/nx,ny,nlvl1/)

        print*, coords

        !  offset = 1105228800 ! this is the time offset with respect to
        !                      ! 01-01-1970 00UTC for the seconds stored in the time variable

        ! we start with the time


        !  call check(nf90_inq_varid(ncid, 'time', timeID))

        !  call check(nf90_get_var(ncid, timeID, time))

        !  call gmtime(time+offset,gmtvalues)

        !  write(yyyy,'(i4)') gmtvalues(6)+1900
        !  write(mm,'(i2.2)') gmtvalues(5)+1
        !  write(dd,'(i2.2)') gmtvalues(4)+2
        !  write(hhmm,'(i2.2,i2.2)') gmtvalues(3),gmtvalues(2)

        !  write(yyyy,'(a4)') nc_file(len-12:len-9)
        !  write(mm,'(a2)') nc_file(len-8:len-7)
        !  write(dd,'(a2)') nc_file(len-6:len-5)
        !  write(hhmm,'(a4)') nc_file(len-4:len-3)//'00'

        yyyy = '2008'
        mm = '01'
        dd = '17'
        hhmm = '1800'

        if (verbose .gt. 1) print*, "Reading T_G"
        call get_2d_vars(ncid, 'T_G', tsurf, nx,ny, start_2d, count_2d)
        if (verbose .gt. 1) print*, "Reading PS"
        call get_2d_vars(ncid, 'PS', psurf, nx,ny, start_2d, count_2d)
        if (verbose .gt. 1) print*, "Reading T"
        call get_3d_vars(ncid, 'T', temp, nx,ny, nlvl, start_3d, count_fl)
        if (verbose .gt. 1) print*, "Reading P"
        call get_3d_vars(ncid, 'P', pres, nx,ny, nlvl, start_3d, count_fl)
        if (verbose .gt. 1) print*, "Reading U_10M"
        call get_2d_vars(ncid, 'U_10M', u_dum, nx,ny, start_2d, count_2d)
        if (verbose .gt. 1) print*, "Reading V_10M"
        call get_2d_vars(ncid, 'V_10M', v_dum, nx,ny, start_2d, count_2d)
        allocate(u(nx,ny,nlvl),v(nx,ny,nlvl))
        u(:,:,nlvl) = u_dum
        v(:,:,nlvl) = v_dum
        deallocate(u_dum,v_dum)
        if (verbose .gt. 1) print*, "Reading QV"
        call get_3d_vars(ncid, 'QV', vapor, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QC', cloud, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QI', ice, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QR', rain, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QG', graupel, nx,ny, nlvl, start_3d, count_fl)
        call get_3d_vars(ncid, 'QS', snow, nx,ny, nlvl, start_3d, count_fl)
        if (n_moments == 2) then
            call get_3d_vars(ncid, 'QH', hail, nx,ny, nlvl, start_3d, count_fl)
              !	nc file with total number concentration
            call check(nf90_open(trim(crm_data2), nf90_nowrite, ncidn))
            call get_3d_vars(ncidn, 'QNC', cloud_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncidn, 'QNI', ice_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncidn, 'QNR', rain_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncidn, 'QNG', graupel_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncidn, 'QNS', snow_n, nx,ny, nlvl, start_3d, count_fl)
            call get_3d_vars(ncidn, 'QNH', hail_n, nx,ny, nlvl, start_3d, count_fl)

         	call check(nf90_close(ncidn))
        else
            allocate(hail(nx,ny,nlvl))
            hail = 0.
        endif

        call check(nf90_close(ncid))

        call check(nf90_open(trim(crm_constants), nf90_nowrite, ncidc))

        call get_2d_vars(ncid, 'lon', lon, nx,ny, start_2d, count_2d)
        call get_2d_vars(ncid, 'lat', lat, nx,ny, start_2d, count_2d)
        call get_3d_vars(ncidc, 'HHL', hhl, nx,ny, nlvl1, start_3d, count_hl)
        call get_2d_vars(ncidc, 'FR_LAND', lfrac, nx,ny, start_2d, count_2d)

        call check(nf90_close(ncidc))

        allocate(h(nx,ny,nlvl))

        allocate(cloudhl(nx,ny,nlvl1),vaporhl(nx,ny,nlvl1))
        allocate(icehl(nx,ny,nlvl1),snowhl(nx,ny,nlvl1))
        allocate(graupelhl(nx,ny,nlvl1),rainhl(nx,ny,nlvl1))
        if (n_moments == 2) then
            allocate(hailhl(nx,ny,nlvl1))
            allocate(icehl_n(nx,ny,nlvl1),snowhl_n(nx,ny,nlvl1))
            allocate(graupelhl_n(nx,ny,nlvl1),rainhl_n(nx,ny,nlvl1))
            allocate(cloudhl_n(nx,ny,nlvl1),hailhl_n(nx,ny,nlvl1))
        endif

        allocate(preshl(nx,ny,nlvl1),temphl(nx,ny,nlvl1))

        call half_levels(nx,ny,nlvl)

        ! If we got this far, everything worked as expected. Yipee!
        !  print *,"*** SUCCESS reading example file ", trim(nc_file)//".nc", "! "

        if (verbose .gt. 0) print*, 'finished reading cosmo_data from acp'

    end subroutine read_cosmo_netcdf_narval

    subroutine get_dimensions(ncid, lvl_dim_name,nl)
        implicit none

        integer :: nl
        integer, intent(in) :: ncid
        integer :: lvlDimID

        character(len=*) :: lvl_dim_name

        call check(nf90_inq_dimid(ncid,lvl_dim_name,lvlDimID))

        call check(nf90_inquire_dimension(ncid,lvlDimID,len=nl))

        return

    end subroutine get_dimensions

    subroutine get_3d_vars(ncid, var_name, var, &
    dim1, dim2, dim3, start_point, counter)

        implicit none

        integer :: VarID

        integer, intent(in) :: ncid, dim1, dim2, dim3

        integer, dimension(3), intent(in) :: start_point, counter

        real, dimension(:,:,:), allocatable, intent(inout) :: var

        character (len = *), intent(in) :: var_name

        allocate(var(dim1, dim2, dim3))

        call check(nf90_inq_varid(ncid, var_name, VarID))

        call check(nf90_get_var(ncid, VarID, var(1:dim1,1:dim2,:), &
        start = start_point, count=counter))

        return

    end subroutine get_3d_vars

    subroutine get_2d_vars(ncid, var_name, var, dim1, dim2, start_point, counter)

        implicit none

        integer :: VarID

        integer, intent(in) :: ncid, dim1, dim2

        integer, dimension(2), intent(in) :: start_point, counter

        real, dimension(:,:), allocatable, intent(inout) :: var

        character (len = *), intent(in) :: var_name

        allocate(var(dim1, dim2))

        call check(nf90_inq_varid(ncid, var_name, VarID))

        call check(nf90_get_var(ncid, VarID, var(1:dim1,1:dim2), &
        start = start_point, count = counter))

        return

    end subroutine get_2d_vars
    subroutine check(status)

        implicit none

        integer, intent(in) :: status

        if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
        end if
    end subroutine check

    subroutine half_levels(nx,ny,nlvl)

        implicit none

        integer :: i,j,k
        integer :: nx,ny,nlvl, nlvl1

        real :: meant,dz_for_pb,dz,xp
        real :: r_d,r_v

        r_d=287.05        ! [J/(kg*K)] = [m^2/(s^2*K)]
        r_v=461.5
        nlvl1 = nlvl+1

        do i = 1, nx ! lon
            do j = 1, ny ! lat

                ! First we calculate the height of the full levels

                do k = 1, nlvl
                    h(i,j,k) = (hhl(i,j,k+1) + hhl(i,j,k))/2
                end do

                do k = 2, nlvl
                    temphl(i,j,k) = (temp(i,j,k)+temp(i,j,k-1))/2
                    vaporhl(i,j,k) = (vapor(i,j,k)+vapor(i,j,k-1))/2
                    dz = h(i,j,k)-h(i,j,k-1)
                    if (pres(i,j,k-1) == pres(i,j,k)) then
                        preshl(i,j,k) = pres(i,j,k)
                    else
                        xp = -alog(pres(i,j,k)/pres(i,j,k-1))/dz
                        preshl(i,j,k) = -pres(i,j,k-1)/xp*(exp(-xp*dz)-1.)/dz
                    end if
                end do

                temphl(i,j,1) = temp(i,j,1) + 0.25*(temp(i,j,1)-temp(i,j,2))
                vaporhl(i,j,1) = vapor(i,j,1) + 0.25*(vapor(i,j,1)-vapor(i,j,2))
                temphl(i,j,nlvl1) = temp(i,j,nlvl) + 0.25*(temp(i,j,nlvl)-temp(i,j,nlvl-1))
                vaporhl(i,j,nlvl1) = vapor(i,j,nlvl) + 0.25*(vapor(i,j,nlvl)-vapor(i,j,nlvl-1))

                meant = (temphl(i,j,1) + temp(i,j,1))/2
                dz_for_pb = (hhl(i,j,1)-h(i,j,1))
                preshl(i,j,1) = pres(i,j,1) * exp(-9.8062/(meant*287.05)*dz_for_pb)

                meant = (temphl(i,j,nlvl1) + temp(i,j,nlvl))/2
                dz_for_pb = (hhl(i,j,nlvl1)-h(i,j,nlvl))

                preshl(i,j,nlvl1) = pres(i,j,nlvl) * exp(-9.8062/(meant*287.05)*dz_for_pb)

            end do
        end do
        return

    end subroutine half_levels

    subroutine integrate_profiles(hhl,p,t,qc,qr,qi,qs,qg,qh,qv,cwp,iwp,rwp,swp,gwp,hwp,iwv,nlvl)

        use settings, only: n_moments

        implicit none

        integer :: i, nlvl

        real :: z ! [m] layer thickness

        real, parameter ::r_d = 287.05 , &  ! [J/(kg K)]  individual gas constant for dry air
        r_v = 461.5       ! [J/(kg K)]  individual gas constant for water vapor

        real, dimension(nlvl) :: cwcl, iwcl, rwcl, swcl, gwcl, iwvl, hwcl

        real, dimension(nlvl), intent(in) :: p, t, qc, qi, qr, qs ,qg, qv, qh
        real, dimension(nlvl+1), intent(in) :: hhl

        real, dimension(nlvl) :: wet_rho  ! [kg/m**3] density of humid air (incl hydrometeors)
        real, intent(out) :: cwp, iwp, rwp, swp, gwp, iwv, hwp

        do i = nlvl,1,-1
            if (n_moments == 1) wet_rho(i) = p(i)/(r_d*(1+(r_v/r_d-1)*qv(i)-qs(i)-qc(i)-qi(i)-qr(i)-qg(i))*t(i))
            if (n_moments == 2) wet_rho(i) = p(i)/(r_d*(1+(r_v/r_d-1)*qv(i)-qs(i)-qc(i)-qi(i)-qr(i)-qg(i)-qh(i))*t(i))
            z = (hhl(i)-hhl(i+1))
            cwcl(i) = wet_rho(i) * z * qc(i)
            iwcl(i) = wet_rho(i) * z * qi(i)
            rwcl(i) = wet_rho(i) * z * qr(i)
            swcl(i) = wet_rho(i) * z * qs(i)
            gwcl(i) = wet_rho(i) * z * qg(i)
            if (n_moments == 2) hwcl(i) = wet_rho(i) * z * qh(i)
            iwvl(i) = wet_rho(i) * z * qv(i)
        end do

        cwp = sum(cwcl)
        iwp = sum(iwcl)
        rwp = sum(rwcl)
        swp = sum(swcl)
        gwp = sum(gwcl)
        if (n_moments == 2) hwp = sum(hwcl)
        iwv = sum(iwvl)

        return

    end subroutine integrate_profiles
end module cosmo_netcdf
