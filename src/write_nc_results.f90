subroutine write_nc_results

    use kinds
    use vars_output
    use vars_atmosphere, only: atmo_ngridx, atmo_ngridy, atmo_max_nlyrs,&
      atmo_year, atmo_month, atmo_day, atmo_time, atmo_model_i,&
      atmo_model_j, atmo_lfrac, atmo_lon, atmo_lat, atmo_iwv
    use netcdf
    use settings, only: active, passive, creator, radar_mode, &
    n_moments, radar_nfft, radar_mode, nfrq, freqs, nc_out_file
    use report_module

    implicit none

    integer :: ncid
    integer :: dlonID, dlatID, dangID, dfrqID, doutID, dstokesID, dlayerID,dnfftID

    integer :: isVarID, jsVarID, lonVarID, latVarID, lfracVarID, iwvVarID, cwpVarID,&
    iwpVarID, rwpVarID, swpVarID, gwpVarID, hwpVarID, &
    tbVarID, heightVarID, &
    ZeVarID,rQualVarID, &
    AttAtmoVarID, AttHydroVarID,&
    lSloVarID,rSloVarID, lEdgVarID, rEdgVarID, kurtVarID, skewVarID, swVarID, velVarID,&
    frequencyVarID, anglesVarID, RadarVelID, RadarSpecID, RadarSNRID
          
          
    integer :: nang = 32, nout = 2, nstokes = 2

    integer, dimension(2) :: dim2d
    integer, dimension(3) :: dim3d
    integer, dimension(4) :: dim4d
    integer, dimension(5) :: dim5d
    integer, dimension(6) :: dim6d

    integer :: today(3), now(3)

    character(300) :: timestring
    character(40) ::gitVersion,gitHash
    character(10) ::attUnit,zeUnit

    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'write_nc_results'
    
    
    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
    !get git data
    call versionNumber(gitVersion,gitHash)

    if (verbose >= 2) call report(info,'Writing ', nc_out_file)

    call check(nf90_create(path=nc_out_file,cmode=nf90_clobber,ncid=ncid))

    ! for netcdf history get meta data
    call idate(today)   ! today(1)=day, (2)=month, (3)=year
    call itime(now)     ! now(1)=hour, (2)=minute, (3)=second

    write (timestring , "(i2.2, '/', i2.2, '/', i4.4, ' ',  i2.2, ':', i2.2, ':', i2.2)") &
    today(2), today(1), today(3), now
    ! write meta data
    call check(nf90_put_att(ncid,nf90_global, "history", "Created with Pamtra (Version: "//trim(gitVersion)// &
    ", Git Hash: "//trim(gitHash)//")  by "//trim(creator)//" (University of Cologne, IGMK) at "//timestring))
    !TODO: we assume that the times are all the same!
    call check(nf90_put_att(ncid,nf90_global, "data_time",atmo_year(1,1)//"/"//atmo_month(1,1)//"/"// &
          atmo_day(1,1)//"-"//atmo_time(1,1)))

    !make dimensions
    call check(nf90_def_dim(ncid, 'nlon', atmo_ngridx, dlonID))
    call check(nf90_def_dim(ncid, 'nlat', atmo_ngridy, dlatID))
    call check(nf90_def_dim(ncid, 'nfreq', nfrq, dfrqID))
    if (passive) then
        call check(nf90_def_dim(ncid, 'nang', nang, dangID))
        call check(nf90_def_dim(ncid, 'nout', nout, doutID))
        call check(nf90_def_dim(ncid, 'nstokes', nstokes, dstokesID))
    end if
    if (active) then
        call check(nf90_def_dim(ncid, 'atmo_max_nlyrs', atmo_max_nlyrs, dlayerID))
    end if
    if ((active) .and. (radar_mode .eq. "spectrum")) then
        call check(nf90_def_dim(ncid, 'nfft', radar_nfft, dnfftID))
    end if

    !1dim
    if (passive) then
        call check(nf90_def_var(ncid,'angle', nf90_float,(/dangID/), anglesVarID))
        call check(nf90_put_att(ncid, anglesVarID, "units", "deg"))
        call check(nf90_put_att(ncid, anglesVarID, "missing_value", -9999))
    end if

    call check(nf90_def_var(ncid,'frequency', nf90_float,(/dfrqID/), frequencyVarID))
    call check(nf90_put_att(ncid, frequencyVarID, "units", "GHz"))
    call check(nf90_put_att(ncid, frequencyVarID, "missing_value", -9999))

    !create variables and apply meta data
    dim2d = (/dlatID,dlonID/)
    !  call put_2d_var(ncid,'longitude',lons,2,/ngridx,ngridy/)
    call check(nf90_def_var(ncid,'model_i', nf90_int,dim2d, isVarID))
    call check(nf90_put_att(ncid, isVarID, "units", "-"))
    call check(nf90_put_att(ncid, isVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'model_j', nf90_int,dim2d, jsVarID))
    call check(nf90_put_att(ncid, jsVarID, "units", "-"))
    call check(nf90_put_att(ncid, jsVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'longitude', nf90_float,dim2d, lonVarID))
    call check(nf90_put_att(ncid, lonVarID, "units", "deg.dec"))
    call check(nf90_put_att(ncid, lonVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'latitude', nf90_float,dim2d, latVarID))
    call check(nf90_put_att(ncid, latVarID, "units", "deg.dec"))
    call check(nf90_put_att(ncid, latVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'lfrac', nf90_float,dim2d, lfracVarID))
    call check(nf90_put_att(ncid, lfracVarID, "units", "-"))
    call check(nf90_put_att(ncid, lfracVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'iwv', nf90_float,dim2d, iwvVarID))
    call check(nf90_put_att(ncid, iwvVarID, "units", "kg/m^2"))
    call check(nf90_put_att(ncid, iwvVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'cwp', nf90_float,dim2d, cwpVarID))
    call check(nf90_put_att(ncid, cwpVarID, "units", "kg/m^2"))
    call check(nf90_put_att(ncid, cwpVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'iwp', nf90_float,dim2d, iwpVarID))
    call check(nf90_put_att(ncid, iwpVarID, "units", "kg/m^2"))
    call check(nf90_put_att(ncid, iwpVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'rwp', nf90_float,dim2d, rwpVarID))
    call check(nf90_put_att(ncid, rwpVarID, "units", "kg/m^2"))
    call check(nf90_put_att(ncid, rwpVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'swp', nf90_float,dim2d, swpVarID))
    call check(nf90_put_att(ncid, swpVarID, "units", "kg/m^2"))
    call check(nf90_put_att(ncid, swpVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'gwp', nf90_float,dim2d, gwpVarID))
    call check(nf90_put_att(ncid, gwpVarID, "units", "kg/m^2"))
    call check(nf90_put_att(ncid, gwpVarID, "missing_value", -9999))

    call check(nf90_def_var(ncid,'hwp', nf90_float,dim2d, hwpVarID))
    call check(nf90_put_att(ncid, hwpVarID, "units", "kg/m^2"))
    call check(nf90_put_att(ncid, hwpVarID, "missing_value", -9999))


    if (active) then

        dim3d = (/dlayerID,dlatID,dlonID/)
        call check(nf90_def_var(ncid,'height', nf90_float,dim3d, heightVarID))
        call check(nf90_put_att(ncid, heightVarID, "units", "m"))
        call check(nf90_put_att(ncid, heightVarID, "missing_value", -9999))

        dim4d = (/dfrqID,dlayerID,dlatID,dlonID/)

        call check(nf90_def_var(ncid,'Attenuation_Atmosphere', nf90_double,dim4d, AttAtmoVarID))
        call check(nf90_put_att(ncid, AttAtmoVarID, "units", attUnit))
        call check(nf90_put_att(ncid, AttAtmoVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Ze', nf90_double,dim4d, ZeVarID))
        call check(nf90_put_att(ncid, ZeVarID, "units", zeUnit))
        call check(nf90_put_att(ncid, ZeVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Attenuation_Hydrometeors', nf90_double,dim4d, AttHydroVarID))
        call check(nf90_put_att(ncid, AttHydroVarID, "units", attUnit))
        call check(nf90_put_att(ncid, AttHydroVarID, "missing_value", -9999))
        if ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments")) then

            call check(nf90_def_var(ncid,'Radar_SNR', nf90_double,dim4d, RadarSNRID))
            call check(nf90_put_att(ncid, RadarSNRID, "units", "dB"))
            call check(nf90_put_att(ncid, RadarSNRID, "missing_value", -9999))

            call check(nf90_def_var(ncid,'Radar_FallVelocity', nf90_double,dim4d, velVarID))
            call check(nf90_put_att(ncid, velVarID, "units", "m/s"))
            call check(nf90_put_att(ncid, velVarID, "missing_value", -9999))

            call check(nf90_def_var(ncid,'Radar_SpectralWidth', nf90_double,dim4d, swVarID))
            call check(nf90_put_att(ncid, swVarID, "units", "m/s"))
            call check(nf90_put_att(ncid, swVarID, "missing_value", -9999))

            call check(nf90_def_var(ncid,'Radar_Skewness', nf90_double,dim4d, skewVarID))
            call check(nf90_put_att(ncid, skewVarID, "units", "-"))
            call check(nf90_put_att(ncid, skewVarID, "missing_value", -9999))

            call check(nf90_def_var(ncid,'Radar_Kurtosis', nf90_double,dim4d, kurtVarID))
            call check(nf90_put_att(ncid, kurtVarID, "units", "-"))
            call check(nf90_put_att(ncid, kurtVarID, "missing_value", -9999))

            call check(nf90_def_var(ncid,'Radar_LeftSlope', nf90_double,dim4d, lSloVarID))
            call check(nf90_put_att(ncid, lSloVarID, "units", "dB/(m/s)"))
            call check(nf90_put_att(ncid, lSloVarID, "missing_value", -9999))

            call check(nf90_def_var(ncid,'Radar_RightSlope', nf90_double,dim4d, rSloVarID))
            call check(nf90_put_att(ncid, rSloVarID, "units", "dB/(m/s)"))
            call check(nf90_put_att(ncid, rSloVarID, "missing_value", -9999))

            call check(nf90_def_var(ncid,'Radar_LeftEdge', nf90_double,dim4d, lEdgVarID))
            call check(nf90_put_att(ncid, lEdgVarID, "units", "m/s"))
            call check(nf90_put_att(ncid, lEdgVarID, "missing_value", -9999))

            call check(nf90_def_var(ncid,'Radar_RightEdgee', nf90_double,dim4d, rEdgVarID))
            call check(nf90_put_att(ncid, rEdgVarID, "units", "m/s"))
            call check(nf90_put_att(ncid, rEdgVarID, "missing_value", -9999))

            call check(nf90_def_var(ncid,'Radar_Quality', nf90_int,dim4d, rQualVarID))
            call check(nf90_put_att(ncid, rQualVarID, "units", "bytes"))
            call check(nf90_put_att(ncid, rQualVarID, "description", "1st byte: aliasing; "// &
            "2nd byte: 2nd peak present; 7th: no peak found"))
            call check(nf90_put_att(ncid, rQualVarID, "missing_value", -9999))

            if (radar_mode == "spectrum") then
                call check(nf90_def_var(ncid,'Radar_Velocity', nf90_double,(/dnfftID/), RadarVelID))
                call check(nf90_put_att(ncid, RadarVelID, "units", "m/s"))
                call check(nf90_put_att(ncid, RadarVelID, "missing_value", -9999))

                dim5d = (/dnfftID,dfrqID,dlayerID,dlatID,dlonID/)
                call check(nf90_def_var(ncid,'Radar_Spectrum', nf90_double,dim5d, RadarSpecID))
                call check(nf90_put_att(ncid, RadarSpecID, "units", "dBz"))
                call check(nf90_put_att(ncid, RadarSpecID, "missing_value", -9999))
 
            end if
        end if
    end if

    if (passive) then
        dim6d = (/dstokesID,dfrqID,dangID,doutID,dlatID,dlonID/)
        call check(nf90_def_var(ncid,'tb', nf90_double,dim6d, tbVarID))
        call check(nf90_put_att(ncid, tbVarID, "units", "K"))
        call check(nf90_put_att(ncid, tbVarID, "missing_value", -9999))
    end if

    call check(nf90_enddef(ncid))

    call check(nf90_put_var(ncid, frequencyVarID, freqs(1:nfrq)))
    if (passive) then
        call check(nf90_put_var(ncid, anglesVarID, angles_deg))
    end if
    call check(nf90_put_var(ncid, isVarID, &
      RESHAPE( atmo_model_i, (/ atmo_ngridy, atmo_ngridx/), ORDER = (/2,1/))))
    call check(nf90_put_var(ncid, jsVarID, &
      RESHAPE( atmo_model_j, (/ atmo_ngridy, atmo_ngridx/), ORDER = (/2,1/))))
    call check(nf90_put_var(ncid, lonVarID, &
      RESHAPE( atmo_lon, (/ atmo_ngridy, atmo_ngridx/), ORDER = (/2,1/))))
    call check(nf90_put_var(ncid, latVarID, &
      RESHAPE( atmo_lat, (/ atmo_ngridy, atmo_ngridx/), ORDER = (/2,1/))))
    call check(nf90_put_var(ncid, lfracVarID, &
      RESHAPE( atmo_lfrac, (/ atmo_ngridy, atmo_ngridx/), ORDER = (/2,1/))))
    call check(nf90_put_var(ncid, iwvVarID,  &
      RESHAPE( atmo_iwv, (/ atmo_ngridy, atmo_ngridx/), ORDER = (/2,1/))))

    call check(nf90_put_var(ncid, cwpVarID, cwps))
    call check(nf90_put_var(ncid, iwpVarID, iwps))
    call check(nf90_put_var(ncid, rwpVarID, rwps))
    call check(nf90_put_var(ncid, swpVarID, swps))
    call check(nf90_put_var(ncid, gwpVarID, gwps))
    call check(nf90_put_var(ncid, hwpVarID, hwps))
    if (passive) then
        call check(nf90_put_var(ncid, tbVarID, tb))
    end if

    if (active) then                             !reshapeing needed due to Fortran's crazy Netcdf handling...
        call check(nf90_put_var(ncid, heightVarID, &
        RESHAPE( radar_hgt, (/ atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/3,2,1/))))
        call check(nf90_put_var(ncid, ZeVarID, &
        RESHAPE( Ze, (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, AttHydroVarID, &
        RESHAPE( Att_hydro, (/nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, AttAtmoVarID, &
        RESHAPE( Att_atmo, (/nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
        if ((radar_mode == "moments") .or.(radar_mode == "spectrum") ) then
            call check(nf90_put_var(ncid, velVarID, &
            RESHAPE( radar_moments(:,:,:,:,1), (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
            call check(nf90_put_var(ncid, swVarID, &
            RESHAPE( radar_moments(:,:,:,:,2), (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
            call check(nf90_put_var(ncid, skewVarID, &
            RESHAPE( radar_moments(:,:,:,:,3), (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
            call check(nf90_put_var(ncid, kurtVarID, &
            RESHAPE( radar_moments(:,:,:,:,4), (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
            call check(nf90_put_var(ncid, lSloVarID, &
            RESHAPE( radar_slopes(:,:,:,:,1), (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
            call check(nf90_put_var(ncid, rSloVarID, &
            RESHAPE( radar_slopes(:,:,:,:,2), (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))

            call check(nf90_put_var(ncid, lEdgVarID, &
            RESHAPE( radar_edge(:,:,:,:,1), (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
            call check(nf90_put_var(ncid, rEdgVarID, &
            RESHAPE( radar_edge(:,:,:,:,2), (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))

            call check(nf90_put_var(ncid, RadarSNRID, &
            RESHAPE( radar_snr, (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
            call check(nf90_put_var(ncid, rQualVarID, &
            RESHAPE( radar_quality(:,:,:,:), (/ nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/4,3,2,1/))))
            if (radar_mode == "spectrum") then
                call check(nf90_put_var(ncid, RadarVelID, radar_vel))
                call check(nf90_put_var(ncid, RadarSpecID, &
                RESHAPE( radar_spectra, (/ radar_nfft, nfrq, atmo_max_nlyrs, atmo_ngridy, atmo_ngridx/), ORDER = (/5,4,3,2,1/))))
            end if
        end if

    end if

    call check(nf90_close(ncid))

    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    
    return

contains


    subroutine check(status)

        integer, intent(in) :: status

        if(status /= nf90_noerr) then
            print *, "NC ERROR: ",trim(nf90_strerror(status))
            stop "Stopped"
        end if

        return

    end subroutine check

end subroutine write_nc_results
