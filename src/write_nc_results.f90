subroutine write_nc_results(nc_file)

  use kinds
  use vars_output
  use vars_atmosphere, only: ngridx, ngridy,nlyr,freqs,nfrq, year, month, day, time
  use netcdf
  use nml_params, only: active, passive, creator, verbose, zeSplitUp, n_moments
  implicit none

  integer :: ncid
  integer :: dlonID, dlatID, dangID, dfrqID, doutID, dstokesID, dlayerID

  integer :: isVarID, jsVarID, lonVarID, latVarID, lfracVarID, iwvVarID, cwpVarID,&
       iwpVarID, rwpVarID, swpVarID, gwpVarID, hwpVarID, &
       tbVarID, heightVarID, &
       ZeVarID,ZeCwVarID,ZeRrVarID,ZeCiVarID,ZeSnVarID,ZeGrVarID,ZeHaVarID, &
       AttAtmoVarID, AttHydroVarID,&
       AttCwVarID, AttRrVarID, AttCiVarID, AttSnVarID, AttGrVarID, AttHaVarID, &
       frequencyVarID, anglesVarID
          
          
  integer :: nang = 32, nout = 2, nstokes = 2

  integer, dimension(2) :: dim2d
  integer, dimension(3) :: dim3d
  integer, dimension(4) :: dim4d
  integer, dimension(6) :: dim6d

  integer :: today(3), now(3)

  character(300) :: nc_file, timestring, user
  character(40) ::gitVersion,gitHash

  if (verbose .gt. 0) print*,"writing: ", nc_file
  !get git data
  call versionNumber(gitVersion,gitHash)

  call check(nf90_create(path=nc_file,cmode=nf90_clobber,ncid=ncid))

  ! for netcdf history get meta data
  call idate(today)   ! today(1)=day, (2)=month, (3)=year
  call itime(now)     ! now(1)=hour, (2)=minute, (3)=second

  write (timestring , "(i2.2, '/', i2.2, '/', i4.4, ' ',  i2.2, ':', i2.2, ':', i2.2)") &
       today(2), today(1), today(3), now
  ! write meta data
  call check(nf90_put_att(ncid,nf90_global, "history", "Created with Pamtra (Version: "//trim(gitVersion)// &
       ", Git Hash: "//trim(gitHash)//")  by "//trim(creator)//" (University of Cologne, IGMK) at "//timestring))
  call check(nf90_put_att(ncid,nf90_global, "data_time",year//"/"//month//"/"//day//"-"//time(1:2)//":"//time(3:4)))

  !make dimensions
  call check(nf90_def_dim(ncid, 'nlon', ngridx, dlonID))
  call check(nf90_def_dim(ncid, 'nlat', ngridy, dlatID))
  call check(nf90_def_dim(ncid, 'nfreq', nfrq, dfrqID))
  if (passive) then
     call check(nf90_def_dim(ncid, 'nang', nang, dangID))
     call check(nf90_def_dim(ncid, 'nout', nout, doutID))
     call check(nf90_def_dim(ncid, 'nstokes', nstokes, dstokesID))
  end if
  if (active) then
     call check(nf90_def_dim(ncid, 'nlyr', nlyr, dlayerID))
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

     if (zeSplitUp) then

        call check(nf90_def_var(ncid,'Ze_cloud water', nf90_double,dim4d, ZeCwVarID))
        call check(nf90_put_att(ncid, ZeCwVarID, "units", "dBz"))
        call check(nf90_put_att(ncid, ZeCwVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Ze_rain_water', nf90_double,dim4d, ZeRrVarID))
        call check(nf90_put_att(ncid, ZeRrVarID, "units", "dBz"))
        call check(nf90_put_att(ncid, ZeRrVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Ze_cloud_ice', nf90_double,dim4d, ZeCiVarID))
        call check(nf90_put_att(ncid, ZeCiVarID, "units", "dBz"))
        call check(nf90_put_att(ncid, ZeCiVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Ze_snow', nf90_double,dim4d, ZeSnVarID))
        call check(nf90_put_att(ncid, ZeSnVarID, "units", "dBz"))
        call check(nf90_put_att(ncid, ZeSnVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Ze_graupel', nf90_double,dim4d, ZeGrVarID))
        call check(nf90_put_att(ncid, ZeGrVarID, "units", "dBz"))
        call check(nf90_put_att(ncid, ZeGrVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Attenuation_cloud_water', nf90_float,dim4d, AttCwVarID))
        call check(nf90_put_att(ncid, AttCwVarID, "units", "dB"))
        call check(nf90_put_att(ncid, AttCwVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Attenuation_rain', nf90_float,dim4d, AttRrVarID))
        call check(nf90_put_att(ncid, AttRrVarID, "units", "dB"))
        call check(nf90_put_att(ncid, AttRrVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Attenuation_cloud_ice', nf90_float,dim4d, AttCiVarID))
        call check(nf90_put_att(ncid, AttCiVarID, "units", "dB"))
        call check(nf90_put_att(ncid, AttCiVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Attenuation_snow', nf90_float,dim4d, AttSnVarID))
        call check(nf90_put_att(ncid, AttSnVarID, "units", "dB"))
        call check(nf90_put_att(ncid, AttSnVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Attenuation_graupel', nf90_float,dim4d, AttGrVarID))
        call check(nf90_put_att(ncid, AttGrVarID, "units", "dB"))
        call check(nf90_put_att(ncid, AttGrVarID, "missing_value", -9999))

        if (n_moments .eq. 2) then
           call check(nf90_def_var(ncid,'Ze_hail', nf90_double,dim4d, ZeHaVarID))
           call check(nf90_put_att(ncid, ZeHaVarID, "units", "dBz"))
           call check(nf90_put_att(ncid, ZeHaVarID, "missing_value", -9999))
       
           call check(nf90_def_var(ncid,'Attenuation_hail', nf90_float,dim4d, AttHaVarID))
           call check(nf90_put_att(ncid, AttHaVarID, "units", "dB"))
           call check(nf90_put_att(ncid, AttHaVarID, "missing_value", -9999))
        end if
        
     else
     
        call check(nf90_def_var(ncid,'Ze', nf90_double,dim4d, ZeVarID))
        call check(nf90_put_att(ncid, ZeVarID, "units", "dBz"))
        call check(nf90_put_att(ncid, ZeVarID, "missing_value", -9999))

        call check(nf90_def_var(ncid,'Attenuation_Hydrometeors', nf90_float,dim4d, AttHydroVarID))
        call check(nf90_put_att(ncid, AttHydroVarID, "units", "dB"))
        call check(nf90_put_att(ncid, AttHydroVarID, "missing_value", -9999))

     end if


     call check(nf90_def_var(ncid,'Attenuation_Atmosphere', nf90_float,dim4d, AttAtmoVarID))
     call check(nf90_put_att(ncid, AttAtmoVarID, "units", "dB"))
     call check(nf90_put_att(ncid, AttAtmoVarID, "missing_value", -9999))


  end if

  if (passive) then
     dim6d = (/dstokesID,dfrqID,dangID,doutID,dlatID,dlonID/)
     call check(nf90_def_var(ncid,'tb', nf90_double,dim6d, tbVarID))
     call check(nf90_put_att(ncid, tbVarID, "units", "K"))
     call check(nf90_put_att(ncid, tbVarID, "missing_value", -9999))
  end if

  call check(nf90_enddef(ncid))
  !  call check(nf90_inq_varid(ncid, 'longitude', VarId))



  call check(nf90_put_var(ncid, frequencyVarID, freqs(1:nfrq)))
  if (passive) then
	  call check(nf90_put_var(ncid, anglesVarID, angles_deg))
  end if
  call check(nf90_put_var(ncid, isVarID, is))
  call check(nf90_put_var(ncid, jsVarID, js))
  call check(nf90_put_var(ncid, lonVarID, lons))
  call check(nf90_put_var(ncid, latVarID, lats))
  call check(nf90_put_var(ncid, lfracVarID, lfracs))
  call check(nf90_put_var(ncid, iwvVarID, iwvs))
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
          RESHAPE( hgt, (/ nlyr, ngridy, ngridx/), ORDER = (/3,2,1/))))

     if (zeSplitUp) then
        call check(nf90_put_var(ncid, ZeCwVarID, &
          RESHAPE( Ze_cw, (/ nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, ZeRrVarID, &
          RESHAPE( Ze_rr, (/ nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, ZeCiVarID, &
          RESHAPE( Ze_ci, (/ nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, ZeSnVarID, &
          RESHAPE( Ze_sn, (/ nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, ZeGrVarID, &
          RESHAPE( Ze_gr, (/ nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
          
        call check(nf90_put_var(ncid, AttCwVarID, &
          RESHAPE( Att_cw, (/nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, AttRrVarID, &
          RESHAPE( Att_rr, (/nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, AttCiVarID, &
          RESHAPE( Att_ci, (/nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, AttSnVarID, &
          RESHAPE( Att_sn, (/nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, AttGrVarID, &
          RESHAPE( Att_gr, (/nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        if (n_moments .eq. 2) then
           call check(nf90_put_var(ncid, ZeHaVarID, &
             RESHAPE( Ze_ha, (/ nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
           call check(nf90_put_var(ncid, AttHaVarID, &
             RESHAPE( Att_ha, (/nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        end if           
     else
        call check(nf90_put_var(ncid, ZeVarID, &
          RESHAPE( Ze, (/ nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
        call check(nf90_put_var(ncid, AttHydroVarID, &
          RESHAPE( Att_hydro, (/nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
     end if




     call check(nf90_put_var(ncid, AttAtmoVarID, &
          RESHAPE( Att_atmo, (/nfrq, nlyr, ngridy, ngridx/), ORDER = (/4,3,2,1/))))
  end if

  call check(nf90_close(ncid))

  return

contains

  !   subroutine put_2d_var(ncid,varname,var,ndims,dims)
  ! 
  !   use kinds 
  !   implicit none
  ! 
  !   integer :: ncid, ndims, VarID
  !   integer, dimension(:) :: dims
  !   real(kind=dbl) , dimension(:) :: var
  !   character(:) :: varname
  ! 
  ! 
  !   return
  ! 
  !   end subroutine put_2d_var

  subroutine check(status)

    integer, intent(in) :: status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if

    return

  end subroutine check

end subroutine write_nc_results
