subroutine write_nc_results(nc_file)

  use kinds
  use vars_output
  use vars_atmosphere, only: ngridx, ngridy,nlyr
  use netcdf
  use nml_params, only: active, passive, creator
  implicit none

  integer :: ncid
  integer :: dlonID, dlatID, dangID, dfreID, doutID, dstokesID, dlayerID
  integer :: isVarID, jsVarID, lonVarID, latVarID, lfracVarID, t_gVarID, wind10uVarID, wind10vVarID, iwvVarID, cwpVarID,&
	     iwpVarID, rwpVarID, swpVarID, gwpVarID, flux_upVarID, flux_downVarID, &
	     tbVarID, heightVarID, ZeVarID, PiaAtmoVarID, PiaHydroVarID
  

  integer :: nang = 32, nfre = 1, nout = 2, nstokes = 2

  integer, dimension(2) :: dim2d
  integer, dimension(3) :: dim3d
  integer, dimension(4) :: dim4d
  integer, dimension(5) :: dim5d

  integer :: today(3), now(3)

  character(100) :: nc_file, timestring, user


  call check(nf90_create(path=nc_file,cmode=nf90_noclobber,ncid=ncid))


  ! for netcdf history get meta data
  call idate(today)   ! today(1)=day, (2)=month, (3)=year
  call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
  write (timestring , "(i2.2, '/', i2.2, '/', i4.4, ' ',  i2.2, ':', i2.2, ':', i2.2)") &
	today(2), today(1), today(3), now
  ! write meta data
  call check(nf90_put_att(ncid,nf90_global, "history", "Created with Fortran by "//trim(creator)//&
	" (University of Cologne, IGMK) at "//timestring))

  !make dimensions
  call check(nf90_def_dim(ncid, 'nlon', ngridx, dlonID))
  call check(nf90_def_dim(ncid, 'nlat', ngridy, dlatID))
  call check(nf90_def_dim(ncid, 'nang', nang, dangID))
  call check(nf90_def_dim(ncid, 'nfre', nfre, dfreID))
  call check(nf90_def_dim(ncid, 'nout', nout, doutID))
  call check(nf90_def_dim(ncid, 'nstokes', nstokes, dstokesID))
if (active) then
  call check(nf90_def_dim(ncid, 'nlyr', nlyr, dlayerID))
end if

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

  call check(nf90_def_var(ncid,'t_g', nf90_float,dim2d, t_gVarID))
  call check(nf90_put_att(ncid, t_gVarID, "units", "K"))
  call check(nf90_put_att(ncid, t_gVarID, "missing_value", -9999))


  call check(nf90_def_var(ncid,'wind10u', nf90_float,dim2d, wind10uVarID))
  call check(nf90_put_att(ncid, wind10uVarID, "units", "m/s"))
  call check(nf90_put_att(ncid, wind10uVarID, "missing_value", -9999))

  call check(nf90_def_var(ncid,'wind10v', nf90_float,dim2d, wind10vVarID))
  call check(nf90_put_att(ncid, wind10vVarID, "units", "m/s"))
  call check(nf90_put_att(ncid, wind10vVarID, "missing_value", -9999))


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


if (active) then
  dim3d = (/dlayerID,dlatID,dlonID/)
  call check(nf90_def_var(ncid,'height', nf90_float,dim3d, heightVarID))
  call check(nf90_put_att(ncid, heightVarID, "units", "m"))
  call check(nf90_put_att(ncid, heightVarID, "missing_value", -9999))

  call check(nf90_def_var(ncid,'Ze', nf90_float,dim3d, ZeVarID))
  call check(nf90_put_att(ncid, ZeVarID, "units", "dBz"))
  call check(nf90_put_att(ncid, ZeVarID, "missing_value", -9999))

  call check(nf90_def_var(ncid,'PIA Hydrometeors', nf90_float,dim3d, PiaHydroVarID))
  call check(nf90_put_att(ncid, PiaHydroVarID, "units", "dB"))
  call check(nf90_put_att(ncid, PiaHydroVarID, "missing_value", -9999))

  call check(nf90_def_var(ncid,'PIA Atmosphere', nf90_float,dim3d, PiaAtmoVarID))
  call check(nf90_put_att(ncid, PiaAtmoVarID, "units", "dB"))
  call check(nf90_put_att(ncid, PiaAtmoVarID, "missing_value", -9999))

end if

  dim4d = (/dstokesID,doutID,dlatID,dlonID/)
  call check(nf90_def_var(ncid,'flux_up', nf90_double,dim4d, flux_upVarID))
  call check(nf90_put_att(ncid, flux_upVarID, "units", "J s^−1 m^−2 sr^−1 m^−1"))
  call check(nf90_put_att(ncid, flux_upVarID, "missing_value", -9999))

  call check(nf90_def_var(ncid,'flux_down', nf90_double,dim4d, flux_downVarID))
  call check(nf90_put_att(ncid, flux_downVarID, "units", "J s^−1 m^−2 sr^−1 m^−1"))
  call check(nf90_put_att(ncid, flux_downVarID, "missing_value", -9999))

  dim5d = (/dstokesID,dangID,doutID,dlatID,dlonID/)
  call check(nf90_def_var(ncid,'tb', nf90_double,dim5d, tbVarID))
  call check(nf90_put_att(ncid, tbVarID, "units", "K"))
  call check(nf90_put_att(ncid, tbVarID, "missing_value", -9999))

  call check(nf90_enddef(ncid))
!  call check(nf90_inq_varid(ncid, 'longitude', VarId))


  call check(nf90_put_var(ncid, isVarID, is))
  call check(nf90_put_var(ncid, jsVarID, js))
  call check(nf90_put_var(ncid, lonVarID, lons))
  call check(nf90_put_var(ncid, latVarID, lats))
  call check(nf90_put_var(ncid, lfracVarID, lfracs))
  call check(nf90_put_var(ncid, t_gVarID, t_g))
  call check(nf90_put_var(ncid, wind10uVarID, w10u))
  call check(nf90_put_var(ncid, wind10vVarID, w10v))
  call check(nf90_put_var(ncid, iwvVarID, iwvs))
  call check(nf90_put_var(ncid, cwpVarID, cwps))
  call check(nf90_put_var(ncid, iwpVarID, iwps))
  call check(nf90_put_var(ncid, rwpVarID, rwps))
  call check(nf90_put_var(ncid, swpVarID, swps))
  call check(nf90_put_var(ncid, gwpVarID, gwps))
  call check(nf90_put_var(ncid, flux_upVarID, flux_up))
  call check(nf90_put_var(ncid, flux_downVarID, flux_down))
  call check(nf90_put_var(ncid, tbVarID, tb))

if (active) then                             !reshapeing needed due to Fortran's crazy Netcdf handling...
  call check(nf90_put_var(ncid, heightVarID, RESHAPE( hgt, (/ nlyr, ngridy, ngridx/), ORDER = (/3,2,1/))))
  call check(nf90_put_var(ncid, ZeVarID, RESHAPE( Ze, (/ nlyr, ngridy, ngridx/), ORDER = (/3,2,1/))))
  call check(nf90_put_var(ncid, PiaHydroVarID, RESHAPE( PIA_hydro, (/ nlyr, ngridy, ngridx/), ORDER = (/3,2,1/))))
  call check(nf90_put_var(ncid, PiaAtmoVarID, RESHAPE( PIA_atmo, (/ nlyr, ngridy, ngridx/), ORDER = (/3,2,1/))))
end if

  call check(nf90_close(ncid))

  deallocate(lons,lats,lfracs,t_g,w10u,w10v,iwvs,cwps,iwps,rwps,swps,gwps,flux_up,flux_down,tb,hgt,Ze,PIA_atmo, PIA_hydro)

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
