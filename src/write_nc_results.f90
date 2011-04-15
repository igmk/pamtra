subroutine write_nc_results(nc_file)

  use kinds
  use vars_output
  use vars_atmosphere, only: ngridx, ngridy
  use netcdf

  implicit none

  integer :: ncid
  integer :: dlonID, dlatID, dangID, dfreID, doutID, dstokesID
  integer :: isVarID, jsVarID, lonVarID, latVarID, lfracVarID, wind10VarID, iwvVarID, cwpVarID,&
	     iwpVarID, rwpVarID, swpVarID, gwpVarID, flux_upVarID, flux_downVarID, &
	     tb_upVarID, tb_downVarID
  

  integer :: nang = 16, nfre =1, nout = 2, nstokes = 2

  integer, dimension(2) :: dim2d
  integer, dimension(4) :: dim4d
  integer, dimension(5) :: dim5d

  character(100) :: nc_file

  call check(nf90_create(path='/net/roumet/mech/pamtra/output/'//nc_file,cmode=nf90_noclobber,ncid=ncid))

  call check(nf90_def_dim(ncid, 'nlon', ngridx, dlonID))
  call check(nf90_def_dim(ncid, 'nlat', ngridy, dlatID))
  call check(nf90_def_dim(ncid, 'nang', nang, dangID))
  call check(nf90_def_dim(ncid, 'nfre', nfre, dfreID))
  call check(nf90_def_dim(ncid, 'nout', nout, doutID))
  call check(nf90_def_dim(ncid, 'nstokes', nstokes, dstokesID))

  dim2d = (/dlatID,dlonID/)
!  call put_2d_var(ncid,'longitude',lons,2,/ngridx,ngridy/)
  call check(nf90_def_var(ncid,'model_i', nf90_int,dim2d, isVarID))
  call check(nf90_def_var(ncid,'model_j', nf90_int,dim2d, jsVarID))
  call check(nf90_def_var(ncid,'longitude', nf90_float,dim2d, lonVarID))
  call check(nf90_def_var(ncid,'latitude', nf90_float,dim2d, latVarID))
  call check(nf90_def_var(ncid,'lfrac', nf90_float,dim2d, lfracVarID))
  call check(nf90_def_var(ncid,'wind10', nf90_float,dim2d, wind10VarID))
  call check(nf90_def_var(ncid,'iwv', nf90_float,dim2d, iwvVarID))
  call check(nf90_def_var(ncid,'cwp', nf90_float,dim2d, cwpVarID))
  call check(nf90_def_var(ncid,'iwp', nf90_float,dim2d, iwpVarID))
  call check(nf90_def_var(ncid,'rwp', nf90_float,dim2d, rwpVarID))
  call check(nf90_def_var(ncid,'swp', nf90_float,dim2d, swpVarID))
  call check(nf90_def_var(ncid,'gwp', nf90_float,dim2d, gwpVarID))

  dim4d = (/dstokesID,doutID,dlatID,dlonID/)
  call check(nf90_def_var(ncid,'flux_up', nf90_double,dim4d, flux_upVarID))
  call check(nf90_def_var(ncid,'flux_down', nf90_double,dim4d, flux_downVarID))

  dim5d = (/dstokesID,dangID,doutID,dlatID,dlonID/)
  call check(nf90_def_var(ncid,'tb_up', nf90_double,dim5d, tb_upVarID))
  call check(nf90_def_var(ncid,'tb_down', nf90_double,dim5d, tb_downVarID))
  call check(nf90_enddef(ncid))
!  call check(nf90_inq_varid(ncid, 'longitude', VarId))
  call check(nf90_put_var(ncid, isVarID, is))
  call check(nf90_put_var(ncid, jsVarID, js))
  call check(nf90_put_var(ncid, lonVarID, lons))
  call check(nf90_put_var(ncid, latVarID, lats))
  call check(nf90_put_var(ncid, lfracVarID, lfracs))
  call check(nf90_put_var(ncid, wind10VarID, w10s))
  call check(nf90_put_var(ncid, iwvVarID, iwvs))
  call check(nf90_put_var(ncid, cwpVarID, cwps))
  call check(nf90_put_var(ncid, iwpVarID, iwps))
  call check(nf90_put_var(ncid, rwpVarID, rwps))
  call check(nf90_put_var(ncid, swpVarID, swps))
  call check(nf90_put_var(ncid, gwpVarID, gwps))
  call check(nf90_put_var(ncid, flux_upVarID, flux_up))
  call check(nf90_put_var(ncid, flux_downVarID, flux_down))
  call check(nf90_put_var(ncid, tb_upVarID, tb_up))
  call check(nf90_put_var(ncid, tb_downVarID, tb_down))

  call check(nf90_close(ncid))

  deallocate(lons,lats,lfracs,w10s,iwvs,cwps,iwps,rwps,swps,gwps,flux_up,flux_down,tb_up,tb_down)

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