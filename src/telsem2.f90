subroutine telsem2(errorstatus, month, lon, lat, freq, emissivity)
  use mod_mwatlas_nt_bin!, only: rttov_readmw_atlas,&
  !rttov_closemw_atlas, test_inputs, emis_interp_ind_sing
  use kinds
  use report_module
  use constants, only: rad2deg
  use settings, only: data_path, nummu, mu_values
  
  implicit none

  integer :: i
  real(kind=sgl), intent(in) :: lon, lat
  real(kind=dbl), intent(in) :: freq
  real(kind=dbl), dimension(2,nummu), intent(out) :: emissivity 
  !DECLARATIONS ==============================================

  !INPUT PARAMETERS
  integer(kind=long) :: imonth  !(1->12)
  character(len=2), intent(in) :: month
  character(len=80) :: dir  !directory of emis database
  real(kind=dbl) :: theta  !(0->60deg)
  real(kind=sgl) :: lon_tmp, lat_tmp
  real(kind=sgl) :: a1,a2

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg    
  character(22) :: nameOfRoutine = 'telsem2'

  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)  

  !====================================================
  !===read the atlas
  !====================================================
    
  dir=data_path(:len_trim(data_path))//'/emissivity/'

  call rttov_readmw_atlas(err,trim(dir),month)

  if (err /= 0) then
     msg = 'error in rttov_readmw_atlas'
     call report(err,msg, nameOfRoutine)
     errorstatus = err
     return
  end if

  read(month,'(I2)') imonth
  ! the data in the database as a range in longitude from 0 to 360. therefore we need to transform 
  ! negative longitudes to the range 180 to 359
  if (lon < 0.) then
    lon_tmp = lon + 360.
  else
    lon_tmp = lon
  end if
  if (lat > 87.5) then
    lat_tmp = 87.5
  else
    lat_tmp = lat
  end if
  do i = 1,nummu
    theta = acos(mu_values(i))*rad2deg
    !====================================================
    !===sample on a single pixel
    !====================================================
    !===test the inputs
    call test_inputs(err,imonth,lat_tmp,lon_tmp,theta,freq)
    if (err /= 0) then
        msg = 'error in test_inputs'
        call report(err,msg, nameOfRoutine)
        call rttov_closemw_atlas()
        errorstatus = err
        return
    end if
    call emis_interp_ind_sing(err,lat_tmp,lon_tmp,theta,freq,emissivity(1,i),emissivity(2,i))
    if (err /= 0) then
        msg = 'error in emis_interp_ind_sing'
        call report(err,msg, nameOfRoutine)
        call rttov_closemw_atlas()
        errorstatus = err
        return
    end if
  end do
    
  call rttov_closemw_atlas()

  errorstatus = err

  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

end subroutine telsem2

