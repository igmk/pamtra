module nml_params
  ! Description:
  ! Definition of all name list paramters for pamtra
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  0.1   17/11/2009    creation of file 
  use kinds
 
  implicit none

  integer :: verbose, n_moments

  real(kind=dbl) :: obs_height     ! upper level output height [m] (> 100000. for satellite)

  real(kind=dbl) :: emissivity

  real(kind=dbl) :: N_0snowDsnow, N_0grauDgrau, N_0rainD, SP

  logical :: lphase_flag, &        ! flag for phase function calculation
	     lgas_extinction, &    ! gas extinction desired
	     lhyd_extinction, &    ! hydrometeor extinction desired
	     grid_calc, &	   ! calculate with dispatch or not
	     write_nc		   ! write netcdf output

  character(5) :: EM_snow, EM_grau
  character(3) :: SD_snow, SD_grau, SD_rain

  character(3) :: gas_mod

end module nml_params
