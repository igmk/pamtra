subroutine hydrometeor_extinction(errorstatus,f,nx,ny,fi)

  use kinds
  use vars_atmosphere, only: nlyr, temp, q_hydro,&
      cwc_q, iwc_q, rwc_q, swc_q, gwc_q
  use settings, only: verbose, hydro_threshold
  use constants
  use descriptor_file
  use drop_size_dist
  use report_module
  implicit none

!   use mod_io_strings
!   use conversions
!   use tmat_snow_db
!   use tmat_rain_db
!   use vars_output, only: radar_spectra, radar_snr, radar_moments,&
! 	radar_quality, radar_slope, Ze, Att_hydro !output of the radar simulator for jacobian mode
!         use report_module
! 

  real(kind=dbl), intent(in) :: f
  integer, intent(in) ::  nx
  integer, intent(in) ::  ny
  integer, intent(in) ::  fi



  integer ::  nz
  integer :: ih
  
  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=40) :: nameOfRoutine = 'hydrometeor_extinction'
  
  if (verbose .gt. 1) print*, nx, ny, 'Entering hydrometeor_extinction'



!TMP
    q_hydro(1,:) = cwc_q(:)
    q_hydro(2,:) = iwc_q(:)
    q_hydro(3,:) = rwc_q(:)
    q_hydro(4,:) = swc_q(:)
    q_hydro(5,:) = gwc_q(:)
  


  grid_z: do nz = 1, nlyr  ! loop over all layers

     if (verbose .gt. 1) print*, 'Layer: ', nz

      hydros: do ih = 1,n_hydro

	if (q_hydro(ih,nz) >= hydro_threshold) then

      ! fill 0-D variable for the run_drop_size routine
	  hydro_name = hydro_name_arr(ih)
	  as_ratio   = as_ratio_arr(ih)
	  liq_ice    = liq_ice_arr(ih)
	  rho_ms     = rho_ms_arr(ih)
	  a_ms       = a_ms_arr(ih)
	  b_ms       = b_ms_arr(ih)
	  moment_in  = moment_in_arr(ih)
	  nbin       = nbin_arr(ih)
	  dist_name  = dist_name_arr(ih)
	  p_1        = p_1_arr(ih)
	  p_2        = p_2_arr(ih)
	  p_3        = p_3_arr(ih)
	  p_4        = p_4_arr(ih)
	  d_1        = d_1_arr(ih)
	  d_2        = d_2_arr(ih)

	  q_h        = q_hydro(ih,nz)
	  n_tot      = 0.
	  r_eff      = 0.
	  t          = temp(nz)
	  if (verbose >= 2) print*, ih, hydro_name
	  call run_drop_size_dist(err)


	if (err == 2) then
	  msg = 'Error in run_drop_size_dist'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
	end if

	end if
    end do hydros

  end do grid_z !end of cycle over the vertical layers



  if (verbose .gt. 1) print*, 'Exiting hydrometeor_extinction'
  return

end subroutine hydrometeor_extinction

