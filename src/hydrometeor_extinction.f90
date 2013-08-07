subroutine hydrometeor_extinction(errorstatus,f,nx,ny,fi)

  use kinds
  use vars_atmosphere, only: nlyr, temp, q_hydro, q_hum,&
      cwc_q, iwc_q, rwc_q, swc_q, gwc_q, press,&
      delta_hgt_lev, hydros_present
  use settings, only: verbose, hydro_threshold
  use constants
  use descriptor_file
  use drop_size_dist
  use report_module
  use scatProperties
  implicit none

!   use mod_io_strings
!   use conversions
!   use tmat_snow_db
!   use tmat_rain_db
!   use vars_output, only: radar_spectra, radar_snr, radar_moments,&
!        radar_quality, radar_slope, Ze, Att_hydro !output of the radar simulator for jacobian mode
!         use report_module
! 

  real(kind=dbl), intent(in) :: f
  real(kind=dbl) ::    scatter_matrix_scatcnv(nstokes,nummu,nstokes,nummu,4)
  real(kind=dbl) ::    extinct_matrix_scatcnv(nstokes,nstokes,nummu,2)
  real(kind=dbl) ::    emis_vector_scatcnv(nstokes,nummu,2)

  integer, intent(in) ::  nx
  integer, intent(in) ::  ny
  integer, intent(in) ::  fi

  integer ::  nz
  integer :: ih
  CHARACTER(len=64) :: scatfiles(nlyr)
  
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
  
  call allocate_scatProperties()

  grid_z: do nz = 1, nlyr  ! loop over all layers

    call prepare_rt3_scatProperties(nz)
    call prepare_rt4_scatProperties(nz)
    hydros_present(nz) = .false.
    
    if (verbose .gt. 1) print*, 'Layer: ', nz

    hydros: do ih = 1,n_hydro



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
      scat_name  = scat_name_arr(ih)
      p_1        = p_1_arr(ih)
      p_2        = p_2_arr(ih)
      p_3        = p_3_arr(ih)
      p_4        = p_4_arr(ih)
      d_1        = d_1_arr(ih)
      d_2        = d_2_arr(ih)

! Convert specific quantities [kg/kg] in absolute ones [kg/m3]
      q_h        = q2abs(q_hydro(ih,nz),temp(nz),press(nz),q_hum(nz),&
                   q_hydro(1,nz),q_hydro(2,nz),q_hydro(3,nz),q_hydro(4,nz),q_hydro(5,nz))
      n_tot      = 0.
      r_eff      = 0.
      t          = temp(nz)
      pressure   = press(nz)

      if (verbose >= 2) print*, ih, hydro_name

      if (q_h < hydro_threshold) then
	if (verbose >=3) print*, nx,ny,nz,ih,hydro_name, "q_h below threshold", q_h
	CYCLE
      end if

      hydros_present(nz) = .true.

      call allocateVars_drop_size_dist()

      call run_drop_size_dist(err)
      if (err == 2) then
	msg = 'Error in run_drop_size_dist'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	return
      end if

      call calc_scatProperties(err,f,nz)
      if (err == 2) then
	msg = 'Error in calc_scatProperties'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	return
      end if

      call deallocateVars_drop_size_dist()

    end do hydros

    call finalize_rt3_scatProperties(nz)

    !convert rt3 to rt4 input
    if (nlegen_coef>0) then

      call scatcnv(err,scatfiles(nz),nlegen_coef,legen_coef,kexttot(nz),salbedo,&
	scatter_matrix_scatcnv,extinct_matrix_scatcnv,emis_vector_scatcnv)
      if (err /= 0) then
	  msg = 'error in scatcnv!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
      end if   
      scattermatrix(nz,:,:,:,:,:) = scattermatrix(nz,:,:,:,:,:) + scatter_matrix_scatcnv
      extmatrix(nz,:,:,:,:) = extmatrix(nz,:,:,:,:) + extinct_matrix_scatcnv
      emisvec(nz,:,:,:) = emisvec(nz,:,:,:) + emis_vector_scatcnv
    end if

    if (active .and. hydros_present(nz)) then
      call radar_simulator(err,radar_spec, back(nz), kexttot(nz), f,&
	delta_hgt_lev(nz),nz,nx,ny,fi)
      if (err /= 0) then
	  msg = 'error in radar_simulator!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
      end if   
    end if


  end do grid_z !end of cycle over the vertical layers

  call deallocate_scatProperties()


  if (verbose .gt. 1) print*, 'Exiting hydrometeor_extinction'
  return

end subroutine hydrometeor_extinction

