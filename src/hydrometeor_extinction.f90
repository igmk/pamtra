subroutine hydrometeor_extinction(errorstatus,f,nx,ny,fi)

  use kinds
  use vars_atmosphere, only: nlyr, temp, q_hydro, q_hum,&
      cwc_q, iwc_q, rwc_q, swc_q, gwc_q, press,&
      delta_hgt_lev, hydros_present
  use settings, only: verbose, hydro_threshold, save_psd
  use constants
  use descriptor_file
  use drop_size_dist
  use report_module
  use scatProperties
  use vars_output, only: psd_area, psd_d_bound, psd_f, psd_mass

  implicit none

!   use mod_io_strings
!   use conversions
!   use tmat_snow_db
!   use tmat_rain_db
!   use vars_output, only: radar_spectra, radar_snr, radar_moments,&
!        radar_quality, radar_slopes, Ze, Att_hydro !output of the radar simulator for jacobian mode
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
      ! start with the ones which might come from a 4D python field.
      if (PRODUCT(SHAPE(as_ratio_arr)) == n_hydro) then
        as_ratio   = as_ratio_arr(1,1,1,ih)
      else
        as_ratio   = as_ratio_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(rho_ms_arr)) == n_hydro) then
        rho_ms   = rho_ms_arr(1,1,1,ih)
      else
        rho_ms   = rho_ms_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(a_ms_arr)) == n_hydro) then
        a_ms   = a_ms_arr(1,1,1,ih)
      else
        a_ms   = a_ms_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(b_ms_arr)) == n_hydro) then
        b_ms   = b_ms_arr(1,1,1,ih)
      else
        b_ms   = b_ms_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(alpha_as_arr)) == n_hydro) then
        alpha_as   = alpha_as_arr(1,1,1,ih)
      else
        alpha_as   = alpha_as_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(beta_as_arr)) == n_hydro) then
        beta_as   = beta_as_arr(1,1,1,ih)
      else
        beta_as   = beta_as_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(nbin_arr)) == n_hydro) then
        nbin   = nbin_arr(1,1,1,ih)
      else
        nbin   = nbin_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(p_1_arr)) == n_hydro) then
        p_1   = p_1_arr(1,1,1,ih)
      else
        p_1   = p_1_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(p_2_arr)) == n_hydro) then
        p_2   = p_2_arr(1,1,1,ih)
      else
        p_2   = p_2_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(p_3_arr)) == n_hydro) then
        p_3   = p_3_arr(1,1,1,ih)
      else
        p_3   = p_3_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(p_4_arr)) == n_hydro) then
        p_4   = p_4_arr(1,1,1,ih)
      else
        p_4   = p_4_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(d_1_arr)) == n_hydro) then
        d_1   = d_1_arr(1,1,1,ih)
      else
        d_1   = d_1_arr(nx,ny,nz,ih)
      end if 

      if (PRODUCT(SHAPE(d_2_arr)) == n_hydro) then
        d_2   = d_2_arr(1,1,1,ih)
      else
        d_2   = d_2_arr(nx,ny,nz,ih)
      end if 

      !these ones are fore sure 1D
      hydro_name = hydro_name_arr(ih)
      liq_ice    = liq_ice_arr(ih)
      moment_in  = moment_in_arr(ih)
      dist_name  = dist_name_arr(ih)
      scat_name  = scat_name_arr(ih)
      vel_size_mod  = vel_size_mod_arr(ih)


! Convert specific quantities [kg/kg] in absolute ones [kg/m3]
      q_h        = q2abs(q_hydro(ih,nz),temp(nz),press(nz),q_hum(nz),&
                   q_hydro(1,nz),q_hydro(2,nz),q_hydro(3,nz),q_hydro(4,nz),q_hydro(5,nz))
      n_tot      = 0.
      r_eff      = 0.
      layer_t    = temp(nz)
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

      if (save_psd) then
        psd_d_bound(nx,ny,nz,ih,1:nbin+1) =  d_bound_ds
        psd_f(nx,ny,nz,ih,1:nbin+1) = f_ds
        psd_mass(nx,ny,nz,ih,1:nbin+1) = mass_ds
        psd_area(nx,ny,nz,ih,1:nbin+1) = area_ds
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

