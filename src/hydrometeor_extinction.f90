subroutine hydrometeor_extinction(errorstatus)

  use kinds
  use vars_atmosphere, only: atmo_nlyrs, atmo_temp, atmo_q_hum,&
      atmo_press,&
      atmo_delta_hgt_lev, atmo_hydro_q, atmo_hydro_reff, atmo_hydro_n
  use vars_rt, only:rt_hydros_present
  use settings, only: verbose, hydro_threshold, save_psd
  use constants
  use descriptor_file
  use drop_size_dist
  use report_module
  use scatProperties
  use vars_output, only: psd_area, psd_d_bound, psd_f, psd_mass
  use vars_index, only: i_x,i_y, i_z

  implicit none

!   use mod_io_strings
!   use conversions
!   use tmat_snow_db
!   use tmat_rain_db
!   use vars_output, only: radar_spectra, radar_snr, radar_moments,&
!        radar_quality, radar_slopes, Ze, Att_hydro !output of the radar simulator for jacobian mode
!         use report_module
! 

  real(kind=dbl) ::    scatter_matrix_scatcnv(nstokes,nummu,nstokes,nummu,4)
  real(kind=dbl) ::    extinct_matrix_scatcnv(nstokes,nstokes,nummu,2)
  real(kind=dbl) ::    emis_vector_scatcnv(nstokes,nummu,2)

  integer :: ih !index hydrometeor
  CHARACTER(len=64), dimension(atmo_nlyrs(i_x,i_y)) :: scatfiles
  
  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=40) :: nameOfRoutine = 'hydrometeor_extinction'
  
  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)



  call allocate_scatProperties()

  grid_z: do i_z = 1, atmo_nlyrs(i_x,i_y)  ! loop over all layers

    call prepare_rt3_scatProperties()
    call prepare_rt4_scatProperties()
    rt_hydros_present(i_z) = .false.
    
    if (verbose .gt. 1) print*, 'Layer: ', i_z

    hydros: do ih = 1,n_hydro



      ! fill 0-D variable for the run_drop_size routine
      ! start with the ones which might come from a 4D python field.
      if (PRODUCT(SHAPE(as_ratio_arr)) == n_hydro) then
        as_ratio   = as_ratio_arr(1,1,1,ih)
      else
        as_ratio   = as_ratio_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(rho_ms_arr)) == n_hydro) then
        rho_ms   = rho_ms_arr(1,1,1,ih)
      else
        rho_ms   = rho_ms_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(a_ms_arr)) == n_hydro) then
        a_ms   = a_ms_arr(1,1,1,ih)
      else
        a_ms   = a_ms_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(b_ms_arr)) == n_hydro) then
        b_ms   = b_ms_arr(1,1,1,ih)
      else
        b_ms   = b_ms_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(alpha_as_arr)) == n_hydro) then
        alpha_as   = alpha_as_arr(1,1,1,ih)
      else
        alpha_as   = alpha_as_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(beta_as_arr)) == n_hydro) then
        beta_as   = beta_as_arr(1,1,1,ih)
      else
        beta_as   = beta_as_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(nbin_arr)) == n_hydro) then
        nbin   = nbin_arr(1,1,1,ih)
      else
        nbin   = nbin_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(p_1_arr)) == n_hydro) then
        p_1   = p_1_arr(1,1,1,ih)
      else
        p_1   = p_1_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(p_2_arr)) == n_hydro) then
        p_2   = p_2_arr(1,1,1,ih)
      else
        p_2   = p_2_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(p_3_arr)) == n_hydro) then
        p_3   = p_3_arr(1,1,1,ih)
      else
        p_3   = p_3_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(p_4_arr)) == n_hydro) then
        p_4   = p_4_arr(1,1,1,ih)
      else
        p_4   = p_4_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(d_1_arr)) == n_hydro) then
        d_1   = d_1_arr(1,1,1,ih)
      else
        d_1   = d_1_arr(i_x,i_y,i_z,ih)
      end if 

      if (PRODUCT(SHAPE(d_2_arr)) == n_hydro) then
        d_2   = d_2_arr(1,1,1,ih)
      else
        d_2   = d_2_arr(i_x,i_y,i_z,ih)
      end if 

      !these ones are fore sure 1D
      hydro_name = hydro_name_arr(ih)
      liq_ice    = liq_ice_arr(ih)
      moment_in  = moment_in_arr(ih)
      dist_name  = dist_name_arr(ih)
      scat_name  = scat_name_arr(ih)
      vel_size_mod  = vel_size_mod_arr(ih)


! Convert specific quantities [kg/kg] in absolute ones [kg/m3]
!       q_h        = q2abs(q_hydro(ih,i_z),atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z),q_hum(i_z),&
!                    q_hydro(1,i_z),q_hydro(2,i_z),q_hydro(3,i_z),q_hydro(4,i_z),q_hydro(5,i_z))
      q_h        = q2abs(atmo_hydro_q(i_x,i_y,i_z, ih),atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z),&
                  atmo_q_hum(i_x,i_y,i_z),0._dbl)
      n_tot      = q2abs(atmo_hydro_n(i_x,i_y,i_z, ih),atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z),&
                  atmo_q_hum(i_x,i_y,i_z),0._dbl)
      r_eff      = atmo_hydro_reff(i_x,i_y,i_z, ih)
      layer_t    = atmo_temp(i_x,i_y,i_z)
      pressure   = atmo_press(i_x,i_y,i_z)

      if (verbose >= 2) print*, ih, hydro_name

      if (q_h < hydro_threshold) then
	if (verbose >=3) print*, i_x,i_y,i_z,ih,hydro_name, "q_h below threshold", q_h
	CYCLE
      end if

      rt_hydros_present(i_z) = .true.

      call allocateVars_drop_size_dist()

      call run_drop_size_dist(err)
      if (err == 2) then
	msg = 'Error in run_drop_size_dist'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	return
      end if

      if (save_psd) then
        psd_d_bound(i_x,i_y,i_z,ih,1:nbin+1) =  d_bound_ds
        psd_f(i_x,i_y,i_z,ih,1:nbin+1) = f_ds
        psd_mass(i_x,i_y,i_z,ih,1:nbin+1) = mass_ds
        psd_area(i_x,i_y,i_z,ih,1:nbin+1) = area_ds
      end if

      call calc_scatProperties(err)
      if (err == 2) then
	msg = 'Error in calc_scatProperties'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	return
      end if

      call deallocateVars_drop_size_dist()

    end do hydros

    call finalize_rt3_scatProperties()

    !convert rt3 to rt4 input
    if (nlegen_coef>0) then

      call scatcnv(err,scatfiles(i_z),nlegen_coef,legen_coef,rt_kexttot(i_z),salbedo,&
	scatter_matrix_scatcnv,extinct_matrix_scatcnv,emis_vector_scatcnv)
      if (err /= 0) then
	  msg = 'error in scatcnv!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
      end if   
      rt_scattermatrix(i_z,:,:,:,:,:) = rt_scattermatrix(i_z,:,:,:,:,:) + scatter_matrix_scatcnv
      rt_extmatrix(i_z,:,:,:,:) = rt_extmatrix(i_z,:,:,:,:) + extinct_matrix_scatcnv
      rt_emisvec(i_z,:,:,:) = rt_emisvec(i_z,:,:,:) + emis_vector_scatcnv
    end if

    if (active .and. rt_hydros_present(i_z)) then
      call radar_simulator(err,radar_spec, rt_back(i_z), rt_kexttot(i_z),&
	atmo_delta_hgt_lev(i_x,i_y,i_z))
      if (err /= 0) then
	  msg = 'error in radar_simulator!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
      end if   
    end if


  end do grid_z !end of cycle over the vertical layers

  call deallocate_scatProperties()


  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
  return

end subroutine hydrometeor_extinction

