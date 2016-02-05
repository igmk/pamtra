subroutine hydrometeor_extinction(errorstatus)

  use kinds
  use vars_atmosphere, only: atmo_nlyrs, atmo_temp, atmo_q_hum,&
      atmo_press,&
      atmo_delta_hgt_lev, atmo_hydro_q, atmo_hydro_reff, atmo_hydro_n
  use vars_rt, only: &
      rt_hydros_present, &
      rt_kexttot,&
      rt_back, &
      rt_scattermatrix, &
      rt_extmatrix, &
      rt_emisvec
  use vars_hydroFullSpec, only: &
      hydrofs_rho_ds, &
      hydrofs_d_ds, &
      hydrofs_d_bound_ds, &
      hydrofs_n_ds, &
      hydrofs_mass_ds, &
      hydrofs_area_ds, &
      hydrofs_nbins
  use settings, only: &
      verbose, &
      hydro_threshold, &
      save_psd, &
      save_ssp, &
      hydro_fullSpec, &
      passive, &
      active, &
      nstokes, &
      nummu, &
      radar_attenuation, &
      hydro_includeHydroInRhoAir, &
      lhyd_scattering, &
      lhyd_emission, &
      lhyd_absorption
  use constants
  use descriptor_file
  use drop_size_dist
  use report_module
  use scatProperties, only: &
    radar_spec, &
    vel_size_mod, &
    scat_name, &
    salbedo, &
    legen_coef, &
    nlegen_coef, &
    allocate_scatproperties, &
    deallocate_scatproperties, &
    prepare_rt3_scatproperties, &
    prepare_rt4_scatproperties, &
    calc_scatproperties, &
    finalize_rt3_scatproperties
  use vars_output, only: &
    out_psd_area, &
    out_psd_d, &
    out_psd_n, &
    out_psd_mass, &
    out_scatter_matrix, &
    out_extinct_matrix, &
    out_emis_vector

  use vars_index, only: i_x,i_y, i_z, i_h

  implicit none

  real(kind=dbl) ::    scatter_matrix_scatcnv(nstokes,nummu,nstokes,nummu,4)
  real(kind=dbl) ::    extinct_matrix_scatcnv(nstokes,nstokes,nummu,2)
  real(kind=dbl) ::    emis_vector_scatcnv(nstokes,nummu,2)

  integer(kind=long) :: increment, start, stop

  
  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=40) :: nameOfRoutine = 'hydrometeor_extinction'
  
  interface
    subroutine radar_simulator(errorstatus,particle_spectrum,back,kexthydro,&
    delta_h)
        use kinds
        use settings, only: radar_npol, radar_nfft_aliased
        implicit none
        real(kind=dbl), dimension(radar_npol),intent(in) ::  back
        real(kind=dbl),intent(in) ::  delta_h,kexthydro
        real(kind=dbl), dimension(radar_npol,radar_nfft_aliased),intent(in):: particle_spectrum
        integer(kind=long), intent(out) :: errorstatus
      end subroutine
  end interface

  if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)



  call allocate_scatProperties()

  !in case we want to calczulate the attenuation top-down we have to reverse the order
  if (TRIM(radar_attenuation) == "top-down") then
    increment = -1
    start = atmo_nlyrs(i_x,i_y) 
    stop =  1
  else
    increment = 1
    start = 1
    stop = atmo_nlyrs(i_x,i_y)
  end if

  grid_z: do i_z = start, stop,increment  ! loop over all layers

    call prepare_rt3_scatProperties()
    call prepare_rt4_scatProperties()
    rt_hydros_present(i_z) = .false.
    
    if (verbose .gt. 1) print*, 'Layer: ', i_z

    hydros: do i_h = 1,n_hydro

      !these ones are fore sure 1D
      hydro_name = hydro_name_arr(i_h)
      liq_ice    = liq_ice_arr(i_h)
      scat_name  = TRIM(scat_name_arr(i_h))
      vel_size_mod  = TRIM(vel_size_mod_arr(i_h))
      !shortcut if the full spectrum is provided
      if (hydro_fullSpec) then

        nbin = hydrofs_nbins
        call allocateVars_drop_size_dist()

        density2scat(:) = hydrofs_rho_ds(i_x,i_y,i_z,i_h,:)
        diameter2scat(:) = hydrofs_d_ds(i_x,i_y,i_z,i_h,:)
        d_ds(:) = hydrofs_d_ds(i_x,i_y,i_z,i_h,:)
        d_bound_ds(:) = hydrofs_d_bound_ds(i_x,i_y,i_z,i_h,:)
        n_ds(:) = hydrofs_n_ds(i_x,i_y,i_z,i_h,:)
        mass_ds(:) = hydrofs_mass_ds(i_x,i_y,i_z,i_h,:)
        area_ds(:) = hydrofs_area_ds(i_x,i_y,i_z,i_h,:)
        delta_d_ds(:) = d_bound_ds(2:) - d_bound_ds(:-1)


        pressure = atmo_press(i_x,i_y,i_z)
        layer_t = atmo_temp(i_x,i_y,i_z)
        
        if (PRODUCT(SHAPE(canting_arr)) == n_hydro) then
          dsd_canting   = canting_arr(1,1,1,i_h)
        else
          dsd_canting   = canting_arr(i_x,i_y,i_z,i_h)
        end if         

        if (SUM(n_ds) > 0.d0) then
          rt_hydros_present(i_z) = .true.
        else
          call deallocateVars_drop_size_dist()
          ! not neccesary to set rt_hydros_present to false - its default value is false!
          CYCLE
        end if

      else
        moment_in  = moment_in_arr(i_h)
        dist_name  = dist_name_arr(i_h)

        ! fill 0-D variable for the run_drop_size routine
        ! start with the ones which might come from a 4D python field.
        if (PRODUCT(SHAPE(as_ratio_arr)) == n_hydro) then
          as_ratio   = as_ratio_arr(1,1,1,i_h)
        else
          as_ratio   = as_ratio_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(rho_ms_arr)) == n_hydro) then
          rho_ms   = rho_ms_arr(1,1,1,i_h)
        else
          rho_ms   = rho_ms_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(a_ms_arr)) == n_hydro) then
          a_ms   = a_ms_arr(1,1,1,i_h)
        else
          a_ms   = a_ms_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(b_ms_arr)) == n_hydro) then
          b_ms   = b_ms_arr(1,1,1,i_h)
        else
          b_ms   = b_ms_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(alpha_as_arr)) == n_hydro) then
          alpha_as   = alpha_as_arr(1,1,1,i_h)
        else
          alpha_as   = alpha_as_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(beta_as_arr)) == n_hydro) then
          beta_as   = beta_as_arr(1,1,1,i_h)
        else
          beta_as   = beta_as_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(nbin_arr)) == n_hydro) then
          nbin   = nbin_arr(1,1,1,i_h)
        else
          nbin   = nbin_arr(i_x,i_y,i_z,i_h)
        end if 

!   Force nbin = 2 when monodisperse distribution is used. Needed by radar simulator.
        if ( (trim(dist_name) == 'mono') .or. (trim(dist_name) == 'mono_cosmo_ice')) &
          nbin = 2


        if (PRODUCT(SHAPE(p_1_arr)) == n_hydro) then
          p_1   = p_1_arr(1,1,1,i_h)
        else
          p_1   = p_1_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(p_2_arr)) == n_hydro) then
          p_2   = p_2_arr(1,1,1,i_h)
        else
          p_2   = p_2_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(p_3_arr)) == n_hydro) then
          p_3   = p_3_arr(1,1,1,i_h)
        else
          p_3   = p_3_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(p_4_arr)) == n_hydro) then
          p_4   = p_4_arr(1,1,1,i_h)
        else
          p_4   = p_4_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(d_1_arr)) == n_hydro) then
          d_1   = d_1_arr(1,1,1,i_h)
        else
          d_1   = d_1_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(d_2_arr)) == n_hydro) then
          d_2   = d_2_arr(1,1,1,i_h)
        else
          d_2   = d_2_arr(i_x,i_y,i_z,i_h)
        end if 

        if (PRODUCT(SHAPE(canting_arr)) == n_hydro) then
          dsd_canting   = canting_arr(1,1,1,i_h)
        else
          dsd_canting   = canting_arr(i_x,i_y,i_z,i_h)
        end if 

        !short cut in case we disabled the particle
        if (dist_name == "disabled") then
          if (verbose >=3) print*, i_x,i_y,i_z,i_h,hydro_name, dist_name, "DISABLED"
          CYCLE
        end if

  ! Convert specific quantities [kg/kg] in absolute ones [kg/m3]
	if (hydro_includeHydroInRhoAir .and. .not.(isnan(sum(atmo_hydro_q(i_x,i_y,i_z, :))))) then
	  q_h        = q2abs(atmo_hydro_q(i_x,i_y,i_z, i_h),atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z),&
		      atmo_q_hum(i_x,i_y,i_z),sum(atmo_hydro_q(i_x,i_y,i_z, :)))
	  n_tot      = q2abs(atmo_hydro_n(i_x,i_y,i_z, i_h),atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z),&
		    atmo_q_hum(i_x,i_y,i_z),sum(atmo_hydro_q(i_x,i_y,i_z, :)))
	else
	  q_h        = q2abs(atmo_hydro_q(i_x,i_y,i_z, i_h),atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z),&
		    atmo_q_hum(i_x,i_y,i_z),0._dbl)
	  n_tot      = q2abs(atmo_hydro_n(i_x,i_y,i_z, i_h),atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z),&
		    atmo_q_hum(i_x,i_y,i_z),0._dbl)
	end if
        r_eff      = atmo_hydro_reff(i_x,i_y,i_z, i_h)
        layer_t    = atmo_temp(i_x,i_y,i_z)
        pressure   = atmo_press(i_x,i_y,i_z)


! In case moment_in = 13 (number_concentration (#/kg) and mass_mixing ratio (kg/kg))
! if number_concentration is less than 1 particle per kg, cycle.
        if (moment_in .eq. 13 .and. (q_h < hydro_threshold .or. n_tot <= 1.)) then
          if (verbose >=3) print*, i_x,i_y,i_z,i_h,hydro_name, "q_h or n_tot below threshold", q_h, n_tot
          CYCLE
        endif
        if (verbose >= 2) print*, i_h, hydro_name
        if (verbose >= 3) print*, "i_h, q_h, n_tot, r_eff", i_h, q_h, n_tot, r_eff
        if ((q_h < hydro_threshold .or. isnan(q_h)) .and. &
            (n_tot <=0 .or. isnan(n_tot)) .and. &
            (r_eff <=0 .or. isnan(r_eff)) .and. &
            (moment_in /= 0)) then
          if (verbose >=3) print*, i_x,i_y,i_z,i_h,hydro_name, "q_h below threshold", q_h
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

        if (all(n_ds <= 0)) then
          if (verbose >=3) print*, i_x,i_y,i_z,i_h,hydro_name, "all n_ds < 0", MAXVAL(n_ds)
          call deallocateVars_drop_size_dist()
          CYCLE
        end if

        if (save_psd) then
          out_psd_d(i_x,i_y,i_z,i_h,1:nbin) =  d_ds
          out_psd_n(i_x,i_y,i_z,i_h,1:nbin) = n_ds
          out_psd_mass(i_x,i_y,i_z,i_h,1:nbin) = mass_ds
          out_psd_area(i_x,i_y,i_z,i_h,1:nbin) = area_ds
        end if

      end if !hydro_fullSpec

      if (active .or. passive) then
        call calc_scatProperties(err)
        if (err == 2) then
          msg = 'Error in calc_scatProperties'
          call report(err, msg, nameOfRoutine)
          errorstatus = err
          return
        end if
      end if

      call deallocateVars_drop_size_dist()
    end do hydros

    call finalize_rt3_scatProperties()

    !convert rt3 to rt4 input
    if (nlegen_coef>0) then
      call scatcnv(err,nlegen_coef,legen_coef,rt_kexttot(i_z),salbedo,&
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

    if (.not. lhyd_scattering) rt_scattermatrix(i_z,:,:,:,:,:) = 0.d0
    if (.not. lhyd_emission) rt_emisvec(i_z,:,:,:) = 0.d0
    
    
    if (active .and. rt_hydros_present(i_z)) then
      call radar_simulator(err,radar_spec, rt_back(i_z,:), rt_kexttot(i_z),&
	atmo_delta_hgt_lev(i_x,i_y,i_z))
      if (err /= 0) then
	  msg = 'error in radar_simulator!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  return
      end if   
    else
      if (verbose >= 5) print*, "Skipped radar simulator because of active .and. rt_hydros_present(i_z)", &
            active, rt_hydros_present(i_z)

    end if


  end do grid_z !end of cycle over the vertical layers

  call deallocate_scatProperties()


  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
  return

end subroutine hydrometeor_extinction

