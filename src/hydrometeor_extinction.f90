subroutine hydrometeor_extinction(errorstatus)

  use kinds
  use vars_atmosphere, only: atmo_nlyrs, atmo_temp, atmo_q_hum,&
      atmo_press,&
      atmo_delta_hgt_lev, atmo_hydro_q, atmo_hydro_reff, atmo_hydro_n
  use vars_rt, only:rt_hydros_present
  use vars_hydroFullSpec, only: &
      hydrofs_delta_d_ds, &
      hydrofs_density2scat, &
      hydrofs_diameter2scat, &
      hydrofs_d_bound_ds, &
      hydrofs_f_ds, &
      hydrofs_mass_ds, &
      hydrofs_area_ds, &
      hydrofs_nbins
  use settings, only: &
      verbose, &
      hydro_threshold, &
      save_psd, &
      hydro_fullSpec, &
      passive, &
      active          
  use constants
  use descriptor_file
  use drop_size_dist
  use report_module
  use scatProperties
  use vars_output, only: out_psd_area, out_psd_d_bound, out_psd_f, out_psd_mass
  use vars_index, only: i_x,i_y, i_z, i_h

  implicit none

!   use mod_io_strings
!   use conversions
!   use tmat_snow_db
!   use tmat_rain_db
!   use vars_output, only: out_radar_spectra, out_radar_snr, out_radar_moments,&
!        out_radar_quality, radar_slopes, Ze, out_att_hydro !output of the radar simulator for jacobian mode
!         use report_module
! 

  real(kind=dbl) ::    scatter_matrix_scatcnv(nstokes,nummu,nstokes,nummu,4)
  real(kind=dbl) ::    extinct_matrix_scatcnv(nstokes,nstokes,nummu,2)
  real(kind=dbl) ::    emis_vector_scatcnv(nstokes,nummu,2)

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

    hydros: do i_h = 1,n_hydro


      !these ones are fore sure 1D
      hydro_name = hydro_name_arr(i_h)
      liq_ice    = liq_ice_arr(i_h)
      scat_name  = scat_name_arr(i_h)
      vel_size_mod  = vel_size_mod_arr(i_h)
      !shortcut if the full spectrum is provided
      if (hydro_fullSpec) then

        nbin = hydrofs_nbins
        call allocateVars_drop_size_dist()

        delta_d_ds(:) = hydrofs_delta_d_ds(i_x,i_y,i_z,i_h,:)
        density2scat(:) = hydrofs_density2scat(i_x,i_y,i_z,i_h,:)
        diameter2scat(:) = hydrofs_diameter2scat(i_x,i_y,i_z,i_h,:)
        d_bound_ds(:) = hydrofs_d_bound_ds(i_x,i_y,i_z,i_h,:)
        f_ds(:) = hydrofs_f_ds(i_x,i_y,i_z,i_h,:)
        mass_ds(:) = hydrofs_mass_ds(i_x,i_y,i_z,i_h,:)
        area_ds(:) = hydrofs_area_ds(i_x,i_y,i_z,i_h,:)

        pressure = atmo_press(i_x,i_y,i_z)
        layer_t = atmo_temp(i_x,i_y,i_z)
        rt_hydros_present(i_z) = .true.

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


        !short cut in case we disabled the particle
        if (dist_name == "disabled") then
          if (verbose >=3) print*, i_x,i_y,i_z,i_h,hydro_name, dist_name, "DISABLED"
          CYCLE
        end if

  ! Convert specific quantities [kg/kg] in absolute ones [kg/m3]
        q_h        = q2abs(atmo_hydro_q(i_x,i_y,i_z, i_h),atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z),&
                    atmo_q_hum(i_x,i_y,i_z),sum(atmo_hydro_q(i_x,i_y,i_z, :)))
        n_tot      = q2abs(atmo_hydro_n(i_x,i_y,i_z, ih),atmo_temp(i_x,i_y,i_z),atmo_press(i_x,i_y,i_z),&
                  atmo_q_hum(i_x,i_y,i_z),sum(atmo_hydro_q(i_x,i_y,i_z, :)))
        r_eff      = atmo_hydro_reff(i_x,i_y,i_z, i_h)
        layer_t    = atmo_temp(i_x,i_y,i_z)
        pressure   = atmo_press(i_x,i_y,i_z)

        if (verbose >= 2) print*, i_h, hydro_name

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

        if (save_psd) then
          out_psd_d_bound(i_x,i_y,i_z,i_h,1:nbin+1) =  d_bound_ds
          out_psd_f(i_x,i_y,i_z,i_h,1:nbin+1) = f_ds
          out_psd_mass(i_x,i_y,i_z,i_h,1:nbin+1) = mass_ds
          out_psd_area(i_x,i_y,i_z,i_h,1:nbin+1) = area_ds
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

