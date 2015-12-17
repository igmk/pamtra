module scatProperties


  use kinds
  use settings, only: maxnleg
  implicit none

  !rt3 style
  real(kind=dbl) :: salbedo
  real(kind=dbl), dimension(6,maxnleg) :: legen_coef
  integer :: nlegen_coef

  !rt4 style

  !needed by rt3 and rt4
  character(len=30) :: scat_name
  character(len=30) :: vel_size_mod
  real(kind=dbl), allocatable, dimension(:,:) :: radar_spec

contains

  subroutine allocate_scatProperties()
    use kinds
    use settings, only: radar_nfft_aliased, radar_npol
    implicit none
    allocate(radar_spec(radar_npol,radar_nfft_aliased)) 
    return
  end subroutine allocate_scatProperties

  subroutine prepare_rt3_scatProperties()
    use kinds
    implicit none



    salbedo     = 0.d0
    legen_coef(:,:)   = 0.d0
    nlegen_coef     = 0
    return
  end subroutine prepare_rt3_scatProperties

  subroutine prepare_rt4_scatProperties()
    use kinds
    use vars_rt, only: rt_kexttot,&
         rt_back,&
         rt_scattermatrix, &
         rt_extmatrix, &
         rt_emisvec
    use vars_index, only: i_z

    implicit none

    rt_kexttot(i_z) = 0.d0  
    rt_back(i_z,:) = 0.d0

    rt_scattermatrix(i_z,:,:,:,:,:)=0.d0
    rt_extmatrix(i_z,:,:,:,:)=0.d0
    rt_emisvec(i_z,:,:,:)=0.d0

    radar_spec(:,:) = 0.d0
    return

  end subroutine prepare_rt4_scatProperties

  subroutine calc_scatProperties(errorstatus)

    use kinds
    use settings, only: nstokes,&
         nummu,&
         radar_nfft_aliased, &
         active, &
         radar_mode, &
         maxnleg, &
         freqs, &
         hydro_fullSpec, &
         radar_pol, &
         radar_npol
    use mie_spheres, only: calc_mie_spheres
    use tmatrix, only: calc_tmatrix
    use rayleigh_gans, only: calc_self_similar_rayleigh_gans, &
         calc_rayleigh_gans, &
         calc_rayleigh
    use vars_rt, only: rt_kexttot,&
         rt_back,&
         rt_scattermatrix, &
         rt_extmatrix, &
         rt_emisvec, &
         rt_hydros_present
    use report_module
    use drop_size_dist, only: liq_ice,&
         nbin,&
         diameter2scat, &
         delta_d_ds, &
         n_ds,&
         density2scat,&
         as_ratio, &
         d_bound_ds, &
         pressure, &
         layer_t, &
         mass_ds, &
         area_ds, &
         dsd_canting, &
         soft_rho_eff, &
         d_ds
    use constants, only: pi, Im
    use vars_index, only: i_x, i_y, i_z, i_f, i_h, i_p
    use vars_hydroFullSpec, only: hydrofs_as_ratio, hydrofs_canting

    implicit none

    real(kind=dbl) :: freq

    real(kind=dbl), dimension(maxnleg) :: legen_coef1_hydro 
    real(kind=dbl), dimension(maxnleg) :: legen_coef2_hydro
    real(kind=dbl), dimension(maxnleg) :: legen_coef3_hydro 
    real(kind=dbl), dimension(maxnleg) :: legen_coef4_hydro

    real(kind=dbl), dimension(nstokes,nummu,nstokes,nummu,4) :: scatter_matrix_hydro
    real(kind=dbl), dimension(nstokes,nstokes,nummu,2) :: extinct_matrix_hydro
    real(kind=dbl), dimension(nstokes,nummu,2) :: emis_vector_hydro
    real(kind=dbl), dimension(radar_nfft_aliased) :: radar_spec_hydro
    real(kind=dbl), dimension(nbin) :: num_density
    real(kind=dbl), dimension(radar_npol) :: back_hydro
    real(kind=dbl), dimension(radar_npol,radar_nfft_aliased) :: back_spec_dia
    real(kind=dbl), dimension(radar_nfft_aliased) :: back_spec_mie, back_spec_liu
    real(kind=dbl), dimension(radar_nfft_aliased) :: back_spec_rg
    real(kind=dbl), allocatable, dimension(:) :: as_ratio_list, canting_list
    real(kind=dbl) :: kext_hydro
    real(kind=dbl) :: salb_hydro
    real(kind=dbl) :: back_hydro_mie, back_hydro_liu
    real(kind=dbl) :: back_hydro_rg
    real(kind=dbl) :: refre
    real(kind=dbl) :: refim
    real(kind=dbl) :: absind
    real(kind=dbl) :: abscof
    real(kind=dbl) :: rg_kappa
    real(kind=dbl) :: rg_beta
    real(kind=dbl) :: rg_gamma

    character(30) :: tokenized(3)

    integer(kind=long) :: pos1, pos2, nn 

    complex(kind=dbl) :: refIndex  

    integer :: nlegen_coef_hydro
    integer :: jj
    integer(kind=long) :: liu_type

    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=40) :: nameOfRoutine = 'calc_scatProperties'

    interface
       subroutine radar_spectrum(&
            errorstatus, &
            nbins,&             !in
            diameter_spec,&     !in
            back,&              !in
            back_spec,&         !in
            temp,&              !in
            press,&             !in
            frequency,&         !in
            rho_particle,&      !in
            vel_size_mod,&      !in
            mass,&              !in
            area,&              !in
            particle_spec)      !out

         use kinds
         use settings, only: radar_nfft_aliased
         implicit none

         integer,intent(in) ::  nbins 

         real(kind=dbl), dimension(nbins),intent(in):: diameter_spec, back_spec
         real(kind=dbl), dimension(nbins),intent(in):: mass, area, rho_particle
         character(len=30),intent(in) :: vel_size_mod
         real(kind=dbl), intent(in):: temp, frequency, press,back
         real(kind=dbl), intent(out), dimension(radar_nfft_aliased):: particle_spec
         integer(kind=long), intent(out) :: errorstatus
       end subroutine radar_spectrum
    end interface

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    err = 0

    if ((scat_name == "disabled") .or. (.not. rt_hydros_present(i_z))) then
       if (verbose >= 3) print*, "OK, we are done here"
       return
    end if
    freq = freqs(i_f)

    ! initialize empty results
    emis_vector_hydro(:,:,:)        = 0.d0
    extinct_matrix_hydro(:,:,:,:)   = 0.d0
    scatter_matrix_hydro(:,:,:,:,:) = 0.d0
    radar_spec_hydro(:) = 0.d0

    legen_coef1_hydro(:)   = 0.d0
    legen_coef2_hydro(:)  = 0.d0
    legen_coef3_hydro(:)  = 0.d0
    legen_coef4_hydro(:)  = 0.d0
    nlegen_coef_hydro     = 0

    ! normalize particle density
    num_density = n_ds / delta_d_ds


    !get the refractive index
    if (liq_ice == 1) then
       call ref_water(err,0.d0, layer_t-273.15, freq, refre, refim, absind, abscof)
       if (err > 0) then
          errorstatus = fatal
          msg = 'Error in ref_water'
          call report(errorstatus, msg, nameOfRoutine)
          return
       end if
    else if (liq_ice == -1) then
       call ref_ice(layer_t, freq, refre, refim)
    else
       errorstatus = fatal
       print*,"liq_ice=", liq_ice
       msg = 'Did not understand variable liq_ice'
       call report(errorstatus, msg, nameOfRoutine)
       return
    end if
    refIndex = refre-Im*refim  ! mimicking a


    !some fixed settings for Tmatrix and rg
    allocate(as_ratio_list(nbin))
    allocate(canting_list(nbin))
    if (hydro_fullSpec) then
       as_ratio_list(:) = hydrofs_as_ratio(i_x,i_y,i_z,i_h,:)
       canting_list(:) = hydrofs_canting(i_x,i_y,i_z,i_h,:)
    else
       as_ratio_list(:) =  as_ratio
       canting_list(:) =  dsd_canting
    end if

    where (canting_list < 0) canting_list = 0.d0
    where (isnan(canting_list)) canting_list = 0.d0
    where (as_ratio_list < 0) as_ratio_list = 0.d0
    where (isnan(as_ratio_list)) as_ratio_list = 0.d0
!!!!modern RT4 routines !!!
    if (TRIM(scat_name) == "tmatrix") then

       call calc_tmatrix(err,&
            freq*1.d9,&
            refIndex,&
            liq_ice,&
            nbin,&
            diameter2scat, &
            delta_d_ds, &
            num_density,&
            density2scat,&
            as_ratio_list,& 
            canting_list, &
            layer_t, &
            scatter_matrix_hydro(:,:,:,:,1:2),&
            extinct_matrix_hydro(:,:,:,1),&
            emis_vector_hydro(:,:,1),&
            back_spec_dia)

       if (allocated(as_ratio_list)) deallocate(as_ratio_list)
       if (allocated(canting_list)) deallocate(canting_list)

       if (err /= 0) then
          msg = 'error in calc_tmatrix!'
          call report(err, msg, nameOfRoutine)
          errorstatus = err
          return
       end if

       !fill up the matrices
       scatter_matrix_hydro(:,:,:,:,4) = scatter_matrix_hydro(:,:,:,:,1) 
       scatter_matrix_hydro(:,:,:,:,3) = scatter_matrix_hydro(:,:,:,:,2)
       extinct_matrix_hydro(:,:,:,2) = extinct_matrix_hydro(:,:,:,1)
       emis_vector_hydro(:,:,2) = emis_vector_hydro(:,:,1)

       do i_p= 1, radar_npol
          if (radar_pol(i_p) == "NN") then
             back_hydro(i_p) = scatter_matrix_hydro(1,16,1,16,2) !scatter_matrix(A,B;C;D;E) backscattering is M11 of Mueller or Scattering Matrix (A;C=1), in quadrature 2 (E) first 16 (B) is 180deg (upwelling), 2nd 16 (D) 0deg (downwelling). this definition is looking from BELOW, scatter_matrix(1,16,1,16,3) would be from above!
          else if (radar_pol(i_p) == "HH") then
             !1.Vivekanandan, J., Adams, W. M. & Bringi, V. N. Rigorous Approach to Polarimetric Radar Modeling of Hydrometeor Orientation Distributions. Journal of Applied Meteorology 30, 1053–1063 (1991).
             back_hydro(i_p) = + scatter_matrix_hydro(1,16,1,16,2) &
                  - scatter_matrix_hydro(1,16,2,16,2) & 
                  - scatter_matrix_hydro(2,16,1,16,2) & 
                  + scatter_matrix_hydro(2,16,2,16,2) 
          else if (radar_pol(i_p) == "VV") then
             !1.Vivekanandan, J., Adams, W. M. & Bringi, V. N. Rigorous Approach to Polarimetric Radar Modeling of Hydrometeor Orientation Distributions. Journal of Applied Meteorology 30, 1053–1063 (1991).
             back_hydro(i_p) = + scatter_matrix_hydro(1,16,1,16,2) &
                  + scatter_matrix_hydro(1,16,2,16,2) & 
                  + scatter_matrix_hydro(2,16,1,16,2) & 
                  + scatter_matrix_hydro(2,16,2,16,2) 
          else if (radar_pol(i_p) == "HV") then
             !1.Vivekanandan, J., Adams, W. M. & Bringi, V. N. Rigorous Approach to Polarimetric Radar Modeling of Hydrometeor Orientation Distributions. Journal of Applied Meteorology 30, 1053–1063 (1991).
             back_hydro(i_p) = + scatter_matrix_hydro(1,16,1,16,2) &
                  - scatter_matrix_hydro(1,16,2,16,2) & 
                  + scatter_matrix_hydro(2,16,1,16,2) & 
                  - scatter_matrix_hydro(2,16,2,16,2) 
          else
             errorstatus = fatal
             msg = 'do not understand radar_pol(i_p): '//radar_pol(i_p)
             call report(errorstatus, msg, nameOfRoutine)
             return
          end if
       end do
       if (verbose >= 5) then
          print*, "S11",scatter_matrix_hydro(1,16,1,16,2)
          print*, "S12",scatter_matrix_hydro(1,16,2,16,2) 
          print*, "S21",scatter_matrix_hydro(2,16,1,16,2) 
          print*, "S22",scatter_matrix_hydro(2,16,2,16,2) 
       end if

       back_hydro(:) = 4*pi*back_hydro(:)!/k**2 !eq 4.82 Bohren&Huffman without k**2 (because of different definition of Mueller matrix according to Mishenko AO 2000). note that scatter_matrix contains already squard entries!
       kext_hydro = extinct_matrix_hydro(1,1,16,1) !11 of extinction matrix (=not polarized), at 0°, first quadrature. equal to extinct_matrix(1,1,16,2)

       !       back_hydro(1) = scatter_matrix_hydro(1,16,1,16,2) !scatter_matrix(A,B;C;D;E) backscattering is M11 of Mueller or Scattering Matrix (A;C=1), in quadrature 2 (E) first 16 (B) is 180deg (upwelling), 2nd 16 (D) 0deg (downwelling). this definition is lokkiing from BELOW, scatter_matrix(1,16,1,16,3) would be from above!
       !       !hh
       !       back_hydro(2) = (scatter_matrix_hydro(1,16,1,16,2) &
       !            - scatter_matrix_hydro(1,16,2,16,2) &
       !            - scatter_matrix_hydro(2,16,1,16,2) &
       !            + scatter_matrix_hydro(2,16,2,16,2) ) /2.d0
       !       !vv
       !       back_hydro(3) = (scatter_matrix_hydro(1,16,1,16,2) &
       !            + scatter_matrix_hydro(1,16,2,16,2) &
       !            + scatter_matrix_hydro(2,16,1,16,2) &
       !            + scatter_matrix_hydro(2,16,2,16,2) ) /2.d0
       !       !vh 
       !       back_hydro(4) = (scatter_matrix_hydro(1,16,1,16,2) &
       !            - scatter_matrix_hydro(1,16,2,16,2) &
       !            + scatter_matrix_hydro(2,16,1,16,2) &
       !            - scatter_matrix_hydro(2,16,2,16,2) ) /2.d0

       ! rayleigh gans only for active!
    else if (scat_name(:16) == "ss-rayleigh-gans") then


       if (len(TRIM(scat_name)) > 16) then
          pos1 = 1
          nn = 0
          DO
             pos2 = INDEX(scat_name(pos1:), "_")
             IF (pos2 == 0) THEN
                nn = nn + 1
                tokenized(nn) = scat_name(pos1:)
                EXIT
             END IF
             nn = nn + 1
             tokenized(nn) = scat_name(pos1:pos1+pos2-2)
             pos1 = pos2+pos1
          END DO
          read(tokenized(2),*) rg_kappa
          read(tokenized(3),*) rg_beta
       else
          !take default values for aggregates of bullet rosettes or columns from Hogan and Westbrook 2014
          rg_kappa = 0.19d0
          rg_beta = 0.23d0
       end if
       rg_gamma = 5.d0/3.d0 

       if (verbose >= 5) print*,scat_name,  rg_gamma, rg_kappa, rg_beta

       call calc_self_similar_rayleigh_gans(err,&
            freq*1d9,&
            liq_ice, &
            nbin,&
            diameter2scat,&
            delta_d_ds, &
            num_density,&
            mass_ds, &
            as_ratio_list, &
            canting_list, &
            refre, &
            refim, & !positive(?)
            rg_kappa, &
            rg_beta, &
            rg_gamma, &
            !OUT
            back_spec_rg, &
            back_hydro_rg )

       kext_hydro = 0.d0
       salb_hydro = 0.d0

       if (allocated(as_ratio_list)) deallocate(as_ratio_list)
       if (allocated(canting_list)) deallocate(canting_list)

       ! for mie polarisatuion does not matter
       do i_p=1, radar_npol
          back_spec_dia(i_p,:) = back_spec_rg(:)
          back_hydro(i_p) = back_hydro_rg
       end do

       if (err /= 0) then
          msg = 'error in calc_self_similar_rayleigh_gans!'
          call report(err, msg, nameOfRoutine)
          errorstatus = err
          return
       end if

       ! rayleigh gans only for active!
    else if (TRIM(scat_name) == "rayleigh-gans") then

       call calc_rayleigh_gans(err,&
            freq*1d9,&
            liq_ice, &
            nbin,&
            diameter2scat,&
            delta_d_ds, &
            num_density,&
            density2scat, &
            as_ratio_list, &
            canting_list, &
            refre, &
            refim, & !positive(?)
            !OUT
            back_spec_rg, &
            back_hydro_rg )

       kext_hydro = 0.d0
       salb_hydro = 0.d0

       if (allocated(as_ratio_list)) deallocate(as_ratio_list)
       if (allocated(canting_list)) deallocate(canting_list)

       ! for mie polarisatuion does not matter
       do i_p=1, radar_npol
          back_spec_dia(i_p,:) = back_spec_rg(:)
          back_hydro(i_p) = back_hydro_rg
       end do

       if (err /= 0) then
          msg = 'error in calc_rayleigh_gans!'
          call report(err, msg, nameOfRoutine)
          errorstatus = err
          return
       end if

       ! pure rayleigh d**6 only for active!
    else if (TRIM(scat_name) == "rayleigh") then

       call calc_rayleigh(err,&
            freq*1d9,&
            liq_ice, &
            nbin,&
            diameter2scat,&
            delta_d_ds, &
            num_density,&
            density2scat, &
            refre, &
            refim, & !positive(?)
            !OUT
            back_spec_rg, &
            back_hydro_rg )

       kext_hydro = 0.d0
       salb_hydro = 0.d0

       if (allocated(as_ratio_list)) deallocate(as_ratio_list)
       if (allocated(canting_list)) deallocate(canting_list)

       ! for mie polarisatuion does not matter
       do i_p=1, radar_npol
          back_spec_dia(i_p,:) = back_spec_rg(:)
          back_hydro(i_p) = back_hydro_rg
       end do

       if (err /= 0) then
          msg = 'error in calc_rayleigh!'
          call report(err, msg, nameOfRoutine)
          errorstatus = err
          return
       end if


!!!!old style RT3 routines !!!

    else if (TRIM(scat_name) == "mie-sphere") then

       call calc_mie_spheres(err,&
            freq*1d9,&
            layer_t,&
            liq_ice,&
            nbin,&
            diameter2scat,&
            delta_d_ds, &
            num_density,&
            density2scat, &
            refre, &
            refim, & !positive(?)
            !OUT
            kext_hydro,&
            salb_hydro,&
            back_hydro_mie,&
            nlegen_coef_hydro,&
            legen_coef1_hydro,&
            legen_coef2_hydro,&
            legen_coef3_hydro,&
            legen_coef4_hydro,&
            back_spec_mie)    

       ! for mie polarisation does not matter
       do i_p=1, radar_npol
          back_spec_dia(i_p,:) = back_spec_mie(:)
          back_hydro(i_p) = back_hydro_mie
       end do
       nlegen_coef = max(nlegen_coef,nlegen_coef_hydro)
       if (err /= 0) then
          msg = 'error in calc_mie_spheres!'
          call report(err, msg, nameOfRoutine)
          errorstatus = err
          return
       end if
    else if (scat_name(:5) == "liudb") then
      read(scat_name(7:8),*) liu_type
      call dda_db_liu(err,freq, layer_t, liu_type, refre-Im*refim,nbin,&
            diameter2scat,&
            delta_d_ds, &
            num_density,&
      kext_hydro, salb_hydro, back_hydro_liu,  &
     nlegen_coef_hydro, legen_coef1_hydro, legen_coef2_hydro, legen_coef3_hydro, &
     legen_coef4_hydro, back_spec_liu)
       ! for random oriented particles in Liu DB polarisation does not matter
       do i_p=1, radar_npol
          back_spec_dia(i_p,:) = back_spec_liu(:)
          back_hydro(i_p) = back_hydro_liu
       end do
       nlegen_coef = max(nlegen_coef,nlegen_coef_hydro)
      if (err /= 0) then
        msg = 'error in dda_db_liu!'
        call report(err, msg, nameOfRoutine)
        errorstatus = err
        return
      end if
    else
       msg = 'do not understand scat_name: '//scat_name
       errorstatus = fatal
       call report(errorstatus, msg, nameOfRoutine)
       return
    end if

    !sum up
    rt_kexttot(i_z) = rt_kexttot(i_z) + kext_hydro
    rt_back(i_z,:) = rt_back(i_z,:) + back_hydro(:)

    !sum up rt4 style
    rt_scattermatrix(i_z,:,:,:,:,:) = rt_scattermatrix(i_z,:,:,:,:,:) + scatter_matrix_hydro
    rt_extmatrix(i_z,:,:,:,:) = rt_extmatrix(i_z,:,:,:,:) + extinct_matrix_hydro
    rt_emisvec(i_z,:,:,:) = rt_emisvec(i_z,:,:,:) + emis_vector_hydro

    !sum up rt3 style
    if (rt_kexttot(i_z) == 0.d0) then 
       salbedo = 0.d0
    else 
       salbedo = salbedo + (salb_hydro*kext_hydro)
    end if

    if (nlegen_coef > 0) then
       nlegen_loop: do jj = 1, nlegen_coef
          legen_coef(1,jj) = legen_coef(1,jj) + (legen_coef1_hydro(jj) * salb_hydro * kext_hydro)
          legen_coef(2,jj) = legen_coef(2,jj) + (legen_coef2_hydro(jj) * salb_hydro * kext_hydro)
          legen_coef(3,jj) = legen_coef(3,jj) + (legen_coef3_hydro(jj) * salb_hydro * kext_hydro)
          legen_coef(4,jj) = legen_coef(4,jj) + (legen_coef4_hydro(jj) * salb_hydro * kext_hydro)

       end do nlegen_loop
       legen_coef(5,:) = legen_coef(1,:)
       legen_coef(6,:) = legen_coef(3,:)
    end if

    !do checks
    if ((rt_kexttot(i_z) .lt. 0.d0) .or. isnan(rt_kexttot(i_z))) then
       print*, "rt_kexttot(i_z)",rt_kexttot(i_z)
       msg = 'rt_kexttot(i_z) smaller than zero or nan!'
       errorstatus = fatal
       call report(errorstatus, msg, nameOfRoutine)
       return
    end if

    if (ANY(rt_back(i_z,:) .lt. 0.d0) .or. ANY(isnan(rt_back(i_z,:)))) then
       print*, "rt_back(i_z)",rt_back(i_z,:)
       msg = 'rt_back(i_z,:) smaller than zero or nan!'
       errorstatus = fatal
       call report(errorstatus, msg, nameOfRoutine)
       return
    end if

    if (nlegen_coef .gt. maxnleg-1) then
       print*, "nlegen_coef",nlegen_coef
       msg = 'nlegen_coef greater maxnleg'
       errorstatus = fatal
       call report(errorstatus, msg, nameOfRoutine)
       return
    end if

    if ((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments"))) then

       do  i_p= 1, radar_npol
          call radar_spectrum(err,nbin,d_ds, back_hydro(i_p),  back_spec_dia(i_p,:),layer_t,pressure,freq,&
               soft_rho_eff,vel_size_mod,mass_ds,area_ds,radar_spec_hydro)
          if (err /= 0) then
             msg = 'error in radar_spectrum!'
             call report(err, msg, nameOfRoutine)
             errorstatus = err
             return
          end if
          radar_spec(i_p,:) = radar_spec(i_p,:)+ radar_spec_hydro(:)
       end do
    end if




    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return 



  end subroutine calc_scatProperties

  subroutine finalize_rt3_scatProperties()
    use vars_index, only: i_x, i_y, i_z, i_f, i_h
    use vars_rt, only: rt_kexttot
    implicit none


    salbedo = salbedo / rt_kexttot(i_z)
    legen_coef(:,:) = legen_coef(:,:) / (salbedo * rt_kexttot(i_z))

    return
  end subroutine finalize_rt3_scatProperties

  subroutine deallocate_scatProperties()
    implicit none
    if (allocated(radar_spec)) deallocate(radar_spec)
    return
  end subroutine deallocate_scatProperties




end module scatProperties
