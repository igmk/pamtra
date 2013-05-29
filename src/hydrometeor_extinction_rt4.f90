subroutine hydrometeor_extinction_rt4(f,nx,ny,fi)

  use kinds
  use vars_atmosphere
  use settings, only: tmp_path, active, passive, dump_to_file, &
                       n_moments, quad_type, nummu, EM_snow, EM_grau, &
		       EM_hail, EM_ice, EM_rain, EM_cloud, as_ratio, &
                       use_rain_db, use_snow_db, data_path, &
		       jacobian_mode, radar_nfft_aliased, radar_mode
  use constants
  use mod_io_strings
  use conversions
  use tmat_snow_db
  use tmat_rain_db
  use vars_output, only: radar_spectra, radar_snr, radar_moments,&
	radar_quality, radar_slope, Ze, Att_hydro !output of the radar simulator for jacobian mode
        use report_module

  implicit none

  integer, intent(in) :: nx,ny,fi
  integer, parameter :: maxleg = 200
  integer, parameter :: nstokes = 2

  integer :: jj, nz

  integer :: nlegencw, nlegenci, nlegenrr, nlegensn, nlegengr, nlegenha

  real(kind=dbl) :: f, qwc, nc, cwc

  real(kind=dbl) :: salbcw, salbrr, salbci, salbsn, salbgr, salbha


  real(kind=dbl), dimension(200) ::  LEGENcw, LEGENrr, LEGENci, LEGENgr, LEGENsn, LEGENha,      &
       LEGEN2cw, LEGEN2rr, LEGEN2ci, LEGEN2gr, LEGEN2sn, LEGEN2ha,  &
       LEGEN3cw, LEGEN3rr, LEGEN3ci, LEGEN3gr, LEGEN3sn, LEGEN3ha,  &
       LEGEN4cw, LEGEN4rr, LEGEN4ci, LEGEN4gr, LEGEN4sn, LEGEN4ha


!   real(kind=dbl), dimension(2) :: P11, ang

  real(kind=dbl) :: threshold ! threshold value for hydrometeor extinction as mass mixing ratio

  real(kind=dbl), dimension(6,100) :: coef



  real(kind=dbl), dimension(nstokes,nummu,nstokes,nummu,4) :: rain_scat,snow_scat,&
      ice_scat,graupel_scat,hail_scat,cloud_scat
  real(kind=dbl), dimension(nstokes,nstokes,nummu,2) :: rain_ext,snow_ext,&
      ice_ext,graupel_ext,hail_ext,cloud_ext
  real(kind=dbl), dimension(nstokes,nummu,2) :: rain_emis,snow_emis,&
      ice_emis,graupel_emis,hail_emis,cloud_emis


  real(kind=dbl), dimension(radar_nfft_aliased) :: cloud_spec, rain_spec, snow_spec,&
      ice_spec, graupel_spec, hail_spec, full_spec

  CHARACTER*64 SCATFILES(nlyr)


  integer, parameter :: nquad = 16

  logical :: didNotChange 

    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'hydrometeor_extinction_rt4'  
  
  if (verbose .gt. 1) print*, nx, ny, 'Entering hydrometeor_extinction_rt4'

  ! INITIALIZATION OF LEGENDRE COEFFICIENTS

 

  nlegen = 0
  legen   = 0.d0
  legen2  = 0.d0
  legen3  = 0.d0
  legen4  = 0.d0
  coef = 0.d0

  hydros_present = .false.

  threshold = 1.d-10   ! [kg/kg]

  if (use_rain_db) then
      rdb_file = data_path(:len_trim(data_path))//'/tmatrix/tmatrix_rain.dat'
      call initialize_rain_db
  end if
  if (use_snow_db) then
      write(as_str,'(f3.1)') as_ratio
      sdb_file = data_path(:len_trim(data_path))//'/tmatrix/tmatrix_s_'//as_str//'.dat'
      call initialize_snow_db
  end if

  if (verbose .gt. 1) print*, 'start loop over layer'

  grid_z: do nz = 1, nlyr  ! loop over all layers
! print *,temp(nz)
     if (verbose .gt. 1) print*, 'Layer: ', nz

  !jacobian mode take profile 1,1 as a reference, all other are compared to this one

  if (jacobian_mode .and. ((nx .ne. 1) .or. (ny .ne. 1))) then
    !check whether profil is the same as in the one of ny=1, nx=1
    !make boolean
    didNotChange= ( &
      (jac_temp_lev(nz) .eq. temp_lev(nz)) .and. &
      (jac_relhum_lev(nz) .eq. relhum_lev(nz)) .and. &
      (jac_temp_lev(nz-1) .eq. temp_lev(nz-1)) .and. &
      (jac_relhum_lev(nz-1) .eq. relhum_lev(nz-1)) .and. &
      (jac_cwc_q(nz) .eq. cwc_q(nz)) .and. &
      (jac_iwc_q(nz) .eq. iwc_q(nz)) .and. &
      (jac_rwc_q(nz) .eq. rwc_q(nz)) .and. &
      (jac_swc_q(nz) .eq. swc_q(nz)) .and. &
      (jac_gwc_q(nz) .eq. gwc_q(nz)))

    if (n_moments .eq. 2) then
      didNotChange = (didNotChange .and. &
	(jac_hwc_q(nz) .eq. hwc_q(nz)) .and. &
	(jac_cwc_n(nz) .eq. cwc_n(nz)) .and. &
	(jac_iwc_n(nz) .eq. iwc_n(nz)) .and. &
	(jac_rwc_n(nz) .eq. rwc_n(nz)) .and. &
	(jac_swc_n(nz) .eq. swc_n(nz)) .and. &
	(jac_gwc_n(nz) .eq. gwc_n(nz)) .and. &
	(jac_hwc_n(nz) .eq. hwc_n(nz)))
      end if

    if (verbose .gt. 1) print*,"jacobian_mode:",nx,ny,nz,didNotChange
    !if layer is identical, then reference jac_xx is used
    if (didNotChange) then
      scattermatrix(nz,:,:,:,:,:)=jac_scattermatrix(nz,:,:,:,:,:)
      extmatrix(nz,:,:,:,:)=jac_extmatrix(nz,:,:,:,:)
      emisvec(nz,:,:,:)=jac_emisvec(nz,:,:,:)
      hydros_present(nz)=jac_hydros_present(nz)
            
      kextsn(nz) = jac_kextsn(nz)
      backsn(nz) = jac_backsn(nz)
      kextcw(nz) = jac_kextcw(nz)
      backcw(nz) = jac_backcw(nz)
      kextrr(nz) = jac_kextrr(nz)
      backrr(nz) = jac_backrr(nz)
      kextgr(nz) = jac_kextgr(nz)
      backgr(nz) = jac_backgr(nz)
      kextci(nz) = jac_kextci(nz)
      backci(nz) = jac_backci(nz)
      kextha(nz) = jac_kextha(nz)
      backha(nz) = jac_backha(nz)

      kexttot(nz) = kextsn(nz) + kextcw(nz) + kextrr(nz) + kextgr(nz) + kextci(nz) + kextha(nz)
      back(nz) = backcw(nz) + backrr(nz) + backci(nz) + backsn(nz) + backgr(nz) + backha(nz)

  !short cut for the radar simulator: the ouput is taken directly from the final arrays at position 1,1.
  if (active) then
    Ze(nx,ny,nz,fi) = Ze(1,1,nz,fi)
    Att_hydro(nx,ny,nz,fi) = Att_hydro(1,1,nz,fi)

    if ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments")) then
      radar_spectra(nx,ny,nz,fi,:) = radar_spectra(1,1,nz,fi,:)
      radar_snr(nx,ny,nz,fi) =   radar_snr(1,1,nz,fi)
      radar_moments(nx,ny,nz,fi,:) = radar_moments(1,1,nz,fi,:)
      radar_slope(nx,ny,nz,fi,:) =   radar_slope(1,1,nz,fi,:)
      radar_quality(nx,ny,nz,fi) =   radar_quality(1,1,nz,fi)
    end if
  end if

      CYCLE

      end if
    end if



     !---------------------------salinity------------------------------
     ! calculation of the single scattering properties
     ! of hydrometeors. cloud water and cloud ice are 
     ! with respect to radius. whereas the distribution 
     ! of precipitating particles is with respect to diameter.
     !---------------------------------------------------------

     !---------------------------------------------------------
     !        single scattering properties of cloud water
     !---------------------------------------------------------

     nlegencw = 0 
     legencw  = 0.d0
     legen2cw = 0.d0
     legen3cw = 0.d0
     legen4cw = 0.d0
     cloud_scat = 0.d0
     cloud_emis= 0.d0
     cloud_ext = 0.d0
     kextcw(nz) = 0.d0 
     salbcw = 0.d0 
     backcw(nz) = 0.d0 
     cloud_spec(:) = 0.d0
     if ((cwc_q(nz) .ge. threshold) .and. (EM_cloud .ne. 'disab')) then
        hydros_present(nz) = .true.
     	if (n_moments .eq. 1) then
	     	qwc = q2abs(cwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
		nc = 0.d0
     	else if (n_moments .eq. 2) then
	     	qwc = q2abs(cwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	        nc =  q2abs(cwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
       	end if

    	call cloud_ssp(f,qwc,temp(nz),press(nz),&
             	maxleg, nc, kextcw(nz), salbcw, backcw(nz),  &
             	nlegencw, legencw, legen2cw, legen3cw, legen4cw, &
		cloud_scat, cloud_ext, cloud_emis,cloud_spec)
     end if


     !---------------------------------------------------------
     !       single scattering properties of ice crystals
     !---------------------------------------------------------

     nlegenci = 0 
     legenci  = 0.d0
     legen2ci = 0.d0
     legen3ci = 0.d0
     legen4ci = 0.d0
     ice_scat = 0.d0
     ice_emis= 0.d0
     ice_ext = 0.d0
     kextci(nz) = 0.0d0 
     salbci = 0.0d0 
     backci(nz) = 0.0d0 
     ice_spec(:) = 0.d0
     if ((iwc_q(nz) .ge. threshold) .and. (EM_ice .ne. 'disab')) then
        hydros_present(nz) = .true.
	     if (n_moments .eq. 1) then
	     	qwc = q2abs(iwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
		nc = 0.d0
	     else if (n_moments .eq. 2) then
	     	qwc = q2abs(iwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	        nc = q2abs(iwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	     end if
	     call ice_ssp(f,qwc,temp(nz),press(nz),&
	             maxleg,nc,kextci(nz), salbci, backci(nz),  &
	             nlegenci, legenci, legen2ci, legen3ci, legen4ci, &
			ice_scat, ice_ext, ice_emis,ice_spec)

     end if

     !---------------------------------------------------------
     !       single scattering properties of rain
     !---------------------------------------------------------

     nlegenrr = 0
     legenrr  = 0.d0
     legen2rr = 0.d0
     legen3rr = 0.d0
     legen4rr = 0.d0
     rain_scat = 0.0d0
     rain_ext = 0.0d0
     rain_emis = 0.0d0

     kextrr(nz) = 0.d0
     salbrr = 0.d0
     backrr(nz) = 0.d0
     rain_spec(:) = 0.d0

     if ((rwc_q(nz) .ge. threshold) .and. EM_rain .ne. "disab") then
        hydros_present(nz) = .true.
     	if (n_moments .eq. 1) then
	       qwc = q2abs(rwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
	       cwc = q2abs(cwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
	       nc = 0.d0
	else if (n_moments .eq. 2) then
            qwc = q2abs(rwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
            cwc = q2abs(cwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	    nc = q2abs(rwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
         end if
    	     	call rain_ssp(f,qwc,cwc,temp(nz),press(nz),&
	                maxleg,nc,kextrr(nz), salbrr, backrr(nz),  &
	                nlegenrr, legenrr, legen2rr, legen3rr, legen4rr, &
			rain_scat, rain_ext, rain_emis,rain_spec)

     end if

     !---------------------------------------------------------
     !       single scattering properties of snow
     !---------------------------------------------------------

     nlegensn = 0 
     legensn = 0.0d0
     legen2sn = 0.0d0
     legen3sn = 0.0d0
     legen4sn = 0.0d0
     snow_scat = 0.0d0
     snow_ext = 0.0d0
     snow_emis = 0.0d0
    kextsn(nz) = 0.0d0
    salbsn = 0.0d0
    backsn(nz) = 0.0d0
     snow_spec(:) = 0.d0
     if ((swc_q(nz) .ge. threshold) .and. (EM_snow .ne. 'disab')) then
        hydros_present(nz) = .true.
     	 if (n_moments .eq. 1) then
	     	qwc = q2abs(swc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
		nc = 0.d0
	     else if (n_moments .eq. 2) then
	     	qwc = q2abs(swc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	        nc = q2abs(swc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	  end if
    	  call snow_ssp(f,qwc,temp(nz),press(nz),&
	                 maxleg,nc,kextsn(nz), salbsn, backsn(nz),  &
	                 nlegensn, legensn, legen2sn, legen3sn, legen4sn,&
			 snow_scat, snow_ext, snow_emis,snow_spec)

     endif

     !---------------------------------------------------------
     !       single scattering properties of graupel
     !---------------------------------------------------------

     nlegengr = 0 
     legengr = 0.0d0
     legen2gr = 0.0d0
     legen3gr = 0.0d0
     legen4gr = 0.0d0
     graupel_scat = 0.0d0
     graupel_ext = 0.0d0
     graupel_emis = 0.0d0
        kextgr(nz) = 0.0d0
        salbgr = 0.0d0
        backgr(nz) = 0.0d0
     graupel_spec(:) = 0.d0
     if ((gwc_q(nz) .ge. threshold) .and. (EM_grau .ne. 'disab'))then
        hydros_present(nz) = .true.
	     if (n_moments .eq. 1) then
	     	qwc = q2abs(gwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
		nc = 0.d0
	     else if (n_moments .eq. 2) then
	     	qwc = q2abs(gwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	        nc = q2abs(gwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	     end if
	      	call grau_ssp(f,qwc,temp(nz),press(nz),&
	             maxleg,nc, kextgr(nz), salbgr, backgr(nz),  &
	             nlegengr, legengr, legen2gr, legen3gr, legen4gr,&
		     graupel_scat, graupel_ext, graupel_emis,graupel_spec)

     endif

     !---------------------------------------------------------
     !       single scattering properties of hail
     !---------------------------------------------------------

     nlegenha = 0
     legenha = 0.0d0
     legen2ha = 0.0d0
     legen3ha = 0.0d0
     legen4ha = 0.0d0
     hail_scat = 0.0d0
     hail_ext = 0.0d0
     hail_emis = 0.0d0

      kextha(nz) = 0.0d0
      salbha = 0.0d0
      backha(nz) = 0.0d0
     hail_spec(:) = 0.d0
     if (n_moments .eq. 2) then
        if ((hwc_q(nz) .ge. threshold) .and. (EM_hail .ne. 'disab')) then
           hydros_present(nz) = .true.
	       qwc = q2abs(hwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	       nc = q2abs(hwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
           call hail_ssp(f,qwc,temp(nz),press(nz),&
                maxleg,nc,kextha(nz), salbha, backha(nz),  &
                nlegenha, legenha, legen2ha, legen3ha, legen4ha,&
		     hail_scat, hail_ext, hail_emis,hail_spec)
        endif
     endif

     nlegen(nz) = max(nlegen(nz),nlegencw,nlegenci,nlegenrr,nlegensn,nlegengr,nlegenha)

     if (verbose .gt. 1) print*, 'End of scattering calc for layer: ', nz

     !CCCCCCCCCCCCC   END OF SINGLE SCATTERING PROPERTY COMPUTATIONS  CCCCCCC

     !                                                                       
     !           Summing up the scattering parameters and writing the
     !           input file of the scattering properties of each layer
     !                                                                       


     kexttot(nz) = kextsn(nz) + kextcw(nz) + kextrr(nz) + kextgr(nz) + kextci(nz) + kextha(nz)
     back(nz) = backcw(nz) + backrr(nz) + backci(nz) + backsn(nz) + backgr(nz) + backha(nz)

     if (kexttot(nz) .lt. 0.) write(*,*) 'something wrong'
     if (kexttot(nz) .le. 0.) then 
        salbtot(nz) = 0.d0
     else 
        salbtot(nz) = (salbcw * kextcw(nz) + salbrr *       &
             kextrr(nz) + salbci * kextci(nz) + salbsn * kextsn(nz) + salbgr *    &
             kextgr(nz) + salbha * kextha(nz)) / kexttot(nz)
     endif
     !
     !      absorp(nz) = (1.0 - salbtot(nz) ) * kexttot(nz)

!!!!!!!!!!!!!!!!! check whether hgt_lev needs to be km or m !!!!!!!!!!!!!!!!!

     !   summing up the Legendre coefficient                                 

!     if (kexttot(nz) .gt. 0.0 .or. salbtot(nz) .gt. 0.0) then ! there are hydrometeors present
     if (hydros_present(nz)) then ! there are hydrometeors present
      if (nlegen(nz) .gt. 0) then

        do jj = 1, Nlegen(nz)
           legen(nz,jj) = (legencw (jj) * salbcw * kextcw(nz) + legenrr ( &
                jj) * salbrr * kextrr(nz) + legenci (jj) * salbci * kextci(nz) + &
                legensn (jj) * salbsn * kextsn(nz) + legengr (jj) * salbgr * &
                kextgr(nz) + legenha (jj) * salbha * kextha(nz)) / (salbtot (nz) &
                * kexttot (nz) )

           legen2(nz,jj) = (legen2cw (jj) * salbcw * kextcw(nz) +         &
                legen2rr (jj) * salbrr * kextrr(nz) + legen2ci (jj) * salbci &
                * kextci(nz) + legen2sn (jj) * salbsn * kextsn(nz) + legen2gr (  &
                jj) * salbgr * kextgr(nz) + legen2ha (jj) * salbha * kextha(nz) ) &
                / (salbtot (nz) * kexttot (nz) )

           legen3(nz,jj) = (legen3cw (jj) * salbcw * kextcw(nz) +         &
                legen3rr (jj) * salbrr * kextrr(nz) + legen3ci (jj) * salbci &
                * kextci(nz) + legen3sn (jj) * salbsn * kextsn(nz) + legen3gr (  &
                jj) * salbgr * kextgr(nz) + legen3ha (jj) * salbha * kextha(nz)) &
                / (salbtot (nz) * kexttot (nz) )

           legen4(nz,jj) = (legen4cw(jj) * salbcw * kextcw(nz) +         &
                legen4rr(jj) * salbrr * kextrr(nz) + legen4ci(jj) * salbci &
                * kextci(nz) + legen4sn(jj) * salbsn * kextsn(nz) + legen4gr (  &
                jj) * salbgr * kextgr(nz) + legen4ha(jj) * salbha * kextha(nz)) &
                / (salbtot(nz) * kexttot(nz))
           g_coeff(nz) = legen (nz,2) / 3.0d0
           coef(1,jj) = legen(nz,jj)
           coef(2,jj) = legen2(nz,jj)
           coef(3,jj) = legen3(nz,jj)
           coef(4,jj) = legen4(nz,jj)
           coef(5,jj) = legen(nz,jj)
           coef(6,jj) = legen3(nz,jj)
           if (dump_to_file) then
              write (22, 1005) jj - 1, legen (nz,jj), legen2 (nz,jj),        &
                   legen3(nz,jj), legen4(nz,jj), legen(nz,jj), legen3(nz,jj)
           end if
	    end do ! end of cycle over Legendre coefficient



     call scatcnv(scatfiles(nz),nlegen(nz),coef,kexttot(nz),salbtot(nz),&
     scattermatrix(nz,:,:,:,:,:),extmatrix(nz,:,:,:,:),emisvec(nz,:,:,:))

    else
     scattermatrix(nz,:,:,:,:,:)=0.d0
     extmatrix(nz,:,:,:,:)=0.d0
     emisvec(nz,:,:,:)=0.d0


    end if

     scattermatrix(nz,:,:,:,:,:)=scattermatrix(nz,:,:,:,:,:)+&
	rain_scat+snow_scat+ice_scat+graupel_scat+hail_scat+cloud_scat
     extmatrix(nz,:,:,:,:)=extmatrix(nz,:,:,:,:)+&
	rain_ext+snow_ext+ice_ext+graupel_ext+hail_ext+cloud_ext
     emisvec(nz,:,:,:)=emisvec(nz,:,:,:)+&
	rain_emis+snow_emis+ice_emis+graupel_emis+hail_emis+cloud_emis

     full_spec = cloud_spec+rain_spec+snow_spec+ &
      ice_spec+graupel_spec+hail_spec



  if (active) then
     call radar_simulator(err,full_spec, back(nz), kexttot(nz), f,&
      temp(nz),delta_hgt_lev(nz),nz,nx,ny,fi)
    if (err /= 0) then
	msg = 'error in radar_simulator!'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	stop !return
    end if   
  end if

 end if !end if hydrometeors present

  end do grid_z !end of cycle over the vertical layers



  if ((nx .eq. 1 ) .and. (ny .eq. 1 ) .and. jacobian_mode) then
    !for jacobian mode safe results of 1,1 grid
    jac_scattermatrix=scattermatrix
    jac_extmatrix=extmatrix
    jac_emisvec=emisvec
    jac_hydros_present = hydros_present
    
    jac_kextsn = kextsn
    jac_backsn = backsn
    jac_kextcw = kextcw
    jac_backcw = backcw
    jac_kextrr = kextrr
    jac_backrr = backrr
    jac_kextgr = kextgr
    jac_backgr = backgr
    jac_kextci = kextci
    jac_backci = backci
    jac_kextha = kextha
    jac_backha = backha

    !to see what was changed we need also the profile
    jac_temp_lev=temp_lev
    jac_relhum_lev=relhum_lev
    jac_cwc_q=cwc_q
    jac_iwc_q=iwc_q
    jac_rwc_q=rwc_q
    jac_swc_q=swc_q
    jac_gwc_q=gwc_q
    if (n_moments .eq. 2) then
      jac_hwc_q=hwc_q
      jac_cwc_n=cwc_n
      jac_iwc_n=iwc_n
      jac_rwc_n=rwc_n
      jac_swc_n=swc_n
      jac_gwc_n=gwc_n
      jac_hwc_n=hwc_n
    end if
  end if



  if (use_snow_db) call close_snow_db()
  if (use_rain_db) call close_rain_db()

  if (verbose .gt. 1) print*, 'Exiting hydrometeor_extinction_rt4'
1005 format  (i3,6(1x,f10.7))
  return

end subroutine hydrometeor_extinction_rt4
