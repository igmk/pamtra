module scatProperties

  use kinds
  use settings, only: nstokes,&
      nummu,&
      radar_nfft_aliased, &
      active, &
      radar_mode, &
      maxnleg
  use mie_spheres, only: calc_mie_spheres
  use tmatrix, only: calc_tmatrix
  use vars_atmosphere, only: kexttot,&
      back,&
      scattermatrix, &
      extmatrix, &
      emisvec
  use report_module
  use drop_size_dist, only:
        liq_ice,&
        nbin,&
        soft_d_eff, &
        delta_d_ds, &
        n_ds,&
        soft_rho_eff,&
        as_ratio
  use constants, only: pi, Im

  implicit none

    !rt3 style
    real(kind=dbl) :: salbedo
    real(kind=dbl), dimension(6,200) :: legen_coef
    integer :: nlegen_coef

    !rt4 style

    !needed by rt3 and rt4
    character(len=15) :: scat_name
    real(kind=dbl), allocatable, dimension(:) :: radar_spec

  contains
  
  subroutine allocate_scatProperties()
   implicit none
   allocate(radar_spec(radar_nfft_aliased)) 
   return
  end subroutine allocate_scatProperties

    subroutine prepare_rt3_scatProperties(iz)
      implicit none

      integer, intent(in) :: iz


      kexttot(iz) = 0.d0
      salbedo     = 0.d0
      back(iz) = 0.d0
      legen_coef(:,:)   = 0.d0
      nlegen_coef     = 0
      return
    end subroutine prepare_rt3_scatProperties

    subroutine prepare_rt4_scatProperties(iz)
      implicit none
      integer, intent(in) :: iz

      

      scattermatrix(iz,:,:,:,:,:)=0.d0
      extmatrix(iz,:,:,:,:)=0.d0
      emisvec(iz,:,:,:)=0.d0

      return

    end subroutine prepare_rt4_scatProperties

  subroutine calc_scatProperties(errorstatus, &
    freq,&
    iz)
    implicit none

    real(kind=dbl), intent(in) :: freq
    integer, intent(in) :: iz

    real(kind=dbl), dimension(200) :: legen_coef1_hydro 
    real(kind=dbl), dimension(200) :: legen_coef2_hydro
    real(kind=dbl), dimension(200) :: legen_coef3_hydro 
    real(kind=dbl), dimension(200) :: legen_coef4_hydro

    real(kind=dbl), dimension(nstokes,nummu,nstokes,nummu,4) :: scatter_matrix_hydro
    real(kind=dbl), dimension(nstokes,nstokes,nummu,2) :: extinct_matrix_hydro
    real(kind=dbl), dimension(nstokes,nummu,2) :: emis_vector_hydro
    real(kind=dbl), dimension(radar_nfft_aliased) :: radar_spec_hydro
    real(kind=dbl), dimension(nbin) :: back_spec_dia
    real(kind=dbl) :: kext_hydro
    real(kind=dbl) :: salb_hydro
    real(kind=dbl) :: back_hydro
    real(kind=dbl) :: refre
    real(kind=dbl) :: refim
    real(kind=dbl) :: absind
    real(kind=dbl) :: abscof
      
    complex(kind=dbl) :: refIndex  
      
    !tmp:
    real(kind=dbl) :: a_mice 
    real(kind=dbl) :: b_mice 
    real(kind=dbl) :: a_as_ice 
    real(kind=dbl) :: b_as_ice


    integer :: nlegen_coef_hydro
    integer :: jj


    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=40) :: nameOfRoutine = 'calc_scatProperties'


    !initilaize empyt results
    emis_vector_hydro(:,:,:)        = 0.d0
    extinct_matrix_hydro(:,:,:,:)   = 0.d0
    scatter_matrix_hydro(:,:,:,:,:) = 0.d0
    radar_spec_hydro(:) = 0.d0

    legen_coef1_hydro(:)   = 0.d0
    legen_coef2_hydro(:)  = 0.d0
    legen_coef3_hydro(:)  = 0.d0
    legen_coef4_hydro(:)  = 0.d0
    nlegen_coef_hydro     = 0


    !get the refractive index
     if (liq_ice == 1) then
        call ref_water(0.d0, t-273.15, freq, refre, refim, absind, abscof)
      else if (liq_ice == -1) then
        call ref_ice(t, freq, refre, refim)
      else
        errorstatus = fatal
        print*,"liq_ice=", liq_ice
        msg = 'Did not understand variable liq_ice'
        call report(errorstatus, msg, nameOfRoutine)
        return
      end if
      refIndex = refre-Im*refim  ! mimicking a


!!!!modern RT4 routines !!!

    if (scat_name == "TMatrix") then
    
      !some fixed settings for Tmatrix

  
      call calc_tmatrix(err,&
        freq*1.d9,&
        refIndex,&
        liq_ice,&
        nbin,&
        soft_d_eff, &
        delta_d_ds, &
        n_ds,&
        soft_rho_eff,&
        as_ratio,& 
        scatter_matrix_hydro(:,:,:,:,1:2),&
        extinct_matrix_hydro(:,:,:,1),&
        emis_vector_hydro(:,:,1),&
        back_spec_dia)

        
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

      back_hydro = scatter_matrix_hydro(1,16,1,16,2) !scatter_matrix(A,B;C;D;E) backscattering is M11 of Mueller or Scattering Matrix (A;C=1), in quadrature 2 (E) first 16 (B) is 180deg (upwelling), 2nd 16 (D) 0deg (downwelling). this definition is lokkiing from BELOW, scatter_matrix(1,16,1,16,3) would be from above!
      back_hydro = 4*pi*back_hydro!/k**2 !eq 4.82 Bohren&Huffman without k**2 (because of different definition of Mueller matrix according to Mishenko AO 2000). note that scatter_matrix contains already squard entries!
      kext_hydro = extinct_matrix_hydro(1,1,16,1) !11 of extinction matrix (=not polarized), at 0°, first quadrature. equal to extinct_matrix(1,1,16,2)

!!!!old style RT3 routines !!!

    else if (scat_name == "mie-sphere") then
      call calc_mie_spheres(err,&
            freq*1d9,&
            t,&
            liq_ice,&
            nbin,&
            soft_d_eff,&
            delta_d_ds, &
            n_ds,&
            soft_rho_eff, &
            refre, &
            refim, & !positive(?)
!OUT
            kext_hydro,&
            salb_hydro,&
            back_hydro,&
            nlegen_coef_hydro,&
            legen_coef1_hydro,&
            legen_coef2_hydro,&
            legen_coef3_hydro,&
            legen_coef4_hydro,&
            back_spec_dia)    
          
      nlegen_coef = max(nlegen_coef,nlegen_coef_hydro)
      if (err /= 0) then
          msg = 'error in calc_mie_spheres!'
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
    kexttot(iz) = kexttot(iz) + kext_hydro
    back(iz) = back(iz) + back_hydro

    !sum up rt4 style
    scattermatrix(iz,:,:,:,:,:) = scattermatrix(iz,:,:,:,:,:) + 
scatter_matrix_hydro
    extmatrix(iz,:,:,:,:) = extmatrix(iz,:,:,:,:) + extinct_matrix_hydro
    emisvec(iz,:,:,:) = emisvec(iz,:,:,:) + emis_vector_hydro


    !sum up rt3 style
    if (kexttot(iz) == 0.d0) then 
      salbedo = 0.d0
    else 
      salbedo = salbedo + (salb_hydro*kext_hydro)
    end if

    if (nlegen_coef > 0) then
      nlegen_loop: do jj = 1, nlegen_coef
        legen_coef(1,jj) = legen_coef(1,jj) + (legen_coef1_hydro(jj) * 
salb_hydro * kext_hydro)
        legen_coef(2,jj) = legen_coef(2,jj) + (legen_coef2_hydro(jj) * 
salb_hydro * kext_hydro)
        legen_coef(3,jj) = legen_coef(3,jj) + (legen_coef3_hydro(jj) * 
salb_hydro * kext_hydro)
        legen_coef(4,jj) = legen_coef(4,jj) + (legen_coef4_hydro(jj) * 
salb_hydro * kext_hydro)

      end do nlegen_loop
      legen_coef(5,:) = legen_coef(1,:)
      legen_coef(6,:) = legen_coef(3,:)
    end if





    !do checks
    if (kexttot(iz) .lt. 0. .or. isnan(kexttot(iz))) then
      print*, "kexttot(iz)",kexttot(iz)
      msg = 'kexttot(iz) smaller than zero or nan!'
      errorstatus = fatal
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if        

    if (back(iz) .lt. 0. .or. isnan(back(iz))) then
      print*, "back(iz)",back(iz)
      msg = 'back(iz) smaller than zero or nan!'
      errorstatus = fatal
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if      

    if (nlegen_coef .gt. maxnleg) then
      print*, "nlegen_coef",nlegen_coef
      msg = 'nlegen_coef greater maxnleg'
      errorstatus = fatal
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if   


    if ((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. 
"moments"))) then
     print*, "TODO transfer correct mass size and area size relation ro radar 
simulator"
     a_mice = 0.82d0
     b_mice = 2.5d0
     !area-size relation in SI
     a_as_ice = 0.12028493607054538 !0.24 in CGS
     b_as_ice = 1.85d0 !from mitchell 1996 similar to a_msnow&b_snow


      call radar_spectrum(nbin,d_bound_ds, back(iz),  
back_spec_dia,t,pressure,freq,&
        "ice",a_mice,b_mice,a_as_ice,b_as_ice,radar_spec_hydro)
      
      radar_spec(:) = radar_spec(:)+ radar_spec_hydro(:)
    end if
    
    


    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return 



  end subroutine calc_scatProperties

  subroutine finalize_rt3_scatProperties(iz)
    implicit none
    integer, intent(in) :: iz


    salbedo = salbedo / kexttot(iz)
    legen_coef(:,:) = legen_coef(:,:) / (salbedo * kexttot(iz))

    return
  end subroutine finalize_rt3_scatProperties

  subroutine deallocate_scatProperties()
   implicit none
   if (allocated(radar_spec)) deallocate(radar_spec)
    return
  end subroutine deallocate_scatProperties

 


end module 