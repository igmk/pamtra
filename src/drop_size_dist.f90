!+ driver for the pamtra hydrometeor module
module drop_size_dist

! Modules used:

  use kinds, only: dbl, &      ! integer parameter specifying double precision
                   long        ! integer parameter specifying long integer

  use report_module

  use constants, only: rho_water
  use settings, only: hydro_limit_density_area

  use conversions, only: q2abs

  implicit none
  character(len=10)   :: hydro_name           ! hydrometeor name
  real(kind=dbl)      :: as_ratio             ! aspect ratio
  integer(kind=long)  :: liq_ice              ! liquid = 1; ice = -1

  integer(kind=long)  :: moment_in            ! Moments input via "PAMTRA input file" (1=q; 2=Ntot; 3=reff; 12=q&ntot; 13=q&reff; 23=ntot&reff)
  integer(kind=long)  :: nbin                 ! Number of bins for the drop-size distribution
  character(len=15)   :: dist_name            ! name of the distribution
  real(kind=dbl)      :: p_1, p_2, p_3, p_4   ! Drop-size parameters from hydrometeor descriptor file
  real(kind=dbl)      :: d_1, d_2             ! Minimum and maximum particle diameters

  real(kind=dbl)      :: q_h                  ! hydrometeor absolute concentration [kg/m3]
  real(kind=dbl)      :: n_tot                ! Total hydrometeor number concentration [#/m3]
  real(kind=dbl)      :: r_eff                ! Effective radius [m]

  real(kind=dbl)      :: layer_t              ! Layer temperature [K]
  real(kind=dbl)      :: pressure             ! Layer Pressure [Pa]
  real(kind=dbl)      :: dsd_canting          ! canting angle in deg of the particles

! make_mass_size IN/OUT
  real(kind=dbl)      :: rho_ms               ! density of the particle [kg/m^3]
  real(kind=dbl)      :: a_ms, b_ms           ! Mass-size parameter a [kg/m^(1/b)] and b [#] 

  real(kind=dbl)      :: alpha_as, beta_as    ! Area-size parameter alpha [m^(2-b)] and beta [#] 

! make_dist_params OUT
  real(kind=dbl)                 ::  n_0, lambda, gam, mu    ! parameter of the modified Gamma distribution
  real(kind=dbl)                 ::  n_t, sig, d_ln          ! parameter of the log-normal distribution
  real(kind=dbl)                 ::  d_mono                  ! diameter for the monodisperse distribution
  real(kind=dbl)                 ::  d_m, n_0_star           ! additional parameter for the normalized gamam distribution (with mu)

! make_dist OUT
  real(kind=dbl), dimension(:), allocatable  :: d_ds         ! particle diameter                        [m]
  real(kind=dbl), dimension(:), allocatable  :: d_bound_ds   ! boundaries of the particle dimeter bins  [m]
  real(kind=dbl), dimension(:), allocatable  :: delta_d_ds   ! width of the particle dimeter bins       [m]
  real(kind=dbl), dimension(:), allocatable  :: n_ds         ! particle number concentration            [#/m3]
  real(kind=dbl), dimension(:), allocatable  :: f_ds         ! drop-size distribution at the borders    [#/m4]

! calc_moment OUT
  real(kind=dbl)                 :: m_0                      ! 0th moment, total number                 [#/m3]
  real(kind=dbl)                 :: m_32                     ! 3nd/2rd moment, effective radius         [m]
  real(kind=dbl)                 :: am_b                     ! total mass, a_ms * bth moment            [kg/m3]

! make_soft_spheroid OUT
  real(kind=dbl), dimension(:), allocatable  :: soft_d_eff   ! particle diameter of soft spheroids      [m]
  real(kind=dbl), dimension(:), allocatable  :: soft_rho_eff ! particle density of soft spheroids       [kg/m^3]

! Particles density & diameter
  real(kind=dbl), dimension(:), allocatable  :: density2scat ! particle density for scattering routines[kg/m^3]
  real(kind=dbl), dimension(:), allocatable  :: diameter2scat! particle diameter for scattering routines [m]
! Particle Mass & area
  real(kind=dbl), dimension(:), allocatable  :: mass_ds ! particle mass for radar simulator [kg]
  real(kind=dbl), dimension(:), allocatable  :: area_ds ! particle cross section area for radar simulator [m^2]

! SSRGA parameter
  real(kind=dbl), dimension(:), allocatable  :: rg_kappa_ds ! kappa parameter
  real(kind=dbl), dimension(:), allocatable  :: rg_beta_ds ! kappa parameter
  real(kind=dbl), dimension(:), allocatable  :: rg_gamma_ds ! kappa parameter 
  real(kind=dbl), dimension(:), allocatable  :: rg_zeta_ds ! kappa parameter
  
 contains

subroutine allocateVars_drop_size_dist

  implicit none

  allocate(d_ds(nbin))
  allocate(n_ds(nbin))
  allocate(delta_d_ds(nbin))
  allocate(density2scat(nbin))
  allocate(diameter2scat(nbin))
  allocate(d_bound_ds(nbin+1))
  allocate(f_ds(nbin+1))
  allocate(mass_ds(nbin))
  allocate(area_ds(nbin))
  allocate(rg_kappa_ds(nbin))
  allocate(rg_beta_ds(nbin))
  allocate(rg_gamma_ds(nbin))
  allocate(rg_zeta_ds(nbin))

end subroutine allocateVars_drop_size_dist

subroutine deallocateVars_drop_size_dist

  implicit none

  if (allocated(d_ds)) deallocate(d_ds)
  if (allocated(n_ds)) deallocate(n_ds)
  if (allocated(density2scat)) deallocate(density2scat)
  if (allocated(diameter2scat)) deallocate(diameter2scat)
  if (allocated(d_bound_ds)) deallocate(d_bound_ds)
  if (allocated(delta_d_ds)) deallocate(delta_d_ds)  
  if (allocated(f_ds)) deallocate(f_ds)
  if (allocated(soft_d_eff)) deallocate(soft_d_eff)
  if (allocated(soft_rho_eff)) deallocate(soft_rho_eff)
  if (allocated(area_ds)) deallocate(area_ds)
  if (allocated(mass_ds)) deallocate(mass_ds)
  if (allocated(rg_kappa_ds)) deallocate(rg_kappa_ds)
  if (allocated(rg_beta_ds)) deallocate(rg_beta_ds)
  if (allocated(rg_gamma_ds)) deallocate(rg_gamma_ds)
  if (allocated(rg_zeta_ds)) deallocate(rg_zeta_ds)

end subroutine deallocateVars_drop_size_dist

subroutine run_drop_size_dist(errorstatus)

! Error handling

  integer(kind=long)  :: errorstatus, ibin
  integer(kind=long)  :: err
  character(len=80)   :: msg
  character(len=14)   :: nameOfRoutine = 'run_drop_size_dist'

  err = 0

  call make_mass_size(errorstatus)

  if (errorstatus == 2) then
    msg = 'Error in make_mass_size'
    call report(errorstatus, msg, nameOfRoutine)
    return
  end if

  call make_dist_params(errorstatus)

  if (verbose >= 4) then
    write(6,'(2(a15),8(a20))') 'hydro_name','dist_name','n_0','lambda','mu','gam','n_t','sig','d_ln','d_mono', 'd_m', 'n_0_star'
    write(6,'(2(a15),8(e20.10))') trim(hydro_name),trim(dist_name),n_0, lambda, mu, gam, n_t, sig, d_ln, d_mono, d_m, n_0_star
  endif

  if (errorstatus == 2) then
    msg = 'Error in make_dist_params'
    call report(errorstatus, msg, nameOfRoutine)
    return
  end if

  call make_dist(errorstatus)

  if (errorstatus == 2) then
    msg = 'Error in make_dist'
    call report(err, msg, nameOfRoutine)
    return
  end if

! Print out the distribution
!   if (verbose >= 4) then
!     do ibin=1,nbin+1
!       print*,'   distribution: bin boundaries[m], drop_size_dist[1/m4]',d_bound_ds(ibin),f_ds(ibin)
!     enddo
!   endif
  if (verbose >= 4) then
    do ibin=1,nbin
      print*,'   distribution: bin center[m], drop_size_dist[1/m4], % of total number, % of total mass',&
                 d_ds(ibin),n_ds(ibin),n_ds(ibin)/n_tot*100.,n_ds(ibin)*mass_ds(ibin)/q_h*100.
    enddo
  print*,'sum(n_ds)',sum(n_ds)
  endif


! Calculate particle MASS at bin boundaries
  do ibin=1,nbin
    mass_ds(ibin) = a_ms * d_ds(ibin)**b_ms
    !if mass density is larger than ice it is corrected, see drop_size_dist.f90
  enddo

! Calculate particle AREA at bin boundaries
  area_ds(:) = -99.
  if ((alpha_as > 0.) .and. (beta_as > 0.)) then
    do ibin=1,nbin
      area_ds(ibin) = alpha_as * d_ds(ibin)**beta_as
      !if area is larger than a square:
      if ((liq_ice == -1) .and. &
          (hydro_limit_density_area) .and. &
          (area_ds(ibin) > d_ds(ibin)**2)) then
          Write( msg, '("area too large:", e10.2, e10.2)' )  area_ds(ibin), d_ds(ibin)**2
        if (verbose >= 1)  call report(warning, msg, nameOfRoutine)
        area_ds(ibin) =  d_ds(ibin)**2
      end if
    enddo
  endif

  call calc_moment(errorstatus)

  if (errorstatus == 2) then
    msg = 'Error in calc_moment'
    call report(err, msg, nameOfRoutine)
    return
  end if

  if (liq_ice == -1) then
    call make_soft_spheroid(errorstatus)
  endif

! fill in the density and diameter array for the scattering routines
  if (liq_ice == -1) then  ! ice
    density2scat = soft_rho_eff
    diameter2scat = soft_d_eff
  endif
  if (liq_ice == 1) then  ! liquid
    density2scat(:) = rho_water
    diameter2scat = d_ds
  endif

  if (errorstatus == 2) then
    msg = 'Error in calc_moment'
    call report(err, msg, nameOfRoutine)
    return
  end if

  ! print results of dist if verbosity is high
  if (verbose >= 2) then
    call check_print(dist_name)
  end if 


end subroutine run_drop_size_dist

end module drop_size_dist