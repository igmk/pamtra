!+ soft spheroid diameter and density calculation
subroutine make_soft_spheroid(errorstatus)

! Description:
!  The routine calculate diameter and density (soft_d_eff, 
!  soft_rho_eff) of soft spheres for each of the size bins.
!  The relations used are:
!  m_i = a * D^b 
!  m_i = rho_i * volume_i
!  volume_i = pi/6 * rho_i * soft_d_eff^3 * (1/As_ratio)^2  for prolate (as_ratio > 1.)
!  volume_i = pi/6 * rho_i * soft_d_eff^3 * As_ratio        for oblate  (as_ratio < 1.)
!  From these relation we get soft_rho_eff, if rho_ms is set to -99, or 
!  soft_d_eff otherwise.
!  as_ratio = dimension_along_axis_of_rotation / dimension_perpendicular_to_axis_of_rotation
! 
! Owner: IGMK
!
! History:
!
! Version Date         Comment
! ------- ----         -------
! 0.      21/06/2013   Creation - E. Orlandi
!
! Code Description:
!  Language: Fortran 90.
!  Software Standards: "European Standards for Writing and
!   Documenting Exchangeable Fortran 90 Code".
!
! Parent Module: 
!
! Declarations:

! Modules used:

  use kinds, only: dbl, &      ! integer parameter specifying double precision
                    long       ! integer parameter specifying long integer

  use report_module

  use constants, only: pi, rho_ice

  use drop_size_dist, only: rho_ms, as_ratio, d_ds, nbin, mass_ds,  &    ! IN
     soft_rho_eff, soft_d_eff,liq_ice                                     ! OUT
  use settings, only: hydro_limit_density_area, hydro_softsphere_min_density, freqs
  use vars_index, only: i_f

  implicit none

!- End of header ---------------------------------------------------------------

! Local array:

  real(kind=dbl), dimension(nbin) :: mass         ! particle mass  [kg]

! Local scalar:

  integer(kind=long) :: i

! Error handling

  integer(kind=long), intent(out) :: errorstatus
!   integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=18) :: nameOfRoutine = 'make_soft_spheroid'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  allocate(soft_rho_eff(nbin)) ! Emil, why isn't this done in allocateVars_drop_size_dist ?
  allocate(soft_d_eff(nbin)) ! Emil, why isn't this done in allocateVars_drop_size_dist ?

  soft_rho_eff(:) = 0.
  soft_d_eff(:) = 0.
! Particle mass
  mass = mass_ds

! density set to a fixed value by the user
! Calculate the diameter of the soft spheroids
  if (rho_ms > 0.) then
    soft_rho_eff(:) = rho_ms
! oblate spheroid or sphere
    if (as_ratio <= 1.) then
      do i=1,nbin
        if (as_ratio < 0.) soft_d_eff(i) = ((6._dbl * mass(i)) / (pi *  rho_ms           ))**(1._dbl/3._dbl) ! spheres
        if (as_ratio > 0.) soft_d_eff(i) = ((6._dbl * mass(i)) / (pi *  rho_ms * as_ratio))**(1._dbl/3._dbl)
      enddo
    endif
! prolate spheroid
    if (as_ratio > 1.) then
      do i=1,nbin
        soft_d_eff(i) = ((6._dbl * mass(i) * as_ratio**2._dbl) / (pi * rho_ms))**(1._dbl/3._dbl)
      enddo
    endif
! Calculate the density of the soft spheroids
  elseif (rho_ms < 0.) then
    soft_d_eff = d_ds
    ! oblate spheroid or sphere
    if (as_ratio <= 1.) then
      do i=1,nbin
!         if (as_ratio < 0.) soft_rho_eff(i) = (6._dbl * mass(i)) / (pi *  d_ds(i)**3._dbl)
        if (as_ratio == -99.) soft_rho_eff(i) = (6._dbl * mass(i)) / (pi *  d_ds(i)**3._dbl) ! assumes spheres for the default value -99
        if (as_ratio == -10.) soft_rho_eff(i) = 917._dbl !CLOUD ICE 
        if (as_ratio == -20.) soft_rho_eff(i) = 0.863_dbl * freqs(i_f) + 115._dbl !SNOW
        if (as_ratio == -30.) soft_rho_eff(i) = 0.815_dbl * freqs(i_f) + 11.2_dbl !GRAUPEL
        if (as_ratio > 0.) soft_rho_eff(i) = (6._dbl * mass(i)) / (pi *  d_ds(i)**3._dbl * as_ratio) ! oblate spheroid
        if (soft_rho_eff(i) < hydro_softsphere_min_density) soft_rho_eff(i) = hydro_softsphere_min_density  ! used to be 5._dbl
        if (as_ratio == -10. .or. as_ratio == -20. .or. as_ratio == -30.) &
             soft_d_eff(i) = ((6._dbl * mass(i)) / (pi *  soft_rho_eff(i) ))**(1._dbl/3._dbl)
      enddo
    endif
    ! prolate spheroid
    if (as_ratio > 1.) then
      do i=1,nbin
        soft_rho_eff(i) = (6._dbl * mass(i) * as_ratio**2._dbl) / (pi *  d_ds(i)**3._dbl)
      enddo
    endif
  else
    msg = 'rho_ms = 0 in '
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  endif

  !   change mass and density in case density is larger than 917 or below hydro_softsphere_min_density
  if ((liq_ice == -1) .and. hydro_limit_density_area) then 
    do i=1,nbin
      if ((soft_rho_eff(i) > rho_ice) .or. (soft_rho_eff(i) < hydro_softsphere_min_density)) then
        if (soft_rho_eff(i) > rho_ice) then
          Write( msg, '("density too high:", f10.2)' )  soft_rho_eff(i)   
          if (verbose >= 1) call report(warning, msg, nameOfRoutine)
          soft_rho_eff(i) = rho_ice
        else
          Write( msg, '("density too low:", f10.2)' )  soft_rho_eff(i)   
          if (verbose >= 1) call report(warning, msg, nameOfRoutine)
          soft_rho_eff(i) = hydro_softsphere_min_density
        end if
        if (as_ratio <= 0.) then
          mass(i) =  (pi *  d_ds(i)**3._dbl) * soft_rho_eff(i) / 6._dbl
        else if ((as_ratio > 0.d0) .and. (as_ratio <= 1.d0)) then
          mass(i) =  (pi *  d_ds(i)**3._dbl * as_ratio) * soft_rho_eff(i) / 6._dbl
        else if (as_ratio > 1.d0) then
          mass(i) =  (pi *  d_ds(i)**3._dbl) * soft_rho_eff(i)  /  (6._dbl * as_ratio**2._dbl) 
        end if
      end if
    end do
  end if

  if ((minval(soft_rho_eff) <= 0.) .or. (minval(soft_d_eff) <= 0.)) then
    msg = 'something wrong in make_soft_spheroid!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  endif

  errorstatus = success
  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

  if (verbose >= 4) then
    print*,'soft spheroids diameter [m]'
    print*,soft_d_eff
    print*,'soft spheroids density [kg/m3]'
    print*,soft_rho_eff
  endif

  return

end subroutine make_soft_spheroid