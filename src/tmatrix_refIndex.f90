subroutine tmatrix_refIndex(freq,t,as_ratio_in,diameter,particle_mass,ptype,scatter_matrix,extinct_matrix,emis_vector)

  use kinds
  use report_module
  use constants, only: pi, c, im
  use settings, only : nummu, nstokes
  implicit none

  real(kind=dbl), intent(in) :: freq
  real(kind=dbl), intent(in) :: t
  real(kind=dbl), intent(in) :: as_ratio_in
  real(kind=dbl), intent(in) :: particle_mass
  real(kind=dbl), intent(in) :: diameter
  character(len=4), intent(in) :: ptype

  real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector

  real(kind=dbl) :: wavelength
  real(kind=dbl) :: wave_num
  real(kind=dbl) :: equiv_radius
  real(kind=dbl) :: density_eff
  real(kind=dbl) :: refre, refim
  real(kind=dbl) :: eu_alpha, eu_beta
  real(kind=dbl) :: bin_wgt
  real(kind=dbl) :: del_d
  real(kind=dbl) :: as_ratio

  integer :: azimuth_num, azimuth0_num

  complex(kind=ext) :: mindex
  complex(kind=dbl) :: msphere, eps_mix,m_ice

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'tmatrix_refIndex'

  interface

    subroutine tmatrix_calc(quad,qua_num,frequency,wave_num,ref_index,axi, nstokes,&
    as_ratio, alpha, beta, azimuth_num, azimuth0_num,&
    scatter_matrix,extinct_matrix,emis_vector)
      use kinds
      implicit none
      character(1), intent(in) :: quad
      integer, intent(in) :: qua_num
      real(kind=dbl), intent(in) :: frequency
      real(kind=dbl), intent(in) :: wave_num
      complex(kind=ext) :: ref_index
      real(kind=dbl), intent(in) :: axi
      integer, intent(in) :: nstokes
      real(kind=dbl), intent(in) :: as_ratio
      real(kind=dbl), intent(in) :: alpha
      real(kind=dbl), intent(in) :: beta
      integer, intent(in) :: azimuth_num
      integer, intent(in) :: azimuth0_num
      real(kind=dbl), intent(out), dimension(nstokes,qua_num,nstokes,qua_num,2) :: scatter_matrix
      real(kind=dbl), intent(out), dimension(nstokes,nstokes,qua_num) :: extinct_matrix
      real(kind=dbl), intent(out), dimension(nstokes,qua_num) :: emis_vector
    end subroutine tmatrix_calc
  end interface


  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  if (verbose >= 4) print*, "freq,t,as_ratio_in,diameter,particle_mass"
  if (verbose >= 4) print*, freq,t,as_ratio_in,diameter,particle_mass

  wavelength = c/freq !
  wave_num = 2.0_dbl*pi/wavelength
  eu_alpha = 0.0_dbl    ! orientation of the particle [°]
  eu_beta = 0.0_dbl!orientation of the particle [°]
  azimuth_num = 30
  azimuth0_num = 1


  if (ptype == "snow" .or. ptype == "ice") then
    call ref_ice(t, freq*1d-9, refre, refim)
    m_ice = refre-Im*refim  ! mimicking a

      !XINXINs code
  ! 	    CALL CAL_REFRACTIVE_INDEX('S',t,freq, diameter, as_ratio, particle_mass*1.d3, equiv_radius, mindex)

      !diameter of sphere with same mass
      equiv_radius = 0.5_dbl*diameter*as_ratio_in**(1.0_dbl/3.0_dbl)
      density_eff = (3.d0/4.d0 * particle_mass) / (pi * equiv_radius**3)


      if (density_eff > 917.d0) then
	if (verbose >= 0) print*, "WARNING changed density from ", density_eff, "kg/m3 to 917 kg/m3 for d=", diameter
	density_eff = 917.d0
      end if
      if (verbose >= 4) print*, "density_eff, equiv_radius, diameter, particle_mass,as_ratio, "
      if (verbose >= 4) print*, density_eff, equiv_radius, diameter, particle_mass, as_ratio_in

      msphere = eps_mix((1.d0,0.d0),m_ice,density_eff)
      mindex =conjg(msphere) !different convention
    else
      print*, nameOfRoutine, " do not know ptype:", ptype
      stop
    end if

    !make it numerically more stable
    if (as_ratio_in .eq. 1.d0) then
      as_ratio = 0.9999
    else
      as_ratio = as_ratio_in
    end if

    call tmatrix_calc('L',nummu,freq,wave_num,mindex,equiv_radius,nstokes,&
            as_ratio, eu_alpha, eu_beta, azimuth_num, azimuth0_num, &
            scatter_matrix,extinct_matrix,emis_vector)



  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return

end subroutine tmatrix_refIndex



