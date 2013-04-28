subroutine tmatrix_sql(freq,t,as_ratio,diameter,particle_mass,ptype,par_name,scatter_matrix,extinct_matrix,emis_vector)

  use kinds
  use report_module
  use constants, only: pi, c, im
  use settings, only : nummu, nstokes
  use report_module
  use sql_tools

  implicit none

  real(kind=dbl), intent(in) :: freq !in Hz
  real(kind=dbl), intent(in) :: t
  real(kind=dbl), intent(in) :: as_ratio
  real(kind=dbl), intent(in) :: particle_mass
  real(kind=dbl), intent(in) :: diameter
  character(len=4), intent(in) :: ptype
  character(len=20), intent(in) :: par_name


  real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector

  real(kind=dbl) :: freq_round
  real(kind=dbl) :: t_round
  real(kind=dbl) :: as_ratio_round
  real(kind=dbl) :: particle_mass_round
  real(kind=dbl) :: diameter_round
  real(kind=dbl) :: particle_mass_magnitude
  real(kind=dbl) :: diameter_magnitude


  logical :: found

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'tmatrix_sql'

  interface
    subroutine tmatrix_refIndex(freq,t,as_ratio,diameter,particle_mass,ptype,scatter_matrix,extinct_matrix,emis_vector)
      use kinds
      use settings, only : nummu, nstokes
      implicit none

      real(kind=dbl), intent(in) :: freq
      real(kind=dbl), intent(in) :: t
      real(kind=dbl), intent(in) :: as_ratio
      real(kind=dbl), intent(in) :: particle_mass
      character(len=4), intent(in) :: ptype
      real(kind=dbl), intent(in) :: diameter

      real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
      real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
      real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector
    end subroutine tmatrix_refIndex
  end interface

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)


  freq_round = NINT(freq*1.d-8)*1.d8
  t_round = NINT(t/1.d0)*1.d0
  as_ratio_round = NINT(as_ratio*20.d0)/20.d0

  diameter_magnitude = 10.d0**(INT(log10(diameter)-3))! change 3 to 4 for one digit more accuracy
  particle_mass_magnitude = 10.d0**(INT(log10(particle_mass)-3))

  diameter_round = DBLE(NINT(diameter/diameter_magnitude)&
		  *diameter_magnitude)
  particle_mass_round = DBLE(NINT(particle_mass/particle_mass_magnitude)&
		   *particle_mass_magnitude)

		   
  if (verbose >= 4) print*, "freq,t,as_ratio,diameter, particle_mass, ptype"
  if (verbose >= 4) print*, freq,t,as_ratio,diameter, particle_mass, ptype

  if (verbose >= 4) print*, "freq_round,t_round,as_ratio_round,diameter_round, particle_mass_round, ptype"
  if (verbose >= 4) print*, freq_round,t_round,as_ratio_round,diameter_round, particle_mass_round, ptype


 call sql_tools_get_entry(par_name,freq_round,t_round,as_ratio_round,diameter_round,&
    particle_mass_round,ptype,&
    scatter_matrix,extinct_matrix,emis_vector,found)

  if (.not. found) then
    if (verbose >= 3) call report(info,'Did not find particle in databse ', nameOfRoutine)
    call tmatrix_refIndex(freq_round,t_round,as_ratio_round,diameter_round,&
      particle_mass_round,ptype,&
      scatter_matrix,extinct_matrix,emis_vector)

    
    call sql_tools_write_entry(par_name,freq_round,t_round,as_ratio_round,diameter_round,&
      particle_mass_round,ptype,&
      scatter_matrix,extinct_matrix,emis_vector)

  end if



  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return
end subroutine tmatrix_sql

