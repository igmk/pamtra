subroutine run_rt3(nx,ny,fi,freq,frq_str)
  use kinds
  use constants
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
  use double_moments_module 
!  use mod_io_strings

  implicit none

  integer, intent(in) :: nx,ny,fi 
  real(kind=dbl), intent(in) :: freq ! frequency [GHz]
  character(6), intent(in) :: frq_str !from commandline

  integer, dimension(maxlay) :: OUTLEVELS
  integer :: ise, imonth ! filehandle for the emissivity data
  integer :: nz

  real(kind=dbl), dimension(maxv) :: MU_VALUES
  real(kind=dbl) :: wavelength       ! microns
  real(kind=dbl) :: GROUND_TEMP, ground_albedo
  real(kind=sgl) :: lat, lon, lfrac
  real(kind=dbl) :: wind10u, wind10v
  real(kind=dbl) :: land_emissivity

  complex(kind=dbl) :: eps_water, & ! function to calculate the dielectic properties of (salt)water
       epsi         ! result of function eps_water
  complex(kind=dbl) :: GROUND_INDEX

  character(300) :: OUT_FILE_PAS, OUT_FILE_ACT !file names if no nc
  character(80) :: femis ! filename for the emissivity databases
  character(3) :: xstr, ystr

  wavelength = c / (freq*1.d3)   ! microns

  write(xstr, '(i3.3)') profiles(nx,ny)%isamp
  write(ystr, '(i3.3)') profiles(nx,ny)%jsamp

  if (verbose .gt. 0) print*, "calculating: ", frq_str, " Y:",ny, " of ", ngridy, "X:", nx, " of ", ngridx

  ground_temp = profiles(nx,ny)%temp_lev(0)       ! K
  lat = profiles(nx,ny)%latitude                  ! °
  lon = profiles(nx,ny)%longitude                 ! °
  lfrac = profiles(nx,ny)%land_fraction
  relhum_lev = profiles(nx,ny)%relhum_lev         ! %
  press_lev = profiles(nx,ny)%press_lev           ! Pa
  temp_lev = profiles(nx,ny)%temp_lev             ! K
  hgt_lev = profiles(nx,ny)%hgt_lev               ! m

  cwc_q = profiles(nx,ny)%cloud_water_q           ! kg/kg
  iwc_q = profiles(nx,ny)%cloud_ice_q             ! kg/kg
  rwc_q = profiles(nx,ny)%rain_q                  ! kg/kg
  swc_q = profiles(nx,ny)%snow_q                  ! kg/kg
  gwc_q = profiles(nx,ny)%graupel_q               ! kg/kg

  if (n_moments .eq. 2) then
     hwc_q = profiles(nx,ny)%hail_q              ! kg/kg
     cwc_n = profiles(nx,ny)%cloud_water_n       ! #/kg
     iwc_n = profiles(nx,ny)%cloud_ice_n         ! #/kg
     rwc_n = profiles(nx,ny)%rain_n              ! #/kg
     swc_n = profiles(nx,ny)%snow_n              ! #/kg
     gwc_n = profiles(nx,ny)%graupel_n           ! #/kg
     hwc_n = profiles(nx,ny)%hail_n              ! #/kg
  end if

  press = profiles(nx,ny)%press                   ! Pa
  temp = profiles(nx,ny)%temp                     ! K
  relhum = profiles(nx,ny)%relhum                 ! %
  vapor_pressure = profiles(nx,ny)%vapor_pressure ! Pa
  rho_vap = profiles(nx,ny)%rho_vap               ! kg/m^3
  q_hum = profiles(nx,ny)%q_hum                   ! kg/kg


  if (verbose .gt. 1) print*, nx,ny, 'type to local variables done' 

  ! Determine surface properties

  if (lfrac .ge. 0.5 .and. lfrac .le. 1.0) then
     ground_type = 'S' ! changed to specular after advice of cathrine prigent
     ise=13
     read(month,'(i2)') imonth
     if (imonth .ge. 7 .and. imonth .le. 12) then
        femis = data_path(:len_trim(data_path))//'/emissivity/ssmi_mean_emis_92'//month//'_direct'
     else if (imonth .ge. 1 .and. imonth .lt. 7) then
        femis = data_path(:len_trim(data_path))//'/emissivity/ssmi_mean_emis_93'//month//'_direct'
     else
        print*, nx,ny, "Warning: No emissivity data found!"
        stop
     end if
     open(ise,file=trim(femis),status='old',form='unformatted',&
          access='direct',recl=28)
     ! land_emis could give polarized reflectivities

     call land_emis(ise,lon,lat,freq,land_emissivity)
     close(ise)
     if (land_emissivity .lt. 0.01) then
	     land_emissivity = 0.94d0
     end if
     ground_albedo = 1.d0 - land_emissivity

  else if (lfrac .ge. 0.0 .and. lfrac .lt. 0.5) then
     ! computing the refractive index of the sea (Fresnel) surface
     ground_type = 'O'
     ground_albedo = 1.0d0
     epsi = eps_water(salinity, ground_temp - 273.15d0, freq)
     ground_index = dconjg(sqrt(epsi))
  else
     ! this is for ground_type specified in run_params.nml
     ground_albedo = 1.d0 - emissivity
  end if

  if (verbose .gt. 1) print*, nx,ny, 'Surface emissivity calculated!'

  ! gaseous absorption
  ! 
  ! kextatmo   extinction by moist air [Np/m]
  !
  if (lgas_extinction) then
     !returns kextatmo!
     call get_atmosg(freq)
  else
     kextatmo = 0.0D0 ! for the whole column
  end if

  if (verbose .gt. 1) print*, nx,ny, 'Gas absorption calculated'


  ! hydrometeor extinction desired


  if (lhyd_extinction) call hydrometeor_extinction(freq,frq_str)

  !
  if (dump_to_file) call dump_profile(frq_str)

  !&&&&&&&&   I/O FILE NAMES   &&&&&&&&&&&&&&&&&&

  OUT_FILE_PAS = output_path(:len_trim(output_path))//"/"//&
       date_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_passive"

  OUT_FILE_ACT = output_path(:len_trim(output_path))//"/"//&
       date_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_active"


  if (active) then
     call calculate_active(OUT_FILE_ACT,freq,hgt(nx,ny,:),Ze(nx,ny,:,fi),Attenuation_atmo(nx,ny,:,fi),&
          Attenuation_hydro(nx,ny,:,fi))
     if (verbose .gt. 1) print*, nx,ny, 'calculate_active done'
  end if

  if (write_nc) then
     !      Output integrated quantities
     call collect_boundary_output(lon,lat,lfrac,&
          profiles(nx,ny)%iwv, profiles(nx,ny)%cwp,profiles(nx,ny)%iwp,profiles(nx,ny)%rwp,profiles(nx,ny)%swp, &
          profiles(nx,ny)%gwp,profiles(nx,ny)%hwp,profiles(nx,ny)%isamp,profiles(nx,ny)%jsamp,nx,ny)
     if (verbose .gt. 1) print*, nx,ny, 'collect_boundary_output done'
  end if

  ! find the output level
  ! in rt3 layers are reversed

  if (obs_height .gt. 99999. .or. obs_height .gt. hgt_lev(nlyr)) then
     outlevels(1) = 1
  else if (obs_height .lt. 0.1 .or. obs_height .lt. hgt_lev(1)) then
     outlevels(1) = nlyr + 1
  else
     out_search: do nz = 1, nlyr
        if (hgt_lev(nz) .ge. obs_height) then
           if (abs(hgt_lev(nz) - obs_height) .lt. abs(hgt_lev(nz-1) - obs_height)) then
              outlevels(1) = nlyr-nz+1
           else
              outlevels(1) = nlyr-nz+2
           end if
           exit out_search
        end if
     end do out_search
  end if

  OUTLEVELS(2) = nlyr+1    ! this is the bottom

  if (passive .eqv. .true.) then

     if (verbose .gt. 1) print*, nx,ny, "Entering rt3 ...."

     call RT3(NSTOKES, NUMMU, AZIORDER, MU_VALUES, src_code,     &
          out_file_pas, QUAD_TYPE, deltam, DIRECT_FLUX,     &
          DIRECT_MU, GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO,         &
          GROUND_INDEX, SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,          &
          NOUTLEVELS, OUTLEVELS, NUMAZIMUTHS,&
          nx,ny,fi)

     !calculate human readable angles!
     angles_deg(1:NUMMU) = 180-(180.*acos(MU_VALUES(NUMMU:1:-1))/pi)
     angles_deg(1+NUMMU:2*NUMMU) = (180.*acos(MU_VALUES(1:NUMMU))/pi)

     if (verbose .gt. 1) print*, nx,ny, "....rt3 finished"

  end if

end subroutine run_rt3
