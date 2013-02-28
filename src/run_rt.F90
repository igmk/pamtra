subroutine run_rt(errorstatus, nx,ny,fi,freq,frq_str)

  use kinds, only: long, dbl
  use constants, only: c,&
                        pi,&
                        errorstatus_fatal
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
  use double_moments_module 
!  use mod_io_strings

  implicit none

#include "error_report.interface"

  integer(kind=long), intent(in) :: nx,ny,fi
  real(kind=dbl), intent(in) :: freq ! frequency [GHz]
  character(8), intent(in) :: frq_str !from commandline

  integer(kind=long), dimension(maxlay) :: OUTLEVELS
  integer(kind=long) :: ise, imonth ! filehandle for the emissivity data
  integer(kind=long) :: nz

  real(kind=dbl), dimension(maxv) :: MU_VALUES
  real(kind=dbl) :: wavelength       ! microns
  real(kind=dbl) :: GROUND_TEMP, ground_albedo

  real(kind=dbl) :: land_emissivity

  complex(kind=dbl) :: eps_water, & ! function to calculate the dielectic properties of (salt)water
       epsi         ! result of function eps_water
  complex(kind=dbl) :: GROUND_INDEX

  character(300) :: OUT_FILE_PAS, OUT_FILE_ACT !file names if no nc
  character(80) :: femis ! filename for the emissivity databases
  character(3) :: xstr, ystr

  ! i/o-length test for the emissivity file
  integer(kind=long) :: iolsgl

  ! Error handling

  integer(kind=long), intent(out) :: errorstatus
  character(len=80) :: ErrMsg
  character(len=14) :: NameOfRoutine = 'run_rt'

  inquire(iolength=iolsgl) 1._sgl

  wavelength = c / (freq*1.d3)   ! microns
  GROUND_TEMP = temp_lev(0)

  if (verbose .gt. 0) print*, "calculating: ", frq_str, " Y:",ny, " of ", ngridy, "X:", nx, " of ", ngridx


  write(xstr, '(i3.3)') model_i
  write(ystr, '(i3.3)') model_j

  ! This GCE model format does not have all the fields expected by    
  ! the radiative transfer code (i.e. total pressure, and water vapor 
  ! pressure for this model).  Assign/compute the missing fields first
  ! make layer averages

  call get_atmosG0

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
          access='direct',recl=iolsgl*7)
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
    call get_gasabs(errorstatus,freq)
    if (errorstatus == errorstatus_fatal) Then
       errmsg = 'error in get_gasabs'
       call error_report(errorstatus, errmsg, NameOfRoutine)
       return
    end if
  else
     kextatmo = 0.0D0 ! for the whole column
  end if
  !save atmospheric attenuation and height for radar
  if (active) then
    Att_atmo(nx,ny,:,fi)  = 10*log10(exp(kextatmo*delta_hgt_lev))
    radar_hgt(nx,ny,:) = hgt(:)
  end if


  if (verbose .gt. 1) print*, nx,ny, 'Gas absorption calculated'


  ! hydrometeor extinction desired
  if (lhyd_extinction) then
    if (rt_mode .eq. 'rt3') then
        call hydrometeor_extinction_rt3(freq,frq_str)
    elseif (rt_mode .eq. 'rt4') then
        call hydrometeor_extinction_rt4(freq,frq_str,nx,ny,fi)!hier nx, ny
    end if
  end if



  !
  if (dump_to_file) call dump_profile(frq_str)

  !&&&&&&&&   I/O FILE NAMES   &&&&&&&&&&&&&&&&&&

  OUT_FILE_PAS = output_path(:len_trim(output_path))//"/"//&
       date_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_passive"

  OUT_FILE_ACT = output_path(:len_trim(output_path))//"/"//&
       date_str//'x'//xstr//'y'//ystr//'f'//frq_str//"_active"
! 
!   if ((active) .and. ((radar_mode == "simple") .or. (radar_mode == "splitted")))  then
!      call calculate_active(OUT_FILE_ACT,freq,&
!           Ze(nx,ny,:,fi),Ze_cw(nx,ny,:,fi),Ze_rr(nx,ny,:,fi),Ze_ci(nx,ny,:,fi),&
!           Ze_sn(nx,ny,:,fi),Ze_gr(nx,ny,:,fi),Ze_ha(nx,ny,:,fi),&
!           Att_atmo(nx,ny,:,fi),Att_hydro(nx,ny,:,fi),Att_cw(nx,ny,:,fi),Att_rr(nx,ny,:,fi),&
!           Att_ci(nx,ny,:,fi),Att_sn(nx,ny,:,fi),Att_gr(nx,ny,:,fi),Att_ha(nx,ny,:,fi))
!      if (verbose .gt. 1) print*, nx,ny, 'calculate_active done'
!      
!   end if

  !save active to ASCII
  if (active .and. (write_nc .eqv. .false.) .and. (in_python .eqv. .false.)) then
    call save_active(OUT_FILE_ACT,nx,ny,fi)
  end if


  if (write_nc) then
     !      Output integrated quantities
     call collect_boundary_output(lon,lat,lfrac,&
          iwv, cwp,iwp,rwp,swp, &
          gwp,hwp,model_i,model_j,nx,ny)
     if (verbose .gt. 1) print*, nx,ny, 'collect_boundary_output done'
  end if

  ! find the output level
  ! in rt3 and rt4 layers are reversed

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

    if (rt_mode .eq. 'rt3') then
        if (verbose .gt. 1) print*, nx,ny, "Entering rt3 ...."

        call RT3(NSTOKES, NUMMU, AZIORDER, MU_VALUES, src_code, &
          out_file_pas, QUAD_TYPE, deltam, DIRECT_FLUX,     &
          DIRECT_MU, GROUND_TEMP, GROUND_TYPE, GROUND_ALBEDO,  &
          GROUND_INDEX, SKY_TEMP, WAVELENGTH, UNITS, OUTPOL,  &
          NOUTLEVELS, OUTLEVELS, nx,ny,fi)

        if (verbose .gt. 1) print*, nx,ny, "....rt3 finished"
    elseif (rt_mode .eq. 'rt4') then
        if (verbose .gt. 1) print*, nx,ny, "Entering rt4 ...."

        call rt4(nstokes,nummu,mu_values,out_file_pas,quad_type,ground_temp,&
        ground_type,ground_albedo,ground_index,sky_temp,&
        wavelength,units,outpol,noutlevels,outlevels,nx,ny,fi)

        if (verbose .gt. 1) print*, nx,ny, "....rt4 finished"

    else
        print*, 'no rt_mode selected'
        stop
    end if
    !calculate human readable angles!
    angles_deg(1:NUMMU) = 180-(180.*acos(MU_VALUES(NUMMU:1:-1))/pi)
    angles_deg(1+NUMMU:2*NUMMU) = (180.*acos(MU_VALUES(1:NUMMU))/pi)

  end if


end subroutine run_rt
