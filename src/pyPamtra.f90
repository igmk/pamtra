subroutine pyPamtra(input_file,namelist_file,frequency)

! ,&
! in_year,in_month,in_day,in_time,in_ngridx,in_ngridy,in_nlyr,in_deltax,in_deltay,in_&
! in_lat,in_lon,in_lfrac,in_relhum_lev,in_press_lev,in_temp_lev,in_hgt_lev,in_&
! in_model_i,in_model_j,in_wind10u,in_wind10v,in_iwv,in_cwp,in_&
! in_iwp,in_rwp,in_swp,in_gwp,in_hwp,in_cwc_q,in_iwc_q,in_rwc_q,in_swc_q,in_gwc_q)

  use kinds
  use constants !physical constants live here
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
  use vars_profile
  use double_moments_module !double moments variables are stored here
  use mod_io_strings !some strings for nice filenames




  !     The code reads a full (e.g. COSMO) grid and computes for each 
  !     profile the radiative transfer for the given frequencies                       
  !                                                                      
  !     By convention, the quantities followed  by "_lev"
  !     are given at the layer heights while the quantitites w/o
  !     "_lev" are layer average quantities                           
  !                                                                       

  implicit none


!!! internal "handle command line parameters" !!! 

  integer :: inarg, ff
  character(40) :: gitHash, gitVersion
  character(6) :: formatted_frqstr !function call

!!! set by "handle command line parameters" !!! 

  character(99),intent(in)  :: input_file !name of profile
  character(300),intent(in) :: namelist_file
  real, intent(in) :: frequency

!f2py intent(in) :: input_file,namelist_file,frequency

  character(6), dimension(maxfreq) :: frqs_str !from commandline
  character(7) :: frq_str_s,frq_str_e

!!!loop variables
  integer ::  fi,nx, ny

!!!output variables
  character(300) ::nc_out_file

print *,"Hello"
print *,input_file,namelist_file,frequency

  !get git data
  call versionNumber(gitVersion,gitHash)

print *,gitVersion,gitHash


nfrq = 1
allocate(freqs(nfrq))
frqs_str(1) = "024.00"
freqs(1) = frequency

print *,frqs_str, freqs


! !!! read variables from namelist file
!   call nml_params_read(namelist_file) !from nml_params.f90

    !set namelist defaults!
    verbose=2

!     write_nc=.true.
    dump_to_file=.false.
!     input_path='profile/'
!     output_path='output/'
    tmp_path='/tmp/'
!     data_path='data/'
write_nc=.false.
input_path='../test/referenceProfile'
output_path='../test/tmp'
data_path='/home/mech/models/pamtra/data/'

    obs_height=833000.
    units='T'
    outpol='VH'
    freq_str=''
    file_desc=''
    creator='Pamtrauser'

    active=.true.
    passive=.true.

    ground_type='S'
    salinity=33.0
    emissivity=0.6

    lgas_extinction=.true.
    gas_mod='R98'

    lhyd_extinction=.true.
    lphase_flag = .true.

    SD_snow='Exp' 
    N_0snowDsnow=7.628 
    EM_snow='icesf' 
    SP=0.2 
    isnow_n0=1
    liu_type=8

    SD_grau='Exp' 
    N_0grauDgrau=4.0 
    EM_grau='surus'

    EM_ice='mieic'

    SD_rain='Exp' 
    N_0rainD=8.0

    n_moments=1
    moments_file='snowCRYSTAL'




print *,"read namelist"

  ! create frequency string of not set in pamtra
  if (freq_str .eq. "") then
     ! get integer and character frequencies
    frq_str_s = "_"//frqs_str(1)
     if (nfrq .eq. 1) then
       frq_str_e = ""
     else
       frq_str_e = "-"//frqs_str(nfrq)
     end if
     freq_str = frq_str_s//frq_str_e
  end if
!      frq_str_list = frq_str_list(:len_trim(frq_str_list)) // "_" //  frqs_str(ff)

  if (verbose .gt. 1) print *,"input_file: ",input_file(:len_trim(input_file)),&
       " namelist file: ",namelist_file," freq: ",freq_str

!!! read n-moments file
  if (n_moments .eq. 2) call double_moments_module_read(moments_file) !from double_moments_module.f90

!!! read the data
  call vars_profile_read_profile(input_file) !from vars_atmosphere.f90

year = profiles_year
month = profiles_month
day = profiles_day
time = profiles_time
ngridx = profiles_ngridx
ngridy = profiles_ngridy
nlyr = profiles_nlyr
deltax = profiles_deltax
deltay = profiles_deltay
  ! now allocate variables
  call allocate_vars


  if (write_nc .eqv. .false.) call mod_io_strings_get_filename()


  if (verbose .gt. 1) print*, 'Start loop over frequencies & profiles!'




  grid_f: do fi =1, nfrq
     grid_y: do ny = 1, ngridy !ny_in, ny_fin  
        grid_x: do nx = 1, ngridx !nx_in, nx_fin   

         !   ground_temp = profiles(nx,ny)%temp_lev(0)       ! K
         lat = profiles(nx,ny)%latitude                  ! °
         lon = profiles(nx,ny)%longitude                 ! °
         lfrac = profiles(nx,ny)%land_fraction
         relhum_lev = profiles(nx,ny)%relhum_lev         ! %
         press_lev = profiles(nx,ny)%press_lev           ! Pa
         temp_lev = profiles(nx,ny)%temp_lev             ! K
         hgt_lev = profiles(nx,ny)%hgt_lev               ! m

         model_i = profiles(nx,ny)%isamp
         model_j = profiles(nx,ny)%jsamp
         wind10u = profiles(nx,ny)%wind_10u
         wind10v = profiles(nx,ny)%wind_10v

         iwv = profiles(nx,ny)%iwv
         cwp = profiles(nx,ny)%cwp
         iwp = profiles(nx,ny)%iwp
         rwp = profiles(nx,ny)%rwp
         swp = profiles(nx,ny)%swp
         gwp = profiles(nx,ny)%gwp
         hwp = profiles(nx,ny)%hwp


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

           !run the model
           call run_rt3(nx,ny,fi,freqs(fi),frqs_str(fi))

        end do grid_x
     end do grid_y
  end do grid_f
end subroutine pyPamtra