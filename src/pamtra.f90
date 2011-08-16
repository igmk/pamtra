program pamtra

  use kinds
  use constants !physical constants live here
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
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

  character(99)  :: input_file !name of profile
  character(300) :: namelist_file
  character(6), dimension(maxfreq) :: frqs_str !from commandline
  character(5*7) :: frq_str_list ! for the filename only!

!!!loop variables
  integer ::  fi,nx, ny

!!!output variables
  character(300) ::nc_out_file

  !get git data
  call versionNumber(gitVersion,gitHash)

  !get command line parameters
  inarg = iargc()

  if (inarg .lt. 3) then
     print *,'Usage: pamtra profile_file namelist_file (list of frequencies)'
     print *,'Example: ./pamtra rt_comp_single.dat run_params.nml 35 94'
     print *,'See namelist file for further pamtra options'
     print *,''
     print *,'Version:  '//gitVersion
     print *,'Git Hash: '//gitHash
     stop
  else if (inarg .gt. maxfreq + 2) then
     print *,'Too many frequencies! Increase maxfreq!'
     stop
  end if
  call getarg(1,input_file)
  call getarg(2,namelist_file)

  nfrq = inarg - 2
  allocate(freqs(nfrq))

  frq_str_list = "" 
  !get integer and character frequencies
  do ff = 1, inarg-2
     call getarg(ff+2,frqs_str(ff))
     read(frqs_str(ff),*) freqs(ff)
     frqs_str(ff) = formatted_frqstr(frqs_str(ff))
     frq_str_list = frq_str_list(:len_trim(frq_str_list)) // "_" //  frqs_str(ff)
  end do

  if (verbose .gt. 1) print *,"input_file: ",input_file(:len_trim(input_file)),&
       " namelist file: ",namelist_file," freq: ",frqs_str

!!! read variables from namelist file
  call nml_params_read(namelist_file) !from nml_params.f90

!!! read n-moments file
  if (n_moments .eq. 2) call double_moments_module_read(moments_file) !from double_moments_module.f90

!!! read the data
  call vars_atmosphere_read_profile(input_file) !from vars_atmosphere.f90

  ! now allocate variables
  call allocate_vars


  ! This GCE model format does not have all the fields expected by    
  ! the radiative transfer code (i.e. total pressure, and water vapor 
  ! pressure for this model).  Assign/compute the missing fields first
  ! make layer averages
  call get_atmosG0

  if (write_nc .eqv. .false.) call mod_io_strings_get_filename()


  if (verbose .gt. 1) print*, 'Start loop over frequencies & profiles!'

  grid_f: do fi =1, nfrq
     grid_y: do ny = 1, ngridy !ny_in, ny_fin  
        grid_x: do nx = 1, ngridx !nx_in, nx_fin   

           !run the model
           call run_rt3(nx,ny,fi,frqs_str)

        end do grid_x
     end do grid_y
  end do grid_f

  if (write_nc) then
     nc_out_file = trim(output_path)//"/"//trim(input_file(1:len_trim(input_file)-4))//&
          trim(frq_str_list)//'_res.nc'
     call write_nc_results(nc_out_file)
  end if

end program pamtra