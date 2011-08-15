program pamtra

  use kinds
  use constants
  use nml_params !all settings go here
  use vars_atmosphere !input variables and reading routine
  use vars_output !output variables
  use double_moments_module 
  use mod_io_strings

  !     Radiative transfer code to process COSMO-model derived profiles   
  !     The code reads a full COSMO grid and computes for each profile the  
  !     radiative transfer for the given frequency                        
  !                                                                      
  !     This code is completely self contained and needs no databases     
  !     or lookup tables to run.  By convention, the quantities followed  
  !     by "_lev" are given at the layer heights while the quantitites
  !     w/o "_lev" are layer average quantities                           
  !                                                                       
  !                    **  RADTRAN I/O SPECIFICATIONS  **

  implicit none





!!! NOT YET SORTED STUFF

!   integer :: jj, nf,  nz, 

integer ::  fi,nx, ny!, &nlev,
!              offset1,offset2, length1, length2,&
!              ise, imonth ! filehandle for the emissivity data
! 
! !   integer :: NUMRAD, NLEGENcw, NLEGENci, NLEGENgr, &
! !        NLEGENsn, NLEGENrr
! 
! 
!   integer :: NP_LUT, N_temp
! 
!   integer :: j_temp, ind_temp
! 
!   integer :: i_bot, i_top, nnz, N_layer_new, ss, tt, length
! 





!    real(kind=dbl) ::  gammln
! 


! 
!   real(kind=dbl) :: E1, E2
! 
!   real(kind=dbl) :: RAD1, RAD2, refre, refim,&
!        N_0sr, Coeff_corr, AD, BD,  &
!        n0S, lambda_D, tmp, N_0snowD, N_0grauD,              &
!        Coeff_snow, a_mgraup, b_g, b_snow,        &
!        a_msnow, Coeff_grau, den_liq, drop_mass, del_r, den_ice,&
!        densnow, fvol_ice
! 
!   real(kind=dbl) ::  ABSIND, ABSCOF
! 
!   real(kind=dbl) :: atm_ext, kextcw, salbcw, asymcw, kextrr, salbrr, asymrr,  &
!        kextci, salbci, asymci, kextsn, salbsn, asymsn, kextgr, salbgr,   &
!        asymgr, salbhl, asymhl, backcw, backrr, backci, backsn, backgr
!                                                       
! 
!   real(kind=dbl) :: mu, D_0, D_max, D_min, A1, A2


! 
! 
! !  real(kind=dbl) :: emissivity     ! land surface reflectivity
! 
!   real(kind=dbl) :: Angle_view, Angle_zenith, I1, I2,     &
!        Angle_viewdeg, Upar, Vpar, Ds_cloud_slant, Ds_obs_point,&
!        DH_bot_intersect, H_bot, H_top, DH_top_intersect, K_extbot,&
!        K_exttop, T_bot, T_top                                            
! 
!   real(kind=dbl), dimension(maxleg) :: LEGENcw, LEGENrr, LEGENci, LEGENgr, LEGENsn,       &
!        LEGEN2cw, LEGEN2rr, LEGEN2ci, LEGEN2gr, LEGEN2sn,  &
!        LEGEN3cw, LEGEN3rr, LEGEN3ci, LEGEN3gr, LEGEN3sn,  &
!        LEGEN4cw, LEGEN4rr, LEGEN4ci, LEGEN4gr, LEGEN4sn
! 
! 
! 
! 
! !   complex(kind=dbl) :: GROUND_INDEX
!   complex(kind=dbl) :: MINDEX, m_air, m_MG, m_ice

!   character :: rLWC_str*4                                                      

!   character(300) :: OUT_FILE_PAS, OUT_FILE_ACT, tmp_file1, 

!   character :: ssstr*1, ttstr*1, Anglestr*4, FILEOUT3D*65
! 
! 
! 


  character(300) ::nc_out_file



! temporary variables


!   real(kind=dbl) :: lwc, iwc, rwc, gwc, swc



!!! INTERNAL "HANDLE COMMAND LINE PARAMETERS" !!! 

  integer :: inarg, ff
  character(40) :: gitHash, gitVersion
  character(6) :: formatted_frqstr !function call

!!! SET BY "HANDLE COMMAND LINE PARAMETERS" !!! 

  character(99)  :: input_file !name of profile
  character(300) :: namelist_file
  character(6), dimension(maxfreq) :: frqs_str !from commandline
  character(5*7) :: frq_str_list ! for the filename only!



!!! HANDLE COMMAND LINE PARAMETERS !!!

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

!!! READ NAMELIST FILE !!!
call nml_params_read(namelist_file) !from nml_params.f90

!!! READ MOMENTS FILE !!!
if (n_moments .eq. 2) call double_moments_module_read(moments_file) !from double_moments_module.f90

!!! READ PROFILES !!!
call vars_atmosphere_read_profile(input_file)

! now allocate variables
call allocate_vars


  !                                                                       
  !     This GCE model format does not have all the fields expected by    
  !     the radiative transfer code (i.e. total pressure, and water vapor 
  !     pressure for this model).  Assign/compute the missing fields first
  !make layer averages
  call get_atmosG0

  if (write_nc .eqv. .false.) call mod_io_strings_get_filename()

  if (verbose .gt. 1) print*, 'Start loop over frequencies!'

grid_f: do fi =1, nfrq


  if (verbose .gt. 1) print*, 'Start loop over profiles!'

  grid_y: do ny = 1, ngridy !ny_in, ny_fin  
    grid_x: do nx = 1, ngridx !nx_in, nx_fin   



call run_rt3(nx,ny,fi,frqs_str)




    end do grid_x
  end do grid_y
end do grid_f

  if (write_nc) then
    nc_out_file = output_path(1:len_trim(output_path))//"/"//trim(input_file(1:len_trim(input_file)-4))//&
                    frq_str_list(1:len_trim(frq_str_list))//'_res.nc'
    call write_nc_results(nc_out_file)
  end if

end program pamtra
