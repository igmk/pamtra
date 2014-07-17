module descriptor_file

  use kinds, only: dbl, &      ! integer parameter specifying double precision
                   long        ! integer parameter specifying long integer
  use settings, only: descriptor_file_name
  use report_module

  implicit none

  character(len=15), dimension(:),allocatable   :: hydro_name_arr           ! hydrometeor name
  real(kind=dbl), dimension(:,:,:,:),allocatable      :: as_ratio_arr             ! aspect ratio
  integer(kind=long), dimension(:),allocatable  :: liq_ice_arr              ! liquid = 1; ice = -1
  real(kind=dbl), dimension(:,:,:,:),allocatable      :: rho_ms_arr               ! density of the particle [kg/m^3]
  real(kind=dbl), dimension(:,:,:,:),allocatable      :: canting_arr    ! canting angle of the particle [kg/m^3]
  real(kind=dbl), dimension(:,:,:,:),allocatable      :: a_ms_arr, b_ms_arr       ! Mass-size parameter a [kg/m^(1/b)] and b [#] 
  real(kind=dbl), dimension(:,:,:,:),allocatable      :: alpha_as_arr, beta_as_arr! Area-size parameter alpha [m^(2-b)] and beta [#] 

  integer(kind=long), dimension(:),allocatable  :: moment_in_arr            ! Moments input via "PAMTRA input file" (1=Ntot; 2reff; 3=q; 12=ntot&reff; 13=ntot&q; 23=reff&q)
  integer(kind=long), dimension(:,:,:,:),allocatable  :: nbin_arr                 ! Number of bins for the drop-size distribution
  character(len=15), dimension(:),allocatable   :: dist_name_arr            ! name of the distribution
  character(len=15), dimension(:),allocatable   :: scat_name_arr            ! name of the scattering model
  character(len=30), dimension(:),allocatable   :: vel_size_mod_arr         ! name of the velocity-size model to be used
  real(kind=dbl), dimension(:,:,:,:),allocatable      :: p_1_arr, p_2_arr         ! Drop-size parameters from hydrometeor descriptor file
  real(kind=dbl), dimension(:,:,:,:),allocatable      :: p_3_arr, p_4_arr         ! Drop-size parameters from hydrometeor descriptor file
  real(kind=dbl), dimension(:,:,:,:),allocatable      :: d_1_arr, d_2_arr         ! Minimum and maximum particle diameters
  integer(kind=long)  :: n_hydro                                            ! Number of hydrometeors

  real(kind=dbl), dimension(:),allocatable      :: q_h_arr                  ! Specific hydrometeor concentration [kg/kg]
  real(kind=dbl), dimension(:),allocatable      :: n_tot_arr                ! Total hydrometeor number concentration [#/kg]
  real(kind=dbl), dimension(:),allocatable      :: r_eff_arr                ! Effective radius [m]


 contains

subroutine read_descriptor_file(errorstatus)


  implicit none

! Error handling

  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'read_descriptor_file'

! work variables
  logical             :: file_exists
  character(len=1)    :: work2
  integer             :: work1, i

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

! Check file existence
  INQUIRE(FILE=trim(descriptor_file_name), EXIST=file_exists)
  if (.not.(file_exists)) then
    msg = 'descriptor file does not exist!'
    errorstatus = fatal
    call report(errorstatus, msg, nameOfRoutine)
    return
  endif


! Get number of hydrometeors
  n_hydro = 0


  open(unit=111,file=trim(descriptor_file_name))
  do
    read(111,*,IOSTAT=work1)  work2
    if (work1 /= 0) exit
    if (work2 /= '!')  n_hydro = n_hydro + 1
  end do
  rewind(111)

  call allocate_descriptor_file(err)
  if (err /= 0) then
      msg = 'error in allocate_descriptor_file!'
      call report(err, msg, nameOfRoutine)
      errorstatus = err
      return
  end if  
! Read the hydrometeor descriptors
  i=1
  do 
    read(111,*,IOSTAT=work1) work2
    if (work1 /= 0) exit
    if (work2 /= '!') then
      backspace(111)
      read(111,*) hydro_name_arr(i), as_ratio_arr(1,1,1,i), liq_ice_arr(i),rho_ms_arr(1,1,1,i), a_ms_arr(1,1,1,i), &
                b_ms_arr(1,1,1,i), alpha_as_arr(1,1,1,i), beta_as_arr(1,1,1,i), moment_in_arr(i), nbin_arr(1,1,1,i),&
                dist_name_arr(i), p_1_arr(1,1,1,i), p_2_arr(1,1,1,i), p_3_arr(1,1,1,i),  &
                p_4_arr(1,1,1,i), d_1_arr(1,1,1,i), d_2_arr(1,1,1,i), scat_name_arr(i), vel_size_mod_arr(i)
      if (verbose >= 4 .and. i == 1) write(6,'(a15,7(a12),a15,6(a12),2(a15))'),'hydro_name','as_ratio','liq_ice','rho_ms','a_ms',&
                                     'b_ms','alpha_as','beta_as','moment_in','nbin','dist_name','p_1',&
                                      'p_2','p_3','p_4','d_1','d_2','scat_name', 'vel_size_mod'
      if (verbose >= 4) write(6,'(a15,e12.3,i12,5(e12.3),2(i12),a15,6(e12.3),2(a15))') &
                trim(hydro_name_arr(i)), as_ratio_arr(1,1,1,i), liq_ice_arr(i),rho_ms_arr(1,1,1,i), a_ms_arr(1,1,1,i), &
                b_ms_arr(1,1,1,i), alpha_as_arr(1,1,1,i), beta_as_arr(1,1,1,i), moment_in_arr(i), nbin_arr(1,1,1,i),&
                trim(dist_name_arr(i)), p_1_arr(1,1,1,i), p_2_arr(1,1,1,i), p_3_arr(1,1,1,i),&
                p_4_arr(1,1,1,i), d_1_arr(1,1,1,i), d_2_arr(1,1,1,i), trim(scat_name_arr(i)), trim(vel_size_mod_arr(i))
      i = i+1
    endif
  end do

  close(unit=111)

  errorstatus = err
  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)


  return
end subroutine read_descriptor_file

subroutine allocate_descriptor_file(errorstatus)
  integer(kind=long), intent(out) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=30) :: nameOfRoutine = 'allocate_descriptor_file'

  if (verbose >= 5) call report(info,'Start of ', nameOfRoutine)
  call assert_true(err,n_hydro>0,&
      "n_hydro must be greater zero")  
  if (err > 0) then
      errorstatus = fatal
      msg = "assertation error"
      print*, "n_hydro", n_hydro, err
      call report(errorstatus, msg, nameOfRoutine)
      return
  end if    

! Allocate variables
  allocate(hydro_name_arr(n_hydro))
  allocate(as_ratio_arr(1,1,1,n_hydro))
  allocate(liq_ice_arr(n_hydro))
  allocate(rho_ms_arr(1,1,1,n_hydro))
  allocate(canting_arr(1,1,1,n_hydro))
  allocate(a_ms_arr(1,1,1,n_hydro))
  allocate(b_ms_arr(1,1,1,n_hydro))
  allocate(alpha_as_arr(1,1,1,n_hydro))
  allocate(beta_as_arr(1,1,1,n_hydro))
  allocate(moment_in_arr(n_hydro))
  allocate(nbin_arr(1,1,1,n_hydro))
  allocate(dist_name_arr(n_hydro))
  allocate(p_1_arr(1,1,1,n_hydro))
  allocate(p_2_arr(1,1,1,n_hydro))
  allocate(p_3_arr(1,1,1,n_hydro))
  allocate(p_4_arr(1,1,1,n_hydro))
  allocate(d_1_arr(1,1,1,n_hydro))
  allocate(d_2_arr(1,1,1,n_hydro))
  allocate(scat_name_arr(n_hydro))
  allocate(vel_size_mod_arr(n_hydro))

  !fill dummy values
  hydro_name_arr = "-9999"
  as_ratio_arr = -9999.d0
  liq_ice_arr = -9999
  rho_ms_arr = -9999.d0
  canting_arr = -9999.d0
  a_ms_arr = -9999.d0
  b_ms_arr = -9999.d0
  alpha_as_arr = -9999.d0
  beta_as_arr = -9999.d0
  moment_in_arr = -9999
  nbin_arr = -9999.d0
  dist_name_arr = "-9999"
  p_1_arr = -9999.d0
  p_2_arr = -9999.d0
  p_3_arr = -9999.d0
  p_4_arr = -9999.d0
  d_1_arr = -9999.d0
  d_2_arr = -9999.d0
  scat_name_arr = "-9999"
  vel_size_mod_arr = "-9999"


  errorstatus = err
  if (verbose >= 5) call report(info,'End of ', nameOfRoutine)
  return
end subroutine allocate_descriptor_file



subroutine deallocate_descriptor_file()

  character(len=30) :: nameOfRoutine = 'deallocate_descriptor_file'
  if (verbose >= 5) call report(info,'Start of ', nameOfRoutine)

  if (allocated(hydro_name_arr)) deallocate(hydro_name_arr)
  if (allocated(as_ratio_arr)) deallocate(as_ratio_arr)
  if (allocated(liq_ice_arr)) deallocate(liq_ice_arr)
  if (allocated(rho_ms_arr)) deallocate(rho_ms_arr)
  if (allocated(canting_arr)) deallocate(canting_arr)
  if (allocated(a_ms_arr)) deallocate(a_ms_arr)
  if (allocated(b_ms_arr)) deallocate(b_ms_arr)
  if (allocated(moment_in_arr)) deallocate(moment_in_arr)
  if (allocated(nbin_arr)) deallocate(nbin_arr)
  if (allocated(dist_name_arr)) deallocate(dist_name_arr)
  if (allocated(p_1_arr)) deallocate(p_1_arr)
  if (allocated(p_2_arr)) deallocate(p_2_arr)
  if (allocated(p_3_arr)) deallocate(p_3_arr)
  if (allocated(p_4_arr)) deallocate(p_4_arr)
  if (allocated(d_1_arr)) deallocate(d_1_arr)
  if (allocated(d_2_arr)) deallocate(d_2_arr)
  if (allocated(scat_name_arr)) deallocate(scat_name_arr)
  if (allocated(alpha_as_arr)) deallocate(alpha_as_arr)
  if (allocated(beta_as_arr)) deallocate(beta_as_arr)
  if (allocated(vel_size_mod_arr)) deallocate(vel_size_mod_arr)

  if (verbose >= 5) call report(info,'End of ', nameOfRoutine)

end subroutine deallocate_descriptor_file

!only for debugging purposes
subroutine printDescriptorVars()
  
  print*, "hydro_name_arr: ", SHAPE(hydro_name_arr), ";  ", hydro_name_arr
  print*, "as_ratio_arr: ", SHAPE(as_ratio_arr), ";  ", as_ratio_arr
  print*, "canting_arr: ", SHAPE(canting_arr), ";  ", canting_arr
  print*, "liq_ice_arr: ", SHAPE(liq_ice_arr), ";  ", liq_ice_arr
  print*, "rho_ms_arr: ", SHAPE(rho_ms_arr), ";  ", rho_ms_arr
  print*, "a_ms_arr: ", SHAPE(a_ms_arr), ";  ", a_ms_arr
  print*, "b_ms_arr: ", SHAPE(b_ms_arr), ";  ", b_ms_arr
  print*, "moment_in_arr: ", SHAPE(moment_in_arr), ";  ", moment_in_arr
  print*, "nbin_arr: ", SHAPE(nbin_arr), ";  ", nbin_arr
  print*, "dist_name_arr: ", SHAPE(dist_name_arr), ";  ", dist_name_arr
  print*, "p_1_arr: ", SHAPE(p_1_arr), ";  ", p_1_arr
  print*, "p_2_arr: ", SHAPE(p_2_arr), ";  ", p_2_arr
  print*, "p_3_arr: ", SHAPE(p_3_arr), ";  ", p_3_arr
  print*, "p_4_arr: ", SHAPE(p_4_arr), ";  ", p_4_arr
  print*, "d_1_arr: ", SHAPE(d_1_arr), ";  ", d_1_arr
  print*, "d_2_arr: ", SHAPE(d_2_arr), ";  ", d_2_arr
  print*, "scat_name_arr: ", SHAPE(scat_name_arr), ";  ", scat_name_arr
  print*, "alpha_as_arr: ", SHAPE(alpha_as_arr), ";  ", alpha_as_arr
  print*, "beta_as_arr: ", SHAPE(beta_as_arr), ";  ", beta_as_arr
  print*, "vel_size_mod_arr: ", SHAPE(vel_size_mod_arr), ";  ", vel_size_mod_arr

end subroutine printDescriptorVars

end module descriptor_file