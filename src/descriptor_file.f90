module descriptor_file

  use kinds, only: dbl, &      ! integer parameter specifying double precision
                   long        ! integer parameter specifying long integer
  use settings, only: descriptor_file_name

  implicit none
  character(len=10), dimension(:),allocatable   :: hydro_name_arr           ! hydrometeor name
  real(kind=dbl), dimension(:),allocatable      :: as_ratio_arr             ! aspect ratio
  integer(kind=long), dimension(:),allocatable  :: liq_ice_arr              ! liquid = 1; ice = -1
  real(kind=dbl), dimension(:),allocatable      :: rho_ms_arr               ! density of the particle [kg/m^3]
  real(kind=dbl), dimension(:),allocatable      :: a_ms_arr, b_ms_arr       ! Mass-size parameter a [kg/m^(1/b)] and b [#] 

  integer(kind=long), dimension(:),allocatable  :: moment_in_arr            ! Moments input via "PAMTRA input file" (1=q; 2=Ntot; 3=reff; 12=q&ntot; 13=q&reff; 23=ntot&reff)
  integer(kind=long), dimension(:),allocatable  :: nbin_arr                 ! Number of bins for the drop-size distribution
  character(len=15), dimension(:),allocatable   :: dist_name_arr            ! name of the distribution
  character(len=15), dimension(:),allocatable   :: scat_name_arr            ! name of the distribution
  real(kind=dbl), dimension(:),allocatable      :: p_1_arr, p_2_arr         ! Drop-size parameters from hydrometeor descriptor file
  real(kind=dbl), dimension(:),allocatable      :: p_3_arr, p_4_arr         ! Drop-size parameters from hydrometeor descriptor file
  real(kind=dbl), dimension(:),allocatable      :: d_1_arr, d_2_arr         ! Minimum and maximum particle diameters
  integer(kind=long)  :: n_hydro                  ! Number of hydrometeors

  real(kind=dbl), dimension(:),allocatable      :: q_h_arr                  ! Specific hydrometeor concentration [kg/kg]
  real(kind=dbl), dimension(:),allocatable      :: n_tot_arr                ! Total hydrometeor number concentration [#/kg]
  real(kind=dbl), dimension(:),allocatable      :: r_eff_arr                ! Effective radius [m]

  real(kind=dbl), dimension(:),allocatable      :: t_arr                    ! Layer temperature [K]

 contains

subroutine read_descriptor_file(errorstatus)

  use report_module

  implicit none

! Error handling

  integer(kind=long), intent(out) :: errorstatus
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
    msg = 'file does not exist!'
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
print*,n_hydro,work2,work1
  end do
  rewind(111)

! Allocate variables
  allocate(hydro_name_arr(n_hydro))
  allocate(as_ratio_arr(n_hydro))
  allocate(liq_ice_arr(n_hydro))
  allocate(rho_ms_arr(n_hydro))
  allocate(a_ms_arr(n_hydro))
  allocate(b_ms_arr(n_hydro))
  allocate(moment_in_arr(n_hydro))
  allocate(nbin_arr(n_hydro))
  allocate(dist_name_arr(n_hydro))
  allocate(p_1_arr(n_hydro))
  allocate(p_2_arr(n_hydro))
  allocate(p_3_arr(n_hydro))
  allocate(p_4_arr(n_hydro))
  allocate(d_1_arr(n_hydro))
  allocate(d_2_arr(n_hydro))
  allocate(scat_name_arr(n_hydro))


! Read the hydrometeor descriptors
  i=1
  do 
    read(111,*,IOSTAT=work1) work2
    if (work1 /= 0) exit
    if (work2 /= '!') then
      backspace(111)
      read(111,*) hydro_name_arr(i), as_ratio_arr(i), liq_ice_arr(i),rho_ms_arr(i), a_ms_arr(i), &
                b_ms_arr(i), moment_in_arr(i), nbin_arr(i),dist_name_arr(i), p_1_arr(i),         &
                p_2_arr(i), p_3_arr(i), p_4_arr(i), d_1_arr(i), d_2_arr(i), scat_name_arr(i)
      i = i+1
    endif
  end do



!   write(6,*)'hydro_name',hydro_name_arr
!   write(6,*)'as_ratio ',as_ratio_arr
!   write(6,*)'liq_ice ',liq_ice_arr
!   write(6,*)'rho_ms ',rho_ms_arr
!   write(6,*)'a_ms ',a_ms_arr
!   write(6,*)'b_ms ',b_ms_arr
!   write(6,*)'moment_in ',moment_in_arr
!   write(6,*)'nbin ',nbin_arr
!   write(6,*)'dist_name ',dist_name_arr
!   write(6,*)'p_1 ',p_1_arr
!   write(6,*)'p_2 ',p_2_arr
!   write(6,*)'p_3 ',p_3_arr
!   write(6,*)'p_4 ',p_4_arr
!   write(6,*)'d_1 ',d_1_arr
!   write(6,*)'d_2 ',d_2_arr


  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)


  return
end subroutine read_descriptor_file

subroutine deallocate_descriptor_file()

  implicit none

  if (allocated(hydro_name_arr)) deallocate(hydro_name_arr)
  if (allocated(as_ratio_arr)) deallocate(as_ratio_arr)
  if (allocated(liq_ice_arr)) deallocate(liq_ice_arr)
  if (allocated(rho_ms_arr)) deallocate(rho_ms_arr)
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

end subroutine deallocate_descriptor_file

end module descriptor_file