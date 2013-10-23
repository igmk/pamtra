module vars_jacobian
  use kinds
  use settings, only : n_moments, nstokes,nummu
  use report_module
  implicit none

  integer(kind=long) :: alloc_status

  !for jacobain mode
  real(kind=dbl), allocatable, dimension(:,:,:,:,:,:) :: jac_scattermatrix
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: jac_extmatrix
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: jac_emisvec
  real(kind=dbl), allocatable, dimension(:) :: jac_kextsn
  real(kind=dbl), allocatable, dimension(:) :: jac_backsn
  real(kind=dbl), allocatable, dimension(:) :: jac_kextcw
  real(kind=dbl), allocatable, dimension(:) :: jac_backcw
  real(kind=dbl), allocatable, dimension(:) :: jac_kextrr
  real(kind=dbl), allocatable, dimension(:) :: jac_backrr
  real(kind=dbl), allocatable, dimension(:) :: jac_kextgr
  real(kind=dbl), allocatable, dimension(:) :: jac_backgr
  real(kind=dbl), allocatable, dimension(:) :: jac_kextci
  real(kind=dbl), allocatable, dimension(:) :: jac_backci
  real(kind=dbl), allocatable, dimension(:) :: jac_kextha
  real(kind=dbl), allocatable, dimension(:) :: jac_backha
  logical, allocatable, dimension(:) :: jac_hydros_present

  !jacobian mode
  real(kind=dbl), allocatable, dimension(:) :: jac_cwc_q, &
       jac_iwc_q, &
       jac_rwc_q, &
       jac_swc_q, &
       jac_gwc_q, &
       jac_hwc_q

  real(kind=dbl), allocatable, dimension(:) :: jac_cwc_n, &
       jac_iwc_n, &
       jac_rwc_n, &
       jac_swc_n, &
       jac_gwc_n, &
       jac_hwc_n

  real(kind=dbl), allocatable, dimension(:) :: jac_relhum_lev,&
       jac_temp_lev

  contains
  subroutine allocate_jacobian_vars(nlyr)
    
    integer(kind=long), intent(in) :: nlyr

    allocate(jac_scattermatrix(nlyr,nstokes,nummu,nstokes,nummu,4),stat=alloc_status)
    allocate(jac_extmatrix(nlyr,nstokes,nstokes,nummu,4),stat=alloc_status)
    allocate(jac_emisvec(nlyr,nstokes,nummu,4),stat=alloc_status)
    allocate(jac_kextsn(nlyr),stat=alloc_status)
    allocate(jac_backsn(nlyr),stat=alloc_status)
    allocate(jac_kextcw(nlyr),stat=alloc_status)
    allocate(jac_backcw(nlyr),stat=alloc_status)
    allocate(jac_kextrr(nlyr),stat=alloc_status)
    allocate(jac_backrr(nlyr),stat=alloc_status)
    allocate(jac_kextgr(nlyr),stat=alloc_status)
    allocate(jac_backgr(nlyr),stat=alloc_status)
    allocate(jac_kextci(nlyr),stat=alloc_status)
    allocate(jac_backci(nlyr),stat=alloc_status)
    allocate(jac_kextha(nlyr),stat=alloc_status)
    allocate(jac_backha(nlyr),stat=alloc_status)
    allocate(jac_hydros_present(nlyr),stat=alloc_status)

    allocate(jac_temp_lev(0:nlyr),stat=alloc_status)
    allocate(jac_relhum_lev(0:nlyr),stat=alloc_status)
    allocate(jac_cwc_q(nlyr),stat=alloc_status)
    allocate(jac_iwc_q(nlyr),stat=alloc_status)
    allocate(jac_rwc_q(nlyr),stat=alloc_status)
    allocate(jac_swc_q(nlyr),stat=alloc_status)
    allocate(jac_gwc_q(nlyr),stat=alloc_status)
    if (n_moments .eq. 2) then
      allocate(jac_hwc_q(nlyr),stat=alloc_status)
      allocate(jac_cwc_n(nlyr),stat=alloc_status)
      allocate(jac_iwc_n(nlyr),stat=alloc_status)
      allocate(jac_rwc_n(nlyr),stat=alloc_status)
      allocate(jac_swc_n(nlyr),stat=alloc_status)
      allocate(jac_gwc_n(nlyr),stat=alloc_status)
      allocate(jac_hwc_n(nlyr),stat=alloc_status)
    end if

  end subroutine allocate_jacobian_vars

  subroutine deallocate_jacobian_vars


    implicit none
    !   integer(kind=long), intent(out) :: errorstatus
    !   integer(kind=long) :: err = 0
    character(len=30) :: nameOfRoutine = 'deallocate_jacobian_vars'

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    if (allocated(jac_scattermatrix)) deallocate(jac_scattermatrix)
    if (allocated(jac_extmatrix)) deallocate(jac_extmatrix)
    if (allocated(jac_emisvec)) deallocate(jac_emisvec)
    if (allocated(jac_hydros_present)) deallocate(jac_hydros_present)
    if (allocated(jac_kextsn)) deallocate(jac_kextsn)
    if (allocated(jac_backsn)) deallocate(jac_backsn)
    if (allocated(jac_kextcw)) deallocate(jac_kextcw)
    if (allocated(jac_backcw)) deallocate(jac_backcw)
    if (allocated(jac_kextrr)) deallocate(jac_kextrr)
    if (allocated(jac_backrr)) deallocate(jac_backrr)
    if (allocated(jac_kextgr)) deallocate(jac_kextgr)
    if (allocated(jac_backgr)) deallocate(jac_backgr)
    if (allocated(jac_kextci)) deallocate(jac_kextci)
    if (allocated(jac_backci)) deallocate(jac_backci)
    if (allocated(jac_kextha)) deallocate(jac_kextha)
    if (allocated(jac_backha)) deallocate(jac_backha)
    if (allocated(jac_relhum_lev)) deallocate(jac_relhum_lev)
    if (allocated(jac_temp_lev)) deallocate(jac_temp_lev)
    if (allocated(jac_cwc_q)) deallocate(jac_cwc_q)
    if (allocated(jac_iwc_q)) deallocate(jac_iwc_q)
    if (allocated(jac_rwc_q)) deallocate(jac_rwc_q)
    if (allocated(jac_swc_q)) deallocate(jac_swc_q)
    if (allocated(jac_gwc_q)) deallocate(jac_gwc_q)
    if (allocated(jac_hwc_q)) deallocate(jac_hwc_q)
    if (allocated(jac_cwc_n)) deallocate(jac_cwc_n)
    if (allocated(jac_iwc_n)) deallocate(jac_iwc_n)
    if (allocated(jac_rwc_n)) deallocate(jac_rwc_n)
    if (allocated(jac_swc_n)) deallocate(jac_swc_n)
    if (allocated(jac_gwc_n)) deallocate(jac_gwc_n)
    if (allocated(jac_hwc_n)) deallocate(jac_hwc_n)

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

  end subroutine deallocate_jacobian_vars
end module vars_jacobian