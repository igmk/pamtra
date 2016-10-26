module vars_rt

  use kinds
  use report_module
  use settings, only : radar_npol
  implicit none

  real(kind=dbl), allocatable, dimension(:,:) :: rt_sfc_emissivity, rt_sfc_reflectivity
  real(kind=dbl), allocatable, dimension(:) :: rt_kextatmo, rt_kexttot
  real(kind=dbl), allocatable, dimension(:,:) :: rt_back
  real(kind=dbl), allocatable, dimension(:,:,:,:,:,:) :: rt_scattermatrix_reverse,rt_scattermatrix
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: rt_extmatrix_reverse,rt_extmatrix
  real(kind=dbl), allocatable, dimension(:,:,:,:) :: rt_emisvec_reverse,rt_emisvec

  logical, allocatable, dimension(:) :: rt_hydros_present, rt_hydros_present_reverse

  contains

  subroutine allocate_rt_vars(errorstatus)
    

    use settings, only: nstokes, nummu
    use vars_atmosphere, only: atmo_nlyrs
    use vars_index, only: i_x, i_y

    integer(kind=long) :: nlyr

    integer(kind=long) :: alloc_status
    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=200) :: msg
    character(len=30) :: nameOfRoutine = 'allocate_rt_vars'

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    nlyr = atmo_nlyrs(i_x,i_y)

    call assert_true(err,(nlyr>0),&
        "nlyr must be greater zero")   
    call assert_true(err,(nstokes>0),&
        "nstokes must be greater zero")   
    call assert_true(err,(nummu>0),&
        "nummu must be greater zero")   
    if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
    end if  

    
    allocate(rt_sfc_emissivity(nstokes,nummu), stat=alloc_status)
    allocate(rt_sfc_reflectivity(nstokes,nummu), stat=alloc_status)
    allocate(rt_kextatmo(nlyr), stat=alloc_status)
    allocate(rt_kexttot(nlyr), stat=alloc_status)
    allocate(rt_back(nlyr, radar_npol), stat=alloc_status)
    allocate(rt_scattermatrix_reverse(nlyr,nstokes,nummu,nstokes,nummu,4),stat=alloc_status)
    allocate(rt_scattermatrix(nlyr,nstokes,nummu,nstokes,nummu,4),stat=alloc_status)
    allocate(rt_extmatrix_reverse(nlyr,nstokes,nstokes,nummu,2),stat=alloc_status)
    allocate(rt_extmatrix(nlyr,nstokes,nstokes,nummu,2),stat=alloc_status)
    allocate(rt_emisvec_reverse(nlyr,nstokes,nummu,2),stat=alloc_status)
    allocate(rt_emisvec(nlyr,nstokes,nummu,2),stat=alloc_status)
    allocate(rt_hydros_present(nlyr),stat=alloc_status)
    allocate(rt_hydros_present_reverse(nlyr),stat=alloc_status)

    rt_sfc_emissivity(:,:) = 0._dbl
    rt_sfc_reflectivity(:,:) = 0._dbl
    
    ! set them to zero, just in case they are not calculated but used for Ze/PIA calculation
    rt_kexttot(:) = 0._dbl
    rt_kextatmo(:) = 0._dbl
    rt_back(:,:) = 0._dbl  

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

    errorstatus = err

  end subroutine allocate_rt_vars

  subroutine deallocate_rt_vars()
    if (allocated(rt_sfc_emissivity)) deallocate(rt_sfc_emissivity)
    if (allocated(rt_sfc_reflectivity)) deallocate(rt_sfc_reflectivity)  
    if (allocated(rt_kextatmo)) deallocate(rt_kextatmo)
    if (allocated(rt_kexttot)) deallocate(rt_kexttot)
    if (allocated(rt_back)) deallocate(rt_back)

    if (allocated(rt_scattermatrix_reverse)) deallocate(rt_scattermatrix_reverse)
    if (allocated(rt_scattermatrix)) deallocate(rt_scattermatrix)
    if (allocated(rt_extmatrix_reverse)) deallocate(rt_extmatrix_reverse)
    if (allocated(rt_extmatrix)) deallocate(rt_extmatrix)
    if (allocated(rt_emisvec_reverse)) deallocate(rt_emisvec_reverse)
    if (allocated(rt_emisvec)) deallocate(rt_emisvec)

    if (allocated(rt_hydros_present_reverse)) deallocate(rt_hydros_present_reverse)
    if (allocated(rt_hydros_present)) deallocate(rt_hydros_present)

  end subroutine deallocate_rt_vars

end module vars_rt