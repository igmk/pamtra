module vars_hydroFullSpec

  use kinds
  use report_module
  implicit none

  save
  
  integer(kind=long) :: hydrofs_nbins
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_rho_ds
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_d_ds
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_d_bound_ds
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_n_ds
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_mass_ds
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_area_ds
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_canting
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_as_ratio
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_fallvelocity
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_rg_kappa_ds
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_rg_beta_ds
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_rg_gamma_ds
  real(kind=dbl), allocatable, dimension(:,:,:,:,:) :: hydrofs_rg_zeta_ds

  contains

  subroutine allocate_hydrofs_vars(errorstatus,nbins)
             
    use vars_atmosphere, only: atmo_ngridx,atmo_ngridy, atmo_max_nlyrs
    use descriptor_file, only: n_hydro

    integer, intent(in) :: nbins
    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=40) :: nameOfRoutine = 'allocate_hydrofs_vars' 

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    hydrofs_nbins = nbins

    call assert_true(err,(atmo_ngridx>0),&
        "atmo_ngridx must be greater zero")   
    call assert_true(err,(atmo_ngridy>0),&
        "atmo_ngridy must be greater zero")   
    call assert_true(err,(atmo_max_nlyrs>0),&
        "atmo_max_nlyrs must be greater zero")
    call assert_true(err,(hydrofs_nbins>0),&
        "hydrofs_nbins must be greater zero")   
    call assert_true(err,(n_hydro>0),&
        "n_hydro must be greater zero")
    if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
    end if  


    allocate(hydrofs_rho_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_d_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_d_bound_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins+1))
    allocate(hydrofs_n_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_mass_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_area_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_canting(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_as_ratio(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_fallvelocity(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_rg_kappa_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_rg_beta_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_rg_gamma_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))
    allocate(hydrofs_rg_zeta_ds(atmo_ngridx,atmo_ngridy, atmo_max_nlyrs,n_hydro,hydrofs_nbins))


    errorstatus = err
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    return
  end subroutine allocate_hydrofs_vars

  subroutine deallocate_hydrofs_vars

    implicit none
    character(len=40) :: nameOfRoutine = 'deallocate_hydrofs_vars'

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    if (allocated(hydrofs_rho_ds)) deallocate(hydrofs_rho_ds)
    if (allocated(hydrofs_d_ds)) deallocate(hydrofs_d_ds)
    if (allocated(hydrofs_d_bound_ds)) deallocate(hydrofs_d_bound_ds)
    if (allocated(hydrofs_n_ds)) deallocate(hydrofs_n_ds)
    if (allocated(hydrofs_mass_ds)) deallocate(hydrofs_mass_ds)
    if (allocated(hydrofs_area_ds)) deallocate(hydrofs_area_ds)
    if (allocated(hydrofs_canting)) deallocate(hydrofs_canting)
    if (allocated(hydrofs_as_ratio)) deallocate(hydrofs_as_ratio)
    if (allocated(hydrofs_fallvelocity)) deallocate(hydrofs_fallvelocity)
    if (allocated(hydrofs_rg_kappa_ds)) deallocate(hydrofs_rg_kappa_ds)
    if (allocated(hydrofs_rg_beta_ds)) deallocate(hydrofs_rg_beta_ds)
    if (allocated(hydrofs_rg_gamma_ds)) deallocate(hydrofs_rg_gamma_ds)
    if (allocated(hydrofs_rg_zeta_ds)) deallocate(hydrofs_rg_zeta_ds)
 
    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)


  end subroutine deallocate_hydrofs_vars

  subroutine print_hydrofs_vars()

    use vars_atmosphere, only: atmo_ngridx,atmo_ngridy, atmo_max_nlyrs
    use descriptor_file, only: n_hydro
    integer(kind=long) :: xx,yy,zz,hh

    do xx=1,atmo_ngridx
      do yy=1,atmo_ngridy
        do zz=1,atmo_max_nlyrs
          do hh=1,n_hydro
            print*, "xx", xx, "yy", yy, "zz", zz, "hh", hh
            print*, "hydrofs_rho_ds", hydrofs_rho_ds(xx,yy,zz,hh, :)
            print*, "hydrofs_d_ds", hydrofs_d_ds(xx,yy,zz,hh, :)
            print*, "hydrofs_d_bound_ds", hydrofs_d_bound_ds(xx,yy,zz,hh, :)
            print*, "hydrofs_n_ds", hydrofs_n_ds(xx,yy,zz,hh, :)
            print*, "hydrofs_mass_ds", hydrofs_mass_ds(xx,yy,zz,hh, :)
            print*, "hydrofs_area_ds", hydrofs_area_ds(xx,yy,zz,hh, :)
            print*, "hydrofs_canting", hydrofs_canting(xx,yy,zz,hh, :)
            print*, "hydrofs_as_ratio", hydrofs_as_ratio(xx,yy,zz,hh, :)
            print*, "hydrofs_fallvelocity", hydrofs_fallvelocity(xx,yy,zz,hh, :)
            print*, "hydrofs_rg_kappa_ds", hydrofs_rg_kappa_ds(xx,yy,zz,hh, :)
            print*, "hydrofs_rg_beta_ds", hydrofs_rg_beta_ds(xx,yy,zz,hh, :)
            print*, "hydrofs_rg_zeta_ds", hydrofs_rg_zeta_ds(xx,yy,zz,hh, :)
            print*, "hydrofs_rg_gamma_ds", hydrofs_rg_gamma_ds(xx,yy,zz,hh, :)

          end do
        end do
      end do
    end do

  end subroutine print_hydrofs_vars


end module vars_hydroFullSpec

