module libRefIndex

  contains
  
  subroutine get_refIndex_water(s,T,f,model,refre,refim)
  use settings, only: liq_mod
  implicit none

  double precision, intent(in) :: s,& ! salinity [0/00]
       T,& ! temperature [°C]
       f   ! frequency [GHz]
  character(3), intent(in) :: model

  double precision, intent(out) :: refre,& ! real part of complex index of refraction n_r []
       refim ! imaginary part of complex index of refraction n_i []

  double precision::   absind,&! absind absorptive index (n_i/n_r) []
       abscoef ! abscoef absorption coefficient (4*pi*n_i*f/c) [1/m]
  integer :: err

  !input
  !f2py intent(in) :: s,T,f, model
  !meta out
  !f2py intent(out) :: refre,refim

  liq_mod = model
  call ref_water(err,s,T,f,refre,refim,absind,abscoef)

  end subroutine get_refIndex_water

  subroutine get_dielecConst_water(s,T,f,K2)

  implicit none

  double precision, intent(in) :: s,& ! salinity [0/00]
       T,& ! temperature [°C]
       f   ! frequency [GHz]

  double precision, intent(out) :: K2  ! dielectric constant |K|^2
 
  double precision :: dielec_water
  !input
  !f2py intent(in) :: s,T,f
  !meta out
  !f2py intent(out) :: K2

  K2 = dielec_water(s,T,f)

  end subroutine get_dielecConst_water

end module libRefIndex