module conversions

contains

real function q2spec(spec_var,t,p,qv,qc,qi,qr,qs,qg,qh)

  use kinds

  implicit none

  real(kind=dbl), parameter :: r_d = 287.05d0,& ! gas constant of dry air
							   r_v = 461.5d0    ! gas constant of water vapor

  real(kind=dbl), intent(in) :: t,&
  								p

  real(kind=dbl), intent(in) :: spec_var,& ! specific variable to convert
  								qv,&
  								qc,&
  								qi,&
  								qr,&
  								qs,&
  								qg

  real(kind=dbl), intent(in), optional :: qh

  if (present(qh)) then
	q2spec = spec_var*p/(r_d*(1.d0+(r_v/r_d-1.d0)*qv-qs-qc-qi-qr-qg-qh)*t)
  else
	q2spec = spec_var*p/(r_d*(1.d0+(r_v/r_d-1.d0)*qv-qs-qc-qi-qr-qg)*t)
  end if


end function q2spec
end module conversions
