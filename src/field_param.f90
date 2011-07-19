subroutine field_param(temp,qs, n0s)

use kinds

implicit none

integer :: nn

real(kind=dbl), intent(in) :: temp, qs
real(kind=dbl), intent(out) :: n0s
real(kind=dbl), dimension(10) :: mma, mmb ! Coeffs for moment relation based on 2nd moment (Field 2005)
real(kind=dbl) :: alf,bet,hlp,zams,m2s,m3s,n0s1,n0s2,sp


n0s1 = 13.5 * 5.65e5 ! parameter in N0S(T)
n0s2 = -0.107        ! parameter in N0S(T), Field et al

zams  = 0.069  ! Formfactor in the mass-size relation of snow particles corresponding to the radius

mma = (/   5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
           0.312550,  0.000204,  0.003199, 0.000000, -0.015952 /)
mmb = (/   0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
           0.060366,  0.000079,  0.000594, 0.000000, -0.003577 /)

  ! Calculate n0s using the temperature-dependent moment
  ! relations of Field et al. (2005)
tm = MAX(MIN(temp-273.15,0.0),-40.0)

nn  = 3
hlp = mma(1)      +mma(2)*tm      +mma(3)*nn       +mma(4)*tm*nn+mma(5)*tm**2 &
    + mma(6)*nn**2+mma(7)*tm**2*nn+mma(8)*tm*nn**2+mma(9)*tm**3+mma(10)*nn**3
alf = 10.0**hlp
bet = mmb(1)      +mmb(2)*tm      +mmb(3)*nn       +mmb(4)*tm*nn+mmb(5)*tm**2 &
    + mmb(6)*nn**2+mmb(7)*tm**2*nn+mmb(8)*tm*nn**2+mmb(9)*tm**3+mmb(10)*nn**3
m2s = qs*zams
m3s = alf*EXP(bet*LOG(m2s))

hlp  = n0s1*EXP(n0s2*tm)
n0s = 13.50 * m2s**4 / m3s**3
n0s = MAX(n0s,0.5*hlp)
n0s = MIN(n0s,1e2*hlp)
n0s = MIN(n0s,1e9)
n0s = MAX(n0s,1e6)

return

end subroutine field_param
