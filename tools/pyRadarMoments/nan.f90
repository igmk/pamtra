function nan()

  use kinds, only: dbl

  implicit none

  real(kind=dbl) :: b, nan

  b = -1.
  nan = set_nan(b)

 contains

function set_nan(a)

  implicit none

  real(kind=dbl) :: a, set_nan

  set_nan = sqrt(a)

  return
end function set_nan

end function nan