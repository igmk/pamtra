module physical_constants

use kinds
implicit none
save

real(kind=dbl), parameter :: speed_of_light = 299792458 ! m/s
real(kind=dbl), parameter :: pi = 3.141592653589793

complex(kind=dbl), parameter :: Im = (0.0, 1.0)

end module physical_constants