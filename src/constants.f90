module constants
  ! Description:
  ! Definition of all constants for pamtra
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  0.1   17/11/2009    creation of file 
  use kinds
 
  implicit none

  real(kind=dbl), parameter :: pi      = 3.141592653589793

  complex(kind=dbl), parameter :: im = (0.0, 1.0)

end module constants