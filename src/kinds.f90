module kinds
  !
  !     *** Define usual kinds for strong typing ***
  !
  implicit none
  save
  !
  !     Integer Kinds
  !     -------------
  !
  integer, parameter :: short = selected_int_kind(3)
  integer, parameter :: long  = selected_int_kind(9)
  !
  !     Real Kinds
  !     ----------
  !
  integer, parameter :: sgl = selected_real_kind(6,37)
  integer, parameter :: dbl = selected_real_kind(13,200)
  integer, parameter :: ext = selected_real_kind(25)
  !
end module kinds
