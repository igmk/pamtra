module nan
contains
  function nan_dbl()

    use kinds, only: dbl

    implicit none

    real(kind=dbl) :: b, nan_dbl

    b = -1.
    nan_dbl = set_nan(b)

   contains

    function set_nan(a)

      implicit none

      real(kind=dbl) :: a, set_nan

      set_nan = sqrt(a)

      return
    end function set_nan

  end function nan_dbl
end module nan