! Unit test for the genuinely stateless eps_ice() function (depends only on
! kinds, no shared module state) -- see CLAUDE.md on why most of src/ is not
! amenable to this kind of test. Run via `meson test`.
program test_eps_ice
  use kinds, only: dbl

  implicit none

  interface
    function eps_ice(T, f)
      use kinds, only: dbl
      real(kind=dbl), intent(in) :: T, f
      complex(kind=dbl) :: eps_ice
    end function eps_ice
  end interface

  complex(kind=dbl) :: got
  real(kind=dbl), parameter :: expected_real = 3.16747_dbl
  real(kind=dbl), parameter :: expected_imag = 2.0919968775736137d-3
  real(kind=dbl), parameter :: tol = 1.0d-6

  got = eps_ice(250.0_dbl, 35.0_dbl)

  if (abs(real(got) - expected_real) > tol .or. abs(aimag(got) - expected_imag) > tol) then
    print *, "eps_ice(250K, 35GHz) mismatch. got=", got, &
             " expected=", expected_real, expected_imag
    stop 1
  end if

  print *, "eps_ice: OK", got
end program test_eps_ice
