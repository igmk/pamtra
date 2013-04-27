subroutine get_scat_coefs(nz,DELTAM, NUMMU, NLEG, COEF, EXTINCTION, SCATTER)
  !        get_scat_coefs returns
  !      the degree of the Legendre series (NLEGEN), the extinction
  !      coefficient, the single scatter albedo, and the Legendre
  !      coefficients.  Only coefficients for the six unique elements
  !      (for randomly oriented particles with a plane of symmetry)
  !      are returned.  If delta-M scaling is desired, the extinction,
  !      single scattering albedo, and diagonal elements of the phase
  !      matrix are scaled according to the M'th Legendre coefficient
  !      (M=NUMMU).
  use kinds
  use vars_atmosphere

  implicit none

  integer :: nz
  INTEGER NLEG, NUMMU
  REAL(kind=dbl) COEF (6, rt3nlegen(nz)+1 ), EXTINCTION, SCATTER
  CHARACTER(1) DELTAM
  INTEGER L
  REAL(kind=dbl) ALBEDO, F

  extinction = rt3kexttot(nz)
  albedo = rt3salbtot(nz)
  scatter = extinction*albedo
  nleg = rt3nlegen(nz)

  DO L = 1, NLEG + 1
     coef(1,l) = rt3legen(nz,l)
     coef(2,l) = rt3legen2(nz,l)
     coef(3,l) = rt3legen3(nz,l)
     coef(4,l) = rt3legen4(nz,l)
     coef(5,l) = rt3legen(nz,l)
     coef(6,l) = rt3legen3(nz,l)
  ENDDO

  IF (DELTAM.EQ.'Y') THEN
     IF (NUMMU + 1.LE.NLEG + 1) THEN
        F = COEF (1, NUMMU + 1) / (2 * NUMMU + 1)
     ELSE
        F = 0.0
     ENDIF
     ALBEDO = SCATTER / EXTINCTION
     EXTINCTION = (1 - ALBEDO * F) * EXTINCTION
     ALBEDO = (1 - F) * ALBEDO / (1 - ALBEDO * F)
     SCATTER = ALBEDO * EXTINCTION
     NLEGEN = NUMMU - 1
     !  Scale only the diagonal phase matrix elements
     DO L = 0, NLEG
        COEF (1, L + 1) = (2 * L + 1) * MAX (0.0d0, (COEF (1, L + 1) / (2 * L + 1) - F) / (1 - F) )
        COEF (3, L + 1) = (2 * L + 1) * MAX (0.0d0, (COEF (3, L + 1) / (2 * L + 1) - F) / (1 - F) )
        COEF (5, L + 1) = (2 * L + 1) * MAX (0.0d0, (COEF (5, L + 1) / (2 * L + 1) - F) / (1 - F) )
        COEF (6, L + 1) = (2 * L + 1) * MAX (0.0d0, (COEF (6, L + 1) / (2 * L + 1) - F) / (1 - F) )
     ENDDO
  ENDIF
  RETURN
END SUBROUTINE get_scat_coefs
