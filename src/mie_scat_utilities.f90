module mie_scat_utilities
  use kinds
  use constants, only: pi, c
  use report_module

  contains

  subroutine miecalc (errorstatus,nterms, x, mn, a, b) 
      !  miecalc calculates the complex mie coefficients an and bn
      !  given the dimensionless size parameter x and the complex
      !  index of refraction (mre,mim).  the number of terms calculated
      !  is given by nterms unless nterms <= 0 or in which case the
      !  appropriate number is calculated and returned in nterms.


      implicit none

      integer :: nterms
      real(kind=dbl) :: x
      complex(kind=dbl) :: mn, a(*), b(*)
      integer :: nstop, n, nn
      integer, parameter :: maxterms = 10000
      real(kind=dbl) :: psin, psim, chin, chim, tmp
      complex(kind=dbl) :: m, y, d(maxterms + 15), xin, xim, ctmp


      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=14) :: nameOfRoutine = 'miecalc'
        
      
      !  if nterms is not specified calculate it
      nstop = int(x + 4.0 * x**0.3334 + 2)
      if (nterms .le. 0) nterms = nstop
      if (nterms .gt. maxterms) then
        errorstatus = fatal
        msg = 'mie calculation requires more terms than available.'
        call report(errorstatus, msg, nameOfRoutine)
        return
      end if

      !  generate the dn's by down recurrence  d = d(log(psi(y)))/dy
      m = conjg(mn)
      y = m * x
      nn = nterms + 15
      d(nn) = cmplx(0.0d0, 0.0d0)
      do n = nn, 2, - 1
        d(n - 1) = n / y - 1.0 / (d(n) + n / y)
      end do

      !  generate the psin's and xin's by upward recurrence
      !  and calculate the an's and bn's from them.
      !  (psin = psi(n), psim = psi(n-1), same for chi)
      psim = cos(x)
      psin = sin(x)
      chim = - sin(x)
      chin = cos(x)
      do n = 1, nterms
        tmp = psin
        psin = (2 * n - 1) / x * psin - psim
        psim = tmp
        tmp = chin
        chin = (2 * n - 1) / x * chin - chim
        chim = tmp
        xin = cmplx(psin, - chin)
        xim = cmplx(psim, - chim)
        ctmp = d(n) / m + n / x
        a(n) = (ctmp * psin - psim) / (ctmp * xin - xim)
        ctmp = m * d (n) + n / x
        b(n) = (ctmp * psin - psim) / (ctmp * xin - xim)
      end do
      errorstatus = err    
      return
  end subroutine miecalc

  subroutine miecross (nterms, x, a, b, qext, qscat, qbackscat)
      !
      !  miecross calculates the extinction, scattering, and
      !  backscatter efficiencies given the mie coefficients an and bn
      !  and the size parameter x.
      !
      ! Output:
      !
      !    qext          extinction efficiency
      !    qscat         scattering efficiency
      !    qbackscat     backscattering efficiency
      !
      ! Note: qbackscat*pi*r2 = back scattering crossection INCLUDING the
      ! akward 4*pi, i.e. identical to sigma_b of eq 4.82 in Bohren & Huffmann
      ! (Max, 10/12)

      implicit none

      integer :: nterms
      real(kind=dbl) :: x, qext, qscat, qbackscat
      complex(kind=dbl) :: a(*), b(*), sum3
      integer :: n
      real(kind=dbl) :: sum1, sum2
       
      
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      do n = 1, nterms
        sum1 = sum1 + (2 * n + 1) * (real(a(n)) + real(b(n)))
        sum2 = sum2 + (2 * n + 1) * (real(a(n) * conjg(a(n)))     &
        + real(b(n) * conjg(b(n))))
        sum3 = sum3 + (2 * n + 1) * ( - 1.0) **n *(a(n) - b(n))
      end do
      qext = 2.0d0 / x**2 * sum1
      qscat = 2.0d0 / x**2 * sum2
      qbackscat = 1.0d0 / x**2 * sum3 * conjg(sum3)

      return
  end subroutine miecross


  subroutine mieangle (nterms, a, b, mu, p1, p2, p3, p4) 
      !  mieangle calculates the intensity scattering matrix elements
      !  (p1,p2,p3,p4) for a particular value of mu (cos(theta)) from the
      !  mie coefficients an's and bn's. the matrix elements are for the
      !  stokes intensity vector (i,q,u,v) and are calculated from the
      !  complex scattering amplitudes s1 and s2.

      implicit none
      integer :: nterms
      real(kind=dbl) :: mu, p1, p2, p3, p4
      complex(kind=dbl) :: a(*), b(*)
      integer :: n
      real(kind=dbl) :: tmp, pin, pim, taun, c
      complex(kind=dbl) :: s1, s2


      s1 = cmplx(0.0, 0.0)
      s2 = cmplx(0.0, 0.0)
      ! sum up the series using the an's and bn's
      pin = 1.0d0
      pim = 0.0d0
      do n = 1, nterms
        taun = n * mu * pin - (n + 1) * pim
        ! calculate the scattering functions at +mu and -mu
        ! using the pin's and the taun's.
        c = (2 * n + 1) / real(n * (n + 1) )
        s1 = s1 + c * (a(n) * pin + b(n) * taun)
        s2 = s2 + c * (b(n) * pin + a(n) * taun)
        ! calculate the angular function pin by up recurrence
        tmp = pin
        pin = ((2 * n + 1) * mu * pin - (n + 1) * pim) / n
        pim = tmp
      end do
      !print*,"s1 ",s1,"    s2 ",s2
      ! calculate the first stokes parameter scattering matrix element
      p1 = 0.5 * (abs(s2) **2 + abs(s1) **2) ! =M11 = M22
      p2 = 0.5 * (abs(s2) **2 - abs(s1) **2)  !=M12 = M21
      p3 = real(conjg(s2) * s1) !=M33 = M44
      p4 = imag(conjg(s1) * s2) !=M34 =-M43
      ! looks like S2=F_vv and S1=F_hh according to Xinxin's Thesis (p.18/19) [Max]
      return
  end subroutine mieangle

  subroutine mieangle_amplScatMat (nterms, a, b, mu, s1, s2) 
      !  mieangle calculates the intensity scattering matrix elements
      !  (p1,p2,p3,p4) for a particular value of mu (cos(theta)) from the
      !  mie coefficients an's and bn's. the matrix elements are for the
      !  stokes intensity vector (i,q,u,v) and are calculated from the
      !  complex scattering amplitudes s1 and s2.


      implicit none
      integer :: nterms
      real(kind=dbl) :: mu
      complex(kind=dbl) :: a(*), b(*)
      integer :: n
      real(kind=dbl) :: tmp, pin, pim, taun, c
      complex(kind=dbl) :: s1, s2


      s1 = cmplx(0.0, 0.0)
      s2 = cmplx(0.0, 0.0)
      ! sum up the series using the an's and bn's
      pin = 1.0d0
      pim = 0.0d0
      do n = 1, nterms
        taun = n * mu * pin - (n + 1) * pim
        ! calculate the scattering functions at +mu and -mu
        ! using the pin's and the taun's.
        c = (2 * n + 1) / real(n * (n + 1) )
        s1 = s1 + c * (a(n) * pin + b(n) * taun)
        s2 = s2 + c * (b(n) * pin + a(n) * taun)
        ! calculate the angular function pin by up recurrence
        tmp = pin
        pin = ((2 * n + 1) * mu * pin - (n + 1) * pim) / n
        pim = tmp
      end do
      !   ! calculate the first stokes parameter scattering matrix element
      !   p1 = 0.5 * (abs(s2) **2 + abs(s1) **2) ! =M11 = M22
      !   p2 = 0.5 * (abs(s2) **2 - abs(s1) **2)  !=M12 = M21
      !   p3 = real(conjg(s2) * s1) !=M33 = M44
      !   p4 = imag(conjg(s1) * s2) !=M34 =-M43
      ! looks like S2=F_vv and S1=F_hh according to Xinxin's Thesis (p.18/19) [Max]
      return
  end subroutine mieangle_amplScatMat



  subroutine amplScatMat_to_scatMat(s1,s2,scatMat)
      ! convertes the amplitude scattering matrix elements S1 and S2 to
      ! scattering or Mueller Matrix!

      implicit none

      integer, parameter :: nstokes = 2

      complex(kind=dbl),intent(in) :: s1, s2
      real(kind=dbl) :: p1,p2!,p3,p4
      real(kind=dbl),intent(out), dimension(nstokes,nstokes) :: scatMat


      p1 = 0.5d0 * (abs(s2) **2 + abs(s1) **2)
      p2 = 0.5d0 * (abs(s2) **2 - abs(s1) **2)
      !p3 = real(conjg(s2) * s1)
      !p4 = imag(conjg(s1) * s2)

      scatMat(:,:) = 0d0

      !see Evans and Stephans 1991, eq6
      scatMat(1,1) = p1
      scatMat(2,2) = p1 !for spheres
      scatMat(1,2) = p2
      scatMat(2,1) = p2
    
      !   scatMat(3,3) = p3
      !   scatMat(4,4) = p3 !for spheres
      !   scatMat(3,4) = p4 !check plus or minus!
      !   scatMat(4,3) = -p4 !check plus or minus!
      !
      return
  end subroutine amplScatMat_to_scatMat


  ! subroutine amplScatMat_to_emissionVector(s1,s2,emissionVector)
  ! 
  !   use kinds
  ! 
  !   implicit none 
  ! 
  !   integer, parameter :: nstokes = 2
  ! 
  ! 
  !   return
  ! end subroutine amplScatMat_to_emissionVector

  subroutine amplScatMat_to_extinctionMatrix(s1,s2,f,extinctionMatrix)

      implicit none

      integer, parameter :: nstokes = 2


      real(kind=dbl), intent(in) :: f
      complex(kind=dbl) :: s1,s2
      real(kind=dbl),intent(out), dimension(nstokes,nstokes) :: extinctionMatrix
      real(kind=dbl) :: k


      k = 2* pi*f*1d9/c

      extinctionMatrix(:,:) = 0d0
      !extinctionMatrix(row,colum)
      extinctionMatrix(1,1) = imag(s1+s2)*2*pi/k
      extinctionMatrix(1,2) = imag(s1-s2)*2*pi/k
      extinctionMatrix(2,2) = imag(s1+s2)*2*pi/k
      extinctionMatrix(2,1) = imag(s1-s2)*2*pi/k

      return
  end subroutine amplScatMat_to_extinctionMatrix

end module mie_scat_utilities
