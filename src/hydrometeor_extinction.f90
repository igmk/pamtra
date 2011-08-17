subroutine hydrometeor_extinction(f,frq_str)

  use kinds
  use vars_atmosphere
  use nml_params, only: verbose, tmp_path, active, passive, dump_to_file, n_moments
  use constants
  use mod_io_strings

  implicit none

  integer, parameter :: maxleg = 200

  integer :: jj, nz

  integer :: nlegencw, nlegenci, nlegenrr, nlegensn, nlegengr, nlegenha

  real(kind=dbl) :: f, wavelength

  real(kind=dbl) :: kextcw, salbcw, kextrr, salbrr,  &
       kextci, salbci, kextsn, salbsn, kextgr, salbgr,  kextha, salbha, &
       backcw, backrr, backci, backsn, backgr, backha

  real(kind=dbl), dimension(200) ::  LEGENcw, LEGENrr, LEGENci, LEGENgr, LEGENsn, LEGENha,      &
       LEGEN2cw, LEGEN2rr, LEGEN2ci, LEGEN2gr, LEGEN2sn, LEGEN2ha,  &
       LEGEN3cw, LEGEN3rr, LEGEN3ci, LEGEN3gr, LEGEN3sn, LEGEN3ha,  &
       LEGEN4cw, LEGEN4rr, LEGEN4ci, LEGEN4gr, LEGEN4sn, LEGEN4ha


  real(kind=dbl), dimension(2) :: P11, ang

  real(kind=dbl) :: threshold ! threshold value for hydrometeor extinction as mass mixing ratio

  character(6), intent(in) :: frq_str !from commandline

  if (verbose .gt. 1) print*, 'Entering hydrometeor_extinction'

  if (dump_to_file) file_ph = ''
  ! INITIALIZATION OF LEGENDRE COEFFICIENTS

  nlegen = 0
  legen   = 0.d0
  legen2  = 0.d0
  legen3  = 0.d0
  legen4  = 0.d0

  threshold = 1.d-5   ! [kg/kg]

  if (verbose .gt. 1) print*, 'start loop over layer'

  grid_z: do nz = 1, nlyr  ! loop over all layers

     if (verbose .gt. 1) print*, 'Layer: ', nz

     !---------------------------salinity------------------------------
     ! calculation of the single scattering properties
     ! of hydrometeors. cloud water and cloud ice are 
     ! with respect to radius. whereas the distribution 
     ! of precipitating particles is with respect to diameter.
     !---------------------------------------------------------

     !---------------------------------------------------------
     !        single scattering properties of cloud water
     !---------------------------------------------------------

     nlegencw = 0 
     legencw  = 0.d0
     legen2cw = 0.d0
     legen3cw = 0.d0
     legen4cw = 0.d0

     kextcw = 0.d0 
     salbcw = 0.d0 
     backcw = 0.d0 

     if (cwc_q(nz) .ge. threshold .and. n_moments .eq. 1) then
    	call cloud_ssp(f,cwc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg, kextcw, salbcw, backcw,  &
             nlegencw, legencw, legen2cw, legen3cw, legen4cw)
     else if (cwc_q(nz) .ge. threshold .and. n_moments .eq. 2) then
    	call cloud_ssp(f,cwc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg, kextcw, salbcw, backcw,  &
             nlegencw, legencw, legen2cw, legen3cw, legen4cw, cwc_n(nz))
     else
        kextcw = 0.0d0
        salbcw = 0.0d0
        backcw = 0.0d0
     end if
     !---------------------------------------------------------
     !       single scattering properties of rain
     !---------------------------------------------------------

     nlegenrr = 0 
     legenrr  = 0.d0
     legen2rr = 0.d0
     legen3rr = 0.d0
     legen4rr = 0.d0

     kextrr = 0.d0 
     salbrr = 0.d0 
     backrr = 0.d0 

     if (rwc_q(nz) .ge. threshold .and. n_moments .eq. 1) then
        call rain_ssp(f,rwc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg,kextrr, salbrr, backrr,  &
             nlegenrr, legenrr, legen2rr, legen3rr, legen4rr)
     else if (rwc_q(nz) .ge. threshold .and. n_moments .eq. 2) then
      	call rain_ssp(f,rwc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg,kextrr, salbrr, backrr,  &
             nlegenrr, legenrr, legen2rr, legen3rr, legen4rr, rwc_n(nz))
     else
        kextrr = 0.0d0
        salbrr = 0.0d0
        backrr = 0.0d0
     end if

     !---------------------------------------------------------
     !       single scattering properties of ice crystals
     !---------------------------------------------------------

     nlegenci = 0 
     legenci  = 0.d0
     legen2ci = 0.d0
     legen3ci = 0.d0
     legen4ci = 0.d0

     kextci = 0.0d0 
     salbci = 0.0d0 
     backci = 0.0d0 

     if (iwc_q(nz) .ge. threshold*1.e-2 .and. n_moments .eq. 1) then
        call ice_ssp(f,iwc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg,kextci, salbci, backci,  &
             nlegenci, legenci, legen2ci, legen3ci, legen4ci)
     else if (iwc_q(nz) .ge. threshold .and. n_moments .eq. 2) then
      	call ice_ssp(f,iwc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg,kextci, salbci, backci,  &
             nlegenci, legenci, legen2ci, legen3ci, legen4ci, iwc_n(nz))
     else
        kextci = 0.0d0
        salbci = 0.0d0
        backci = 0.0d0
     end if

     !---------------------------------------------------------
     !       single scattering properties of snow
     !---------------------------------------------------------

     nlegensn = 0 
     legensn = 0.0d0
     legen2sn = 0.0d0
     legen3sn = 0.0d0
     legen4sn = 0.0d0

     if (swc_q(nz) .ge. threshold .and. n_moments .eq. 1) then
        call snow_ssp(f,swc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg,kextsn, salbsn, backsn,  &
             nlegensn, legensn, legen2sn, legen3sn, legen4sn)
        call legendre2phasefunction(legensn, nlegensn, 2, 200,p11, ang)
        backsn = kextsn * salbsn * P11 (2)
     else if (swc_q(nz) .ge. threshold .and. n_moments .eq. 2) then
      	call snow_ssp(f,swc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg,kextsn, salbsn, backsn,  &
             nlegensn, legensn, legen2sn, legen3sn, legen4sn, swc_n(nz))
        call legendre2phasefunction(legensn, nlegensn, 2, 200,p11, ang)
        backsn = kextsn * salbsn * P11 (2)
     else
        kextsn = 0.0d0
        salbsn = 0.0d0
        backsn = 0.0d0
     endif

     !---------------------------------------------------------
     !       single scattering properties of graupel
     !---------------------------------------------------------

     nlegengr = 0 
     legengr = 0.0d0
     legen2gr = 0.0d0
     legen3gr = 0.0d0
     legen4gr = 0.0d0

     if (gwc_q(nz) .ge. threshold .and. n_moments .eq. 1) then
        call grau_ssp(f,gwc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg,kextgr, salbgr, backgr,  &
             nlegengr, legengr, legen2gr, legen3gr, legen4gr)
        call legendre2phasefunction(legengr, nlegengr, 2, 200, p11, ang)
        backgr = kextgr * salbgr * p11 (2)
     else if (gwc_q(nz) .ge. threshold .and. n_moments .eq. 2) then
      	call grau_ssp(f,gwc_q(nz),temp(nz),press(nz),q_hum(nz),&
             maxleg,kextgr, salbgr, backgr,  &
             nlegengr, legengr, legen2gr, legen3gr, legen4gr, gwc_n(nz))
      	call legendre2phasefunction(legengr, nlegengr, 2, 200, p11, ang)
        backgr = kextgr * salbgr * p11 (2)
     else
        kextgr = 0.0d0
        salbgr = 0.0d0
        backgr = 0.0d0
     endif

     !---------------------------------------------------------
     !       single scattering properties of hail
     !---------------------------------------------------------

     nlegenha = 0
     legenha = 0.0d0
     legen2ha = 0.0d0
     legen3ha = 0.0d0
     legen4ha = 0.0d0

     if (n_moments .eq. 2) then
        if (hwc_q(nz) .ge. 1000.) then
           call hail_ssp(f,hwc_q(nz),temp(nz),press(nz),q_hum(nz),&
                maxleg,kextha, salbha, backha,  &
                nlegenha, legenha, legen2ha, legen3ha, legen4ha, hwc_n(nz))
           call legendre2phasefunction(legenha, nlegenha, 2, 200, p11, ang)
           backha = kextha * salbha * p11 (2)
        else
           kextha = 0.0d0
           salbha = 0.0d0
           backha = 0.0d0
        endif
     else
        kextha = 0.0d0
        salbha = 0.0d0
        backha = 0.0d0
     endif

     nlegen(nz) = max(nlegen(nz),nlegencw,nlegenci,nlegenrr,nlegensn,nlegengr,nlegenha)

     if (verbose .gt. 1) print*, 'End of scattering calc for layer: ', nz

     !CCCCCCCCCCCCC   END OF SINGLE SCATTERING PROPERTY COMPUTATIONS  CCCCCCC

     !                                                                       
     !           Summing up the scattering parameters and writing the
     !           input file of the scattering properties of each layer
     !                                                                       


     kexttot(nz) = kextsn + kextcw + kextrr + kextgr + kextci + kextha
     back(nz) = backcw + backrr + backci + backsn + backgr + backha


     if (kexttot(nz) .lt. 0.) write(*,*) 'something wrong'
     if (kexttot(nz) .le. 0.) then 
        salbtot(nz) = 0.0
     else 
        salbtot(nz) = (salbcw * kextcw + salbrr *       &
             kextrr + salbci * kextci + salbsn * kextsn + salbgr *    &
             kextgr + salbha * kextha) / kexttot(nz)
     endif
     !
     !      absorp(nz) = (1.0 - salbtot(nz) ) * kexttot(nz)

!!!!!!!!!!!!!!!!! check whether hgt_lev needs to be km or m !!!!!!!!!!!!!!!!!

     !   summing up the Legendre coefficient                                 

     if (kexttot(nz) .gt. 0.0 .or.    &
          salbtot(nz) .gt. 0.0) then ! there are hydrometeor present : a PH file is needed
        if (dump_to_file) then
           write(nzstr, '(i2.2)') nz
           FILE_PH(nz) = tmp_path(:len_trim(tmp_path))//'/PHx'//xstr//'y'//ystr//'lev'//Nzstr//'f'//frq_str

           open(unit=21, file=file_PH(nz), STATUS='unknown', &
	        form='FORMATTED')
           write(21,*) kexttot(nz), '   EXINCTION'
           write(21,*) kexttot(nz) * salbtot(nz), '   SCATTERING'
           write(21,*) salbtot(nz), '   SINGLE SCATTERING ALBEDO'
           write(21,*) Nlegen(nz) - 1, '      DEGREE OF LEGENDRE SERIES'
        end if

        do jj = 1, Nlegen(nz)
           legen (nz,jj) = (legencw (jj) * salbcw * kextcw + legenrr ( &
                jj) * salbrr * kextrr + legenci (jj) * salbci * kextci + &
                legensn (jj) * salbsn * kextsn + legengr (jj) * salbgr * &
                kextgr + legenha (jj) * salbha * kextha) / (salbtot (nz) &
                * kexttot (nz) )

           legen2 (nz,jj) = (legen2cw (jj) * salbcw * kextcw +         &
                legen2rr (jj) * salbrr * kextrr + legen2ci (jj) * salbci &
                * kextci + legen2sn (jj) * salbsn * kextsn + legen2gr (  &
                jj) * salbgr * kextgr + legen2ha (jj) * salbha * kextha ) &
                / (salbtot (nz) * kexttot (nz) )

           legen3 (nz,jj) = (legen3cw (jj) * salbcw * kextcw +         &
                legen3rr (jj) * salbrr * kextrr + legen3ci (jj) * salbci &
                * kextci + legen3sn (jj) * salbsn * kextsn + legen3gr (  &
                jj) * salbgr * kextgr + legen3ha (jj) * salbha * kextha) &
                / (salbtot (nz) * kexttot (nz) )

           legen4 (nz,jj) = (legen4cw(jj) * salbcw * kextcw +         &
                legen4rr (jj) * salbrr * kextrr + legen4ci (jj) * salbci &
                * kextci + legen4sn (jj) * salbsn * kextsn + legen4gr (  &
                jj) * salbgr * kextgr + legen4ha (jj) * salbha * kextha) &
                / (salbtot(nz) * kexttot &
                (nz))
           if (dump_to_file) then
              write (21, 1005) jj - 1, legen (nz,jj), legen2 (nz,jj),        &
                   legen3(nz,jj), legen4(nz,jj), legen(nz,jj), legen3(nz,jj)
           end if
           g_coeff (nz) = legen (nz,2) / 3.0d0

	end do ! end of cycle over Legendre coefficient
        if (dump_to_file) close(21)
     end if
  end do grid_z !end of cycle over the vertical layers

  if (verbose .gt. 1) print*, 'Exiting hydrometeor_extinction'

1005 format  (i3,6(1x,f10.7))

  return

end subroutine hydrometeor_extinction
