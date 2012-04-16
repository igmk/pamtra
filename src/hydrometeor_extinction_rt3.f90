subroutine hydrometeor_extinction_rt3(f,frq_str)

  use kinds
  use vars_atmosphere
  use nml_params, only: verbose, tmp_path, active, passive, dump_to_file, n_moments
  use constants
  use mod_io_strings
  use conversions

  implicit none

  integer, parameter :: maxleg = 200

  integer :: jj, nz

  integer :: nlegencw, nlegenci, nlegenrr, nlegensn, nlegengr, nlegenha

  real(kind=dbl) :: f, qwc, nc

  real(kind=dbl) :: salbcw, salbrr, salbci, salbsn, salbgr, salbha


  real(kind=dbl), dimension(200) ::  LEGENcw, LEGENrr, LEGENci, LEGENgr, LEGENsn, LEGENha,      &
       LEGEN2cw, LEGEN2rr, LEGEN2ci, LEGEN2gr, LEGEN2sn, LEGEN2ha,  &
       LEGEN3cw, LEGEN3rr, LEGEN3ci, LEGEN3gr, LEGEN3sn, LEGEN3ha,  &
       LEGEN4cw, LEGEN4rr, LEGEN4ci, LEGEN4gr, LEGEN4sn, LEGEN4ha


  real(kind=dbl), dimension(2) :: P11, ang

  real(kind=dbl) :: threshold ! threshold value for hydrometeor extinction as mass mixing ratio

  character(6), intent(in) :: frq_str !from commandline

  if (verbose .gt. 1) print*, 'Entering hydrometeor_extinction_rt3'

  if (dump_to_file) file_ph = ''
  ! INITIALIZATION OF LEGENDRE COEFFICIENTS

  nlegen = 0
  legen   = 0.d0
  legen2  = 0.d0
  legen3  = 0.d0
  legen4  = 0.d0

  threshold = 1.d-10   ! [kg/kg]

  if (verbose .gt. 1) print*, 'start loop over layer'

  grid_z: do nz = 1, nlyr  ! loop over all layers
! print *,temp(nz)
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

     kextcw(nz) = 0.d0 
     salbcw = 0.d0 
     backcw(nz) = 0.d0 

     if (cwc_q(nz) .ge. threshold) then
     	if (n_moments .eq. 1) then
	     	qwc = q2abs(cwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
    		call cloud_ssp(f,qwc,temp(nz),&
             	maxleg, kextcw(nz), salbcw, backcw(nz),  &
             	nlegencw, legencw, legen2cw, legen3cw, legen4cw)
     	else if (n_moments .eq. 2) then
	     	qwc = q2abs(cwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	        nc = q2abs(cwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
    		call cloud_ssp(f,qwc,temp(nz),&
             	maxleg, kextcw(nz), salbcw, backcw(nz),  &
             	nlegencw, legencw, legen2cw, legen3cw, legen4cw, nc)
       	end if
     else
        kextcw(nz) = 0.0d0
        salbcw = 0.0d0
        backcw(nz) = 0.0d0
     end if


     !---------------------------------------------------------
     !       single scattering properties of ice crystals
     !---------------------------------------------------------

     nlegenci = 0 
     legenci  = 0.d0
     legen2ci = 0.d0
     legen3ci = 0.d0
     legen4ci = 0.d0

     kextci(nz) = 0.0d0 
     salbci = 0.0d0 
     backci(nz) = 0.0d0 

     if (iwc_q(nz) .ge. threshold) then
	     if (n_moments .eq. 1) then
	     	qwc = q2abs(iwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
	        call ice_ssp(f,qwc,temp(nz),&
	             maxleg,kextci(nz), salbci, backci(nz),  &
	             nlegenci, legenci, legen2ci, legen3ci, legen4ci)
	     else if (n_moments .eq. 2) then
	     	qwc = q2abs(iwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	        nc = q2abs(iwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	      	call ice_ssp(f,qwc,temp(nz),&
	             maxleg,kextci(nz), salbci, backci(nz),  &
	             nlegenci, legenci, legen2ci, legen3ci, legen4ci, nc)
	     end if

     else
        kextci(nz) = 0.0d0
        salbci = 0.0d0
        backci(nz) = 0.0d0
     end if

     !---------------------------------------------------------
     !       single scattering properties of rain
     !---------------------------------------------------------

     nlegenrr = 0
     legenrr  = 0.d0
     legen2rr = 0.d0
     legen3rr = 0.d0
     legen4rr = 0.d0

     kextrr(nz) = 0.d0
     salbrr = 0.d0
     backrr(nz) = 0.d0

     if (rwc_q(nz) .ge. threshold) then
     	if (n_moments .eq. 1) then
	       qwc = q2abs(rwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
	       call rain_ssp(f,qwc,temp(nz),&
	            maxleg,kextrr(nz), salbrr, backrr(nz),  &
	            nlegenrr, legenrr, legen2rr, legen3rr, legen4rr)
	    else if (n_moments .eq. 2) then
	     	qwc = q2abs(rwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	        nc = q2abs(rwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	     	call rain_ssp(f,qwc,temp(nz),&
	            maxleg,kextrr(nz), salbrr, backrr(nz),  &
	            nlegenrr, legenrr, legen2rr, legen3rr, legen4rr, nc)
	    end if
     else
        kextrr(nz) = 0.0d0
        salbrr = 0.0d0
        backrr(nz) = 0.0d0
     end if

     !---------------------------------------------------------
     !       single scattering properties of snow
     !---------------------------------------------------------

     nlegensn = 0 
     legensn = 0.0d0
     legen2sn = 0.0d0
     legen3sn = 0.0d0
     legen4sn = 0.0d0

     if (swc_q(nz) .ge. threshold) then
     	 if (n_moments .eq. 1) then
	     	qwc = q2abs(swc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
	        call snow_ssp(f,qwc,temp(nz),&
	             maxleg,kextsn(nz), salbsn, backsn(nz),  &
	             nlegensn, legensn, legen2sn, legen3sn, legen4sn)
	        call legendre2phasefunction(legensn, nlegensn, 2, 200,p11, ang)
	        backsn(nz) = kextsn(nz) * salbsn * P11 (2)
	     else if (n_moments .eq. 2) then
	     	qwc = q2abs(swc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	        nc = q2abs(swc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	      	call snow_ssp(f,qwc,temp(nz),&
	             maxleg,kextsn(nz), salbsn, backsn(nz),  &
	             nlegensn, legensn, legen2sn, legen3sn, legen4sn, nc)
	        call legendre2phasefunction(legensn, nlegensn, 2, 200,p11, ang)
	        backsn(nz) = kextsn(nz) * salbsn * P11 (2)
	     end if
     else
        kextsn(nz) = 0.0d0
        salbsn = 0.0d0
        backsn(nz) = 0.0d0
     endif

     !---------------------------------------------------------
     !       single scattering properties of graupel
     !---------------------------------------------------------

     nlegengr = 0 
     legengr = 0.0d0
     legen2gr = 0.0d0
     legen3gr = 0.0d0
     legen4gr = 0.0d0

     if (gwc_q(nz) .ge. threshold) then
	     if (n_moments .eq. 1) then
	     	qwc = q2abs(gwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz))
	        call grau_ssp(f,qwc,temp(nz),&
	             maxleg,kextgr(nz), salbgr, backgr(nz),  &
	             nlegengr, legengr, legen2gr, legen3gr, legen4gr)
	        call legendre2phasefunction(legengr, nlegengr, 2, 200, p11, ang)
	        backgr(nz) = kextgr(nz) * salbgr * p11 (2)
	     else if (n_moments .eq. 2) then
	     	qwc = q2abs(gwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	        nc = q2abs(gwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	      	call grau_ssp(f,qwc,temp(nz),&
	             maxleg,kextgr(nz), salbgr, backgr(nz),  &
	             nlegengr, legengr, legen2gr, legen3gr, legen4gr, nc)
	      	call legendre2phasefunction(legengr, nlegengr, 2, 200, p11, ang)
	        backgr(nz) = kextgr(nz) * salbgr * p11 (2)
	     end if
     else
        kextgr(nz) = 0.0d0
        salbgr = 0.0d0
        backgr(nz) = 0.0d0
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
        if (hwc_q(nz) .ge. threshold) then
	       qwc = q2abs(hwc_q(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
	       nc = q2abs(hwc_n(nz),temp(nz),press(nz),q_hum(nz),cwc_q(nz),iwc_q(nz),rwc_q(nz),swc_q(nz),gwc_q(nz),hwc_q(nz))
           call hail_ssp(f,qwc,temp(nz),&
                maxleg,kextha(nz), salbha, backha(nz),  &
                nlegenha, legenha, legen2ha, legen3ha, legen4ha, nc)
           call legendre2phasefunction(legenha, nlegenha, 2, 200, p11, ang)
        else
           kextha(nz) = 0.0d0
           salbha = 0.0d0
           backha(nz) = 0.0d0
        endif
     else
        kextha(nz) = 0.0d0
        salbha = 0.0d0
        backha(nz) = 0.0d0
     endif

     nlegen(nz) = max(nlegen(nz),nlegencw,nlegenci,nlegenrr,nlegensn,nlegengr,nlegenha)

     if (verbose .gt. 1) print*, 'End of scattering calc for layer: ', nz

     !CCCCCCCCCCCCC   END OF SINGLE SCATTERING PROPERTY COMPUTATIONS  CCCCCCC

     !                                                                       
     !           Summing up the scattering parameters and writing the
     !           input file of the scattering properties of each layer
     !                                                                       


     kexttot(nz) = kextsn(nz) + kextcw(nz) + kextrr(nz) + kextgr(nz) + kextci(nz) + kextha(nz)
     back(nz) = backcw(nz) + backrr(nz) + backci(nz) + backsn(nz) + backgr(nz) + backha(nz)

     if (kexttot(nz) .lt. 0.) write(*,*) 'something wrong'
     if (kexttot(nz) .le. 0.) then 
        salbtot(nz) = 0.d0
     else 
        salbtot(nz) = (salbcw * kextcw(nz) + salbrr *       &
             kextrr(nz) + salbci * kextci(nz) + salbsn * kextsn(nz) + salbgr *    &
             kextgr(nz) + salbha * kextha(nz)) / kexttot(nz)
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
           legen (nz,jj) = (legencw (jj) * salbcw * kextcw(nz) + legenrr ( &
                jj) * salbrr * kextrr(nz) + legenci (jj) * salbci * kextci(nz) + &
                legensn (jj) * salbsn * kextsn(nz) + legengr (jj) * salbgr * &
                kextgr(nz) + legenha (jj) * salbha * kextha(nz)) / (salbtot (nz) &
                * kexttot (nz) )

           legen2 (nz,jj) = (legen2cw (jj) * salbcw * kextcw(nz) +         &
                legen2rr (jj) * salbrr * kextrr(nz) + legen2ci (jj) * salbci &
                * kextci(nz) + legen2sn (jj) * salbsn * kextsn(nz) + legen2gr (  &
                jj) * salbgr * kextgr(nz) + legen2ha (jj) * salbha * kextha(nz) ) &
                / (salbtot (nz) * kexttot (nz) )

           legen3 (nz,jj) = (legen3cw (jj) * salbcw * kextcw(nz) +         &
                legen3rr (jj) * salbrr * kextrr(nz) + legen3ci (jj) * salbci &
                * kextci(nz) + legen3sn (jj) * salbsn * kextsn(nz) + legen3gr (  &
                jj) * salbgr * kextgr(nz) + legen3ha (jj) * salbha * kextha(nz)) &
                / (salbtot (nz) * kexttot (nz) )

           legen4 (nz,jj) = (legen4cw(jj) * salbcw * kextcw(nz) +         &
                legen4rr (jj) * salbrr * kextrr(nz) + legen4ci (jj) * salbci &
                * kextci(nz) + legen4sn (jj) * salbsn * kextsn(nz) + legen4gr (  &
                jj) * salbgr * kextgr(nz) + legen4ha (jj) * salbha * kextha(nz)) &
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

  if (verbose .gt. 1) print*, 'Exiting hydrometeor_extinction_rt3'

1005 format  (i3,6(1x,f10.7))

  return

end subroutine hydrometeor_extinction_rt3
