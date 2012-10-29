subroutine calculate_active(OUT_FILE_ACT,freq,hgt,&
  Ze,Ze_cw,Ze_rr,Ze_ci,Ze_sn,Ze_gr,Ze_ha,&
  Att_atmo, Att_hydro, Att_cw,Att_rr,Att_ci,Att_sn,Att_gr,Att_ha)
  ! This function computes and writes Ze and PIA
 
  use kinds 
  use vars_atmosphere !gives kextatmo, kexttot, temp, nlyr
  use constants
  use nml_params

  implicit none

  integer :: nz
  real(kind=dbl) :: K2, dielec_water, tau_hydro, tau_atmo, d_hgt,wavelength
  character(300), intent(in) ::OUT_FILE_ACT 
  real(kind=dbl), intent(in) :: freq
  real(kind=dbl), dimension(nlyr), intent(out) :: hgt, &
          Ze,Ze_cw,Ze_rr,Ze_ci,Ze_sn,Ze_gr,Ze_ha, &
          Att_atmo, Att_hydro, Att_cw,Att_rr,Att_ci,Att_sn,Att_gr,Att_ha

  wavelength = c / (freq*1.d9)   ! m
  tau_hydro = 0.d0
  tau_atmo = 0.d0
  do nz = 1, nlyr
     hgt(nz) = (hgt_lev(nz-1)+hgt_lev(nz))*0.5d0
     d_hgt = hgt_lev(nz) - hgt_lev(nz-1)
     K2 = dielec_water(0.D0,temp(nz)-t_abs,freq)


     Ze(nz) = 1d18* (1d0/ (K2*pi**5) ) * back(nz) * (wavelength)**4
     Ze_cw(nz) = 1d18* (1d0/ (K2*pi**5) ) * backcw(nz) * (wavelength)**4
     Ze_rr(nz) = 1d18* (1d0/ (K2*pi**5) ) * backrr(nz) * (wavelength)**4
     Ze_ci(nz) = 1d18* (1d0/ (K2*pi**5) ) * backci(nz) * (wavelength)**4
     Ze_sn(nz) = 1d18* (1d0/ (K2*pi**5) ) * backsn(nz) * (wavelength)**4
     Ze_gr(nz) = 1d18* (1d0/ (K2*pi**5) ) * backgr(nz) * (wavelength)**4
     Ze_ha(nz) = 1d18* (1d0/ (K2*pi**5) ) * backha(nz) * (wavelength)**4

     Att_atmo(nz) = exp(kextatmo(nz)*d_hgt)
     Att_hydro(nz) = exp(kexttot(nz)*d_hgt)
     Att_cw(nz) = exp(kextcw(nz)*d_hgt)
     Att_rr(nz) = exp(kextrr(nz)*d_hgt)
     Att_ci(nz) = exp(kextci(nz)*d_hgt)
     Att_sn(nz) = exp(kextsn(nz)*d_hgt)
     Att_gr(nz) = exp(kextgr(nz)*d_hgt)
     Att_ha(nz) = exp(kextha(nz)*d_hgt)


    if (activeLogScale) then

      Ze(nz) = 10*log10(Ze(nz))
      if (abs(Ze(nz)) .ge. huge(Ze(nz))) Ze(nz) = -9999.d0
      Ze_cw(nz) = 10*log10(Ze_cw(nz))
      if (abs(Ze_cw(nz)) .ge. huge(Ze_cw(nz))) Ze_cw(nz) = -9999.d0
      Ze_rr(nz) = 10*log10(Ze_rr(nz))
      if (abs(Ze_rr(nz)) .ge. huge(Ze_rr(nz))) Ze_rr(nz) = -9999.d0
      Ze_ci(nz) = 10*log10(Ze_ci(nz))
      if (abs(Ze_ci(nz)) .ge. huge(Ze_ci(nz))) Ze_ci(nz) = -9999.d0
      Ze_sn(nz) = 10*log10(Ze_sn(nz))
      if (abs(Ze_sn(nz)) .ge. huge(Ze_sn(nz))) Ze_sn(nz) = -9999.d0
      Ze_gr(nz) = 10*log10(Ze_gr(nz))
      if (abs(Ze_gr(nz)) .ge. huge(Ze_gr(nz))) Ze_gr(nz) = -9999.d0
      Ze_ha(nz) = 10*log10(Ze_ha(nz))
      if (abs(Ze_ha(nz)) .ge. huge(Ze_ha(nz))) Ze_ha(nz) = -9999.d0

      Att_atmo(nz) = 10*log10(Att_atmo(nz))
      Att_hydro(nz) = 10*log10(Att_hydro(nz))
      Att_cw(nz) = 10*log10(Att_cw(nz))
      Att_rr(nz) = 10*log10(Att_rr(nz))
      Att_ci(nz) = 10*log10(Att_ci(nz))
      Att_sn(nz) = 10*log10(Att_sn(nz))
      Att_gr(nz) = 10*log10(Att_gr(nz))
      Att_ha(nz) = 10*log10(Att_ha(nz))
  end if


  end do


  ! 		!PIA is calculated like:
  !  		tau_hydro = 0.d0
  ! 		tau_atmo = 0.d0
  ! 		do nz = nlyr,1,-1
  ! 			d_hgt = hgt_lev(nz) - hgt_lev(nz-1)
  ! 			tau_hydro = tau_hydro + (kexttot(nz)*d_hgt)
  ! 			tau_atmo = tau_atmo + (kextatmo(nz)*d_hgt)
  ! 			PIA_atmo_TD(nz) = 10*log10(exp(2 * tau_atmo))
  ! 			PIA_hydro_TD(nz) =  10*log10(exp(2 * tau_hydro))
  ! 
  ! 		end do

  if ((write_nc .eqv. .false.) .and. (in_python .eqv. .false.)) then
     open (unit=22, file=OUT_FILE_ACT, status='unknown')
     write (22,*) "C           z[m]           Ze[dBz] Attenuation_hydro[dB] Attenuation_atmo[dB]"
     do nz = 1, nlyr
        write (22,2222) hgt(nz), Ze(nz), Att_hydro(nz), Att_atmo(nz)
2222    format(1x, f16.4,1x, f16.4,1x, f16.4,1x, f16.4)
     end do
     close(22)
  end if



end subroutine calculate_active
