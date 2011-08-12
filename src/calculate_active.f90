subroutine calculate_active(OUT_FILE_ACT,freq,hgt,Ze, Attenuation_atmo, Attenuation_hydro)
! This function computes and writes Ze and PIA
! PIA is calculated bottom (BU) up AND top down (TD) for hydrometeors and gaseous extinction seperately.


use kinds 
use vars_atmosphere !gives kextatmo, kexttot, temp, nlyr
use constants
use nml_params

implicit none

integer :: nz
real(kind=dbl) :: K2, dielec_water, tau_hydro, tau_atmo, d_hgt,wavelength
character(300), intent(in) ::OUT_FILE_ACT 
real(kind=dbl), intent(in) :: freq
real(kind=dbl), dimension(nlyr), intent(out) ::      hgt, Ze, Attenuation_atmo, Attenuation_hydro

		wavelength = c / (freq*1.d9)   ! m
		tau_hydro = 0.d0
		tau_atmo = 0.d0
		do nz = 1, nlyr
			hgt(nz) = (hgt_lev(nz-1)+hgt_lev(nz))*0.5d0
			d_hgt = hgt_lev(nz) - hgt_lev(nz-1)
			K2 = dielec_water(0.D0,temp(nz)-t_abs,freq)
			Ze(nz) = 10*log10(1d18* (1d0/ (K2*pi**5) ) * back(nz) * (wavelength)**4)
			if (abs(Ze(nz)) .gt. huge(Ze(nz))) Ze(nz) = -9999.d0
            Attenuation_atmo(nz) = 10*log10(exp(kextatmo(nz)*d_hgt))
            Attenuation_hydro(nz) = 10*log10(exp(kexttot(nz)*d_hgt))

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

		if (write_nc .eqv. .false.) then
			open (unit=22, file=OUT_FILE_ACT, status='unknown')
			write (22,*) "C           z[m]           Ze[dBz] Attenuation_hydro[dB] Attenuation_atmo[dB]"
			do nz = 1, nlyr
				write (22,2222) hgt(nz), Ze(nz), Attenuation_hydro(nz), Attenuation_atmo(nz)
				2222 format(1x, f16.4,1x, f16.4,1x, f16.4,1x, f16.4)
			end do
			close(22)
		end if



end subroutine calculate_active
