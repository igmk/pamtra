subroutine calculate_active(OUT_FILE_ACT,freq,hgt,Ze, PIA_atmo_BU, PIA_hydro_BU, PIA_atmo_TD, PIA_hydro_TD)
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
real(kind=dbl), dimension(nlyr), intent(out) ::      hgt, Ze, PIA_atmo_BU, PIA_hydro_BU, PIA_atmo_TD, PIA_hydro_TD



		wavelength = c / (freq*1.d3)   ! microns
		tau_hydro = 0.d0
		tau_atmo = 0.d0
		do nz = 1, nlyr
			hgt(nz) = (hgt_lev(nz-1)+hgt_lev(nz))*0.5d0
			d_hgt = hgt_lev(nz) - hgt_lev(nz-1)
			K2 = dielec_water(0.D0,temp(nz),freq)
			Ze(nz) = 10*log10(1d18* (1/ (pi*K2) ) * back(nz) * (wavelength*1d-6)**4)

			if (abs(Ze(nz)) .gt. huge(Ze(nz))) Ze(nz) = -9999.d0
			tau_hydro = tau_hydro + (kexttot(nz)*d_hgt)
			tau_atmo = tau_atmo + (kextatmo(nz)*d_hgt)
			PIA_atmo_BU(nz) = 10*log10(exp(2 * tau_atmo))
			PIA_hydro_BU(nz) =  10*log10(exp(2 * tau_hydro))

		end do

 		tau_hydro = 0.d0
		tau_atmo = 0.d0
		do nz = nlyr,1,-1
			d_hgt = hgt_lev(nz) - hgt_lev(nz-1)
			tau_hydro = tau_hydro + (kexttot(nz)*d_hgt)
			tau_atmo = tau_atmo + (kextatmo(nz)*d_hgt)
			PIA_atmo_TD(nz) = 10*log10(exp(2 * tau_atmo))
			PIA_hydro_TD(nz) =  10*log10(exp(2 * tau_hydro))

		end do

		if (write_nc .eqv. .false.) then
			open (unit=22, file=OUT_FILE_ACT, status='unknown')
! 			write (22,*) "C SD snow   : ", SD_snow
! 			write (22,*) "C N0 snow   : ", N0snowstr
! 			write (22,*) "C EM snow   : ", EM_snow
! 			write (22,*) "C SP        : ", SP_str
! 			write (22,*) "C SD graupel: ", SD_grau
! 			write (22,*) "C N0 graupel: ", N0graustr
! 			write (22,*) "C EM graupel: ", EM_grau
! 			write (22,*) "C SD rain   : ", SD_rain
! 			write (22,*) "C N0 rain   : ", N0rainstr
! 			write (22,*) "C"
			write (22,*) "C           z[m]           Ze[dBz] PIA_hydro_BU[dB]  PIA_atmo_BU[dB] PIA_hydro_TD[dB]  PIA_atmo_TD[dB]"
			do nz = 1, nlyr
				write (22,2222) hgt(nz), Ze(nz), PIA_hydro_BU(nz), PIA_atmo_BU(nz), PIA_hydro_TD(nz), PIA_atmo_TD(nz)
				2222 format(1x, f16.4,1x, f16.4,1x, f16.4,1x, f16.4,1x, f16.4,1x, f16.4)
			end do
			close(22)
		end if



end subroutine calculate_active