subroutine save_active(OUT_FILE_ACT,nx,ny,fi)
  ! This function writes Ze and PIA
 
  use kinds 
  use vars_atmosphere, only: atmo_nlyrs
  use vars_output
!   use constants
  use settings
        use report_module
!   use vars_index, only: i_x, i_y
  implicit none

  integer :: nz
  character(300), intent(in) ::OUT_FILE_ACT 
  integer, intent(in) :: nx,ny,fi

  if (verbose .gt. 1) print*, 'Saving active results to ',OUT_FILE_ACT

  if (radar_mode == "simple") then
     open (unit=22, file=OUT_FILE_ACT, status='unknown')
     write (22,*) "C           z[m]           Ze[dBz] Attenuation_hydro[dB] Attenuation_atmo[dB]"
     do nz = 1, atmo_nlyrs(nx,ny)
        write (22,2222) out_radar_hgt(nx,ny,nz), out_Ze(nx,ny,nz,fi), &
              out_att_hydro(nx,ny,nz,fi), out_att_atmo(nx,ny,nz,fi)
2222    format(1x, f16.4,1x, f16.4,1x, f16.4,1x, f16.4)
     end do
     close(22)

  else if ((radar_mode == "moments") .or. (radar_mode == "spectrum")) then
     open (unit=22, file=OUT_FILE_ACT, status='unknown')
     write (22,*) "C           z[m]          Ze[dBz] DopplerVel [m/s]",&
      " Spec Width [m/s]    Skewness [-]     Kurtosis [-] LeftSlope[dBs/m]",&
      " RightSlop[dBs/m] Atten_hydro[dB]    Atten_atmo[dB]"
     do nz = 1, atmo_nlyrs(nx,ny)
        write (22,3333) out_radar_hgt(nx,ny,nz), out_Ze(nx,ny,nz,fi), &
              out_radar_moments(nx,ny,nz,fi,1), out_radar_moments(nx,ny,nz,fi,2), &
              out_radar_moments(nx,ny,nz,fi,3), out_radar_moments(nx,ny,nz,fi,4), &
              out_radar_slopes(nx,ny,nz,fi,1), out_radar_slopes(nx,ny,nz,fi,2), &
              out_att_hydro(nx,ny,nz,fi), out_att_atmo(nx,ny,nz,fi)
3333    format(1x, f16.4,1x, f16.4,1x, f16.4,1x, f16.4,&
              f16.4,1x, f16.4,1x, f16.4,1x, f16.4,&
              f16.4,1x, f16.4,1x)
     end do
     close(22)

  else
    print*,"did not understand radar_mode", radar_mode
    stop
  end if

  if (verbose .gt. 1) print*, 'Active Results saved to ',OUT_FILE_ACT

end subroutine save_active
