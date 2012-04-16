subroutine deallocate_profile_vars

  use vars_atmosphere
  use vars_output
  use nml_params
  use mod_io_strings

  implicit none


  deallocate(hgt_lev,press_lev,press,temp_lev,temp,relhum_lev,relhum,&
    vapor_pressure,rho_vap,q_hum,cwc_q,iwc_q,rwc_q,swc_q,gwc_q)

  if (n_moments .eq. 2) then
     deallocate(hwc_q,cwc_n,iwc_n,rwc_n,swc_n,gwc_n,hwc_n)
  end if

  deallocate(nlegen,rt3nlegen,kextatmo,kexttot,salbtot,rt3kexttot,rt3salbtot,g_coeff,&
    back,legen,legen2,legen3,legen4,rt3legen,rt3legen2,rt3legen3,rt3legen4,&
    rt4hydros_present,hydros_present)

  if (dump_to_file) deallocate(file_ph)



!   deallocate(angles_deg)
!   deallocate(profiles)


end subroutine deallocate_profile_vars
