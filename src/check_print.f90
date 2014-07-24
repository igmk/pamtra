subroutine check_print()

! Modules used:

  use kinds, only: dbl, &      ! integer parameter specifying double precision
                   long        ! integer parameter specifying long integer

  use drop_size_dist, only: dist_name, p_1, p_2, p_3, p_4, moment_in, n_0, lambda, gam, mu,  &
                            n_t, sig, d_ln, d_mono, nbin, m_0, m_32, am_b, q_h, n_tot, r_eff,&
                            d_ds, d_bound_ds, n_ds, f_ds, hydro_name
  use vars_index, only: i_x,i_y, i_z, i_f, i_h

  implicit none

  open(112,file='check_moments',position='append')
  

     write(112,'(a9,i5,i5,i5,i5,i5)')"x y z f h",i_x,i_y, i_z, i_f, i_h
     write(112,'(a10,a10,5x,i5,4(3x,f20.10))')hydro_name,dist_name,moment_in,p_1,p_2,p_3,p_4

! total number concentration
  if (moment_in == 1 .or. moment_in == 12 .or. moment_in == 13) then
    if (ABS(m_0-n_tot)/n_tot*100. >= 3) &
     write(112,'(2(a10,5x,f20.10,5x),a12,5x,f20.10,a2)')'m_0',m_0,'n_tot [m-3]',n_tot,'diff[%]',(m_0-n_tot)/n_tot*100.,' %'
  endif
! effective radius
  if (moment_in == 2 .or. moment_in == 12 .or. moment_in == 23) then
     if (ABS(m_32-r_eff)/r_eff*100. >= 3) &
      write(112,'(2(a10,5x,f20.10,5x),a12,5x,f20.10,a2)')'m_32',m_32,'r_eff [m]',r_eff,'diff[%]',(m_32-r_eff)/r_eff*100.,' %'
  endif
! total mass concentration
  if (moment_in == 3 .or. moment_in == 13 .or. moment_in == 23) then
     if (ABS(am_b-q_h)/q_h*100. .gt. 3) &
      write(112,'(2(a10,5x,f20.10,5x),a12,5x,f20.10,a2)')'m_3',am_b,'mass[kg/m3]',q_h,'diff[%]',(am_b-q_h)/q_h*100.,' %'
  endif

  close(112)
  return
end