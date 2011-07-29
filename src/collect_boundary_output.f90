subroutine collect_boundary_output(t_surf,lon,lat,lfrac,wind10u,wind10v,iwv, &
	cwp,iwp,rwp,swp,gwp,hwp,model_i,model_j,a,b)

  use kinds
  use vars_output

implicit none

  integer, intent(in) :: a,b
                                                                        
      integer, intent(in) :: model_i, model_j
      real(kind=dbl), intent(in) :: t_surf,lon,lat,lfrac,wind10u,wind10v,iwv,cwp,iwp,rwp,swp,gwp,hwp
                                                                        
!      Output integrated quantities
      is(b,a) = model_i
      js(b,a) = model_j
      lons(b,a) = lon
      lats(b,a) = lat
      lfracs(b,a) = lfrac
      t_g(b,a) = t_surf
      w10u(b,a) = wind10u
      w10v(b,a) = wind10v
      iwvs(b,a) = iwv
      cwps(b,a) = cwp
      iwps(b,a) = iwp
      rwps(b,a) = rwp
      swps(b,a) = swp
      gwps(b,a) = gwp
      hwps(b,a) = hwp

  return

end subroutine collect_boundary_output
