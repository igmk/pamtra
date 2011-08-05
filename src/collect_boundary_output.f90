subroutine collect_boundary_output(lon,lat,lfrac,iwv, &
	cwp,iwp,rwp,swp,gwp,hwp,model_i,model_j,a,b)

  use kinds
  use vars_output

implicit none

  integer, intent(in) :: a,b
                                                                        
      integer, intent(in) :: model_i, model_j
      real(kind=sgl), intent(in) :: lon,lat,lfrac,iwv,cwp,iwp,rwp,swp,gwp,hwp
!      Output integrated quantities
      is(b,a) = model_i
      js(b,a) = model_j
      lons(b,a) = lon
      lats(b,a) = lat
      lfracs(b,a) = lfrac
      iwvs(b,a) = iwv
      cwps(b,a) = cwp
      iwps(b,a) = iwp
      rwps(b,a) = rwp
      swps(b,a) = swp
      gwps(b,a) = gwp
      hwps(b,a) = hwp

  return

end subroutine collect_boundary_output
