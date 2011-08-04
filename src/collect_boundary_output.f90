subroutine collect_boundary_output(GROUND_TEMP,lon,lat,lfrac,wind10u,wind10v,iwv, &
	cwp,iwp,rwp,swp,gwp,model_i,model_j,a,b)

  use kinds
  use vars_output


  integer, intent(in) :: a,b
                                                                        
      integer :: model_i, model_j
      real lon,lat,lfrac,wind10u,wind10v,iwv,cwp,iwp,rwp,swp,gwp
                                                                        


!      Output integrated quantities

      is(b,a) = model_i
      js(b,a) = model_j
      lons(b,a) = lon
      lats(b,a) = lat
      lfracs(b,a) = lfrac
      t_g(b,a) = GROUND_TEMP
      w10u(b,a) = wind10u
      w10v(b,a) = wind10v
      iwvs(b,a) = iwv
      cwps(b,a) = cwp
      iwps(b,a) = iwp
      rwps(b,a) = rwp
      swps(b,a) = swp
      gwps(b,a) = gwp

  return

end subroutine collect_boundary_output
