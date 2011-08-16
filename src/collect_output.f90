subroutine collect_output(NSTOKES, NUMMU, AZIORDER,   &
     WAVELENGTH, UNITS, OUTPOL,NOUTLEVELS,        &
     OUTLEVELS, NUMAZIMUTHS,UP_RAD,    &
     DOWN_RAD,a,b,fi)


  use kinds
  use constants, only: pi
  use vars_output
  !   use vars_atmosphere, only: ngridx, ngridy

  integer, intent(in) :: a,b,fi

  INTEGER NSTOKES, NUMMU, AZIORDER  !,NUMAZI
  INTEGER NOUTLEVELS, OUTLEVELS ( * ), NUMAZIMUTHS
  REAL(kind=dbl) WAVELENGTH 
  REAL(kind=dbl) UP_RAD (NSTOKES, NUMMU, AZIORDER + 1, NOUTLEVELS) 
  REAL(kind=dbl) DOWN_RAD (NSTOKES, NUMMU, AZIORDER + 1, NOUTLEVELS)

  INTEGER I, J, K, L, LI, M, N 
  REAL OUT (4), PHI, PHID 


  !                                                             
  !       integer :: model_i, model_j
  !       real lon,lat,lfrac,wind10u,wind10v,iwv,cwp,iwp,rwp,swp,gwp

  N = NUMMU * (AZIORDER + 1) * NOUTLEVELS 
  CALL CONVERT_OUTPUT(UNITS, OUTPOL, NSTOKES, N, WAVELENGTH, 0,UP_RAD)
  CALL CONVERT_OUTPUT(UNITS, OUTPOL, NSTOKES, N, WAVELENGTH, 0,DOWN_RAD)

  NUMAZI = 2 * AZIORDER + 1 

  do l = 1, noutlevels 
     li = outlevels (l) 
     !    For each azimuth and zenith at this level sum the Fourier
     !    azimuth series appropriate for the particular Stokes parameter
     !    and output the radiance.                                
     do k = 1, numazimuths 
        if (numazimuths.eq.1) then 
           phid = 0.0 
        else 
           phid = 180.0 * float (k - 1) / (numazimuths - 1) 
        end if
        phi = pi * phid / 180.0 
        !               output upwelling radiance: -1 < mu < 0                  
        do j = nummu, 1, - 1 
           do i = 1, nstokes 
	      out (i) = 0.0 
	      do m = 0, aziorder 
                 if (i.le.2) then 
                    out(i) = out(i) + cos(m * phi) * up_rad(i, j, m + 1, l) 
                 else 
                    out(i) = out(i) + sin(m * phi) * up_rad(i, j, m + 1, l) 
                 end if
	      end do
           end do
           do i = 1, nstokes
              !	      tb_up(a,b,l,j,i) = out(i)
	      tb(i,fi,(nummu+1)-j,l,b,a) = out(i)
           end do
           !	    write (3, form1) height (li), phid, - mu_values (j), (out (i),i = 1, nstokes)                                                   
        end do
        !               output downwelling radiance: 0 < mu < 1                 
        do j = 1, nummu 
           do i = 1, nstokes 
	      out(i) = 0.0 
	      do m = 0, aziorder 
                 if (i.le.2) then 
                    out(i) = out(i) + cos(m * phi) * down_rad(i, j, m + 1, l) 
                 else 
                    out(i) = out(i) + sin(m * phi) * down_rad(i, j, m + 1, l) 
                 end if
	      end do
           end do
           do i = 1, nstokes
              !	      tb_down(a,b,l,j,i) = out(i)
	      tb(i,fi,nummu+j,l,b,a) = out(i)
           end do
           !	    write (3, form1) height (li), phid, mu_values(j), (out(i),i = 1, nstokes)                                                   
        end do
     end do
  end do
  return

end subroutine collect_output
