subroutine hongdb(errorstatus,freq,temp,hong_type,dia,abs_hong,ext_hong,bsc_hong,g,p_hong)
!
! This routine gives access to the Hong database described in Hong et al  2009 JGR. The database contains 
! single scattering properties (absorption and scattering cross section and phase function elements) for six
! particles calculated with the DDA method.
!
! The database is limited to:
! - frequency range 90 to 874 GHz
! - diameters 2 to 2000 microns
! - temperature 223.15 to 273.15 K
!
! The original database used an old refractive index model after Warren and included only a single temperature
! at 243.15 K. Therefore a correction with the Maetzler refractive index to the absorption coefficient and an extrapolation in temperature of the
! this coefficient after Eriksson et al 2015 AMT has been performed.
!
! Due to the limited range in diameters, a loss of mass occurs.
! 
! Up to now, only the first element of the phase function is used.!!!

  use kinds
  use report_module
  use settings, only: data_path
  
  implicit none
  
  integer :: i, sf,st,sd,floc, tloc, dloc, pos1,pos2,pos3,pos4,pos5,pos6,pos7,pos8
  
  integer, intent(in) :: hong_type
  real(kind=dbl), intent(in) :: freq,temp,dia
  real(kind=dbl), intent(out) :: abs_hong,ext_hong,bsc_hong,g
  real(kind=dbl), dimension(181), intent(out) :: p_hong
  real(kind=dbl) :: area, abs_eff, ext_eff
  
  ! the intermediate datasets
  real(kind=dbl) :: a111,e111
  real(kind=dbl) :: a112,e112
  real(kind=dbl) :: a121,e121
  real(kind=dbl) :: a122,e122
  real(kind=dbl) :: a211,e211
  real(kind=dbl) :: a212,e212
  real(kind=dbl) :: a221,e221
  real(kind=dbl) :: a222,e222
  real(kind=dbl), dimension(181) :: p111,p112,p121,p122,p211,p212,p221,p222,dum_ar
  
  real(kind=dbl) :: freq_l, freq1, freq2, temp_l, temp1, temp2, dia_l, dia1, dia2
  real(kind=dbl), dimension(38) :: dias
  real(kind=dbl), dimension(21) :: freqs
  real(kind=dbl), dimension(5) :: temps
  
  real(kind=dbl) :: dum
!  real(kind=dbl), dimension(11) :: dummy
  
  integer(kind=long) :: errorstatus
  integer(kind=long) :: err
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'hongdb'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  freq_l = freq
  temp_l = temp
  dia_l = dia*1.e6 ! diameter in database is in micron
  err = 0
  abs_eff = 0._dbl
  abs_hong = 0._dbl
  ext_eff = 0._dbl
  ext_hong = 0._dbl
  bsc_hong = 0._dbl
  g = 0._dbl
  p_hong(:) = 0._dbl
  
  freqs = (/90.0_dbl, 118.0_dbl, 157.0_dbl, 166.0_dbl, 183.3_dbl, 190.0_dbl, 203.0_dbl, 220.0_dbl, 243.0_dbl, 325.0_dbl, 340.0_dbl,&
    380.0_dbl, 425.0_dbl, 448.0_dbl, 463.0_dbl, 487.0_dbl, 500.0_dbl, 640.0_dbl, 664.0_dbl, 683.0_dbl, 874.0_dbl/)
  temps =(/223.15_dbl,238.15_dbl,243.15_dbl,258.15_dbl,273.15_dbl/)
  dias = (/2.0_dbl, 4.0_dbl, 6.0_dbl, 8.0_dbl, 10.0_dbl, 12.5_dbl, 15.0_dbl, 20.0_dbl, 25.0_dbl, 30.0_dbl, &
    40.0_dbl, 50.0_dbl, 60.0_dbl, 70.0_dbl, 80.0_dbl, 90.0_dbl, 100.0_dbl, 125.0_dbl, 150.0_dbl, 175.0_dbl, &
    200.0_dbl, 250.0_dbl, 300.0_dbl, 350.0_dbl, 400.0_dbl, 500.0_dbl, 600.0_dbl, 700.0_dbl, 800.0_dbl, 900.0_dbl,&
    1000.0_dbl, 1100.0_dbl, 1200.0_dbl, 1300.0_dbl, 1400.0_dbl, 1600.0_dbl, 1800.0_dbl, 2000.0_dbl/)

  ! we do not do any extrapolation in frequency, temperature and diameter
  
  if ((freq_l < freqs(1)) .or. (freq_l > freqs(21))) then
     errorstatus = fatal
     msg =  'No extraplation for frequency' 
     call report(errorstatus, msg, nameOfRoutine)
     return    
  end if

!   if ((dia_l < dias(1)) .or. (dia_l > dias(38))) then
!      errorstatus = fatal
!      msg =  'No extraplation for diameter' 
!      call report(errorstatus, msg, nameOfRoutine)
!      return    
!     
!   end if
! 
!   if ((temp_l < temps(1)) .or. (temp_l > temps(5))) then
!      errorstatus = fatal
!      msg =  'No extraplation for temperature!' 
!      call report(errorstatus, msg, nameOfRoutine)
!      return    
!   end if

  if (abs(freq_l - freqs(1)) < 0.00001) freq_l = freq_l + 0.00001
  if (abs(freq_l - freqs(21)) < 0.00001) freq_l = freq_l - 0.00001
  if (abs(temp_l - temps(1)) < 0.00001) temp_l = temp_l + 0.00001
  if (abs(temp_l - temps(5)) < 0.00001) temp_l = temp_l - 0.00001
  if (abs(dia_l - dias(1)) < 0.00001) dia_l = dia_l + 0.00001
  if (abs(dia_l - dias(38)) < 0.00001) dia_l = dia_l - 0.00001

  ! find the bounding frequencies

  floc = minloc(abs(freqs - freq_l), 1) ! index of nearest frequency
  if (freqs(floc) > freq_l) then
    floc = floc - 1
    freq2 = freqs(floc + 1)
    freq1 = freqs(floc)
  elseif (freqs(floc) < freq_l) then
    freq1 = freqs(floc)
    freq2 = freqs(floc+1)
  else
    freq1 = freq_l
    freq2 = freq1
  end if
  
  ! find the bounding termperatures

  tloc = minloc(abs(temps - temp_l), 1) ! index of nearest temperature
  if (temps(tloc) > temp_l) then
    tloc = tloc - 1
    temp2 = temps(tloc + 1)
    temp1 = temps(tloc)
  elseif (temps(tloc) < temp_l) then
    temp1 = temps(tloc)
    temp2 = temps(tloc+1)
  else
    temp1 = temp_l
    temp2 = temp1
  end if

  ! find the bounding diameters

  dloc = minloc(abs(dias - dia_l), 1) ! index if nearest diameter
  if (dias(dloc) > dia_l) then
    dloc = dloc - 1
    dia2 = dias(dloc + 1)
    dia1 = dias(dloc)
  elseif (dias(dloc) < dia_l) then
    dia1 = dias(dloc)
    dia2 = dias(dloc+1)
  else
    dia1 = dia_l
    dia2 = dia1
  end if
  
  ! calculate positions for all the datasets we need to perform the interpolation
  ! floc, tloc, and dloc are the indices of the lower boundary
  sf = size(freqs) ! 21
  st = size(temps) ! 5
  sd = size(dias)  ! 38
  pos1 = 1 + hong_type * sf * st * sd + (floc - 1) * st * sd + (tloc - 1) * sd + (dloc - 1)
  pos2 = 1 + hong_type * sf * st * sd + (floc - 1) * st * sd + (tloc - 1) * sd + dloc
  pos3 = 1 + hong_type * sf * st * sd + (floc - 1) * st * sd + (tloc ) * sd + (dloc - 1)
  pos4 = 1 + hong_type * sf * st * sd + (floc - 1) * st * sd + (tloc) * sd + dloc
  pos5 = 1 + hong_type * sf * st * sd + (floc) * st * sd + (tloc - 1) * sd + (dloc - 1)
  pos6 = 1 + hong_type * sf * st * sd + (floc) * st * sd + (tloc - 1) * sd + dloc
  pos7 = 1 + hong_type * sf * st * sd + (floc) * st * sd + (tloc) * sd + (dloc - 1)
  pos8 = 1 + hong_type * sf * st * sd + (floc) * st * sd + (tloc) * sd + dloc
  if (verbose > 3) print*, pos1, pos2, pos3, pos4, pos5, pos6, pos7, pos8
  ! The record length is i4+10*f4+(8*181)*f4 = 5836
!  open(666,file='/home/mech/workspace/HongLiuDB/scripts/test.dda',form='unformatted',access='direct',recl=5836)
  ! The record length is i8+10*f8+(8*181)*f8 = 11672
  open(666,file=trim(data_path)//'/hongdb/hong.dda',form='unformatted',access='direct',recl=11672)
  ! we have to get datasets for all combination
  ! freq1,temp1,dia1
  read(666,rec=pos1) dum, dum, dum, dum, dum, dum, area, e111,a111,dum,dum,p111,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar
  if (verbose > 3) print*,a111,e111
  ! freq1,temp1,dia2
  read(666,rec=pos2) dum, dum, dum, dum, dum, dum, dum, e112,a112,dum,dum,p112,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar
  if (verbose > 3) print*,a112,e112
  ! freq1,temp2,dia1
  read(666,rec=pos3) dum, dum, dum, dum, dum, dum, dum, e121,a121,dum,dum,p121,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar
  if (verbose > 3) print*,a121,e121
  ! freq1,temp2,dia2
  read(666,rec=pos4) dum, dum, dum, dum, dum, dum, dum, e122,a122,dum,dum,p122,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar
  if (verbose > 3) print*,a122,e122
  ! freq2,temp1,dia1
  read(666,rec=pos5) dum, dum, dum, dum, dum, dum, dum, e211,a211,dum,dum,p211,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar
  if (verbose > 3) print*,a211,e211
  ! freq2,temp1,dia2
  read(666,rec=pos6) dum, dum, dum, dum, dum, dum, dum, e212,a212,dum,dum,p212,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar
  if (verbose > 3) print*,a212,e212
  ! freq2,temp2,dia1
  read(666,rec=pos7) dum, dum, dum, dum, dum, dum, dum, e221,a221,dum,dum,p221,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar
  if (verbose > 3) print*,a221,e221
  ! freq2,temp2,dia2
  read(666,rec=pos8) dum, dum, dum, dum, dum, dum, dum, e222,a222,dum,dum,p222,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar,dum_ar
  if (verbose > 3) print*,a222,e222
  close(666)

  ! interpolation
!  print*, dia_l,dia1,dia2,temp_l,temp1,temp2,freq_l,freq1,freq2,a111,a112,a121,a122,a211,a212,a221,a222,abs_hong
  ! absorption cross section
  call intpol(dia_l,dia1,dia2,temp_l,temp1,temp2,freq_l,freq1,freq2,a111,a112,a121,a122,a211,a212,a221,a222,abs_eff)
  abs_hong = abs_eff*area*1.e-12 ! rescale from 1/micron**2 to 1/m**2
  if (verbose > 3) print*, dia_l,dia1,dia2,temp_l,temp1,temp2,freq_l,freq1,freq2,a111,a112,a121,a122,a211,a212,a221,a222,abs_hong
  ! scattering cross section
  call intpol(dia_l,dia1,dia2,temp_l,temp1,temp2,freq_l,freq1,freq2,e111,e112,e121,e122,e211,e212,e221,e222,ext_eff)
  ext_hong = ext_eff * area*1.e-12 ! rescale from 1/micron**2 to 1/m**2
  if (verbose > 3) print*, dia_l,dia1,dia2,temp_l,temp1,temp2,freq_l,freq1,freq2,e111,e112,e121,e122,e211,e212,e221,e222,ext_hong
  ! back scattering cross section
!   call intpol(dia_l,dia1,dia2,temp_l,temp1,temp2,freq_l,freq1,freq2,b111,b112,b121,b122,b211,b212,b221,b222,bsc_hong)
!   print*, dia_l,dia1,dia2,temp_l,temp1,temp2,freq_l,freq1,freq2,b111,b112,b121,b122,b211,b212,b221,b222,bsc_hong
  ! phase function
  do i = 1,181
    call intpol(dia_l,dia1,dia2,temp_l,temp1,temp2,freq_l,freq1,freq2,&
      p111(i),p112(i),p121(i),p122(i),p211(i),p212(i),p221(i),p222(i),p_hong(i))
  end do
  
  bsc_hong = p_hong(181)*area*1.e-12 ! rescale from 1/micron**2 to 1/m**2

  errorstatus = err    
  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

  return
  
end subroutine hongdb

subroutine intpol(dia,dia1,dia2,temp,temp1,temp2,freq,freq1,freq2,t111,t112,t121,t122,t211,t212,t221,t222,res)

  use kinds
  
  implicit none
  
  real(kind=dbl) :: dia,dia1,dia2,temp,temp1,temp2,freq,freq1,freq2
  real(kind=dbl) :: t111,t112,t121,t122,t211,t212,t221,t222,res
  real(kind=dbl) :: tmp1, tmp2,tmp3,tmp4,tmp5,tmp6
  
  if (dia1 == dia2) then
    tmp1 = t111
    tmp2 = t121
    tmp3 = t211
    tmp4 = t221
  else
    call interpolation(2,1,(/dia1,dia2/),(/t111,t112/),dia,tmp1)
    call interpolation(2,1,(/dia1,dia2/),(/t121,t122/),dia,tmp2)
    call interpolation(2,1,(/dia1,dia2/),(/t211,t212/),dia,tmp3)
    call interpolation(2,1,(/dia1,dia2/),(/t221,t222/),dia,tmp4)
  end if
  if (temp1 == temp2) then
    tmp5 = tmp1
    tmp6 = tmp3
  else
    call interpolation(2,1,(/temp1,temp2/),(/tmp1,tmp2/),temp,tmp5)
    call interpolation(2,1,(/temp1,temp2/),(/tmp3,tmp4/),temp,tmp6)
  end if
!   call interpolation(1,2,(/temp1,temp2/),(//),temp,tmp)
!   call interpolation(1,2,(/temp1,temp2/),(//),temp,tmp)
  if (freq1 == freq2) then
    res = tmp5
  else
    call interpolation(2,1,(/freq1,freq2/),(/tmp5,tmp6/),freq,res)
  end if
!   call interpolation(1,2,(/freq1,freq2/),(//),freq,tmp)
!   call interpolation(1,2,(/freq1,freq2/),(//),freq,tmp)
!   call interpolation(1,2,(/freq1,freq2/),(//),freq,tmp)
!   call interpolation(1,2,(//),,,res)
  
  return
  
end subroutine intpol