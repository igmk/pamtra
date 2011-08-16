! Get and interpolate land surface emissivities
!
subroutine land_emis (&
     ise, & ! in
     lon, & ! in
     lat, & ! in
     freqdbl,& ! in
     emissivity)  ! out
  !
  ! Description:
  ! To read and interpolate surface emissivities provided
  ! by Catherine Prigent and received from the SSMI satellite
  ! for all frequencies
  ! 
  ! Method:
  ! Catherine Prigent AMSU-A land surface emissivity estimation for numerical
  ! weather prediction assimilation schemes 
  !
  ! History:
  ! Version  Date        Comment
  ! -------  ----        -------
  !  0.1     14/01/2005  Adopted to mesonh2mwmod (M Mech)
  ! 
  ! Modules used:

  use equcom
  use kinds

  implicit none
  !subroutine arguments:
  integer, intent(in) :: ise

  real(kind=sgl), intent(in) :: lon,lat

  !local constants:
  real(kind=sgl), parameter :: dlat=0.25

  !local variables:
  integer :: ilon,ilat,cellnum,i,ii,ilon3,ilat3
  integer :: ilat4,ilon4,nb

  integer, dimension(4) :: cell

  real(kind=sgl) :: alat,alon,alat3,alon3,alat4,alon4
  real(kind=sgl) :: emiv,emih,theta,ff
  real(kind=dbl), intent(in) :: freqdbl ! frequency in [GHz]
  real(kind=sgl) :: freq ! frequency in [GHz]

  real(kind=sgl), dimension(4) :: a
  real(kind=sgl), dimension(7) :: emis
  real(kind=sgl), dimension(4,7) :: emis_p

  real(kind=dbl), intent(out) :: emissivity

  ! Equal area computations

freq = real(freqdbl)

call equare


  !-------------------------------------------
  ! 1. Get and interpolate emissivities
  !-------------------------------------------

  !-------------------------------------------
  ! 1.1 Calculate number of dataset to request
  !-------------------------------------------
  theta = 50

  if (lon .lt. 0.) then
     alon=lon+360.
  else
     alon=lon
  end if
  
  alat = lat

  cellnum=0
  ilat=int((alat+90.)/dlat)+1
  ilon=int(alon/(360./ncells(ilat)))+1
  do i=1,ilat-1
     cellnum=ncells(i)+cellnum
  end do

  ! dataset number
  cellnum=cellnum+ilon

  ! get emissivities

  read(ise,rec=cellnum) (emis(i), i= 1,7)
  !----------------------------------
  ! 1.2 Check if we have emissivities
  !----------------------------------
  if (emis(1) .lt. 0.1) then

     do ii=1,7
        emis(ii)=0.
     end do

     ! * bouchage de trous sur les 4 les plus proches attention passage par 0
     ! * cell(1)=cellnul+1, cell(2)=cellnum-1, et les 2 autres celll(3) et 
     ! * cell(4)

     cell(1)=cellnum+1
     cell(2)=cellnum-1
     alat3=alat+dlat
     alon3=alon
     if (alon3.le.0) alon3=alon3+360.
     ilat3=int((alat3+90.)/dlat)+1
     ilon3=int(alon3/(360./ncells(ilat3)))+1
     cell(3)=0
     do ii=1,ilat3-1
        cell(3)=ncells(ii)+cell(3)
     end do
     cell(3)=cell(3)+ilon3 

     alat4=alat-dlat
     alon4=alon
     if (alon4.le.0) alon4=alon4+360.
     ilat4=int((alat4+90.)/dlat)+1
     ilon4=int(alon4/(360./ncells(ilat4)))+1
     cell(4)=0
     do ii=1,ilat4-1
        cell(4)=ncells(ii)+cell(4)
     end do
     cell(4)=cell(4)+ilon4 

     nb=0
     do i=1,4
        read(ise,rec=cell(i)) (emis_p(i,ii),ii=1,7)
        if (emis_p(i,1).gt.0.1) then
           nb=nb+1
           do ii=1,7
              emis(ii)=emis_p(i,ii)+emis(ii)
           end do
        end if
     end do
     ! If all 4 gridpoints of the interpolation don't have data
     if (nb == 0) then
        cell(1)=cellnum+2
        cell(2)=cellnum-2
        alat3=alat+2*dlat
        alon3=alon
        if (alon3.le.0) alon3=alon3+360.
        ilat3=int((alat3+90.)/dlat)+1
        ilon3=int(alon3/(360./ncells(ilat3)))+1
        cell(3)=0
        do ii=1,ilat3-1
           cell(3)=ncells(ii)+cell(3)
        end do
        cell(3)=cell(3)+ilon3 
        
        alat4=alat-2*dlat
        alon4=alon
        if (alon4.le.0) alon4=alon4+360.
        ilat4=int((alat4+90.)/dlat)+1
        ilon4=int(alon4/(360./ncells(ilat4)))+1
        cell(4)=0
        do ii=1,ilat4-1
           cell(4)=ncells(ii)+cell(4)
        end do
        cell(4)=cell(4)+ilon4 
        
        nb=0
        do i=1,4
           read(ise,rec=cell(i)) (emis_p(i,ii),ii=1,7)
           if (emis_p(i,1).gt.0.1) then
              nb=nb+1
              do ii=1,7
                 emis(ii)=emis_p(i,ii)+emis(ii)
              end do
           end if
        end do
     end if
     if (nb.ge.2) then
        do ii=1,7
           emis(ii)=emis(ii)/(nb*1.)
        end do
     end if
  end if
  
  !------------------------------------
  ! 1.3 Interpolate to used frequencies
  !------------------------------------
     if (freq .lt. 19.35) then
        emiv=emis(1)
        emih=emis(2)
     end if
     if (freq .ge. 19.35 .and. freq .lt. 37.) then
        emiv=emis(1)+(freq-19.35)*(emis(1)-emis(4))/(19.35-37.)
        emih=emis(2)+(freq-19.35)*(emis(2)-emis(5))/(19.35-37.)
     end if
     if (freq .ge. 37. .and. freq .lt. 85.5) then
        emiv=emis(4)+(freq-37.)*(emis(4)-emis(6))/(37.-85.5)
        emih=emis(5)+(freq-37.)*(emis(5)-emis(7))/(37.-85.5)
     end if
     if (freq .ge. 85.5) then 
        emiv=emis(6)
        emih=emis(7)
     end if

     if (freq .lt. 85.5) then
        ff=freq
     else
        ff=85.5
     end if

! We are not dealing with that, do we?
!
!     if (ipol(freq).eq.2) then
!        ! * case H at nadir
        a(1)=.00327*ff+0.08
        a(2)=-3.90e-5*ff-0.00421
        a(3)=3.02e-6*ff-0.000046
        a(4)=-6.6e-8*ff+1.8e-6
!     else
        ! * case V at nadir
!        a(1)=.00327*ff+0.08
!        a(2)=-4.74e-5*ff-0.00529
!        a(3)=3.26e-6*ff+0.000475
!        a(4)=-6.6e-8*ff-7.7e-6
!     end if

     emiv=(emiv+emih)/2.+(emiv-emih)*(a(1)+a(2)*theta+&
         a(3)*theta*theta+a(4)*theta*theta*theta)
     emih=emiv

     emissivity = dble((emiv+emih)/2.d0)

  return

end subroutine land_emis
