subroutine equare(ncells)

!   use equcom
  use kinds

  implicit none

  real(kind=sgl), dimension(720), intent(out) :: ncells
  real(kind=sgl), parameter :: dlat=0.25

  integer, parameter :: maxlat=180./dlat,maxlon=360./dlat

  integer :: totcel,maxlt2,lat,icellr,lat1,lat2

  real(kind=dbl) :: rcells,pi,aezon,hezon,aecell,xlatb,&
       rlatb,rlate,xlate,htb,hte,htzone,dlongr,rcelat,azone,&
       rcellr,acell,asq

  real(kind=sgl), parameter :: rearth = 6371.2

  pi = 4.0 * atan(1.0)
  rcelat = (dlat*pi)/180.
  totcel = 0                                                          
  hezon = rearth*sin(rcelat)                                          
  aezon = 2*pi*rearth*hezon                                          
  aecell = (aezon*dlat)/360.                                          
  maxlt2 = maxlat/2                                                   
  do lat = 1, maxlt2                                               
     xlatb = (lat-1)*dlat                                                
     xlate = xlatb+dlat                                                  
     rlatb = (2.d0*pi*xlatb)/360.d0
     rlate = (2.d0*pi*xlate)/360.d0
     htb = rearth*sin(rlatb)                                             
     hte = rearth*sin(rlate)                                             
     htzone = hte-htb                                                    
     azone = 2*pi*rearth*htzone                                         
     rcells = azone/aecell                                               
     icellr = (rcells+.5d0)                                                
     totcel = totcel+2*icellr                                            
     rcellr = icellr                                                     
     dlongr = 360.d0/rcellr                                                
     acell = azone/rcellr                                                
     asq = azone/maxlon                                                  
     lat1 = lat+maxlt2                                                   
     lat2 = maxlt2+1-lat                                                 
     ncells(lat1)=icellr                                               
     ncells(lat2)=icellr                                               
  end do
  
  return                                                            
  
end subroutine equare
