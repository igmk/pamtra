MODULE mod_mwatlas_nt_bin
!   This module has been written by
!   Catherine Prigent
!   Filipe Aires  
!  
!    November 2015
!
!   and modified for PAMTRA by
!
!    Mario Mech
!
!    November 2016

  use kinds
  use report_module

  implicit none

    ! number of lines in the atlas
    integer(kind=long) :: atlas_ndat
    ! number of channels in the atlas
    integer(kind=long) :: atlas_nchan=7
    ! name of the atlas (including version number)
    CHARACTER(len=22) :: atlas_name
    ! month of the atlas
    integer(kind=long) :: atlas_month                  
    ! resolution of the atlas (equal-area)
    real(kind=dbl) :: atlas_dlat=0.25
    ! number of cells per lat band
    integer(kind=long), dimension(720) :: atlas_ncells
    ! the first cellnumber of lat band
    integer(kind=long), dimension(720) :: atlas_firstcell
    ! Emissivities
    real(kind=dbl), allocatable,dimension(:,:) :: atlas_emis!emis(ndat,nchan)
    ! Surface classes
    integer(kind=long), allocatable,dimension(:) :: atlas_class1
    integer(kind=long), allocatable,dimension(:) :: atlas_class2
    ! cellnumber of each of the pixels in the atlas
    integer(kind=long), allocatable,dimension(:) :: atlas_cellnum
    ! "Correspondance" vector indicating that for the ith element, the j so that EMIS(j,...) is the emissivity of cellnumber i.
    integer(kind=long) :: atlas_correspondance(660066) 

!=======================================================
!ROUTINES ==============================================
!=======================================================
CONTAINS

  !------------------------------------------------------------------
  SUBROUTINE test_inputs(errorstatus, imonth,lat,lon,theta,freq)
    !======to test if inputs are correct
    integer(kind=long), INTENT (IN) :: imonth
    real(kind=dbl), INTENT (IN)    :: theta, freq
    real(kind=sgl), INTENT (IN)    :: lat
    real(kind=sgl) :: lon
    
    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg    
    character(22) :: nameOfRoutine = 'test_input'

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine) 
    
    call assert_false(err,(imonth<1 .OR. imonth>12),nameOfRoutine//': Month must be between 1 and 12.')
    call assert_false(err,(lat<-90 .OR. lat>90),nameOfRoutine//': Latitude must be between -90. and 90.')
    call assert_false(err,(lon<0 .OR. lon>360),nameOfRoutine//': Longitude must be between 0 and 360')
    call assert_false(err,(theta<0 .OR. theta>90),nameOfRoutine//': Zenith angle must be between 0 and 90')
    call assert_false(err,(freq<10 .OR. freq>700),nameOfRoutine//': Frequency must be between 10 and 700.')
    if (lat == 90.) then
        ! Otherwise the data reports an error
        lon = 0.
    end if
    
    if (err /= 0) then
      print*, imonth, lat, lon, theta, freq
      msg = 'error in value check in test_inputs for telsem2'
      call report(err, msg, nameOfRoutine)
      errorstatus = err
      return
    end if
    errorstatus = err
    
  if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
  
  END SUBROUTINE test_inputs


  !------------------------------------------------------------------
  SUBROUTINE rttov_readmw_atlas(errorstatus,dir,month)
    !======Read a monthly atlas
    integer(kind=long) :: imonth
    CHARACTER(len=*), INTENT (IN)       :: dir
    !----TRANSITORY VARIABLES
    integer(kind=long)   :: ipos
    integer(kind=long)   :: j,i
    integer(kind=long)   :: iiin=21    ! unit for input
    integer(kind=long)   :: irec
    CHARACTER(len=2), intent(in) :: month
    integer(kind=long)   :: cellnum
    real(kind=sgl)      :: ssmi(7)
    integer(kind=long)   :: cur_class1,cur_class2
    integer(kind=long)   :: take ! flag to take or not the pixel atlas for location constrains
    ! data records within each monthly file
    integer(kind=long), dimension(12), parameter :: ndata = (/224556,221099,221682,224855,227966,&
        229641,231213,231146,231298,232414,233123,230372/)
    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg    
    character(22) :: nameOfRoutine = 'rttov_readmw_atlas'
  
    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)    
    !initialisation
    err = 0
    
    inquire (iolength=irec) cellnum, ssmi(1:7),cur_class1,cur_class2

    !
    read(month,'(I2.2)') imonth
    atlas_month=imonth
    CALL equare_telsem2(atlas_dlat,atlas_ncells,atlas_firstcell)
    atlas_nchan=7
    atlas_name='ssmi_mean_emis_climato'
    !----INITIALISATIONS
    DO j=1,660066
       atlas_correspondance(j)=-777
    END DO
    open(unit=iiin, file=dir//'ssmi_mean_emis_climato_'//&
    month//'_cov_interpol_M2'//'.dat',access='direct',form='unformatted',recl=irec,iostat=err)
    IF( err /= 0 ) THEN
       msg = 'Error opening '//dir//'ssmi_mean_emis_climato_'//month//'_cov_interpol_M2'//'.dat'
       call report(err,msg, nameOfRoutine)
       errorstatus = err
     return
    END IF
    atlas_ndat=ndata(imonth)
    if (verbose >= 3 ) WRITE(0,*) 'Nb data=',atlas_ndat

    !----ALLOCATE VARIABLES
    ALLOCATE(atlas_emis(atlas_ndat,atlas_nchan)) 
    ALLOCATE(atlas_class1(atlas_ndat)) 
    ALLOCATE(atlas_class2(atlas_ndat)) 
    ALLOCATE(atlas_cellnum(atlas_ndat)) 

    ipos=0
    DO j=1,atlas_ndat
       READ(iiin,rec=j) cellnum, (ssmi(i),i=1,atlas_nchan),cur_class1,cur_class2

       take=1
       IF ((cur_class1 > 0).AND.(cur_class2 > 0).AND.(ipos<atlas_ndat).AND.(take==1)) THEN
          ipos=ipos+1
          DO i=1,atlas_nchan
             atlas_emis(ipos,i)=dble(ssmi(i))
          END DO
          atlas_cellnum(ipos)=cellnum 
          atlas_class1(ipos)=cur_class1
          atlas_class2(ipos)=cur_class2 
          atlas_correspondance(cellnum)=ipos 
       END IF          
    END DO
    atlas_ndat=ipos;
    CLOSE(iiin)

    errorstatus = err

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)
    
  END SUBROUTINE rttov_readmw_atlas

  !------------------------------------------------------------------
  SUBROUTINE rttov_closemw_atlas()
    !======free a monthly atlas
    DEALLOCATE(atlas_emis) 
    DEALLOCATE(atlas_class1) 
    DEALLOCATE(atlas_class2) 
    DEALLOCATE(atlas_cellnum) 
  END SUBROUTINE rttov_closemw_atlas
! 
!   !------------------------------------------------------------------
  SUBROUTINE equare_telsem2(DLAT,NCELLS,FIRSTCELL)
    !======computes the number of cells in a lattitude band
    !======and the first cellnumber of a lattitude band
    IMPLICIT none
    ! EQUAL-AREA COMPUTATIONS                                           
    ! I  TOTCEL             TOTAL NUMBER OF E.A. BOXES                  
    !  NCELLS(720)        NUMBER OF E.A. BOXES PER LAT ZONE           
    !  TOCELL(1440,720)   BOX NUMBER OF E.A. BOX (1-659600)         
    real(kind=dbl), INTENT (IN)      :: dlat
    integer(kind=long) :: NCELLS(:)
    integer(kind=long) :: FIRSTCELL(:)
    integer(kind=long) :: maxlat,maxlon
    integer(kind=long) :: TOTCEL
    integer(kind=long), ALLOCATABLE   :: TOCELL(:,:)
    integer(kind=long) :: maxlt2,lat,icellr,lat1,lat2,numcel,numcls,lon
    DOUBLE PRECISION    :: rcells,rearth,pi,aezon,hezon,aecell,xlatb
    DOUBLE PRECISION    :: rlatb,rlate,xlate,htb,hte,htzone,rcelat,azone
    !DOUBLE PRECISION    :: dlongr,acell,asq,twopi,halfpi,rcellr
    integer(kind=long) :: I

    maxlat=FLOOR(180./DLAT)
    maxlon=FLOOR(360./DLAT)

    ALLOCATE(TOCELL(maxlon,maxlat))

    !COMMON /EQUCOM/TOCEll,NCELLS
    REARTH = 6371.2d0                                                   
    PI     = 2.0d0 * ASIN(1.d0)   
    !HALFPI=PI/2.d0
    !TWOPI=2.d0*PI                                                       
    RCELAT=(DLAT*PI)/180.d0
    TOTCEL=0d0                                                          
    ! CALCULATE HEIGHT AND AREA OF EQUATORIAL ZONE                          
    HEZON=REARTH*SIN(RCELAT)                                          
    AEZON=2.d0*PI*REARTH*HEZON                                          
    ! CALCULATE AREA OF EQUATORIAL CELL                                    
    AECELL=(AEZON*DLAT)/360.d0                                          
    ! print*,'aecell',aecell
    ! COMPUTE LONGITUDE ZONES FOR EACH LATITUDE ZONE                       
    MAXLT2=MAXLAT/2                                                   
    DO LAT=1,MAXLT2                                               
       XLATB=(LAT-1)*DLAT                                                
       XLATE=XLATB+DLAT                                                  
       RLATB=(2.0d0*PI*XLATB)/360.0d0
       RLATE=(2.0d0*PI*XLATE)/360.0d0
       !CALCULATE HEIGHTS OF LATB,LATE,ZONE                                  
       HTB=REARTH*SIN(RLATB)                                             
       HTE=REARTH*SIN(RLATE)                                             
       HTZONE=HTE-HTB                                                    
       AZONE=2.d0*PI*REARTH*HTZONE                                         
       !CALCULATE NUMBER OF CELLS 
       RCELLS=AZONE/AECELL                                               
       ICELLR=FLOOR(RCELLS+.50d0)                                                
       !AUGMENT TOTAL # GRID CELLS (BOTH HEMISPHERES)                     
       TOTCEL=TOTCEL+2*ICELLR                                            
       !RCELLR=ICELLR                                                     
       !DLONGR=360.0d0/RCELLR                                                
       !ACELL=AZONE/RCELLR                                                
       !ASQ=AZONE/MAXLON                                                  
       !CREATE TABLE OF LONGITUDES
       LAT1=LAT+MAXLT2                                                   
       LAT2=MAXLT2+1-LAT                                                 
       NCELLS(LAT1)=ICELLR                                               
       NCELLS(LAT2)=ICELLR                                               
    END DO
    NUMCEL = 0                                                     
    ! THROUGH EACH LAT ZONE                                        
    DO LAT=1,MAXLAT                                              
       NUMCLS = NCELLS(LAT)                                           
       ! if (lat.eq.9) print*,lat,numcls,360./numcls
       ! LOOP THROUGH EACH LON FOR THIS LAT ZONE                           
       DO LON=1,NUMCLS                                           
          NUMCEL = NUMCEL + 1                                         
          ! FILL TOCELL ARRAY WITH STARTING CELL NUMBER
          TOCELL(LON,LAT) = NUMCEL
       END DO
    END DO
    DEALLOCATE(TOCELL)

    ! search for the first cellnumber in each lat band
    FIRSTCELL(1)=1
    DO I=2,maxlat
       FIRSTCELL(I)=FIRSTCELL(I-1)+NCELLS(I-1)
    END DO
  END SUBROUTINE equare_telsem2
! 
! 
!   !------------------------------------------------------------------
  FUNCTION calc_cellnum(lat,lon)
    !======computes the cellnumber given the lat and long
    !======using the ncells included in the atlas
    IMPLICIT none
    real(kind=sgl), INTENT (IN) :: lat, lon
!    type(atlas_emis_mw), INTENT (IN) :: atlas
    integer(kind=long) :: calc_cellnum
    integer(kind=long) :: ilat,ilon, i
    ! search for the cellnum in the older file
    calc_cellnum=0
    ilat=int((lat+90.)/atlas_dlat)+1
    ilon=int(lon/(360./atlas_ncells(ilat)))+1
    do i=1,ilat-1
       calc_cellnum=atlas_ncells(i)+calc_cellnum
    end do
    calc_cellnum=calc_cellnum+ilon
  END FUNCTION calc_cellnum
! 
!   !------------------------------------------------------------------
  SUBROUTINE get_coordinates(cellnum,atlas_dlat,atlas_firstcell,atlas_ncells,lat,lon)
    !======computes lat and lon given the cellnumber
    IMPLICIT none
    integer(kind=long), INTENT (IN) :: cellnum
    ! resolution of the atlas (equal-area)
    real(kind=dbl) :: atlas_dlat
    ! number of cells per lat band
    integer(kind=long), dimension(720) :: atlas_ncells
    ! the first cellnumber of lat band
    integer(kind=long), dimension(720) :: atlas_firstcell
!    type(atlas_emis_mw), INTENT (IN) :: atlas
    real(kind=sgl), INTENT (OUT) :: lat, lon
    !integer(kind=long) :: ilat,ilon
    integer(kind=long) :: i
    real(kind=dbl) :: res_lat ! latitude resolution
    integer(kind=long) :: index_lat_max,index_lat,index_lon

    res_lat = atlas_dlat 

    index_lat_max  = int(180/res_lat)
 
    IF (cellnum >= atlas_firstcell(index_lat_max)) THEN
       index_lat = index_lat_max
       lat =(index_lat - 0.5)*res_lat - 90
       index_lon = cellnum - atlas_firstcell(index_lat_max)+1 
       lon = (index_lon - 0.5)*(360.0/atlas_ncells(index_lat))
    ELSE
       DO i=1,index_lat_max-1 
          IF ( (cellnum>=atlas_firstcell(i)) .AND. (cellnum<atlas_firstcell(i+1)) ) THEN
	     index_lat = i
	     lat = (index_lat - 0.5)*res_lat- 90
	     !	cout << i << "  " << *lat <<endl;
	     index_lon = cellnum - atlas_firstcell(i)+1
	     lon = (index_lon - 0.5)*(360.0/atlas_ncells(index_lat))
          END IF
       END DO
    END IF
  END SUBROUTINE get_coordinates
! 
!   !------------------------------------------------------------------
!   SUBROUTINE emis_interp_ind_sing(lat,lon,theta,freq,&
!   atlas_correspondance,atlas_emis,atlas_emis_err,atlas_class1,atlas_class2,ev,eh,stdv,stdh,verb)
  SUBROUTINE emis_interp_ind_sing(errorstatus,lat,lon,theta,freq,ev,eh)
    !======interpolates emissivities 
    !          IND: individual cellnumber atlas
    !          SING: singular channel
    Implicit none     
    real(kind=sgl), INTENT (IN)    :: lat, lon
    real(kind=dbl), INTENT (IN)    :: theta, freq

    real(kind=dbl), INTENT (OUT)   :: ev, eh   !resultat de l'interpolation
!     integer(kind=long), INTENT (IN) :: verb
    integer :: verbosity = 0
    integer(kind=long) :: ipos
    integer(kind=long) :: cellnum
    real(kind=dbl)    :: ev_a(3),eh_a(3)      !emissivity in the atlas

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(22) :: nameOfRoutine = 'emis_interp_ind_sing'

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)  
    
    !Initialisations
    ev=0
    eh=0

    cellnum=calc_cellnum(lat,lon)
    ipos=atlas_correspondance(cellnum)
    IF (ipos>0) THEN
       ev_a(1)=atlas_emis(ipos,1)
       eh_a(1)=atlas_emis(ipos,2)
       ev_a(2)=atlas_emis(ipos,4)
       eh_a(2)=atlas_emis(ipos,5)
       ev_a(3)=atlas_emis(ipos,6)
       eh_a(3)=atlas_emis(ipos,7)
       CALL emis_interp(lat,lon,theta,freq,atlas_class1(ipos),atlas_class2(ipos),ev_a,eh_a,ev,eh)
    END IF
    IF (verbosity >= 3) THEN
       WRITE(0,*) 'Cellnum=',cellnum,' lat=',lat,' lon=',lon,' ipos=',ipos
    END IF
    
    errorstatus = err

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)    
  END SUBROUTINE emis_interp_ind_sing
! 
!   !------------------------------------------------------------------
  SUBROUTINE interp_freq2(emiss19,emiss37,emiss85,f,class2,emiss)  
    !======Linear interpolation of emissivity up to 85 GHz.
    !======Above 85 GHz, for most surfaces constant emissivity  
    !======except when emiss85>emiss37 for 3 classes of sea ice (class2 11 12 13) 
    !======and for surface water (class2 10) (ie when water like behavior)
    !======For most surface types, very limited frequency dependence or not enough 
    !======convincing evidence for a reliable frequency dependence.   
    IMPLICIT none
    real(kind=dbl), INTENT (IN) :: emiss19,emiss37,emiss85,f
    integer(kind=long), INTENT (IN) :: class2
    real(kind=dbl), INTENT (OUT) :: emiss
    real(kind=dbl) :: a, b, c,rapport43_32(4),rapport54_43(4)
     
    DATA rapport43_32 / 0.62, 0.37, 0.46, 0.63 / 
    DATA rapport54_43 / 0.30, 0.60, 0.47, 0.35 / 
    
    ! class2 1 to 5 idem as TELSEM
    ! class2 10=water
    ! class2 11 to 16 sea ice
    ! class2 17 to 22 snow and continental ice
    ! Frequency interpolation above 85 GHz only for classes 10 to 13 (water like classes)

    IF (f<=19.35) THEN
       a=1
       b=0
       c=0
       emiss = emiss19
    ELSE IF ((19.35<f).AND.(f<=37.)) THEN
       a=(37.-f  )/(37.-19.35)
       b=(f-19.35)/(37.-19.35)
       c=0
       emiss = a*emiss19+b*emiss37
    ELSE IF ((f>37.).AND.(f<85.5)) THEN
       a=0
       b=(85.5-f )/(85.5-37)
       c=(f-37   )/(85.5-37)
       emiss = b*emiss37+c*emiss85
    ELSE IF (85.5<=f) THEN 
        a=0
        b=0
        c=1  
        emiss = emiss85 
        IF ((class2>9).AND.(class2<14).AND.(emiss85>emiss37)) THEN
           IF (f<=150) THEN 
               emiss=emiss85+(f-85.5)*((emiss85-emiss37)/(85.5-37.))*rapport43_32(class2-9)
               IF (emiss>1) emiss=1
           END IF
           IF ((f>150).AND.(f<=190)) THEN
               emiss=emiss85+(150.-85.5)*((emiss85-emiss37)/(85.5-37.))*rapport43_32(class2-9)
               emiss=emiss+(f-150.)*((emiss-emiss85)/(150.-85.5))*rapport54_43(class2-9)
               IF (emiss>1) emiss=1     
           END IF
           IF (f>190) THEN
               emiss=emiss85+(150.-85.5)*((emiss85-emiss37)/(85.5-37.))*rapport43_32(class2-9)
               emiss=emiss+(190.-150.)*((emiss-emiss85)/(150.-85.5))*rapport54_43(class2-9)
               IF (emiss>1) emiss=1     
           END IF          
        END IF
    END IF
  END SUBROUTINE interp_freq2
! 
!  !------------------------------------------------------------------
   SUBROUTINE emis_interp(lat,lon,theta,freq,class1,class2,ev,eh,emiss_interp_v,emiss_interp_h)
    !======interpolation of emissivities for angle, freq
    IMPLICIT none     
    integer(kind=long), INTENT (IN) :: class1,class2
    real(kind=sgl), INTENT (IN)    :: lat, lon
    real(kind=dbl), INTENT (IN)    :: theta, freq, ev(3), eh(3)
    real(kind=dbl), INTENT (OUT)   :: emiss_interp_v, emiss_interp_h
    real(kind=dbl)    :: e0, theta0, theta53 , emiss_scal_v(3), emiss_scal_h(3)
    real(kind=dbl)    :: S1_v, S1_h, S2_v, S2_h, S_v, S_h, a0, a1, a2, a3, b0, b1
    real(kind=dbl)    :: b2, b3, em53_v, em53_h, emtheta_v, emtheta_h
    real(kind=dbl)    :: a0_k0(3,10),a0_k1(3,10),a0_k2(3,10)
    real(kind=dbl)    :: a0_eveh(3,10),a1_eveh(3,10),a2_eveh(3,10),a3_eveh(3,10)
    real(kind=dbl)    :: b0_eveh(3,10),b1_eveh(3,10),b2_eveh(3,10),b3_eveh(3,10)
    integer(kind=long) :: j
    !COMMON /EMISSIVITE/emiss_interp_v,emiss_interp_h
    data a0_k0/0.11509,0.091535,0.34796,0.10525,0.16627,0.24434, &
         & 0.29217,0.23809,0.28954,0.17516,0.19459,0.28697, &
         & 0.10521,0.12126,0.30278,0.18212,0.19625,0.14551, &
         & -0.19202,0.5411,0.03739,0.10292,0.5486,-0.058937, &
         & -0.022672,0.44492,-0.058448,-0.33894,-0.17621,0.14742/
    data a0_k1/0.61168,0.59095,0.7918,0.60271,0.69213,0.62218, &
         &  0.32728,0.34334,0.37062,0.51217,0.4491,0.50101, &
         & 0.48913,0.41932,0.29734,0.64474,0.30637,0.031107, &
         & 1.0405,0.17538,1.3215,0.61819,0.31298,1.7218, &
         & 0.87761,0.47583,1.2583,1.0959,0.92842,0.51033/
    data a0_k2/0.26726,0.32033,-0.14778,0.28547,0.13592,0.13193, &
         & 0.37178,0.41813,0.33875,0.30203,0.35479,0.20189, &
         & 0.40663,0.47493,0.40668,0.14811,0.52382,0.86634, &
         & 0.14286,0.27164,-0.37947,0.2737,0.12001,-0.67315, &
         & 0.13492,0.065463,-0.19316,0.24905,0.25475,0.34637/
    data a0_eveh/0.9592599869E+00,0.9565299749E+00,0.9511899948E+00, &
         & 0.9560700059E+00,0.9541199803E+00,0.9483199716E+00, &
         & 0.9461100101E+00,0.9439799786E+00,0.9387800097E+00, &
         & 0.9317600131E+00,0.9289000034E+00,0.9236800075E+00, &
         & 0.9208700061E+00,0.9190599918E+00,0.9105200171E+00, &
         & 0.9162799716E+00,0.8937299848E+00,0.8014699817E+00, &
         & 0.9570500255E+00,0.9213600159E+00,0.7893999815E+00, &
         & 0.9639400244E+00,0.9530599713E+00,0.8850200176E+00, &
         & 0.9685299993E+00,0.9622600079E+00,0.9118800163E+00, &
         & 0.8997200131E+00,0.9012699723E+00,0.9107499719E+00/
    data a1_eveh/0.3627802414E-07,-0.7778328204E-08,0.4396108011E-07, &
         & 0.2503205394E-06,0.1996262995E-06,0.2929977541E-06, &
         & 0.4190530660E-06,0.3655744649E-06,0.3519195673E-06, &
         & 0.5574374313E-06,0.5273076340E-06,0.5376484182E-06, &
         & 0.1026844529E-05,0.9679998811E-06,0.8616486866E-06, &
         & 0.3180800832E-06,0.2886778532E-06,0.2310362675E-06, &
         & -0.1118036366E-06,-0.1502856577E-06,0.4842232926E-07, &
         & -0.8410978580E-08,-0.3478669441E-07,0.2209441590E-06, &
         & 0.2485776633E-06,0.1800235907E-06,0.2510202251E-06, &
         & 0.2687000915E-06,0.1740325644E-06,0.3562134339E-06/
    data a2_eveh/0.3067140824E-05,0.2520012231E-05,0.4831396382E-05, &
         & 0.8213598448E-05,0.7378375358E-05,0.1022081960E-04, &
         & 0.1225889173E-04,0.1165553113E-04,0.1188659007E-04, &
         & 0.1693615741E-04,0.1648317448E-04,0.1715818144E-04, &
         & 0.2744720041E-04,0.2642072104E-04,0.2671847506E-04, &
         & 0.1349592094E-04,0.1261523357E-04,0.5447756394E-05, &
         & 0.2064244654E-05,0.1919016057E-06,0.5940860319E-06, &
         & 0.5334760772E-05,0.4130339221E-05,0.4104662821E-05, &
         & 0.6530796327E-05,0.5727014013E-05,0.7451782039E-05, &
         & 0.1071246970E-04,0.9539280654E-05,0.1034286015E-04/
    data a3_eveh/-0.2004991551E-07,-0.6895366056E-07, &
         & -0.2047409282E-06, &
         & -0.7322448425E-07,-0.1273002681E-06,-0.2729916844E-06, &
         & -0.9421125213E-07,-0.1683332300E-06,-0.2726891637E-06, &
         & -0.1317753799E-06,-0.2107972250E-06,-0.3556060904E-06, &
         & -0.1889465580E-06,-0.2757958271E-06,-0.4909850304E-06, &
         & 0.7339644004E-08,-0.4058669560E-06,-0.4146343997E-06, &
         & 0.6170279931E-07,-0.1998567996E-06,-0.4713119139E-07, &
         & -0.1361754887E-07,-0.1765622955E-06,-0.2348146637E-06, &
         & -0.3901189061E-07,-0.1305666189E-06,-0.1533838798E-06, &
         & -0.2679148992E-07,-0.4441960044E-07,-0.1815613899E-06/
    data b0_eveh/0.9592599869E+00,0.9565299749E+00,0.9511899948E+00, &
         & 0.9560700059E+00,0.9541199803E+00,0.9483199716E+00, &
         & 0.9461100101E+00,0.9439799786E+00,0.9387800097E+00, &
         & 0.9317600131E+00,0.9289000034E+00,0.9236800075E+00, &
         & 0.9208700061E+00,0.9190599918E+00,0.9105200171E+00, &
         & 0.9162799716E+00,0.8937299848E+00,0.8014699817E+00, &
         & 0.9570500255E+00,0.9213600159E+00,0.7893999815E+00, &
         & 0.9639400244E+00,0.9530599713E+00,0.8850200176E+00, &
         & 0.9685299993E+00,0.9622600079E+00,0.9118800163E+00, &
         & 0.8997200131E+00,0.9012699723E+00,0.9107499719E+00/
    data b1_eveh/0.3626608347E-07,-0.7786279177E-08,0.4393379172E-07, &
         & 0.2502746099E-06,0.1995944388E-06,0.2929554341E-06, &
         & 0.4189516289E-06,0.3655020180E-06,0.3518483140E-06, &
         & 0.5572838404E-06,0.5271903092E-06,0.5375342766E-06, &
         & 0.1026605219E-05,0.9677979733E-06,0.8614680951E-06, &
         & 0.3179358714E-06,0.2884899004E-06,0.2308632219E-06, &
         & -0.1118781370E-06,-0.1503948681E-06,0.4834672396E-07, &
         & -0.8455684153E-08,-0.3485171618E-07,0.2208606134E-06, &
         & 0.2485595019E-06,0.1799959364E-06,0.2509846695E-06, &
         & 0.2686167306E-06,0.1739760478E-06,0.3561317214E-06/
    data b2_eveh/0.3065537157E-05,0.2518960400E-05,0.4829731552E-05, &
         & 0.8209894986E-05,0.7375769655E-05,0.1021809931E-04, &
         & 0.1225203869E-04,0.1165053800E-04,0.1188218721E-04, &
         & 0.1692612022E-04,0.1647546378E-04,0.1715117833E-04, &
         & 0.2743142431E-04,0.2640772436E-04,0.2670711910E-04, &
         & 0.1348545720E-04,0.1260529825E-04,0.5439695997E-05, &
         & 0.2058213340E-05,0.1860650656E-06,0.5898303925E-06, &
         & 0.5330772183E-05,0.4126528893E-05,0.4100859314E-05, &
         & 0.6528573977E-05,0.5725009032E-05,0.7449450095E-05, &
         & 0.1070590315E-04,0.9534271157E-05,0.1033751869E-04/
    data b3_eveh/-0.1370247134E-06,-0.1436897747E-06, &
         & -0.2954870411E-06, &
         & -0.3118435643E-06,-0.2916583242E-06,-0.4311032171E-06, &
         & -0.5048401022E-06,-0.4662823869E-06,-0.5206445053E-06, &
         & -0.7210980471E-06,-0.6662896794E-06,-0.7548637200E-06, &
         & -0.1110204039E-05,-0.1030801400E-05,-0.1140921199E-05, &
         & -0.6330818110E-06,-0.9186441048E-06,-0.7947813856E-06, &
         & -0.3242539890E-06,-0.5027602583E-06,-0.2777987334E-06, &
         & -0.2747250676E-06,-0.3811997260E-06,-0.4102405455E-06, &
         & -0.1994112324E-06,-0.2555484855E-06,-0.2842682534E-06, &
         & -0.4413041665E-06,-0.3717419474E-06,-0.4975536854E-06/
    ! Interpolation en angle
    DO j = 1,3
       ! Calcul par regression multilineaire de la valeur e0 en theta=0°
       e0 = a0_k0(j,class1)+a0_k1(j,class1)*ev(j)+a0_k2(j,class1)*eh(j)
       ! Lecture des coefficients des polynomes ev et eh
       a0 = a0_eveh(j,class1)
       a1 = a1_eveh(j,class1)
       a2 = a2_eveh(j,class1)
       a3 = a3_eveh(j,class1)
       b0 = b0_eveh(j,class1)
       b1 = b1_eveh(j,class1)
       b2 = b2_eveh(j,class1)
       b3 = b3_eveh(j,class1)
       theta0 = 0.
       theta53 = 53.
       ! Polarisation verticale
       S1_v = ((theta-theta53)/(theta0-theta53)) * ((e0-a0)/a0)
       em53_v = a3*(theta53**3) + a2*(theta53**2) + a1*theta53 + a0
       S2_v =((theta-theta0)/(theta53-theta0))*((ev(j)-em53_v)/em53_v)
       S_v = 1 + S1_v + S2_v
       emtheta_v = a3*(theta**3) + a2*(theta**2) + a1*theta + a0
       emiss_scal_v(j) = S_v * emtheta_v     
       ! Polarisation horizontale
       S1_h = ((theta-theta53)/(theta0-theta53)) * ((e0-b0)/b0)
       em53_h = b3*(theta53**3) + b2*(theta53**2) + b1*theta53 + b0
       S2_h =((theta-theta0)/(theta53-theta0))*((eh(j)-em53_h)/em53_h)
       S_h = 1 + S1_h + S2_h
       emtheta_h = b3*(theta**3) + b2*(theta**2) + b1*theta + b0
       emiss_scal_h(j) = S_h * emtheta_h     
    END DO
    ! Interpolation en frequence
    !emiss_interp_v=interp_freq(emiss_scal_v(1),emiss_scal_v(2),emiss_scal_v(3),freq)
    !emiss_interp_h=interp_freq(emiss_scal_h(1),emiss_scal_h(2),emiss_scal_h(3),freq)  
    CALL interp_freq2(emiss_scal_v(1),emiss_scal_v(2),emiss_scal_v(3),freq,class2,emiss_interp_v)
    CALL interp_freq2(emiss_scal_h(1),emiss_scal_h(2),emiss_scal_h(3),freq,class2,emiss_interp_h)  
    ! Cas ev<eh: on fait la moyenne entre les deux
    IF (emiss_interp_v < emiss_interp_h) THEN
       emiss_interp_v = (emiss_interp_v + emiss_interp_h)/2.
       emiss_interp_h =  emiss_interp_v
    END IF
  END SUBROUTINE emis_interp


END MODULE mod_mwatlas_nt_bin


