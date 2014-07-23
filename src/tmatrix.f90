module tmatrix
  use kinds
  use constants, only: pi, c
  use settings, only: nummu, nstokes, verbose, active, passive, &
    tmatrix_db, tmatrix_db_path, radar_npol, radar_pol, liq_mod
  use rt_utilities, only: lobatto_quadrature
  use report_module

  implicit none

  contains
  
  subroutine calc_tmatrix(errorstatus,&
    frequency,&
    ref_index,&
    phase,&
    nbins,&
    dmax,&
    del_d,&
    ndens,&
    density,&
    as_ratio,&
    canting,&
    temp, &
    scatter_matrix,&
    extinct_matrix,&
    emis_vector,&
    back_spec)
    

      !  calculate the matrix and vectors, for spectrum
      !
      !   input:
      !       frequency       double  frequency [Hz]
      !       ref_index       double refractive index of pure ice or water
      !       phase           integer
      !       nbins           integer no of bins
      !       dmax            dbl (nbins)  dmax [m]
      !       del_d           width of dmax bins [m]    
      !       ndens           number density (not normed) [1/m³]
      !       density         dbl (nbins) density of softspheres
      !       as_ratio        double  aspect ratio
      !       canting        double  canting angle (deg) -> beta in tmatrix code
      !       temp           double temperature [K]
      !
      !   output:
      !       scatter_matrix  double  scattering matrix []
      !       extinct_matrix  double  extinction matrix []
      !       emis_vector     double  emission vector []
      !       back_spec       double  backscattering spectrum [nbins]

      use vars_index, only: i_p
      implicit none

      real(kind=dbl), intent(in) :: frequency
      complex(kind=dbl), intent(in) :: ref_index
      integer, intent(in) :: phase
      integer, intent(in) :: nbins
      real(kind=dbl), dimension(nbins), intent(in) :: dmax
      real(kind=dbl), dimension(nbins), intent(in) :: del_d
      real(kind=dbl), dimension(nbins), intent(in) :: ndens    
      real(kind=dbl), dimension(nbins), intent(in) :: density
      real(kind=dbl), dimension(nbins), intent(in) :: as_ratio
      real(kind=dbl), dimension(nbins), intent(in) :: canting
      real(kind=dbl), intent(in) :: temp

      real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
      real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
      real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector
      real(kind=dbl), intent(out), dimension(radar_npol,nbins) :: back_spec

      complex(kind=dbl) :: mMix
      complex(kind=dbl) :: eps_mix 
      complex(kind=ext) :: mindex
      real(kind=dbl) :: alpha
      real(kind=dbl) :: beta
      real(kind=dbl) :: axi
      real(kind=dbl) :: ndens_eff
      real(kind=dbl) :: del_d_eff
      integer :: azimuth_num
      integer :: azimuth0_num   
      integer :: ir
      integer :: n_lines, work1
      character(len=1)    :: work2
      character(1) :: quad
      character(12) ::db_file
      character(600) ::db_path
      character(5) ::format_str
      logical ::file_exists, file_OK
  
      real(kind=dbl), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix_part
      real(kind=dbl), dimension(nstokes,nstokes,nummu) :: extinct_matrix_part
      real(kind=dbl), dimension(nstokes,nummu) :: emis_vector_part
  
      real(kind=dbl), dimension(nstokes*nummu*nstokes*nummu*2) :: scatter_matrix_flat
      real(kind=dbl), dimension(nstokes*nstokes*nummu) :: extinct_matrix_flat
      real(kind=dbl), dimension(nstokes*nummu) :: emis_vector_flat

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=14) :: nameOfRoutine = 'tmatrix_clac'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)    
      
      if (verbose >= 4) print*,"frequency,ref_index,phase,nbins,dmax,del_d,ndens,density,as_ratio"
      if (verbose >= 4) print*,frequency,ref_index,phase,nbins,dmax,del_d,"N ",&
        ndens,"rho ", density,"AR ",as_ratio

      err = 0

      call assert_false(err,(isnan(frequency) .or. frequency < 0.d0),&
	  "nan or negative frequency")
      call assert_false(err,(isnan(real(ref_index)) .or. isnan(imag(ref_index))),&
	  "nan ref_index")
      call assert_true(err,(nbins>0),&
	  "nbins>0")
      call assert_false(err,(any(isnan(dmax)) .or. any(dmax <= 0.d0)),&
	  "nan or negative dmax")
      call assert_false(err,(any(isnan(del_d)) .or. any(del_d <= 0.d0)),&
	  "nan del_d")
      call assert_false(err,(any(isnan(ndens)) .or. any(ndens < 0.d0)),&
	  "nan or negative ndens")
      call assert_true(err,SUM(ndens)>0,&
          "sum(ndens) must be greater zero")    
      call assert_false(err,(any(isnan(density)) .or. any(density < 0.d0)),&
	  "nan or negative density")
      call assert_false(err,any(isnan(as_ratio)) .or. any(as_ratio < 0.d0),&
          "nan or negative as_ratio")
      call assert_false(err,any(isnan(canting)) .or. any(canting < 0.d0),&
          "nan or negative canting")
      if (err > 0) then
	  errorstatus = fatal
	  msg = "assertation error"
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      end if   

      !T Matrix settings
      alpha = 0.0_dbl    ! orientation of the particle [°]
      azimuth_num = 30
      azimuth0_num = 1   
      quad ="L" !quadratur
	
    !initialize
      back_spec(:,:) = 0.d0
      scatter_matrix = 0.d0
      extinct_matrix = 0.d0
      emis_vector = 0.d0

    do ir = 1, nbins


        ndens_eff = ndens(ir)
        del_d_eff = del_d(ir)

      beta = canting(ir)

      !in case we have no hydrometeors, we need no tmatrix calculations!
      if (ndens_eff == 0.d0) then
        if (verbose >= 4) print*, "Skipped iteration", ir, "because ndens_eff", ndens_eff
        CYCLE
      end if
      call assert_true(err,ndens_eff>0,&
          "nan or negative ndens_eff")
      call assert_true(err,del_d_eff>0,&
          "nan or negative del_d_eff")
      call assert_true(err,density(ir)>0,&
          "nan or negative density(ir)")
      if (err > 0) then
          errorstatus = fatal
          msg = "assertation error"
          call report(errorstatus, msg, nameOfRoutine)
          return
      end if   
      if (phase == -1 .and. density(ir) /= 917.d0) then
          mMix = eps_mix((1.d0,0.d0),ref_index,density(ir))
      else
                mMix = ref_index
      end if      
      mindex =conjg(mMix) !different convention

      !we want the volume equivalent radius
      if (as_ratio(ir) <= 1) then
        !oblate, with axis of rotation vertically
        axi = 0.5_dbl*dmax(ir)*as_ratio(ir)**(1.0_dbl/3.0_dbl) 
      else 
        !prolate, with axis of rotation vertically
        axi = 0.5_dbl*dmax(ir)/as_ratio(ir)**(2.0_dbl/3.0_dbl) 
      end if

      if (tmatrix_db == "none") then
        call calc_single_tmatrix(err,quad,nummu,frequency,mindex,axi, nstokes,&
            as_ratio(ir), alpha, beta, azimuth_num, azimuth0_num,&
            scatter_matrix_part,extinct_matrix_part,emis_vector_part)
        if (err /= 0) then
            msg = 'error in calc_single_tmatrix!'
            call report(err, msg, nameOfRoutine)
            errorstatus = err
            return
        end if          
      else if (tmatrix_db == "file") then
        
        db_path =""
!         write(db_path,'(A4,A6,A1,4(A6,I3.3),A6,E12.6,2(A6,ES36.30),4(A6,ES14.8),A1)'),&
!                   "/v01","/quad_", quad, &
!                   "/numu_",nummu,"/azno_",azimuth_num, "/a0no_", azimuth0_num, "/nsto_",nstokes,&
!                   "/freq_",frequency, &
!                   "/min1_", REAL(mindex), "/min2_",IMAG(mindex), &
!                   "/axxi_",axi, "/asra_",as_ratio(ir), "/alph_",alpha, "/beta_",beta,"/"

        if (phase == 1) then !ref_index name depends on phase
          write(db_path,'(A4,A6,A1,4(A6,I3.3),A6,ES12.6,A6,SP,I3.2,SS,2(A6,ES10.4),A6,A3,4(A6,ES10.4),A1)'),&
                  "/v01","/quad_", quad, &
                  "/numu_",nummu,"/azno_",azimuth_num, "/a0no_", azimuth0_num, "/nsto_",nstokes,&
                  "/freq_",frequency, &
                  "/phas_",phase, "/temp_",temp, "/dens_",density(ir), "/refi_",liq_mod,&
                  "/diam_",axi, "/asra_",as_ratio(ir), "/alph_",alpha, "/beta_",beta,"/"
        else if (phase == -1) then
          write(db_path,'(A4,A6,A1,4(A6,I3.3),A6,ES12.6,A6,SP,I3.2,SS,2(A6,ES10.4),A6,A3,4(A6,ES10.4),A1)'),&
                  "/v01","/quad_", quad, &
                  "/numu_",nummu,"/azno_",azimuth_num, "/a0no_", azimuth0_num, "/nsto_",nstokes,&
                  "/freq_",frequency, &
                  "/phas_",phase, "/temp_",temp, "/dens_",density(ir), "/refi_","Mae",&
                  "/diam_",axi, "/asra_",as_ratio(ir), "/alph_",alpha, "/beta_",beta,"/"
        else
            errorstatus = fatal
            msg = "Do not understand phase"
            call report(errorstatus, msg, nameOfRoutine)
            return
        end if  

        db_file = "spheroid.dat"
        INQUIRE(FILE=TRIM(tmatrix_db_path)//TRIM(db_path)//TRIM(db_file), EXIST=file_exists)

        !check whether file is not empty
        if (file_exists) then
          open(112,file=TRIM(tmatrix_db_path)//TRIM(db_path)//TRIM(db_file),action="READ")
          n_lines = 0
          do
            read(112,*,IOSTAT=work1)  work2
            if (work1 /= 0) exit
            n_lines = n_lines + 1
          end do
          rewind(112)
          if (n_lines == 3) then 
            file_OK = .true.
          else
            file_OK = .false.
          end if

        end if
        if (file_exists .and. file_OK) then

          if (verbose > 0) print * , TRIM(db_path)//TRIM(db_file), " exists, opening"

          open(112,file=TRIM(tmatrix_db_path)//TRIM(db_path)//TRIM(db_file),action="READ")

          write(format_str,"(I5.5)") SHAPE(scatter_matrix_flat)
          read(112,"("//format_str//"(ES25.17, 2x))")scatter_matrix_flat

          write(format_str,"(I5.5)") SHAPE(extinct_matrix_flat)
          read(112,"("//format_str//"(ES25.17, 2x))")extinct_matrix_flat

          write(format_str,"(I5.5)") SHAPE(emis_vector_flat)
          read(112,"("//format_str//"(ES25.17, 2x))")emis_vector_flat

          close(112)
          if (verbose > 1) print * , TRIM(db_path)//TRIM(db_file), " closed"

          err = 0
          call assert_true(err,PRODUCT(SHAPE(scatter_matrix_part)) == PRODUCT(SHAPE(scatter_matrix_flat)),&
              "shape of scatter_matrix_flat does not match")
          call assert_true(err,PRODUCT(SHAPE(extinct_matrix_part)) == PRODUCT(SHAPE(extinct_matrix_flat)),&
              "shape of extinct_matrix_flat does not match")
          call assert_true(err,PRODUCT(SHAPE(emis_vector_part)) == PRODUCT(SHAPE(emis_vector_flat)),&
              "shape of emis_vector_flat does not match")
          if (err > 0) then
              errorstatus = fatal
              msg = "assertation error"
              call report(errorstatus, msg, nameOfRoutine)
              return
          end if   

          scatter_matrix_part =reshape(scatter_matrix_flat,SHAPE(scatter_matrix_part))
          extinct_matrix_part =reshape(extinct_matrix_flat,SHAPE(extinct_matrix_part))
          emis_vector_part =reshape(emis_vector_flat,SHAPE(emis_vector_part))

          ! for debugging only: calculate scatter matrix and compare with file
          if (verbose .gt. 20) then
            call calc_single_tmatrix(err,quad,nummu,frequency,mindex,axi, nstokes,&
                as_ratio(ir), alpha, beta, azimuth_num, azimuth0_num,&
                scatter_matrix_part,extinct_matrix_part,emis_vector_part)
            if (err /= 0) then
                msg = 'error in calc_single_tmatrix!'
                call report(err, msg, nameOfRoutine)
                errorstatus = err
                return
            end if 
            print*, MAXVAL((ABS(reshape(emis_vector_flat,SHAPE(emis_vector_part)) &
                - emis_vector_part) /emis_vector_part))
            print*, MAXVAL((ABS(reshape(extinct_matrix_flat,SHAPE(extinct_matrix_part)) &
                - extinct_matrix_part) /extinct_matrix_part))
            print*, MAXVAL((ABS(reshape(scatter_matrix_flat,SHAPE(scatter_matrix_part)) &
                - scatter_matrix_part) /scatter_matrix_part))

          end if

        else
          !file does not exist or is not OK
          
          if (file_exists) close(112)

          if (verbose > 0) print * , TRIM(tmatrix_db_path)//TRIM(db_path)//TRIM(db_file), " NOT FOUND. calculating..."
          CALL EXECUTE_COMMAND_LINE("mkdir -p "//TRIM(tmatrix_db_path)//TRIM(db_path))

          call calc_single_tmatrix(err,quad,nummu,frequency,mindex,axi, nstokes,&
              as_ratio(ir), alpha, beta, azimuth_num, azimuth0_num,&
              scatter_matrix_part,extinct_matrix_part,emis_vector_part)
          if (err /= 0) then
              msg = 'error in calc_single_tmatrix!'
              call report(err, msg, nameOfRoutine)
              errorstatus = err
              return
          end if 

          scatter_matrix_flat =reshape(scatter_matrix_part,SHAPE(scatter_matrix_flat))
          extinct_matrix_flat =reshape(extinct_matrix_part,SHAPE(extinct_matrix_flat))
          emis_vector_flat =reshape(emis_vector_part,SHAPE(emis_vector_flat))
          if (verbose > 1) print * ,TRIM(db_path)//TRIM(db_file), " open..."

          open(113,file=TRIM(tmatrix_db_path)//TRIM(db_path)//TRIM(db_file),ACTION="WRITE")

          write(format_str,"(I5.5)") SHAPE(scatter_matrix_flat)
          write(113,"("//format_str//"(ES25.17, 2x))")scatter_matrix_flat

          write(format_str,"(I5.5)") SHAPE(extinct_matrix_flat)
          write(113,"("//format_str//"(ES25.17, 2x))")extinct_matrix_flat

          write(format_str,"(I5.5)") SHAPE(emis_vector_flat)
          write(113,"("//format_str//"(ES25.17, 2x))")emis_vector_flat

          close(113)
          if (verbose > 1) print * , TRIM(db_path)//TRIM(db_file), " closed"

        end if

      else
            msg = 'do not understand tmatrix_db: '//tmatrix_db
            call report(err, msg, nameOfRoutine)
            errorstatus = err
            return
      end if

      do i_p= 1, radar_npol
        if (radar_pol(i_p) == "NN") then
          !scatter_matrix(A,B;C;D;E) backscattering is M11 of Mueller or Scattering Matrix (A;C=1), in quadrature 2 (E) first 16 (B) is 180deg (upwelling), 2nd 16 (D) 0deg (downwelling). this definition is lokkiing from BELOW, sc
          back_spec(i_p,ir) = 4*pi*ndens_eff*scatter_matrix_part(1,16,1,16,2)
      else if (radar_pol(i_p) == "HH") then
        !1.Vivekanandan, J., Adams, W. M. & Bringi, V. N. Rigorous Approach to Polarimetric Radar Modeling of Hydrometeor Orientation Distributions. Journal of Applied Meteorology 30, 1053–1063 (1991).
        back_spec(i_p,ir) = + scatter_matrix_part(1,16,1,16,2) &
                            - scatter_matrix_part(1,16,2,16,2) & 
                            - scatter_matrix_part(2,16,1,16,2) & 
                            + scatter_matrix_part(2,16,2,16,2) 
      else if (radar_pol(i_p) == "VV") then
        back_spec(i_p,ir) = + scatter_matrix_part(1,16,1,16,2) &
                            + scatter_matrix_part(1,16,2,16,2) & 
                            + scatter_matrix_part(2,16,1,16,2) & 
                            + scatter_matrix_part(2,16,2,16,2) 
!       else if (radar_pol(i_p) == "HV") then
!         back_spec(i_p,ir) = + scatter_matrix_part(1,16,1,16,2) &
!                             - scatter_matrix_part(1,16,2,16,2) & 
!                             + scatter_matrix_part(2,16,1,16,2) & 
!                             - scatter_matrix_part(2,16,2,16,2) 
      else
          msg = 'do not understand radar_pol(i_p): '//radar_pol(i_p)
          err = fatal
          call report(err, msg, nameOfRoutine)
          errorstatus = err
          return
        end if
      end do 
      scatter_matrix = scatter_matrix + scatter_matrix_part * ndens_eff * del_d_eff
      extinct_matrix = extinct_matrix + extinct_matrix_part * ndens_eff * del_d_eff
      emis_vector = emis_vector + emis_vector_part * ndens_eff * del_d_eff
    
    end do !nbins

    call assert_false(err,any(isnan(scatter_matrix)),&
        "nan in scatter matrix")
    call assert_false(err,any(isnan(extinct_matrix)),&
        "nan in extinct_matrix")
    call assert_false(err,any(isnan(emis_vector)),&
        "nan in emis_vector")
    call assert_false(err,any(isnan(back_spec)),&
        "nan in back_spec")	  
    if (err > 0) then
        errorstatus = fatal
        msg = "assertation error"
        call report(errorstatus, msg, nameOfRoutine)
        return
    end if   	

    errorstatus = err
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine) 
    return
    
      
  end subroutine calc_tmatrix


  subroutine calc_single_tmatrix(errorstatus,quad,qua_num,frequency,ref_index,&
      axi, nstokes,&
      as_ratio, alpha, beta, azimuth_num, azimuth0_num,&
      scatter_matrix,extinct_matrix,emis_vector)


      !  calculate the matrix and vectors, for a single particle with a single orientation
      !
      !   input:
      !       quad            char    name of quadrature
      !       quad_num        int     number of quadrature angles
      !       frequency       double  frequency [Hz]
      !       ref_index       complex refractive index (positve im part)
      !       axi             double  equivalent sphere radius [m]
      !       nstokes         int     number of stokes components
      !       as_ratio        double  aspect ratio
      !       alpha           double  orientation of the particle [°]
      !       beta            double  orientation of the particle [°]
      !       azimuth_num     int     number of azimuth angles
      !       azimuth0_num    int     number of azimuth angles
      !
      !   output:
      !       scatter_matrix  double  scattering matrix []
      !       extinct_matrix  double  extinction matrix []
      !       emis_vector     double  emission vector []



      implicit none


      character(1), intent(in) :: quad
      integer, intent(in) :: qua_num
      real(kind=dbl), intent(in) :: frequency
      complex(kind=ext) :: ref_index
      real(kind=dbl), intent(in) :: axi
      integer, intent(in) :: nstokes
      real(kind=dbl), intent(in) :: as_ratio
      real(kind=dbl), intent(in) :: alpha
      real(kind=dbl), intent(in) :: beta
      integer, intent(in) :: azimuth_num
      integer, intent(in) :: azimuth0_num
      real(kind=dbl), intent(out), dimension(nstokes,qua_num,nstokes,qua_num,2) :: scatter_matrix
      real(kind=dbl), intent(out), dimension(nstokes,nstokes,qua_num) :: extinct_matrix
      real(kind=dbl), intent(out), dimension(nstokes,qua_num) :: emis_vector

      real(kind=dbl) :: wave_num
      real(kind=dbl) ::thet0, thet, phi, phi0,&
	  rat, phi_weights, phi0_weights
      integer ::  m, n, &
	ii,jj,kk,ll, kkk1

      !       parameter(qua_num=8, nstokes = 4)


      complex(kind=dbl) :: s11, s12, s21, s22
      real(kind=dbl) :: qua_angle(qua_num), qua_weights(qua_num)
      real(kind=dbl) :: &
	  scatt_matrix_tmp1_11, scatt_matrix_tmp1_12,&
	  scatt_matrix_tmp1_21, scatt_matrix_tmp1_22
      real(kind=dbl) ::&
	  emis_vector_tmp1_11(2*qua_num), emis_vector_tmp1_12(2*qua_num),&
	  emis_vector_tmp2_11,emis_vector_tmp2_12
      real(kind=dbl) :: thet0_weights, thet_weights

      integer :: nmax, np, qua_start

      real(kind=ext) :: mrr, mri, lam
      ! some factors that stay constant during calculations
      real(kind=dbl) :: fact_sca
      complex(kind=dbl) :: fact_ext

      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=80) :: msg
      character(len=30) :: nameOfRoutine = 'tmatrix_clac_single'

	if (verbose >= 3) call report(info,'Start of ', nameOfRoutine) 
	
      if (verbose >= 5) print*, "quad,qua_num,frequency,ref_index,axi, nstokes,as_ratio, alpha, beta, azimuth_num, azimuth0_num"
      if (verbose >= 5) print*, quad,qua_num,frequency,ref_index,axi, nstokes,as_ratio, alpha, beta, azimuth_num, azimuth0_num

      np = -1
      rat = 1.e0

      
      LAM = c/(frequency)

      wave_num = 2.0_dbl*pi/LAM
      
      err = 0
      call assert_true(err,as_ratio> 0.d0,"nan or negative in as_ratio")    
      call assert_true(err,axi> 0.d0,"nan or negative in axi")    
      call assert_true(err,frequency> 0.d0,"nan or negative in frequency")   
      call assert_true(err,(wave_num > 0.d0),"nan or <= 0 in wave-num")
      call assert_true(err,(LAM > 0.d0),"nan or <= 0 in LAM")
      call assert_true(err,(AS_RATIO > 0.d0),"nan or <= 0 in AS_RATIO")

      if (err > 0) then
	  errorstatus = fatal
	  msg = "assertation error"
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      end if   	
      
  !     axi = axi!*1.e6
      mrr = REAL(ref_index)
      mri = abs(IMAG(ref_index))


      if (active .eqv. .true. .and. passive .eqv. .false.) then
        qua_start = 16
      else
        qua_start = 1
      end if

      if (verbose >= 5) print*,"lam, mrr,mri, AXI, AS_RATIO, RAT, NP"
      if (verbose >= 5) print*,lam, mrr,mri, AXI, AS_RATIO, RAT, NP
  ! call the tmatrix routine amplq -> fills common block /TMAT/
      call tmatrix_amplq(lam, mrr,mri, AXI, AS_RATIO, RAT, NP,nmax)

      extinct_matrix = 0.d0
      scatter_matrix = 0.d0
      emis_vector = 0.d0
      emis_vector_tmp2_11 = 0.d0
      emis_vector_tmp2_12 = 0.d0
      ! if the particle is rotationally-symmetric, reduce calculation time for orientation-averaging
      ! if not, do orientation averaging for incident and scatterred directions

      !      write(*,*)wave_num

      fact_sca = 0.5e0/(wave_num**2)
      fact_ext = 2.0e0*pi*cmplx(0.,1.)/wave_num**2.e0
      ! calculate the quadrature angle, number and weight according to quadrature method
      ! subroutine lobatto_quadrature and gauss_legendre_quadrature in the file named 'refractive_index.f'
      if (quad(1:1).eq.'l'.or.quad(1:1).eq.'L') then
	  call lobatto_quadrature(qua_num,qua_angle,qua_weights)
      else
	  errorstatus = fatal
	  msg = "did not understand 'quad'="//quad
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      end if 		  

      ! for each quadrature angle
      ii = 1
      do 1241 jj = qua_start, qua_num
	  thet0=acos(qua_angle(jj)*(-1.)**(real(ii)-1))*180.d0/pi
	  thet0_weights = qua_weights(jj)
	  if(thet0.gt.179.9999)thet0=180.0d0
	  ! initializing the emis vector summation
	  emis_vector_tmp1_11 = 0.d0
	  emis_vector_tmp1_12 = 0.d0

	  do 1242 kk = 1, 2
	      kkk1 = kk
	      do 1243 ll = qua_start, qua_num
		  thet=acos(qua_angle(ll)*(-1.)**(real(kk)-1))*180.d0/pi
		  thet_weights=qua_weights(ll)
		  if(thet.gt.179.9999)thet=180.d0

		  do 1244 m = 1, azimuth0_num ! 1
		      phi0 = 360.0d0/(real(azimuth0_num))*(real(m)-1.d0)
		      phi0_weights = 1.d0/360.d0*(360.d0/azimuth0_num)
		      !		     if(azimuth0_num.eq.1)phi0 = 0.0

		      scatt_matrix_tmp1_11 = 0.d0
		      scatt_matrix_tmp1_12 = 0.d0

		      scatt_matrix_tmp1_21 = 0.d0
		      scatt_matrix_tmp1_22 = 0.d0

		      do 1245 n = 1, azimuth_num ! 30
			  phi = 360.d0/real(azimuth_num)*(real(n)-1.d0)
			  phi_weights = 1.d0/360.d0*(360.d0/azimuth_num)

		  CALL tmatrix_AMPL(NMAX,dble(LAM),THET0,THET,PHI0,PHI,ALPHA,BETA,&
		    S11,S12,S21,S22)

			  s11 = s11*wave_num
			  s12 = s12*wave_num
			  s21 = s21*wave_num
			  s22 = s22*wave_num

			  scatt_matrix_tmp1_11 = scatt_matrix_tmp1_11 + (fact_sca*&
			      (s11*dconjg(s11)+s12*dconjg(s12)+s21*dconjg(s21)+s22*dconjg(s22)))*phi_weights

			  scatt_matrix_tmp1_12 = scatt_matrix_tmp1_12 + (fact_sca*&
			      (s11*dconjg(s11)-s12*dconjg(s12)+s21*dconjg(s21)-s22*dconjg(s22)))*phi_weights

			  scatt_matrix_tmp1_21 = scatt_matrix_tmp1_21 + (fact_sca*&
			      (s11*dconjg(s11)+s12*dconjg(s12)-s21*dconjg(s21)-s22*dconjg(s22)))*phi_weights

			  scatt_matrix_tmp1_22 = scatt_matrix_tmp1_22 + (fact_sca*&
			      (s11*dconjg(s11)-s12*dconjg(s12)-s21*dconjg(s21)+s22*dconjg(s22)))*phi_weights

			  if (phi0 .eq. phi .and. thet0 .eq. thet) then
			      extinct_matrix(1,1,jj) = extinct_matrix(1,1,jj)+phi0_weights*(-real((s11 + s22)*fact_ext))
			      extinct_matrix(1,2,jj) = extinct_matrix(1,2,jj)+phi0_weights*(-real((s11 - s22)*fact_ext))
			      extinct_matrix(2,1,jj) = extinct_matrix(2,1,jj)+phi0_weights*(-real((s11 - s22)*fact_ext))
			      extinct_matrix(2,2,jj) = extinct_matrix(2,2,jj)+phi0_weights*(-real((s11 + s22)*fact_ext))
			  end if
  1245                    continue   ! phi

  
		      scatter_matrix(1,ll,1,jj,kkk1) = scatter_matrix(1,ll,1,jj,kkk1) + scatt_matrix_tmp1_11*phi0_weights
		      scatter_matrix(1,ll,2,jj,kkk1) = scatter_matrix(1,ll,2,jj,kkk1) + scatt_matrix_tmp1_12*phi0_weights
		      scatter_matrix(2,ll,1,jj,kkk1) = scatter_matrix(2,ll,1,jj,kkk1) + scatt_matrix_tmp1_21*phi0_weights
		      scatter_matrix(2,ll,2,jj,kkk1) = scatter_matrix(2,ll,2,jj,kkk1) + scatt_matrix_tmp1_22*phi0_weights

  1244                continue  ! phi0

		  ! calculate the summation of the scattering matrix in the whole sphere
		  emis_vector_tmp1_11(ll+(kk-1)*qua_num) = scatter_matrix(1,ll,1,jj,kkk1)*thet_weights*2.*pi
		  emis_vector_tmp1_12(ll+(kk-1)*qua_num) = scatter_matrix(1,ll,2,jj,kkk1)*thet_weights*2.*pi

  1243            continue ! thet ll
  1242        continue

	  emis_vector(1,jj) = extinct_matrix(1,1,jj) - sum(emis_vector_tmp1_11)
	  emis_vector(2,jj) = extinct_matrix(1,2,jj) - sum(emis_vector_tmp1_12)

  1241    continue ! thet0 jj


      errorstatus = err
      if (verbose >= 3) call report(info,'End of ', nameOfRoutine) 
      return

  end subroutine calc_single_tmatrix

end module tmatrix