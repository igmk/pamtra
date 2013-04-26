subroutine tmatrix_sql(freq,t,as_ratio,diameter,particle_mass,ptype,scatter_matrix,extinct_matrix,emis_vector)

  use kinds
  use report_module
  use constants, only: pi, c, im
  use settings, only : nummu, nstokes
  use report_module
  use sqlite
  use sql_tools

  implicit none

  real(kind=dbl), intent(in) :: freq !in Hz
  real(kind=dbl), intent(in) :: t
  real(kind=dbl), intent(in) :: as_ratio
  real(kind=dbl), intent(in) :: particle_mass
  real(kind=dbl), intent(in) :: diameter
  character(len=4), intent(in) :: ptype

  real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector

  real(kind=dbl) :: freq_round
  real(kind=dbl) :: t_round
  real(kind=dbl) :: as_ratio_round
  real(kind=dbl) :: particle_mass_round
  real(kind=dbl) :: diameter_round
  real(kind=dbl) :: particle_mass_magnitude
  real(kind=dbl) :: diameter_magnitude

character(len=200) :: fname
character(len=20) :: par_name

  logical :: found

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'tmatrix_sql'

  interface
    subroutine tmatrix_refIndex(freq,t,as_ratio,diameter,particle_mass,ptype,scatter_matrix,extinct_matrix,emis_vector)
      use kinds
      use settings, only : nummu, nstokes
      implicit none

      real(kind=dbl), intent(in) :: freq
      real(kind=dbl), intent(in) :: t
      real(kind=dbl), intent(in) :: as_ratio
      real(kind=dbl), intent(in) :: particle_mass
      character(len=4), intent(in) :: ptype
      real(kind=dbl), intent(in) :: diameter

      real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
      real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
      real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector
    end subroutine tmatrix_refIndex
  end interface

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

fname="/work/mmaahn/test.sqlite"
par_name="test"

call sql_open_table(fname,par_name)




!   if (verbose >= 4) print*, "freq,t,as_ratio,diameter, particle_mass"
!   if (verbose >= 4) print*, freq,t,as_ratio,diameter, particle_mass

  freq_round = NINT(freq*1.d-8)*1.d8
  t_round = NINT(t/10.d0)*10.d0
  as_ratio_round = NINT(as_ratio*20.d0)/20.d0

  diameter_magnitude = 10.d0**(INT(log10(diameter)-2))! change 2 to 3 for one digit more accuracy
  particle_mass_magnitude = 10.d0**(INT(log10(particle_mass)-2))

  diameter_round = DBLE(NINT(diameter/diameter_magnitude)&
		  *diameter_magnitude)
  particle_mass_round = DBLE(NINT(particle_mass/particle_mass_magnitude)&
		   *particle_mass_magnitude)

  print*, "D", diameter, diameter_round
  print*, "M", particle_mass, particle_mass_round
  

!   freq_round = freq
!   t_round = t
!   as_ratio_round = as_ratio

  if (verbose >= 4) print*, "freq_round,t_round,as_ratio_round,diameter_round, particle_mass_round, ptype"
  if (verbose >= 4) print*, freq_round,t_round,as_ratio_round,diameter_round, particle_mass_round, ptype


 call sql_get_entry(par_name,freq_round,t_round,as_ratio_round,diameter_round,&
    particle_mass_round,ptype,&
    scatter_matrix,extinct_matrix,emis_vector,found)
! 




  if (.not. found) then

    call tmatrix_refIndex(freq_round,t_round,as_ratio_round,diameter_round,&
      particle_mass_round,ptype,&
      scatter_matrix,extinct_matrix,emis_vector)

    
    call sql_write_entry(par_name,freq_round,t_round,as_ratio_round,diameter_round,&
      particle_mass_round,ptype,&
      scatter_matrix,extinct_matrix,emis_vector)

  end if


  call sql_close_table()

  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return
end subroutine tmatrix_sql


subroutine tmatrix_refIndex(freq,t,as_ratio_in,diameter,particle_mass,ptype,scatter_matrix,extinct_matrix,emis_vector)

  use kinds
  use report_module
  use constants, only: pi, c, im
  use settings, only : nummu, nstokes
  implicit none

  real(kind=dbl), intent(in) :: freq
  real(kind=dbl), intent(in) :: t
  real(kind=dbl), intent(in) :: as_ratio_in
  real(kind=dbl), intent(in) :: particle_mass
  real(kind=dbl), intent(in) :: diameter
  character(len=4), intent(in) :: ptype

  real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector

  real(kind=dbl) :: wavelength
  real(kind=dbl) :: wave_num
  real(kind=dbl) :: equiv_radius
  real(kind=dbl) :: density_eff
  real(kind=dbl) :: refre, refim
  real(kind=dbl) :: eu_alpha, eu_beta
  real(kind=dbl) :: bin_wgt
  real(kind=dbl) :: del_d
  real(kind=dbl) :: as_ratio

  integer :: azimuth_num, azimuth0_num

  complex(kind=ext) :: mindex
  complex(kind=dbl) :: msphere, eps_mix,m_ice

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=14) :: nameOfRoutine = 'tmatrix_refIndex'

  interface

    subroutine tmatrix_calc(quad,qua_num,frequency,wave_num,ref_index,axi, nstokes,&
    as_ratio, alpha, beta, azimuth_num, azimuth0_num,&
    scatter_matrix,extinct_matrix,emis_vector)
      use kinds
      implicit none
      character(1), intent(in) :: quad
      integer, intent(in) :: qua_num
      real(kind=dbl), intent(in) :: frequency
      real(kind=dbl), intent(in) :: wave_num
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
    end subroutine tmatrix_calc
  end interface


  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  if (verbose >= 4) print*, "freq,t,as_ratio_in,diameter,particle_mass"
  if (verbose >= 4) print*, freq,t,as_ratio_in,diameter,particle_mass

  wavelength = c/freq !
  wave_num = 2.0_dbl*pi/wavelength
  eu_alpha = 0.0_dbl    ! orientation of the particle [°]
  eu_beta = 0.0_dbl!orientation of the particle [°]
  azimuth_num = 30
  azimuth0_num = 1


  if (ptype == "snow" .or. ptype == "ice") then
    call ref_ice(t, freq*1d-9, refre, refim)
    m_ice = refre-Im*refim  ! mimicking a

      !XINXINs code
  ! 	    CALL CAL_REFRACTIVE_INDEX('S',t,freq, diameter, as_ratio, particle_mass*1.d3, equiv_radius, mindex)

      !diameter of sphere with same mass
      equiv_radius = 0.5_dbl*diameter*as_ratio_in**(1.0_dbl/3.0_dbl)
      density_eff = (3.d0/4.d0 * particle_mass) / (pi * equiv_radius**3)


      if (density_eff > 917.d0) then
	print*, "WANRING changed density from ", density_eff, "kg/m3 to 917 kg/m3 for d=", diameter
	density_eff = 917.d0
      end if
      if (verbose >= 4) print*, "density_eff, equiv_radius, diameter, particle_mass,as_ratio, "
      if (verbose >= 4) print*, density_eff, equiv_radius, diameter, particle_mass, as_ratio_in

      msphere = eps_mix((1.d0,0.d0),m_ice,density_eff)
      mindex =conjg(msphere) !different convention
    else
      print*, nameOfRoutine, " do not know ptype:", ptype
      stop
    end if

    !make it numerically more stable
    if (as_ratio_in .eq. 1.d0) then
      as_ratio = 0.9999
    else
      as_ratio = as_ratio_in
    end if

    call tmatrix_calc('L',nummu,freq,wave_num,mindex,equiv_radius,nstokes,&
            as_ratio, eu_alpha, eu_beta, azimuth_num, azimuth0_num, &
            scatter_matrix,extinct_matrix,emis_vector)



  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return

end subroutine tmatrix_refIndex





subroutine tmatrix_calc(quad,qua_num,frequency,wave_num,ref_index,axi, nstokes,&
    as_ratio, alpha, beta, azimuth_num, azimuth0_num,&
    scatter_matrix,extinct_matrix,emis_vector)

!DB VERSION

    !  calculate the matrix and vectors, for a single particle with a single orientation
    !
    !   input:
    !       quad            char    name of quadrature
    !       quad_num        int     number of quadrature angles
    !       frequency       double  frequency [Hz]
    !       wave_num        double  wave number [1/m]
    !       ref_index       complex refractive index
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

    use kinds
    use rt_utilities, only: lobatto_quadrature

    implicit none


    character(1), intent(in) :: quad
    integer, intent(in) :: qua_num
    real(kind=dbl), intent(in) :: frequency
    real(kind=dbl), intent(in) :: wave_num
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


    real(kind=dbl) :: maximum_size, minimum_size, &
        num_0, lambda_0, temperature, &
        particle_size, particle_mass, pi
    real(kind=dbl) ::thet0, thet, phi, phi0,&
        rat, sum_tmp, phi_weights, phi0_weights
    integer :: i, j, k, l, m, n, &
       particle_num,&
       ii,jj,kk,ll,kkk,kkk1,kkk2

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


    integer :: nmax, np

    real(kind=ext) :: mrr, mri, lam
    ! some factors that stay constant during calculations
    real(kind=dbl) :: fact_sca
    complex(kind=dbl) :: fact_ext


    np = -1
    rat = 1.e0

    LAM = 3.0e8_ext/(frequency)!*1e6
!     axi = axi!*1.e6
    mrr = REAL(ref_index)
    mri = abs(IMAG(ref_index))


!      RAT=.963711
!      LAM=6.283185   
!      MRR=1.5D0
!      MRI=0.02D0
!      NP=-1

! !ICHOICE=1  NCHECK=1
! RAT= .963711
! !PROLATE SPHEROIDS, A/B=   .5000000
! 
! LAM=  6.283185   
! MRR= .1500D+01   
! MRI= .2000D-01


!ACCURACY OF COMPUTATIONS DDELT =  .10D-02
!EQUAL-SURFACE-AREA-SPHERE RADIUS= 10.0000
!thet0= 56.00  thet= 65.00  phi0=114.00  phi=128.00  alpha=145.00  beta= 52.00
!AMPLITUDE MATRIX
!S11=-.50941D+01 + i* .24402D+02
!S12=-.19425D+01 + i* .19971D+01
!S21=-.11521D+01 + i*-.30977D+01
!S22=-.69323D+01 + i* .24748D+02
!PHASE MATRIX
!  650.3172  -17.9846   10.0498  -12.7580
!  -21.1462  631.6323 -127.3059   87.2144
!    6.8321  132.6132  635.2768  -34.7730
!   -9.6629  -78.1229   51.4094  643.1738
! time =     .35 min
!


! call the tmatrix routine amplq -> fills common block /TMAT/
    call amplq(lam, mrr,mri, AXI, AS_RATIO, RAT, NP,nmax)

    extinct_matrix = 0.e0
    scatter_matrix = 0.e0
    emis_vector = 0.e0
    emis_vector_tmp2_11 = 0.e0
    emis_vector_tmp2_12 = 0.e0
    ! if the particle is rotationally-symmetric, reduce calculation time for orientation-averaging
    ! if not, do orientation averaging for incident and scatterred directions

    !      write(*,*)wave_num
    pi = acos(-1.d0)

    fact_sca = 0.5e0/(wave_num**2)
    fact_ext = 2.0e0*pi*cmplx(0.,1.)/wave_num**2.e0
    ! calculate the quadrature angle, number and weight according to quadrature method
    ! subroutine lobatto_quadrature and gauss_legendre_quadrature in the file named 'refractive_index.f'
    if (quad(1:1).eq.'l'.or.quad(1:1).eq.'L') then
        call lobatto_quadrature(qua_num,qua_angle,qua_weights)
    endif

    ! for each quadrature angle
    ii = 1
    do 1241 jj = 1, qua_num
        thet0=acos(qua_angle(jj)*(-1.)**(real(ii)-1))*180.d0/pi
        thet0_weights = qua_weights(jj)
        if(thet0.gt.179.9999)thet0=180.0d0
        ! initializing the emis vector summation
        emis_vector_tmp1_11 = 0.d0
        emis_vector_tmp1_12 = 0.d0

        do 1242 kk = 1, 2
            kkk1 = kk
            do 1243 ll = 1, qua_num
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



		CALL AMPL(NMAX,dble(LAM),THET0,THET,PHI0,PHI,ALPHA,BETA,&
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



     return
 end subroutine tmatrix_calc

