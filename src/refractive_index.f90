subroutine refractive_index(errorstatus,type, t, f, particle_size, as_ratio, particle_mass, ref)

    use kinds
    use constants, only: pi
    use report_module
    use eps_water, only: get_eps_water
    implicit none

    real(kind=dbl) :: t, f,particle_size, as_ratio, volume
    real(kind=dbl) :: particle_mass, min_dim
    real(kind=dbl) ::  sal, mix_de
    complex(kind=dbl) :: eps_ice, eps_mix
    complex(kind=ext) ::  ref
    complex(kind=dbl) ::  ref_dbl
    ! input:
    !  tempr: temperautre, in k
    !  frq: frequency, in hz,
    !  freq: in refrective index calculation, convert freq to ghz
    !  particle_size: maximum dimension of a single particle, in meter
    !  as_ratio: aspect ratio, minimun_dimension/maximum_dimension
    !  mass: mass of a single particle, in gram
    ! output:
    !  axi: equivalent volume sphere radius, in meter
    !  ref: refractive index of snow/ice
  
    character(1) :: type

    !      if(phase.eq.'S'.or.phase.eq.'s')then
    !	    call refsnow(tempr, frq, particle_size,as_ratio, particle_mass, ref)
    !          print*, tempr, frq, particle_size,as_ratio, particle_mass, axi, ref
    !      elseif(phase.eq.'L'.or. phase.eq.'l'.or.phase .eq. 'r') then
    !        sal = 0
    !        call refwater(sal,frq,tempr,ref)
    !      endif

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'refractive_index'


    err = 0
    if ((type .eq. 'c') .or. (type .eq. 'r')) then
        sal = 0._dbl
         call get_eps_water(err,sal,t-273.15_dbl,f,ref_dbl)
          if (err > 0) then
              errorstatus = fatal
              msg = "error in get_eps_water"
              call report(errorstatus, msg, nameOfRoutine)
              return
          end if   
         ref = sqrt(ref_dbl)
    else if (type .eq. 'i') then
        ref = sqrt(eps_ice(t,f))
    else if ((type .eq. 's') .or. (type .eq. 'g') .or. (type .eq. 'h')) then
        ! the minimum dimension of oblates, particle_size in meter
        min_dim = particle_size*as_ratio
        ! particle volume,
        volume= 4.0_dbl*pi/3.0_dbl*(particle_size/2.0_dbl)**2*(min_dim/2.0_dbl)
        mix_de = particle_mass/volume
!        if (mix_de .lt. 100._dbl) mix_de = 100._dbl
        print*, particle_mass, volume, mix_de
        ref = sqrt(eps_mix((1._dbl,0._dbl),eps_ice(t,f),mix_de))
    else
              errorstatus = fatal
              msg = "No appropriate dielectric model found"
              call report(errorstatus, msg, nameOfRoutine)
              return
    end if

    return
end subroutine refractive_index

! xinxins version
subroutine cal_refractive_index(phase, tempr, frq, particle_size,&
as_ratio, particle_mass, axi, ref)

    use kinds

    implicit none

    real(kind=dbl) :: tempr, frq, particle_size, as_ratio, axi
    real(kind=dbl) :: particle_mass, sal
    complex(kind=ext) :: ref
    character(1) :: phase

    ! input:
    !  tempr: temperautre, in k
    !  frq: frequency, in hz,
    !  freq: in refrective index calculation, convert freq to ghz
    !  particle_size: maximum dimension of a single particle, in meter
    !  as_ratio: aspect ratio, minimun_dimension/maximum_dimension
    !  mass: mass of a single particle, in gram
    ! output:
    !  axi: equivalent volume sphere radius, in meter
    !  ref: refractive index of snow/ice

    if((phase.eq.'S').or.(phase.eq.'s'))then
        call refsnow(tempr, frq, particle_size,as_ratio, particle_mass, axi, ref)
    else if((phase.eq.'L').or. (phase.eq.'l').or.(phase .eq. 'r')) then
        sal = 0._dbl
        call refwater(sal,frq,tempr,ref)
        axi = particle_size/2._dbl
    !   write(*,*)'refractive_index.f','axi = ',axi
    end if

    return
end subroutine cal_refractive_index

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     snow refractive index

subroutine refsnow(tempr, frq, particle_size, as_ratio, particle_mass, axi, snow_ref)

    use kinds
    use constants, only: pi
    implicit none

    real(kind=dbl) :: tempr, frq, freq, particle_size, as_ratio, axi, mass
    real(kind=dbl) :: particle_mass
    ! input:
    !  tempr: temperautre, in k
    !  frq: frequency, in hz,
    !  freq: in refrective index calculation, convert freq to ghz
    !  particle_size: maximum dimension of a single particle, in meter
    !  as_ratio: aspect ratio, minimun_dimension/maximum_dimension
    !  mass: mass of a single particle, in gram
    ! output:
    !  axi: equivalent volume sphere radius, in meter
    !  snow_ref: refractive index of snow

    ! internal variables
    real(kind=dbl) :: min_dim, max_dim, volm
    real(kind=dbl) :: air_dens, ice_dens
    real(kind=dbl) :: realp, alpha1, b1, b, b2, betam, delta_beta, beta1, imagp
    complex(kind=ext) :: ice_ref, snow_ref

    ! min_dim: minimum dimension, in m
    ! max_dim: maximum dimension, in m
    ! volm: particle volume, in cm3
    ! realp,...imagp, temp value for snow permittivity calculation
    ! air_dens: the density of air, g/cm3 (mixing ice with air -- snow)
    ! ice_dens: the density of ice, g/cm3
    !      write(*,*)particle_size, as_ratio
    air_dens = 1.25e-3_dbl
    ice_dens = 0.917_dbl
    freq = frq/1.e9_dbl
    mass = particle_mass
    !      write(*,*)mass
    ! the minimum dimension of oblates, cm. particle_size in meter
    min_dim = particle_size*as_ratio*1.0e2_dbl
    ! the maximum dimension of oblates, cm. particle_size in meter
    max_dim = particle_size*1.0e2_dbl
    ! particle volume,
    volm = 4.0_dbl*pi/3.0_dbl*(max_dim/2.0_dbl)*(max_dim/2.0_dbl)*(min_dim/2.0_dbl)

    ! equivalent-volume sphere's radius, in meter
    axi = 0.5_dbl*(particle_size**3.0_dbl*as_ratio)**(1.0_dbl/3.0_dbl)

    !====================================================================
    !      write(*,*)max_dim, min_dim,as_ratio

    !=========== calculate refractive index of ice
    !=== real part of ice dielectric
    realp = 3.1884_dbl+9.1_dbl*1.e-4_dbl*(tempr-273.0_dbl)
    !=== imaginary part of ice dielectric%%%%%%%%
    alpha1 = (0.00504_dbl+0.0062_dbl*(300._dbl/tempr-1._dbl))*exp(-22.1_dbl*(300._dbl/tempr-1._dbl))
    b1 = 0.0207_dbl
    b = 335.0_dbl
    b2 = 1.16_dbl*1.0e-11_dbl
    betam = b1/tempr*exp(b/tempr)/(exp(b/tempr)-1._dbl)**2+b2*freq*freq
    delta_beta = exp(-9.963_dbl+0.0372_dbl*(tempr-273.16_dbl))
    beta1 = betam+delta_beta
    imagp = alpha1/freq+beta1*freq

    ice_ref = realp*(1._ext,0._ext)+(alpha1/freq+beta1*freq)*(0._ext,1._ext)
    ice_ref = sqrt(ice_ref)

    !======================================================================
    !        relative volume of ice in particle
    snow_ref = (2._ext*(air_dens-mass/volm)/(air_dens-ice_dens)&
    *ice_ref*ice_ref-2*(air_dens-mass/volm)/(air_dens-ice_dens)&
    +ice_ref*ice_ref+2)/&
    (ice_ref*ice_ref+2-(air_dens-mass/volm)/(air_dens-ice_dens)&
    *ice_ref*ice_ref+(air_dens-mass/volm)/(air_dens-ice_dens))
    snow_ref = sqrt(snow_ref)
    !      write(*,*)air_dens, mass, volm, ice_dens, ice_ref
    !      write(*,*)'ice_ref',ice_ref
    !      write(*,*)snow_ref
    !**************************************************************************

    return
end

!     water refractive index
subroutine refwater(sal1,freq,tempr,ref)

    use kinds
    use constants, only: pi

    implicit none

    real(kind=dbl) :: sal1, freq, temp, mrr, mri, sal,tempr
    complex(kind=ext) :: first_term,second_term,third_term,eps,nn,term1,term2
    real(kind=dbl) :: a_1,a_2,a_3,a_4,a_5,a_6,a_7,a_8
    real(kind=dbl) :: a_9,a_10,a_11,a_12,a_13,a_14,a_15,a_16,a_17,a_18
    real(kind=dbl) :: eps_s, eps_1, tau_1, tau_2, eps_inf
    real(kind=dbl) :: c_alpha_0, d_alpha_0, alpha_0, d_p, p, sigma_35
    real(kind=dbl) :: sigma
    complex(kind=ext) :: ref
    real(kind=dbl) :: q,c_p,alpha_1

    sal = sal1 * 1.e-3_dbl
    !      write(*,*)freq,sal
    !*** check the input ranges:
    ! if freq gt 1d12 and verbose eq 1 then print, '!!!frequency range: 0-1000 ghz!!!, extrapolation'
    ! if t gt 303.15 or t lt 273.15 and verbose eq 1 then print, '!!!temperature range: 0-30 degc!!!, extrapolation'
    ! if sal gt 40. or sal lt 0. and verbose eq 1 then print, '!!!salinity range: 0-40 ppt!!!, extrapolation'

    !*** convert temperature from kelvin to degree celsius
    temp = tempr - 273.15_dbl
    !;--------------------------------------------------------------------------------------------------------
    !;coeffs and calculation of eps(freq, temp, sal) according to (5.21, p.445)
    !;--------------------------------------------------------------------------------------------------------

    !;*** coefficients a_i (table 5.5 or p. 454):

    a_1  =  0.46606917e-2_dbl
    a_2  = -0.26087876e-4_dbl
    a_3  = -0.63926782e-5_dbl
    a_4  =  0.63000075e1_dbl
    a_5  =  0.26242021e-2_dbl
    a_6  = -0.42984155e-2_dbl
    a_7  =  0.34414691e-4_dbl
    a_8  =  0.17667420e-3_dbl
    a_9  = -0.20491560e-6_dbl
    a_10 =  0.58366888e3_dbl
    a_11 =  0.12634992e3
    a_12 =  0.69227972e-4_dbl
    a_13 =  0.38957681e-6_dbl
    a_14 =  0.30742330e3_dbl
    a_15 =  0.12634992e3_dbl
    a_16 =  0.37245044e1_dbl
    a_17 =  0.92609781e-2_dbl
    a_18 = -0.26093754e-1_dbl


    !;*** calculate parameter functions (5.24)-(5.28), p.447

    eps_s   = 87.85306_dbl * exp(-0.00456992_dbl * temp - a_1*sal - &
    a_2*sal*sal - a_3*sal*temp)
    eps_1   = a_4 * exp( -a_5*temp - a_6*sal - a_7*sal*temp)
    tau_1   = (a_8 + a_9*sal) * exp( a_10 / (temp + a_11)) * 1e-9_dbl
    tau_2   = (a_12 + a_13*sal) * exp( a_14 / (temp + a_15)) * 1e-9_dbl
    eps_inf = a_16 + a_17*temp + a_18*sal
    !      write(*,*)eps_s,eps_1,tau_1,tau_2,eps_inf

    !;*** calculate seawater conductivity (5.20), p.437

    if (sal.gt.0.) then
        c_alpha_0 =  (6.9431_dbl + 3.2841_dbl * sal - 0.099486_dbl * sal**2.)
        d_alpha_0 =  (84.85_dbl + 69.024_dbl * sal + sal**2.)
        alpha_0   =  c_alpha_0 / d_alpha_0
        alpha_1   = 49.843_dbl - 0.2276_dbl * sal + 0.00198_dbl * sal**2.
        q = 1.000_dbl + alpha_0*(temp - 15.0_dbl) / (temp + alpha_1)
        c_p = (37.5109_dbl + 5.45216_dbl * sal + 0.014409_dbl * sal**2.)
        d_p = (1004.75_dbl + 182.283_dbl * sal + sal**2.)
        p = sal * c_p / d_p
        sigma_35  = 2.903602_dbl + 8.607e-2_dbl * temp+4.738817e-4_dbl*temp**2. &
        -2.991e-6_dbl * temp**3. + 4.3041e-9_dbl * temp**4.
        sigma = sigma_35 * p * q
    else
        sigma = 0._dbl
    endif
    !;just reduce pc time.... calculation would give the same!

    !;*** finally apply the interpolation formula (5.21)
    term1 = 1.*(1.0_ext,0._ext) -2._ext*pi*freq*tau_1*(0._ext,1.0_ext)
    first_term  = (eps_s - eps_1) /  term1
    term2 = 1._ext*(1.0_ext,0._ext) -2._ext*pi*freq*tau_2*(0._ext,1._ext)
    second_term = (eps_1 - eps_inf) / term2
    !;third_term  = dcomplex(eps_inf, (17.9751d * sigma / freq ))
    third_term = eps_inf

    eps = first_term + second_term + third_term
    !      write(*,*)  term1,term2
    !      write(*,*)  first_term,second_term,third_term
    !;calculate refractivity and mass/volume absorption coefficient
    !;frequency in hz, lwc in g/m³
    !      re = (eps-1)/(eps+2)
    !      mass_abscof = 6.d*3.1415926*imag(re)*freq*1d-3/cl
    !      vol_abscof = mass_abscof * lwc
    !;*** convert to refractive index

    nn  = sqrt(eps)
    mrr = real(nn)
    mri = imag(nn)
    ref = nn

    if(mri.lt.0)ref=conjg(nn)

    return
end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5




subroutine vector_sum(dim1, vector, sum1)

    integer dim1
    real*8 sum1, vector(1:dim1)
    integer i
    sum1 = 0    ! initialize the value
    do i = 1, dim1, 1
        sum1 = sum1 + vector(i)
    enddo

    return
end

subroutine gauss_distribution(beta, canting_num, canting_std, &
canting_mean, canting_weights)
    real*8 canting_std, canting_mean, canting_weights
    real*8 sum_p, canting_angle, beta
    integer canting_num

    pi = acos(-1.0)
    sum_p = 0.0

    do 1300 i = 1, canting_num
        canting_angle = 90.0/(real(canting_num)-1)*(real(i)-1)
        sum_p = sum_p + &
        exp(-(canting_angle-canting_mean)**2.0/(2.0*canting_std**2))
1300 continue

     if (canting_std.eq.0) then
         if (beta.eq.0) then
             canting_weights = 1
         else
             canting_weights = 0
         endif
     else
         canting_weights = 1.0/sum_p*exp(-((beta-canting_mean)**2.0)/(2.0*canting_std**2.0))
     !      write(*,*)canting_weights
     endif

     if (canting_num.eq.1) then
         canting_weights = 1d0
     endif



     return
 end
