module eps_water
  use kinds 
  implicit none
  contains

  subroutine get_eps_water(errorstatus, s,T,f,eps_water)

    ! This function calculates the complex permittivity
    ! of natural water including soluted salt.
    !
    ! Input:
    !	s  salinity [ppt]
    !	t  temperature [°C]
    !	f  frequency [GHz]
    !
    ! Result:
    !	 eps_water complex permittivity of natural water
    ! 


    use kinds 
    use settings, only: liq_mod
    use report_module
    implicit none

    real(kind=dbl), intent(in) :: s,& ! salinity [0/00]
        T,& ! temperature [°C]
        f   ! frequency [GHz]

    complex(kind=dbl), intent(out) :: eps_water
    integer(kind=long), intent(out) :: errorstatus

    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'get_eps_water'

    if (verbose >= 4) call report(info,'Start of ', nameOfRoutine)    


    if (liq_mod .eq. 'Ell') then
      if (verbose >= 4) print*, 'Take Ellison model for refractive index'
      call eps_water_ellison(s,T,f,eps_water)

    else if (liq_mod .eq. 'Lie') then
      if (verbose >= 4) print*, 'Take Liebe model for refractive index'
      call eps_water_liebe(s,T,f,eps_water)

    else if (liq_mod .eq. 'Ray') then
      if (verbose >= 4) print*, 'Take Ray model for refractive index'
      call eps_water_ray(s,T,f,eps_water)

    else if (liq_mod .eq. 'Sto') then
      if (verbose >= 4) print*, 'Take Stogryn model for refractive index'
      call eps_water_stogryn(s,T,f,eps_water)

    else if (liq_mod .eq. 'TKC') then
      if (verbose >= 4) print*, 'Take Turner, Kneifel, and Cadeddu model for refractive index'
      call eps_water_tkc(T,f,eps_water)

    else 
      errorstatus = fatal
      msg = "Do not know liq_mod: "//liq_mod
      call report(errorstatus, msg, nameOfRoutine)
      return
    endif

    if (verbose >= 10) print*, "K2, s,T,f", abs((eps_water-1)/(eps_water+2))**2, liq_mod, s,T,f
    if (verbose >= 4) call report(info,'End of ', nameOfRoutine)    
    
    return
  end subroutine get_eps_water

  !! refractive index
  !  ref_wat = sqrt(eps_water)
  !   refre = real(ref_wat)
  !   refim = aimag(ref_wat)

  ! absind  absorptive index 
  ! abscof  absorption coefficient    [1/m]

  !   absind = refim/refre
  !   abscoef = (4*pi*refim*f*1.e9/c)



      subroutine eps_water_ellison(s,T,f,eps)

      ! This function calculates the complex permittivity
      ! of natural water including soluted salt.
      ! Valid parameter range: s 0-40, T 0-30, f 0-500 (1000)
      !
      ! Input:
      ! s  salinity [ppt]
      ! t  temperature [°C]
      ! f  frequency [GHz]
      !
      ! Result:
      !  eps_water complex permittivity of natural water
      ! 
      ! References:
      !      Maetzler 2006: Thermal microwave radiation: Application for remote sensing
      !      section 5.2.5.4

      use kinds 

      implicit none

      real(kind=dbl), intent(in) :: s,& ! salinity [0/00]
            T,& ! temperature [°C]
            f   ! frequency [GHz]
      complex(kind=dbl), intent(out) :: eps

      real(kind=dbl), dimension(18), parameter :: a = (/&
            0.46606917e-02, & ! a1
            -0.26087876e-04, & ! a2
            -0.63926782e-05, & ! a3
            0.63000075e+01, & ! a4
            0.26242021e-02, & ! a5
            -0.42984155e-02, & ! a6
            0.34414691e-04, & ! a7
            0.17667420e-03, & ! a8
            -0.20491560e-06, & ! a9
            0.58366888e+03, & ! a10
            0.12634992e+03, & ! a11
            0.69227972e-04, & ! a12
            0.38957681e-06, & ! a13
            0.30742330e+03, & ! a14
            0.12634992e+03, & ! a15
            0.37245044e+01, & ! a16
            0.92609781e-02, & ! a17
            -0.26093754e-01  & ! a18
            /)

      real(kind=dbl), parameter :: pi = 3.141592653589793

      real(kind=dbl) :: tau_1, tau_2, q, sig

      real(kind=dbl) :: p, alpha_0, alpha_1
      real(kind=dbl) :: sig_s35

      complex(kind=dbl) :: eps_s, eps_1, eps_inf


      p = s*(37.5109+5.45216*s+0.014409*s**2)/(1004.75+182.283*s+s**2)
      alpha_0 = (6.9431+3.2841*s-0.099486*s**2)/(84.85+69.024*s+s**2)
      alpha_1 = 49.843-0.2276*s+0.00198*s**2
      sig_s35 = 2.903602+8.607e-02*T+4.738817e-04*T**2-2.991e-06*T**3+4.3041e-09*T**4

      eps_s = 87.85306*exp(-0.00456992*T-a(1)*s-a(2)*s**2-a(3)*s*T)
      eps_1 = a(4)*exp(-a(5)*T-a(6)*s-a(7)*s*T)
      tau_1 = (a(8)+a(9)*s)*exp(a(10)/(T+a(11)))
      tau_2 = (a(12)+a(13)*s)*exp(a(14)/(T+a(15)))
      eps_inf = a(16)+a(17)*T+a(18)*s
      q = 1+(alpha_0*(T-15))/(T+alpha_1)
      sig = sig_s35*p*q

      eps = (eps_s-eps_1)/(1-(0.,1)*2*pi*f*tau_1)+&
            (eps_1-eps_inf)/(1-(0.,1.)*2*pi*f*tau_2)+eps_inf+(0.,1.)*17.9751*sig/f

      return

      end subroutine eps_water_ellison

      subroutine eps_water_liebe(s,T,f,eps)

      ! This function calculates the complex permittivity
      ! of natural water including soluted salt.
      ! Valid parameter range: f 0-1000
      !
      ! Input:
      ! s  salinity [ppt]
      ! T  temperature [°C]
      ! f  frequency [GHz]
      !
      ! Result:
      !  eps_water complex permittivity of natural water
      ! 
      ! References:
      !      Liebe, Hufford and Manabe, Int. J. IR & MM Waves, Vol. 12 (1991), pp. 659-675
      !      Liebe et al., AGARD Conf. Proc. 542, May 1993

      use kinds 

      implicit none

      real(kind=dbl), intent(in) :: s,& ! salinity [0/00]
            T,& ! temperature [°C]
            f   ! frequency [GHz]
      complex(kind=dbl), intent(out) :: eps

      real(kind=dbl), parameter :: pi = 3.141592653589793

      real(kind=dbl) :: theta_1, fp, fs, T_K

      complex(kind=dbl) :: eps_0, eps_1, eps_2


      !convert temperature to Kelvin
      T_K = T + 273.15
      theta_1 = 1.-300./T_K

      fp = 20.1*exp(7.88*theta_1) !from eq. 2b
      fs = 39.8*fp

      eps_0 = 77.66-103.3*theta_1
      eps_1 = 0.0671*eps_0
      eps_2 = 3.52 !from MPM93

      eps = (eps_0-eps_1)/cmplx(1.,-f/fp)+&
            (eps_1-eps_2)/cmplx(1.,-f/fs) + eps_2

      return

      end subroutine eps_water_liebe

      subroutine eps_water_ray(s,T,f,eps)

      ! This function calculates the complex permittivity
      ! of natural water including soluted salt.
      ! Valid parameter range: 0.2 microns to 10 cm (temperature dependence only considered beyond 0.1 cm)
      !
      ! Input:
      ! s  salinity [ppt]
      ! T  temperature [°C]
      ! f  frequency [GHz]
      !
      ! Result:
      !  eps_water complex permittivity of natural water
      ! 
      ! References:
      ! 
      ! 0.2 um - 0.69 um
      ! 
      ! Hale,G., and M. Querry,1972.
      ! Optical constants of water in the 200 nm to 200 um wavelength regi
      ! Applied Optics,12,3,555-563.
      ! 
      ! 0.69 um - 2.0 um
      ! 
      ! Palmer,K.F., and D. Williams,1974.
      ! Optical properties of water in the near infrared.
      ! Journal of the Optical Society of America,64,8,1107-1110.
      ! 
      ! 2.0 um - 1000.0 um
      ! 
      ! Downing,H.D., and D. Williams,1975.
      ! Optical constants of water in the infrared.
      ! Journal of Geophysical Review,80,12,1656-1661.
      ! 
      ! 1.0 mm - 10.0 cm
      ! 
      ! Ray,P.S.,1972.
      ! Broadband complex refractive indices of ice and water.
      ! Applied Optics,11,8,1836-1844.
      ! 

      use kinds 

      implicit none

      real(kind=dbl), intent(in) :: s,& ! salinity [0/00]
            T,& ! temperature [°C]
            f   ! frequency [GHz]
      complex(kind=dbl), intent(out) :: eps

      integer :: i, i1, i2

      real (kind=dbl) :: wl

      integer, parameter :: numwat = 518
      real (kind=dbl), parameter :: wlmin = 0.2
      real (kind=dbl), parameter :: wlmax = 100000.0
      real (kind=dbl), parameter :: cutwat = 1000.0

      real(kind=dbl), dimension(518), parameter :: wltabw = (/&
               .20000,    .22500,    .25000,    .27500,    .30000,    .32500,&
               .35001,    .37500,    .40000,    .42501,    .45000,    .47499,&
               .50000,    .52499,    .54999,    .57501,    .59999,    .62500,&
               .64998,    .67499,    .68966,    .70175,    .71429,    .72464,&
               .73529,    .74627,    .75188,    .75758,    .76923,    .78125,&
               .79365,    .80645,    .81301,    .81967,    .83333,    .84746,&
               .86207,    .87719,    .89286,    .90909,    .92593,    .93458,&
               .94340,    .95238,    .96154,    .97276,    .98039,    .99010,&
               1.00000,   1.01010,   1.02041,   1.03093,   1.04167,   1.05263,&
               1.06952,   1.08696,   1.09890,   1.11111,   1.12360,   1.13636,&
               1.14943,   1.16279,   1.17647,   1.19048,   1.20482,   1.21951,&
               1.23457,   1.25000,   1.26582,   1.28205,   1.29870,   1.31579,&
               1.33333,   1.35135,   1.36986,   1.38889,   1.40845,   1.42857,&
               1.44300,   1.47059,   1.49254,   1.51515,   1.53846,   1.56250,&
               1.58730,   1.61290,   1.63934,   1.66667,   1.69492,   1.72414,&
               1.75439,   1.78571,   1.80180,   1.81818,   1.85185,   1.88679,&
               1.92678,   1.96078,   2.00000,   2.02020,   2.04082,   2.06186,&
               2.08333,   2.10526,   2.12766,   2.15054,   2.17391,   2.19780,&
               2.22222,   2.24719,   2.27273,   2.29885,   2.32558,   2.35294,&
               2.38095,   2.40964,   2.43902,   2.46914,   2.50000,   2.50627,&
               2.51256,   2.51889,   2.52525,   2.53165,   2.53807,   2.54453,&
               2.55102,   2.55754,   2.56410,   2.57069,   2.57732,   2.58398,&
               2.59067,   2.59740,   2.60417,   2.61097,   2.61780,   2.62467,&
               2.63158,   2.63852,   2.64550,   2.65252,   2.65957,   2.66667,&
               2.67380,   2.68097,   2.68817,   2.69542,   2.70270,   2.71003,&
               2.71739,   2.72480,   2.73224,   2.73973,   2.74725,   2.75482,&
               2.76243,   2.77008,   2.77778,   2.78552,   2.79330,   2.80112,&
               2.80899,   2.81690,   2.82486,   2.83286,   2.84091,   2.84900,&
               2.85714,   2.86533,   2.87356,   2.88184,   2.89017,   2.89855,&
               2.90698,   2.91545,   2.92398,   2.93255,   2.94118,   2.94985,&
               2.95858,   2.96736,   2.97619,   2.98507,   2.99401,   3.00300,&
               3.01205,   3.02115,   3.03030,   3.03951,   3.04878,   3.05810,&
               3.06748,   3.07692,   3.08642,   3.09598,   3.10559,   3.11526,&
            3.12500,   3.13480,   3.14465,   3.15457,   3.16456,   3.17460,&
            3.18471,   3.19489,   3.20513,   3.21543,   3.22581,   3.23625,&
            3.24675,   3.25733,   3.26797,   3.27869,   3.28947,   3.30033,&
            3.31126,   3.32226,   3.33333,   3.34448,   3.35570,   3.36700,&
            3.37838,   3.38983,   3.40136,   3.41297,   3.42466,   3.43643,&
            3.44828,   3.46021,   3.47222,   3.48432,   3.49650,   3.50877,&
            3.52113,   3.53357,   3.54610,   3.55872,   3.57143,   3.58423,&
            3.59712,   3.61011,   3.62319,   3.63636,   3.64964,   3.66300,&
            3.67647,   3.69004,   3.70370,   3.71747,   3.73134,   3.74532,&
            3.75940,   3.77358,   3.78788,   3.80228,   3.81679,   3.83142,&
            3.84615,   3.86100,   3.87597,   3.89105,   3.90625,   3.92157,&
            3.93701,   3.95257,   3.96825,   3.98406,   4.00000,   4.01606,&
            4.03226,   4.04858,   4.06504,   4.08163,   4.09836,   4.11523,&
            4.13223,   4.14938,   4.16667,   4.18410,   4.20168,   4.21941,&
            4.23729,   4.25532,   4.27350,   4.29185,   4.31034,   4.32900,&
            4.34783,   4.36681,   4.38596,   4.40529,   4.42478,   4.44444,&
            4.46429,   4.48430,   4.50450,   4.52489,   4.54545,   4.56621,&
            4.58716,   4.60829,   4.62963,   4.65116,   4.67290,   4.69484,&
            4.71698,   4.73934,   4.76190,   4.78469,   4.80769,   4.83092,&
            4.85437,   4.87805,   4.90196,   4.92611,   4.95050,   4.97512,&
            5.00000,   5.02513,   5.05051,   5.07614,   5.10204,   5.12821,&
            5.15464,   5.18135,   5.20833,   5.23560,   5.26316,   5.29101,&
            5.31915,   5.34759,   5.37634,   5.40541,   5.43478,   5.46448,&
            5.49451,   5.52486,   5.55556,   5.58659,   5.61798,   5.64972,&
            5.68182,   5.71429,   5.74713,   5.78035,   5.81395,   5.84795,&
            5.88235,   5.91716,   5.95238,   5.98802,   6.02410,   6.06061,&
            6.09756,   6.13497,   6.17284,   6.21118,   6.25000,   6.28931,&
            6.32911,   6.36943,   6.41026,   6.45161,   6.49351,   6.53595,&
            6.57895,   6.62252,   6.66667,   6.71141,   6.75676,   6.80272,&
            6.84932,   6.89655,   6.94444,   6.99301,   7.04225,   7.09220,&
            7.14286,   7.19424,   7.24638,   7.29927,   7.35294,   7.40741,&
            7.46269,   7.51880,   7.57576,   7.63359,   7.69231,   7.75194,&
            7.81250,   7.87402,   7.93651,   8.00000,   8.06452,   8.13008,&
            8.19672,   8.26446,   8.33333,   8.40336,   8.47458,   8.54701,&
            8.62069,   8.69565,   8.77193,   8.84956,   8.92857,   9.00901,&
            9.09091,   9.17431,   9.25926,   9.34579,   9.43396,   9.52381,&
            9.61538,   9.70874,   9.80392,   9.90099,  10.00000,  10.10101,&
            10.20408,  10.30928,  10.41667,  10.52632,  10.63830,  10.75269,&
            10.86957,  10.98901,  11.11111,  11.23596,  11.36364,  11.49425,&
            11.62791,  11.76471,  11.90476,  12.04819,  12.19512,  12.34568,&
            12.50000,  12.65823,  12.82051,  12.98701,  13.15789,  13.33333,&
            13.51351,  13.69863,  13.88889,  14.08451,  14.28571,  14.49275,&
            14.70588,  14.92537,  15.15152,  15.38462,  15.62500,  15.87302,&
            16.12903,  16.39344,  16.66667,  16.94915,  17.24138,  17.54386,&
            17.85714,  18.18182,  18.51852,  18.86792,  19.23077,  19.60784,&
            20.00000,  20.40816,  20.83333,  21.27660,  21.73913,  22.22222,&
            22.72727,  23.25581,  23.80952,  24.39024,  25.00000,  25.64103,&
            26.31579,  27.02703,  27.77778,  28.57143,  29.41176,  30.30303,&
            31.25000,  32.25806,  33.33333,  34.48276,  35.71429,  37.03704,&
            38.46154,  40.00000,  41.66667,  43.47826,  45.45455,  47.61905,&
            50.00000,  52.63158,  55.55556,  58.82353,  62.50000,  66.66667,&
            71.42857,  76.92308,  83.33333,  90.90909, 100.00000, 111.11111,&
            125.00000, 142.85714, 166.66667, 200.00000, 250.00000, 333.33333,&
            500.00000,1000.00000/)

      real(kind=dbl), dimension(518), parameter :: rntabw = (/&
         1.396,1.373,1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,&
         1.336,1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.332,1.332,&
         1.332,1.332,1.332,1.332,1.332,1.332,1.331,1.331,1.331,1.331,1.331,&
         1.330,1.330,1.330,1.330,1.330,1.329,1.329,1.329,1.329,1.329,1.328,&
         1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,1.328,&
         1.327,1.327,1.327,1.327,1.326,1.326,1.326,1.326,1.325,1.325,1.325,&
         1.325,1.325,1.324,1.324,1.324,1.324,1.323,1.323,1.323,1.322,1.322,&
         1.321,1.321,1.321,1.320,1.320,1.319,1.319,1.318,1.318,1.317,1.316,&
         1.315,1.314,1.314,1.313,1.312,1.312,1.311,1.310,1.309,1.307,1.306,&
         1.301,1.301,1.300,1.298,1.298,1.296,1.295,1.294,1.293,1.291,1.289,&
         1.287,1.285,1.282,1.280,1.277,1.274,1.270,1.265,1.261,1.260,1.259,&
         1.257,1.256,1.255,1.254,1.252,1.250,1.249,1.247,1.246,1.243,1.241,&
         1.240,1.238,1.235,1.232,1.230,1.227,1.224,1.221,1.218,1.214,1.210,&
         1.205,1.200,1.195,1.191,1.185,1.179,1.172,1.166,1.157,1.149,1.144,&
         1.139,1.138,1.138,1.139,1.141,1.144,1.149,1.154,1.158,1.161,1.165,&
         1.171,1.177,1.183,1.191,1.199,1.212,1.220,1.233,1.246,1.258,1.271,&
         1.282,1.293,1.305,1.317,1.329,1.342,1.353,1.364,1.376,1.386,1.398,&
         1.407,1.417,1.426,1.434,1.442,1.450,1.457,1.465,1.471,1.476,1.480,&
         1.483,1.486,1.487,1.487,1.487,1.486,1.485,1.482,1.479,1.477,1.474,&
         1.472,1.467,1.464,1.461,1.457,1.454,1.451,1.448,1.444,1.441,1.437,&
         1.434,1.431,1.427,1.425,1.421,1.418,1.415,1.413,1.410,1.407,1.405,&
         1.403,1.400,1.398,1.396,1.394,1.392,1.390,1.388,1.387,1.385,1.383,&
         1.382,1.379,1.378,1.377,1.375,1.374,1.372,1.371,1.370,1.369,1.367,&
         1.366,1.365,1.363,1.361,1.361,1.360,1.358,1.358,1.357,1.355,1.354,&
         1.353,1.352,1.351,1.350,1.349,1.348,1.348,1.347,1.346,1.345,1.344,&
         1.344,1.343,1.342,1.341,1.340,1.340,1.338,1.337,1.337,1.335,1.334,&
         1.334,1.333,1.332,1.332,1.331,1.330,1.330,1.330,1.329,1.329,1.329,&
         1.328,1.328,1.327,1.327,1.327,1.327,1.327,1.326,1.326,1.326,1.325,&
         1.325,1.325,1.325,1.325,1.325,1.324,1.324,1.323,1.322,1.322,1.321,&
         1.320,1.319,1.318,1.318,1.317,1.316,1.314,1.313,1.311,1.310,1.308,&
         1.306,1.304,1.302,1.299,1.297,1.294,1.291,1.288,1.285,1.282,1.278,&
         1.275,1.271,1.267,1.262,1.256,1.251,1.247,1.242,1.241,1.241,1.247,&
         1.265,1.289,1.311,1.332,1.349,1.354,1.356,1.354,1.350,1.345,1.341,&
         1.337,1.333,1.330,1.326,1.324,1.322,1.320,1.319,1.318,1.317,1.316,&
         1.315,1.314,1.313,1.311,1.310,1.309,1.308,1.307,1.306,1.305,1.303,&
         1.302,1.301,1.300,1.298,1.296,1.295,1.294,1.293,1.291,1.288,1.286,&
         1.285,1.283,1.281,1.279,1.276,1.274,1.271,1.269,1.267,1.264,1.261,&
         1.259,1.256,1.253,1.249,1.246,1.242,1.238,1.234,1.230,1.224,1.220,&
         1.214,1.208,1.202,1.194,1.189,1.181,1.174,1.168,1.162,1.156,1.149,&
         1.143,1.139,1.135,1.132,1.132,1.131,1.132,1.130,1.130,1.134,1.138,&
         1.142,1.157,1.171,1.182,1.189,1.201,1.213,1.223,1.236,1.249,1.264,&
         1.277,1.289,1.303,1.313,1.324,1.335,1.348,1.361,1.372,1.385,1.396,&
         1.407,1.419,1.431,1.441,1.451,1.462,1.470,1.480,1.488,1.496,1.504,&
         1.510,1.515,1.521,1.527,1.532,1.537,1.541,1.545,1.549,1.552,1.552,&
         1.552,1.550,1.546,1.543,1.541,1.539,1.537,1.534,1.532,1.529,1.525,&
         1.528,1.542,1.567,1.600,1.640,1.689,1.746,1.801,1.848,1.890,1.929,&
         1.960,1.982,1.997,2.000,2.010,2.020,2.040,2.070,2.110,2.150,2.225,&
         2.481/)

      real(kind=dbl), dimension(518), parameter :: cntabw = (/&
         1.1000e-07,4.9000e-08,3.4000e-08,2.4000e-08,1.6000e-08,1.1000e-08,&
         6.5000e-09,3.5000e-09,1.9000e-09,1.3000e-09,1.0000e-09,9.4000e-10,&
         1.0000e-09,1.3000e-09,2.0000e-09,3.6000e-09,1.1000e-08,1.4000e-08,&
         1.6000e-08,2.2000e-08,2.7000e-08,3.8000e-08,5.6000e-08,7.7300e-08,&
         1.3900e-07,1.6300e-07,1.6800e-07,1.6400e-07,1.5400e-07,1.4300e-07,&
         1.3300e-07,1.2500e-07,1.2400e-07,1.3000e-07,2.0400e-07,2.6100e-07,&
         2.9400e-07,3.5300e-07,4.3300e-07,5.4300e-07,8.7700e-07,1.1800e-06,&
         1.6100e-06,2.4400e-06,3.6000e-06,3.9800e-06,3.9200e-06,3.7000e-06,&
         3.3100e-06,2.8200e-06,2.3100e-06,1.9000e-06,1.5700e-06,1.3700e-06,&
         1.2600e-06,1.4400e-06,1.6800e-06,2.0500e-06,2.8900e-06,4.9600e-06,&
         8.8700e-06,1.0900e-05,1.1500e-05,1.1800e-05,1.2000e-05,1.1800e-05,&
         1.1500e-05,1.1000e-05,1.0800e-05,1.1500e-05,1.3800e-05,1.7500e-05,&
         2.3900e-05,4.1600e-05,5.9400e-05,1.0100e-04,2.4100e-04,3.5200e-04,&
         3.6400e-04,3.3400e-04,2.5800e-04,1.8800e-04,1.4800e-04,1.2000e-04,&
         1.0200e-04,8.7300e-05,7.9200e-05,7.4900e-05,7.6200e-05,8.5500e-05,&
         1.0600e-04,1.3000e-04,1.3600e-04,1.3700e-04,1.5900e-04,8.6300e-04,&
         1.9000e-03,1.7000e-03,1.1000e-03,9.0000e-04,7.3100e-04,6.1700e-04,&
         5.1400e-04,4.5200e-04,4.0000e-04,3.5900e-04,3.4100e-04,3.3800e-04,&
         3.4500e-04,3.7600e-04,4.1600e-04,4.6500e-04,5.4200e-04,6.5200e-04,&
         7.9200e-04,9.6800e-04,1.2300e-03,1.5600e-03,1.9000e-03,1.9500e-03,&
         2.0000e-03,2.0500e-03,2.0700e-03,2.1000e-03,2.1200e-03,2.1500e-03,&
         2.1900e-03,2.2400e-03,2.2700e-03,2.3100e-03,2.3400e-03,2.3900e-03,&
         2.4300e-03,2.4800e-03,2.5700e-03,2.7000e-03,2.9800e-03,3.3000e-03,&
         4.0200e-03,4.3700e-03,4.8200e-03,5.3600e-03,6.2700e-03,7.3200e-03,&
         8.5500e-03,1.0500e-02,1.2700e-02,1.4500e-02,1.6400e-02,1.8600e-02,&
         2.0500e-02,2.8200e-02,3.8000e-02,4.6200e-02,5.4800e-02,6.4900e-02,&
         7.4400e-02,8.3600e-02,9.2700e-02,1.0200e-01,1.1200e-01,1.2100e-01,&
         1.3100e-01,1.4200e-01,1.5400e-01,1.6700e-01,1.8000e-01,1.9400e-01,&
         2.0600e-01,2.1800e-01,2.2900e-01,2.3900e-01,2.4900e-01,2.5800e-01,&
         2.6500e-01,2.7100e-01,2.7600e-01,2.8000e-01,2.8100e-01,2.8200e-01,&
         2.8200e-01,2.7900e-01,2.7600e-01,2.7200e-01,2.6700e-01,2.6200e-01,&
         2.5500e-01,2.5000e-01,2.4300e-01,2.3600e-01,2.2800e-01,2.2000e-01,&
         2.1200e-01,2.0400e-01,1.9500e-01,1.8300e-01,1.7300e-01,1.6300e-01,&
         1.5300e-01,1.4400e-01,1.3400e-01,1.2500e-01,1.1700e-01,1.1000e-01,&
         9.9400e-02,9.2000e-02,8.5500e-02,7.8500e-02,7.1600e-02,6.5300e-02,&
         6.0000e-02,5.5000e-02,5.0400e-02,4.6200e-02,4.2200e-02,3.8500e-02,&
         3.4800e-02,3.1500e-02,2.9700e-02,2.7900e-02,2.6200e-02,2.5000e-02,&
         2.2900e-02,2.1000e-02,1.9300e-02,1.7700e-02,1.6300e-02,1.5100e-02,&
         1.3800e-02,1.2800e-02,1.1800e-02,1.1000e-02,1.0100e-02,9.4100e-03,&
         8.6600e-03,8.0700e-03,7.3700e-03,6.8300e-03,6.2500e-03,5.7900e-03,&
         5.3800e-03,5.0600e-03,4.7300e-03,4.4900e-03,4.2400e-03,4.0500e-03,&
         3.8900e-03,3.7600e-03,3.6300e-03,3.5500e-03,3.4700e-03,3.4000e-03,&
         3.3500e-03,3.3600e-03,3.3500e-03,3.3900e-03,3.4000e-03,3.4800e-03,&
         3.5200e-03,3.6300e-03,3.7000e-03,3.7800e-03,3.8900e-03,3.9900e-03,&
         4.1000e-03,4.2200e-03,4.3300e-03,4.5000e-03,4.6500e-03,4.7900e-03,&
         4.9400e-03,5.1200e-03,5.3100e-03,5.4900e-03,5.6800e-03,5.8600e-03,&
         6.0800e-03,6.3100e-03,6.5300e-03,6.7300e-03,6.9600e-03,7.2200e-03,&
         7.4900e-03,7.7900e-03,8.0600e-03,8.3300e-03,8.6400e-03,8.9600e-03,&
         9.2700e-03,9.6600e-03,1.0000e-02,1.0400e-02,1.0800e-02,1.1200e-02,&
         1.1700e-02,1.2200e-02,1.2600e-02,1.3100e-02,1.3600e-02,1.4000e-02,&
         1.4500e-02,1.4900e-02,1.5200e-02,1.5400e-02,1.5600e-02,1.5700e-02,&
         1.5700e-02,1.5700e-02,1.5500e-02,1.5300e-02,1.5100e-02,1.4800e-02,&
         1.4600e-02,1.4300e-02,1.4000e-02,1.3700e-02,1.3300e-02,1.2900e-02,&
         1.2600e-02,1.2200e-02,1.1800e-02,1.1500e-02,1.1000e-02,1.0800e-02,&
         1.0500e-02,1.0300e-02,1.0100e-02,1.0000e-02,9.9300e-03,9.9000e-03,&
         9.9500e-03,1.0000e-02,1.0200e-02,1.0400e-02,1.0700e-02,1.1000e-02,&
         1.1500e-02,1.2000e-02,1.2800e-02,1.3800e-02,1.5000e-02,1.6600e-02,&
         1.8500e-02,2.0500e-02,2.4200e-02,2.9300e-02,3.3200e-02,4.2900e-02,&
         5.4400e-02,6.8800e-02,8.4000e-02,1.0210e-01,1.1700e-01,1.3000e-01,&
         1.3200e-01,1.2400e-01,1.0600e-01,8.8000e-02,7.4000e-02,6.1800e-02,&
         5.3500e-02,4.8400e-02,4.4700e-02,4.2000e-02,3.9800e-02,3.8300e-02,&
         3.7300e-02,3.7000e-02,3.6600e-02,3.6300e-02,3.6000e-02,3.5700e-02,&
         3.5500e-02,3.5200e-02,3.5000e-02,3.4700e-02,3.4600e-02,3.4300e-02,&
         3.4200e-02,3.4200e-02,3.4200e-02,3.4300e-02,3.4200e-02,3.4200e-02,&
         3.4200e-02,3.4200e-02,3.4200e-02,3.4400e-02,3.4500e-02,3.4600e-02,&
         3.4900e-02,3.5100e-02,3.5100e-02,3.5100e-02,3.5200e-02,3.5600e-02,&
         3.5900e-02,3.6100e-02,3.6200e-02,3.6600e-02,3.7000e-02,3.7400e-02,&
         3.7800e-02,3.8300e-02,3.8700e-02,3.9200e-02,3.9800e-02,4.0500e-02,&
         4.1100e-02,4.1700e-02,4.2400e-02,4.3400e-02,4.4300e-02,4.5300e-02,&
         4.6700e-02,4.8100e-02,4.9700e-02,5.1500e-02,5.3400e-02,5.5700e-02,&
         5.8900e-02,6.2200e-02,6.6100e-02,7.0700e-02,7.6400e-02,8.2800e-02,&
         8.9800e-02,9.7300e-02,1.0700e-01,1.1800e-01,1.3000e-01,1.4400e-01,&
         1.5900e-01,1.7600e-01,1.9200e-01,2.0800e-01,2.2600e-01,2.4300e-01,&
         2.6000e-01,2.7700e-01,2.9200e-01,3.0500e-01,3.1700e-01,3.2800e-01,&
         3.3800e-01,3.4700e-01,3.5600e-01,3.6500e-01,3.7300e-01,3.7900e-01,&
         3.8600e-01,3.9200e-01,3.9700e-01,4.0300e-01,4.0800e-01,4.1200e-01,&
         4.1700e-01,4.2000e-01,4.2300e-01,4.2500e-01,4.2700e-01,4.2800e-01,&
         4.2700e-01,4.2700e-01,4.2600e-01,4.2500e-01,4.2300e-01,4.2100e-01,&
         4.1800e-01,4.1500e-01,4.1100e-01,4.0800e-01,4.0400e-01,4.0100e-01,&
         3.9700e-01,3.9400e-01,3.9000e-01,3.8600e-01,3.8200e-01,3.7700e-01,&
         3.7200e-01,3.6800e-01,3.6300e-01,3.5900e-01,3.5600e-01,3.5200e-01,&
         3.5300e-01,3.5700e-01,3.6100e-01,3.6800e-01,3.7500e-01,3.8500e-01,&
         3.9800e-01,4.1400e-01,4.3600e-01,4.6900e-01,5.0500e-01,5.3900e-01,&
         5.7100e-01,5.9700e-01,6.1800e-01,6.2900e-01,6.2200e-01,6.0800e-01,&
         5.9300e-01,5.7700e-01,5.5700e-01,5.3200e-01,5.0700e-01,4.8700e-01,&
         4.6600e-01,4.5000e-01,4.4400e-01,4.3800e-01,4.6000e-01,5.2700e-01,&
         7.1800e-01,8.4657e-01/)

      real(kind=dbl), parameter :: pi = 3.141592653589793
      
      real (kind=dbl) :: fac, rn, cn
      real (kind=dbl) :: T1, T2, xl

      real(kind=dbl) :: sigma, alpha

      real(kind=dbl) :: eps_s, eps_inf, xlams
      real(kind=dbl) :: powtrm, denom
      real(kind=dbl) :: eps_r, eps_i

      complex(kind=dbl) ::  n

      !convert frequency to microns
      wl = 1.e6*2.99792458e8/(f*1.e9)
      if ((wl .lt. wlmin) .or. (wl .gt. wlmax)) then
        print*,wl, wlmin, wlmax, "out of wavelength specifications"
        stop 
        end if 
      !0.2 micron to 1000.0 microns - table lookup
      if (wl .gt. cutwat) goto 356
      do i=2, numwat
         if (wl .gt. wltabw(i)) goto 343
         i1 = i-1
         i2 = i
         goto 348
      343 CONTINUE
      end do 

      i1 = numwat - 2 
      i2 = numwat - 1
      348 CONTINUE
      fac = (wl-wltabw(i1))/(wltabw(i2)-wltabw(i1))
      rn = rntabw(i1)+fac*(rntabw(i2)-rntabw(i1))
      cn = cntabw(i1)+fac*(cntabw(i2)-cntabw(i1))
      goto 402
      
      !0.1 cm to 10 cm
      !define temperature terms and wavelength in cm
      356 CONTINUE
      T1=T+273.0
      T2=T-25.0
      xl=wl/10000.0

      !define frequency independent conductivity (sigma) and spread parameter (alpha)
      !sigma given by Saxton,J.A.,1949,Wireless Engineer,26,p.288
      !alpha given by Ray (eq. 7B)
      sigma=12.5664e8
      alpha=-16.8129/T1+0.0609265

      !define static dielectric constant (eps_s) - Ray eq. 4
      !define high frequency dielectric constant (eps_inf) - Ray eq. 7A
      !define relaxation wavelength in cm (xlams) - Ray eq. 7C
      
      !temperature dependence of eps_s given by
      !Wyman,J., and E.N.Ingalls,1938,Jour.Am.Chem.Soc.,60,p.1182
      eps_s = 78.54*(1.0-4.579e-3*T2+1.19e-5*T2**2-2.8e-8*T2**3)
      eps_inf = 5.27137+0.0216474*T-0.00131198*T**2
      xlams = 0.00033836*exp(2513.98/T1)

      !calculate expressions used for dielectric constant
      powtrm=(xlams/xl)**(1.-alpha)
      denom=1.+2*powtrm*sin(pi*alpha/2)+(xlams/xl)**(2.*(1.-alpha))

      !calculation of dielectric constant
      !real part - Ray eq. 5
      eps_r = eps_inf+(eps_s-eps_inf)*(1.0+powtrm*sin(pi*alpha/2))/denom

      !imaginary part - Ray eq. 6
      eps_i = (sigma*xl/18.8496e10)+(eps_s-eps_inf)*powtrm*cos(pi*alpha/2)/denom

      !complex permittivity
      eps = dcmplx(eps_r,eps_i)
!       print*,eps_water_ray
      n = sqrt(eps)
      rn = real(n)
      cn = aimag(n)

      !correction to imaginary index to account for the
      !remaining absorption bands - Ray eq. 8 (table 2)
      if (wl .gt. 3000.0) goto 402
      cn = cn + 0.39*exp(-abs(log10(wl/17.0)/0.45)**1.3)&
               + 0.41*exp(-abs(log10(wl/62.0)/0.35)**1.7)&
               + 0.25*exp(-abs(log10(wl/300.0)/0.47)**3.0)

      402 CONTINUE
      n = dcmplx(rn,cn)
      eps = n**2
      return

      end subroutine eps_water_ray

      subroutine eps_water_stogryn(s,T,f,eps)

      ! This function calculates the complex permittivity
      ! of natural water including soluted salt.
      !
      ! Input:
      ! s  salinity [ppt]
      ! t  temperature [°C]
      ! f  frequency [GHz]
      !
      ! Result:
      !  eps_water complex permittivity of natural water
      ! 
      ! References:
      !      Maetzler 2006: Thermal microwave radiation: Application for remote sensing
      !      section 5.2.4.2 (after Stogryn, A. P., Bull, H. T., Rubayi, K., and Iravanchy, S.: 
      !      The microwave permittivity of sea and fresh water. Aerojet internal report, 1996)

      use kinds 

      implicit none

      real(kind=dbl), intent(in) :: s,& ! salinity [0/00]
            T,& ! temperature [°C]
            f   ! frequency [GHz]
      complex(kind=dbl), intent(out) :: eps

      real(kind=dbl) :: a, b

      real(kind=dbl) :: two_pi_tau_1, two_pi_tau_2, q, sig

      real(kind=dbl) :: p, alpha_0, alpha_1
      real(kind=dbl) :: sig_s35

      complex(kind=dbl) :: eps_s, eps_1, eps_inf

      p = s*(37.5109+5.45216*s+0.014409*s**2)/(1004.75+182.283*s+s**2)
      alpha_0 = (6.9431+3.2841*s-0.099486*s**2)/(84.85+69.024*s+s**2)
      alpha_1 = 49.843-0.2276*s+0.00198*s**2
      sig_s35 = 2.903602+8.607e-02*T+4.738817e-04*T**2-2.991e-06*T**3+4.3041e-09*T**4

      a = 1.-(s*(0.03838+0.00218*s)*(79.88+T))/((12.01+s)*(52.53+T)) 
      b = 1.-s*((0.03409+0.002817*s)/(7.69+s))-s*T*((0.00246+0.00141*T)/(188.-7.57*T+T**2))

      eps_s = ((37088.6-82.168*T)/(421.854+T))*a
      eps_1 = (0.0787*((37088.6-82.168*T)/(421.854+T)))*a
      two_pi_tau_1 = ((255.04+0.7246*T)/((49.25+T)*(45.+T)))*b
      two_pi_tau_2 = 0.00628
      eps_inf = 4.05+0.0186*T
      q = 1+(alpha_0*(T-15))/(T+alpha_1)
      sig = sig_s35*p*q

      eps = (eps_s-eps_1)/(1-(0.,1.)*f*two_pi_tau_1)+&
            (eps_1-eps_inf)/(1-(0.,1.)*f*two_pi_tau_2)+eps_inf+(0.,1.)*17.9751*sig/f

      return

      end subroutine eps_water_stogryn

      subroutine eps_water_tkc(T,f,eps)   

      ! Abstract:
      !  The "Turner-Kneifel-Cadeddu" liquid water absorption model (submitted to JTECH 2015).
      !  It was built using both laboratory observations (primarily at warm temperatures) and 
      !  field data observed by MWRs at multiple frequencies at supercool temperatures. The field
      !  data were published in Kneifel et al. JAMC 2014.  The strength of the TKC model is the 
      !  use of an optimal estimation framework to determine the empirical coefficients of the 
      !  double-Debye model.  A full description of this model is given in 
      !
      !  Turner, D.D., S. Kneifel, and M.P. Cadeddu, 2015: An improved liquid
      !  water absorption model in the microwave for supercooled liquid clouds.
      !  J. Atmos. Oceanic Technol., submitted April 2015.
      !
      ! Note that the model is designed to operate over the frequency range 
      ! from 0.5 to 500 GHz, and temperatures from -40 degC to +50 degC.  
      !
      ! Authors:
      !  Dave Turner, National Severe Storms Laborotory / NOAA
      !  Stefan Kneifel, McGill University and the University of Cologne
      !  Maria Cadeddu, Argonne National Laboratory
      !
      ! Input:
      ! T  temperature [°C]
      ! f  frequency [GHz]
      !
      ! Result:
      !  eps_water_tkc complex permittivity of natural water
      !       
      use kinds
      use constants, only: c, pi

      implicit none
      
      real(kind=dbl), intent(in) :: &
	    T,& ! temperature [°C]
      f   ! frequency [GHz]
      complex(kind=dbl), intent(out) :: eps

      real(kind=dbl) :: frq, a_1, b_1, c_1, d_1, a_2, b_2, c_2, d_2, t_c, &
	      eps_s, delta_1, tau_1, delta_2, tau_2, term1_p1, term2_p1, &
	      eps1, eps2
      

      ! Empirical coefficients for the TKC model. The first 4 are a1, b1, c1, and d1, 
      ! the next four are a2, b2, c2, and d2, and the last one is tc.
      real(kind=dbl), dimension(9), parameter :: coef = (/&
	  8.169396d+01, 4.410555d-03, 1.208992d-13, 6.768869d+02, &
	  1.597733d+00, 1.060228d-02, 9.982113d-15, 5.720517d+02, &
	  1.351758d+02/)

      ! Convert the frequency from GHz to Hz

      frq = f * 1d9


      ! This helps to understand how things work below
      a_1 = coef(1)
      b_1 = coef(2)
      c_1 = coef(3)
      d_1 = coef(4)

      a_2 = coef(5)
      b_2 = coef(6)
      c_2 = coef(7)
      d_2 = coef(8)

      t_c = coef(9)


      ! Compute the static dielectric permittivity (Eq 6)
      eps_s = 87.9144d0 - 0.404399d0 * T + 9.58726d-4 * T**2. - 1.32802d-6 * T**3.

      ! Compute the components of the relaxation terms (Eqs 9 and 10)
      ! First Debye component
      delta_1 = a_1 * exp(-b_1 * T)
      tau_1   = c_1 * exp(d_1 / (T + t_c))
      ! Second Debye component
      delta_2 = a_2 * exp(-b_2 * T)
      tau_2   = c_2 * exp(d_2 / (T + t_c))

      ! Compute the relaxation terms (Eq 7) for the two Debye components
      term1_p1 = (tau_1**2.*delta_1) / (1.d0 + (2.d0*pi*frq*tau_1)**2.)
      term2_p1 = (tau_2**2.*delta_2) / (1.d0 + (2.d0*pi*frq*tau_2)**2.)

      ! Compute the real permittivitity coefficient (Eq 4)
      eps1 = eps_s - ((2.d0*pi*frq)**2.)*(term1_p1 + term2_p1) 
   

      ! Compute the relaxation terms (Eq 8) for the two Debye components
      term1_p1 = (tau_1 * delta_1) / (1.d0 + (2.d0*pi*frq*tau_1)**2.)
      term2_p1 = (tau_2 * delta_2) / (1.d0 + (2.d0*pi*frq*tau_2)**2.)

      ! Compute the imaginary permittivitity coefficient (Eq 5)
      eps2 = 2.d0*pi*frq * (term1_p1 + term2_p1)
 
      eps = complex(eps1, eps2)

!       ! Compute the mass absorption coefficient (Eq 1)
!       RE = (epsilon-1)/(epsilon+2)
!       alpha = 6.d*!dpi*IMAGINARY(RE)*frq*1d-3/cl

      return
  
  end subroutine eps_water_tkc
end module eps_water

