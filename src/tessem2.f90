subroutine tessem2(errorstatus,freq,zenith_angle,wind_speed,temperature,salinity,emissivity,reflectivity)

    ! inputs: - freq (GHz) between 10 and 700 GHz
    !         - angle from nadir (degree) 
    !         - surface wind speed at 10 m (m/s) between 0 and 25m/s
    !         - sea surface temperature (K) between 270 and 310 K
    !         - salinity (deg/degdeg) / typically 36 (between 0 and 40deg/degdeg)
    ! ouputs: - sea surface emissivity for the given orthogonal polarization

    use kinds
    use report_module
    use nan, only: nan_dbl
    
    implicit none

    integer(kind=long) :: i,j,pol
    
    real(kind=dbl), intent(in) :: freq, zenith_angle, wind_speed, temperature, salinity
    real(kind=dbl), intent(out), dimension(2) :: emissivity, reflectivity 
    integer(kind=long), parameter :: size_mmax=15

    real(kind=dbl), dimension(2,6,size_mmax) :: nn_params ! (pol, nodes, number) b1, b2, x_min, x_max, y_min, y_max
    real(kind=dbl), dimension(2,2,size_mmax,size_mmax) :: w_params ! (pol, 2, nodes, number))
    real(kind=dbl), dimension(15,5) :: w_tmp
    
    real(kind=dbl), dimension(5) :: x
    real(kind=dbl), dimension(2) :: y
    real(kind=dbl), dimension(2,15) :: trans, new_x, new_y

    ! Error handling

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=14) :: nameOfRoutine = 'tessem2'

    nn_params = nan_dbl()
    w_params = nan_dbl()
    
    nn_params(1,1,1:15) = (/-4.425637,-9.833792,1.338131,-1.688713,0.121844,&
                        -0.020075,0.983048,-5.365536,0.759206,1.705090,&
                        -0.861439,2.639012,4.438483,-1.376479,0.920075/)
    nn_params(1,2,1) = -3.126669
    w_tmp = reshape((/0.988378,-0.053936,-1.166637, 0.213664,-0.264162,&
        -0.946456,-0.084722,-2.169090, 1.038926, 1.814769,&
        -0.588050,1.553885, 3.672554, 0.198168, 0.520655,&
        -0.865033, 6.091850, 0.514930, 1.203841, 1.349278,&
        -0.875449,-1.582926, 1.683341, 1.679280, 0.533539,&
        1.047984, 0.534832, 0.429764, 1.978477,-0.982343,&
        0.024583,-2.832310,-0.249698,-0.264623,-0.385522,&
        0.197084, 0.071463,-1.025112,-0.365062, 0.099878,&
        0.308190,-0.034068, 0.011647, 0.158408,-0.512171,&
        -1.243063,-0.008161, 0.323026,-0.040724, 0.526443,&
        -1.044800, 0.035568, 0.057517,-0.559605,-0.086835,&
        0.058220, 0.284328, 0.235683,-0.050183,-0.074034,&
        2.514267, 0.002472,-0.249759,-0.017623, 2.265044,&
        1.006211,-0.000011,-0.010390, 0.285759,-0.024043,&
        -0.010873,-0.063544,-0.043210,-0.000275, 0.017342/),(/15,5/))
    w_params(1,1,1:15,1:5) = w_tmp
    w_params(1,2,1,1:15) = (/-0.105095,-5.132794,-0.068984,-1.133355,-0.016834,&
                        0.021879,-1.791876,2.727342,0.070811,0.556869,&
                        -1.172738,-1.725109,1.650515,-1.123763,-0.929365/)
    nn_params(1,3,1:5) = (/10.000000,0.000000,0.000000,273.149990,0.000000/)
    nn_params(1,4,1:5) = (/700.000000,88.000000,24.000000,309.149994,40.000000/)
    nn_params(1,5,1) = 0.272829
    nn_params(1,6,1) = 0.997013
    
    nn_params(2,1,1:15) = (/-0.492534,3.549479,1.209267,-1.201308,0.657236,&
                        0.522424,0.439041,1.070375,0.917859,2.999629,&
                        0.150532,-0.225454,-1.418973,-5.616032,3.077098/)
    nn_params(2,2,1) = -3.191476
    w_tmp = reshape((/0.914705,-0.244028, 0.217031, 0.022331,-0.757337,&
                    -0.514908,-0.682990, 1.319635, 1.052532, 0.327665,&
                    -0.197895, 0.529253,-1.846079,-5.075406, 3.362876,&
                    -0.907370,-1.694778,-0.628096, 2.230164, 0.003549,&
                    -0.290958,-0.013354, 0.122967, 0.001593,-0.061660,&
                    1.609950,-0.056635, 1.059906,-0.224954, 0.942301,&
                    1.806560, 0.821685,-1.057938, 0.082932, 0.045666,&
                    0.119414, 0.032902,-0.178454,-0.058426,-0.098223,&
                    0.291156,-0.036004,-0.524861, 0.089880, 0.268074,&
                    -0.237687, 0.040081,-0.019350, 0.028382, 0.831685,&
                    -1.168739, 0.715660, 0.498660,-0.716106, 0.916730,&
                    0.060666,-0.605959, 0.214584,-0.233417,-0.207224,&
                    0.013675,-0.000319, 0.008815,-0.004765, 0.823811,&
                    -0.244888, 0.265740, 1.596743,-2.614260,-2.705295,&
                    -0.003804, 0.105391,-0.008504, 0.003713,-0.018812/),(/15,5/))
    w_params(2,1,1:15,1:5) = w_tmp
    w_params(2,2,1,1:15) = (/0.019943,2.713515,-0.147582,-0.305444,-0.275328,&
                                    -0.054163,0.628523,0.028561,0.044241,0.034885,&
                                    -0.168618,0.580301,-0.157557,-0.586725,0.179877/)
    nn_params(2,3,1:5) = (/10.000000,0.000000,0.000000,273.149990,0.000000/)
    nn_params(2,4,1:5) = (/700.000000,88.000000,24.000000,309.149994,40.000000/)
    nn_params(2,5,1) = 0.018036
    nn_params(2,6,1) = 0.884545

    if (verbose >= 1) call report(info,'Start of ', nameOfRoutine)
    err = 0

    x = (/freq, zenith_angle, wind_speed,temperature,salinity/)
    
    do pol = 1,2
    ! preprocessing      
    do i=1,5
       new_x(pol,i)=-1+(x(i)-nn_params(pol,3,i))/(nn_params(pol,4,i)-nn_params(pol,3,i))*2
    enddo
    ! propagation
    do i=1,15
       trans(pol,i)=nn_params(pol,1,i);
       do j=1,5
        trans(pol,i)=trans(pol,i)+w_params(pol,1,i,j)*new_x(pol,j)
       enddo
       trans(pol,i)= 2/(1+exp(-2*trans(pol,i)))-1
    enddo
    do i=1,1
       new_y(pol,i)=nn_params(pol,2,i)
       do j=1,15
          new_y(pol,i)=new_y(pol,i)+w_params(pol,2,i,j)*trans(pol,j)
       enddo
    enddo
    ! postprocessing      
       y(pol)=nn_params(pol,5,1)+(new_y(pol,1)+1)/2*(nn_params(pol,6,1)-nn_params(pol,5,1))
    end do

    emissivity = y
    reflectivity = 1. - y
    
    if (err > 0) then
      errorstatus = fatal
      msg = "error in calculations from tessem2"
      call report(errorstatus, msg, nameOfRoutine)
      return
    end if 

    errorstatus = err
    if (verbose >= 1) call report(info,'End of ', nameOfRoutine)

end subroutine tessem2
