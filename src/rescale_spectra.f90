module rescale_spec
    !
    ! Description:
    !rescale a particle spectrum to another size interval.
    !
    ! Current Code Owner: IGMK
    !
    ! History:
    !
    ! Version   Date       Comment
    ! -------   ----       -------
    ! 0.1       28/11/2012 created M.Maahn
    ! 0.2       15/04/2013 Application of European Standards for Writing and
    !                      Documenting Exchangeable Fortran 90 Code - M. Maahn

    ! Code Description:
    !   Language:		Fortran 90.
    !   Software Standards: "European Standards for Writing and  
    !     Documenting Exchangeable Fortran 90 Code". 
    !


contains


  subroutine rescale_spectra(&
    errorstatus,&     	! out
    nx1,& 	!in
    nx2,& 	!in
    sort,& 	!in
    x1,& 	!in
    y1,& 	!in
    x2,& 	!in
    y2)  	!out
    !(c) M.Maahn, IGMK, 11/2012

    !nx1,in: length of x1, y1
    !nx2,in: length of x2, y2
    !sort: logical: sort x1 and y1 before interpolation?
    !x1, in: x-values original data
    !y1, in: y-values original values
    !x2, in: center of new bins
    !y2, out: result

    !this routine averages and interpolates for every interval depending what is needed
    !in praxis, at the edges of the interval, values are interpolated  
    ! (thus we can be sure, that there are at least 2 values in every interval)
    !interpolated and original values are combined and sorted
    !then we average each interval.

  !  use settings, only: verbose
    use kinds
    use report_module
    implicit none


    integer, intent(in) :: nx1,nx2
    real(kind=dbl), intent(in), dimension(nx1) :: x1,y1
    real(kind=dbl), intent(in), dimension(nx2) :: x2
    real(kind=dbl), intent(out), dimension(nx2) :: y2
    logical, intent(in) :: sort
    real(kind=dbl), dimension(nx1) :: x1_sorted,y1_sorted
    real(kind=dbl), dimension(nx2+1) :: x2_shift

    real(kind=dbl), dimension(nx1+nx2+1) :: x12,x12_sorted
    real(kind=dbl), dimension(nx1+nx2+1) :: y12,y12_sorted
    real(kind=dbl), dimension(nx2+1) :: y2_interp

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=15) :: nameOfRoutine = 'rescale_spectra'

    interface
      subroutine dsort (errorstatus, dx, dy, n, kflag)
	use kinds
	implicit none
	integer(kind=long), intent(out) :: errorstatus	
	real(kind=dbl), dimension(n), intent(inout) :: dx, dy
	integer, intent(in) :: n, kflag
      end subroutine dsort
    end interface

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

    call assert_true(err,nx1>1,&
        "nx1 must be greater 1") 
    call assert_true(err,nx2>1,&
        "nx2 must be greater 1")   
    call assert_true(err,maxval(x1)<=maxval(x2),&
        "max x1 must be smaller than x2")      
    call assert_true(err,minval(x1)>=minval(x2),&
        "min x1 must be larger than x2")    
    if (err /= 0) then
      msg = 'error in tranforming the spectrum to velocity space...'

      print*, "nx1, nx2, maxval(x1), maxval(x2), minval(x1), minval(x2)"
      print*, nx1, nx2, maxval(x1), maxval(x2), minval(x1), minval(x2) 
      call report(err, msg, nameOfRoutine)
      errorstatus = err
      return
    end if


    if (all(y1==0)) then
      msg = "all values zero"
      call report(warning, msg, nameOfRoutine)
    end if

    if (verbose >= 7) then
      print*, "IN X", x1
      print*, "IN Y", y1
    end if
    
    x1_sorted = x1
    y1_sorted = y1

    if (sort) then
      call dsort(err,x1_sorted, y1_sorted, nx1, 2)
      if (err /= 0) then
	  msg = 'error in dsort!'
	  call report(err, msg, nameOfRoutine)
	  errorstatus = err
	  stop !return
      end if       
    end if


    x2_shift(2:nx2) = x2(2:) - 0.5d0*(x2(2:) - x2(1:nx2-1))
    x2_shift(1) = x2(1) -  0.5d0*(x2(2) - x2(1))
    x2_shift(nx2+1) = x2(nx2) +  0.5d0*(x2(nx2) - x2(nx2-1))  


    call interpolate_spectra(err,nx1,nx2+1,x1_sorted,y1_sorted,x2_shift,y2_interp)
    if (err /= 0) then
	msg = 'error in interpolate_spectra!'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	return
    end if

    !join interpolated and original array
    x12(1:nx1)=x1_sorted
    x12(nx1+1:nx1+nx2+1)=x2_shift
    y12(1:nx1)=y1_sorted
    y12(nx1+1:nx1+nx2+1)=y2_interp

    !make order right
    x12_sorted = x12
    y12_sorted=y12
    call dsort(err,x12_sorted, y12_sorted, nx1+nx2+1, 2)
    if (err /= 0) then
	msg = 'error in dsort!'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	stop !return
    end if     

    call average_spectra(err,nx1+nx2+1,nx2+1,x12_sorted,y12_sorted,x2_shift,y2)
    if (err /= 0) then
	msg = 'error in average_spectra!'
	call report(err, msg, nameOfRoutine)
	errorstatus = err
	return
    end if

    if (all(y2==0)) then
      msg = "all calculated values zero"
      call report(warning, msg, nameOfRoutine)
    end if

    if (verbose >= 7) then
      print*, "OUT X", x2
      print*, "OUT Y", y2
    end if    
    
    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    errorstatus = err
    return

  end subroutine rescale_spectra


  subroutine average_spectra&
    (errorstatus,&     	! out
    nx12,& !in
    nx2,& !in
    x12_sorted,& !in
    y12_sorted,& !in
    x2,& !in
    y_result) !out
    !averages the spectrum
    !borders of the averaged intervalls must be already present in in x12_sorted!
    !works only in combination with average_spectra

    !(c) M.Maahn, IGMK, 11/2012

    use settings, only: verbose
    use kinds
    use report_module
    implicit none

    integer, intent(in) :: nx12,nx2
    real(kind=dbl), intent(in), dimension(nx12) :: x12_sorted,y12_sorted
    real(kind=dbl), intent(in), dimension(nx2) :: x2
    real(kind=dbl), intent(out), dimension(nx2-1) :: y_result

    integer :: ii,jj1,jj2

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=15) :: nameOfRoutine = 'average_spectra'

    interface
      SUBROUTINE locate (xx, n, x, j) 
	use kinds                                                                 
	INTEGER j, n 
	REAL(kind=dbl) x, xx (n) 
      end SUBROUTINE locate 
    end interface

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    if (all(y12_sorted==0)) then
      msg = "all input values zero"
      call report(warning, msg, nameOfRoutine)
    end if

    !step zero
    call locate (x12_sorted, nx12, x2(1), jj1)
    do ii=1,nx2-1
      !find indices
      call locate (x12_sorted, nx12, x2(ii+1), jj2) 
      !locate does not work properly if an entry is searched EQUAL to the last
      if (jj1 .eq. jj2) then
	jj2 = jj1 +1 
      end if

      if ((jj2 .gt. nx12) .or. (jj1 .gt. nx12)) then
	  errorstatus = fatal
	  msg = 'Exiting averaging loop'
	  call report(errorstatus, msg, nameOfRoutine)
	  return
      else
	  err = success
      end if

      !make the averaging, first width half of the weights on the left side
      y_result(ii) = SUM(y12_sorted(jj1:jj2-1) * 0.5d0*(x12_sorted(jj1+1:jj2)-x12_sorted(jj1:jj2-1)  ) ) 
      !now right-side weights
      y_result(ii) = y_result(ii) + SUM(y12_sorted(jj1+1:jj2) * 0.5d0*(x12_sorted(jj1+1:jj2)-x12_sorted(jj1:jj2-1)  ) ) 
      !devide by weights
      y_result(ii) = y_result(ii) / SUM((x12_sorted(jj1+1:jj2)-x12_sorted(jj1:jj2-1)))
      !save idnex for next iteration
      jj1=jj2

    end do

    if (all(y_result==0)) then
      msg = "all calculated values zero"
      call report(warning, msg, nameOfRoutine)
    end if

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

    errorstatus = err
    return

  end subroutine average_spectra

  subroutine interpolate_spectra(&
    errorstatus,&     	! out
    nx1,& !in
    nx2,& !in
    x1,& !in
    y1,& !in
    x2,& !in
    y2)  !out

    !interpolation which gives back ZERO if x2 is more than 1bin out of range of x1!
    !values MUST be sorted

    use settings, only: verbose
    use kinds
    use report_module
    implicit none

    integer :: i
    integer :: nx1,nx2

    integer:: ix2

    real(kind=dbl), intent(in), dimension(nx1) :: x1,y1
    real(kind=dbl), dimension(nx1+2) :: x1_ext,y1_ext
    real(kind=dbl), intent(in), dimension(nx2) :: x2
    real(kind=dbl), intent(out), dimension(nx2) :: y2

    integer(kind=long), intent(out) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=15) :: nameOfRoutine = 'interpolate_spectra'

    if (verbose >= 3) call report(info,'Start of ', nameOfRoutine)

    !extend x1 and y1 to give a "0" reference point
    x1_ext(1)     = x1(1) - (x1(2)-x1(1))
    x1_ext(nx1+2) = x1(nx1) +  (x1(2)-x1(1))
    y1_ext(:) = 0.d0 

    x1_ext(2:nx1+1) = x1
    y1_ext(2:nx1+1) = y1

    ix2 = 0
    y2(:) = 0.d0

    do i = 1, nx2
      call locate(x1_ext,nx1+2,x2(i),ix2)
      !points our of range?
      if ((ix2 .eq. 0) .or. (ix2 .eq. nx1+2)) then
	y2(i) = 0.d0
      else
	y2(i) = (x2(i)-x1_ext(ix2))*(y1_ext(ix2+1)-y1_ext(ix2))/(x1_ext(ix2+1)-x1_ext(ix2))+y1_ext(ix2)
      end if
    if (verbose .gt. 5) print*, "y2(i),i",y2(i),i
    end do

    if (verbose >= 3) call report(info,'End of ', nameOfRoutine)

    errorstatus = err
    return

  end subroutine interpolate_spectra
end module rescale_spec