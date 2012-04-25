module tmat_db

    use kinds

    implicit none

    integer :: i

    integer, parameter :: nf = 10, & ! number of frequencies in db
                          nt = 21, & ! number of temperatures in db
                          na = 3     ! number of aspect ratios in db
!    integer, parameter :: nf = 1, & ! number of frequencies in db
!                          nt = 23, & ! number of temperatures in db
!                          na = 1     ! number of aspect ratios in db

    integer, parameter :: nbins = 50 ! number of size bins

    integer, parameter :: nscat_records = 2048, &
                          next_records = 64, &
                          nemis_records = 32

    integer, parameter :: lbin_data = (nscat_records+next_records+nemis_records)

    integer :: ioldble ! length of a double precission record in a processor dependent manner
    integer :: lrec    ! length of a record we want to write

    integer :: rec_poslt, &
               rec_posut, &
               db_unit, &
               db_status

    real(kind=sgl), dimension(nf), parameter :: &
        db_freqs = (/176.311, 178.811, 180.311, 181.511, 182.311, 184.311, 185.111, 186.311, 187.811, 190.311/)
!        db_freqs = (/150.0/)

    real(kind=sgl), dimension(na), parameter :: &
        db_as = (/0.3, 0.5, 0.7/)
!        db_as = (/0.4/)

    real(kind=sgl), dimension(nt) :: &
        db_temps = (/(real(i), i=235,275,2)/)
!        db_temps = (/(real(i), i=251,273)/)


    character(40) :: db_file

    contains

! initialize the data base (file)

        subroutine initialize

        implicit none

        real(kind=dbl) :: test_iolength   ! double precission real to test the iolength
        db_unit = 25

        inquire(iolength=ioldble) test_iolength
        lrec = ioldble*lbin_data

        open(db_unit, file=trim(db_file), form='unformatted', access='direct', recl=lrec, iostat=db_status)

        if (db_status == 0) print*, 'database ', db_file, ' opened successfully!'

        return

        end subroutine initialize

! get the data

        subroutine get_data(type,f,t,ys)

        implicit none

        integer :: if,it,ir

        real(kind=sgl), intent(in) :: f,t

        real(kind=dbl), dimension(lbin_data) :: bound_lt,&
                                                bound_ut

        real(kind=dbl) :: xs, xl, xu
        real(kind=dbl), dimension(lbin_data) ::  yu, yl
        real(kind=dbl), dimension(lbin_data), intent(out) :: ys

        character(1), intent(in) :: type

        print*, "getting data for: ", type,f,t

        call check_input(f,t)

        call locate(db_freqs,nf,f,if)
        call locate(db_temps,nt,t,it)

        ! it lower temperature position
        ! nbins number of size bins
        ! lrec length of one record containing scattering and extinction matrices and emission vector
        ! for on frequency, temperature, aspectio ratio, and size bin

        ir = 1
!        rec_pos = (if-1)*(na*nt*nbins)+(ia-1)*(nt*nbins)+(it-1)*nbins+ir
        rec_poslt = (if-1)*(nt*nbins)+(it-1)*nbins+ir
        rec_posut = (if-1)*(nt*nbins)+(it)*nbins+ir

        print*, if, db_freqs(if)
        print*, it, db_temps(it), db_temps(it+1)
        print*, rec_poslt, rec_posut, lbin_data
        read(db_unit,rec=rec_poslt) bound_lt
        print*, bound_lt(1)
        read(db_unit,rec=rec_posut) bound_ut
        print*, bound_ut(1)

        xs = t
        xl = db_temps(it)
        xu = db_temps(it+1)
        yl = bound_lt
        yu = bound_ut
        ys = (xs-xl)*(yu-yl)/(xu-xl)+yl

        return

        end subroutine get_data

        subroutine check_input(f,t)

        implicit none

        real(kind=sgl), intent(in) :: f,t

        if (f .lt. db_freqs(1) .or. f .gt. db_freqs(nf)) then
            print*, 'Frequency ', f, ' out of range!'
            stop
        elseif (t .lt. db_temps(1) .or. t .gt. db_temps(nt)) then
            print*, 'Temperature ', t, ' out of range!'
            stop
        end if

        return

        end subroutine check_input

!        subroutine interpolation(nx1,nx2,x1,y1,x2,y2)
!
!          implicit none
!
!          integer :: i
!          integer :: nx1,nx2
!
!          integer:: ix2
!
!          real, intent(in), dimension(nx1) :: x1,y1
!          real, intent(in), dimension(nx2) :: x2
!          real, intent(out), dimension(nx2) :: y2
!
!          ix2 = 0
!
!          do i = 1, nx2
!             call locate(x1,nx1,x2(i),ix2)
!             y2(i) = (x2(i)-x1(ix2))*(y1(ix2+1)-y1(ix2))/(x1(ix2+1)-x1(ix2))+y1(ix2)
!          end do
!
!          return
!
!        end subroutine interpolation

! close the data base (file)

        subroutine close_db

        implicit none

        close(db_unit, iostat=db_status)

        if (db_status == 0) print*, "database closed successfully!"

        return

        end subroutine close_db
end module tmat_db
