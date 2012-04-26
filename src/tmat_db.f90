module tmat_db

    use kinds

    implicit none

    integer :: itemp

    integer, parameter :: nf = 11, & ! number of frequencies in db
                          nt = 21, & ! number of temperatures in db
                          na = 3     ! number of aspect ratios in db
!    integer, parameter :: nf = 1, & ! number of frequencies in db
!                          nt = 23, & ! number of temperatures in db
!                          na = 1     ! number of aspect ratios in db

    integer, parameter :: dbbins = 50 ! number of size bins

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
        db_freqs = (/150.00,176.311, 178.811, 180.311, 181.511, 182.311, 184.311, 185.111, 186.311, 187.811, 190.311/)
!        db_freqs = (/150.0/)

    real(kind=sgl), dimension(na), parameter :: &
        db_as = (/0.3, 0.5, 0.7/)
!        db_as = (/0.4/)

    real(kind=sgl), dimension(nt) :: &
        db_temps = (/(real(itemp), itemp=235,275,2)/)
!        db_temps = (/(real(i), i=251,273)/)


    character(3) :: as_str
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

        subroutine get_data(f,t,ir,s_mat,e_mat,e_vec)

        implicit none

        integer :: if,it,ir,id,l,l1,j1,j2,j

        integer, parameter :: nquad = 16,&
                              nstokes = 2

        real(kind=sgl), intent(in) :: f,t

        real(kind=dbl), dimension(lbin_data) :: bound_lt,&
                                                bound_ut

        real(kind=dbl) :: xs, xl, xu
        real(kind=dbl), dimension(lbin_data) ::  yu, yl
        real(kind=dbl), dimension(lbin_data) :: ys

        real(kind=dbl), dimension(nstokes,nquad,nstokes,nquad,4), intent(out) :: s_mat
        real(kind=dbl), dimension(nstokes,nstokes,nquad,2), intent(out) :: e_mat
        real(kind=dbl), dimension(nstokes,nquad,2), intent(out) :: e_vec

!        print*, "getting data for: ", f,t

        call check_input(f,t)

        call locate_db(db_freqs,nf,f,if)
        call locate_db(db_temps,nt,t,it)

        ! it lower temperature position
        ! dbbins number of size bins
        ! lrec length of one record containing scattering and extinction matrices and emission vector
        ! for on frequency, temperature, aspectio ratio, and size bin

!        rec_pos = (if-1)*(na*nt*dbbins)+(ia-1)*(nt*dbbins)+(it-1)*dbbins+ir
        rec_poslt = (if-1)*(nt*dbbins)+(it-1)*dbbins+ir
        rec_posut = (if-1)*(nt*dbbins)+(it)*dbbins+ir

!        print*, if, db_freqs(if)
!        print*, it, db_temps(it), db_temps(it+1)
!        print*, rec_poslt, rec_posut, lbin_data
        read(db_unit,rec=rec_poslt) bound_lt
!        print*, bound_lt(1)
        read(db_unit,rec=rec_posut) bound_ut
!        print*, bound_ut(1)

        xs = t
        xl = db_temps(it)
        xu = db_temps(it+1)
        yl = bound_lt
        yu = bound_ut
        ys = (xs-xl)*(yu-yl)/(xu-xl)+yl

        id = 1
        do l = 1, 2
            l1 = (l-1)*2 + 1
            do j1 = 1, nquad
                do j2 = 1, nquad
                    s_mat(1,j2,1,j1,l1) = ys(id)
                    s_mat(1,j2,2,j1,l1) = ys(id+1)
                    s_mat(2,j2,1,j1,l1) = ys(id+2)
                    s_mat(2,j2,2,j1,l1) = ys(id+3)
                    id = id + 4
                end do
            end do
        end do
        s_mat(:,:,:,:,4) = s_mat(:,:,:,:,1)
        s_mat(:,:,:,:,2) = s_mat(:,:,:,:,3)

        do j = 1, nquad
             e_mat(1,1,j,1) = ys(id)
             e_mat(1,2,j,1) = ys(id+1)
             e_mat(2,1,j,1) = ys(id+2)
             e_mat(2,2,j,1) = ys(id+3)
            id = id + 4
        enddo
        e_mat(:,:,:,2) = e_mat(:,:,:,1)

        do j = 1, nquad
            e_vec(1,j,1) = ys(id)
            e_vec(2,j,1) = ys(id+1)
            id = id + 2
        enddo
        e_vec(:,:,2) = e_vec(:,:,1)
        print*, s_mat(1,1,1,1,1), e_mat(1,1,1,1), e_vec(2,16,2)

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

!**********************************************************************
!    calculation location for interpolation                           *
!    copr. 1986-92 numerical recipes software +>k-5v1`                *
!                                                                     *
!**********************************************************************
!
subroutine locate_db(xx, n, x, j)
    !    given an array xx(1:n) and given a value x, returns a value j such
    !    that x is between xx(j) and xx(j+1) xx(1:n) must be monotonic, either
    !    increasing or decreasing. j=0 or j=n is returned to indicate
    !    that x is out of range

    integer j, n
    real x, xx (n)
    integer jl, jm, ju
              !initialise lower and
    jl = 0
              ! upper boundaries
    ju = n + 1
                            !if we are nmot yet done
10  if (ju - jl.gt.1) then
                          !compute a midpoint
        jm = (ju + jl) / 2
        if ( (xx (n) .ge.xx (1) ) .eqv. (x.ge.xx (jm) ) ) then
                   !and replace either the lower
            jl = jm
        else
                    !or the upper limit
            ju = jm
        endif
        goto 10
    endif
    if (x.eq.xx (1) ) then
        j = 1
    elseif (x.eq.xx (n) ) then
        j = n - 1
    else
        j = jl
    endif
    return
end subroutine locate_db

! close the data base (file)

        subroutine close_db

        implicit none

        close(db_unit, iostat=db_status)

        if (db_status == 0) print*, "database closed successfully!"

        return

        end subroutine close_db
end module tmat_db
