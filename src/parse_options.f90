subroutine parse_options(gitVersion,gitHash,frqs_str,nfrq)

    use kinds
    use getopt_m
    use nml_params, only: maxfreq
    use file_mod, only: namelist_file, input_file
    use vars_profile, only: coords

    implicit none

    integer :: nfrq
    !    real(kind=dbl), allocatable,dimension(:) :: freqs
    character:: ch
    character(40) :: gitHash, gitVersion
    character(8), dimension(maxfreq) :: frqs_str
    type(option_s):: opts(5)
    opts(1) = option_s( "namelist", .true.,  'n' )
    opts(2) = option_s( "profile", .true.,  'p' )
    opts(3) = option_s( "grid", .true.,  'g' )
    opts(4) = option_s( "freqs", .true., 'f' )
    opts(5) = option_s( "help", .false., 'h' )

    namelist_file = 'run_params.nml'
    input_file = 'standard.dat'
    coords = (/1,1,1,1/)
    frqs_str = ''
    frqs_str(1) = '89.0'
    nfrq = 1

    do
        select case( getopt( "n:cp:cg:cf:ch", opts ))
            case( char(0))
                exit
            case( 'n' )
                namelist_file = optarg
            case( 'p' )
                input_file = optarg
            case( 'g' )
                if (optarg(len_trim(optarg):) .ne. ',') &
                optarg = trim(optarg)//','
                call process_grid(optarg,coords)
            case( 'f' )
                if (optarg(len_trim(optarg):) .ne. ',') &
                optarg = trim(optarg)//','
                !                nf = countsubstring(optarg,',')
                !                allocate(freqs(nf))
                call process_freq(optarg,frqs_str,nfrq)
            case( '?' )
                print *, 'unknown option ', optopt
                stop
            case('h')
                print*,'Usage: pamtra [options]'
                print*,''
                print*,'Available options:'
                print*,'   -n|--namelist     namelist file (default run_params.nml)'
                print*,'   -p|--profile      profile file  (default standard.dat)'
                print*,'   -g|--grid         start_lon,end_lon,start_lat,end_lat (number of grid point)'
                print*,'   -f|--freqs        comma seperated list of frequencies (no blanks) (default 89.0)'
                print*,'   -h|--help         print this help'
                print*,''
                print*,'Example: ./pamtra -p rt_comp_single.dat -n run_params.nml -f 35,94,'
                print*,'See namelist file for further pamtra options'
                print*,''
                print*,'Version:  '//gitVersion
                print*,'Git Hash: '//gitHash
                stop
            case default
                print *, 'unhandled option ', optopt, ' (this is a bug)'
        end select
    end do
contains

    subroutine process_freq(arg,frqs_str,nfrq)

        use nml_params, only: maxfreq

        implicit none

        integer :: nfrq, ind
        character(len=*), intent(in) :: arg
        character(900) :: arg_loc
        character(8), dimension(maxfreq) :: frqs_str

        nfrq = 0
        arg_loc = arg
        do while (arg_loc .ne. "")
            !    print*, arg_loc
            nfrq = nfrq+1
            ind = index(arg_loc,",")
            !    print*, ind, arg_loc(1:ind-1)
            !    read(arg_loc(1:ind-1),*) freqs(ff)
            frqs_str(nfrq) = arg_loc(1:ind-1)
            !    print*, freqs(ff)
            arg_loc = trim(arg_loc(ind+1:))
        !    print*, arg_loc
        end do

        return

    end subroutine process_freq

    subroutine process_grid(arg,coords)

        integer :: cc, ind
        character(len=*), intent(in) :: arg
        character(150) :: arg_loc

        integer, dimension(4), intent(out) :: coords

        cc = 1
        arg_loc = arg
        do while (arg_loc .ne. "")
            ind = index(arg_loc,",")
            read(arg_loc(1:ind-1),*) coords(cc)
            cc = cc+1
            arg_loc = trim(arg_loc(ind+1:))
        end do

        return

    end subroutine process_grid
end subroutine parse_options
