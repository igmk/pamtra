subroutine parse_options(gitVersion,gitHash)

    use kinds, only: long
    use getopt_m
    use settings, only: nfrq, freqs,frqs_str,namelist_file,&
      input_pathfile,input_file,verbose, descriptor_file_name, output_path
!     use vars_profile, only: coords
    use mod_io_strings, only: formatted_frqstr

    implicit none

    integer(kind=long) :: pos1, pos2, n
    integer(kind=long) :: ff
    character(40) :: gitHash, gitVersion
    type(option_s):: opts(7)
    opts(1) = option_s( "namelist", .true.,  'n' )
    opts(2) = option_s( "profile", .true.,  'p' )
    opts(3) = option_s( "descriptor", .true.,  'd' )    
!     opts(4) = option_s( "grid", .true.,  'g' )
    opts(4) = option_s( "output", .true., 'o' )
    opts(5) = option_s( "freqs", .true., 'f' )
    opts(6) = option_s( "verbose", .true., 'v' )
    opts(7) = option_s( "help", .false., 'h' )

    namelist_file = 'None'
    input_pathfile = 'profile/example_input.lay'
    output_path = "."
!     coords = (/1,1,1,1/)
    frqs_str = ''
    frqs_str(1) = '89.0'
    nfrq = 1
    descriptor_file_name = "descriptorfiles/descriptor_file.txt"
    
    do
        select case( getopt( "n:cp:cd:co:cf:cv:ch", opts ))
            case( char(0))
                exit
            case( 'n' )
                namelist_file = trim(optarg)
            case( 'p' )
                input_pathfile = trim(optarg)
            case( 'o' )
                output_path = trim(optarg)
            case( 'd' )
                descriptor_file_name = trim(optarg)                
!             case( 'g' )
!                 if (optarg(len_trim(optarg):) .ne. ',') &
!                 optarg = trim(optarg)//','
!                 call process_grid(optarg,coords)
            case( 'f' )
                if (optarg(len_trim(optarg):) .ne. ',') &
                optarg = trim(optarg)//','
                !                nf = countsubstring(optarg,',')
                !                allocate(freqs(nf))
                call process_freq(optarg)
            case( 'v' )
                read(optarg,'(i2)') verbose
            case( '?' )
                print *, 'unknown option ', optopt
                stop
            case('h')
                print*,'Usage: pamtra [options]'
                print*,''
                print*,'Available options:'
                print*,'   -n|--namelist     namelist file (default None)'
                print*,'   -p|--profile      profile file  (default profile/standard.dat)'
                print*,'   -o|--output       output directory  (default .)'
                print*,'   -d|--descriptor   descriptor file  (default descriptor_file.txt)'
!                 print*,'   -g|--grid         start_lon,end_lon,start_lat,end_lat (number of grid point)'
                print*,'   -f|--freqs        comma seperated list of frequencies (no blanks) (default 89.0)'
                print*,'   -v|--verbose      integer specifying verbose level between -1 (required by parallel python)'
                print*,'                     and 4 (default 0)'
                print*,'   -h|--help         print this help'
                print*,''
                print*,'Example: ./pamtra -v 1 -p rt_comp_single.dat -n run_params.nml -f 35,94,'
                print*,'See namelist file for further pamtra options'
                print*,''
                print*,'Version:  '//gitVersion
                print*,'Git Hash: '//gitHash
                stop
            case default
                print *, 'unhandled option ', optopt, ' (this is a bug)'
        end select

    end do
    ! the frequency string needs to be converted to a real array

    do ff = 1, nfrq
        read(frqs_str(ff),*) freqs(ff)
        frqs_str(ff) = formatted_frqstr(frqs_str(ff))
    end do
    
 
    !get filename
    pos1 = 1
    n = 0
    DO
      pos2 = INDEX(input_pathfile(pos1:), "/")
      IF (pos2 == 0) THEN
        n = n + 1
        input_file = input_pathfile(pos1:)
        EXIT
      END IF
      n = n + 1
      input_file = input_pathfile(pos1:pos1+pos2-2)
      pos1 = pos2+pos1
  END DO

contains

    subroutine process_freq(arg)

        use kinds, only: long

        implicit none

        integer(kind=long) :: ind
        character(len=*), intent(in) :: arg
        character(500) :: arg_loc


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

!     subroutine process_grid(arg,coords)
! 
!         use kinds, only: long
! 
!         implicit none
! 
!         integer(kind=long) :: cc, ind
!         character(len=*), intent(in) :: arg
!         character(150) :: arg_loc
! 
!         integer(kind=long), dimension(4), intent(out) :: coords
! 
!         cc = 1
!         arg_loc = arg
!         do while (arg_loc .ne. "")
!             ind = index(arg_loc,",")
!             read(arg_loc(1:ind-1),*) coords(cc)
!             cc = cc+1
!             arg_loc = trim(arg_loc(ind+1:))
!         end do
! 
!         return

!     end subroutine process_grid
end subroutine parse_options
