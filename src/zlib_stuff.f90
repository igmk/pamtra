!
! Fortran interfaces for zlib and for POSIX pipes (popen and friends)
! Author:  David L Duffy 2009-2014
! 
! Copyright (c) 2014 David L Duffy
! 
! This software is provided 'as-is', without any express or implied
! warranty. In no event will the authors be held liable for any damages
! arising from the use of this software.
! 
! Permission is granted to anyone to use this software for any purpose,
! including commercial applications, and to alter it and redistribute it
! freely, subject to the following restrictions:
! 
!    1. The origin of this software must not be misrepresented; you must not
!    claim that you wrote the original software. If you use this software
!    in a product, an acknowledgment in the product documentation would be
!    appreciated but is not required.
! 
!    2. Altered source versions must be plainly marked as such, and must not be
!    misrepresented as being the original software.
! 
!    3. This notice may not be removed or altered from any source
!    distribution.
!
! AS182 was published in the journal Applied Statistics, under the copyright of the Royal
! Statistical Society, and I understand was thus made available subject 
! to the restriction that no fee is charged for redistribution
!
! C preprocessor directives support various Fortran compilers. Where
! iso_c_binding is not available, gzipped files can be decompressed
! using the gzip executable (providing the Fortran compiler supports the
! SYSTEM command)
!
! The readline subroutine may fail reading very large files
!
! 20140106: initial release
! 20140110: added in missing iocodes module, fixed typos 
!           pointed out by John David and yeonpil
!
! iostat codes needed for wrinline etc
!

#define ZLIB 1
#define OPEN64 0
#define POPEN 1


module iocodes
#if OPEN64
  integer, parameter :: eofcode = -4001
  integer, parameter :: eolcode = -4006
  character (len=10), parameter :: stream_access = 'sequential'
  character (len=6), parameter :: stream_form = 'binary'
#else
  integer, parameter :: eofcode = -1
  integer, parameter :: eolcode = -2
  character (len=6), parameter :: stream_access = 'stream'
  character (len=11), parameter :: stream_form = 'unformatted'
#endif
end module iocodes
!
! Definition of a port
!   slots: associated file name
!          1=uncompressed 2=gzipped 3=unzipped copy 4=pipe
!          Fortran style logical unit number
!          gzip C-style file handle
!
module ioports
#if !OPEN64
  use, intrinsic :: iso_c_binding
#endif
  integer, parameter :: PORT_STANDARD = 1, PORT_GZIPPED = 2,  &
                        PORT_COPY = 3, PORT_PIPE = 4
  type, public :: ioport
    character (len=256) :: filnam
    character (len=1) :: stat = ' '
    integer :: filtyp
    integer :: fstream
#if !OPEN64
    type (c_ptr) :: handle = c_null_ptr
#endif
  end type ioport
end module ioports
!
! Fortran interface to zlib 
!   based on looking at fgzlib, fgsl and Janus Weil's example
!   on comp.lang.fortran May 2009
!   currently enough functionality to read gzipped text files
!
#if ZLIB
module f95zlib
  use, intrinsic :: iso_c_binding
  use ioports
! buffer for gzread
  integer, parameter :: ZBUFLEN = 65536
  character (len=ZBUFLEN), target :: zbuffer
! current character and end of zbuffer
  integer :: zbufpos=0, zbufend=ZBUFLEN
! gzopen  
  interface
    function gzopen(path, mode) bind(C, name='gzopen')
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: path, mode
      type (c_ptr) :: gzopen
    end function
  end interface
! gzread  
  interface
    function gzread(filehandle, buf, len) bind(C, name='gzread')
      use, intrinsic :: iso_c_binding
      integer (c_int) :: gzread
      type (c_ptr), value :: filehandle
      type (c_ptr), value :: buf
      integer(c_int), value :: len
    end function
  end interface
! gzwrite
  interface
    function gzwrite (filehandle, buf, len) bind(C, name='gzwrite')
      use, intrinsic :: iso_c_binding
      integer (c_int) :: gzwrite
      type (c_ptr), value :: filehandle
      type (c_ptr), value :: buf
      integer(c_int), value :: len
    end function
  end interface
! gzgetc  
  interface
    function gzgetc(filehandle) bind(C, name='gzgetc')
      use, intrinsic :: iso_c_binding
      integer (c_int) :: gzgetc
      type (c_ptr), value :: filehandle
    end function
  end interface
! gzrewind 
  interface
    function gzrewind(filehandle) bind(C, name='gzrewind')
      use, intrinsic :: iso_c_binding
      integer(c_int) :: gzrewind
      type (c_ptr), value :: filehandle
    end function
  end interface
! gzclose  
  interface
    function gzclose(filehandle) bind(C, name='gzclose')
    use, intrinsic :: iso_c_binding
    integer(c_int) :: gzclose
    type (c_ptr), value :: filehandle
    end function
  end interface
contains
!
! Wrapper for gzopen
!   also reinitializes gzread's buffer
!
  subroutine fgz_open(path, mode, fd, ios)
    use, intrinsic :: iso_c_binding
    character(kind=c_char, len=*), intent(in) :: path, mode
    type (ioport) :: fd
#if SUN
    character(kind=c_char, len=len_trim(path)+1) :: cpath
    character(kind=c_char, len=len_trim(mode)+1) :: cmode
    integer :: eos
#endif
    integer :: ios
    ios=0
    fd%filnam=path
    fd%filtyp=PORT_GZIPPED 
    fd%fstream=-1
#if SUN
    eos=len_trim(path)
    cpath=path
    cpath((eos+1):(eos+1))=c_null_char
    eos=len_trim(mode)
    cmode=mode
    cmode((eos+1):(eos+1))=c_null_char
    fd%handle = gzopen(cpath, cmode)
#else
    fd%handle = gzopen(trim(path) // c_null_char, trim(mode) // c_null_char)
#endif
    if (.not.c_associated(fd%handle)) ios=-1
    zbufpos=0
  end subroutine fgz_open
!
! Wrapper for gzrewind
!
  subroutine fgz_rewind(fd, ios)
    use, intrinsic :: iso_c_binding
    type(ioport) :: fd
    integer :: ios
    integer(c_int) :: ir
    ios = 0
    ir = gzrewind(fd%handle)
    if (ir /= 0) ios=ir
    zbufpos=0
  end subroutine fgz_rewind
!
! Wrapper for gzread
!   read one line of text from buffer
!
  subroutine fgz_read(fd, lin, advance, ios)
    use, intrinsic :: iso_c_binding
    use iocodes
    type(ioport) :: fd
    character(len=*) :: lin
    character(len=*), intent(in), optional :: advance
    integer, intent(out) :: ios

    integer :: i, j, linlen, nchar, newzpos, pos
    integer(c_int) :: blen, rlen
!
!  eol morez more
!   F    T    T    read buffer, copy to output
!   F    T    F    read buffer, output full
!   T    F    F    found <NL>
!  advancing
!   no             after output full, exit with buffer pos at end of text
!   yes            after output full, exit with buffer pos at next <NL>
!
    logical :: advancing, eol, more, morez

    type (c_ptr) :: buf = c_null_ptr

    advancing=.true.
    if (present(advance)) advancing=(advance == 'yes') 
    linlen=len(lin)
    ios=0
    lin=' '
    sta=1
    nchar=-1
    pos=0
    j=0
    eol=.false.
    more=.true.
    morez=.true.
    do while (morez)
      j=j+1
! refill buffer if necessary
      if (zbufpos == 0) then
        blen=ZBUFLEN
        buf=c_loc(zbuffer(1:1))
        rlen=gzread(fd%handle, buf, blen)
        if (rlen <= 0) then
          ios=-1
          return
        end if
        zbufpos=1
        zbufend=rlen
      end if
! place buffer index at <NL> or buffer end
! if <NL> will exit after updating output
      newzpos=zbufend+1
      nchar=zbufend-zbufpos+1
      do i=zbufpos, zbufend
        if (zbuffer(i:i) == achar(10)) then
          eol=.true.
          morez=.false.
          newzpos=i+1
          nchar=i-zbufpos
          exit
        end if
      end do
! read in min(buffer, remaining output)
! if not advancing move buffer idx back to last character read and exit
      if (more) then
        if (linlen < pos+nchar) then
          more=.false.
          nchar=linlen-pos
          if (.not.advancing) then
            newzpos=zbufpos+nchar
            morez=.false.
            eol=.false.
          end if
        end if
        lin((pos+1):(pos+nchar))=zbuffer(zbufpos:(zbufpos+nchar-1))
        pos=pos+nchar
      end if
      zbufpos=newzpos
      if (zbufpos > zbufend) then
        zbufpos=0
      end if
    end do
    if (.not.advancing .and. eol) ios=eolcode
  end subroutine fgz_read
!
! write one line of text to a gzipped textfile
!
  subroutine fgz_write(fd, lin, advance, ios)
    use, intrinsic :: iso_c_binding
    use iocodes
    type(ioport) :: fd
    character(len=*) :: lin
    character(len=*), intent(in), optional :: advance
    integer, intent(out) :: ios

    logical :: advancing
    integer :: ioerr, linlen, lsta, lpos
    integer(c_int) :: blen, wlen
    type (c_ptr) :: buf = c_null_ptr

    advancing=.true.
    if (present(advance)) then
      advancing=(advance == 'yes')
    end if
    ios=0
    lpos=0
    lenlin=len_trim(lin)
    do 
      lsta=lpos+1
      lpos=min(lsta+ZBUFLEN-1, lenlin)
      zbuffer=lin(lsta:lpos)
      buf=c_loc(zbuffer(1:1))
      blen=lpos-lsta+1
      wlen=gzwrite(fd%handle, buf, blen)
      ioerr=wlen
      if (ioerr == 0) exit
      if (lpos == lenlin) then
        if (advancing) then
          zbuffer=char(10)
          buf=c_loc(zbuffer(1:1))
          blen=1
          wlen=gzwrite(fd%handle, buf, blen)
          ioerr=wlen
        end if
        exit
      end if
    end do
    if (ioerr == 0) ios=-1
  end subroutine fgz_write
!
! Wrapper for gzclose
!
  subroutine fgz_close(fd, ios)
    use, intrinsic :: iso_c_binding
    type(ioport) :: fd
    integer :: ios
    integer(c_int) :: ic 
    ios = 0
    ic = gzclose(fd%handle)
    if (ic /= 0) ios = ic
  end subroutine fgz_close
end module f95zlib
#endif
!
! Fortran interface to popen
!
#if POPEN
module f95pipes
  use, intrinsic :: iso_c_binding
  use ioports
! popen  
  interface
    function popen(path, mode) bind(C, name='popen')
      use, intrinsic :: iso_c_binding
      character(kind=c_char), dimension(*) :: path, mode
      type (c_ptr) :: popen
    end function
  end interface
! fgets  
  interface
    function fgets(buf, siz, handle) bind(C, name='fgets')
      use, intrinsic :: iso_c_binding
      type (c_ptr) :: fgets
      character(kind=c_char), dimension(*) :: buf
      integer(kind=c_int), value :: siz
      type (c_ptr), value :: handle
    end function
  end interface
! fputs  
  interface
    function fputs(buf, handle) bind(C, name='fputs')
      use, intrinsic :: iso_c_binding
      integer (c_int) :: fputs
      character(kind=c_char), dimension(*) :: buf
      type (c_ptr), value :: handle
    end function
  end interface
! pclose  
  interface
    function pclose(handle) bind(C, name='pclose')
    use, intrinsic :: iso_c_binding
    integer(c_int) :: pclose
    type (c_ptr), value :: handle
    end function
  end interface
contains
! wrapper for popen
!   fd%stat gives mode
  subroutine pipe_open(command, fd, ios)
    use, intrinsic :: iso_c_binding
    character(*), intent(in) :: command
    type (ioport) :: fd
    integer :: ios
    ios=0
    fd%filnam=command
    fd%filtyp=PORT_PIPE
    fd%fstream=-1
    fd%handle =  popen(trim(command) // C_NULL_CHAR, fd%stat // C_NULL_CHAR)
    if (.not.c_associated(fd%handle)) ios=-1
  end subroutine pipe_open
!
! rewind pipe
  subroutine pipe_rewind(fd, ios)
    use, intrinsic :: iso_c_binding
    type (ioport) :: fd
    integer :: ios
    call pipe_close(fd, ios)
    if (ios == 0) call pipe_open(fd%filnam, fd, ios)
  end subroutine pipe_rewind
!
! wrapper for fgets
  subroutine pipe_read(fd, lin, advance, ios)
    use, intrinsic :: iso_c_binding
    type(ioport) :: fd
    character(len=*) :: lin
    character(len=*), intent(in), optional :: advance
    integer, intent(out) :: ios

    integer :: i, eos
    integer(c_int) :: clen
    ios=0
    clen=len(lin)
    lin=' '
    if (.not.c_associated(fgets(lin, clen, fd%handle))) then
      ios=-1
      return
    end if
    eos=2
    do i=1, clen  
      if (lin(i:i) == C_NULL_CHAR) then
        eos=i-2
        exit
      end if
    end do
    lin=lin(1:eos)
  end subroutine pipe_read
!
! wrapper for fputs
  subroutine pipe_write(fd, lin, advance, ios)
    use, intrinsic :: iso_c_binding
    type(ioport) :: fd
    character(len=*) :: lin
    character(len=*), intent(in), optional :: advance
    integer, intent(out) :: ios
    logical :: advancing

    integer(c_int) :: ioerr
    ios=0
    advancing=.true.
    if (present(advance)) then
      advancing=(advance == 'yes')
    end if
    ioerr=fputs(trim(lin) // C_NULL_CHAR, fd%handle)
    ios=ioerr
    if (advancing) then
      ioerr=fputs(char(10) // C_NULL_CHAR, fd%handle)
      ios=ioerr
    end if
  end subroutine pipe_write
!
! wrapper for pclose
  subroutine pipe_close(fd, ios)
    use, intrinsic :: iso_c_binding
    type(ioport) :: fd
    integer :: ios
    integer(c_int) :: ic 
    ios = 0
    ic = pclose(fd%handle)
    if (ic /= 0) ios = ic
  end subroutine pipe_close
end module f95pipes
#endif
#if !(ZLIB)
!
! time/random number generator seeds
!
module rndseed
  integer :: ix = 1, iy = 2, iz = 3
  integer :: initix = 1, initiy = 2, initiz = 3
end module rndseed
!
! random number generators
!
module rngs 
  contains
!  
! Algorithm AS 183 Appl Stat 1982; 31:188
! Returns a pseudo-random number from U(0,1)
!  
! ix,iy,iz should be "randomly" initialised to 1-30000
! eg via time
!  
  function random()
    use rndseed
    real :: random
    ix=171*mod(ix,177)-2*(ix/177)
    iy=172*mod(iy,176)-35*(iy/176)
    iz=170*mod(iz,178)-63*(iz/178)
    if (ix < 0) ix=ix+30269
    if (iy < 0) iy=iy+30307
    if (iz < 0) iz=iz+30323
    random=amod(float(ix)/30269.0+float(iy)/30307.0 + float(iz)/30323.0,1.0)
  end function random
!  
! Return a pseudo-random integer from integer U(lo..hi)
!  
  integer function irandom(lo,hi)
    integer, intent(in) :: lo
    integer, intent(in) :: hi
    real :: x
    x=random()
    irandom=lo+int(float(hi-lo+1)*x)
    if (irandom > hi) irandom=hi
  end function irandom
!
! Random character string
!
  subroutine uniqnam(nchar, str)
    integer, intent(in) :: nchar
    character (len=*) :: str
    integer i
    str=' '
    do i=1, nchar
      str(i:i)=char(96+irandom(1,26))
    end do  
  end subroutine uniqnam
end module rngs
#endif
!
! Output stream, formatting
!
module outstream
  integer :: outstr       ! stream
  integer :: logstr       ! logging stream
  character (len=1) :: tabsep = ' ' ! character to separate output words

contains
!  
! format free output, adding newline if line exceeding specified length
!
  subroutine newlin(sol, eol, pos, newpos)
    integer, intent(in) :: sol, eol
    integer, intent(inout) :: pos
    integer, intent(in) :: newpos
    
    if (pos > eol) then
      pos=newpos
      write(outstr,*)
      write(outstr,'(a)', advance='no') repeat(' ', sol-1)
    end if
    return
  end subroutine newlin
end module outstream
!
! Readline subroutine for either plain or gzipped files -- 
!
module fileio
  use iocodes
  use ioports
#if ZLIB
  use f95zlib
#endif
#if POPEN
  use f95pipes 
#endif
  use outstream
  public :: close_port, newlun, open_port, readline, rewind_port
contains
!
! Find a free Fortran style unit
!
  subroutine newlun(strm)
    integer, intent(out) :: strm
    integer, parameter :: MAXUNITS = 99
    integer :: iport
    logical :: ios
    do iport=8, MAXUNITS
      inquire(iport, opened=ios)
      if (.not.ios) then
        strm=iport
        return
      end if
    end do
    write(*,'(a)') 'No available i/o streams!'
    stop
  end subroutine newlun 
!
! Test if gzipped file, reading magic number 31,139
!
  function isgzipped(filnam)
    logical :: isgzipped
    character (len=*), intent(in) :: filnam
    integer :: s
    character (len=1) :: ch1, ch2
    isgzipped=.false.
    call newlun(s)

    open(s, file=filnam, access=stream_access, form=stream_form,  &
         status='old', iostat=ios)
    if (ios /= 0) then
      write(outstr, '(3a)') 'ERROR: Cannot open "', trim(filnam), '".'
      return
    end if
    read(s, iostat=ios) ch1, ch2
    if (ios /= 0) then
      write(outstr, '(3a)') 'ERROR: Could not read header of "', trim(filnam), '".'
    else
      isgzipped=(ichar(ch1) == 31 .and. ichar(ch2) == 139)
    end if
    close(s, status='keep')
    return
  end function isgzipped
!
! Open a (plain or gzipped) file or pipe for reading or writing
!
  subroutine open_port(filnam, port, mode, ios)
#if !ZLIB
    use rngs
#endif
    character (len=*), intent(in) :: filnam
    character (len=1), intent(in) :: mode
    type (ioport) :: port
    integer, intent(out) :: ios

    integer :: eon, strm
    logical :: gzipped, apipe
    character(len=3) :: fileage
#if !(ZLIB)
    character(len=len(filnam)) :: wrkfil
#endif
    ios=0
    apipe=.false.
    gzipped=.true. !ugly patch here
    fileage='old'
    if (mode == 'w') fileage='new'
    eon=len_trim(filnam)
    if (eon == 0) then
      write(outstr,'(a)') 'ERROR: No file name given.'
      ios=1
      return
    end if
    port%stat=mode
#if POPEN
    apipe=((mode == 'r' .and. filnam(eon:eon) == '|') .or.  &
           (mode == 'w' .and. filnam(1:1) == '|'))
#endif
!     if (.not.apipe) gzipped=isgzipped(filnam) !ugly patch here
    if (gzipped) then
#if ZLIB
      call fgz_open(filnam, mode // 'b', port, ios)
#else
      if (mode == 'r') then
        call uniqnam(5, wrkfil)
        wrkfil='sp-' // trim(wrkfil) // '.txt'
        call system('gzip -cd "' // trim(filnam) // '" > ' // wrkfil)
        port%filnam=wrkfil
      else
        if (filnam(max(1,(eon-2)):eon) == '.gz') then
          port%filnam=filnam(1:(eon-3))
        else
          port%filnam=filnam
        end if
      end if
      call newlun(strm)
      open(strm, file=port%filnam, status=fileage, iostat=ios)
      port%filtyp=PORT_COPY
      port%fstream=strm
#endif
#if POPEN
    else if (apipe) then
      if (mode == 'r') then
        call pipe_open(filnam(1:(eon-1)), port, ios)
      else
        call pipe_open(filnam(2:eon), port, ios)
      end if
#endif
    else
      call newlun(strm)
      open(strm, file=filnam, status=fileage, iostat=ios)
      port%filnam=filnam
      port%filtyp=PORT_STANDARD
      port%fstream=strm
    end if
  end subroutine open_port
!
! Reopen a file for reading or writing
!
  subroutine rewind_port(port, ios)
    type(ioport), intent(inout) :: port
    integer, intent(out) :: ios

    ios=0
    if (port%filtyp == PORT_STANDARD .or. port%filtyp == PORT_COPY) then
      rewind(port%fstream)
#if ZLIB
    else if (port%filtyp == PORT_GZIPPED) then
      call fgz_rewind(port, ios)
#endif
#if POPEN
    else if (port%filtyp == PORT_PIPE) then
      call pipe_rewind(port, ios)
#endif
    end if
  end subroutine rewind_port
!
! Read one record from file
!
  subroutine readline(port, lin, advance, ios)
    type (ioport), intent(in) :: port
    character(len=*) :: lin
    character(len=*), optional :: advance
    integer, intent(out) :: ios  
    character (len=3) :: advancing
    ios=0
    advancing='yes'
    if (present(advance)) then
      advancing=advance
    end if
    if (port%filtyp == PORT_STANDARD .or. port%filtyp == PORT_COPY) then
      read(port%fstream,'(a)', advance=advancing, iostat=ios) lin
#if ZLIB
    else if (port%filtyp == PORT_GZIPPED) then
      call fgz_read(port, lin, advance=advancing, ios=ios)
#endif
#if POPEN
    else if (port%filtyp == PORT_PIPE) then
      call pipe_read(port, lin, advance=advancing, ios=ios)
#endif
    end if
  end subroutine readline
!
! Write one record to file
!
  subroutine writeline(port, lin, advance, ios)
    type (ioport), intent(in) :: port
    character(len=*) :: lin
    character(len=*), optional :: advance
    integer, intent(out) :: ios  
    character (len=3) :: advancing
    ios=0
    advancing='yes'
    if (present(advance)) then
      advancing=advance
    end if
    if (port%filtyp == PORT_STANDARD .or. port%filtyp == PORT_COPY) then
      write(port%fstream,'(a)', advance=advancing, iostat=ios) lin
#if ZLIB
    else if (port%filtyp == PORT_GZIPPED) then
      call fgz_write(port, lin, advance=advancing, ios=ios)
#endif
#if POPEN
    else if (port%filtyp == PORT_PIPE) then
      call pipe_write(port, lin, advance=advancing, ios=ios)
#endif
    end if
  end subroutine writeline
!
! Close file for reading -  
!   if gzipped but ZLIB not available, then delete temporary file
!
  subroutine close_port(port, ios)
    type (ioport), intent(in) :: port
    integer, intent(out) :: ios
    ios=0
    if (port%filtyp == PORT_STANDARD) then
      close(port%fstream)
#if ZLIB
    else if (port%filtyp == PORT_GZIPPED) then
      call fgz_close(port, ios)
#else
    else if (port%filtyp == PORT_COPY) then
      close(port%fstream)
      if (port%stat == 'w') then
        call system('gzip  -f "' // trim(port%filnam) // '"')
      else
        call delfile(port%filnam, -2)
      end if
#endif
#if POPEN
    else if (port%filtyp == PORT_PIPE) then
      call pipe_close(port, ios)
#endif
    end if
  end subroutine close_port
!
! Delete a file
!
  subroutine delfile(filnam, plevel)
    use outstream
    character (len=*), intent(in) :: filnam
    integer, intent(in) :: plevel
    integer :: ioerr
#if IFORT || SUN
    integer :: unlink
    ioerr=unlink(filnam)
#elif OPEN64
    ioerr=0
    call unlink(filnam)
#else
    call unlink(filnam, ioerr)
#endif
    if (ioerr == 0 .and. plevel >= 0) then
      write(outstr, '(3a)') 'Deleted file "', trim(filnam),'".'
    end if
  end subroutine delfile
end module fileio
