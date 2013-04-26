module sql_tools
   use sqlite
   use kinds
   use report_module
   use settings, only: nummu, nstokes

type(SQLITE_DATABASE)                      :: db
type(SQLITE_STATEMENT)                     :: stmt, stmt2
type(SQLITE_COLUMN), dimension(:,:), pointer :: writeCols
type(SQLITE_COLUMN), dimension(:), pointer :: readColTmp, readColsH, readColsV
contains


subroutine sql_open_table(fname,par_name)
  implicit none
  character(len=200), intent(in) :: fname
  character(len=20), intent(in) :: par_name

  logical:: db_exists

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=30) :: nameOfRoutine = 'sql_open_table'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  INQUIRE(FILE=TRIM(fname), EXIST=db_exists)
  call sqlite3_open( TRIM(fname), db )

  call sql_create_readCols()


  if (db_exists .eqv. .false.) then

    call sql_create_writeCols()

    call sqlite3_create_table( db,TRIM(par_name)//"_H" , writeCols(1,:),primary="id" )
    call sqlite3_create_table( db,TRIM(par_name)//"_V" , writeCols(2,:),primary="id" )

    deallocate(writeCols)


  end if

  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return

end subroutine sql_open_table


subroutine sql_close_table()
  implicit none

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=30) :: nameOfRoutine = 'sql_close_table'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)


  call sqlite3_close( db )


  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return

end subroutine sql_close_table

subroutine sql_create_writeCols()

   implicit none

  character(len=2) :: i1s, i2s, i3s, i4s, i5s
  character(len=8) :: fmt
  character(len=17) ::colNameSM
  character(len=11) ::colNameEM
  character(len=8) ::colNameEV

  integer :: i1, i2, i3, i4, i5
  integer, dimension(2) :: lls
  integer :: lenCols

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=30) :: nameOfRoutine = 'sql_create_writeCols'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)


  lenCols=7 +((nstokes * nummu) + (nstokes**2 * nummu) + (nstokes**2 *nummu**2 * 2))/2

  allocate(writeCols(nstokes,lenCols))


    do i1=1, nstokes
      call sqlite3_column_props( writeCols(i1,1), "id", SQLITE_INT )
      call sqlite3_column_props( writeCols(i1,2), "frequency", SQLITE_CHAR, 16  )
      call sqlite3_column_props( writeCols(i1,3), "temperature", SQLITE_CHAR, 16  )
      call sqlite3_column_props( writeCols(i1,4), "diameter", SQLITE_CHAR, 16  )
      call sqlite3_column_props( writeCols(i1,5), "mass", SQLITE_CHAR, 16 )
      call sqlite3_column_props( writeCols(i1,6), "as_ratio", SQLITE_CHAR, 16  )
      call sqlite3_column_props( writeCols(i1,7), "type", SQLITE_CHAR, 4 )

    end do

    fmt = '(I2.2)'
    lls = (/8,8/)
  ! 
    do i1=1, nstokes
      write (i1s,fmt) i1
      do i2=1, nummu
	write (i2s,fmt) i2
	do i3=1, nstokes
	  write (i3s,fmt) i3
	  do i4=1, nummu
	    write (i4s,fmt) i4
	    do i5=1, 2
	      write (i5s,fmt) i5
	      colNameSM = "SM_"//i1s//"_"//i2s//"_"//i3s//"_"//i4s//"_"//i5s
	      call sqlite3_column_props( writeCols(i1,lls(i1)), colNameSM, SQLITE_REAL )
	      lls(i1) = lls(i1)+1
	    end do
	  end do
	end do
      end do    
    end do

    do i1=1, nstokes
      write (i1s,fmt) i1
      do i2=1, nstokes
	write (i2s,fmt) i2
	do i3=1, nummu
	  write (i3s,fmt) i3
	  colNameEM = "EM_"//i1s//"_"//i2s//"_"//i3s
	  call sqlite3_column_props( writeCols(i1,lls(i1)), colNameEM, SQLITE_REAL )
	      lls(i1) = lls(i1)+1
	end do
      end do 
    end do

    do i1=1, nstokes
      write (i1s,fmt) i1
      do i2=1, nummu
	write (i2s,fmt) i2
	colNameEV = "EV_"//i1s//"_"//i2s
	call sqlite3_column_props( writeCols(i1,lls(i1)), colNameEV, SQLITE_REAL )
	lls(i1) = lls(i1)+1
      end do 
    end do

  if (verbose >= 4) print*, "created write_cols, len=*", lls(1)-1, lls(2)-1

  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return
end subroutine sql_create_writeCols

subroutine sql_create_readCols()

   implicit none

  character(len=2) :: i1s, i2s, i3s, i4s, i5s
  character(len=8) :: fmt
  character(len=17) ::colNameSM
  character(len=11) ::colNameEM
  character(len=8) ::colNameEV

  integer :: i1, i2, i3, i4, i5
  integer, dimension(2) :: lls
  integer :: lenCols

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=30) :: nameOfRoutine = 'sql_create_readCols'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)


  lenCols=7 +((nstokes * nummu) + (nstokes**2 * nummu) + (nstokes**2 *nummu* nummu * 2))/2

  allocate(readColsH(lenCols))
  allocate(readColsV(lenCols))


      call sqlite3_column_props( readColsH(1), "id", SQLITE_INT )
      call sqlite3_column_props( readColsH(2), "frequency", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsH(3), "temperature", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsH(4), "diameter", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsH(5), "mass", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsH(6), "as_ratio", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsH(7), "type", SQLITE_CHAR, 4 )

    fmt = '(I2.2)'
    lls = (/8,8/)
  ! 
    i1=1
    write (i1s,fmt) i1
      do i2=1, nummu
	write (i2s,fmt) i2
	do i3=1, nstokes
	  write (i3s,fmt) i3
	  do i4=1, nummu
	    write (i4s,fmt) i4
	    do i5=1, 2
	      write (i5s,fmt) i5
	      colNameSM = "SM_"//i1s//"_"//i2s//"_"//i3s//"_"//i4s//"_"//i5s
	      call sqlite3_column_props( readColsH(lls(i1)), colNameSM, SQLITE_REAL )
	      lls(i1) = lls(i1)+1
	    end do
	  end do
	end do
      end do    

      do i2=1, nstokes
	write (i2s,fmt) i2
	do i3=1, nummu
	  write (i3s,fmt) i3
	  colNameEM = "EM_"//i1s//"_"//i2s//"_"//i3s
	  call sqlite3_column_props( readColsH(lls(i1)), colNameEM, SQLITE_REAL )
	      lls(i1) = lls(i1)+1
	end do
      end do 

      do i2=1, nummu
	write (i2s,fmt) i2
	colNameEV = "EV_"//i1s//"_"//i2s
	call sqlite3_column_props( readColsH(lls(i1)), colNameEV, SQLITE_REAL )
	lls(i1) = lls(i1)+1
      end do 

      call sqlite3_column_props( readColsV(1), "id", SQLITE_INT )
      call sqlite3_column_props( readColsV(2), "frequency", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsV(3), "temperature", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsV(4), "diameter", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsV(5), "mass", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsV(6), "as_ratio", SQLITE_CHAR, 16 )
      call sqlite3_column_props( readColsV(7), "type", SQLITE_CHAR, 4 )


    fmt = '(I2.2)'
    lls = (/8,8/)
  ! 

    i1=2
    write (i1s,fmt) i1
      do i2=1, nummu
	write (i2s,fmt) i2
	do i3=1, nstokes
	  write (i3s,fmt) i3
	  do i4=1, nummu
	    write (i4s,fmt) i4
	    do i5=1, 2
	      write (i5s,fmt) i5
	      colNameSM = "SM_"//i1s//"_"//i2s//"_"//i3s//"_"//i4s//"_"//i5s
	      call sqlite3_column_props( readColsV(lls(i1)), colNameSM, SQLITE_REAL )
	      lls(i1) = lls(i1)+1
	    end do
	  end do
	end do
      end do    

      do i2=1, nstokes
	write (i2s,fmt) i2
	do i3=1, nummu
	  write (i3s,fmt) i3
	  colNameEM = "EM_"//i1s//"_"//i2s//"_"//i3s
	  call sqlite3_column_props( readColsV(lls(i1)), colNameEM, SQLITE_REAL )
	      lls(i1) = lls(i1)+1
	end do
      end do 

      do i2=1, nummu
	write (i2s,fmt) i2
	colNameEV = "EV_"//i1s//"_"//i2s
	call sqlite3_column_props( readColsV(lls(i1)), colNameEV, SQLITE_REAL )
	lls(i1) = lls(i1)+1
      end do 


  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return

end subroutine sql_create_readCols

subroutine sql_get_entry(par_name,freq,t,as_ratio,diameter,particle_mass,ptype,scatter_matrix,extinct_matrix,emis_vector,found)

  implicit none

  character(len=20), intent(in) :: par_name
  real(kind=dbl), intent(in) :: freq !in Hz
  real(kind=dbl), intent(in) :: t
  real(kind=dbl), intent(in) :: as_ratio
  real(kind=dbl), intent(in) :: particle_mass
  real(kind=dbl), intent(in) :: diameter
  character(len=4), intent(in) :: ptype

  character(len=23) :: freq_s !in Hz
  character(len=23) :: t_s
  character(len=23) :: as_ratio_s
  character(len=23) :: particle_mass_s
  character(len=23) :: diameter_s

  character(len=10) :: fmt
  character(len=300) :: whereString
  integer :: ll, jj
  integer :: i1, i2, i3, i4, i5
  character(len=2) :: i1s, i2s, i3s, i4s, i5s
  character(len=8) :: jjs
  integer,dimension(2) :: lls

  real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector
  logical, intent(out) :: found
  real(kind=dbl) :: realresult
  logical :: finished
  character(len=22) :: par_nameHV

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=30) :: nameOfRoutine = 'sql_get_entry'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  call sql_create_readCols()


  par_nameHV ="HV"
  scatter_matrix(:,:,:,:,:) = 0.d0
  extinct_matrix(:,:,:) = 0.d0
  emis_vector(:,:) = 0.d0

  fmt="(ES16.10)"

  write(freq_s,fmt), freq
  write(t_s,fmt), t
  write(as_ratio_s,fmt), as_ratio
  write(diameter_s,fmt), diameter
  write(particle_mass_s,fmt), particle_mass

  if (verbose >= 4) print*, freq_s,t_s,as_ratio_s ,diameter_s , particle_mass_s


    fmt = '(I2.2)'
    lls = (/8,8/)

  do i1=1, nstokes
    write (i1s,fmt) i1

  if (i1 == 1) then
  whereString = ' WHERE (frequency="'//trim(freq_s)//'"'//&
    ' and temperature="'//trim(t_s)//'"'//&
    ' and diameter="'//trim(diameter_s)//'"'//&
    ' and mass="'//trim(particle_mass_s)//'"'//&
    ' and as_ratio="'//trim(as_ratio_s)//'"'//&
    ' and type="'//TRIM(ptype)//'" )'
  else
    write (jjs,'(I7.7)') jj
    whereString = ' WHERE id='//jjs
  end if


  if (verbose >= 4) print*, "SELECT * FROM "//trim(par_name)//'_'//par_nameHV(i1:i1)//' '//trim(whereString)
  call sqlite3_reset( stmt )

  call sqlite3_prepare( db, "SELECT * FROM "//trim(par_name)//'_'//par_nameHV(i1:i1)//' '//trim(whereString), stmt, readColsH )


!   call sqlite3_prepare_select( db, trim(par_name)//'_H', readColsH, stmt,whereString)

  call sqlite3_next_row( db,stmt, readColsH, finished )

  if (finished .eqv. .false.) then 


  if (verbose >= 4) print*, "found entry!"

  if (i1 == 1) call sqlite3_get_column( readColsH(1), jj )

      do i2=1, nummu
	write (i2s,fmt) i2
	do i3=1, nstokes
	  write (i3s,fmt) i3
	  do i4=1, nummu
	    write (i4s,fmt) i4
	    do i5=1, 2
	      write (i5s,fmt) i5
	      call sqlite3_get_column( readColsH(lls(i1)), realresult )
	      scatter_matrix(i1,i2,i3,i4,i5) = realresult
	      lls(i1) = lls(i1)+1
	    end do
	  end do
	end do
      end do    


      do i2=1, nstokes
	write (i2s,fmt) i2
	do i3=1, nummu
	  write (i3s,fmt) i3
	  call sqlite3_get_column( readColsH(lls(i1)), realresult )
	  extinct_matrix(i1,i2,i3) = realresult
	  lls(i1) = lls(i1)+1
	end do
      end do 

      do i2=1, nummu
	write (i2s,fmt) i2
	call sqlite3_get_column( readColsH(lls(i1)), realresult )
	emis_vector(i1,i2) = realresult
	lls(i1) = lls(i1)+1
      end do 

  if (verbose >= 4) print*, "read", lls(i1) -1, "of", size(readColsH)

    deallocate(readColsH)
    found=.true.
  else
  if (verbose >= 4) print*, "NOTHING FOUND!"
    found=.false.
  end if
   


  end do





!   deallocate(readColsH,readColsV)

  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

  return

end subroutine sql_get_entry

subroutine sql_write_entry(par_name,freq,t,as_ratio,diameter,&
    particle_mass,ptype,&
    scatter_matrix,extinct_matrix,emis_vector)

  implicit none

  character(len=20), intent(in) :: par_name
  character(len=22) :: par_nameHV

  real(kind=dbl), intent(in) :: freq
  real(kind=dbl), intent(in) :: t
  real(kind=dbl), intent(in) :: as_ratio
  real(kind=dbl), intent(in) :: particle_mass
  real(kind=dbl), intent(in) :: diameter
  character(len=4), intent(in) :: ptype

  real(kind=dbl), intent(out), dimension(nstokes,nummu,nstokes,nummu,2) :: scatter_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nstokes,nummu) :: extinct_matrix
  real(kind=dbl), intent(out), dimension(nstokes,nummu) :: emis_vector

  logical:: finished

  integer :: jj
  integer :: i1, i2, i3, i4, i5
  integer, dimension(2) :: lls, llsOffset

  character(len=10) ::jjs
  character(len=13) ::whereString
  character(len=2) :: polChar

  character(len=23) :: freq_s !in Hz
  character(len=23) :: t_s
  character(len=23) :: as_ratio_s
  character(len=23) :: particle_mass_s
  character(len=23) :: diameter_s
  character(len=10) :: fmt

  integer(kind=long) :: errorstatus
  integer(kind=long) :: err = 0
  character(len=80) :: msg
  character(len=30) :: nameOfRoutine = 'sql_write_entry'

  if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

  !first get the index we want to write to

  jj=0
  polChar="HV"

  allocate(readColTmp(1))
  call sqlite3_column_query( readColTmp(1), "id", SQLITE_INT )
  call sqlite3_reset( stmt )
  call sqlite3_prepare_select( db, trim(par_name)//"_V", readColTmp, stmt,'order by id desc limit 1')

  call sqlite3_next_row(db,  stmt, readColTmp, finished )

  if (finished .eqv. .true.) then
    if (verbose >= 4) call report(info,'table empty', nameOfRoutine)
    jj = 1
  else
    call sqlite3_get_column( readColTmp(1), jj )
    if (verbose >= 4) print*, 'found entries', jj
    jj = jj+1
  end if
  deallocate(readColTmp)

! call sqlite3_reset( stmt )

  write (jjs,'(I7.7)') jj
  whereString = "id="//jjs

  if (verbose >= 4) print*, 'writing index no', jj, jjs


  call sql_create_writeCols()
  

  fmt="(ES16.10)"

  write(freq_s,fmt), freq
  write(t_s,fmt), t
  write(as_ratio_s,fmt), as_ratio
  write(diameter_s,fmt), diameter
  write(particle_mass_s,fmt), particle_mass

  lls = (/1,1/)
  llsOffset = (/7,7/)

  !write the index and infos
  do i1=1, nstokes
  call sqlite3_begin( db )


  if (verbose >= 4) print*, "JJ", jj
    call sqlite3_set_column( writeCols(i1,1), jj )
    call sqlite3_set_column( writeCols(i1,2), trim(freq_s) )
    call sqlite3_set_column( writeCols(i1,3), trim(t_s ))
    call sqlite3_set_column( writeCols(i1,4), trim(diameter_s ))
    call sqlite3_set_column( writeCols(i1,5), trim(particle_mass_s ))
    call sqlite3_set_column( writeCols(i1,6), trim(as_ratio_s) )
    call sqlite3_set_column( writeCols(i1,7), trim(ptype) )

    par_nameHV = TRIM(par_name)//"_"//polChar(i1:i1)
    call sqlite3_insert( db, TRIM(par_nameHV), writeCols(i1,1:7) )


      do i2=1, nummu
	do i3=1, nstokes
	  do i4=1, nummu
	    do i5=1, 2
	      call sqlite3_set_column( writeCols(i1,lls(i1)+llsOffset(i1)),&
		scatter_matrix(i1,i2,i3,i4,i5) )
	        lls(i1) = lls(i1)+1
	    end do
	  end do
	end do
      !we cannot change more than 999 variables at once...
      if (lls(i1) >= 500) then
	  par_nameHV = TRIM(par_name)//"_"//polChar(i1:i1)
	  call sqlite3_update( db,TRIM(par_nameHV) ,&
	    writeCols(i1,llsOffset(i1)+1:llsOffset(i1)+lls(i1)-1), TRIM(whereString) )

	if (verbose >= 4) print*, "UPDATING ", par_nameHV, i1, lls(i1), llsOffset(i1)+1, llsOffset(i1)+lls(i1)-1
	llsOffset(i1) = llsOffset(i1) + lls(i1)-1
        lls(i1) = 1

      end if

      end do    

      do i2=1, nstokes
	do i3=1, nummu
	  call sqlite3_set_column( writeCols(i1,lls(i1)+llsOffset(i1)),&
	    extinct_matrix(i1,i2,i3) )
	      lls(i1) = lls(i1)+1
	end do
      end do 

      do i2=1, nummu
	call sqlite3_set_column( writeCols(i1,lls(i1)+llsOffset(i1)),&
	   emis_vector(i1,i2) )
	lls(i1) = lls(i1)+1
      end do 

    par_nameHV = TRIM(par_name)//"_"//polChar(i1:i1)
    if (verbose >= 4) print*, "UPDATING ", par_nameHV, i1, lls(i1), llsOffset(i1)+1, llsOffset(i1)+lls(i1)-1
    call sqlite3_update( db, TRIM(par_nameHV),&
      writeCols(i1,llsOffset(i1)+1:llsOffset(i1)+lls(i1)-1), TRIM(whereString) )
  call sqlite3_commit( db )  

    end do    


  deallocate(writeCols)
  if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
  return

end subroutine sql_write_entry


subroutine sql_open_table2(fname,par_name)

   implicit none

character(len=200), intent(in) :: fname
character(len=20), intent(in) :: par_name


character(len=2) :: i1s, i2s, i3s, i4s, i5s
character(len=8) :: fmt
character(len=17) ::colNameSM
character(len=11) ::colNameEM
character(len=8) ::colNameEV
character(len=10) ::jjs
character(len=13) ::whereString
integer :: i1, i2, i3, i4, i5, ll, jj
integer, dimension(2) :: lls
integer :: lenCols, completion

integer, parameter ::  sqlMaxValues=998

type(SQLITE_DATABASE)                      :: db
type(SQLITE_STATEMENT)                     :: stmt, stmt2
type(SQLITE_COLUMN), dimension(:,:), pointer :: columns
type(SQLITE_COLUMN), dimension(:), pointer :: res1, results

LOGICAL :: db_exists, finished
real(kind=dbl) :: rand, realresult

lenCols=7 +((nstokes * nummu) + (nstokes**2 * nummu) + (nstokes**2 *nummu* nummu * 2))/2

!+((nstokes**2 *nummu**2 * 2))/2 ! + ((nstokes**2 * nummu) + (nstokes * nummu) + (nstokes**2 *nummu**2 * 2))/2
allocate(columns(nstokes,lenCols))
allocate(results(lenCols))

INQUIRE(FILE=TRIM(fname), EXIST=db_exists)
call sqlite3_open( TRIM(fname), db )

  if (verbose >= 4) print*, 1, TRIM(fname), "db_exists", db_exists, lenCols



  do i1=1, nstokes
  call sqlite3_column_props( columns(i1,1), "id", SQLITE_INT )
  call sqlite3_column_props( columns(i1,2), "frequency", SQLITE_REAL )
  call sqlite3_column_props( columns(i1,3), "temperature", SQLITE_REAL )
  call sqlite3_column_props( columns(i1,4), "diameter", SQLITE_REAL )
  call sqlite3_column_props( columns(i1,5), "mass", SQLITE_REAL )
  call sqlite3_column_props( columns(i1,6), "as_ratio", SQLITE_REAL )
  call sqlite3_column_props( columns(i1,7), "type", SQLITE_CHAR, 4 )

  call sqlite3_column_props( results(1), "id", SQLITE_INT )
  call sqlite3_column_props( results(2), "frequency", SQLITE_REAL )
  call sqlite3_column_props( results(3), "temperature", SQLITE_REAL )
  call sqlite3_column_props( results(4), "diameter", SQLITE_REAL )
  call sqlite3_column_props( results(5), "mass", SQLITE_REAL )
  call sqlite3_column_props( results(6), "as_ratio", SQLITE_REAL )
  call sqlite3_column_props( results(7), "type", SQLITE_CHAR, 4 )

end do
! call sqlite3_create_table( db,TRIM(par_name) , columns(1:6) )

  fmt = '(I2.2)'
  lls = (/8,8/)
! 
  do i1=1, nstokes
    write (i1s,fmt) i1
    do i2=1, nummu
      write (i2s,fmt) i2
      do i3=1, nstokes
	write (i3s,fmt) i3
	do i4=1, nummu
	  write (i4s,fmt) i4
	  do i5=1, 2
	    write (i5s,fmt) i5
	    colNameSM = "SM_"//i1s//"_"//i2s//"_"//i3s//"_"//i4s//"_"//i5s
  if (verbose >= 4) print*, colNameSM
	    call sqlite3_column_props( columns(i1,lls(i1)), colNameSM, SQLITE_REAL )
	    call sqlite3_column_props( results(lls(i1)), colNameSM, SQLITE_REAL )

!   call sqlite3_add_column( db,TRIM(par_name) , columns(ll) )

	    lls(i1) = lls(i1)+1
	  end do
	end do
      end do
    end do    
  end do
  do i1=1, nstokes
    write (i1s,fmt) i1
    do i2=1, nstokes
      write (i2s,fmt) i2
      do i3=1, nummu
	write (i3s,fmt) i3
	colNameEM = "EM_"//i1s//"_"//i2s//"_"//i3s
	call sqlite3_column_props( columns(i1,lls(i1)), colNameEM, SQLITE_REAL )
	call sqlite3_column_props( results(lls(i1)), colNameEM, SQLITE_REAL )

!   call sqlite3_add_column( db,TRIM(par_name) , columns(ll) )
	    lls(i1) = lls(i1)+1
      end do
    end do 
  end do
  do i1=1, nstokes
    write (i1s,fmt) i1
    do i2=1, nummu
      write (i2s,fmt) i2
      colNameEV = "EV_"//i1s//"_"//i2s
      call sqlite3_column_props( columns(i1,lls(i1)), colNameEV, SQLITE_REAL )
      call sqlite3_column_props( results(lls(i1)), colNameEV, SQLITE_REAL )

!   call sqlite3_add_column( db,TRIM(par_name) , columns(ll) )

	    lls(i1) = lls(i1)+1
    end do 
  end do

if (db_exists .eqv. .false.) then


call sqlite3_create_table( db,TRIM(par_name)//"_H" , columns(1,:),primary="id" )
call sqlite3_create_table( db,TRIM(par_name)//"_V" , columns(2,:),primary="id" )


print*, 5
end if


print*, lenCols


! call sqlite3_reset( stmt )
print*, "Calling"

allocate(res1(1))
call sqlite3_column_query( res1(1), "id", SQLITE_INT )

call sqlite3_prepare_select( db, 'test_H', res1, stmt,'order by id desc limit 1')
print*, "h"

call sqlite3_next_row(db,  stmt, res1, finished )
call sqlite3_get_column( res1(1), jj )


print*, "RES"

print*, jj
! call sqlite3_reset( stmt )


jj = jj+1
      write (jjs,'(I7.7)') jj
      whereString = "id="//jjs

      call sqlite3_set_column( columns(1,1), jj )
      call sqlite3_set_column( columns(2,1), jj )
   call sqlite3_begin( db )

    call sqlite3_insert( db, TRIM(par_name)//"_H", columns(1,1:1) )
    call sqlite3_insert( db, TRIM(par_name)//"_V", columns(2,1:1) )
! print*, ll, columns(1,1)%name, columns(1,1)%int_value
! print*, ll, columns(2,1)%name, columns(2,1)%int_value



print*, jjs, whereString

  do i1=1, nstokes

      do ll=2, 6
      call random(1, .false., rand)
rand=2.d0
      call sqlite3_set_column( columns(i1,ll), rand )
      end do
      call sqlite3_set_column( columns(i1,7), "snow" )

      do ll=8, 900
      call sqlite3_set_column( columns(i1,ll), rand )
! print*, ll, columns(i1,ll)%name, columns(i1,ll)%double_value
      end do

  if (i1 == 1) then
    call sqlite3_update( db, TRIM(par_name)//"_H", columns(1,1:900), TRIM(whereString) )
  else
    call sqlite3_update( db, TRIM(par_name)//"_V", columns(2,1:900), TRIM(whereString) )
  end if  

      do ll=901, lenCols
      call sqlite3_set_column( columns(i1,ll), rand )
! print*, ll, columns(i1,ll)%name, columns(i1,ll)%double_value
      end do

  if (i1 == 1) then
    call sqlite3_update( db, TRIM(par_name)//"_H", columns(1,901:lenCols), TRIM(whereString) )
  else
    call sqlite3_update( db, TRIM(par_name)//"_V", columns(2,901:lenCols), TRIM(whereString) )
  end if  
  end do

    call sqlite3_commit( db )  


print*, "Calling"

call sqlite3_reset( stmt )


call sqlite3_prepare_select( db, 'test_V', results, stmt,' WHERE (id=1 and temperature=2. and type="snow")')
print*, "h"

call sqlite3_next_row( db,stmt, results, finished )

if (finished .eqv. .false.) then 
  call sqlite3_get_column( results(10), realresult )
print*, realresult
else
  print*, "NOTHING FOUND!"
end if


call sqlite3_close( db )
print*, 6
deallocate(columns)
deallocate(results)
end subroutine

end module sql_tools