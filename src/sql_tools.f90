module sql_tools
  use sqlite
  use kinds
  use report_module
  use settings, only: nummu, nstokes

  type(SQLITE_DATABASE)                      :: db
  type(SQLITE_STATEMENT)                     :: stmt, stmt2
  type(SQLITE_COLUMN), dimension(:,:), pointer :: writeCols
  type(SQLITE_COLUMN), dimension(:), pointer :: readColTmp, readCols
  contains

  subroutine sql_tools_open_table(fname,par_name)
    implicit none
    character(len=200), intent(in) :: fname
    character(len=20), intent(in) :: par_name

    logical:: db_exists

    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'sql_tools_open_table'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

    INQUIRE(FILE=TRIM(fname), EXIST=db_exists)
    call sqlite3_open( TRIM(fname), db )

  !   call sql_tools_create_readCols()


    if (db_exists .eqv. .false.) then

      call sql_tools_create_writeCols()

      call sqlite3_create_table( db,TRIM(par_name)//"_H" , writeCols(1,:),primary="id" )
      call sqlite3_create_table( db,TRIM(par_name)//"_V" , writeCols(2,:),primary="id" )

      deallocate(writeCols)


    end if

    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return

  end subroutine sql_tools_open_table


  subroutine sql_tools_close_table()
    implicit none

    integer(kind=long) :: errorstatus
    integer(kind=long) :: err = 0
    character(len=80) :: msg
    character(len=30) :: nameOfRoutine = 'sql_tools_close_table'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)


    call sqlite3_close( db )


    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
    return

  end subroutine sql_tools_close_table

  subroutine sql_tools_create_writeCols()

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
    character(len=30) :: nameOfRoutine = 'sql_tools_create_writeCols'

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
  end subroutine sql_tools_create_writeCols
! 
!   subroutine sql_tools_create_readCols()
! 
!     implicit none
! 
!     character(len=2) :: i1s, i2s, i3s, i4s, i5s
!     character(len=8) :: fmt
!     character(len=17) ::colNameSM
!     character(len=11) ::colNameEM
!     character(len=8) ::colNameEV
! 
!     integer :: i1, i2, i3, i4, i5
!     integer, dimension(2) :: lls
!     integer :: lenCols
! 
!     integer(kind=long) :: errorstatus
!     integer(kind=long) :: err = 0
!     character(len=80) :: msg
!     character(len=30) :: nameOfRoutine = 'sql_tools_create_readCols'
! 
!     if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)
! 
! 
!     lenCols=7 +((nstokes * nummu) + (nstokes**2 * nummu) + (nstokes**2 *nummu* nummu * 2))/2
! 
!     allocate(readCols(lenCols))
! 
! 
! 	call sqlite3_column_props( readCols(1), "id", SQLITE_INT )
! 	call sqlite3_column_props( readCols(2), "frequency", SQLITE_CHAR, 16 )
! 	call sqlite3_column_props( readCols(3), "temperature", SQLITE_CHAR, 16 )
! 	call sqlite3_column_props( readCols(4), "diameter", SQLITE_CHAR, 16 )
! 	call sqlite3_column_props( readCols(5), "mass", SQLITE_CHAR, 16 )
! 	call sqlite3_column_props( readCols(6), "as_ratio", SQLITE_CHAR, 16 )
! 	call sqlite3_column_props( readCols(7), "type", SQLITE_CHAR, 4 )
! 
!       fmt = '(I2.2)'
!       lls = (/8,8/)
!     ! 
!       i1=1
!       write (i1s,fmt) i1
! 	do i2=1, nummu
! 	  write (i2s,fmt) i2
! 	  do i3=1, nstokes
! 	    write (i3s,fmt) i3
! 	    do i4=1, nummu
! 	      write (i4s,fmt) i4
! 	      do i5=1, 2
! 		write (i5s,fmt) i5
! 		colNameSM = "SM_"//i1s//"_"//i2s//"_"//i3s//"_"//i4s//"_"//i5s
! 		call sqlite3_column_props( readCols(lls(i1)), colNameSM, SQLITE_REAL )
! 		lls(i1) = lls(i1)+1
! 	      end do
! 	    end do
! 	  end do
! 	end do    
! 
! 	do i2=1, nstokes
! 	  write (i2s,fmt) i2
! 	  do i3=1, nummu
! 	    write (i3s,fmt) i3
! 	    colNameEM = "EM_"//i1s//"_"//i2s//"_"//i3s
! 	    call sqlite3_column_props( readCols(lls(i1)), colNameEM, SQLITE_REAL )
! 		lls(i1) = lls(i1)+1
! 	  end do
! 	end do 
! 
! 	do i2=1, nummu
! 	  write (i2s,fmt) i2
! 	  colNameEV = "EV_"//i1s//"_"//i2s
! 	  call sqlite3_column_props( readCols(lls(i1)), colNameEV, SQLITE_REAL )
! 	  lls(i1) = lls(i1)+1
! 	end do 
! 
! 
!     if (verbose >= 2) call report(info,'End of ', nameOfRoutine)
!     return
! 
!   end subroutine sql_tools_create_readCols

  subroutine sql_tools_get_entry(par_name,freq,t,as_ratio,diameter,particle_mass,ptype,&
      scatter_matrix,extinct_matrix,emis_vector,found)

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
    character(len=30) :: nameOfRoutine = 'sql_tools_get_entry'

    if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)



    par_nameHV ="HV"
    scatter_matrix(:,:,:,:,:) = 0.d0
    extinct_matrix(:,:,:) = 0.d0
    emis_vector(:,:) = 0.d0

!     fmt="(ES16.10)"
! 
!     write(freq_s,fmt), freq
!     write(t_s,fmt), t
!     write(as_ratio_s,fmt), as_ratio
!     write(diameter_s,fmt), diameter
!     write(particle_mass_s,fmt), particle_mass
! 
!     if (verbose >= -4) print*, freq_s,t_s,as_ratio_s ,diameter_s , particle_mass_s

    
    fmt="(ES12.6)"

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

    call sqlite3_prepare( db, "SELECT * FROM "//trim(par_name)//'_'//par_nameHV(i1:i1)//' '//trim(whereString), stmt, readCols )


  !   call sqlite3_prepare_select( db, trim(par_name)//'_H', readCols, stmt,whereString)

    call sqlite3_next_row( db,stmt, readCols, finished )

    if (finished .eqv. .true.) then 

      if (verbose >= 4) print*, "NOTHING FOUND!"
      found=.false.
      deallocate(readCols)
      return
      
    else
    if (verbose >= 4) print*, "found entry!"

    if (i1 == 1) call sqlite3_get_column( readCols(1), jj )
	do i2=1, nummu
	  write (i2s,fmt) i2
	  do i3=1, nstokes
	    write (i3s,fmt) i3
	    do i4=1, nummu
	      write (i4s,fmt) i4
	      do i5=1, 2
		write (i5s,fmt) i5
		call sqlite3_get_column( readCols(lls(i1)), realresult )
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
	    call sqlite3_get_column( readCols(lls(i1)), realresult )
	    extinct_matrix(i1,i2,i3) = realresult
	    lls(i1) = lls(i1)+1
	  end do
	end do 

	do i2=1, nummu
	  write (i2s,fmt) i2
	  call sqlite3_get_column( readCols(lls(i1)), realresult )
	  emis_vector(i1,i2) = realresult
	  lls(i1) = lls(i1)+1
	end do 

      if (verbose >= 4) print*, "read", lls(i1) -1, "of", size(readCols)

      deallocate(readCols)
      found=.true.
    end if
    


    end do





  !   deallocate(readCols,readColsV)

    if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

    return

  end subroutine sql_tools_get_entry

  subroutine sql_tools_write_entry(par_name,freq,t,as_ratio,diameter,&
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
    character(len=30) :: nameOfRoutine = 'sql_tools_write_entry'

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


    call sql_tools_create_writeCols()
    

    fmt="(ES12.6)"

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

  end subroutine sql_tools_write_entry

end module sql_tools