subroutine error_msg(infile,i,j)

    implicit none

	integer, intent(in),optional :: i,j

	character(25) :: hostname

	character(100), intent(in) :: infile

    call hostnm(hostname)

	print*, "Error in file action occured on ", trim(infile), " at ", i,j
	print*, "hostname: ", hostname

	return
end subroutine error_msg
