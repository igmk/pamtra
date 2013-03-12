module file_mod

    implicit none
    save

    character(99)  :: input_file        ! name of profile
    character(300) :: namelist_file     ! name of nml_file
    character(300) :: nc_out_file       ! name of netcdf output file

end module file_mod
