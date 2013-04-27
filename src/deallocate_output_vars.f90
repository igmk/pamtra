subroutine deallocate_output_vars

    use vars_atmosphere
    use vars_output
    use settings
    use mod_io_strings

    implicit none

    if (write_nc)  deallocate(is,js,lons,lats,lfracs,iwvs,cwps,iwps,rwps,swps,gwps,hwps)
    if (write_nc .or. in_python) deallocate(tb)

    if (active) then
        deallocate(Ze, radar_hgt, Att_hydro, Att_atmo)
    end if

    if ((active) .and. ((radar_mode .eq. "simple") .or. (radar_mode .eq. "splitted"))) then
        deallocate(Ze_cw,Ze_rr,Ze_ci,Ze_sn,Ze_gr,Ze_ha, &
        Att_cw,Att_rr,Att_ci,Att_sn,Att_gr,Att_ha, &
        )
    end if
    if((active) .and. ((radar_mode .eq. "spectrum") .or. (radar_mode .eq. "moments"))) then
        deallocate(radar_spectra, radar_snr, radar_vel,radar_moments, radar_slope, radar_quality)
    end if
end subroutine deallocate_output_vars
