subroutine mass_size_rel(phase,ms_re_flag, p_size,particle_mass,as_ratio)

    use kinds
    use constants, only: pi

    implicit none

    real(kind=dbl) :: p_size, particle_size, particle_mass
    real(kind=dbl) :: as_ratio
    integer :: ms_re_flag
    character*1 phase

    ! input variables:
    !  p_size: in meter, m
    !  ms_re_flag: flag for mass size relationship,
    !  as_ratio: aspect ratio of the particles
    ! output variable:
    !  particle_mass: in gram
    !
    ! flag for mass size relationship
    ! ms_re_flag    1: matrosov, 2007, jas
    !               2: fixed density, 0.1 g/cm3
    !               3: fixed density, 0.2 g/cm
    !               4: brown and francis, 1995, j. atmos. oceanic technol.
    !               5: heymsfield et al., 2010, j. atmos. sci.
    ! internal variables
    !     particle_size, in micron-meter, um
    ! in mass-size relationship equations, for matrosov model, only the maximum dimension
    ! is needed, the aspect ratio is used to calculate particle with fixed densities
    ! or liquid water droplets

    particle_size = p_size*1.0d6
    if ((phase(1:1) .eq. 's') .or. (phase(1:1) .eq. 'S')) then
        if (ms_re_flag .eq. 1) then
            if ((particle_size .ge. 0.005d4) .and. (particle_size .le. 0.2d4)) then ! 50um~200um
                particle_mass = 0.003d0*(particle_size/1.0d4)**2.0d0
            else if ((particle_size .gt. 0.2d4) .and. (particle_size .le. 2.0d4)) then ! 200um~2cm
                particle_mass = 0.0067d0*(particle_size/1.0d4)**2.5d0
            else if (particle_size .gt. 2.0d4) then ! >2cm
                particle_mass = 0.0047d0*(particle_size/1.0d4)**3.0d0
            endif
        else if (ms_re_flag .eq. 2) then
            particle_mass=0.1d0*(particle_size/1d4)**3d0*as_ratio*pi/6d0
        else if (ms_re_flag .eq. 3) then
            particle_mass=0.2d0*(particle_size/1d4)**3d0*as_ratio*pi/6d0
        else if (ms_re_flag .eq. 4) then
            if (particle_size .lt. 66d0) then !<66um
                particle_mass = 480d0*(particle_size/1d6)**3d0*1d3
            else if (particle_size .ge. 66d0)then !66um
                particle_mass = 0.0121d0*(particle_size/1d6)**1.9d0*1d3
            end if
        else if (ms_re_flag .eq. 5) then
            if (particle_size .lt. 67d0) then ! <67um
                particle_mass = 480d0*(particle_size/1d6)**3d0*1d3
            else if (particle_size .ge. 67d0) then ! >67um
                particle_mass = 0.0837d0*(particle_size/1d6)**2.1d0*1d3
            end if
        end if
    !     if the phase is liquid, use the water density instead of ice density
    !     for water dropltets, mass size relationship is useless
    else if ((phase(1:1) .eq. 'l') .or. (phase(1:1) .eq. 'L')) then
        particle_mass = 1.d6*4.d0/3.d0*pi*(p_size/2.d0)**3.d0
    end if

    return
end subroutine mass_size_rel
