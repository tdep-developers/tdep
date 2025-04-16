
!> Generate the self-energy at a single q-point
subroutine generate(se, qpoint, qdir, wp, uc, fc, fct, fcf, ise, isf, qp, dr, opts, tmr, mw, mem, verbosity)
    !> self energy
    class(lo_phonon_selfenergy), intent(out) :: se
    !> qpoint of interest
    class(lo_qpoint), intent(in) :: qpoint
    !> direction of probe?
    real(r8), dimension(3), intent(in) :: qdir
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> interpolated self-energy
    type(lo_interpolated_selfenergy), intent(in) :: ise
    !> tabulated and smeared spectral functions
    type(lo_spectralfunction_helper), intent(in) :: isf
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> timer (start and stop outside this routine)
    type(lo_timer), intent(inout) :: tmr
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot
    integer, intent(in) :: verbosity

    ! scattering rates container
    type(lo_listofscatteringrates) :: sr
    !> harmonic properties at Gamma
    type(lo_phonon_dispersions_qpoint) :: gp

    call mem%tick()
    call tmr%tick()
    ! First thing to do is to initialize all arrays, and set all options
    setopts: block
        integer :: i

        se%integrationtype = opts%integrationtype   ! How to integrate
        ! Number of points on the energy axis
        if (se%integrationtype .eq. 4) then
            se%n_energy = ise%n_energy
        else
            se%n_energy = opts%nf
        end if
        se%n_mode = uc%na*3                ! Number of modes
        se%smearing_prefactor = opts%sigma             ! Adjustment of default smearing parameter
        se%isotope_scattering = opts%isotopescattering ! include isotope scattering?
        se%thirdorder_scattering = opts%thirdorder        ! include threephonon term?
        se%fourthorder_scattering = opts%fourthorder       ! include fourphonon term?
        se%diagonal = opts%diagonal          ! calculate only the diagonal part of the self-energy
        se%qdir = qdir                   ! Make note on the probe direction
        ! Make sure we know exactly when to skip symmetry things
        se%skipsym = .false.
        if (opts%qpointpath) se%skipsym = .true.
        if (opts%grid) se%skipsym = .true.
        ! Make space for all the necessary arrays
        allocate (se%energy_axis(se%n_energy))
        allocate (se%im_3ph(se%n_energy, se%n_mode))
        allocate (se%im_iso(se%n_energy, se%n_mode))
        allocate (se%re_3ph(se%n_energy, se%n_mode))
        allocate (se%re_4ph(se%n_energy, se%n_mode))
        ! Get the energy range
        if (se%integrationtype .eq. 4) then
            se%energy_axis = ise%energy
        else
            se%energy_axis = 0.0_r8
            call lo_linspace(0.0_r8, opts%maxf, se%energy_axis)
        end if

        se%im_3ph = 0.0_r8
        se%im_iso = 0.0_r8
        se%re_3ph = 0.0_r8
        se%re_4ph = 0.0_r8

        ! Now things should be clean and nice. Maybe say what we are about to do?
        if (verbosity .gt. 1) then
            write (*, *) '         liso:', se%isotope_scattering
            write (*, *) '  lthirdorder:', se%thirdorder_scattering
            write (*, *) ' lfourthorder:', se%fourthorder_scattering
            write (*, *) '        sigma:', se%smearing_prefactor*lo_mean(dr%default_smearing)*lo_frequency_hartree_to_thz, ' Thz'
            write (*, *) '  temperature: ', tochar(opts%temperature), ' K'
            write (*, *) '  frequencies:'
            do i = 1, dr%n_mode
                write (*, *) '    mode ', tochar(i, -3), ', omega: ', tochar(wp%omega(i)*lo_frequency_Hartree_to_THz, 6), ' THz'
            end do
        end if

        ! Get harmonic properties at Gamma with the proper direction for this q-point.
        call gp%generate(fc, uc, mem, qvec=[0.0_r8, 0.0_r8, 0.0_r8], qdirection=qdir)

        ! If a closed grid, get Gamma consistent with the grid.
        if (opts%grid) then
            do i = 1, qp%n_irr_point
                if (norm2(qp%ip(i)%r) .gt. lo_sqtol) cycle
                gp%omega = dr%iq(i)%omega
                gp%egv = dr%iq(i)%egv
                gp%vel = dr%iq(i)%vel
                exit
            end do
        end if

        ! And some space for auxiliary information
        allocate (se%xmid(se%n_mode))
        allocate (se%xlo(se%n_mode))
        allocate (se%xhi(se%n_mode))
        allocate (se%scalingfactor(se%n_mode))
        allocate (se%tau(se%n_mode))
        se%xmid = 0.0_r8
        se%xlo = 0.0_r8
        se%xhi = 0.0_r8
        se%scalingfactor = 0.0_r8
        se%tau = 0.0_r8
    end block setopts

    call tmr%tock('self-energy initialization')

    ! Now calculate all the matrix elements
    if (se%integrationtype .eq. 4) then
        ! Keep the eigenvectors after
        call sr%generate(qpoint, wp, gp, qp, dr, uc, fc, fct, se%isotope_scattering, se%thirdorder_scattering, se%skipsym, .true., opts%grid, tmr, mw, mem, verbosity)
    else
        ! Do not keep eigenvectors
        call sr%generate(qpoint, wp, gp, qp, dr, uc, fc, fct, se%isotope_scattering, se%thirdorder_scattering, se%skipsym, .false., opts%grid, tmr, mw, mem, verbosity)
    end if

    ! Start to actually calculate stuff, isotope things first
    if (se%isotope_scattering) then
        select case (se%integrationtype)
        case (1:2)
            call isotope_imaginary_selfenergy_gaussian(wp, qp, dr, se, sr, mw, mem, verbosity)
        case (3)
            call isotope_imaginary_selfenergy_tetrahedron(wp, qp, dr, se, sr, mw, mem, verbosity)
        case (4:5)
            ! Nothing, taken care of with the other integrals
        case default
            call lo_stop_gracefully(['Unknown integration type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end select
    end if

    call tmr%tock('isotope integrals')

    ! Then three-phonon things
    if (se%thirdorder_scattering) then
        select case (se%integrationtype)
        case (1:2)
            call threephonon_imaginary_selfenergy_gaussian(wp, se, sr, qp, dr, opts%temperature, mw, mem, verbosity)
        case (3)
            call threephonon_imaginary_selfenergy_tetrahedron(wp, qp, sr, dr, opts%temperature, se, mw, mem, verbosity)
        case (4)
            call convolution_imaginary_selfenergy(se, wp, qp, dr, sr, ise, isf, uc, opts%temperature, mw, mem, verbosity)
        case (5)
            call convolution_imaginary_selfenergy(se, wp, qp, dr, sr, ise, isf, uc, opts%temperature, mw, mem, verbosity)
        case default
            call lo_stop_gracefully(['Unknown integration type'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end select
    end if

    call tmr%tock('three-phonon integrals')

    ! Maybe fourthorder things
    if (se%fourthorder_scattering) then
        fourthorder: block
            real(r8), dimension(:), allocatable :: delta
            integer :: j

            allocate (delta(dr%n_mode))
            call fourphonon_selfenergy(qpoint, wp, gp, qp, uc, opts%temperature, dr, fcf, delta, se%skipsym, mw, mem, verbosity)
            do j = 1, se%n_mode
                se%re_4ph(:, j) = delta(j)
            end do
            deallocate (delta)
        end block fourthorder
    end if

    call tmr%tock('four-phonon self-energy')

    ! finalize to ensure that it's reasonable.
    sanity: block
        ! Make sure the selfenergy is zero where it's supposed to be. First
        ! make sure it's at least not negative.
        se%im_3ph = max(se%im_3ph, 0.0_r8)
        se%im_iso = max(se%im_iso, 0.0_r8)
        ! zero at zero
        se%im_3ph(1, :) = 0.0_r8
        se%im_iso(1, :) = 0.0_r8
        ! zero at the end
        se%im_3ph(se%n_energy, :) = 0.0_r8
        se%im_iso(se%n_energy, :) = 0.0_r8
    end block sanity

    ! Kramers-Kronig-transform the imaginary part to get the real.
    if (se%thirdorder_scattering) then
        kktransform: block
            real(r8) :: pref = 2.0_r8/lo_pi
            complex(r8), dimension(:), allocatable :: z0
            complex(r8) :: eta
            real(r8), dimension(:), allocatable :: x, xs, y0
            real(r8) :: xp, dlx
            integer :: imode, ie, ctr

            call mem%allocate(x, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(xs, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(z0, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(y0, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            x = se%energy_axis                  ! copy of x-axis
            xs = x**2                           ! x^2, precalculated
            dlx = x(2) - x(1)                     ! prefactor for riemann integral
            eta = lo_imag*(x(2) - x(1))*1E-8_r8   ! small imaginary thing to avoid divergence

            ctr = 0
            se%re_3ph = 0.0_r8
            do imode = 1, se%n_mode
            do ie = 1, se%n_energy
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                xp = x(ie)
                y0 = se%im_3ph(:, imode)*x
                ! To anyone reading this code:
                ! The correct way to compute the principal value would be
                ! z0 = (xp + eta)**2 - xs
                ! However, it doesn't matter for the Im -> Re transform because we
                ! don't really hit 0.
                z0 = xp**2 - xs + eta
                y0 = real(y0/z0, r8)
                y0(1) = y0(1)*0.5_r8
                y0(se%n_energy) = y0(se%n_energy)*0.5_r8
                se%re_3ph(ie, imode) = sum(y0)*dlx
            end do
            end do
            call mw%allreduce('sum', se%re_3ph)
            se%re_3ph = se%re_3ph*pref

            call mem%deallocate(x, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(xs, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(z0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(y0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block kktransform
        call tmr%tock('Kramers-Kronig transformation')
    end if

    ! Normalize the spectral functions
    normalize: block
        real(r8), parameter :: integraltol = 1E-13_r8
        real(r8), dimension(:), allocatable :: yim, yre, ysf
        real(r8) :: f0, f1, f2, f3
        integer :: imode

        call mem%allocate(yim, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(yre, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ysf, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        yim = 0.0_r8
        yre = 0.0_r8
        ysf = 0.0_r8

        se%scalingfactor = 0.0_r8
        se%xmid = 0.0_r8
        se%xlo = 0.0_r8
        se%xhi = 0.0_r8
        se%tau= 0.0_r8
        do imode = 1, dr%n_mode
            if (wp%omega(imode) .lt. lo_freqtol) cycle
            if (mod(imode, mw%n) .ne. mw%r) cycle

            ! Evaluate rough spectral function
            yim = se%im_3ph(:, imode) + se%im_iso(:, imode)
            yre = se%re_3ph(:, imode) + se%re_4ph(:, imode)
            call evaluate_spectral_function(se%energy_axis, yim, yre, wp%omega(imode), ysf)
            ! Find peaks in the spectral function
            call find_spectral_function_max_and_fwhm(wp%omega(imode), dr%omega_max, se%energy_axis, yim, yre, se%xmid(imode), se%xlo(imode), se%xhi(imode))
            ! Integrate every which way to get normalization.
            call integrate_spectral_function(se%energy_axis, wp%omega(imode), yim, yre, se%xmid(imode), se%xlo(imode), se%xhi(imode), &
                                             1.0_r8, opts%temperature, integraltol, 1.0_r8, f0, f1, f2, f3)
            se%scalingfactor(imode) = 1.0_r8/f0
            se%tau(imode) = f1
        end do
        call mw%allreduce('sum', se%scalingfactor)
        call mw%allreduce('sum', se%xmid)
        call mw%allreduce('sum', se%xhi)
        call mw%allreduce('sum', se%xlo)
        call mw%allreduce('sum', se%tau)

        call mem%deallocate(yim, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(yre, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ysf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block normalize
    call tmr%tock('spectral function normalization')

    ! Make sure I did not mess up with memory management.
    call mem%tock(__FILE__, __LINE__, mw%comm)
end subroutine

!> Generate the self-energy at a single q-point
subroutine generate_interp(se, qpoint, qdir, wp, uc, ise, opts, mw, mem)
    !> self energy
    class(lo_phonon_selfenergy), intent(out) :: se
    !> qpoint of interest
    class(lo_qpoint), intent(in) :: qpoint
    !> direction of probe?
    real(r8), dimension(3), intent(in) :: qdir
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> self-energy
    type(lo_interpolated_selfenergy), intent(in) :: ise
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> fancy interpolation thing
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    call mem%tick()

    se%n_energy = ise%n_energy           ! Number of points on the energy axis
    se%n_mode = uc%na*3                ! Number of bands
    se%smearing_prefactor = opts%sigma             ! Adjustment of default smearing parameter
    se%integrationtype = opts%integrationtype   ! How to integrate
    se%isotope_scattering = opts%isotopescattering ! include isotope scattering?
    se%thirdorder_scattering = opts%thirdorder        ! include threephonon term?
    se%fourthorder_scattering = opts%fourthorder       ! include fourphonon term?
    se%diagonal = opts%diagonal          ! calculate only the diagonal part of the self-energy
    se%skipsym = opts%qpointpath        ! If on a path, skip using symmetry. For now, should reintroduce it later.
    se%qdir = qdir                   ! Make note on the probe direction
    ! Make space for all the necessary arrays
    allocate (se%energy_axis(se%n_energy))
    allocate (se%im_3ph(se%n_energy, se%n_mode))
    allocate (se%im_iso(se%n_energy, se%n_mode))
    allocate (se%re_3ph(se%n_energy, se%n_mode))
    allocate (se%re_4ph(se%n_energy, se%n_mode))
    se%energy_axis = ise%energy
    se%im_3ph = 0.0_r8
    se%im_iso = 0.0_r8
    se%re_3ph = 0.0_r8
    se%re_4ph = 0.0_r8

    call ise%diagonal_selfenergy(uc, qpoint%r, wp%omega, wp%egv, se%re_3ph, se%im_3ph)

    ! Make sure I did not mess up with memory management.
    call mem%tock(__FILE__, __LINE__, mw%comm)
end subroutine
