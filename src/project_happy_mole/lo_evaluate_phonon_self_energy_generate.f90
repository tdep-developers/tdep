submodule(lo_evaluate_phonon_self_energy) lo_evaluate_phonon_self_energy_generate
implicit none
contains

!> Generate the self-energy at a single q-point
module subroutine generate(se, qpoint, qdir, wp, uc, fc, fct, fcf, ise, qp, dr, &
    temperature, max_energy, n_energy, integrationtype, sigma,&
    use_isotope, use_thirdorder, use_fourthorder, offdiagonal, &
    tmr, mw, mem, verbosity)
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
    !> pre-calculated self-energies
    type(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8), intent(in) :: temperature
    !> max energy to consider
    real(r8), intent(in) :: max_energy
    !> number of energies
    integer, intent(in) :: n_energy
    !> how are we going to integrate
    integer, intent(in) :: integrationtype
    !> adjustment of smearing parameter
    real(r8), intent(in) :: sigma
    !> include isotope scattering
    logical, intent(in) :: use_isotope
    !> include thirdorder scattering
    logical, intent(in) :: use_thirdorder
    !> include fourthorder scattering
    logical, intent(in) :: use_fourthorder
    !> include off-diagonal terms in the self-energy
    logical, intent(in) :: offdiagonal
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
    ! harmonic properties at Gamma
    type(lo_phonon_dispersions_qpoint) :: gp
    ! is the input q-point on the grid
    logical :: ongrid

    call mem%tick()
    call tmr%tick()
    ! First thing to do is to initialize all arrays, and set all options
    setopts: block
        integer :: i

        se%integrationtype = integrationtype  ! How to integrate

        ! Number of points on the energy axis
        if (se%integrationtype .eq. 5) then
            se%n_energy = ise%n_energy
        else
            se%n_energy = n_energy
        end if
        se%n_mode = uc%na*3                          ! Number of modes
        se%smearing_prefactor = sigma                ! Adjustment of default smearing parameter
        se%isotope_scattering = use_isotope          ! include isotope scattering?
        se%thirdorder_scattering = use_thirdorder    ! include threephonon term?
        se%fourthorder_scattering = use_fourthorder  ! include fourphonon term?
        se%diagonal = .not.offdiagonal               ! calculate only the diagonal part of the self-energy
        se%qdir = qdir                               ! Make note on the probe direction
        ! Make sure we know exactly when to skip symmetry things
        se%skipsym = .true.
        ! Make space for all the necessary arrays
        allocate(se%energy_axis(se%n_energy))
        allocate(se%im_3ph(se%n_energy, se%n_mode))
        allocate(se%im_iso(se%n_energy, se%n_mode))
        allocate(se%re_3ph(se%n_energy, se%n_mode))
        allocate(se%re_4ph(se%n_energy, se%n_mode))
        ! Get the energy range
        if (se%integrationtype .eq. 5) then
            se%energy_axis = ise%omega
        else
            se%energy_axis = 0.0_r8
            call lo_linspace(0.0_r8, max_energy, se%energy_axis)
        end if
        ! Set self-energies to nothing, for now.
        se%im_3ph = 0.0_r8
        se%im_iso = 0.0_r8
        se%re_3ph = 0.0_r8
        se%re_4ph = 0.0_r8


        ! Get harmonic properties at Gamma with the proper direction for this q-point.
        call gp%generate(fc, uc, mem, qvec=[0.0_r8, 0.0_r8, 0.0_r8], qdirection=qdir)

        ! Figure out if the desired q-point resides on the grid.
        select type(qp)
        type is(lo_fft_mesh)
            ongrid=qp%is_point_on_grid( qpoint%r )
        class default
            ongrid=.false.
        end select

        ! If a closed grid, get Gamma consistent with the grid.
        if ( ongrid ) then
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
        se%xmid = 0.0_r8
        se%xlo = 0.0_r8
        se%xhi = 0.0_r8
        se%scalingfactor = 0.0_r8

        ! Now things should be clean and nice. Maybe say what we are about to do?
        if (verbosity .gt. 1) then
            write (*, *) '         liso:', se%isotope_scattering
            write (*, *) '  lthirdorder:', se%thirdorder_scattering
            write (*, *) ' lfourthorder:', se%fourthorder_scattering
            write (*, *) '        sigma:', se%smearing_prefactor*lo_mean(dr%default_smearing)*lo_frequency_hartree_to_thz, ' Thz'
            write (*, *) '  temperature: ', tochar(temperature), ' K'
            write (*, *) '  closed grid: ',ongrid
            write (*, *) '  frequencies:'
            do i = 1, dr%n_mode
                write (*, *) '    mode ', tochar(i, -3), ', omega: ', tochar(wp%omega(i)*lo_frequency_Hartree_to_THz, 6), ' THz'
            end do
        end if

    end block setopts

    call tmr%tock('self-energy initialization')

    ! First we evaluate all the matrix elements
    call sr%generate(qpoint, wp, gp, qp, dr, uc, fc, fct, se%isotope_scattering, se%thirdorder_scattering, se%skipsym, .false., ongrid, tmr, mw, mem, verbosity)

    ! Start to actually calculate stuff, isotope things first
    if (se%isotope_scattering) then
        call isotope_imaginary_selfenergy(wp,se,qp,dr,sr,mw,mem,verbosity)
        call tmr%tock('isotope integrals')
    end if

    ! ! Then three-phonon things
    if (se%thirdorder_scattering) then
        call threephonon_imaginary_selfenergy(wp,se,qp,dr,sr,ise,uc,temperature,mw,mem,verbosity)
        call tmr%tock('three-phonon integrals')
    end if

    ! Maybe fourthorder things
    if (se%fourthorder_scattering) then
        fourthorder: block
            real(r8), dimension(:), allocatable :: delta
            integer :: j

            allocate (delta(dr%n_mode))
            call fourphonon_real_selfenergy(qpoint, wp, gp, qp, uc, temperature, dr, fcf, delta, se%skipsym, mw, mem, verbosity)
            do j = 1, se%n_mode
                se%re_4ph(:, j) = delta(j)
            end do
            deallocate (delta)
        end block fourthorder
        call tmr%tock('four-phonon self-energy')
    end if

    ! finalize to ensure that it's reasonable.
    sanity: block
        ! Make sure the selfenergy is zero where it's supposed to be.
        ! First make sure it's at least not negative.
        se%im_3ph = max(se%im_3ph, 0.0_r8)
        se%im_iso = max(se%im_iso, 0.0_r8)
        ! zero at zero
        se%im_3ph(1, :) = 0.0_r8
        se%im_iso(1, :) = 0.0_r8
        ! zero at the end
        se%im_3ph(se%n_energy, :) = 0.0_r8
        se%im_iso(se%n_energy, :) = 0.0_r8

        ! Here I suppose we could add some kind of tapering function to make
        ! sure it goes smoothly to zero at large Omega.
    end block sanity

    ! Kramers-Kronig-transform the imaginary part to get the real.
    if (se%thirdorder_scattering) then
        kktransform: block
            real(r8), parameter :: pref = 2.0_r8/lo_pi
            complex(r8), dimension(:), allocatable :: z0
            complex(r8) :: eta
            real(r8), dimension(:), allocatable :: x, xs, y0
            real(r8) :: xp, dlx
            integer :: imode, ie, ctr

            call mem%allocate(x, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(xs, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(z0, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(y0, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            x = se%energy_axis                    ! copy of x-axis
            xs = x**2                             ! x^2, precalculated
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

            ! Here we can add a shift so that Re(0)=0. Maybe a good idea?
            do imode=1,se%n_mode
                xp=se%re_3ph(1,imode)
                se%re_3ph(:,imode) = se%re_3ph(:,imode)-xp
            enddo

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
        real(r8) :: f0, f1, f2
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
                                             1.0_r8, temperature, integraltol, f0, f1, f2)
            se%scalingfactor(imode) = 1.0_r8/f0
        end do
        call mw%allreduce('sum', se%scalingfactor)
        call mw%allreduce('sum', se%xmid)
        call mw%allreduce('sum', se%xhi)
        call mw%allreduce('sum', se%xlo)

        call mem%deallocate(yim, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(yre, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ysf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block normalize
    call tmr%tock('spectral function normalization')

    ! Make sure I did not mess up with memory management.
    call mem%tock(__FILE__, __LINE__, mw%comm)
end subroutine

end submodule