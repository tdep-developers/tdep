submodule(phonondamping) phonondamping_dos
!! Helper routines that did not fit elsewhere. Maybe deprecated.
implicit none
contains

!> calculate things as a density of states
module subroutine get_intensity_as_dos(pd, qpd, drd, uc, fc, fct, fcf, ise, sf, qp, dr, opts, tmr, mw, mem, verbosity)
    !> phonon density of states
    type(lo_phonon_dos), intent(out) :: pd
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qpd
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(inout) :: drd
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> self-energy
    type(lo_interpolated_selfenergy), intent(in) :: ise
    !> spectral functions on grid
    type(lo_spectralfunction_helper), intent(out) :: sf
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> timer
    type(lo_timer), intent(inout) :: tmr
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! stuff for the DOS mesh
    real(r8) :: timer, t0, t1
    real(r8), dimension(3), parameter :: qdir = [1.0_r8, 0.0_r8, 0.0_r8]
    integer :: solrnk

    timer = walltime()
    t0 = timer
    t1 = timer

    ! Start the timer
    init: block
        integer :: i

        ! make space in the dos
        pd%n_atom = uc%na
        pd%n_mode = uc%na*3
        pd%n_dos_point = opts%nf
        pd%dosmin = 0.0_r8
        pd%dosmax = -1.0_r8
        pd%integrationtype = -1   ! no choice here
        pd%smearing_prefactor = opts%sigma
        allocate (pd%omega(pd%n_dos_point))
        allocate (pd%dos(pd%n_dos_point))
        allocate (pd%pdos_site(pd%n_dos_point, pd%n_atom))
        allocate (pd%pdos_mode(pd%n_dos_point, pd%n_mode))
        pd%omega = 0.0_r8
        pd%dos = 0.0_r8
        pd%pdos_site = 0.0_r8
        pd%pdos_mode = 0.0_r8

        ! While I'm at it, make space for shifts and linewidths.
        do i = 1, qpd%n_irr_point
            allocate (drd%iq(i)%shift3(drd%n_mode))
            allocate (drd%iq(i)%shift4(drd%n_mode))
            allocate (drd%iq(i)%linewidth(drd%n_mode))
            drd%iq(i)%shift3 = 0.0_r8
            drd%iq(i)%shift4 = 0.0_r8
            drd%iq(i)%linewidth = 0.0_r8
        end do

        ! Read/write rank
        solrnk = mw%n - 1

        ! Initialize the helper (only on one rank, to save a little space)
        if (mw%r .eq. solrnk) then
            sf%n_irr_qpoint = qpd%n_irr_point
            sf%n_energy = opts%nf
            sf%max_energy = opts%maxf
            sf%temperature = opts%temperature
            allocate (sf%qv(3, sf%n_irr_qpoint))
            allocate (sf%spectralfunction(sf%n_energy, dr%n_mode, sf%n_irr_qpoint))
            sf%qv = 0.0_r8
            sf%spectralfunction = 0.0_r8
        end if

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'CALCULATING BROADENED PHONON DOS'
        end if

    end block init

    ! We can start by calculating the lineshape at all the points we need.
    lshp: block
        type(lo_spectralfunction_helper) :: isf
        type(lo_phonon_selfenergy) :: se
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(:, :), allocatable :: buf_shift_3rd, buf_shift_4th, buf_linewidth, buf_spectralfunction
        real(r8), dimension(:), allocatable :: buf0, buf1, siteproj
        real(r8) :: sigma
        integer :: imode, iq, iatom

        t0 = walltime()
        ! Make some space
        call mem%allocate(buf0, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf1, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(siteproj, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_shift_3rd, [dr%n_mode, qpd%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_shift_4th, [dr%n_mode, qpd%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_linewidth, [dr%n_mode, qpd%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_spectralfunction, [opts%nf, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf0 = 0.0_r8
        buf1 = 0.0_r8
        siteproj = 0.0_r8
        buf_shift_3rd = 0.0_r8
        buf_shift_4th = 0.0_r8
        buf_linewidth = 0.0_r8
        buf_spectralfunction = 0.0_r8

        !if ( verbosity .gt. 0 ) call lo_progressbar_init()
        do iq = 1, qpd%n_irr_point

            if (verbosity .gt. 0) then
                write (*, *) '... starting q-point', iq, 'out of', qpd%n_irr_point
            end if

            ! get the self-energy
            call se%generate(qpd%ip(iq), qdir, drd%iq(iq), uc, fc, fct, fcf, ise, isf, qp, dr, opts, tmr, mw, mem, verbosity=-1)

            buf_spectralfunction = 0.0_r8
            ! Add together to get the DOS
            do imode = 1, dr%n_mode
                ! Make it parallel
                if (mod(imode, mw%n) .ne. mw%r) cycle
                ! Skip acoustic at Gamma
                if (drd%iq(iq)%omega(imode) .lt. lo_freqtol) cycle
                ! Get the spectral function
                buf1 = 0.0_r8
                call evaluate_spectral_function(se%energy_axis, se%im_3ph(:, imode) + se%im_iso(:, imode), se%re_3ph(:, imode) + se%re_4ph(:, imode), drd%iq(iq)%omega(imode), buf1)
                buf1 = buf1/lo_trapezoid_integration(se%energy_axis, buf1)

                ! Stash the unsmeared spectral function
                buf_spectralfunction(:, imode) = buf1

                ! Smear spectral function
                sigma = qpd%adaptive_sigma(qpd%ip(iq)%radius, drd%iq(iq)%vel(:, imode), drd%default_smearing(imode), pd%smearing_prefactor)
                call gaussian_smear_spectral_function(se%energy_axis, sigma, buf1)
                ! Normalize it after smearing? Probably a good idea.
                buf1 = buf1/lo_trapezoid_integration(se%energy_axis, buf1)

                ! site-projections
                do iatom = 1, uc%na
                    cv0 = drd%iq(iq)%egv((iatom - 1)*3 + 1:iatom*3, imode)
                    siteproj(iatom) = abs(dot_product(cv0, conjg(cv0)))
                end do

                ! Add together in the right place
                pd%pdos_mode(:, imode) = pd%pdos_mode(:, imode) + buf1*qpd%ip(iq)%integration_weight
                do iatom = 1, uc%na
                    pd%pdos_site(:, iatom) = pd%pdos_site(:, iatom) + buf1*siteproj(iatom)*qpd%ip(iq)%integration_weight
                end do

                ! Store shifts and linewidth
                buf_shift_3rd(imode, iq) = lo_linear_interpolation(se%energy_axis, se%re_3ph(:, imode), drd%iq(iq)%omega(imode))
                buf_shift_4th(imode, iq) = lo_linear_interpolation(se%energy_axis, se%re_4ph(:, imode), drd%iq(iq)%omega(imode))
                buf_linewidth(imode, iq) = lo_linear_interpolation(se%energy_axis, se%im_3ph(:, imode), drd%iq(iq)%omega(imode)) + lo_linear_interpolation(se%energy_axis, se%im_iso(:, imode), drd%iq(iq)%omega(imode))
            end do

            ! Store spectral function
            call mw%allreduce('sum', buf_spectralfunction)

            if (mw%r .eq. solrnk) then
                ! Store q-coordinate
                sf%qv(:, iq) = qpd%ip(iq)%r
                ! Store spectral function
                sf%spectralfunction(:, :, iq) = buf_spectralfunction
            end if

            !if ( verbosity .gt. 0 .and. iq .lt. qpd%n_irr_point ) then
            !    call lo_progressbar(' ... lineshapes across q-mesh',iq,qpd%n_irr_point,walltime()-t0)
            !endif
        end do
        call mw%allreduce('sum', pd%pdos_site)
        call mw%allreduce('sum', pd%pdos_mode)
        call mw%allreduce('sum', buf_shift_3rd)
        call mw%allreduce('sum', buf_shift_4th)
        call mw%allreduce('sum', buf_linewidth)

        ! Store shifts and widths
        do iq = 1, qpd%n_irr_point
            drd%iq(iq)%shift3 = buf_shift_3rd(:, iq)
            drd%iq(iq)%shift4 = buf_shift_4th(:, iq)
            drd%iq(iq)%linewidth = buf_linewidth(:, iq)
        end do

        ! Make sure to store the energy axis
        pd%omega = se%energy_axis
        pd%dosmax = maxval(se%energy_axis)

        call mem%deallocate(buf0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(siteproj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_shift_3rd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_shift_4th, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_linewidth, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        !if ( verbosity .gt. 0 ) then
        !    t1=walltime()
        !    call lo_progressbar(' ... lineshapes across q-mesh',qpd%n_irr_point,qpd%n_irr_point,t1-t0)
        !    t0=t1
        !endif

    end block lshp

    ! Now fix the DOS so that it looks nice.
    makenice: block
        integer :: i, j, k
        real(r8), dimension(:, :), allocatable :: sitebuf
        real(r8) :: f0, f1, dostol

        ! I might have gotten a contribution at zero frequency due to the smearing. That has to
        ! be removed. I subtract the line that goes from the dos at 0 to the max frequency.
        do j = 1, pd%n_mode
            f1 = pd%pdos_mode(1, j)
            do i = 1, pd%n_dos_point
                f0 = (pd%n_dos_point - i)*1.0_r8/((pd%n_dos_point - 1)*1.0_r8)
                f0 = f0**2
                pd%pdos_mode(i, j) = pd%pdos_mode(i, j) - f0*f1
                if (pd%pdos_mode(i, j) .lt. lo_freqtol) pd%pdos_mode(i, j) = 0.0_r8
            end do
            pd%pdos_mode(1, j) = 0.0_r8
        end do
        do j = 1, pd%n_atom
            f1 = pd%pdos_site(1, j)
            do i = 1, pd%n_dos_point
                f0 = (pd%n_dos_point - i)*1.0_r8/((pd%n_dos_point - 1)*1.0_r8)
                f0 = f0**2
                pd%pdos_site(i, j) = pd%pdos_site(i, j) - f0*f1
                if (pd%pdos_site(i, j) .lt. lo_freqtol) pd%pdos_site(i, j) = 0.0_r8
            end do
            pd%pdos_site(1, j) = 0.0_r8
        end do

        ! Sum up contributions and normalize things
        pd%dos = 0.0_r8
        do j = 1, pd%n_mode
            f0 = 1.0_r8/lo_trapezoid_integration(pd%omega, pd%pdos_mode(:, j))
            pd%dos = pd%dos + pd%pdos_mode(:, j)
            pd%pdos_mode(:, j) = pd%pdos_mode(:, j)*f0
        end do
        pd%dos(1) = 0.0_r8
        ! Get a tolerance where to slice off the phonon dos.
        dostol = lo_mean(pd%dos)*lo_sqtol
        pd%dos = lo_chop(pd%dos, dostol)
        ! Normalize the total
        f0 = lo_trapezoid_integration(pd%omega, pd%dos)
        if (verbosity .gt. 0) write (*, *) '... raw normalization (should be 1):', f0/dr%n_mode
        pd%dos = pd%dos*dr%n_mode/f0
        ! Adjust the mode projected so that they sum up to the total
        do i = 1, pd%n_dos_point
            f0 = sum(pd%pdos_mode(i, :))
            if (f0 .gt. dostol) then
                pd%pdos_mode(i, :) = pd%pdos_mode(i, :)*pd%dos(i)/f0
            else
                pd%pdos_mode(i, :) = 0.0_r8
            end if
        end do

        ! Fix the degeneracy of the site-projected
        call mem%allocate(sitebuf, [pd%n_dos_point, pd%n_atom], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        sitebuf = 0.0_r8
        do i = 1, pd%n_atom
            ! enfore site degeneracy
            do j = 1, uc%sym%degeneracy(i)
                k = uc%sym%degenerate_atom(j, i)
                sitebuf(:, i) = sitebuf(:, i) + pd%pdos_site(:, k)/(1.0_r8*uc%sym%degeneracy(j))
            end do
            f0 = lo_trapezoid_integration(pd%omega, sitebuf(:, i))
            sitebuf(:, i) = sitebuf(:, i)*3.0_r8/f0
        end do
        pd%pdos_site = sitebuf
        call mem%deallocate(sitebuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        ! normalize the projections so that they add up to the total
        do i = 1, pd%n_dos_point
            f0 = sum(pd%pdos_site(i, :))
            if (f0 .gt. dostol) then
                pd%pdos_site(i, :) = pd%pdos_site(i, :)*pd%dos(i)/f0
            else
                pd%pdos_site(i, :) = 0.0_r8
            end if
        end do
    end block makenice
end subroutine

!> calculate things as a density of states
module subroutine get_intensity_as_dos_interp(pd, tc, qp, dr, uc, ise, opts, mw, mem, verbosity)
    !> phonon density of states
    type(lo_phonon_dos), intent(out) :: pd
    !> thermal conductivity container
    type(lo_thermal_conductivity), intent(out) :: tc
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> self-energy
    type(lo_interpolated_selfenergy), intent(in) :: ise
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! stuff for the DOS mesh
    real(r8) :: timer, t0, t1
    real(r8), dimension(3), parameter :: qdir = [1.0_r8, 0.0_r8, 0.0_r8]

    timer = walltime()
    t0 = timer
    t1 = timer

    ! Make some space for the DOS
    init: block
        integer :: i

        ! make space in the dos
        pd%n_atom = uc%na
        pd%n_mode = uc%na*3
        pd%n_dos_point = ise%n_energy
        pd%dosmin = 0.0_r8
        pd%dosmax = -1.0_r8
        pd%integrationtype = -1
        pd%smearing_prefactor = opts%sigma
        allocate (pd%omega(pd%n_dos_point))
        allocate (pd%dos(pd%n_dos_point))
        allocate (pd%pdos_site(pd%n_dos_point, pd%n_atom))
        allocate (pd%pdos_mode(pd%n_dos_point, pd%n_mode))
        pd%omega = ise%energy
        pd%dos = 0.0_r8
        pd%pdos_site = 0.0_r8
        pd%pdos_mode = 0.0_r8

        ! While I'm at it, make space for shifts and linewidths.
        do i = 1, qp%n_irr_point
            allocate (dr%iq(i)%shift3(dr%n_mode))
            allocate (dr%iq(i)%shift4(dr%n_mode))
            allocate (dr%iq(i)%linewidth(dr%n_mode))
            dr%iq(i)%shift3 = 0.0_r8
            dr%iq(i)%shift4 = 0.0_r8
            dr%iq(i)%linewidth = 0.0_r8
        end do

        ! Make some space for the thermal conductivity container?
        ! tc%n_mode=dr%n_mode
        ! tc%n_energy=ise%n_energy
        ! allocate(tc%omega( tc%n_energy ))
        ! allocate(tc%spectral_kappa(ise%n_energy,3,3,dr%n_mode,dr%n_mode))
        ! allocate(tc%kappa(3,3,dr%n_mode,dr%n_mode,qp%n_irr_point))
        ! allocate(tc%lifetime(dr%n_mode,dr%n_mode,qp%n_irr_point))
        ! allocate(tc%mean_free_path(dr%n_mode,dr%n_mode,qp%n_irr_point))
        ! tc%omega = ise%energy
        ! tc%spectral_kappa=0.0_r8
        ! tc%kappa=0.0_r8
        ! tc%lifetime=0.0_r8
        ! tc%mean_free_path=0.0_r8

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'CALCULATING BROADENED PHONON DOS WITH INTERPOLATED SELFENERGY'
        end if
    end block init

    ! We can start by calculating the lineshape at all the points we need.
    lshp: block
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(:, :), allocatable :: buf_shift_3rd, buf_shift_4th, buf_linewidth, buf_tau
        real(r8), dimension(:, :), allocatable :: buf_spectralfunction, buf_Re, buf_Im
        real(r8), dimension(:), allocatable :: buf0, buf1, siteproj, taper
        real(r8) :: sigma
        integer :: imode, iq, iatom

        t0 = walltime()
        ! Make some space
        call mem%allocate(buf0, ise%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf1, ise%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(taper, ise%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(siteproj, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_shift_3rd, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_shift_4th, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_linewidth, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_tau, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_spectralfunction, [ise%n_energy, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_Re, [ise%n_energy, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_Im, [ise%n_energy, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        buf0 = 0.0_r8
        buf1 = 0.0_r8
        taper = 0.0_r8
        siteproj = 0.0_r8
        buf_shift_3rd = 0.0_r8
        buf_shift_4th = 0.0_r8
        buf_linewidth = 0.0_r8
        buf_tau = 0.0_r8
        buf_spectralfunction = 0.0_r8
        buf_Re = 0.0_r8
        buf_Im = 0.0_r8

        call taperfn_im(ise%energy, dr%omega_max, dr%omega_min, taper)

        if (verbosity .gt. 0) call lo_progressbar_init()
        qploop: do iq = 1, qp%n_irr_point
            ! Make it parallel
            if (mod(iq, mw%n) .ne. mw%r) cycle

            ! get the (diagonal) self-energy
            call ise%diagonal_selfenergy(uc, qp%ip(iq)%r, dr%iq(iq)%omega, dr%iq(iq)%egv, buf_Re, buf_Im)
            !buf_Im=max(buf_Im,0.0_r8)

            ! Add together to get the DOS
            do imode = 1, dr%n_mode
                ! Skip acoustic at Gamma
                if (dr%iq(iq)%omega(imode) .lt. lo_freqtol) cycle
                ! Get the lineshape
                call evaluate_spectral_function(pd%omega, buf_Im(:, imode), buf_Re(:, imode), dr%iq(iq)%omega(imode), buf0)
                ! Taper, smear and normalize
                buf1 = buf0*taper
                sigma = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, imode), dr%default_smearing(imode), pd%smearing_prefactor)
                call gaussian_smear_spectral_function(ise%energy, sigma, buf1)
                buf1 = buf1/lo_trapezoid_integration(ise%energy, buf1)

                ! site-projections
                do iatom = 1, uc%na
                    cv0 = dr%iq(iq)%egv((iatom - 1)*3 + 1:iatom*3, imode)
                    siteproj(iatom) = abs(dot_product(cv0, conjg(cv0)))
                end do

                ! Add together in the right place
                pd%pdos_mode(:, imode) = pd%pdos_mode(:, imode) + buf1*qp%ip(iq)%integration_weight
                do iatom = 1, uc%na
                    pd%pdos_site(:, iatom) = pd%pdos_site(:, iatom) + buf1*siteproj(iatom)*qp%ip(iq)%integration_weight
                end do

                ! While I'm at it, get the thermal conductivity?
                ! Normalize the spectral function to be on the safe side.
                buf0 = buf0/lo_trapezoid_integration(ise%energy, buf0)
                !call spectral_tau(ise%energy,buf0,ise%temperature,buf_tau(imode,iq))

                ! Store shifts and linewidth
                buf1 = lo_gauss(ise%energy, dr%iq(iq)%omega(imode), sigma)
                buf1 = buf1/sum(buf1)
                buf1 = buf1*buf_Im(:, imode)
                buf_shift_3rd(imode, iq) = lo_linear_interpolation(ise%energy, buf_Re(:, imode), dr%iq(iq)%omega(imode))
                buf_linewidth(imode, iq) = lo_linear_interpolation(ise%energy, buf_Im(:, imode), dr%iq(iq)%omega(imode))
            end do

            if (verbosity .gt. 0 .and. iq .lt. qp%n_irr_point) then
                call lo_progressbar(' ... lineshapes across q-mesh', iq, qp%n_irr_point, walltime() - t0)
            end if
        end do qploop

        call mw%allreduce('sum', pd%pdos_site)
        call mw%allreduce('sum', pd%pdos_mode)
        call mw%allreduce('sum', buf_shift_3rd)
        call mw%allreduce('sum', buf_shift_4th)
        call mw%allreduce('sum', buf_linewidth)
        call mw%allreduce('sum', buf_tau)

        ! Store shifts and widths
        do iq = 1, qp%n_irr_point
            dr%iq(iq)%shift3 = buf_shift_3rd(:, iq)
            ! DANGER ZONE! REUSING SPACE I SHOULD NOT REUSE.
            dr%iq(iq)%shift4 = buf_tau(:, iq)
            dr%iq(iq)%linewidth = buf_linewidth(:, iq)
        end do

        call mem%deallocate(buf0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(taper, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(siteproj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_shift_3rd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_shift_4th, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_linewidth, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_tau, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_Re, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_Im, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... lineshapes across q-mesh', qp%n_irr_point, qp%n_irr_point, t1 - t0)
            t0 = t1
        end if
    end block lshp

    ! Now fix the DOS so that it looks nice.
    ! @TODO move this to type_phonon_dos
    makenice: block
        integer :: i, j, k
        real(r8), dimension(:, :), allocatable :: sitebuf
        real(r8) :: f0, f1, dostol

        ! I might have gotten a contribution at zero frequency due to the smearing. That has to
        ! be removed. I subtract the line that goes from the dos at 0 to the max frequency.
        do j = 1, pd%n_mode
            f1 = pd%pdos_mode(1, j)
            do i = 1, pd%n_dos_point
                f0 = (pd%n_dos_point - i)*1.0_r8/((pd%n_dos_point - 1)*1.0_r8)
                f0 = f0**2
                pd%pdos_mode(i, j) = pd%pdos_mode(i, j) - f0*f1
                if (pd%pdos_mode(i, j) .lt. lo_freqtol) pd%pdos_mode(i, j) = 0.0_r8
            end do
            pd%pdos_mode(1, j) = 0.0_r8
        end do
        do j = 1, pd%n_atom
            f1 = pd%pdos_site(1, j)
            do i = 1, pd%n_dos_point
                f0 = (pd%n_dos_point - i)*1.0_r8/((pd%n_dos_point - 1)*1.0_r8)
                f0 = f0**2
                pd%pdos_site(i, j) = pd%pdos_site(i, j) - f0*f1
                if (pd%pdos_site(i, j) .lt. lo_freqtol) pd%pdos_site(i, j) = 0.0_r8
            end do
            pd%pdos_site(1, j) = 0.0_r8
        end do

        ! Sum up contributions and normalize things
        pd%dos = 0.0_r8
        do j = 1, pd%n_mode
            f0 = 1.0_r8/lo_trapezoid_integration(pd%omega, pd%pdos_mode(:, j))
            pd%dos = pd%dos + pd%pdos_mode(:, j)
            pd%pdos_mode(:, j) = pd%pdos_mode(:, j)*f0
        end do
        pd%dos(1) = 0.0_r8
        ! Get a tolerance where to slice off the phonon dos.
        dostol = lo_mean(pd%dos)*lo_sqtol
        pd%dos = lo_chop(pd%dos, dostol)
        ! Normalize the total
        f0 = lo_trapezoid_integration(pd%omega, pd%dos)
        if (verbosity .gt. 0) write (*, *) '... raw normalization (should be 1):', f0/dr%n_mode
        pd%dos = pd%dos*dr%n_mode/f0
        ! Adjust the mode projected so that they sum up to the total
        do i = 1, pd%n_dos_point
            f0 = sum(pd%pdos_mode(i, :))
            if (f0 .gt. dostol) then
                pd%pdos_mode(i, :) = pd%pdos_mode(i, :)*pd%dos(i)/f0
            else
                pd%pdos_mode(i, :) = 0.0_r8
            end if
        end do

        ! Fix the degeneracy of the site-projected
        call mem%allocate(sitebuf, [pd%n_dos_point, pd%n_atom], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        sitebuf = 0.0_r8
        do i = 1, pd%n_atom
            ! enfore site degeneracy
            do j = 1, uc%sym%degeneracy(i)
                k = uc%sym%degenerate_atom(j, i)
                sitebuf(:, i) = sitebuf(:, i) + pd%pdos_site(:, k)/(1.0_r8*uc%sym%degeneracy(j))
            end do
            f0 = lo_trapezoid_integration(pd%omega, sitebuf(:, i))
            sitebuf(:, i) = sitebuf(:, i)*3.0_r8/f0
        end do
        pd%pdos_site = sitebuf
        call mem%deallocate(sitebuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        ! normalize the projections so that they add up to the total
        do i = 1, pd%n_dos_point
            f0 = sum(pd%pdos_site(i, :))
            if (f0 .gt. dostol) then
                pd%pdos_site(i, :) = pd%pdos_site(i, :)*pd%dos(i)/f0
            else
                pd%pdos_site(i, :) = 0.0_r8
            end if
        end do
    end block makenice
end subroutine

end submodule
