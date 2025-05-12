submodule(lo_selfenergy_interpolation) lo_selfenergy_interpolation_grid
use gottochblandat, only: lo_lorentz
implicit none

contains

!> Calculate the spectral function along a path in the BZ
module subroutine spectral_function_grid_interp(ise, uc, fc, grid_density, smearing_prefactor, temperature, tc, pd, mw, mem)
    !> interpolated self-energy thing
    class(lo_interpolated_selfenergy_grid), intent(in) :: ise
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> q-grid density?
    integer, dimension(3), intent(in) :: grid_density
    !> smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    !> temperature
    real(r8), intent(in) :: temperature
    !> thermal conductivity
    type(lo_thermal_conductivity), intent(out) :: tc
    !> phonon dos
    type(lo_phonon_dos), intent(out) :: pd
    !> mpi communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    class(lo_qpoint_mesh), allocatable :: qp
    type(lo_phonon_dispersions) :: dr
    real(r8) :: timer, t0, t1
    integer :: solrnk

    ! Start timers
    timer = walltime()
    t0 = timer
    t1 = timer

    ! Set some basic things
    init: block
        integer :: i

        if (mw%talk) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'Calculating spectral function on a grid with fancy interpolation'
        end if

        ! Which rank do we solve serial things on?
        solrnk = 0

        ! Create the q-point mesh
        call lo_generate_qmesh(qp, uc, grid_density, 'fft', timereversal=.true., headrankonly=.false., mw=mw, mem=mem, verbosity=-1)
        ! Harmonic dispersions
        call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=-1)

        ! Then I guess we make some space for the DOS?
        ! make space in the dos
        pd%n_atom = uc%na
        pd%n_mode = uc%na*3
        pd%n_dos_point = ise%n_energy
        pd%dosmin = 0.0_r8
        pd%dosmax = -1.0_r8
        pd%integrationtype = -1   ! no choice here
        pd%smearing_prefactor = smearing_prefactor
        allocate (pd%omega(pd%n_dos_point))
        allocate (pd%dos(pd%n_dos_point))
        allocate (pd%pdos_site(pd%n_dos_point, pd%n_atom))
        allocate (pd%pdos_mode(pd%n_dos_point, pd%n_mode))
        pd%omega = ise%omega
        pd%dosmax = maxval(pd%omega)
        pd%dos = 0.0_r8
        pd%pdos_site = 0.0_r8
        pd%pdos_mode = 0.0_r8

        if (mw%talk) then
            t1 = walltime()
            write (lo_iou, *) '... made space (', tochar(t1 - t0), 's)'
            t0 = t1
        end if

        call tc%init(dr, ise%n_energy, pd%dosmax, temperature, mw)
    end block init

    ! Calculate self-energies / spectral functions
    specfun: block
        real(r8), parameter :: integraltol=1E-13_r8
        complex(r8), dimension(3) :: cv0
        !real(r8), parameter :: largedistance = 1.0E10_r8*lo_A_to_bohr
        real(r8), dimension(:, :), allocatable :: buf_spectral, buf_spectral_smeared
        real(r8), dimension(:, :), allocatable :: buf_sigmaIm, buf_sigmaRe
        real(r8), dimension(:), allocatable :: buf_taper, siteproj, normalizationfactor, xmid, xlo, xhi
        real(r8), dimension(3), parameter :: qdir = [1.0_r8, 0.0_r8, 0.0_r8]
        real(r8) :: f0, f1, f2, sigma
        integer :: iq, imode, iatom

        ! A little space for temporary buffers
        call mem%allocate(buf_spectral, [pd%n_dos_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_spectral_smeared, [pd%n_dos_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_sigmaIm, [pd%n_dos_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_sigmaRe, [pd%n_dos_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_taper, pd%n_dos_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(siteproj, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(normalizationfactor, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(xmid, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(xlo, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(xhi, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf_spectral = 0.0_r8
        buf_spectral_smeared = 0.0_r8
        buf_sigmaIm = 0.0_r8
        buf_sigmaRe = 0.0_r8
        buf_taper = 0.0_r8
        siteproj = 0.0_r8

        do iq = 1, qp%n_irr_point

            if (mw%talk) then
                t1 = walltime()
                write (lo_iou, *) '... q-point '//tochar(iq)//' out of '//tochar(qp%n_irr_point)
                t0 = t1
                !write (*, '(a,i6,a,i6)') 'lineshape/phonondamping running qpt number ', iq, ' of ', qp%n_irr_point
            end if

            !if (mw%talk)
            ! get the self-energy
            !call se%generate(qp%ip(iq), qdir, dr%iq(iq), uc, fc, fct, fcf, ise, isf, qp, dr, opts, tmr, mw, mem, verbosity=-1)



            ! Check the scaling factor?
            ! evaluate the spectral function
            buf_spectral = 0.0_r8
            buf_spectral_smeared = 0.0_r8
            buf_sigmaIm = 0.0_r8
            buf_sigmaRe = 0.0_r8
            normalizationfactor = 0.0_r8
            xmid = 0.0_r8
            xlo = 0.0_r8
            xhi = 0.0_r8

            ! Evaluate the self-energy at this q
            call ise%evaluate(uc,qp%ip(iq)%r,dr%iq(iq),buf_sigmaRe,buf_sigmaIm,mw)

            if (mw%talk) then
                t1 = walltime()
                write (lo_iou, *) '     evaluate interpolation (', tochar(t1 - t0), 's)'
                t0 = t1
                !write (*, '(a,i6,a,i6)') 'lineshape/phonondamping running qpt number ', iq, ' of ', qp%n_irr_point
            end if


            do imode = 1, dr%n_mode
                if (mod(imode, mw%n) .ne. mw%r) cycle
                ! Get the spectral functions
                if (dr%iq(iq)%omega(imode) .gt. lo_freqtol) then
                    ! Taper self-energy
                    call taperfn_im(ise%omega, pd%dosmax, dr%iq(iq)%omega(imode), buf_taper)
                    buf_sigmaIm(:, imode) = buf_sigmaIm(:, imode)*buf_taper

                    ! First up would be to find a normalization factor per mode
                    call find_spectral_function_max_and_fwhm(dr%iq(iq)%omega(imode), dr%omega_max, ise%omega, buf_sigmaIm(:,imode), buf_sigmaRe(:,imode), xmid(imode),xlo(imode),xhi(imode))
                    call integrate_spectral_function(ise%omega, dr%iq(iq)%omega(imode), buf_sigmaIm(:,imode), buf_sigmaRe(:,imode), xmid(imode),xlo(imode),xhi(imode), &
                    1.0_r8, temperature, integraltol, f0, f1, f2)
                    normalizationfactor(imode)=1.0_r8/f0

                    call evaluate_spectral_function(ise%omega, buf_sigmaIm(:, imode), buf_sigmaRe(:, imode), dr%iq(iq)%omega(imode), buf_spectral(:, imode))
                    sigma = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, imode), dr%default_smearing(imode), smearing_prefactor)
                    buf_spectral_smeared(:, imode) = buf_spectral(:, imode)
                    call gaussian_smear_spectral_function(ise%omega, sigma, buf_spectral_smeared(:, imode))
                    ! Normalize both.
                    f0 = lo_trapezoid_integration(ise%omega, buf_spectral(:, imode))
                    buf_spectral(:, imode) = buf_spectral(:, imode)/f0
                    f0 = lo_trapezoid_integration(ise%omega, buf_spectral_smeared(:, imode))
                    buf_spectral_smeared(:, imode) = buf_spectral_smeared(:, imode)/f0

                    ! Accumulate DOS?
                    do iatom = 1, uc%na
                        cv0 = dr%iq(iq)%egv((iatom - 1)*3 + 1:iatom*3, imode)
                        siteproj(iatom) = abs(dot_product(cv0, conjg(cv0)))
                    end do

                    ! Add together in the right place
                    pd%pdos_mode(:, imode) = pd%pdos_mode(:, imode) + buf_spectral_smeared(:, imode)*qp%ip(iq)%integration_weight
                    do iatom = 1, uc%na
                        pd%pdos_site(:, iatom) = pd%pdos_site(:, iatom) + buf_spectral_smeared(:, imode)*siteproj(iatom)*qp%ip(iq)%integration_weight
                    end do

                else
                    buf_spectral(:, imode) = 0.0_r8
                    buf_spectral_smeared(:, imode) = 0.0_r8
                end if
            end do

            if (mw%talk) then
                t1 = walltime()
                write (lo_iou, *) '     integrals and normalization (', tochar(t1 - t0), 's)'
                t0 = t1
                !write (*, '(a,i6,a,i6)') 'lineshape/phonondamping running qpt number ', iq, ' of ', qp%n_irr_point
            end if

            call mw%allreduce('sum', buf_spectral)
            call mw%allreduce('sum', buf_spectral_smeared)
            call mw%allreduce('sum', buf_sigmaIm)
            call mw%allreduce('sum', buf_sigmaRe)
            call mw%allreduce('sum', normalizationfactor)
            call mw%allreduce('sum', xmid)
            call mw%allreduce('sum', xlo)
            call mw%allreduce('sum', xhi)

            ! Accumulate thermal transport data, can get it for free here anyway.
            call tc%accumulate(iq, qp, dr, fc, uc, buf_spectral, buf_spectral_smeared, buf_sigmaIm, buf_sigmaRe, normalizationfactor, xmid, xlo, xhi, mw, mem)

            if (mw%talk) then
                t1 = walltime()
                write (lo_iou, *) '     kappa accumulation (', tochar(t1 - t0), 's)'
                t0 = t1
                !write (*, '(a,i6,a,i6)') 'lineshape/phonondamping running qpt number ', iq, ' of ', qp%n_irr_point
            end if


            ! ! Report
            ! if (verbosity .gt. 0) then
            !     write (*, *) '    did', iq, 'out of', qp%n_irr_point
            ! end if
        end do

        ! Sync across ranks
        call mw%allreduce('sum', pd%pdos_site)
        call mw%allreduce('sum', pd%pdos_mode)

        ! Make sure to keep the frequency axis
        !pd%omega = se%energy_axis
        !pd%dosmax = maxval(se%energy_axis)

        call tc%report(qp, mw)

        ! Cleanup
        call mem%deallocate(buf_spectral, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_spectral_smeared, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_sigmaIm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_sigmaRe, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_taper, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(siteproj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(normalizationfactor, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(xmid, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(xlo, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(xhi, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    end block specfun


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
        !if ( verbosity .gt. 0 ) write (*, *) '... raw normalization (should be 1):', f0/dr%n_mode
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