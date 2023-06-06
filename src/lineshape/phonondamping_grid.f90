submodule(phonondamping) phonondamping_grid
!! Calculate self-energies on a mesh
implicit none
contains

!> Get the self-energy on a closed grid.
module subroutine get_selfenergy_on_closed_grid(sf, tc, pd, qp, dr, uc, fc, fct, fcf, ise, opts, tmr, mw, mem, verbosity)
    !> spectral functions on grid
    type(lo_spectralfunction_helper), intent(out) :: sf
    !> thermal conductivity
    type(lo_thermal_conductivity), intent(out) :: tc
    !> phonon dos
    type(lo_phonon_dos), intent(out) :: pd
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(inout) :: dr
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

    type(lo_spectralfunction_helper) :: isf
    integer :: solrnk

    ! Initialize helper
    init: block
        ! Read/write rank
        solrnk = mw%n - 1

        ! Initialize the helper (only on one rank, to save a little space)
        if (mw%r .eq. solrnk) then
            sf%n_irr_qpoint = qp%n_irr_point
            sf%n_energy = opts%nf
            sf%max_energy = opts%maxf
            sf%temperature = opts%temperature
            allocate (sf%qv(3, sf%n_irr_qpoint))
            allocate (sf%spectralfunction(sf%n_energy, dr%n_mode, sf%n_irr_qpoint))
            sf%qv = 0.0_r8
            sf%spectralfunction = 0.0_r8
        end if

        ! Obtain spectral functions. This is intended as a self-consistent
        ! procedure, which is why you need spectral functions to get spectral functions.
        ! For the first iteration we don't bother, just use a normal integration type.
        select case (opts%integrationtype)
            !case(4)
            ! Use the interpolation guy
            ! Smear the spectral functions?
        case (5)
            ! Read from file
            call isf%read_from_hdf5('infile.grid_spectral_function.hdf5', qp, opts%maxf, opts%temperature, opts%nf, mw, verbosity)
        end select

        ! Initialize the thermal conductivity helper.
        call tc%init(dr, opts%nf, opts%maxf, opts%temperature, mw)
    end block init

    ! Initialize the phonon DOS, cost nothing so might as well.
    initdos: block
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
    end block initdos

    call tmr%tock('initialize grid calculation')

    ! To integrate I need to smear the spectral functions appropriately
    select case (opts%integrationtype)
    case (4:5)
        !call isf%smear(qp,dr,opts%maxf,opts%sigma,mw,mem)
    case default
        ! do nothing
    end select

    ! Calculate self-energies / spectral functions
    specfun: block
        complex(r8), dimension(3) :: cv0
        real(r8), parameter :: largedistance = 1.0E10_r8*lo_A_to_bohr
        real(r8), dimension(:, :), allocatable :: buf_spectral, buf_spectral_smeared
        real(r8), dimension(:, :), allocatable :: buf_sigmaIm, buf_sigmaRe
        real(r8), dimension(:), allocatable :: buf_taper, siteproj
        real(r8), dimension(3), parameter :: qdir = [1.0_r8, 0.0_r8, 0.0_r8]
        real(r8) :: f0, sigma
        type(lo_phonon_selfenergy) :: se
        integer :: iq, imode, iatom

        ! A little space for temporary buffers
        call mem%allocate(buf_spectral, [opts%nf, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_spectral_smeared, [opts%nf, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_sigmaIm, [opts%nf, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_sigmaRe, [opts%nf, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_taper, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(siteproj, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf_spectral = 0.0_r8
        buf_spectral_smeared = 0.0_r8
        buf_sigmaIm = 0.0_r8
        buf_sigmaRe = 0.0_r8
        buf_taper = 0.0_r8
        siteproj = 0.0_r8

        do iq = 1, qp%n_irr_point
            if (mw%talk) write (*, '(a,i6,a,i6)') 'lineshape/phonondamping running qpt number ', iq, ' of ', qp%n_irr_point
            ! get the self-energy
            call se%generate(qp%ip(iq), qdir, dr%iq(iq), uc, fc, fct, fcf, ise, isf, qp, dr, opts, tmr, mw, mem, verbosity=-1)

            ! evaluate the spectral function
            buf_spectral = 0.0_r8
            buf_spectral_smeared = 0.0_r8
            buf_sigmaIm = 0.0_r8
            buf_sigmaRe = 0.0_r8
            do imode = 1, dr%n_mode
                if (mod(imode, mw%n) .ne. mw%r) cycle
                ! Get the spectral functions
                if (dr%iq(iq)%omega(imode) .gt. lo_freqtol) then
                    buf_sigmaIm(:, imode) = se%im_3ph(:, imode) + se%im_iso(:, imode)
                    buf_sigmaRe(:, imode) = se%re_3ph(:, imode) + se%re_4ph(:, imode)
                    call taperfn_im(se%energy_axis, opts%maxf, dr%iq(iq)%omega(imode), buf_taper)
                    buf_sigmaIm(:, imode) = buf_sigmaIm(:, imode)*buf_taper
                    call evaluate_spectral_function(se%energy_axis, buf_sigmaIm(:, imode), buf_sigmaRe(:, imode), dr%iq(iq)%omega(imode), buf_spectral(:, imode))

                    sigma = qp%adaptive_sigma(qp%ip(iq)%radius, dr%iq(iq)%vel(:, imode), dr%default_smearing(imode), opts%sigma)
                    buf_spectral_smeared(:, imode) = buf_spectral(:, imode)
                    call gaussian_smear_spectral_function(se%energy_axis, sigma, buf_spectral_smeared(:, imode))
                    ! Normalize both.
                    f0 = lo_trapezoid_integration(se%energy_axis, buf_spectral(:, imode))
                    buf_spectral(:, imode) = buf_spectral(:, imode)/f0
                    f0 = lo_trapezoid_integration(se%energy_axis, buf_spectral_smeared(:, imode))
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
            call mw%allreduce('sum', buf_spectral)
            call mw%allreduce('sum', buf_spectral_smeared)
            call mw%allreduce('sum', buf_sigmaIm)
            call mw%allreduce('sum', buf_sigmaRe)

            ! The spectral functions can be used for many things, so better use it
            ! while we have it. First would be to just store it away:
            if (mw%r .eq. solrnk) then
                ! Store q-coordinate
                sf%qv(:, iq) = qp%ip(iq)%r
                ! Store spectral function
                sf%spectralfunction(:, :, iq) = buf_spectral
            end if

            call tmr%tock('evaluate spectral function')

            ! Accumulate thermal transport data, can get it for free here anyway.
            call tmr%tick()
            call tc%accumulate(iq, qp, dr, fc, uc, buf_spectral, buf_spectral_smeared, buf_sigmaIm, buf_sigmaRe, se%scalingfactor, se%xmid, se%xlo, se%xhi, mw, mem)
            call tmr%tock('evaluate thermal conductivity')

            ! Report
            if (verbosity .gt. 0) then
                write (*, *) '    did', iq, 'out of', qp%n_irr_point
            end if
        end do

        ! Make sure to keep the frequency axis
        pd%omega = se%energy_axis
        pd%dosmax = maxval(se%energy_axis)

        call tc%report(qp, mw)

        ! Cleanup
        call mem%deallocate(buf_spectral, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_spectral_smeared, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_sigmaIm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_sigmaRe, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_taper, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(siteproj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    end block specfun

    ! Now fix the DOS so that it looks nice.
    makenice: block
        integer :: i, j, k
        real(r8), dimension(:, :), allocatable :: sitebuf
        real(r8) :: f0, f1, dostol

        ! Sync across ranks
        call mw%allreduce('sum', pd%pdos_site)
        call mw%allreduce('sum', pd%pdos_mode)

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

    ! Do some other things maybe? later point.
end subroutine

! !> calculate things as a density of states
! subroutine get_selfenergy_on_grid(pd,qpd,drd,uc,fc,fct,fcf,ise,qp,dr,opts,sigmaRe,sigmaIm,tmr,mw,mem,verbosity)
!     !> phonon density of states
!     type(lo_phonon_dos), intent(out) :: pd
!     !> q-point mesh
!     class(lo_qpoint_mesh), intent(in) :: qpd
!     !> harmonic properties on this mesh
!     type(lo_phonon_dispersions), intent(inout) :: drd
!     !> crystal structure
!     type(lo_crystalstructure), intent(inout) :: uc
!     !> second order force constant
!     type(lo_forceconstant_secondorder), intent(inout) :: fc
!     !> third order force constant
!     type(lo_forceconstant_thirdorder), intent(in) :: fct
!     !> fourth order force constant
!     type(lo_forceconstant_fourthorder), intent(in) :: fcf
!     !> self-energy
!     type(lo_interpolated_selfenergy), intent(in) :: ise
!     !> q-point mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> harmonic properties on this mesh
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> all settings
!     type(lo_opts), intent(in) :: opts
!     !> real part of self-energy
!     real(r8), dimension(:,:,:), intent(inout) :: sigmaRe
!     !> imaginary part of self-energy
!     real(r8), dimension(:,:,:), intent(inout) :: sigmaIm
!     !> timer
!     type(lo_timer), intent(inout) :: tmr
!     !> mpi helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> memory tracker
!     type(lo_mem_helper), intent(inout) :: mem
!     !> talk a lot?
!     integer, intent(in) :: verbosity
!
!     real(r8), dimension(3), parameter :: qdir=[1.0_r8,0.0_r8,0.0_r8]
!     type(lo_spectralfunction_helper) :: isf
!     type(lo_phonon_selfenergy) :: se
!     integer :: iq
!
!     sigmaRe=0.0_r8
!     sigmaIm=0.0_r8
!     do iq=1,qpd%n_irr_point
!         ! get the self-energy
!         call se%generate(qpd%ip(iq),qdir,drd%iq(iq),uc,fc,fct,fcf,ise,isf,qp,dr,opts,tmr,mw,mem,verbosity=-1)
!         ! Store it
!         sigmaRe(:,:,iq)=se%re_3ph+se%re_4ph
!         sigmaIm(:,:,iq)=se%im_3ph+se%im_iso
!         ! Report
!         if ( verbosity .gt. 0 ) then
!             write(*,*) 'did',iq,'out of',qpd%n_irr_point
!         endif
!     enddo
! end subroutine

!> lineshape on just a list of q-vectors
module subroutine get_selfenergy_on_points(qvec, wp, qp, dr, uc, fc, fct, fcf, ise, opts, sigmaRe, sigmaIm, tmr, mw, mem, verbosity)
    !> q-vectors
    real(r8), dimension(:, :), intent(in) :: qvec
    !> dispersions on points
    type(lo_phonon_dispersions_qpoint), dimension(:), intent(inout) :: wp
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(inout) :: dr
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
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> real part of self-energy
    real(r8), dimension(:, :, :), intent(inout) :: sigmaRe
    !> imaginary part of self-energy
    real(r8), dimension(:, :, :), intent(inout) :: sigmaIm
    !> timer
    type(lo_timer), intent(inout) :: tmr
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(3), parameter :: qdir = [1.0_r8, 0.0_r8, 0.0_r8]

    type(lo_qpoint) :: qpoint
    type(lo_spectralfunction_helper) :: isf
    type(lo_phonon_selfenergy) :: se
    integer :: iq

    sigmaRe = 0.0_r8
    sigmaIm = 0.0_r8
    do iq = 1, size(qvec, 2)
        ! get harmonic things at this q-point
        qpoint%r = qvec(:, iq)
        call lo_get_small_group_of_qpoint(qpoint, uc)
        call wp(iq)%generate(fc, uc, mem, qpoint, qdirection=qdir)

        ! get the self-energy
        call se%generate(qpoint, qdir, wp(iq), uc, fc, fct, fcf, ise, isf, qp, dr, opts, tmr, mw, mem, verbosity=-1)
        ! Store it
        sigmaRe(:, :, iq) = se%re_3ph + se%re_4ph
        sigmaIm(:, :, iq) = se%im_3ph + se%im_iso
        ! Report
        if (verbosity .gt. 0) then
            write (*, *) 'did', iq, 'out of', size(qvec, 2)
        end if
    end do
end subroutine

end submodule
