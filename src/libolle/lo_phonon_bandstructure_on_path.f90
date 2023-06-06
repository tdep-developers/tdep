module lo_phonon_bandstructure_on_path
!!
!! Deal with phonon dispersions along a path in the BZ
!!
    use konstanter, only: r8, i8, lo_iou, lo_huge, lo_hugeint, lo_gnuplotterminal, lo_tol, lo_sqtol, lo_freqtol, &
                          lo_frequency_hartree_to_icm, lo_frequency_hartree_to_meV, lo_frequency_hartree_to_thz, &
                          lo_groupvel_hartreebohr_to_ms, lo_exitcode_param, lo_frequency_Hartree_to_Hz, lo_bohr_to_A, &
                          lo_exitcode_symmetry, lo_gitbranch, lo_gitrevision
    use gottochblandat, only: open_file, walltime, lo_chop, tochar, lo_trueNtimes, lo_progressbar_init, lo_progressbar
    use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
    use lo_memtracker, only: lo_mem_helper
    use hdf5_wrappers, only: lo_hdf5_helper
    use type_crystalstructure, only: lo_crystalstructure
    use type_phonon_dispersions, only: lo_phonon_dispersions_qpoint
    use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
    use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
    use type_qpointmesh, only: lo_bandstructure

    implicit none
    private
    public :: lo_phonon_bandstructure

!> metadata for an active q-point
    type lo_active_qpoint
        !> where on the path is this q-point
        integer :: ind_path = -lo_hugeint
        !> how many irreducible components in the Raman tensor
        integer :: nx_raman = -lo_hugeint
        !> how many irreducible components in the IR tensor
        integer :: nx_IR = -lo_hugeint
        !> realspace raman coefficient matrix
        real(r8), dimension(:, :), allocatable :: coeff_Raman
        !> realspace IR coefficient matrix
        real(r8), dimension(:, :), allocatable :: coeff_IR
        !> raman coefficient per mode
        real(r8), dimension(:, :, :), allocatable :: coeff_Raman_mode
        !> IR coefficient per mode
        real(r8), dimension(:, :, :), allocatable :: coeff_IR_mode
    end type

!> Phonon dispersions along a path in the Brillouin zone
    type, extends(lo_bandstructure) :: lo_phonon_bandstructure
        !> number of phonon modes
        integer :: n_mode = -lo_hugeint
        !> points that hold dispersion data
        type(lo_phonon_dispersions_qpoint), dimension(:), allocatable :: p

        !> number of active (IR/Raman) q-points
        integer :: n_active_qpoint = -lo_hugeint
        !> active q-points
        type(lo_active_qpoint), dimension(:), allocatable :: active_qpoint

        !> energy axis for spectral function
        real(r8), dimension(:), allocatable :: energy_axis
        !> spectral function (q,energy)
        real(r8), dimension(:, :), allocatable :: spectral_function
        !> spectral function with thermal prefactor (q,energy)
        real(r8), dimension(:, :), allocatable :: spectral_function_with_prefactor
        !> real part of the self-energy (q,energy,mode)
        real(r8), dimension(:, :, :), allocatable :: selfenergy_real
        !> imaginary part of the self-energy (q,energy,mode)
        real(r8), dimension(:, :, :), allocatable :: selfenergy_imag
    contains
        !> calculate the bandstructure
        procedure :: generate
        !> allgather over mpi
        procedure :: allgather
        !> sort the bands for nicer plots
        procedure :: sort_bands
        !> write files
        procedure :: write_dispersive_property
        !> write files to hdf5
        procedure :: write_to_hdf5
        !> write spectral_function to hdf5
        procedure :: write_spectral_function_to_hdf5
        !> size in memory, in bytes
        procedure :: size_in_mem => bs_size_in_mem
        !> destroy
        procedure :: destroy
    end type

contains

!> Get the dispersions along the path
    subroutine generate(bs, p, fc, timereversal, mw, mem, verbosity, npts, fct, readpathfromfile)
        !> bandstructure
        class(lo_phonon_bandstructure), intent(out) :: bs
        !> crystal structure
        type(lo_crystalstructure), intent(inout) :: p
        !> forceconstant
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        !> Enforce time-reversal symmetry?
        logical, intent(in) :: timereversal
        !> MPI helper
        type(lo_mpi_helper), intent(inout) :: mw
        !> memory tracker
        type(lo_mem_helper), intent(inout) :: mem
        !> how much to talk
        integer, intent(in) :: verbosity
        !> How many points between each high symmetry point?
        integer, intent(in), optional :: npts
        !> thirdorder forceconstant
        type(lo_forceconstant_thirdorder), intent(in), optional :: fct
        !> should the path be read from file?
        logical, intent(in), optional :: readpathfromfile

        real(r8) :: timer
        logical :: readpath, gruneisen

        ! Start timer
        timer = walltime()

        ! set all the options, fetch the path and so on.
        initandheuristics: block
            integer :: i

            if (verbosity .gt. 0) then
                write (lo_iou, *) ''
                write (lo_iou, *) 'Calculating phonon dispersions on a q-point path'
            end if

            ! Make sure all the symmetry stuff is there!
            if (p%info%havewedge .eqv. .false.) then
                call p%classify('wedge', timereversal=timereversal)
            end if
            if (p%info%decidedtimereversal .eqv. .false.) then
                call lo_stop_gracefully(['time-reversal symmetry needs to be decided for dispersions on a path'], &
                                        lo_exitcode_param, __FILE__, __LINE__, mw%comm)
            end if
            if (p%sym%timereversal .neqv. timereversal) then
                call lo_stop_gracefully(['Conflicting instructions regarding time-reversal symmetry for dispersions on path'], &
                                        lo_exitcode_param, __FILE__, __LINE__, mw%comm)
            end if

            ! calculate Grunesien parameters?
            if (present(fct)) then
                gruneisen = .true.
            else
                gruneisen = .false.
            end if
            ! read or generate path
            if (present(readpathfromfile)) then
                readpath = readpathfromfile
            else
                readpath = .false.
            end if
            ! Number of points
            if (present(npts)) then
                bs%n_point_per_path = npts
            else
                bs%n_point_per_path = 100
            end if
            ! Coordinates of the path
            if (readpath) then
                call bs%read_path_from_file(p, mw, verbosity)
            else
                call bs%standardpath(p, mw, verbosity)
            end if

            ! Now the heuristics are done, set some more constants and make space
            bs%n_mode = p%na*3
            allocate (bs%p(bs%n_point))
            do i = 1, bs%n_point
                allocate (bs%p(i)%omega(bs%n_mode))
                allocate (bs%p(i)%egv(bs%n_mode, bs%n_mode))
                allocate (bs%p(i)%vel(3, bs%n_mode))
                allocate (bs%p(i)%degeneracy(bs%n_mode))
                allocate (bs%p(i)%degenmode(bs%n_mode, bs%n_mode))
                if (gruneisen) allocate (bs%p(i)%gruneisen(bs%n_mode))
                bs%p(i)%omega = 0.0_r8
                bs%p(i)%egv = 0.0_r8
                bs%p(i)%vel = 0.0_r8
                bs%p(i)%degeneracy = 0
                bs%p(i)%degenmode = 0
                if (gruneisen) bs%p(i)%gruneisen = 0.0_r8
            end do
        end block initandheuristics

        ! calculate the actual dispersions
        getdispersions: block
            complex(r8), dimension(:, :, :), allocatable :: Dq
            complex(r8), dimension(:, :), allocatable :: D
            integer, dimension(:, :), allocatable :: di
            real(r8), dimension(3) :: qdir
            real(r8) :: f0
            integer :: q, lq, j, nq

            call mem%allocate(D, [bs%n_mode, bs%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(Dq, [bs%n_mode, bs%n_mode, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(di, [bs%n_mode, bs%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            D = 0.0_r8
            Dq = 0.0_r8
            di = 0

            ! Count number of q-points per rank for some reason.
            nq = 0
            do q = 1, bs%n_point
                if (mod(q, mw%n) .eq. mw%r) nq = nq + 1
            end do

            if (verbosity .gt. 0) call lo_progressbar_init()
            lq = 0
            do q = 1, bs%n_point
                if (mod(q, mw%n) .ne. mw%r) cycle
                lq = lq + 1
                j = bs%q(q)%path
                ! get the direction for lo-to splitting
                qdir = bs%segment(j)%r2 - bs%segment(j)%r1
                f0 = norm2(qdir)
                if (f0 .gt. lo_tol) then
                    qdir = qdir/f0
                else
                    call random_number(qdir)
                    qdir = qdir/norm2(qdir)
                end if
                ! the actual dispersions
                call bs%p(q)%generate(fc, p, mem, bs%q(q), qdirection=qdir)
                ! perhaps gruneisen parameters as well?
                if (gruneisen) then
                    call fct%mode_gruneisen_parameter(p, bs%q(q), bs%p(q)%omega, bs%p(q)%egv, bs%p(q)%gruneisen)
                end if
                ! and report?
                if (verbosity .gt. 0) then
                    if (lo_trueNtimes(q, 127, nq)) call lo_progressbar('... calculating frequencies', lq, nq, walltime() - timer)
                end if
            end do
            ! And finally communicate it
            call bs%allgather(mw, mem)

            ! finish reporting
            if (verbosity .gt. 0) call lo_progressbar('... calculating frequencies', nq, nq, walltime() - timer)

            ! cleanup
            call mem%deallocate(D, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(Dq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block getdispersions

        if (verbosity .gt. 0) write (lo_iou, *) 'Calculated dispersion on path in ', tochar(walltime() - timer), 's'
    end subroutine

!> Distribute a phonon dispersion object over MPI, all-to-all
    subroutine allgather(bs, mw, mem)
        !> phonon dispersions
        class(lo_phonon_bandstructure), intent(inout) :: bs
        !> mpi helper
        type(lo_mpi_helper), intent(inout) :: mw
        !> memory tracker
        type(lo_mem_helper), intent(inout) :: mem

        character, dimension(:), allocatable :: sbuf, rbuf
        integer :: i, j, ctr, ctr_tot, pos

        ! Start counting bytes per rank that should be sent.
        ctr = 0
        do i = 1, bs%n_point
            if (mod(i, mw%n) .ne. mw%r) cycle
            ! first want I pack the q-index
            ctr = ctr + storage_size(i)/8
            ! then the rest of the q-point
            ctr = ctr + bs%p(i)%size_packed()
        end do

        ! Get size per rank and offset
        call mw%size_and_offset(ctr, ctr_tot)

        ! Space for communication buffers
        if (ctr .gt. 0) then
            call mem%allocate(sbuf, ctr, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
        else
            call mem%allocate(sbuf, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
        end if
        call mem%allocate(rbuf, ctr_tot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Pack things into send buffer
        pos = 0
        do i = 1, bs%n_point
            if (mod(i, mw%n) .ne. mw%r) cycle
            ! pack the index
            call mw%pack(i, sbuf, pos)
            ! pack the point
            call bs%p(i)%pack_to_buf(sbuf, pos, mw)
        end do
        ! Now I can commence operation actual operation
        call mw%allgatherv(sbuf, rbuf, mw%ctr_size_per_rank, mw%ctr_offset_per_rank, __FILE__, __LINE__)
        ! Clear the send buffer
        call mem%deallocate(sbuf, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
        ! Start unpacking
        pos = 0
        do j = 1, bs%n_point
            ! fetch which q-index is stored here
            call mw%unpack(i, rbuf, pos)
            ! unpack the point
            call bs%p(i)%unpack_from_buf(rbuf, pos, mw)
        end do
        ! Clear the recieve buffer
        call mem%deallocate(rbuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end subroutine

!> Rearrange the bands such that they are continous. Does not quite work.
    subroutine sort_bands(bs, mem, verbosity)
        !> bandstructure
        class(lo_phonon_bandstructure), intent(inout) :: bs
        !> memory tracker
        type(lo_mem_helper), intent(inout) :: mem
        !> talk?
        integer, intent(in) :: verbosity

        complex(r8), dimension(:, :, :), allocatable :: egvbuf
        real(r8), dimension(:, :, :), allocatable :: velbuf
        real(r8), dimension(:, :), allocatable :: ombuf
        real(r8) :: f0, f1, f2, dsum, deltaom
        integer, dimension(:, :, :), allocatable :: bitot
        integer, dimension(:, :), allocatable :: bi
        integer, dimension(:), allocatable :: iperm
        integer :: i, j, l, path, b1, b2, bcurr
        logical, dimension(:, :), allocatable :: matched

        if (verbosity .gt. 0) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'Attemting to sort bands into something continous'
        end if

        call mem%allocate(bi, [bs%n_mode, bs%n_point_per_path], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(bitot, [bs%n_mode, bs%n_point_per_path, bs%n_path], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ombuf, [bs%n_mode, bs%n_point_per_path], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(velbuf, [3, bs%n_mode, bs%n_point_per_path], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(egvbuf, [bs%n_mode, bs%n_mode, bs%n_point_per_path], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(iperm, bs%n_mode, &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(matched, [bs%n_mode, bs%n_point_per_path], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        bi = 0
        bitot = 0
        ombuf = 0.0_r8
        velbuf = 0.0_r8
        egvbuf = 0.0_r8
        iperm = 0
        matched = .false.

        ! First figure out the global deltaomega, as in what is the largest change between two qpoints for a single band.
        deltaom = 0.0_r8
        do path = 1, bs%n_path
        do i = 1, bs%n_point_per_path - 1
            j = (path - 1)*bs%n_point_per_path + i
            do b1 = 1, bs%n_mode
                deltaom = max(deltaom, abs(bs%p(j)%omega(b1) - bs%p(j + 1)%omega(b1)))
            end do
        end do
        end do
        deltaom = deltaom + lo_freqtol

        bitot = 0
        ! Rearrange it per segment
        do path = 1, bs%n_path

            ! Fetch the frequencies and velocities
            do i = 1, bs%n_point_per_path
                j = (path - 1)*bs%n_point_per_path + i
                ombuf(:, i) = bs%p(j)%omega
                velbuf(:, :, i) = bs%p(j)%vel
                egvbuf(:, :, i) = bs%p(j)%egv
            end do

            ! Try to figure out the actual connectivity
            matched = .false.
            matched(:, 1) = .true.
            bi = 0
            do b1 = 1, bs%n_mode
                ! start with something
                bi(b1, 1) = b1
                bcurr = b1
                do i = 2, bs%n_point_per_path
                    ! figure out what the best match is
                    l = 0
                    dsum = lo_huge
                    do b2 = 1, bs%n_mode
                        if (matched(b2, i)) cycle                              ! only match once!
                        f0 = abs(ombuf(bcurr, i - 1) - ombuf(b2, i))                  ! delta in frequency
                        !if ( f0 .gt. 10*deltaom ) cycle                        ! should be reasonable close in frequency
                        f1 = norm2(velbuf(:, bcurr, i - 1) - velbuf(:, b2, i))          ! delta in group velocity
                        f2 = abs(dot_product(egvbuf(:, bcurr, i - 1), egvbuf(:, b2, i)))
                        f2 = abs(1.0_r8 - f2)
                        f1 = (f1*10)**2*f0**2*f2
                        if (f1 .lt. dsum) then
                            l = b2
                            dsum = f1
                        end if
                    end do
                    if (l .eq. 0) then
                        call lo_stop_gracefully(['I am not a smart man.'], lo_exitcode_symmetry, __FILE__, __LINE__)
                    end if
                    ! Now I have used that band, never use it again!
                    matched(l, i) = .true.

                    ! destroy the old one so that it does not get matched again
                    bcurr = l
                    bi(b1, i) = bcurr
                end do
            end do

            ! Store the alignment of this path
            bitot(:, :, path) = bi
        end do

        ! rearrange again so the it becomes continous between segments, sort of
        do path = 2, bs%n_path
            iperm = bitot(:, bs%n_point_per_path, path - 1)
            do i = 1, bs%n_point_per_path
                bitot(:, i, path) = bitot(iperm, i, path)
            end do
        end do

        ! Now rearrange all the things
        do path = 1, bs%n_path
        do j = 1, bs%n_point_per_path
            i = (path - 1)*bs%n_point_per_path + j
            iperm = bitot(:, j, path)
            if (allocated(bs%p(i)%omega)) bs%p(i)%omega = bs%p(i)%omega(iperm)
            if (allocated(bs%p(i)%vel)) bs%p(i)%vel = bs%p(i)%vel(:, iperm)
            if (allocated(bs%p(i)%egv)) bs%p(i)%egv = bs%p(i)%egv(:, iperm)
            if (allocated(bs%p(i)%gruneisen)) bs%p(i)%gruneisen = bs%p(i)%gruneisen(iperm)
            if (allocated(bs%p(i)%linewidth)) bs%p(i)%linewidth = bs%p(i)%linewidth(iperm)
            if (allocated(bs%p(i)%shift3)) bs%p(i)%shift3 = bs%p(i)%shift3(iperm)
            if (allocated(bs%p(i)%shift4)) bs%p(i)%shift4 = bs%p(i)%shift4(iperm)
            if (allocated(bs%p(i)%p_plus)) bs%p(i)%p_plus = bs%p(i)%p_plus(iperm)
            if (allocated(bs%p(i)%p_minus)) bs%p(i)%p_minus = bs%p(i)%p_minus(iperm)
            if (allocated(bs%p(i)%p_iso)) bs%p(i)%p_iso = bs%p(i)%p_iso(iperm)
            if (allocated(bs%p(i)%qs)) bs%p(i)%qs = bs%p(i)%qs(iperm)
            if (allocated(bs%p(i)%F0)) bs%p(i)%F0 = bs%p(i)%F0(:, iperm)
            if (allocated(bs%p(i)%Fn)) bs%p(i)%Fn = bs%p(i)%Fn(:, iperm)
            if (allocated(bs%p(i)%kappa)) bs%p(i)%kappa = bs%p(i)%kappa(:, :, iperm)
            if (allocated(bs%p(i)%F)) bs%p(i)%F = bs%p(i)%F(iperm)
            if (allocated(bs%p(i)%S)) bs%p(i)%S = bs%p(i)%S(iperm)
            if (allocated(bs%p(i)%deltaF3)) bs%p(i)%deltaF3 = bs%p(i)%deltaF3(iperm)
            if (allocated(bs%p(i)%deltaF4)) bs%p(i)%deltaF4 = bs%p(i)%deltaF4(iperm)
            if (allocated(bs%p(i)%deltaS3)) bs%p(i)%deltaS3 = bs%p(i)%deltaS3(iperm)
            if (allocated(bs%p(i)%deltaS4)) bs%p(i)%deltaS4 = bs%p(i)%deltaS4(iperm)
            if (allocated(bs%p(i)%electronphononphasespace)) &
                bs%p(i)%electronphononphasespace = bs%p(i)%electronphononphasespace(iperm)
            if (allocated(bs%p(i)%omega_elastic)) bs%p(i)%omega_elastic = bs%p(i)%omega_elastic(iperm)
            if (allocated(bs%p(i)%threephononphasespace)) bs%p(i)%threephononphasespace = bs%p(i)%threephononphasespace(iperm)
            if (allocated(bs%p(i)%mfp)) bs%p(i)%mfp = bs%p(i)%mfp(:, iperm)
            if (allocated(bs%p(i)%scalar_mfp)) bs%p(i)%scalar_mfp = bs%p(i)%scalar_mfp(iperm)
            if (allocated(bs%p(i)%thermal_prefactor)) bs%p(i)%thermal_prefactor = bs%p(i)%thermal_prefactor(iperm)
            if (allocated(bs%p(i)%pyroelectric_amplitude)) bs%p(i)%pyroelectric_amplitude = bs%p(i)%pyroelectric_amplitude(iperm)
            if (allocated(bs%p(i)%degeneracy)) bs%p(i)%degeneracy = bs%p(i)%degeneracy(iperm)
            if (allocated(bs%p(i)%degenmode)) bs%p(i)%degenmode = bs%p(i)%degenmode(:, iperm)
        end do
        end do

        ! And cleanup
        call mem%deallocate(bi, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(bitot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ombuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(velbuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(egvbuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(iperm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(matched, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end subroutine

!> dump a dispersion plot
    subroutine write_to_hdf5(bs, uc, enhet, filename, mem, hdftag)
        !> the dispersion relations along the path
        class(lo_phonon_bandstructure), intent(in) :: bs
        !> crystal structure
        type(lo_crystalstructure), intent(in) :: uc
        !> unit
        character(len=3), intent(in) :: enhet
        !> filename
        character(len=*), intent(in) :: filename
        !> memory tracker
        type(lo_mem_helper), intent(inout) :: mem
        !> optionally, write to an already open hdf5 file.
        integer(i8), intent(in), optional :: hdftag

        type(lo_hdf5_helper) :: h5
        character(len=10), dimension(:), allocatable :: dumlbl
        character(len=1000) :: lblstr, unitstr
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(:, :), allocatable :: d2
        real(r8), dimension(:, :, :), allocatable :: d3
        real(r8) :: unitfactor, f0
        integer :: i, j, k

        select case (enhet)
        case ('thz')
            unitfactor = lo_frequency_hartree_to_THz
            unitstr = 'THz'
        case ('mev')
            unitfactor = lo_frequency_hartree_to_meV
            unitstr = 'meV'
        case ('icm')
            unitfactor = lo_frequency_hartree_to_icm
            unitstr = 'cm^-^1'
        case default
            call lo_stop_gracefully(['Unknown unit'], lo_exitcode_param, __FILE__, __LINE__)
        end select

        ! create the file, or use the provided tag
        if (present(hdftag)) then
            ! Write to the tag that was provided
            h5%file_id = hdftag
        else
            ! Create a new file.
            call h5%init(__FILE__, __LINE__)
            call h5%open_file('write', trim(filename))
        end if

        ! Dump some metadata
        call h5%open_group('write', 'metadata')
        call h5%store_attribute(lo_gitbranch, h5%group_id, 'git_branch')
        call h5%store_attribute(lo_gitrevision, h5%group_id, 'git_revision')
        call h5%close_group()

        ! Dump the scalar x-axis
        call h5%store_data(bs%q_axis/lo_bohr_to_A, h5%file_id, 'q_values', enhet='A^-1', dimensions='q-point')
        ! The ticks for the x-axis
        call h5%store_data(bs%q_axis_ticks/lo_bohr_to_A, h5%file_id, 'q_ticks', enhet='A^-1')
        ! The labels for the x-ticks, really annoying
        allocate (dumlbl(size(bs%q_axis_tick_labels, 1)))
        do i = 1, size(dumlbl, 1)
            k = 0
            dumlbl(i) = '         '
            do j = 1, len(dumlbl(i))
                if (bs%q_axis_tick_labels(i) (j:j) .ne. '|') then
                    k = k + 1
                    dumlbl(i) (k:k) = trim(bs%q_axis_tick_labels(i) (j:j))
                end if
            end do
            if (trim(dumlbl(i)) .eq. 'Γ') then
                dumlbl(i) = 'G'
            end if
        end do
        lblstr = ''
        do i = 1, size(dumlbl, 1)
            lblstr = trim(lblstr)//' '//trim(dumlbl(i))
        end do
        deallocate (dumlbl)
        call h5%store_attribute(trim(adjustl(lblstr)), h5%file_id, 'q_tick_labels')

        ! Time to dump the dispersive properties:
        if (allocated(bs%p(1)%omega)) then
            call mem%allocate(d2, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            d2 = 0.0_r8
            do i = 1, bs%n_point
                d2(:, i) = bs%p(i)%omega*unitfactor
            end do
            call h5%store_data(d2, h5%file_id, 'frequencies', enhet=trim(unitstr), dimensions='q-point,mode')
            call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        if (allocated(bs%p(1)%egv)) then
            call mem%allocate(d3, [bs%n_mode, bs%n_mode, bs%n_point], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            do i = 1, bs%n_point
                d3(:, :, i) = lo_chop(real(bs%p(i)%egv), lo_sqtol)
            end do
            call h5%store_data(d3, h5%file_id, 'eigenvectors_re', enhet='dimensionless', dimensions='q-point,mode,atomxyz')
            do i = 1, bs%n_point
                d3(:, :, i) = lo_chop(aimag(bs%p(i)%egv), lo_sqtol)
            end do
            call h5%store_data(d3, h5%file_id, 'eigenvectors_im', enhet='dimensionless', dimensions='q-point,mode,atomxyz')
            call mem%deallocate(d3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            ! projection of each frequency onto normal modes
            call mem%allocate(d3, [uc%na, bs%n_mode, bs%n_point], &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            d3 = 0.0_r8
            do i = 1, bs%n_point
                do j = 1, bs%n_mode
                do k = 1, uc%na
                    cv0 = bs%p(i)%egv((k - 1)*3 + 1:k*3, j)
                    d3(k, j, i) = abs(dot_product(conjg(cv0), cv0))**2
                end do
                end do
                ! make sure it's normalized correctly
                do j = 1, bs%n_mode
                    d3(:, j, i) = d3(:, j, i)/sum(d3(:, j, i))
                end do
            end do
            call h5%store_data(d3, h5%file_id, 'site_projection_per_mode', enhet='dimensionless', dimensions='q-point,mode,atom')
            call mem%deallocate(d3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        if (allocated(bs%p(1)%vel)) then
            call mem%allocate(d3, [3, bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            d3 = 0.0_r8
            do i = 1, bs%n_point
                d3(:, :, i) = lo_chop(bs%p(i)%vel, lo_sqtol)*lo_groupvel_hartreebohr_to_ms
            end do
            call h5%store_data(d3, h5%file_id, 'group_velocities', enhet='m/s', dimensions='q-point,mode,xyz')
            call mem%deallocate(d3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        if (allocated(bs%p(1)%vel) .and. allocated(bs%p(1)%linewidth)) then
            call mem%allocate(d2, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(d3, [3, bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            d2 = 0.0_r8
            d3 = 0.0_r8
            do i = 1, bs%n_point
            do j = 1, bs%n_mode
                f0 = bs%p(i)%linewidth(j)
                if (bs%p(i)%linewidth(j) .gt. lo_freqtol) then
                    f0 = f0*lo_frequency_hartree_to_Hz
                    d3(:, j, i) = bs%p(i)%vel(:, j)*lo_groupvel_hartreebohr_to_ms*0.5_r8/f0
                    d2(j, i) = norm2(d3(:, j, i))
                end if
            end do
            end do
            call h5%store_data(d3, h5%file_id, 'mean_free_path', enhet='m', dimensions='q-point,mode,xyz')
            call h5%store_data(d2, h5%file_id, 'scalar_mean_free_path', enhet='m', dimensions='q-point,mode')
            call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(d3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        if (allocated(bs%p(1)%gruneisen)) then
            call mem%allocate(d2, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            d2 = 0.0_r8
            do i = 1, bs%n_point
                d2(:, i) = bs%p(i)%gruneisen
            end do
            call h5%store_data(d2, h5%file_id, 'mode_gruneisen_parameters', enhet='dimensionless', dimensions='q-point,mode')
            call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        if (allocated(bs%p(1)%linewidth)) then
            call mem%allocate(d2, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            d2 = 0.0_r8
            do i = 1, bs%n_point
                d2(:, i) = bs%p(i)%linewidth*unitfactor
            end do
            call h5%store_data(d2, h5%file_id, 'linewidth', enhet=trim(unitstr), dimensions='q-point,mode')
            call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        if (allocated(bs%p(1)%shift3)) then
            call mem%allocate(d2, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            d2 = 0.0_r8
            do i = 1, bs%n_point
                d2(:, i) = bs%p(i)%shift3*unitfactor
            end do
            call h5%store_data(d2, h5%file_id, 'shift_thirdorder', enhet=trim(unitstr), dimensions='q-point,mode')
            call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
        if (allocated(bs%p(1)%shift4)) then
            call mem%allocate(d2, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            d2 = 0.0_r8
            do i = 1, bs%n_point
                d2(:, i) = bs%p(i)%shift4*unitfactor
            end do
            call h5%store_data(d2, h5%file_id, 'shift_fourthorder', enhet=trim(unitstr), dimensions='q-point,mode')
            call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
        if (allocated(bs%p(1)%pyroelectric_amplitude)) then
            call mem%allocate(d2, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            d2 = 0.0_r8
            do i = 1, bs%n_point
                d2(:, i) = bs%p(i)%pyroelectric_amplitude
            end do
            call h5%store_data(d2, h5%file_id, 'pyroelectric_amplitude', enhet='A*sqrt(amu)', dimensions='q-point,mode')
            call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        call mem%allocate(d2, [3, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        d2 = 0.0_r8
        do i = 1, bs%n_point
            d2(:, i) = bs%q(i)%r/lo_Bohr_to_A
        end do
        call h5%store_data(d2, h5%file_id, 'q_vector', enhet='A^-1', dimensions='q-point,xyz')
        call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! If applicable, store which q-points might be active (or not active for that matter).
        if (bs%n_active_qpoint .gt. 0) then
            storeactive: block
                integer :: iq

                ! Make a note how many points
                call h5%store_attribute(bs%n_active_qpoint, h5%file_id, 'n_active_qpoint')

                do iq = 1, bs%n_active_qpoint
                    call h5%open_group('write', 'active_qpoint_'//tochar(iq))

                    ! Store some metadata
                    call h5%store_attribute(bs%active_qpoint(iq)%nx_raman, h5%group_id, 'nx_Raman')
                    call h5%store_attribute(bs%active_qpoint(iq)%nx_IR, h5%group_id, 'nx_IR')
                    call h5%store_attribute(bs%active_qpoint(iq)%ind_path, h5%group_id, 'index_on_path')

                    if (bs%active_qpoint(iq)%nx_raman .gt. 0) then
                        call h5%store_data(bs%active_qpoint(iq)%coeff_Raman, h5%group_id, &
                                           'coeff_Raman', enhet='dimensionless')
                        call h5%store_data(bs%active_qpoint(iq)%coeff_Raman_mode, h5%group_id, &
                                           'coeff_Raman_per_mode', enhet='dimensionless')
                    end if

                    if (bs%active_qpoint(iq)%nx_IR .gt. 0) then
                        call h5%store_data(bs%active_qpoint(iq)%coeff_IR, h5%group_id, &
                                           'coeff_IR', enhet='dimensionless')
                        call h5%store_data(bs%active_qpoint(iq)%coeff_IR_mode, h5%group_id, &
                                           'coeff_IR_per_mode', enhet='dimensionless')
                    end if

                    call h5%close_group()
                end do
            end block storeactive
        else
            ! Make a note that there are no active points
            call h5%store_attribute(0, h5%file_id, 'n_active_qpoint')
        end if

        ! close the file and hdf5, if relevant.
        if (present(hdftag) .eqv. .false.) then
            call h5%close_file()
            call h5%destroy(__FILE__, __LINE__)
        end if
    end subroutine

!> Write something along a path in the BZ.
    subroutine write_dispersive_property(bs, enhet, property, filename, pdfplot)
        !> phonon dispersions with all the data
        class(lo_phonon_bandstructure), intent(in) :: bs
        !> frequency unit
        character(len=3), intent(in) :: enhet
        !> filename
        character(len=*), intent(in) :: filename
        !> produce a .gnuplot file that generates a pdf?
        logical, intent(in) :: pdfplot
        !> what property should be plotted
        character(len=*), intent(in) :: property

        character(len=1000) :: gpfn, opf, unitstr, ylabel
        character(len=10), dimension(:), allocatable :: dumstr
        real(r8), dimension(:), allocatable :: dum
        real(r8) :: unitfactor
        integer :: i, j, k, u

        ! Get the unit conversion factor
        select case (enhet)
        case ('thz')
            unitfactor = lo_frequency_hartree_to_THz
            unitstr = 'THz'
        case ('mev')
            unitfactor = lo_frequency_hartree_to_meV
            unitstr = 'meV'
        case ('icm')
            unitfactor = lo_frequency_hartree_to_icm
            unitstr = 'cm^-^1'
        case default
            unitfactor = 1.0_r8
        end select

        ! the output format for the raw data
        opf = "("//tochar(bs%n_mode + 1)//"(1X,E18.12))"
        select case (trim(property))
        case ('frequency')
            ! write the phonon dispersions
            ylabel = 'Frequency ('//trim(unitstr)//')'
        case ('groupvelocity')
            ! phonon group velocities
            unitstr = 'km/s'
            unitfactor = lo_groupvel_hartreebohr_to_ms*1E-3_r8
            ylabel = 'Group velocity ('//trim(unitstr)//')'
        case ('gruneisen')
            ! mode gruneisen parameters
            unitstr = ''
            unitfactor = 1.0_r8
            ylabel = 'Mode Gruneisen parameter'
        case ('linewidth')
            ! write phonon linewidths
            ylabel = 'Linewidth ('//trim(unitstr)//')'
        case ('shift')
            ! the anharmonic frequency shifts
            ylabel = 'Shift ('//trim(unitstr)//')'
        case ('deltaF')
            ! anharmonic contribution to the free energy
            ! I default this to meV/atom
            unitstr = 'meV'
            unitfactor = 1000.0_r8
            ylabel = 'Free energy shift ('//trim(unitstr)//'/atom)'
        case ('electronphononphasespace')
            ! electron phonon phase space
            unitstr = ''
            unitfactor = 1.0_r8
            ylabel = 'Electron-phonon phase space'
        case ('elasticfrequency')
            ! write the phonon dispersions
            ylabel = 'Frequency ('//trim(unitstr)//')'
            opf = "(7(1X,E18.12))"
        case ('threephononphasespace')
            ! write the phonon dispersions
            ylabel = 'Three-phonon phase space'
        end select

        ! Get the filename
        gpfn = trim(filename)//'.gnuplot'

        ! Dump the raw data
        allocate (dum(bs%n_mode))
        u = open_file('out', trim(filename))
        do i = 1, bs%n_point
            select case (trim(property))
            case ('frequency')
                ! write the phonon dispersions
                write (u, opf) bs%q_axis(i), bs%p(i)%omega*unitfactor
            case ('gruneisen')
                ! write gruneisen parameters
                write (u, opf) bs%q_axis(i), bs%p(i)%gruneisen*unitfactor
            case ('groupvelocity')
                ! here I have to calculate the norms
                do j = 1, bs%n_mode
                    dum(j) = norm2(bs%p(i)%vel(:, j))
                end do
                write (u, opf) bs%q_axis(i), dum*unitfactor
            case ('linewidth')
                ! write phonon linewidths
                write (u, opf) bs%q_axis(i), bs%p(i)%linewidth*unitfactor
            case ('shift')
                ! the anharmonic frequency shifts
                dum = 0.0_r8
                if (allocated(bs%p(i)%shift3)) dum = dum + bs%p(i)%shift3
                if (allocated(bs%p(i)%shift4)) dum = dum + bs%p(i)%shift4
                write (u, opf) bs%q_axis(i), dum*unitfactor
            case ('deltaF')
                ! anharmonic contribution to the free energy
                write (u, opf) bs%q_axis(i), (bs%p(i)%deltaF3 + bs%p(i)%deltaF4)*unitfactor
            case ('electronphononphasespace')
                write (u, opf) bs%q_axis(i), bs%p(i)%electronphononphasespace*unitfactor
            case ('elasticfrequency')
                ! write the phonon dispersions
                write (u, opf) bs%q_axis(i), bs%p(i)%omega(1:3)*unitfactor, bs%p(i)%omega_elastic*unitfactor
            case ('threephononphasespace')
                write (u, opf) bs%q_axis(i), bs%p(i)%threephononphasespace
            end select
        end do
        close (u)
        deallocate (dum)

        ! Dump x-axis ticks in matlab format, for neater plotting
        j = size(bs%q_axis_tick_labels, 1)
        allocate (dumstr(j))
        do i = 1, size(dumstr, 1)
            k = 0
            dumstr(i) = '          '
            do j = 1, len(dumstr(i))
                if (bs%q_axis_tick_labels(i) (j:j) .ne. '|') then
                    k = k + 1
                    dumstr(i) (k:k) = trim(bs%q_axis_tick_labels(i) (j:j))
                end if
            end do
            if (trim(dumstr(i)) .eq. 'Γ') then
                dumstr(i) = 'G'
            end if
        end do
        ! j=open_file('out','disprel_xtck')
        !     opf="("//tochar(bs%n_path+1)//"(1X,A) )"
        !     write(j,opf) dumstr
        !     opf="("//tochar(bs%n_path+1)//"(1X,F12.8) )"
        !     write(j,opf) bs%q_axis_ticks
        ! close(j)
        deallocate (dumstr)
        ! Nice gnuplot file
        if (pdfplot) then
            u = open_file('out', trim(trim(gpfn)//'_pdf'))
        else
            u = open_file('out', gpfn)
        end if
        ! Choose terminal
        if (pdfplot) then
            write (u, *) 'set terminal pdf size 8cm,7cm enhanced font "CMU Serif,10"'
            write (u, *) 'set output "'//trim(filename)//'.pdf"'
        else
            write (u, *) 'set terminal '//lo_gnuplotterminal//' size 500,350 enhanced font "CMU Serif,10"'
        end if

        ! set the ticks
        write (u, *) 'unset xtics'
        write (u, *) 'set xtics ( "'//trim(bs%q_axis_tick_labels(1))//'" 0.0 ) '
        do i = 2, bs%n_path + 1
            write (u, '(A)', advance='no') 'set xtics add ('
            write (u, '(A)', advance='no') '"'//trim(bs%q_axis_tick_labels(i))//'" '
            write (u, '(F9.6)', advance='no') real(bs%q_axis_ticks(i))
            write (u, *) ' )'
        end do

        ! set gridlines at tics, and the y-label
        write (u, *) 'set grid xtics lc rgb "#888888" lw 1 lt 0'
        write (u, *) 'set xzeroaxis linewidth 0.1 linecolor 0 linetype 1'
        write (u, *) 'set ytics scale 0.5'
        write (u, *) 'set xtics scale 0.5'
        write (u, *) 'set mytics 10'
        write (u, *) 'unset key'
        write (u, *) 'set ylabel "', trim(ylabel), '"'

        ! Plot
        if (trim(property) .eq. 'elasticfrequency') then
            write (u, '(A)', advance='no') 'plot'
            do j = 2, 4
                write (u, '(A)', advance='no') ' "'//trim(filename)//'" u 1:'//tochar(j)//' w line lc rgb "#618712"'
                write (u, '(A)') ',\'
            end do
            do j = 5, 7
                write (u, '(A)', advance='no') ' "'//trim(filename)//'" u 1:'//tochar(j)//' w line lc rgb "#612752"'
                if (j .lt. bs%n_mode + 4) then
                    write (u, '(A)') ',\'
                end if
            end do
        else
            write (u, '(A)', advance='no') 'plot'
            do j = 2, bs%n_mode + 1
                write (u, '(A)', advance='no') ' "'//trim(filename)//'" u 1:'//tochar(j)//' w line lc rgb "#618712"'
                if (j .lt. bs%n_mode + 1) then
                    write (u, '(A)') ',\'
                end if
            end do
        end if
        close (u)
    end subroutine

!> write lineshapes along a path as an hdf5 file
    subroutine write_spectral_function_to_hdf5(bs, enhet, filename)
        !> bandstructure and other things
        class(lo_phonon_bandstructure), intent(inout) :: bs
        !> the frequency unit
        character(len=3), intent(in) :: enhet
        !> the filename
        character(len=*), intent(in) :: filename

        type(lo_hdf5_helper) :: h5
        real(r8) :: unitfactor
        integer :: i, j, k
        character(len=30), dimension(:), allocatable :: dumlbl
        character(len=2000) :: lblstr, unitstr

        ! Fix the unit?
        select case (enhet)
        case ('thz')
            unitfactor = lo_frequency_Hartree_to_THz
            unitstr = 'Thz'
        case ('icm')
            unitfactor = lo_frequency_Hartree_to_icm
            unitstr = '1/cm'
        case ('mev')
            unitfactor = lo_frequency_Hartree_to_meV
            unitstr = 'meV'
        case default
            call lo_stop_gracefully(["Unknown unit, try 'thz', 'mev' or 'icm'"], lo_exitcode_param, __FILE__, __LINE__)
        end select

        ! Create a new file.
        call h5%init(__FILE__, __LINE__)
        call h5%open_file('write', trim(filename))

        ! the x-axis ticks
        call h5%store_data(bs%q_axis_ticks/lo_bohr_to_A, h5%file_id, 'q_ticks', enhet='1/A')
        ! and the tick labels, horribly annoying.
        allocate (dumlbl(size(bs%q_axis_tick_labels, 1)))
        do i = 1, size(dumlbl)
            k = 0
            dumlbl(i) = '         '
            do j = 1, len(dumlbl(i))
                if (bs%q_axis_tick_labels(i) (j:j) .ne. '|') then
                    k = k + 1
                    dumlbl(i) (k:k) = trim(bs%q_axis_tick_labels(i) (j:j))
                end if
            end do
            if (trim(dumlbl(i)) .eq. 'Γ') then
                dumlbl(i) = 'G'
            end if
        end do
        lblstr = ''
        do i = 1, size(dumlbl, 1)
            lblstr = trim(lblstr)//' '//trim(dumlbl(i))
        end do
        ! simple stuff
        call h5%store_attribute(trim(adjustl(lblstr)), h5%file_id, 'q_tick_labels')
        call h5%store_data(bs%q_axis/lo_bohr_to_A, h5%file_id, 'q_values', enhet='dimensionless')
        call h5%store_data(bs%energy_axis*unitfactor, h5%file_id, 'energy_values', enhet=trim(unitstr))
        call h5%store_attribute(trim(unitstr), h5%file_id, 'energy_unit')
        ! Store the actual spectralfunction and selfenergy. The unit conversion is a little
        ! annoying, these arrays can be huge, so I don't want to allocate them again. I just
        ! scale them, and then unscale after writing.
        bs%spectral_function = bs%spectral_function/unitfactor
        call h5%store_data(bs%spectral_function, h5%file_id, 'spectral_function', enhet='1/'//trim(unitstr), dimensions='energy,q')
        bs%spectral_function = bs%spectral_function*unitfactor

        if (allocated(bs%spectral_function_with_prefactor)) then
            bs%spectral_function_with_prefactor = bs%spectral_function_with_prefactor/unitfactor
            call h5%store_data(bs%spectral_function_with_prefactor, h5%file_id, &
                               'spectral_function_with_prefactor', enhet='dimensionless', dimensions='energy,q')
            bs%spectral_function_with_prefactor = bs%spectral_function_with_prefactor*unitfactor
        end if

        if (allocated(bs%selfenergy_real)) then
            bs%selfenergy_real = bs%selfenergy_real*unitfactor
            bs%selfenergy_imag = bs%selfenergy_imag*unitfactor
            call h5%store_data(bs%selfenergy_real, h5%file_id, &
                               'real_self_energy', enhet=trim(unitstr), dimensions='mode,energy,q')
            call h5%store_data(bs%selfenergy_imag, h5%file_id, &
                               'imaginary_self_energy', enhet=trim(unitstr), dimensions='mode,energy,q')
            bs%selfenergy_real = bs%selfenergy_real/unitfactor
            bs%selfenergy_imag = bs%selfenergy_imag/unitfactor
        end if

        ! Close the file
        call h5%close_file()
        call h5%destroy(__FILE__, __LINE__)
    end subroutine

!> measure size in memory, in bytes
    function bs_size_in_mem(b) result(mem)
        !> the dispersion relations along the path
        class(lo_phonon_bandstructure), intent(in) :: b
        !> memory in bytes
        integer :: mem

        integer :: i
        mem = 0
        mem = mem + storage_size(b)

        if (allocated(b%energy_axis)) mem = mem + storage_size(b%energy_axis)*size(b%energy_axis)
        if (allocated(b%spectral_function)) mem = mem + storage_size(b%spectral_function)*size(b%spectral_function)
        if (allocated(b%spectral_function_with_prefactor)) then
            mem = mem + storage_size(b%spectral_function_with_prefactor)*size(b%spectral_function_with_prefactor)
        end if
        if (allocated(b%selfenergy_real)) mem = mem + storage_size(b%selfenergy_real)*size(b%selfenergy_real)
        if (allocated(b%selfenergy_imag)) mem = mem + storage_size(b%selfenergy_imag)*size(b%selfenergy_imag)
        if (allocated(b%q_axis)) mem = mem + storage_size(b%q_axis)*size(b%q_axis)
        if (allocated(b%q_axis_ticks)) mem = mem + storage_size(b%q_axis_ticks)*size(b%q_axis_ticks)
        if (allocated(b%q_axis_tick_labels)) mem = mem + storage_size(b%q_axis_tick_labels)*size(b%q_axis_tick_labels)
        if (allocated(b%symb_q_start)) mem = mem + storage_size(b%symb_q_start)*size(b%symb_q_start)
        if (allocated(b%symb_q_end)) mem = mem + storage_size(b%symb_q_end)*size(b%symb_q_end)
        if (allocated(b%segment)) mem = mem + storage_size(b%segment)*size(b%segment)

        mem = mem/8
        if (allocated(b%q)) then
            do i = 1, size(b%q)
                mem = mem + b%q(i)%size_in_mem()
            end do
        end if
        if (allocated(b%p)) then
            do i = 1, size(b%p)
                mem = mem + b%p(i)%size_in_mem()
            end do
        end if
    end function

    !> destroy the structure
    !&<
    subroutine destroy(bs)
        !> the dispersion relations along the path
        class(lo_phonon_bandstructure), intent(inout) :: bs

        if ( allocated(bs%energy_axis)                      ) deallocate(bs%energy_axis)
        if ( allocated(bs%spectral_function)                ) deallocate(bs%spectral_function)
        if ( allocated(bs%spectral_function_with_prefactor) ) deallocate(bs%spectral_function_with_prefactor)
        if ( allocated(bs%selfenergy_real)                  ) deallocate(bs%selfenergy_real)
        if ( allocated(bs%selfenergy_imag)                  ) deallocate(bs%selfenergy_imag)
        if ( allocated(bs%q_axis)                           ) deallocate(bs%q_axis)
        if ( allocated(bs%q_axis_ticks)                     ) deallocate(bs%q_axis_ticks)
        if ( allocated(bs%q_axis_tick_labels)               ) deallocate(bs%q_axis_tick_labels)
        if ( allocated(bs%symb_q_start)                     ) deallocate(bs%symb_q_start)
        if ( allocated(bs%symb_q_end)                       ) deallocate(bs%symb_q_end)
        if ( allocated(bs%segment)                          ) deallocate(bs%segment)
        if ( allocated(bs%q)                                ) deallocate(bs%q)
        if ( allocated(bs%p)                                ) deallocate(bs%p)
        ! Make sure we set everything to nothing?
        bs%n_mode=-lo_hugeint
        bs%n_point=-lo_hugeint
        bs%n_path=-lo_hugeint
        bs%n_point_per_path=-lo_hugeint
    end subroutine
    !&>

! !> generate dispersion relations from elastic constants
! subroutine generate_elastic(bs,p,fc)
!     !> the bandstructure
!     class(lo_phonon_bandstructure), intent(inout) :: bs
!     !> the crystal structure
!     type(lo_crystalstructure), intent(inout) :: p
!     !> the force constant
!     type(lo_forceconstant_secondorder), intent(in) :: fc
!
!     ! complex(r8), dimension(:,:), allocatable :: D,egv
!     ! real(r8), dimension(:,:), allocatable :: Cq
!     ! integer :: i,j
!
!     write(lo_iou,*) 'FIXME GENERATE ELASTIC CONSTANTS DISPERSIONS'
!     stop
!
!     ! ! Make some space
!     ! do i=1,bs%n_point
!     !     if ( allocated(bs%p(i)%omega_elastic) .eqv. .false. ) then
!     !         lo_allocate(bs%p(i)%omega_elastic(3))
!     !     endif
!     ! enddo
!     !
!     ! ! Get the actual frequencies
!     ! lo_allocate(Cq(3,3))
!     ! lo_allocate(D(3,3))
!     ! lo_allocate(egv(3,3))
!     ! do i=1,bs%n_point
!     !     call cq_matrix(fc,p,bs%q(i)%w,Cq)
!     !     D=Cq
!     !     call lo_complex_hermitian_eigenproblem(D,bs%p(i)%omega_elastic,egv)
!     !     do j=1,3
!     !         bs%p(i)%omega_elastic(j)=sqrt(bs%p(i)%omega_elastic(j))
!     !     enddo
!     ! enddo
!     ! lo_deallocate(Cq)
!     ! lo_deallocate(D)
!     ! lo_deallocate(egv)
!     ! contains
!     !
!     ! !> Solve the complex hermition eigenvalue problem, returns them sorted
!     ! subroutine lo_complex_hermitian_eigenproblem(m,eigenvalues,eigenvectors)
!     !     !> The hermitian matrix
!     !     complex(r8), dimension(:,:), intent(in) :: m
!     !     !> The eigenvalues
!     !     real(r8), dimension(:), intent(out) :: eigenvalues
!     !     !> The eigenvectors
!     !     complex(r8), dimension(:,:), intent(out) :: eigenvectors
!     !     !
!     !     integer :: i,j,n
!     !     complex(r8), dimension(:,:), allocatable :: dum,dumegv
!     !     complex(r8), dimension(:), allocatable :: dumom
!     !     real(r8), dimension(:), allocatable :: dumreal
!     !     integer, dimension(:), allocatable :: ind
!     !     !
!     !     n=size(m,1)
!     !     allocate(dum(n,n),dumegv(n,n),ind(n),dumom(n),dumreal(n))
!     !
!     !     dum=m
!     !     call lo_zgeev(dum,dumom,dumegv)
!     !     dumreal=abs(dumom)
!     !     call qsort(dumreal,ind)
!     !
!     !     do i=1,n
!     !         j=ind(i)
!     !         eigenvalues(i)=abs(dumom(j))
!     !         eigenvectors(:,i)=dumegv(:,j)
!     !     enddo
!     !     deallocate(dum,dumegv,ind,dumom,dumreal)
!     ! end subroutine
!     !
!     ! subroutine cq_matrix(fc,uc,q,Cq)
!     !     !> the force constant
!     !     type(lo_forceconstant_secondorder), intent(in) :: fc
!     !     !> the crystal structure
!     !     type(lo_crystalstructure), intent(in) :: uc
!     !     !> the q-vector
!     !     real(r8), dimension(3), intent(in) :: q
!     !     !> the eigenvectors
!     !     real(r8), dimension(3,3), intent(out) :: Cq
!     !     !
!     !     integer :: i,j,k,l
!     !     real(r8), dimension(3) :: qv
!     !     real(r8) :: density
!     !
!     !     ! Get the dynamical matrix thing
!     !     qv=q*lo_twopi
!     !     Cq=0.0_r8
!     !     do i=1,3
!     !     do k=1,3
!     !         do j=1,3
!     !         do l=1,3
!     !             Cq(i,k)=Cq(i,k)+fc%elastic_constants_tensor(i,j,l,k)*qv(j)*qv(l)
!     !         enddo
!     !         enddo
!     !     enddo
!     !     enddo
!     !     density=sum(uc%mass)/uc%volume
!     !     Cq=Cq/density
!     ! end subroutine
! end subroutine

end module
