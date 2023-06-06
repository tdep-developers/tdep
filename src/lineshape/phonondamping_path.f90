submodule(phonondamping) phonondamping_path
!! Calculate self-energies on a path
implicit none
contains

!> Calculate the spectral function along a path in the BZ
module subroutine spectral_function_along_path(bs, uc, fc, fct, fcf, ise, qp, dr, opts, tmr, mw, mem)
    !> the bandstructure
    type(lo_phonon_bandstructure), intent(inout) :: bs
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
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> timer
    type(lo_timer), intent(inout) :: tmr
    !> mpi communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), parameter :: time_report_interval = 5.0_r8 ! how often to report progress
    real(r8), dimension(:, :, :), allocatable :: buf_re, buf_im, thermal_disp
    real(r8), dimension(:, :), allocatable :: buf_shift_3rd, buf_shift_4th, buf_linewidth
    real(r8), dimension(:, :), allocatable :: buf_xval
    real(r8), dimension(:), allocatable :: buf_saxis
    real(r8) :: timer, t0, t1
    integer, dimension(:, :), allocatable :: buf_qind
    integer :: n_pts_per_path, n_pts, solrnk

    ! Start timers
    timer = walltime()
    t0 = timer
    t1 = timer

    ! Set some basic things
    init: block
        integer :: i, j, ipath

        if (mw%talk) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'Calculating spectral function along path'
        end if

        ! Which rank do we solve serial things on?
        solrnk = 0

        ! Figure out the number of 'real' points on each path, not the
        ! eventual number of points we are going to interpolate to later.
        n_pts_per_path = 0
        do i = 1, bs%n_point_per_path, opts%stride
            n_pts_per_path = n_pts_per_path + 1
        end do
        n_pts = bs%n_path*n_pts_per_path

        ! store some information pertain
        call mem%allocate(buf_xval, [n_pts_per_path, bs%n_path], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_qind, [n_pts_per_path, bs%n_path], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        ! Make a little space for buffers
        if (mw%r .eq. solrnk) then
            call mem%allocate(buf_re, [opts%nf, dr%n_mode, n_pts], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_im, [opts%nf, dr%n_mode, n_pts], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_shift_3rd, [dr%n_mode, n_pts], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_shift_4th, [dr%n_mode, n_pts], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_linewidth, [dr%n_mode, n_pts], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_saxis, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            buf_re = 0.0_r8
            buf_im = 0.0_r8
            buf_shift_3rd = 0.0_r8
            buf_shift_4th = 0.0_r8
            buf_linewidth = 0.0_r8
            buf_saxis = 0.0_r8
        end if

        ! Store some indexing things so that I can manage the stride.
        buf_xval = 0.0_r8
        buf_qind = 0
        do ipath = 1, bs%n_path
            j = 0
            do i = 1, bs%n_point_per_path, opts%stride
                j = j + 1
                buf_xval(j, ipath) = bs%q_axis((ipath - 1)*bs%n_point_per_path + i)
                buf_qind(j, ipath) = (ipath - 1)*bs%n_point_per_path + i
            end do
        end do

        ! Also, on the head rank, make space for the final selfenergy.
        if (mw%r .eq. solrnk) then
            allocate (bs%energy_axis(opts%nf))
            allocate (bs%spectral_function(bs%n_point, opts%nf))
            allocate (bs%selfenergy_real(bs%n_point, opts%nf, dr%n_mode))
            allocate (bs%selfenergy_imag(bs%n_point, opts%nf, dr%n_mode))
            bs%energy_axis = 0.0_r8
            bs%spectral_function = 0.0_r8
            bs%selfenergy_real = 0.0_r8
            bs%selfenergy_imag = 0.0_r8
            do i = 1, bs%n_point
                allocate (bs%p(i)%linewidth(bs%n_mode))
                allocate (bs%p(i)%shift3(bs%n_mode))
                allocate (bs%p(i)%shift4(bs%n_mode))
                bs%p(i)%linewidth = 0.0_r8
                bs%p(i)%shift3 = 0.0_r8
                bs%p(i)%shift4 = 0.0_r8
            end do
        end if

        ! I will also need the thermal displacement matrix, eventually.
        call mem%allocate(thermal_disp, [3, 3, uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        thermal_disp = 0.0_r8
        call dr%thermal_displacement_matrix(qp, uc, opts%temperature, thermal_disp, mw, mem)

        if (mw%talk) then
            write (lo_iou, *) '              number of point:', bs%n_point, '(', n_pts, ')'
            write (lo_iou, *) '              number of paths:', bs%n_path
            write (lo_iou, *) '       number points per path:', bs%n_point_per_path, '(', n_pts_per_path, ')'
            write (lo_iou, *) '                       stride:', opts%stride
            write (lo_iou, *) '      number of energy values:', opts%nf
        end if
    end block init

    ! Get the self-energies
    selfenergy: block
        type(lo_spectralfunction_helper) :: isf
        type(lo_phonon_selfenergy) :: se
        real(r8), dimension(3) :: qdir
        real(r8) :: f0
        integer :: ipath, ipt, jpt, iqp, imode

        do ipath = 1, bs%n_path
        do ipt = 1, n_pts_per_path
            jpt = (ipath - 1)*n_pts_per_path + ipt
            iqp = buf_qind(ipt, ipath)

            ! Generate the self-energy
            qdir = bs%segment(ipath)%r2 - bs%segment(ipath)%r1
            f0 = norm2(qdir)
            if (f0 .gt. 1E-10_r8) then
                qdir = qdir/f0
            else
                ! just pick something, I dunno
                qdir = [1.0_r8, 0.0_r8, 0.0_r8]
            end if

            call se%generate(bs%q(iqp), qdir, bs%p(iqp), uc, fc, fct, fcf, ise, isf, qp, dr, opts, tmr, mw, mem, verbosity=-1)

            if (mw%r .eq. solrnk) then
                ! Store self-energy in buffers
                buf_re(:, :, jpt) = se%re_3ph + se%re_4ph
                buf_im(:, :, jpt) = se%im_3ph + se%im_iso
                ! Store the shifts from third and fourth order while they are available
                do imode = 1, dr%n_mode
                    buf_shift_3rd(imode, jpt) = lo_linear_interpolation(se%energy_axis, se%re_3ph(:, imode), bs%p(iqp)%omega(imode))
                    buf_shift_4th(imode, jpt) = lo_linear_interpolation(se%energy_axis, se%re_4ph(:, imode), bs%p(iqp)%omega(imode))
                end do
                ! While we have the self-energy available it makes sense to
                ! stuff away the energy axis for the self-energy. First a sanity
                ! check!
                if (jpt .gt. 1) then
                    if (sum(abs(buf_saxis - se%energy_axis))/opts%nf .gt. 1E-50_r8) then
                        call lo_stop_gracefully(['Inconsistent energy axis, should be impossible'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
                    end if
                end if
                buf_saxis = se%energy_axis
            end if

            if (mw%talk) then
                if (walltime() - t0 .gt. time_report_interval) then
                    call lo_looptimer('... spectralfunction along path', timer, walltime(), jpt, n_pts)
                    t0 = walltime()
                end if
            end if
        end do
        end do
    end block selfenergy

    ! Sensible place to start is to fix the values at small q, they
    ! can look weird otherwise. This is fast, so I do it serially.
    if (mw%r .eq. solrnk) then
        fixsmallq: block
            real(r8), dimension(:), allocatable :: bufr, bufi
            real(r8), dimension(3) :: qv1, qv2
            real(r8) :: minomega, f0
            integer :: i_gamma, i_okpoint, ipt, jpt, iqp, ipath, imode
            logical :: beginning

            call mem%allocate(bufr, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufi, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            bufr = 0.0_r8
            bufi = 0.0_r8

            ! Smallest frequency that could have a proper selfenergy
            minomega = dr%omega_min*0.5_r8

            t0 = walltime()
            call lo_progressbar_init()
            pathloop: do ipath = 1, bs%n_path
                ! Locate Gamma on this path
                qv1 = bs%segment(ipath)%r1 - uc%bz%gshift(bs%segment(ipath)%r1 + lo_degenvector)
                qv2 = bs%segment(ipath)%r2 - uc%bz%gshift(bs%segment(ipath)%r2 + lo_degenvector)
                ! does it contain gamma?
                if (norm2(qv1) .gt. lo_tol .and. norm2(qv2) .gt. lo_tol) cycle pathloop
                ! yup, contains Gamma, locate it.
                if (norm2(qv1) .lt. lo_tol) then
                    i_gamma = 1
                    beginning = .true.
                else
                    i_gamma = n_pts_per_path
                    beginning = .false.
                end if

                ! Now try to fix each band
                modeloop: do imode = 1, bs%n_mode
                    ! Skip modes that don't need fixing, i.e. anything not acoustic.
                    if (bs%p(i_gamma)%omega(imode) .gt. minomega) cycle modeloop
                    ! Locate a neat point not gamma?
                    i_okpoint = -1
                    if (beginning) then
                        do ipt = 1, n_pts_per_path
                            iqp = buf_qind(ipt, ipath)
                            if (bs%p(iqp)%omega(imode) .gt. minomega) then
                                i_okpoint = ipt
                                exit
                            end if
                        end do
                    else
                        do ipt = n_pts_per_path, 1, -1
                            iqp = buf_qind(ipt, ipath)
                            if (bs%p(iqp)%omega(imode) .gt. minomega) then
                                i_okpoint = ipt
                                exit
                            end if
                        end do
                    end if
                    if (i_okpoint .eq. -1) then
                        call lo_stop_gracefully(['I do not know how to index points'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                    end if
                    ! Grab selfenergies for interpolation
                    jpt = (ipath - 1)*n_pts_per_path + i_okpoint
                    bufr = buf_re(:, imode, jpt)
                    bufi = buf_im(:, imode, jpt)
                    if (beginning) then
                        do ipt = i_gamma, i_okpoint
                            ! 0-1 scalar parameter for interpolation? makes sense perhaps.
                            f0 = lo_chop((ipt - i_gamma)/real(i_okpoint - i_gamma, r8), lo_sqtol)
                            jpt = (ipath - 1)*n_pts_per_path + ipt
                            buf_re(:, imode, jpt) = f0*bufr
                            buf_im(:, imode, jpt) = f0*bufi
                        end do
                    else
                        do ipt = i_okpoint, i_gamma
                            ! 0-1 scalar parameter for interpolation
                            f0 = lo_chop((ipt - i_gamma)/real(i_okpoint - i_gamma, r8), lo_sqtol)
                            jpt = (ipath - 1)*n_pts_per_path + ipt
                            buf_re(:, imode, jpt) = f0*bufr
                            buf_im(:, imode, jpt) = f0*bufi
                        end do
                    end if
                end do modeloop
                if (ipath .lt. bs%n_path) call lo_progressbar(' ... fixing tiny q', ipath, bs%n_path, walltime() - t0)
            end do pathloop
            call mem%deallocate(bufr, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(bufi, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            call lo_progressbar(' ... fixing tiny q', bs%n_path, bs%n_path, walltime() - t0)
        end block fixsmallq

        ! Now, interpolate the self-energy to all the points. Still only on the head rank.
        interpolate: block
            real(r8), dimension(:, :, :), allocatable :: bufr, bufi
            real(r8), dimension(:), allocatable :: bufvi, bufvj
            real(r8) :: x, f0, f1, f2
            integer :: imode, jmode, ipath, ie, ipt, jpt, i

            t0 = walltime()
            call mem%allocate(bufr, [n_pts_per_path, opts%nf, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufi, [n_pts_per_path, opts%nf, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufvi, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufvj, n_pts_per_path, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            bufr = 0.0_r8
            bufi = 0.0_r8
            bufvi = 0.0_r8
            bufvj = 0.0_r8

            call lo_progressbar_init()
            do ipath = 1, bs%n_path
                ! prefetch a buffer?
                bufr = 0.0_r8
                bufi = 0.0_r8
                do ipt = 1, n_pts_per_path
                    jpt = (ipath - 1)*n_pts_per_path + ipt
                    do imode = 1, dr%n_mode
                        bufr(ipt, :, imode) = buf_re(:, imode, jpt)
                        bufi(ipt, :, imode) = buf_im(:, imode, jpt)
                    end do
                end do
                ! Interpolate to the right place
                do ipt = 1, bs%n_point_per_path
                    jpt = (ipath - 1)*bs%n_point_per_path + ipt
                    x = bs%q_axis(jpt)
                    ! Interpolate the raw self-energies
                    do imode = 1, dr%n_mode
                    do ie = 1, opts%nf
                        bs%selfenergy_real(jpt, ie, imode) = lo_linear_interpolation(buf_xval(:, ipath), bufr(:, ie, imode), x)
                        bs%selfenergy_imag(jpt, ie, imode) = lo_linear_interpolation(buf_xval(:, ipath), bufi(:, ie, imode), x)
                    end do
                    end do
                    ! Interpolate the shifts and linewidhs along the paths and to the harmonic frequencies
                    do imode = 1, dr%n_mode
                        bufvj = buf_shift_3rd(imode, (ipath - 1)*n_pts_per_path + 1:ipath*n_pts_per_path)
                        bs%p(jpt)%shift3(imode) = lo_linear_interpolation(buf_xval(:, ipath), bufvj, x)

                        bufvj = buf_shift_4th(imode, (ipath - 1)*n_pts_per_path + 1:ipath*n_pts_per_path)
                        bs%p(jpt)%shift4(imode) = lo_linear_interpolation(buf_xval(:, ipath), bufvj, x)

                        bufvi = bs%selfenergy_imag(jpt, :, imode)
                        bs%p(jpt)%linewidth(imode) = lo_linear_interpolation(buf_saxis, bufvi, bs%p(jpt)%omega(imode))
                    end do
                    call lo_progressbar(' ... interpolating selfenergy', jpt, bs%n_point, walltime() - t0)
                end do
            end do
            call mem%deallocate(bufr, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(bufi, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(bufvi, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(bufvj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            ! While here, make sure the shifts and widths are as degenerate as they should
            do ipt = 1, bs%n_point
            do imode = 1, dr%n_mode
                f0 = 0.0_r8
                f1 = 0.0_r8
                do i = 1, bs%p(ipt)%degeneracy(imode)
                    jmode = bs%p(ipt)%degenmode(i, imode)
                    f0 = f0 + bs%p(ipt)%shift3(jmode)
                    f1 = f1 + bs%p(ipt)%shift4(jmode)
                    f2 = f2 + bs%p(ipt)%linewidth(jmode)
                end do
                f0 = f0/real(bs%p(ipt)%degeneracy(imode), r8)
                f1 = f1/real(bs%p(ipt)%degeneracy(imode), r8)
                f2 = f2/real(bs%p(ipt)%degeneracy(imode), r8)
                do i = 1, bs%p(ipt)%degeneracy(imode)
                    jmode = bs%p(ipt)%degenmode(i, imode)
                    bs%p(ipt)%shift3(jmode) = f0
                    bs%p(ipt)%shift4(jmode) = f1
                    bs%p(ipt)%linewidth(jmode) = f2
                end do
            end do
            end do
        end block interpolate

        ! Now that I have the self-energies all over the place, makes sense to create the intensities
        spf: block
            !real(r8), dimension(:,:,:), allocatable :: sigma
            real(r8), dimension(:), allocatable :: bufr, bufi, bufs
            !real(r8) :: f0,f1
            real(r8) :: im_at_gamma, im_lower_limit
            integer :: ipt, imode, ie

            ! put something at gamma to make the intensities not weird.
            im_at_gamma = (buf_saxis(2) - buf_saxis(1))*0.25_r8

            ! Also put a lower limit to the imaginary self-energy. This is to make
            ! the plot look sensible, nothing else. If a peak is sharper than the distance
            ! between two pixels, everything looks wonky.
            if (opts%minsmear .gt. 0) then
                im_lower_limit = (buf_saxis(2) - buf_saxis(1))*0.25_r8*opts%minsmear
            else
                im_lower_limit = (buf_saxis(2) - buf_saxis(1))*0.25_r8
            end if

            call mem%allocate(bufr, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufi, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufs, opts%nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            bufr = 0.0_r8
            bufi = 0.0_r8
            bufs = 0.0_r8

            bs%spectral_function = 0.0_r8
            ! !bs%spectral_function_with_prefactor=0.0_r8
            do ipt = 1, bs%n_point
            do imode = 1, dr%n_mode
                if (bs%p(ipt)%omega(imode) .gt. lo_freqtol) then
                    bufi = bs%selfenergy_imag(ipt, :, imode)
                    bufr = bs%selfenergy_real(ipt, :, imode)
                    ! make sure we at least have some smearing?
                    bufi = max(bufi, im_lower_limit)
                    call evaluate_spectral_function(buf_saxis, bufi, bufr, bs%p(ipt)%omega(imode), bufs)
                    bufs = bufs/lo_trapezoid_integration(buf_saxis, bufs)
                else
                    do ie = 1, opts%nf
                        bufs(ie) = lo_lorentz(buf_saxis(ie), 0.0_r8, im_at_gamma)
                    end do
                    bufs = bufs/lo_trapezoid_integration(buf_saxis, bufs)
                end if
                bs%spectral_function(ipt, :) = bs%spectral_function(ipt, :) + bufs
                ! And add the thermal prefactor thing
                !f1=thermal_pref(bs%q(i)%r,bs%p(i)%egv(:,j),uc,opts%temperature,bs%p(i)%omega(j),sigma)
                !bs%spectral_function_with_prefactor(i,:)=bs%spectral_function_with_prefactor(i,:)+dum*f1
            end do
            end do
            ! And store the energy axis
            bs%energy_axis = buf_saxis
            !
            call mem%deallocate(bufr, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(bufi, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(bufs, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block spf
    end if ! if mw%r == solrnk

    ! And some cleanup
    call mem%deallocate(buf_xval, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_qind, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if (mw%r .eq. solrnk) then
        call mem%deallocate(buf_re, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_im, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_shift_3rd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_shift_4th, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_linewidth, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(buf_saxis, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if
end subroutine

!> Calculate the spectral function along a path in the BZ
module subroutine spectral_function_along_path_interp(bs, uc, ise, opts, mw, mem)
    !> the bandstructure
    type(lo_phonon_bandstructure), intent(inout) :: bs
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> interpolated self-energy thing
    type(lo_interpolated_selfenergy), intent(in) :: ise
    !> all settings
    type(lo_opts), intent(in) :: opts
    !> mpi communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !real(r8), parameter :: time_report_interval=5.0_r8 ! how often to report progress
    real(r8), dimension(:, :, :), allocatable :: buf_re, buf_im
    real(r8), dimension(:, :), allocatable :: buf_shift_3rd, buf_linewidth
    !real(r8), dimension(:,:), allocatable :: buf_xval
    !real(r8), dimension(:), allocatable :: buf_sf
    real(r8) :: timer, t0, t1
    !integer, dimension(:,:), allocatable :: buf_qind
    !integer :: n_pts_per_path,n_pts,solrnk
    !integer :: nb
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
            write (lo_iou, *) 'Calculating spectral function along path with fancy interpolation'
        end if

        ! Which rank do we solve serial things on?
        solrnk = 0

        ! Some space for the spectral function. Will adjust if memory becomes problem.
        call mem%allocate(buf_re, [ise%n_energy, bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_im, [ise%n_energy, bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_shift_3rd, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_linewidth, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf_re = 0.0_r8
        buf_im = 0.0_r8
        buf_shift_3rd = 0.0_r8
        buf_linewidth = 0.0_r8

        ! Also, on the head rank, make space for the final selfenergy.
        if (mw%r .eq. solrnk) then
            allocate (bs%energy_axis(ise%n_energy))
            allocate (bs%spectral_function(bs%n_point, ise%n_energy))
            allocate (bs%selfenergy_real(bs%n_point, ise%n_energy, bs%n_mode))
            allocate (bs%selfenergy_imag(bs%n_point, ise%n_energy, bs%n_mode))
            bs%energy_axis = 0.0_r8
            bs%spectral_function = 0.0_r8
            bs%selfenergy_real = 0.0_r8
            bs%selfenergy_imag = 0.0_r8
            do i = 1, bs%n_point
                allocate (bs%p(i)%linewidth(bs%n_mode))
                allocate (bs%p(i)%shift3(bs%n_mode))
                allocate (bs%p(i)%shift4(bs%n_mode))
                bs%p(i)%linewidth = 0.0_r8
                bs%p(i)%shift3 = 0.0_r8
                bs%p(i)%shift4 = 0.0_r8
            end do
        end if

        if (mw%talk) then
            t1 = walltime()
            write (lo_iou, *) '... made space (', tochar(t1 - t0), 's)'
            t0 = t1
        end if

    end block init

    ! Get the self-energies
    selfenergy: block
        integer :: iq, imode

        if (mw%talk) call lo_progressbar_init()
        do iq = 1, bs%n_point
            if (mod(iq, mw%n) .ne. mw%r) cycle
            ! Get self-energy
            call ise%diagonal_selfenergy(uc, bs%q(iq)%r, bs%p(iq)%omega, bs%p(iq)%egv, buf_re(:, :, iq), buf_im(:, :, iq))
            ! Get shift and width
            do imode = 1, bs%n_mode
                buf_shift_3rd(imode, iq) = lo_linear_interpolation(ise%energy, buf_re(:, imode, iq), bs%p(iq)%omega(imode))
                buf_linewidth(imode, iq) = lo_linear_interpolation(ise%energy, buf_im(:, imode, iq), bs%p(iq)%omega(imode))
            end do
            if (mw%talk .and. iq .lt. bs%n_point) then
                t1 = walltime()
                call lo_progressbar(' ... spectral functions', iq, bs%n_point, t1 - t0)
            end if
        end do
        call mw%allreduce('sum', buf_re)
        call mw%allreduce('sum', buf_im)
        call mw%allreduce('sum', buf_shift_3rd)
        call mw%allreduce('sum', buf_linewidth)

        ! Store in the right place
        if (mw%r .eq. solrnk) then
            do iq = 1, bs%n_point
            do imode = 1, bs%n_mode
                bs%p(iq)%linewidth(imode) = buf_linewidth(imode, iq)
                bs%p(iq)%shift3(imode) = buf_shift_3rd(imode, iq)
                bs%selfenergy_imag(iq, :, imode) = buf_im(:, imode, iq)
                bs%selfenergy_real(iq, :, imode) = buf_re(:, imode, iq)
            end do
            end do
        end if

        if (mw%talk) then
            t1 = walltime()
            call lo_progressbar(' ... spectral functions', bs%n_point, bs%n_point, t1 - t0)
            t0 = t1
        end if
    end block selfenergy

    ! Now that I have the self-energies all over the place, makes sense to create the intensities
    if (mw%r .eq. solrnk) then
        spf: block
            real(r8), dimension(:), allocatable :: bufr, bufi, bufs
            real(r8) :: im_at_gamma, im_lower_limit
            integer :: ipt, imode, ie

            ! put something at gamma to make the intensities not weird.
            im_at_gamma = (ise%energy(2) - ise%energy(1))*0.25_r8

            ! Also put a lower limit to the imaginary self-energy. This is to make
            ! the plot look sensible, nothing else. If a peak is sharper than the distance
            ! between two pixels, everything looks wonky.
            if (opts%minsmear .gt. 0) then
                im_lower_limit = (ise%energy(2) - ise%energy(1))*0.25_r8*opts%minsmear
            else
                im_lower_limit = (ise%energy(2) - ise%energy(1))*0.25_r8
            end if

            call mem%allocate(bufr, ise%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufi, ise%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufs, ise%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            bufr = 0.0_r8
            bufi = 0.0_r8
            bufs = 0.0_r8

            bs%spectral_function = 0.0_r8

            do ipt = 1, bs%n_point
            do imode = 1, bs%n_mode
                if (bs%p(ipt)%omega(imode) .gt. lo_freqtol) then
                    bufi = buf_im(:, imode, ipt)
                    bufr = buf_re(:, imode, ipt)
                    ! make sure we at least have some smearing?
                    bufi = max(bufi, im_lower_limit)
                    call evaluate_spectral_function(ise%energy, bufi, bufr, bs%p(ipt)%omega(imode), bufs)
                    bufs = bufs/lo_trapezoid_integration(ise%energy, bufs)
                else
                    do ie = 1, ise%n_energy
                        bufs(ie) = lo_lorentz(ise%energy(ie), 0.0_r8, im_at_gamma)
                    end do
                    bufs = bufs/lo_trapezoid_integration(ise%energy, bufs)
                end if
                bs%spectral_function(ipt, :) = bs%spectral_function(ipt, :) + bufs
            end do
            end do
            ! And store the energy axis
            bs%energy_axis = ise%energy

            call mem%deallocate(bufr, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(bufi, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(bufs, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block spf
    end if

    call mem%deallocate(buf_re, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_im, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_shift_3rd, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_linewidth, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

! !> Calculate the spectral function along a path in the BZ
! subroutine spectralfunction_along_path(bs,uc,fc,fct,fcf,qp,dr,opts,mw,mem)
!     !> the bandstructure
!     type(lo_phonon_bandstructure), intent(inout) :: bs
!     !> crystal structure
!     type(lo_crystalstructure), intent(inout) :: uc
!     !> second order force constant
!     type(lo_forceconstant_secondorder), intent(inout) :: fc
!     !> third order force constant
!     type(lo_forceconstant_thirdorder), intent(in) :: fct
!     !> fourth order force constant
!     type(lo_forceconstant_fourthorder), intent(in) :: fcf
!     !> q-point mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> harmonic properties on this mesh
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> all settings
!     type(lo_opts), intent(in) :: opts
!     !> mpi communicator
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> memory tracker
!     type(lo_mem_helper), intent(inout) :: mem
!
!     real(r8), parameter :: timereport=10.0_r8
!     real(r8), dimension(:,:,:), allocatable :: imbuf,rebuf
!     real(r8), dimension(:), allocatable :: seax,inax
!     real(r8) :: timer
!
!     ! Make some space and things like that
!     init: block
!         real(r8) :: f0
!         integer :: i
!
!         timer=walltime()
!         ! Figure out how much space is needed to store the buffers, and maybe print a warning
!         ! in case it's some really large amount
!         f0=4*dr%n_mode*opts%nf*bs%n_point*4.0_r8/1024/1024
!         if ( f0 .gt. 200.0_r8 .and. mw%talk ) then
!             write(*,*) '... Will use at least ',tochar(int(f0)),'MB memory per rank, probably a lot more.'
!         endif
!         ! Make space for linewidth, shifts and so on
!         do i=1,bs%n_point
!             allocate(bs%p(i)%linewidth(bs%n_mode))
!             allocate(bs%p(i)%shift3(bs%n_mode))
!             allocate(bs%p(i)%shift4(bs%n_mode))
!             bs%p(i)%linewidth=0.0_r8
!             bs%p(i)%shift3=0.0_r8
!             bs%p(i)%shift4=0.0_r8
!         enddo
!         ! Space for intensity
!         allocate(bs%spectral_function(bs%n_point,opts%nf))
!         allocate(bs%spectral_function_with_prefactor(bs%n_point,opts%nf))
!         allocate(bs%selfenergy_real(bs%n_point,opts%nf,dr%n_mode))
!         allocate(bs%selfenergy_imag(bs%n_point,opts%nf,dr%n_mode))
!         allocate(rebuf(bs%n_point,opts%nf,dr%n_mode))
!         allocate(imbuf(bs%n_point,opts%nf,dr%n_mode))
!         allocate(bs%energy_axis(opts%nf))
!         bs%energy_axis=0.0_r8
!         bs%selfenergy_real=0.0_r8
!         bs%selfenergy_imag=0.0_r8
!         bs%spectral_function=0.0_r8
!         bs%spectral_function_with_prefactor=0.0_r8
!         lo_allocate(seax(opts%nf))
!         lo_allocate(inax(opts%nf))
!         seax=0.0_r8
!         inax=0.0_r8
!     end block init
!
!     ! Get the self-energies
!     selfenergy: block
!         type(lo_phonon_selfenergy) :: se
!         real(r8), dimension(:), allocatable :: x,y,z
!         real(r8) :: f0,f1,t0
!         integer :: i,j,l,ctr,path,q,nq_per_path,band
!
!         t0=walltime()
!
!         ! Count q-points per path, and make space for dummy arrays
!         nq_per_path=0
!         do q=1,bs%n_point_per_path,opts%stride
!             nq_per_path=nq_per_path+1
!         enddo
!         lo_allocate(x(nq_per_path))
!         lo_allocate(y(nq_per_path))
!         lo_allocate(z(nq_per_path))
!
!         rebuf=0.0_r8
!         imbuf=0.0_r8
!         if ( mw%talk ) call lo_progressbar_init()
!         do path=1,bs%n_path
!         do q=1,bs%n_point_per_path,opts%stride
!             i=(path-1)*bs%n_point_per_path+q
!             ! get the self-energies
!             call se%generate(bs%q(i),bs%p(i),uc,fc,fct,fcf,qp,dr,opts,mw,mem)
!
!             ! store them
!             rebuf(i,:,:)=se%re_3ph+se%re_4ph
!             imbuf(i,:,:)=se%im_3ph+se%im_iso
!
!             ! Interpolate the shifts
!             do j=1,se%n_mode
!                 bs%p(i)%shift3(j)=lo_linear_interpolation( se%energy_axis_selfenergy,se%re_3ph(:,j),bs%p(i)%omega(j) )
!                 bs%p(i)%shift4(j)=lo_linear_interpolation( se%energy_axis_selfenergy,se%re_4ph(:,j),bs%p(i)%omega(j) )
!             enddo
!
!             if ( mw%talk ) then
!                 if ( walltime()-t0 .gt. timereport ) then
!                     call lo_looptimer('... spectralfunction along path',timer,walltime(),i,bs%n_point)
!                     t0=walltime()
!                 endif
!             endif
!         enddo
!         enddo
!
!         if ( opts%stride .gt. 1 ) then
!             t0=walltime()
!             bs%selfenergy_real=0.0_r8
!             bs%selfenergy_imag=0.0_r8
!             ! Interpolate the self-energy to all q
!             ctr=0
!             if ( mw%talk ) call lo_progressbar_init()
!             do path=1,bs%n_path
!                 ! fetch the x-values for this path
!                 l=0
!                 do q=1,bs%n_point_per_path,opts%stride
!                     i=(path-1)*bs%n_point_per_path+q
!                     l=l+1
!                     x(l)=bs%q_axis(i)
!                 enddo
!
!                 do band=1,dr%n_mode
!                     do j=1,opts%nf
!                         ! fetch self-energies
!                         l=0
!                         y=0.0_r8
!                         z=0.0_r8
!                         do q=1,bs%n_point_per_path,opts%stride
!                             i=(path-1)*bs%n_point_per_path+q
!                             l=l+1
!                             y(l)=rebuf(i,j,band)
!                             z(l)=imbuf(i,j,band)
!                         enddo
!                         ! interpolate self-energies
!                         do q=1,bs%n_point_per_path
!                             i=(path-1)*bs%n_point_per_path+q
!                             f0=lo_linear_interpolation(x,y,bs%q_axis(i))
!                             f1=lo_linear_interpolation(x,z,bs%q_axis(i))
!                             bs%selfenergy_real(i,j,band)=f0
!                             bs%selfenergy_imag(i,j,band)=max(f1,0.0_r8)
!                         enddo
!                     enddo
!                     ctr=ctr+1
!                     if ( mw%talk ) call lo_progressbar(' ... interpolating selfenergy',ctr,bs%n_path*dr%n_mode,walltime()-t0)
!                 enddo
!             enddo
!             rebuf=bs%selfenergy_real
!             imbuf=bs%selfenergy_imag
!         endif
!         ! and the different intensity axes
!         seax=se%energy_axis_selfenergy
!         inax=se%energy_axis_spectralfunction
!         bs%selfenergy_real=rebuf
!         bs%selfenergy_imag=imbuf
!     end block selfenergy
!
!     ! For tiny q close to Gamma, I linearly interpolate from the next best point, to get neat-looking plots.
!     newsmallq: block
!         real(r8), dimension(:), allocatable :: dumre,dumim,re0,im0
!         real(r8), dimension(3) :: qv1,qv2
!         real(r8) :: minomega,omshift,f0,t0
!         integer :: path,band,i_gamma,i_okpoint
!         integer :: i,ii,str
!         logical :: beginning
!
!         if ( opts%stride .gt. 0 ) then
!             str=opts%stride
!         else
!             str=1
!         endif
!
!         t0=walltime()
!         ! Temporarily store self-energies
!         lo_allocate(dumre(opts%nf))
!         lo_allocate(dumim(opts%nf))
!         lo_allocate(re0(opts%nf))
!         lo_allocate(im0(opts%nf))
!         minomega=dr%omega_min*0.5_r8
!
!         if ( mw%talk ) call lo_progressbar_init()
!
!         do path=1,bs%n_path
!             ! Locate Gamma on this path
!             qv1=bs%segment(path)%r1-uc%bz%gshift( bs%segment(path)%r1 + lo_degenvector )
!             qv2=bs%segment(path)%r2-uc%bz%gshift( bs%segment(path)%r2 + lo_degenvector )
!             ! does it contain gamma?
!             if ( norm2(qv1) .gt. lo_tol .and. norm2(qv2) .gt. lo_tol ) cycle
!             ! Locate gamma?
!             if ( norm2(qv1) .lt. lo_tol ) then
!                 i_gamma=(path-1)*bs%n_point_per_path+1
!                 beginning=.true.
!             else
!                 i_gamma=path*bs%n_point_per_path
!                 beginning=.false.
!             endif
!             ! Now try to fix each band
!             do band=1,bs%n_mode
!                 ! Skip bands that are ok
!                 if ( bs%p(i_gamma)%omega(band) .gt. minomega ) cycle
!                 ! Locate a neat point not gamma?
!                 i_okpoint=-1
!                 if ( beginning ) then
!                     do i=1,bs%n_point_per_path,opts%stride
!                         ii=(path-1)*bs%n_point_per_path+i
!                         if ( bs%p(ii)%omega(band) .gt. minomega ) then
!                             i_okpoint=ii
!                             exit
!                         endif
!                     enddo
!                 else
!                     do i=bs%n_point_per_path,1,-opts%stride
!                         ii=(path-1)*bs%n_point_per_path+i
!                         if ( bs%p(ii)%omega(band) .gt. minomega ) then
!                             i_okpoint=ii
!                             exit
!                         endif
!                     enddo
!                 endif
!                 ! Grab decent self-energy for interpolation
!                 dumre=bs%selfenergy_real(i_okpoint,:,band)
!                 dumim=bs%selfenergy_imag(i_okpoint,:,band)
!                 ! Now try to shift it?
!                 if ( beginning ) then
!                     do i=i_gamma,i_okpoint
!                         ! 0-1 scalar parameter for interpolation
!                         f0=lo_chop((i-i_gamma)/real(i_okpoint-i_gamma,r8),lo_sqtol)
!                         ! How much should I shift omega?
!                         omshift=bs%p(i)%omega(band)-bs%p(i_okpoint)%omega(band)
!                         ! Shift the self-energy to this q-point, rigidly I believe
!                         call lo_put_function_on_new_axis(seax,dumre,seax-omshift,re0)
!                         call lo_put_function_on_new_axis(seax,dumim,seax-omshift,im0)
!                         rebuf(i,:,band)=re0*f0
!                         imbuf(i,:,band)=max(im0,0.0_r8)*f0
!                     enddo
!                 else
!                     do i=i_okpoint,i_gamma
!                         ! 0-1 scalar parameter for interpolation
!                         f0=lo_chop((i-i_gamma)/real(i_okpoint-i_gamma,r8),lo_sqtol)
!                         ! How much should I shift omega?
!                         omshift=bs%p(i)%omega(band)-bs%p(i_okpoint)%omega(band)
!                         ! Shift the self-energy to this q-point, rigidly I believe
!                         call lo_put_function_on_new_axis(seax,dumre,seax-omshift,re0)
!                         call lo_put_function_on_new_axis(seax,dumim,seax-omshift,im0)
!                         rebuf(i,:,band)=re0*f0
!                         imbuf(i,:,band)=max(im0,0.0_r8)*f0
!                     enddo
!                 endif
!             enddo
!             ! Report
!             if ( mw%talk .and. path .lt. bs%n_path ) call lo_progressbar(' ... fixing tiny q',path,bs%n_path,walltime()-t0)
!         enddo
!         ! Store the self-energies
!         bs%selfenergy_real=rebuf
!         bs%selfenergy_imag=imbuf
!         if ( mw%talk ) call lo_progressbar(' ... fixing tiny q',bs%n_path,bs%n_path,walltime()-t0)
!     end block newsmallq
!
!     ! Get the intensities
!     intensities: block
!         real(r8), dimension(:,:,:), allocatable :: sigma
!         real(r8), dimension(:), allocatable :: dum
!         real(r8) :: f0,f1,minsmear
!         integer :: i,j,k
!
!         ! put something at gamma to make the intensities not weird.
!         f0=(seax(2)-seax(1))*0.25_r8
!
!         ! add minimum smearing, if desired
!         if ( opts%minsmear .gt. 0 ) then
!             minsmear=(seax(2)-seax(1))*opts%minsmear*0.25_r8
!             imbuf=max(imbuf,minsmear)
!         else
!             minsmear=(seax(2)-seax(1))*0.25_r8
!         endif
!
!         ! Also get the thermal displacement matrices for DW factors
!         lo_allocate(sigma(3,3,uc%na))
!         call dr%thermal_displacement_matrix(qp,uc,opts%temperature,sigma,mw,mem)
!
!         lo_allocate(dum(opts%nf))
!         do i=1,bs%n_point
!         do j=1,dr%n_mode
!             if ( bs%p(i)%omega(j) .gt. lo_freqtol ) then
!                 call getintensity(seax,imbuf(i,:,j),rebuf(i,:,j),bs%p(i)%omega(j),inax,dum)
!                 dum=dum/lo_trapezoid_integration(inax,dum)
!             else
!                 do k=1,opts%nf
!                     dum(k)=lo_lorentz(inax(k),0.0_r8,f0)
!                 enddo
!                 dum=dum/lo_trapezoid_integration(inax,dum)
!             endif
!             bs%spectral_function(i,:)=bs%spectral_function(i,:)+dum
!
!             f1=thermal_pref(bs%q(i)%r,bs%p(i)%egv(:,j),uc,opts%temperature,bs%p(i)%omega(j),sigma)
!             bs%spectral_function_with_prefactor(i,:)=bs%spectral_function_with_prefactor(i,:)+dum*f1
!             ! And the prefactor thingy
!         enddo
!         enddo
!         bs%energy_axis=inax
!
!         ! Also store the linewidth at the harmonic frequencies, as well as the shifts
!         do i=1,bs%n_point
!         do j=1,dr%n_mode
!             if ( bs%p(i)%omega(j) .gt. lo_freqtol ) then
!                 f0=lo_linear_interpolation(seax,imbuf(i,:,j),bs%p(i)%omega(j))
!                 bs%p(i)%linewidth(j)=max(f0,0.0_r8)
!             else
!                 bs%p(i)%linewidth(j)=0.0_r8
!             endif
!         enddo
!         enddo
!     end block intensities
!
! end subroutine

end submodule
