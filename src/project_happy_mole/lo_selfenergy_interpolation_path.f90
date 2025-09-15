submodule(lo_selfenergy_interpolation) lo_selfenergy_interpolation_path
use gottochblandat, only: lo_lorentz
implicit none

contains

!> Calculate the spectral function along a path in the BZ
module subroutine spectral_function_path_interp(ise, bs, uc, mw, mem)
    !> interpolated self-energy thing
    class(lo_interpolated_selfenergy_grid), intent(inout) :: ise
    !> bandstructure, with harmonic part pre-calculated
    type(lo_phonon_bandstructure), intent(inout) :: bs
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
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
        call mem%allocate(buf_shift_3rd, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_linewidth, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
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

        call mem%allocate(buf_re, [ise%n_energy, bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_im, [ise%n_energy, bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf_re=0.0_r8
        buf_im=0.0_r8

        if (mw%talk) call lo_progressbar_init()
        do iq = 1, bs%n_point
            if (mod(iq, mw%n) .ne. mw%r) cycle
            ! Get self-energy
            call ise%evaluate(uc,bs%q(iq)%r,bs%p(iq),buf_Re(:,:,iq),buf_im(:,:,iq),mem)
            !call ise%diagonal_selfenergy(uc, bs%q(iq)%r, bs%p(iq)%omega, bs%p(iq)%egv, buf_re(:, :, iq), buf_im(:, :, iq))
            ! Get shift and width
            do imode = 1, bs%n_mode
                buf_shift_3rd(imode, iq) = lo_linear_interpolation(ise%omega, buf_re(:, imode, iq), bs%p(iq)%omega(imode))
                buf_linewidth(imode, iq) = lo_linear_interpolation(ise%omega, buf_im(:, imode, iq), bs%p(iq)%omega(imode))
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
            im_at_gamma = (ise%omega(2) - ise%omega(1))*0.25_r8

            ! Also put a lower limit to the imaginary self-energy. This is to make
            ! the plot look sensible, nothing else. If a peak is sharper than the distance
            ! between two pixels, everything looks wonky.
            im_lower_limit = (ise%omega(2) - ise%omega(1))*0.25_r8
            im_lower_limit = (ise%omega(2) - ise%omega(1))*0.25_r8

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
                    call evaluate_spectral_function(ise%omega, bufi, bufr, bs%p(ipt)%omega(imode), bufs)
                    bufs = bufs/lo_trapezoid_integration(ise%omega, bufs)
                else
                    do ie = 1, ise%n_energy
                        bufs(ie) = lo_lorentz(ise%omega(ie), 0.0_r8, im_at_gamma)
                    end do
                    bufs = bufs/lo_trapezoid_integration(ise%omega, bufs)
                end if
                bs%spectral_function(ipt, :) = bs%spectral_function(ipt, :) + bufs
            end do
            end do
            ! And store the energy axis
            bs%energy_axis = ise%omega

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

end submodule