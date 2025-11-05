submodule(lo_evaluate_phonon_self_energy) lo_evaluate_phonon_self_energy_helpers
implicit none
contains

module subroutine add_quadratic_to_fix_normalization(omega,sigmare,sigmaim,x,initial_residual,mw)
    !> harmonic frequencies
    real(r8), dimension(:), intent(in) :: omega
    !> real part of self-energy
    real(r8), dimension(:,:), intent(inout) :: sigmare
    !> imaginary part of self-energy
    real(r8), dimension(:,:), intent(inout) :: sigmaim
    !> energy axis
    real(r8), dimension(:), intent(in) :: x
    !> initial residual
    real(r8), dimension(:), intent(out) :: initial_residual
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw

    real(r8), parameter :: tol_integral = 1E-12_r8
    real(r8), parameter :: derivdelta = 1E-2_r8
    integer, parameter :: derivorder = 3
    integer, parameter :: maxiter = 15
    real(r8), dimension(derivorder*2 + 1, 4) :: sc
    real(r8), dimension(derivorder*2 + 1) :: fval
    real(r8), dimension(:), allocatable :: s_xlo,s_xmax,s_xhi
    real(r8), dimension(:), allocatable :: rbuf0,rbuf1,buf_alpha
    real(r8) :: alpha,y0,y1,dy1,ddy1
    integer :: n_mode,imode,iter,i

    n_mode=size(omega)

    ! First locate some bounds on the spectral functions just so that
    ! integration becomes easier later on.
    allocate(s_xlo(n_mode))
    allocate(s_xmax(n_mode))
    allocate(s_xhi(n_mode))
    s_xlo=0.0_r8
    s_xmax=0.0_r8
    s_xhi=0.0_r8
    do imode=1,n_mode
        if ( omega(imode) .lt. lo_freqtol ) cycle
        if ( mod(imode,mw%n) .ne. mw%r ) cycle
        call lo_find_spectral_function_max_and_fwhm(omega(imode),x,sigmaim(:,imode),sigmare(:,imode),s_xmax(imode),s_xlo(imode),s_xhi(imode))
    enddo
    call mw%allreduce('sum',s_xlo)
    call mw%allreduce('sum',s_xmax)
    call mw%allreduce('sum',s_xhi)

    ! We are going to need a quadrature
    !call lo_

    allocate(rbuf0(size(x)))
    allocate(rbuf1(size(x)))
    allocate(buf_alpha(n_mode))
    rbuf0=0.0_r8
    rbuf1=0.0_r8
    buf_alpha=0.0_r8
    initial_residual=0.0_r8
    do imode=1,n_mode
        ! Skip acoustic modes
        if ( omega(imode) .lt. lo_freqtol ) cycle
        if ( mod(imode,mw%n) .ne. mw%r ) cycle

        alpha=0.0_r8
        rbuf0=sigmare(:,imode)
        do iter=1,maxiter
            ! Safeguard
            if ( iter .eq. maxiter ) then
                ! This did not work out ok, use alpha=0 as a fallback
                alpha=0.0_r8
                exit
            endif
            ! Get the function value at this point?
            rbuf1=rbuf0+alpha*x*x
            call lo_adaptive_reqursive_integral_over_spectral_function(x,omega(imode),sigmaIm(:,imode),rbuf1,s_xmax(imode),s_xlo(imode),s_xhi(imode),tol_integral,y0)
            y1=(1.0_r8-y0)**2

            ! Store the first residual?
            if ( iter .eq. 1 ) then
                initial_residual(imode)=y0
            endif


            ! We might be done?
            if ( abs(y0-1.0_r8) .lt. tol_integral*100 ) then
                buf_alpha(imode)=alpha
                exit
            endif

            ! Evaluate derivatives?
            call lo_centraldifference(3,alpha,derivdelta,sc)
            do i=1,2*derivorder+1
                rbuf1=rbuf0+sc(i,1)*x*x
                call lo_adaptive_reqursive_integral_over_spectral_function(x,omega(imode),sigmaIm(:,imode),rbuf1,s_xmax(imode),s_xlo(imode),s_xhi(imode),tol_integral,fval(i))
                fval(i)=(1.0_r8-fval(i))**2
            enddo
            dy1=sum(fval*sc(:,2))
            ddy1=sum(fval*sc(:,3))
            alpha=alpha - dy1/ddy1
        enddo
    enddo
    call mw%allreduce('sum',buf_alpha)
    call mw%allreduce('sum',initial_residual)
    ! And finally, adjust the real part of the self-energy
    do imode=1,n_mode
        sigmare(:,imode)=sigmare(:,imode)+buf_alpha(imode)*x*x
    enddo

    deallocate(s_xlo)
    deallocate(s_xmax)
    deallocate(s_xhi)
end subroutine

module subroutine tapering_function(x,y)
    !> energy axis
    real(r8), dimension(:), intent(in) :: x
    !> tapering function
    real(r8), dimension(:), intent(out) :: y

    real(r8), parameter :: taperonset = 2.7_r8/3.0_r8

    integer :: n, ii, jj, i, j
    real(r8) :: xmax, f0

    n = size(x) ! number of energy values
    xmax = x(n) ! largest x
    y = 1.0_r8
    y(1) = 0.0_r8

    ! Taper the high end with a cosine
    ii = floor(n*taperonset)
    jj = n-3
    do i = ii, jj
        j = i - ii
        f0 = real(j, r8)/real(jj - ii, r8)
        f0 = (1.0_r8 + cos(f0*lo_pi))*0.5_r8
        y(i) = f0
    end do
    do i = jj, n
        y(i) = 0.0_r8
    end do
end subroutine

!> very careful
subroutine lo_adaptive_reqursive_integral_over_spectral_function(x, omega, sigmaIm, sigmaRe, xmid, xlo, xhi, tolerance, integrated_spectralfunction)
    !> rough x-axis
    real(r8), dimension(:), intent(in) :: x
    !> harmonic frequency
    real(r8), intent(in) :: omega
    !> imaginary part
    real(r8), dimension(:), intent(in) :: sigmaIm
    !> real part
    real(r8), dimension(:), intent(in) :: sigmaRe
    !> peak location
    real(r8), intent(in) :: xmid
    !> fwhm left location
    real(r8), intent(in) :: xlo
    !> fwhm right location
    real(r8), intent(in) :: xhi
    !> relative tolerance to integrate to?
    real(r8), intent(in) :: tolerance
    !> integral of spectral function
    real(r8), intent(out) :: integrated_spectralfunction

    integer, parameter :: maxiter = 50
    integer, parameter :: maxintervals = 10000
    ! Order of quadratures to check.
    integer, parameter :: nlo = 8
    integer, parameter :: nhi = 11
    real(r8), dimension(2, nlo) :: gqlo
    real(r8), dimension(2, nhi) :: gqhi
    real(r8), dimension(nlo) :: lo_sf !, lo_sq, lo_kp
    real(r8), dimension(nhi) :: hi_sf !, hi_sq, hi_kp

    real(r8), dimension(:), allocatable :: current_nodes, next_nodes
    real(r8), dimension(:), allocatable :: current_sf !, current_sq, current_kp
    real(r8), dimension(:), allocatable :: next_sf !, next_sq, next_kp
    logical, dimension(:), allocatable :: current_ok, next_ok
    logical :: segmentok
    real(r8) :: xmax, f0, x0
    real(r8) :: int_sf !, int_sq, int_kp
    real(r8) :: err_sf !, err_sq, err_kp
    real(r8) :: vlo_sf !, vlo_sq, vlo_kp
    real(r8) :: vhi_sf !, vhi_sq, vhi_kp
    integer :: n, n_curr
    integer :: iter, i, j, ctr, ctrok, ctrnew

    n = size(x)
    xmax = x(n)

    ! Maybe try do do this recursively somehow.
    ! Think I will do it like this, I will chop things into a
    ! few intervals. Then do a nlo and nhi point Gaussian quadrature
    ! on these intervals. If within tolerance, fine, if not split
    ! it in half and go again? Can make this faster later, just
    ! proof of concept for now.

    ! Start building the first intervals.
    n_curr = 9
    allocate (current_nodes(n_curr))
    allocate (current_ok(n_curr))
    allocate (current_sf(n_curr))
    !allocate (current_sq(n_curr))
    !allocate (current_kp(n_curr))
    current_nodes = 0.0_r8
    current_ok = .false.
    current_sf = 0.0_r8
    !current_sq = 0.0_r8
    !current_kp = 0.0_r8

    ! Select some sensible starting nodes.
    current_nodes(1) = 0.0_r8
    current_nodes(2) = xlo*0.5_r8
    current_nodes(3) = xlo
    current_nodes(4) = xlo + (xmid - xlo)*0.5_r8
    current_nodes(5) = xmid
    current_nodes(6) = xmid + (xhi - xmid)*0.5_r8
    current_nodes(7) = xhi
    current_nodes(8) = xhi + (xmax - xhi)*0.5_r8
    current_nodes(9) = xmax

    ! Start adaptively integrating
    iterloop: do iter = 1, maxiter
        allocate (next_nodes(2*n_curr))
        allocate (next_ok(2*n_curr))
        allocate (next_sf(2*n_curr))
        !allocate (next_sq(2*n_curr))
        !allocate (next_kp(2*n_curr))
        next_nodes = -1.0_r8
        next_ok = .false.
        next_sf = 0.0_r8
        !next_sq = 0.0_r8
        !next_kp = 0.0_r8

        err_sf = 0.0_r8
        !err_sq = 0.0_r8
        !err_kp = 0.0_r8
        int_sf = 0.0_r8
        !int_sq = 0.0_r8
        !int_kp = 0.0_r8

        ctrok = 0
        ctrnew = 0
        ctr = 0
        do i = 1, n_curr - 1

            ! Check interval.
            if (current_ok(i)) then
                ! This interval is already fine.
                segmentok = .true.
                int_sf = int_sf + current_sf(i)
                !int_sq = int_sq + current_sq(i)
                !int_kp = int_kp + current_kp(i)
            else
                ! We have to actually check. Evaluate functions and integrate.
                call lo_gaussianquadrature(nlo, current_nodes(i), current_nodes(i + 1), gqlo)
                call lo_gaussianquadrature(nhi, current_nodes(i), current_nodes(i + 1), gqhi)
                do j = 1, nlo
                    x0 = gqlo(1, j)
                    f0 = interpolated_spectral_function(x, sigmaIm, sigmaRe, omega, x0)
                    !f1 = lo_planck(temperature, x0)
                    lo_sf(j) = f0
                    !lo_sq(j) = f0*f0*lo_pi
                    !lo_kp(j) = f0*f0*f1*(f1 + 1.0_r8)
                end do
                do j = 1, nhi
                    x0 = gqhi(1, j)
                    f0 = interpolated_spectral_function(x, sigmaIm, sigmaRe, omega, x0)
                    !f1 = lo_planck(temperature, x0)
                    hi_sf(j) = f0
                    !hi_sq(j) = f0*f0*lo_pi
                    !hi_kp(j) = f0*f0*f1*(f1 + 1.0_r8)*x0*x0
                end do
                vlo_sf = sum(lo_sf*gqlo(2, :))
                !vlo_sq = sum(lo_sq*gqlo(2, :))
                !vlo_kp = sum(lo_kp*gqlo(2, :))
                vhi_sf = sum(hi_sf*gqhi(2, :))
                !vhi_sq = sum(hi_sq*gqhi(2, :))
                !vhi_kp = sum(hi_kp*gqhi(2, :))

                ! Check if fine:
                !if ( abs(f0-f1) .lt. tolerance ) then
                if (abs(vlo_sf - vhi_sf) .lt. tolerance) then
                    current_sf(i) = vhi_sf
                    !current_sq(i) = vhi_sq
                    !current_kp(i) = vhi_kp
                    segmentok = .true.
                else
                    segmentok = .false.
                end if
                ! Accumulate things
                int_sf = int_sf + vhi_sf
                !int_sq = int_sq + vhi_sq
                !int_kp = int_kp + vhi_kp

                err_sf = err_sf + abs(vhi_sf - vlo_sf)
                !err_sq = err_sq + abs(vhi_sq - vlo_sq)
                !err_kp = err_kp + abs(vhi_kp - vlo_kp)
            end if

            if (segmentok) then
                ! Just copy
                ctr = ctr + 1
                next_nodes(ctr) = current_nodes(i)
                next_ok(ctr) = .true.
                next_sf(ctr) = current_sf(i)
                !next_sq(ctr) = current_sq(i)
                !next_kp(ctr) = current_kp(i)
                ctrok = ctrok + 1
            else
                ! Split in half? Seems sensible.
                ctr = ctr + 1
                next_nodes(ctr) = current_nodes(i)
                next_ok(ctr) = .false.
                next_sf(ctr) = 0.0_r8
                !next_sq(ctr) = 0.0_r8
                !next_kp(ctr) = 0.0_r8
                ctr = ctr + 1
                next_nodes(ctr) = (current_nodes(i) + current_nodes(i + 1))*0.5_r8
                next_ok(ctr) = .false.
                next_sf(ctr) = 0.0_r8
                !next_sq(ctr) = 0.0_r8
                !next_kp(ctr) = 0.0_r8
                ! make a note I added a segment
                ctrnew = ctrnew + 1
            end if
        end do

        !write(*,*) iter,'int',int_sf,int_kp,err_sf,ctr,n_curr,ctrnew

        ! Check if we are converged
        if (ctr + 1 .eq. n_curr) exit iterloop

        deallocate (current_nodes)
        deallocate (current_ok)
        deallocate (current_sf)
        !deallocate (current_sq)
        !deallocate (current_kp)
        n_curr = ctr + 1
        allocate (current_nodes(n_curr))
        allocate (current_ok(n_curr))
        allocate (current_sf(n_curr))
        !allocate (current_sq(n_curr))
        !allocate (current_kp(n_curr))
        current_nodes = 0.0_r8
        current_ok = .false.
        current_sf = 0.0_r8
        !current_sq = 0.0_r8
        !current_kp = 0.0_r8
        current_nodes(1:ctr) = next_nodes(1:ctr)
        current_nodes(ctr + 1) = xmax
        current_ok(1:ctr) = next_ok(1:ctr)
        current_sf(1:ctr) = next_sf(1:ctr)
        !current_sq(1:ctr) = next_sq(1:ctr)
        !current_kp(1:ctr) = next_kp(1:ctr)
        deallocate (next_nodes)
        deallocate (next_ok)
        deallocate (next_sf)
        !deallocate (next_sq)
        !deallocate (next_kp)
    end do iterloop

    integrated_spectralfunction = int_sf
    !integrated_sf_squared = int_sq
    !integrated_sf_nn = int_kp

    ! A little bit of cleanup?
    !deallocate(current_nodes)
    !deallocate(current_ok)
    !deallocate(current_sf)
    !deallocate(current_sq)
    !deallocate(current_kp)
    !deallocate(next_nodes)
    !deallocate(next_ok)
    !deallocate(next_sf)
    !deallocate(next_sq)
    !deallocate(next_kp)
end subroutine

!> find peak maximum of spectral function
subroutine lo_find_spectral_function_max_and_fwhm(omega, energy, sigmaIm, sigmaRe, energy_max, energy_lo, energy_hi)
    !> harmonic frequency
    real(r8), intent(in) :: omega
    !> energy axis
    real(r8), dimension(:), intent(in) :: energy
    !> imaginary self-energy
    real(r8), dimension(:), intent(in) :: sigmaIm
    !> real self-energy
    real(r8), dimension(:), intent(in) :: sigmaRe
    !> highest peak
    real(r8), intent(out) :: energy_max
    !> fwhm lower
    real(r8), intent(out) :: energy_lo
    !> fwhm higher
    real(r8), intent(out) :: energy_hi

    type(lo_brent_helper) :: bt
    integer, parameter :: maxiter = 200
    integer, parameter :: derivorder = 3
    real(r8), dimension(derivorder*2 + 1, 4) :: sc
    real(r8), dimension(derivorder*2 + 1) :: fval
    real(r8), dimension(:), allocatable :: y0, taperIm
    real(r8) :: x0, xlo, xhi, ylo, yhi, derivdl, peak_max, f0
    integer :: n, i, iter, ind_max, ctr, ind_left, ind_right

    ! First evaluate the spectral function. Make sure to
    ! add a small amount to the imaginary part so that
    ! we are guaranteed to get a resolveable peak.
    n = size(energy)
    allocate (y0(n))
    allocate (taperIm(n))
    y0 = 0.0_r8
    f0 = energy(2) - energy(1)
    y0 = max(sigmaIm, f0)

    call tapering_function(energy,taperIm)
    taperIm = y0*taperIm
    call evaluate_spectral_function(energy, taperIm, sigmaRe, omega, y0)

    ! Rough search to get a starting point
    x0 = -1.0_r8
    f0 = -lo_huge
    ind_max = -1
    do i = 2, n
        if (y0(i) .gt. f0) then
            f0 = y0(i)
            x0 = energy(i)
            ind_max = i
        end if
    end do

    ! Rough search to get fwhm locations to the left and to the right?
    peak_max = y0(ind_max)

    f0 = lo_huge
    ind_left = -1
    do i = 1, ind_max - 1
        if (abs(peak_max*0.5_r8 - y0(i)) .lt. f0) then
            f0 = abs(peak_max*0.5_r8 - y0(i))
            ind_left = i
        end if
    end do

    f0 = lo_huge
    ind_right = -1
    do i = ind_max + 1, n
        if (abs(peak_max*0.5_r8 - y0(i)) .lt. f0) then
            f0 = abs(peak_max*0.5_r8 - y0(i))
            ind_right = i
        end if
    end do

    ! Values that bracket the maximum, to be on the safe side.
    xlo = energy(ind_max - 1)
    xhi = energy(ind_max + 1)
    derivdl = (energy(2) - energy(1))*1E-2_r8

    ! Try a few times to see if we get something that sticks?
    initloop: do iter = 1, 10
        ! Function values (i.e. derivatives)
        call lo_centraldifference(derivorder, xlo, derivdl, sc)
        do i = 1, 2*derivorder + 1
            fval(i) = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, sc(i, 1))
        end do
        ylo = sum(fval*sc(:, 2))

        call lo_centraldifference(derivorder, xhi, derivdl, sc)
        do i = 1, 2*derivorder + 1
            fval(i) = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, sc(i, 1))
        end do
        yhi = sum(fval*sc(:, 2))

        ! If we bracket, we are good!
        if (ylo*yhi .lt. 0.0_r8) then
            exit initloop
        end if

        ! Ok, did not bracket. Try something.
        if (yhi .gt. 0.0_r8) then
            ! Move xhi a little closer?
            xhi = 0.5_r8*xhi + x0*0.5_r8
        end if
        if (ylo .lt. 0.0_r8) then
            ! Move xlo a little closer?
            xlo = 0.5_r8*xlo + x0*0.5_r8
        end if

        ! Give up after a few attempts. This happens when the peak is at zero, I think.
        ! This is a somewhat safe choice anyway.
        if (iter .ge. 10) then
            ! I will just use the rough values instead. Close enough.
            energy_max = energy(ind_max)
            if (ind_left .eq. 1) then
                energy_lo = energy_max*0.5_r8
            else
                energy_lo = energy(ind_left)
            end if
            if (ind_right .eq. n) then
                energy_hi = (energy_max + energy(n))*0.5_r8
            else
                energy_hi = energy(ind_right)
            end if
            ! return early!
            return
        end if
    end do initloop

    call bt%init(xlo, xhi, ylo, yhi, (energy(2) - energy(1))*1E-7_r8)
    x0 = -1.0_r8
    energy_max = -1.0_r8
    peak_max = -1.0_r8
    ctr = 0
    do iter = 1, maxiter
        ! Check for convergence
        if (abs(x0 - bt%next_x()) .lt. bt%tolerance) ctr = ctr + 1
        if (ctr .ge. 3) then
            ! found a zero
            energy_max = bt%current_x()
            peak_max = fval(derivorder + 1)
            exit
        end if
        ! pick next x-value
        x0 = bt%next_x()
        ! evaluate function
        call lo_centraldifference(derivorder, x0, derivdl, sc)
        do i = 1, 2*derivorder + 1
            fval(i) = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, sc(i, 1))
        end do
        f0 = sum(fval*sc(:, 2))
        ! tell the solver what happened
        call bt%update(x0, f0)
    end do

    ! Now do the same thing again, but look for the FWHM to the left.
    energy_lo = -lo_huge
    xhi = energy_max
    yhi = peak_max - peak_max*0.5_r8
    xlo = -1
    ctr = -1
    do i = n, 1, -1
        if (energy(i) .gt. energy_max - (energy(2) - energy(1))) cycle
        if (y0(i) .lt. peak_max*0.45_r8) then
            xlo = energy(i)
            exit
        end if
    end do
    if (xlo .lt. 0.0_r8) xlo = 0.0_r8
    ylo = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, xlo)
    ylo = ylo - peak_max*0.5_r8

    ! Check if something went wrong in the bracketing. Could happen with weird lineshapes.
    if (ylo*yhi .gt. 0.0_r8) then
        energy_max = energy(ind_max)
        if (ind_left .eq. 1) then
            energy_lo = energy_max*0.5_r8
        elseif (energy(ind_left) .ge. energy_max) then
            energy_lo = energy_max*0.5_r8
        else
            energy_lo = energy(ind_left)
        end if
    end if

    if (energy_lo .lt. 0.0_r8) then
        call bt%init(xlo, xhi, ylo, yhi, (energy(2) - energy(1))*1E-9_r8)
        x0 = -1.0_r8
        ctr = 0
        do iter = 1, maxiter
            ! Check for convergence
            if (abs(x0 - bt%next_x()) .lt. bt%tolerance) ctr = ctr + 1
            if (ctr .ge. 3) then
                ! found a zero
                energy_lo = bt%current_x()
                exit
            end if
            ! pick next x-value
            x0 = bt%next_x()
            ! evaluate function
            f0 = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, x0)
            f0 = f0 - peak_max*0.5_r8
            ! tell the solver what happened
            call bt%update(x0, f0)
        end do
    end if

    ! And once again, in the other direction.
    energy_hi = -lo_huge
    xlo = energy_max
    ylo = peak_max - peak_max*0.5_r8
    xhi = -1
    do i = 1, n
        if (energy(i) .lt. energy_max + energy(2) - energy(1)) cycle
        if (y0(i) .lt. peak_max*0.45_r8) then
            xhi = energy(i)
            exit
        end if
    end do
    if (xhi .lt. 0.0_r8) xhi = energy(n - 1)
    yhi = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, xhi)
    yhi = yhi - peak_max*0.5_r8
    if (ylo*yhi .gt. 0.0_r8) then
        ! This did not really work. Give the rough value instead!
        if (ind_right .eq. n) then
            energy_hi = (energy_max + energy(n))*0.5_r8
        else
            energy_hi = energy(ind_right)
        end if
    end if

    ! Only do this if bracketing worked
    if (energy_hi .lt. 0.0_r8) then
        ! Find the exact value
        call bt%init(xlo, xhi, ylo, yhi, (energy(2) - energy(1))*1E-9_r8)
        x0 = -1.0_r8
        ctr = 0
        do iter = 1, maxiter
            ! Check for convergence
            if (abs(x0 - bt%next_x()) .lt. bt%tolerance) ctr = ctr + 1
            if (ctr .ge. 3) then
                ! found a zero
                energy_hi = bt%current_x()
                exit
            end if
            ! pick next x-value
            x0 = bt%next_x()
            ! evaluate function
            f0 = interpolated_spectral_function(energy, sigmaIm, sigmaRe, omega, x0)
            f0 = f0 - peak_max*0.5_r8
            ! tell the solver what happened
            call bt%update(x0, f0)
        end do
    end if

    ! Cleanup
    deallocate (y0)
    deallocate (taperIm)
end subroutine

!> evaluate (diagonal) spectral function at arbitrary point
function interpolated_spectral_function(x, yim, yre, omega, xi) result(yi)
    !> x-axis that defines self-energy
    real(r8), dimension(:), intent(in) :: x
    !> imaginary part of self energy
    real(r8), dimension(:), intent(in) :: yim
    !> real part of self energy
    real(r8), dimension(:), intent(in) :: yre
    !> harmonic frequency
    real(r8), intent(in) :: omega
    !> point to evaluate at
    real(r8), intent(in) :: xi
    !> function value
    real(r8) :: yi

    ! Find the correct interval for interpolation.
    real(r8) :: xhi, f0, f1, f2, sigmaIm, sigmaRe
    integer :: n, ilo, ihi

    n = size(x)
    xhi = x(n)

    ! This should put us in the right interval.
    ilo = floor(n*xi/xhi) + 1
    ! ilo=max(ilo,1)
    ! ilo=min(ilo,n-1)
    ihi = ilo + 1
    ! evaluate imaginary and real self-energy
    f0 = (x(ihi) - xi)/(x(ihi) - x(ilo))
    sigmaIm = f0*yim(ilo) + (1.0_r8 - f0)*yim(ihi)
    sigmaRe = f0*yre(ilo) + (1.0_r8 - f0)*yre(ihi)

    ! get the spectral function
    f0 = 2*omega*sigmaIm*2*omega/lo_pi
    f1 = (xi**2 - omega**2 - 2*omega*sigmaRe)**2
    f2 = 4*(omega**2)*(sigmaIm**2)
    yi = f0/(f1 + f2)
end function

end submodule