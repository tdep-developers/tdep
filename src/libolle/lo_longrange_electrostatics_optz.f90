submodule(lo_longrange_electrostatics) lo_longrange_electrostatics_optz
!!
!! Do some nonlinear magic to make sure that the Born effective charges don't
!! break the Hermiticity of the dynamical matrix when used later on to correct
!! phonons.
!!
use konstanter, only: lo_sqtol
use gottochblandat, only: lo_real_symmetric_eigenvalues_eigenvectors, lo_linear_least_squares, lo_unflatten_2tensor
use quadratures_stencils, only: lo_centraldifference
implicit none

!> helper to keep track of nonlinear minimization. When needed I will break this out to it's own thing, it's kinda neat.
type lo_nonlinear_minimization_helper
    !> which iteration are we on?
    integer :: iter = -lo_hugeint
    !> how many variables
    integer :: nx = -lo_hugeint
    !> how many directions to provide
    integer :: n_direction = -lo_hugeint
    !> which algo did I use this step?
    integer :: algo = -lo_hugeint
    !> Errors (the thing we want to minimize)
    real(r8) :: e_curr, e_prev
    !> Current x
    real(r8), dimension(:), allocatable :: x_curr
    !> Previous x
    real(r8), dimension(:), allocatable :: x_prev
    !> Current gradient
    real(r8), dimension(:), allocatable :: g_curr
    !> Previous gradient
    real(r8), dimension(:), allocatable :: g_prev

    !> Search direction things
    real(r8), dimension(:, :), allocatable :: sdir
    real(r8), dimension(:, :), allocatable :: sx
    real(r8), dimension(:), allocatable :: serr
    real(r8), dimension(:), allocatable :: sstep

    !> BFGS helpers
    real(r8), dimension(:), allocatable :: bfgs_oldx
    real(r8), dimension(:), allocatable :: bfgs_oldg
    real(r8), dimension(:), allocatable :: bfgs_s
    real(r8), dimension(:), allocatable :: bfgs_y
    real(r8), dimension(:, :), allocatable :: bfgs_B
    integer :: bfgs_n_update = -lo_hugeint

    !> CG helpers
    real(r8), dimension(:), allocatable :: cg_oldx
    real(r8), dimension(:), allocatable :: cg_oldg
contains
    !> initialize the solver
    procedure :: init => init_nonlinear_helper
    !> update search directions
    procedure :: searchdir => get_searchdir
    !> update derivative delta
    procedure :: get_delta
    !> update initial step size
    procedure :: get_step
    !> update x (can not type bind this in ifort 2019???? WAT?)
    !procedure :: next_x
end type

contains

!> force Born charges to become Hermitian. Serial for now.
module subroutine force_borncharges_Hermitian(ew, x_Z, coeff_Z, eps, p, verbosity)
    !> Ewald summation helper
    class(lo_ewald_parameters), intent(in) :: ew
    !> irreducible representation of Born charges
    real(r8), dimension(:), intent(inout) :: x_Z
    !> coefficient matrix that can transform the irreducible representation to actual Born charges.
    real(r8), dimension(:, :), intent(in) :: coeff_Z
    !> dielectric constant
    real(r8), dimension(3, 3), intent(in) :: eps
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! Define some magical parameters for the nonlinear solver
    real(r8), parameter :: errtol = 1E-28_r8
    integer, parameter :: derivorder = 3

    ! Normal temporary things
    real(r8), dimension(:, :, :, :), allocatable :: wD, dynmat_reciprocal
    real(r8), dimension(:, :, :), allocatable :: wZ
    real(r8), dimension(:), allocatable :: wf
    real(r8) :: derivdelta
    integer :: nx

    ! Set some sensible things
    init: block
        complex(r8), dimension(:, :, :, :), allocatable :: D0
        real(r8), dimension(:, :, :), allocatable :: w0
        integer :: a1
        ! Number of variables to optimize for
        nx = size(x_Z)

        ! Workspace for born charges, tensor form
        allocate (wZ(3, 3, p%na))
        wZ = 0.0_r8
        ! Workspace for born charges, flattened
        allocate (wf(9*p%na))
        wf = 0.0_r8
        ! Unflatten the initial born charges
        !@TODO think about how annoying this is to change such that it is automagical?
        call dgemv('N', 9*p%na, nx, 1.0_r8, coeff_Z, 9*p%na, x_z, 1, 0.0_r8, wf, 1)
        do a1 = 1, p%na
            wZ(:, :, a1) = lo_unflatten_2tensor(wf((a1 - 1)*9 + 1:a1*9))
        end do
        ! Get the dynamical matrix at Gamma, but without any charges multiplied in.
        ! Also create some workspace for the error and gradient functions.
        allocate (D0(3, 3, p%na, p%na))
        allocate (wD(3, 3, p%na, p%na))
        allocate (dynmat_reciprocal(3, 3, p%na, p%na))
        allocate (w0(3, 3, p%na))
        D0 = 0.0_r8
        wD = 0.0_r8
        w0 = 0.0_r8
        dynmat_reciprocal = 0.0_r8
        call ew%longrange_dynamical_matrix(p, [0.0_r8, 0.0_r8, 0.0_r8], wZ, w0, eps, D0, reconly=.true., chgmult=.false.)
        ! Make sure it's real. Has to be.
        if (sum(abs(aimag(D0))) .gt. lo_sqtol) then
            call lo_stop_gracefully(['Dipole dynamical matrix at Gamma is not real'], lo_exitcode_symmetry, __FILE__, __LINE__)
        else
            ! keep only the real part to work with, slightly faster.
            dynmat_reciprocal = real(D0, r8)
        end if
        ! cleanup
        deallocate (D0)
        deallocate (w0)
    end block init

    ! Figure out settings for the nonlinear solver. Also possibly kill early, in
    ! case things are already fine.
    preopt: block
        real(r8) :: err0

        ! Calculate the current error?
        err0 = Z_error_function(coeff_Z, x_Z, wf, wZ, wD, dynmat_reciprocal)
        ! We could be fine?
        if (err0 .lt. errtol) then
            if (verbosity .gt. 0) then
                write (lo_iou, *) '... Born charges already Hermitian'
            end if
            return
        end if

        ! We were not fine. Too bad.
        write (lo_iou, *) '... Born charges not Hermitian, error:', err0
        write (lo_iou, *) '... DISPERSIONS MIGHT BE IRREGULAR!'
        write (lo_iou, *) '... FIXME!'

        ! fkdev: we currently don't enfore anything until we now how to do this properly
        return
    end block preopt

    ! Apparently we have to optimize for reals. I hate this. 50% of the time, it works 100% of the time.
    opt: block
        ! First some parameters for the minimization
        integer, parameter :: nouter = 1000
        integer, parameter :: ninner = 5000
        integer, parameter :: maxctr = 25
        real(r8), parameter :: gradtol = 1E-20_r8
        real(r8), parameter :: longtime = 1.0_r8 ! if I waited this long, print how it's going.
        ! Things needed to set stuff
        type(lo_nonlinear_minimization_helper) :: mh
        real(r8), dimension(:), allocatable :: x, x1, x2, grad, searchdir
        real(r8) :: err0, err1, err2, step, df0
        real(r8) :: t0, t1
        integer :: outiter, initer, ctr, ictr
        integer :: idirection

        ! start timer
        t0 = walltime()
        t1 = t0

        ! Some workspace
        allocate (x(nx))
        allocate (x1(nx))
        allocate (x2(nx))
        allocate (grad(nx))
        allocate (searchdir(nx))
        x = 0.0_r8
        x1 = 0.0_r8
        x2 = 0.0_r8
        grad = 0.0_r8
        searchdir = 0.0_r8
        ! Start guessing things for the solver
        derivdelta = 1E-6_r8  ! or some other random number
        x = x_Z               ! starting point
        df0 = lo_huge         ! for output
        ! Calculate error and gradient
        err0 = Z_error_function(coeff_Z, x, wf, wZ, wD, dynmat_reciprocal)
        !derivdelta=mh%get_delta(err0,1.0_r8)
        ! Just something to guess with.
        derivdelta = 1E-2_r8
        call Z_error_gradient(coeff_Z, x, derivorder, derivdelta, grad, wf, wZ, wD, dynmat_reciprocal)
        ! Initialize solver
        call mh%init(nx, x, grad, initial_hessian=10.0_r8)

        if (verbosity .gt. 0) then
            write (lo_iou, *) '... optimizing Born charges to remove error'
            write (lo_iou, "(5X,A10,3X,A14,3X,A14)") 'iter', 'err', 'grad'
        end if

        ictr = 0
        outloop: do outiter = 1, nouter
            ! Get the current error
            err0 = Z_error_function(coeff_Z, x, wf, wZ, wD, dynmat_reciprocal)
            ! Update the derivative, if we have enough updates of the Hessian
            derivdelta = min(derivdelta, mh%get_delta(err0, 1E4_r8))
            call Z_error_gradient(coeff_Z, x, derivorder, derivdelta, grad, wf, wZ, wD, dynmat_reciprocal)

            ! Check if we are converged
            if (err0 .lt. errtol) then
                ! Converged
                exit outloop
            end if
            if (norm2(grad) .lt. gradtol) then
                ! Converged, I think
                exit outloop
            end if
            ! Get the search direction
            call mh%searchdir(x, grad, err0)

            ! Do a line search in this direction?
            do idirection = 1, mh%n_direction
                x1 = mh%x_curr
                err0 = mh%e_curr
                err2 = lo_huge
                ctr = 0
                !step=1.0_r8
                searchdir = mh%sdir(:, idirection)
                searchdir = searchdir/norm2(searchdir)
                ! pick a sensible step size?
                step = mh%get_step(searchdir, grad, err0) !/1.5_r8**3
                inloop: do initer = 1, ninner
                    ictr = ictr + 1
                    ! take a step
                    x1 = x1 + step*searchdir
                    ! evaluate error
                    err1 = Z_error_function(coeff_Z, x1, wf, wZ, wD, dynmat_reciprocal)
                    if (err1 .lt. errtol) then
                        ! Could happen that we are really converged here!
                        x = x1
                        err0 = err1
                        exit outloop
                    elseif (err1 .lt. err0) then
                        ! Good step, keep going
                        step = step*1.5_r8
                        ! If it was better than the previous best, keep it!
                        if (err1 .lt. err2) then
                            err2 = err1
                            x2 = x1
                        end if
                    else
                        ! Bad step, backstep and decrease the step size
                        x1 = x1 - step*searchdir
                        step = step*0.25_r8
                        ctr = ctr + 1
                    end if
                    ! Start over with new direction if we tried too many times
                    if (ctr .ge. maxctr) exit inloop

                    if (verbosity .gt. 0) then
                    if (log(err0/df0) .lt. -2.0_r8) then
                        write (lo_iou, "(5X,I10,3X,E14.5,3X,E14.5)") ictr, err0, norm2(grad)
                        df0 = err0
                    end if
                    end if
                end do inloop
                mh%sx(:, idirection) = x2
                mh%serr(idirection) = err2
            end do

            ! Pick the next x
            !call mh%next_x(x,err1)
            call next_x(mh, x, err1)

            ! If we've been fiddling for a long time, report?
            if (verbosity .gt. 0) then
                t1 = walltime()
                if (t1 - t0 .gt. longtime) then
                    t0 = t1
                    write (lo_iou, "(5X,I10,3X,E14.5,3X,E14.5)") ictr, err0, norm2(grad)
                end if
            end if
        end do outloop

        ! Some cleanup
        if (verbosity .gt. 0) then
            write (lo_iou, *) '... this is the best we could do in finite time:'
            write (lo_iou, "(5X,I10,3X,E14.5,3X,E14.5)") ictr, err0, norm2(grad)
        end if

        ! Seems we have fixed the Born charges! Return the neat ones
        x_Z = x
        ! And cleanup
        deallocate (x)
        deallocate (x1)
        deallocate (x2)
        deallocate (grad)
        deallocate (searchdir)
        deallocate (wD)
        deallocate (dynmat_reciprocal)
        deallocate (wZ)
        deallocate (wf)
    end block opt
end subroutine

! below are not exposed

!> gradient for fixing the Born charges
subroutine Z_error_gradient(coeff, x, order, delta, grad, wf, wZ, wD, dm)
    !> coefficient matrix to convert from irreducible
    real(r8), dimension(:, :), intent(in) :: coeff
    !> irreducible
    real(r8), dimension(:), intent(in) :: x
    !> order of stencil
    integer, intent(in) :: order
    !> delta in stencil
    real(r8), intent(in) :: delta
    !> gradient
    real(r8), dimension(:), intent(out) :: grad
    !> workspace
    real(r8), dimension(:), intent(inout) :: wf
    real(r8), dimension(:, :, :), intent(inout) :: wZ
    real(r8), dimension(:, :, :, :), intent(inout) :: wD
    real(r8), dimension(:, :, :, :), intent(in) :: dm

    real(r8), dimension(2*order + 1, 4) :: sc
    real(r8), dimension(2*order + 1) :: vl
    real(r8), dimension(size(x)) :: wx
    integer :: i, j

    do i = 1, size(x)
        call lo_centraldifference(order, x(i), delta, sc)
        do j = 1, 2*order + 1
            wx = x
            wx(i) = sc(j, 1)
            vl(j) = Z_error_function(coeff, wx, wf, wZ, wD, dm)
        end do
        grad(i) = sum(sc(:, 2)*vl)
    end do
end subroutine

!> Calculate the squared violation of hermiticity
function Z_error_function(coeff, x, wf, wZ, wD, dm) result(hermerr)
    !> coefficient matrix to convert from irreducible
    real(r8), dimension(:, :), intent(in) :: coeff
    !> irreducible
    real(r8), dimension(:), intent(in) :: x
    ! work arrays
    real(r8), dimension(:), intent(inout) :: wf
    real(r8), dimension(:, :, :), intent(inout) :: wZ
    real(r8), dimension(:, :, :, :), intent(inout) :: wD
    real(r8), dimension(:, :, :, :), intent(in) :: dm
    ! squared error
    real(r8) :: hermerr

    real(r8), dimension(3, 3) :: m0
    integer :: a1, a2, i, j, ii, jj

    ! Unflatten the Born charges
    wf = matmul(coeff, x)
    do a1 = 1, size(wZ, 3)
        wZ(:, :, a1) = lo_unflatten_2tensor(wf((a1 - 1)*9 + 1:a1*9))
    end do

    ! Multiply in the charges. This can probably be done smarter.
    wD = 0.0_r8
    do a2 = 1, size(wZ, 3)
    do a1 = 1, size(wZ, 3)
        do j = 1, 3
        do i = 1, 3
            do jj = 1, 3
            do ii = 1, 3
                wD(i, j, a1, a2) = wD(i, j, a1, a2) + wZ(i, ii, a1)*wZ(j, jj, a2)*dm(ii, jj, a1, a2)
            end do
            end do
        end do
        end do
    end do
    end do

    hermerr = 0.0_r8
    do a1 = 1, size(wZ, 3)
        m0 = 0.0_r8
        do a2 = 1, size(wZ, 3)
            m0 = m0 + wD(:, :, a1, a2)
        end do
        hermerr = hermerr + sum(abs(m0 - transpose(m0)))**2
    end do
    hermerr = hermerr/((size(wZ, 3)*3)**2)
end function

!> initialize the nonlinear helper
subroutine init_nonlinear_helper(mh, nx, x0, g0, initial_hessian)
    !> helper
    class(lo_nonlinear_minimization_helper), intent(out) :: mh
    !> number of variables
    integer, intent(in) :: nx
    !> starting x
    real(r8), dimension(:), intent(in) :: x0
    !> starting gradient
    real(r8), dimension(:), intent(in) :: g0
    !> starting Hessian
    real(r8), intent(in) :: initial_hessian

    integer :: i

    mh%iter = 0
    mh%nx = nx
    mh%n_direction = 3 ! how many different search directions to provide
    allocate (mh%sdir(nx, mh%n_direction))
    allocate (mh%sx(nx, mh%n_direction))
    allocate (mh%serr(mh%n_direction))
    mh%sdir = 0.0_r8
    mh%sx = 0.0_r8
    mh%serr = 0.0_r8

    allocate (mh%x_curr(mh%nx))
    allocate (mh%x_prev(mh%nx))
    allocate (mh%g_curr(mh%nx))
    allocate (mh%g_prev(mh%nx))
    mh%x_curr = x0
    mh%x_prev = x0
    mh%g_prev = g0
    mh%g_prev = g0

    ! BFGS helper things
    allocate (mh%bfgs_oldg(nx))
    allocate (mh%bfgs_oldx(nx))
    allocate (mh%bfgs_s(nx))
    allocate (mh%bfgs_y(nx))
    allocate (mh%bfgs_B(nx, nx))
    mh%bfgs_oldg = g0
    mh%bfgs_oldx = x0
    mh%bfgs_s = 0.0_r8
    mh%bfgs_y = 0.0_r8
    mh%bfgs_B = 0.0_r8
    do i = 1, nx
        mh%bfgs_B(i, i) = initial_hessian
    end do
    mh%bfgs_n_update = 0

    ! CG helper things
    allocate (mh%cg_oldg(nx))
    allocate (mh%cg_oldx(nx))
    mh%cg_oldg = g0
    mh%cg_oldx = g0
end subroutine

!> get a sensible step size to seed the linear search?
function get_step(mh, searchdir, gradient, y) result(step)
    !> helper
    class(lo_nonlinear_minimization_helper), intent(in) :: mh
    !> direction to look
    real(r8), dimension(:), intent(in) :: searchdir
    !> gradient at this point
    real(r8), dimension(:), intent(in) :: gradient
    !> function value at this point
    real(r8), intent(in) :: y
    !> step size
    real(r8) :: step

    real(r8) :: A, B, f0

    ! Try a few different things?
    A = dot_product(searchdir, matmul(mh%bfgs_B, searchdir))*0.5_r8
    B = dot_product(searchdir, gradient)

    f0 = B**2 - 4*A*y
    ! try to minimize the quadratic equation?
    if (f0 .gt. 1E-30_r8) then
        ! maybe use quadratic equation?
        if (sqrt(f0) .gt. b) then
            ! yup, makes sense
            step = (sqrt(f0) - B)/(2*A)
            return
        end if
    end if

    ! quadratic no good. Linear?
    if (B .lt. -1E-15_r8) then
        step = -y/B
        return
    end if

    ! nothing worked. Hmm. Pick 1.
    step = 1.0_r8
end function

!> get a sensible derivative delta?
function get_delta(mh, y, prefactor) result(delta)
    !> helper
    class(lo_nonlinear_minimization_helper), intent(in) :: mh
    !> current function value
    real(r8), intent(in) :: y
    !> scaling factor
    real(r8), intent(in) :: prefactor
    !> stencil delta
    real(r8) :: delta

    real(r8) :: f0

    ! Norm of the Hessian or something like that
    f0 = norm2(mh%bfgs_B)/(mh%nx**2)
    ! Probably need some tolerance here?
    if (f0 .gt. 1e-15_r8) then
        delta = prefactor*sqrt(y/f0)
    else
        delta = prefactor*sqrt(y/1E-15_r8)
    end if
end function

!> get the next point to look at
subroutine next_x(mh, x1, e1)
    !> helper
    class(lo_nonlinear_minimization_helper), intent(inout) :: mh
    !> next position
    real(r8), dimension(:), allocatable, intent(inout) :: x1
    !> lowest error
    real(r8), intent(out) :: e1

    real(r8) :: f0
    integer :: i

    f0 = lo_huge
    do i = 1, mh%n_direction
        if (mh%serr(i) .lt. f0) then
            f0 = mh%serr(i)
            x1 = mh%sx(:, i)
            e1 = f0
            mh%algo = i
        end if
    end do
end subroutine

!> fetch search directions (yes several candidates)
subroutine get_searchdir(mh, x1, g1, e1)
    !> helper
    class(lo_nonlinear_minimization_helper), intent(inout) :: mh
    !> updated positions
    real(r8), dimension(:), intent(in) :: x1
    !> updated gradient
    real(r8), dimension(:), intent(in) :: g1
    !> updated error
    real(r8), intent(in) :: e1

    real(r8), parameter :: betatol = 1E-20_r8
    real(r8), parameter :: eigvaltol = 1E-10_r8
    real(r8) :: beta1, beta2
    integer :: i, j

    ! Update where we are?
    mh%x_prev = mh%x_curr
    mh%g_prev = mh%g_curr
    mh%e_prev = mh%e_curr
    mh%x_curr = x1
    mh%g_curr = g1
    mh%e_curr = e1

    if (mh%iter .eq. 0) then
        ! If it's the first step, all the directions will be steepest descent
        do i = 1, mh%n_direction
            mh%sdir(:, i) = -g1
        end do
    else
        ! This is a normal step, after the first.
        ! Try to update the search directions in fancy ways. First it's just
        ! steepest descent:
        steepdesc: block
            mh%sdir(:, 1) = -g1
        end block steepdesc

        ! Second one is BFGS
        bfgs: block
            real(r8), dimension(:, :), allocatable :: wB, eigvec
            real(r8), dimension(:), allocatable :: eigval

            allocate (eigvec(mh%nx, mh%nx))
            allocate (eigval(mh%nx))
            allocate (wB(mh%nx, mh%nx))
            eigvec = 0.0_r8
            eigval = 0.0_r8
            wB = 0.0_r8

            mh%bfgs_y = mh%bfgs_oldg - g1
            mh%bfgs_s = mh%bfgs_oldx - x1
            ! Start with a straight update
            beta1 = dot_product(mh%bfgs_y, mh%bfgs_s)
            beta2 = dot_product(mh%bfgs_s, matmul(mh%bfgs_B, mh%bfgs_s))
            if (abs(beta1) .gt. betatol .and. beta2 .gt. betatol) then
                ! Ok numerics seem fine, try BFGS update
                mh%bfgs_s = matmul(mh%bfgs_B, mh%bfgs_s)
                wB = 0.0_r8
                do i = 1, mh%nx
                do j = 1, mh%nx
                    wB(i, j) = mh%bfgs_B(j, i) + mh%bfgs_y(i)*mh%bfgs_y(j)/beta1 - mh%bfgs_s(i)*mh%bfgs_s(j)/beta2
                end do
                end do
                ! Check the eigenvalues
                call lo_real_symmetric_eigenvalues_eigenvectors(wB, eigval, eigvec)
                if (minval(eigval) .gt. eigvaltol) then
                    ! This was a proper update, update everything
                    mh%bfgs_B = wB
                    mh%bfgs_oldg = g1
                    mh%bfgs_oldx = x1
                    mh%bfgs_n_update = mh%bfgs_n_update + 1
                end if
            end if

            ! Anyhow, solve for the search direction
            call lo_linear_least_squares(mh%bfgs_B, -g1, mh%sdir(:, 2))

            deallocate (eigvec)
            deallocate (eigval)
            deallocate (wB)
        end block bfgs

        ! Third one is conjugate gradient
        cg: block
            real(r8) :: f0, f1

            ! Start with Fletcher-Reeves?
            f0 = dot_product(g1, g1)
            f1 = dot_product(mh%cg_oldg, mh%cg_oldg)
            if (f1 .gt. betatol) then
                mh%sdir(:, 3) = -g1 - (f0/f1)*mh%cg_oldg
                mh%cg_oldg = g1
                mh%cg_oldx = x1
            else
                mh%sdir(:, 3) = -g1
            end if
            ! Maybe something else when things are annoying.
        end block cg
    end if

    ! And update the number of iterations
    mh%iter = mh%iter + 1
end subroutine

end submodule
