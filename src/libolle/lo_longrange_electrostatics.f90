module lo_longrange_electrostatics
!!
!! Handle things related to longrange electrostatics. I had to break this out as a separate entity
!! Since I got circular dependencies all over the place.
!!
use konstanter, only: r8, i8, lo_iou, lo_pi, lo_twopi, lo_huge, lo_tol, lo_hugeint, lo_exitcode_param, lo_exitcode_symmetry
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use gottochblandat, only: walltime, tochar, lo_points_on_sphere, lo_determ, lo_invert3x3matrix, lo_sqnorm, lo_progressbar_init, lo_progressbar
use geometryfunctions, only: lo_inscribed_sphere_in_box, lo_bounding_sphere_of_box, lo_increment_dimensions
use lo_brents_method, only: lo_brent_helper
use lo_sorting, only: lo_qsort
implicit none
private

public :: lo_ewald_parameters

!> container for Ewald summation settings
type lo_ewald_parameters
    !> what is the Ewald lambda parameter
    real(r8) :: lambda = -lo_huge
    !> dielectric tensor
    real(r8), dimension(3, 3) :: eps = -lo_huge
    !> inverse dielectric tensor
    real(r8), dimension(3, 3) :: inveps = -lo_huge
    !> how many R-vectors should we sum over
    integer :: n_Rvector = -lo_hugeint
    !> how many G-vectors should we sum over
    integer :: n_Gvector = -lo_hugeint
    !> R-vectors
    real(r8), dimension(:, :), allocatable :: Rvec
    !> G-vectors
    real(r8), dimension(:, :), allocatable :: Gvec
contains
    !> set parameters
    procedure :: set => lo_set_dipole_ewald_parameters
    !> calculate the long-range dynamical matrix
    procedure :: longrange_dynamical_matrix
    !> get the long-range partial elastic constants
    procedure :: longrange_elastic_constant_bracket
    !> calculate the long-range forceconstant for a supercell
    procedure :: supercell_longrange_forceconstant
    !> measure size in memoryx
    procedure :: size_in_mem
    !> destroy
    procedure :: destroy
end type

interface
    module subroutine longrange_dynamical_matrix( &
        ew, p, q, born_effective_charges, born_onsite_correction, eps, D, Dx, Dy, Dz, reconly, chgmult)
        class(lo_ewald_parameters), intent(in) :: ew
        type(lo_crystalstructure), intent(in) :: p
        real(r8), dimension(3), intent(in) :: q
        real(r8), dimension(:, :, :), intent(in) :: born_effective_charges
        real(r8), dimension(:, :, :), intent(in) :: born_onsite_correction
        real(r8), dimension(3, 3), intent(in) :: eps
        complex(r8), dimension(:, :, :, :), intent(out) :: D
        complex(r8), dimension(:, :, :, :), intent(out), optional :: Dx, Dy, Dz
        logical, intent(in), optional :: reconly
        logical, intent(in), optional :: chgmult
    end subroutine
    module subroutine longrange_elastic_constant_bracket( ew, p, born_effective_charges, eps, bracket, rottensor, reconly, mw, verbosity)
        class(lo_ewald_parameters), intent(in) :: ew
        type(lo_crystalstructure), intent(in) :: p
        real(r8), dimension(:, :, :), intent(in) :: born_effective_charges
        real(r8), dimension(3, 3), intent(in) :: eps
        real(r8), dimension(3,3,3,3), intent(out) :: bracket
        real(r8), dimension(:,:,:,:), intent(out) :: rottensor
        logical, intent(in), optional :: reconly
        type(lo_mpi_helper), intent(inout) :: mw
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine supercell_longrange_forceconstant(ew, born_effective_charges, eps, ss, forceconstant, thres)
        class(lo_ewald_parameters), intent(in) :: ew
        real(r8), dimension(:, :, :), intent(in) :: born_effective_charges
        real(r8), dimension(3, 3), intent(in) :: eps
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), dimension(:, :, :, :), intent(out) :: forceconstant
        real(r8), intent(in) :: thres
    end subroutine
    module pure function ewald_H_thingy(x, y, inveps) result(H)
        real(r8), dimension(3), intent(in) :: x
        real(r8), intent(in) :: y
        real(r8), dimension(3, 3), intent(in) :: inveps
        real(r8), dimension(3, 3) :: H
    end function
end interface

contains

!> pick optimum Ewald parameters for dipole terms
subroutine lo_set_dipole_ewald_parameters(ew, p, eps, strategy, tol, verbosity, forced_lambda)
    !> ewald parameters
    class(lo_ewald_parameters), intent(out) :: ew
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> dielectric tensor
    real(r8), dimension(3, 3), intent(in) :: eps
    !> strategy to choose lambda
    integer, intent(in) :: strategy
    !> tolerance
    real(r8), intent(in) :: tol
    !> talk a lot
    integer, intent(in) :: verbosity
    !> force a certain lambda parameter
    real(r8), intent(in), optional :: forced_lambda

    integer, parameter :: npts = 40
    real(r8), dimension(3, npts) :: pts
    real(r8) :: t0, t1

    t0 = walltime()
    t1 = t0

    ! Set some simple things
    init: block
        ! get the points on a sphere
        call lo_points_on_sphere(pts)
        ! store the dielectric tensor, and the inverse
        ew%eps = eps
        ew%inveps = lo_invert3x3matrix(eps)
    end block init

    ! Then pick a lambda parameter
    select case (strategy)
    case (1)
        ! This means pick the optimum lambda for speed
        optlambda: block
            real(r8), parameter :: radtol = 1E-3_r8
            integer, parameter :: niter = 1000
            type(lo_brent_helper) :: bt
            real(r8) :: x0, x1, y0, y1, f0, f1, g0, g1, x, y
            integer :: iter

            ! Start by making sure we bracket the values:
            x0 = 1.0_r8
            do iter = 1, niter
                call ewald_dipole_k_r(x0, pts, ew%eps, ew%inveps, p, tol, f0, f1)
                g0 = 4*((4.0_r8*lo_pi/3.0_r8)*(f0**3))/p%volume  ! # vectors in R
                g1 = ((4.0_r8*lo_pi/3.0_r8)*(f1**3))*p%volume  ! # vectors in K
                y0 = (g0 - g1)/(g0 + g1)
                if (y0 .gt. -1E-3_r8) then
                    x0 = x0*1.5_r8
                else
                    ! done
                    exit
                end if
            end do
            x1 = 0.5_r8
            do iter = 1, niter
                call ewald_dipole_k_r(x1, pts, ew%eps, ew%inveps, p, tol, f0, f1)
                g0 = 4*((4.0_r8*lo_pi/3.0_r8)*(f0**3))/p%volume  ! # vectors in R
                g1 = ((4.0_r8*lo_pi/3.0_r8)*(f1**3))*p%volume  ! # vectors in K
                y1 = (g0 - g1)/(g0 + g1)
                if (y1 .lt. 1E-3_r8) then
                    x1 = x1/1.5_r8
                else
                    ! done
                    exit
                end if
            end do

            ! Start minimizing
            call bt%init(x0, x1, y0, y1, 1E-20_r8)
            y = bt%current_y()
            do iter = 1, niter
                ! check for convergence
                if (abs(bt%current_y()) .lt. radtol) then
                    ! found a zero
                    ew%lambda = bt%current_x()
                    exit
                end if
                ! pick next x-value
                x = bt%next_x()
                ! evaluate function
                call ewald_dipole_k_r(x, pts, ew%eps, ew%inveps, p, tol, f0, f1)
                g0 = 4*((4.0_r8*lo_pi/3.0_r8)*(f0**3))/p%volume  ! # vectors in R
                g1 = ((4.0_r8*lo_pi/3.0_r8)*(f1**3))*p%volume  ! # vectors in K
                y = (g0 - g1)/(g0 + g1)
                ! tell the solver what happened
                call bt%update(x, y)
                if (iter .eq. niter) then
                    call lo_stop_gracefully(['Could not converge Brents mehtod'], lo_exitcode_param, __FILE__, __LINE__)
                end if
            end do
        end block optlambda
    case (2)
        ! This means choose a certain range separation instead
        radlambda: block
            integer, parameter :: niter = 1000
            type(lo_brent_helper) :: bt
            real(r8), parameter :: tol_damp = 1E-5_r8
            real(r8), dimension(3) :: v0
            real(r8) :: l0, l1, y0, y1, x, y, f0, invdete, doublenndist
            integer :: iter, i

            ! Decide on a sensible distance where things should be zero, I heuristically
            ! pick double the nearest neighbour distance. Not a crucial choice, but it works.
            doublenndist = p%nearest_neighbour_distance()*2.0_r8
            ! inverse of the metric
            invdete = 1.0_r8/sqrt(lo_determ(eps))

            ! pick some number to bracket the value safely. This should be extremely safe.
            l0 = 10.0_r8
            l1 = 1E-8_r8
            ! Evaluate the function at these two values
            y0 = 0.0_r8
            y1 = 0.0_r8
            do i = 1, npts
                v0 = matmul(ew%inveps, pts(:, i)*doublenndist)
                f0 = sqrt(dot_product(v0, pts(:, i)*doublenndist))
                y0 = y0 + (l0**3)*sum(abs(ewald_H_thingy(v0*l0, f0*l0, ew%inveps)))*invdete
                y1 = y1 + (l1**3)*sum(abs(ewald_H_thingy(v0*l1, f0*l1, ew%inveps)))*invdete
            end do
            y0 = y0/npts - tol_damp
            y1 = y1/npts - tol_damp
            call bt%init(l0, l1, y0, y1, tol_damp*1E-6_r8)
            ! Iterate to minimize
            do iter = 1, niter
                ! check for convergence
                if (abs(bt%current_y()) .lt. tol_damp*1E-6_r8) then
                    ! found a zero
                    ew%lambda = bt%current_x()
                    exit
                end if
                ! pick next x-value
                x = bt%next_x()
                ! evaluate function
                y = 0.0_r8
                do i = 1, npts
                    v0 = matmul(ew%inveps, pts(:, i)*doublenndist)
                    f0 = sqrt(dot_product(v0, pts(:, i)*doublenndist))
                    y = y + (x**3)*sum(abs(ewald_H_thingy(v0*x, f0*x, ew%inveps)))*invdete
                end do
                y = y/npts - tol_damp
                ! tell the solver what happened
                call bt%update(x, y)
            end do
            ! and we have a lambda!
        end block radlambda
    case (3)
        ! Pick lambda from input
        if (present(forced_lambda)) then
            ew%lambda = forced_lambda
        else
            call lo_stop_gracefully(['Have to provide lambda'], lo_exitcode_param, __FILE__, __LINE__)
        end if
    end select

    if (verbosity .gt. 0) then
        t1 = walltime()
        write (lo_iou, *) '... decided on Ewald lambda=', tochar(ew%lambda), '  (', tochar(t1 - t0), 's)'
        t0 = t1
    end if

    ! Now build the vectors
    buildvecs: block
        integer, parameter :: niter = 1000
        integer, dimension(3) :: bd
        real(r8), dimension(:, :), allocatable :: dr
        real(r8), dimension(:), allocatable :: vn
        real(r8), dimension(3, 3) :: m0
        real(r8), dimension(3) :: v0
        real(r8) :: krad, rrad, f0, f1
        integer, dimension(:), allocatable :: di
        integer :: i, j, k, l

        ! Get the radii in real and reciprocal space
        call ewald_dipole_k_r(ew%lambda, pts, ew%eps, ew%inveps, p, tol, rrad, krad)

        ! Start with realspace, get a supercell large enough to fit the cutoff
        bd = 1
        do
            do i = 1, 3
                m0(:, i) = p%latticevectors(:, i)*(2*bd(i) + 1)
            end do
            f0 = lo_inscribed_sphere_in_box(m0)
            if (f0 .gt. rrad) exit
            bd = lo_increment_dimensions(bd, p%latticevectors)
        end do
        ! count realspace vectors
        f0 = rrad**2
        l = 0
        do i = -bd(1), bd(1)
        do j = -bd(2), bd(2)
        do k = -bd(3), bd(3)
            v0 = p%fractional_to_cartesian(real([i, j, k], r8))
            if (lo_sqnorm(v0) .lt. f0) then
                l = l + 1
            end if
        end do
        end do
        end do
        ew%n_Rvector = l
        allocate (ew%Rvec(3, l))
        allocate (dr(3, l))
        allocate (vn(l))
        allocate (di(l))
        ew%rvec = 0.0_r8
        dr = 0.0_r8
        vn = 0.0_r8
        di = 0
        l = 0
        do i = -bd(1), bd(1)
        do j = -bd(2), bd(2)
        do k = -bd(3), bd(3)
            v0 = p%fractional_to_cartesian(real([i, j, k], r8))
            f1 = lo_sqnorm(v0)
            if (f1 .lt. f0) then
                l = l + 1
                vn(l) = -f1
                dr(:, l) = v0
            end if
        end do
        end do
        end do
        ! Sort backwards by length
        call lo_qsort(vn, di)
        do i = 1, ew%n_Rvector
            ew%Rvec(:, i) = dr(:, di(i))
        end do
        deallocate (dr)
        deallocate (vn)
        deallocate (di)

        ! Then exactly the same thing again for G-vectors
        bd = 1
        do
            do i = 1, 3
                m0(:, i) = p%reciprocal_latticevectors(:, i)*(2*bd(i) + 1)
            end do
            f0 = lo_inscribed_sphere_in_box(m0)
            if (f0 .gt. krad) exit
            bd = lo_increment_dimensions(bd, p%reciprocal_latticevectors)
        end do
        ! count realspace vectors
        f0 = krad**2
        l = 0
        do i = -bd(1), bd(1)
        do j = -bd(2), bd(2)
        do k = -bd(3), bd(3)
            v0 = p%fractional_to_cartesian(real([i, j, k], r8), reciprocal=.true.)
            if (lo_sqnorm(v0) .lt. f0) then
                l = l + 1
            end if
        end do
        end do
        end do
        ew%n_Gvector = l
        allocate (ew%Gvec(3, l))
        allocate (dr(3, l))
        allocate (vn(l))
        allocate (di(l))
        ew%Gvec = 0.0_r8
        dr = 0.0_r8
        vn = 0.0_r8
        di = 0
        l = 0
        do i = -bd(1), bd(1)
        do j = -bd(2), bd(2)
        do k = -bd(3), bd(3)
            v0 = p%fractional_to_cartesian(real([i, j, k], r8), reciprocal=.true.)
            f1 = lo_sqnorm(v0)
            if (f1 .lt. f0) then
                l = l + 1
                vn(l) = -f1
                dr(:, l) = v0
            end if
        end do
        end do
        end do
        ! Sort backwards by length
        call lo_qsort(vn, di)
        do i = 1, ew%n_Gvector
            ew%Gvec(:, i) = dr(:, di(i))
        end do
        deallocate (dr)
        deallocate (vn)
        deallocate (di)
    end block buildvecs

    if (verbosity .gt. 0) then
        t1 = walltime()
        write (lo_iou, *) '... stored ', tochar(ew%n_Rvector), '+', tochar(ew%n_Gvector), &
            ' lattice vectors (', tochar(t1 - t0), 's)'
        t0 = t1
    end if
end subroutine

!> turn the ratio of k-space volume to r-space volume into a function to be minimized
subroutine ewald_dipole_k_r(lambda, pts, eps, inveps, p, tol, rrad, krad)
    !> ewald parameter
    real(r8), intent(in) :: lambda
    !> points to test with
    real(r8), dimension(:, :), intent(in) :: pts
    !> dielectric tensor
    real(r8), dimension(3, 3), intent(in) :: eps, inveps
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> tolerance to target
    real(r8), intent(in) :: tol
    !> radius in r-space
    real(r8), intent(out) :: rrad
    !> radius in k-space
    real(r8), intent(out) :: krad

    integer, parameter :: niter = 1000
    type(lo_brent_helper) :: bt

    ! Init to something
    krad = -1.0_r8
    rrad = -1.0_r8

    ! find radius in K-space
    findkrad: block
        real(r8), dimension(3) :: v0, v1
        real(r8) :: r0, r1, y0, y1, f0, f1, x, y, inv4lambdasq
        integer :: iter, i, j, k

        f0 = lo_inscribed_sphere_in_box(p%reciprocal_latticevectors)
        inv4lambdasq = 1.0_r8/(4.0_r8*lambda*lambda)
        ! starting guesses
        r0 = f0*1E-9_r8
        r1 = f0*1E9_r8
        y0 = 0.0_r8
        y1 = 0.0_r8
        do i = 1, size(pts, 2)
            v0 = pts(:, i)*r0*lo_twopi
            v1 = pts(:, i)*r1*lo_twopi
            f0 = 0.0_r8
            f1 = 0.0_r8
            do j = 1, 3
            do k = 1, 3
                f0 = f0 + v0(k)*eps(k, j)*v0(j)
                f1 = f1 + v1(k)*eps(k, j)*v1(j)
            end do
            end do
            y0 = y0 + exp(-f0*inv4lambdasq)*f0
            y1 = y1 + exp(-f1*inv4lambdasq)*f1
        end do
        y0 = y0 - tol
        y1 = y1 - tol
        ! Make sure we bracket the target
        if (y0*y1 .gt. 0.0_r8) then
            call lo_stop_gracefully(['Could not bracket values'], lo_exitcode_param, __FILE__, __LINE__)
        end if

        call bt%init(r0, r1, y0, y1, tol*1E-3_r8)
        do iter = 1, niter
            ! check for convergence
            if (abs(bt%current_y()) .lt. tol*1E-3_r8) then
                ! found a zero
                krad = bt%current_x()
                exit
            end if
            ! pick next x-value
            x = bt%next_x()
            ! evaluate function
            y = 0.0_r8
            do i = 1, size(pts, 2)
                v0 = pts(:, i)*x*lo_twopi
                f0 = 0.0_r8
                do j = 1, 3
                do k = 1, 3
                    f0 = f0 + v0(k)*eps(k, j)*v0(j)
                end do
                end do
                y = y + exp(-f0*inv4lambdasq)*f0
            end do
            y = y - tol
            ! tell the solver what happened
            call bt%update(x, y)
            if (iter .eq. niter) then
                call lo_stop_gracefully(['Could not converge Brents mehtod'], lo_exitcode_param, __FILE__, __LINE__)
            end if
        end do
        ! then add the bounding sphere to be on the safe side, the longest possible q-vector added later.
        krad = krad + lo_bounding_sphere_of_box(p%reciprocal_latticevectors)
    end block findkrad

    ! Get the radius in R-space, same thing
    findrrad: block
        real(r8), dimension(3) :: v0, v1
        real(r8) :: r0, r1, y0, y1, x, y, f0, f1, invdete
        integer :: iter, i

        invdete = 1.0_r8/sqrt(lo_determ(eps))
        f0 = lo_inscribed_sphere_in_box(p%latticevectors)
        ! This should bracket any sensible choice
        r0 = f0*1E-4_r8
        r1 = f0*1E4_r8
        ! Get the function values
        y0 = 0.0_r8
        y1 = 0.0_r8
        do i = 1, size(pts, 2)
            v0 = matmul(inveps, pts(:, i)*r0)
            v1 = matmul(inveps, pts(:, i)*r1)
            f0 = sqrt(dot_product(v0, pts(:, i)*r0))
            f1 = sqrt(dot_product(v1, pts(:, i)*r1))
            y0 = y0 + (lambda**3)*sum(abs(ewald_H_thingy(v0*lambda, f0*lambda, inveps)))*invdete
            y1 = y1 + (lambda**3)*sum(abs(ewald_H_thingy(v1*lambda, f1*lambda, inveps)))*invdete
        end do
        y0 = y0 - tol
        y1 = y1 - tol
        ! Be sure we bracket the target
        if (y0*y1 .gt. 0.0_r8) then
            call lo_stop_gracefully(['Could not bracket values'], lo_exitcode_param, __FILE__, __LINE__)
        end if
        call bt%init(r0, r1, y0, y1, tol*1E-3_r8)
        do iter = 1, niter
            ! check for convergence
            if (abs(bt%current_y()) .lt. tol*1E-3_r8) then
                ! found a zero
                rrad = bt%current_x()
                exit
            end if
            ! pick next x-value
            x = bt%next_x()
            ! evaluate function
            y = 0.0_r8
            do i = 1, size(pts, 2)
                v0 = matmul(inveps, pts(:, i)*x)
                f0 = sqrt(dot_product(v0, pts(:, i)*x))
                y = y + (lambda**3)*sum(abs(ewald_H_thingy(v0*lambda, f0*lambda, inveps)))*invdete
            end do
            y = y - tol
            ! tell the solver what happened
            call bt%update(x, y)
            if (iter .eq. niter) then
                call lo_stop_gracefully(['Could not converge Brents mehtod'], lo_exitcode_param, __FILE__, __LINE__)
            end if
        end do
        ! then add the bounding sphere, the longest possible rij-vector added later.
        rrad = rrad + lo_bounding_sphere_of_box(p%latticevectors)
    end block findrrad
end subroutine

!> measure size in memory, in bytes
function size_in_mem(ew) result(mem)
    !> ewald settings
    class(lo_ewald_parameters), intent(in) :: ew
    !> size in memory
    integer(i8) :: mem

    mem = 0
    mem = mem + storage_size(ew)
    if (allocated(ew%Rvec)) mem = mem + storage_size(ew%Rvec)*size(ew%Rvec)
    if (allocated(ew%Gvec)) mem = mem + storage_size(ew%Gvec)*size(ew%Gvec)
    mem = mem/8
end function

!> destroy
subroutine destroy(ew)
    !> ewald settings
    class(lo_ewald_parameters), intent(inout) :: ew

    if (allocated(ew%Rvec)) deallocate (ew%Rvec)
    if (allocated(ew%Gvec)) deallocate (ew%Gvec)
    ew%lambda = -lo_huge
    ew%eps = -lo_huge
    ew%inveps = -lo_huge
    ew%n_Rvector = -lo_hugeint
    ew%n_Gvector = -lo_hugeint
end subroutine

end module
