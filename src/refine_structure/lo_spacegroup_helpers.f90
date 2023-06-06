submodule(lo_spacegroup) lo_spacegroup_helpers
!! Helper routines for symmetry stuff
implicit none
contains

!> Delaunay reduction of lattice axistors to canonical form
module subroutine delaunay_reduce_basis(basis, canonical_basis)
    !> original lattice axistors
    real(r8), dimension(3, 3), intent(in) :: basis
    !> canonical lattice axistors
    real(r8), dimension(3, 3), intent(out) :: canonical_basis

    real(r8), dimension(3, 4) :: ext_basis
    real(r8), dimension(3, 7) :: dum
    real(r8), dimension(7) :: d
    real(r8) :: f0
    integer, dimension(7) :: ind
    integer :: i, j, k

    ! First do the Delaunay reduction thing
    ext_basis = 0.0_r8
    ext_basis(:, 1:3) = basis
    do i = 1, 3
        ext_basis(:, 4) = ext_basis(:, 4) - basis(:, i)
    end do
    ! Iterate and stuff
    iterloop: do
        do i = 1, 4
        do j = i + 1, 4
            ! check scalar product
            f0 = dot_product(ext_basis(:, i), ext_basis(:, j))
            if (f0 .gt. lo_sqtol) then
                ! subtract i from all not i and j
                do k = 1, 4
                    if (k .ne. i .and. k .ne. j) then
                        ext_basis(:, k) = ext_basis(:, k) + ext_basis(:, i)
                    end if
                end do
                ! flip sign of i
                ext_basis(:, i) = -ext_basis(:, i)
                ! start over!
                cycle iterloop
            end if
        end do
        end do
        ! If we make it here, we are done!
        exit iterloop
    end do iterloop

    ! Look through all variants to get shortest axis
    dum = 0.0_r8
    dum(:, 1:4) = ext_basis
    dum(:, 5) = ext_basis(:, 1) + ext_basis(:, 2)
    dum(:, 6) = ext_basis(:, 2) + ext_basis(:, 3)
    dum(:, 7) = ext_basis(:, 1) + ext_basis(:, 3)
    do i = 1, 7
        d(i) = norm2(dum(:, i))
    end do
    call lo_qsort(d, ind)

    do i = 3, 7
        canonical_basis(:, 1) = dum(:, ind(1))
        canonical_basis(:, 2) = dum(:, ind(2))
        canonical_basis(:, 3) = dum(:, ind(i))
        if (abs(lo_determ(canonical_basis)) .gt. lo_tol) exit
    end do
    if (lo_determ(canonical_basis) .lt. 0.0_r8) canonical_basis = -canonical_basis
end subroutine

!> Backwards figure out what a rotation matrix does
module subroutine classify_point_operation(m, fold, axis, classification)
    !> operation (Cartesian)
    real(r8), dimension(3, 3), intent(in) :: m
    !> 1,2,3,4,6-fold rotation?
    integer, intent(out) :: fold
    !> axis of rotation, if applicable
    real(r8), dimension(3) :: axis
    !> what kind of operation is it? (1=identiy,2=inversion,3=rotation,4=reflection,5=rotoreflection)
    integer, intent(out) :: classification

    integer :: i, j, k
    real(r8), dimension(3, 3) :: m1, identity
    real(r8), dimension(3) :: v0, v1, welldefinedvalues
    real(r8) :: tr, eta, theta, sintheta, costheta, f0, f1

    ! an identity matrix
    identity = 0.0_r8
    identity(1, 1) = 1.0_r8
    identity(2, 2) = 1.0_r8
    identity(3, 3) = 1.0_r8

    ! Should add some more here.
    welldefinedvalues(1) = 0.0_r8
    welldefinedvalues(2) = 1.0_r8
    welldefinedvalues(3) = 0.5_r8

    ! get some basic things
    eta = lo_determ(m)        ! the determinant
    tr = m(1, 1) + m(2, 2) + m(3, 3) ! the trace
    ! The angle can only be pi,pi/2,pi/3,pi/4,pi/6 or zero
    f0 = (tr - eta)*0.5_r8
    if (abs(f0 + 1.0_r8) .lt. 1E-10_r8) f0 = -1.0_r8
    if (abs(f0 - 1.0_r8) .lt. 1E-10_r8) f0 = 1.0_r8
    f0 = acos(f0)
    if (abs(f0) .lt. lo_radiantol) then
        theta = 0.0_r8
        costheta = cos(theta)
        sintheta = sin(theta)
        fold = 0
    elseif (abs(f0 - lo_twopi/2.0_r8) .lt. lo_radiantol) then
        theta = lo_twopi/2.0_r8
        costheta = cos(theta)
        sintheta = sin(theta)
        fold = 2
    elseif (abs(f0 - lo_twopi/3.0_r8) .lt. lo_radiantol) then
        theta = lo_twopi/3.0_r8
        costheta = cos(theta)
        sintheta = sin(theta)
        fold = 3
    elseif (abs(f0 - lo_twopi/4.0_r8) .lt. lo_radiantol) then
        theta = lo_twopi/4.0_r8
        costheta = cos(theta)
        sintheta = sin(theta)
        fold = 4
    elseif (abs(f0 - lo_twopi/6.0_r8) .lt. lo_radiantol) then
        theta = lo_twopi/6.0_r8
        costheta = cos(theta)
        sintheta = sin(theta)
        fold = 6
    else
        ! small sanity check
        call lo_stop_gracefully(['Strange angle for a rotation matrix: '//tochar(f0*180/lo_pi)], lo_exitcode_symmetry)
    end if

    axis = 0.0_r8
    !whatkind='dunno'
    classification = -1
    ! Ok, got nice angles, try to classify the operations. First the trivial.
    if (norm2(m - identity) .lt. lo_tol) then
        !whatkind='identity'
        theta = 0.0_r8
        fold = -1
        classification = 1
    elseif (norm2(m + identity) .lt. lo_tol) then
        !whatkind='inversion'
        theta = 0.0_r8
        fold = -1
        classification = 2
    elseif (abs(eta*tr + 1.0_r8) .gt. lo_tol .and. abs(eta*tr - 3.0_r8) .gt. lo_tol) then
        ! Easy cases:
        f0 = sqrt((3.0_r8 - eta*tr)*(1.0_r8 + eta*tr)) !0.5_r8
        axis(1) = (m(3, 2) - m(2, 3))/f0
        axis(2) = (m(1, 3) - m(3, 1))/f0
        axis(3) = (m(2, 1) - m(1, 2))/f0

        if (eta .gt. 0.0_r8) then
            ! easy, pure rotation
            !whatkind='rotation'
            classification = 3
        else
            ! pure reflection or rotoreflection
            if (abs(theta) .lt. lo_radiantol) then
                ! I don't think this can happen
                !whatkind='reflection'
                classification = 4
            else
                ! I think this always happens.
                !whatkind='rotoreflection'
                classification = 5
            end if
        end if
    else
        ! Trickier cases
        do i = 1, 3
            f0 = 2*m(i, i) + eta - tr
            v0(i) = f0/(3.0_r8*eta - tr)
            ! check for well defined values. Just to be on the safe side.
            do j = 1, size(welldefinedvalues, 1)
                if (abs(v0(i) - welldefinedvalues(j)) .lt. lo_sqtol) then
                    v0(i) = welldefinedvalues(j)
                end if
            end do
            v0(i) = sqrt(v0(i))
        end do
        ! Figure out the correct signs.
        axis = -2
        f1 = lo_huge
        signloop: do i = -1, 1, 2
        do j = -1, 1, 2
        do k = -1, 1, 2
            v1(1) = i*v0(1)
            v1(2) = j*v0(2)
            v1(3) = k*v0(3)
            if (eta .gt. 0.0_r8) then
                m1 = lo_rotation_matrix_from_axis_and_angle(v1, theta)
            else
                m1 = lo_improper_rotation_matrix_from_axis_and_angle(v1, theta)
            end if
            f0 = norm2(m - m1)
            if (f0 .lt. f1) then
                axis = v1
                f1 = f0
            end if
        end do
        end do
        end do signloop

        if (eta .gt. 0.0_r8) then
            !write(*,*) 'ax',axis
            !whatkind='rotation'
            classification = 3
        else
            if (abs(theta) .lt. lo_radiantol) then
                ! no angle, then it's a pure reflection I suppose
                !whatkind='reflection'
                classification = 4
            else
                !whatkind='rotoreflection'
                classification = 5
            end if
        end if
    end if

    axis = lo_chop(axis, 1E-14_r8)
end subroutine

!> Ensure point operations form a (crystallographic) point group
module subroutine ensure_point_operations_are_a_group(op, tolerance)
    !> set of point operations
    real(r8), dimension(:, :, :), allocatable, intent(inout) :: op
    !> numerical tolerance
    real(r8), intent(in) :: tolerance

    real(r8), dimension(:, :, :), allocatable :: buf0, buf1
    real(r8), dimension(3, 3) :: m0, m1
    integer :: i, j, k, l, n, iter

    ! First keep only the unique operations
    call lo_return_unique(op, buf0, tolerance)
    n = size(buf0, 3)

    ! Decent start. Now make sure it forms a group.
    allocate (buf1(3, 3, 48))
    buf1 = 0.0_r8
    buf1(:, :, 1:n) = buf0

    iterl: do iter = 1, 100
        l = n
        do i = 1, l
        vl2: do j = 1, l
            m0 = matmul(buf1(:, :, i), buf1(:, :, j))
            do k = 1, l
                if (norm2(m0 - buf1(:, :, k)) .lt. tolerance) then
                    cycle vl2
                end if
            end do
            n = n + 1
            if (n .gt. 48) then
                call lo_stop_gracefully(['Point operations are not a crystallographic point group'], lo_exitcode_symmetry, __FILE__, __LINE__)
            end if
            buf1(:, :, n) = m0
            ! it was new!
            cycle iterl
        end do vl2
        end do
        if (n .eq. l) then
            ! If we make it here it's done
            exit iterl
        end if
    end do iterl

    ! Store operations
    if (size(op, 3) .ne. n) then
        deallocate (op)
        allocate (op(3, 3, n))
    end if
    op = buf1(:, :, 1:n)

    ! Put identity first
    m1 = 0.0_r8
    do i = 1, 3
        m1(i, i) = 1.0_r8
    end do

    ! Put identity first
    do i = 1, n
        if (norm2(op(:, :, i) - m1) .lt. tolerance) then
            m0 = op(:, :, 1)
            op(:, :, 1) = m1
            op(:, :, i) = m0
            exit
        end if
    end do
    ! And inversion second, if applicable
    do i = 1, n
        if (norm2(op(:, :, i) + m1) .lt. tolerance) then
            m0 = op(:, :, 2)
            op(:, :, 2) = -m1
            op(:, :, i) = m0
            exit
        end if
    end do
end subroutine

end submodule
