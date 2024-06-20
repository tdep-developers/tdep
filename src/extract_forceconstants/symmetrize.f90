!> Tools to symmetrize things, mostly from
!> https://github.com/ollehellman/aims_wrappers/blob/master/src/relax_with_symmetry_constraints/type_irredcell.f90

module symmetrize

use konstanter, only: r8, lo_sqtol, lo_exitcode_symmetry, lo_status
use type_crystalstructure, only: lo_crystalstructure
use gottochblandat, only: lo_identitymatrix, lo_nullspace_coefficient_matrix, &
                          lo_stop_gracefully, lo_unflatten_2tensor, lo_chop, tochar, &
                          lo_flattentensor
use type_blas_lapack_wrappers, only: lo_dgels
implicit none

contains

!> symmetrize 3x3 stress/strain
subroutine symmetrize_stress(stress, uc, stress_symmetric, verbose)
    !> the matrix that should be symmetrized
    real(r8), dimension(3, 3), intent(in)  :: stress
    real(r8), dimension(3, 3), intent(out) :: stress_symmetric
    logical, intent(in), optional :: verbose
    logical :: talk
    !> the crystal structure with symmetry
    type(lo_crystalstructure), intent(in) :: uc
    real(r8), dimension(:, :), allocatable :: coeffM, wM
    real(r8), dimension(9, 1) :: wV

    real(r8), dimension(:), allocatable :: stress_irred_solution
    real(r8), dimension(9, 9) :: invarM, I9, T, rotM
    real(r8), dimension(3, 3) :: m0
    real(r8), dimension(9) :: dv
    integer :: i, j, n_irred

    talk = .false.
    if (present(verbose)) then
        if (verbose) talk = .true.
    end if

    ! Transposition operator
    T(:, 1) = [1, 0, 0, 0, 0, 0, 0, 0, 0]
    T(:, 2) = [0, 0, 0, 1, 0, 0, 0, 0, 0]
    T(:, 3) = [0, 0, 0, 0, 0, 0, 1, 0, 0]
    T(:, 4) = [0, 1, 0, 0, 0, 0, 0, 0, 0]
    T(:, 5) = [0, 0, 0, 0, 1, 0, 0, 0, 0]
    T(:, 6) = [0, 0, 0, 0, 0, 0, 0, 1, 0]
    T(:, 7) = [0, 0, 1, 0, 0, 0, 0, 0, 0]
    T(:, 8) = [0, 0, 0, 0, 0, 1, 0, 0, 0]
    T(:, 9) = [0, 0, 0, 0, 0, 0, 0, 0, 1]
    ! identity matrix
    call lo_identitymatrix(I9)
    ! now add all the operations of the lattice
    invarM = 0.0_r8
    do i = 1, uc%sym%n
        rotM = expandoperation_pair(uc%sym%op(i)%m)
        invarM = invarM + rotM - I9
    end do
    ! and transposition symmetry
    invarM = invarM + T - I9
    ! and project out the nullspace
    call lo_nullspace_coefficient_matrix(invarM, coeffM, n_irred, tolerance=lo_sqtol)

    ! It's weird if it's not at least 1, or at most 6, since the number of zero singular values
    ! should be the number of degrees of freedom.
    if (n_irred .lt. 1 .or. n_irred .gt. 6) then
        call lo_stop_gracefully(['Strain tensor must have at least 1 and at most 6 degrees of freedom'], &
                                lo_exitcode_symmetry, __FILE__, __LINE__)
    end if
    ! and a little space for the solution
    allocate (stress_irred_solution(n_irred))
    stress_irred_solution = 0.0_r8

    if (talk) then
        write (*, *) '... found ', tochar(n_irred), ' degrees of freedom for the lattice:'
        do i = 1, n_irred
            write (*, *) 'x'//tochar(i)//':'
            dv = 0.0_r8
            dv(i) = 1.0_r8
            m0 = lo_unflatten_2tensor(matmul(coeffM, dv(1:n_irred)))
            do j = 1, 3
                write (*, "(3(2X,F16.9))") m0(:, j)
            end do
        end do
    end if

    ! now solve for the irreducible components
    ! Solve for the irreducible components
    allocate (wM(9, n_irred))
    wM(:, :) = coeffM
    wV(:, 1) = lo_flattentensor(stress)
    call lo_dgels(wM, wV, info=lo_status)
    if (lo_status .ne. 0) then
        write (*, *) 'Lapack failed for the lattice'
        stop
    end if
    stress_irred_solution(:) = wV(1:n_irred, 1)
    ! Store the symetrized stress
    stress_symmetric(1:3, 1:3) = lo_chop(lo_unflatten_2tensor(matmul(coeffM, stress_irred_solution)), 1E-12_r8)
end subroutine

!> symmetrize 3xN forces in unitcell
subroutine symmetrize_forces(forces, uc, forces_sym, verbose)
    !> the matrix that should be symmetrized
    real(r8), dimension(:, :), intent(in)  :: forces
    real(r8), dimension(:, :), intent(out)  :: forces_sym
    logical, intent(in), optional :: verbose
    logical :: talk
    !> the crystal structure with symmetry
    type(lo_crystalstructure), intent(in) :: uc
    real(r8), dimension(:, :), allocatable :: coeffM, dr

    real(r8), dimension(:), allocatable :: forces_irred
    real(r8), dimension(3*uc%na, 3*uc%na) :: invarM, Id, rotM
    real(r8), dimension(:), allocatable :: dv
    integer :: i, j, o, n_irred, dof

    if ((present(verbose)) .and. (verbose)) then
        talk = .true.
    else
        talk = .false.
    end if

    call lo_identitymatrix(Id)

    ! Build a matrix representation how all the internal positions transform under an
    ! operation as large block matrices:
    invarM = 0.0_r8
    do o = 1, uc%sym%n
        rotM = 0.0_r8
        do i = 1, uc%na
            ! the semimagic %fmap here keeps track of which atom is transform to which atom
            ! under each operation. I can mimic this switcharoo by putting the small 3x3
            ! rotation matrix in the right place in the big matrix. That way I don't have to
            ! deal with translations explicitly, and life becomes much easier.
            j = uc%sym%op(o)%fmap(i)
            !rotM((j-1)*3+1:j*3,(i-1)*3+1:3*i)=transpose(p%sym%op(o)%fm)
            ! Mathematica confused me again, this is how it's supposed to be. Hmm.
            rotM((i - 1)*3 + 1:i*3, (j - 1)*3 + 1:3*j) = transpose(uc%sym%op(o)%fm)
        end do
        invarM = invarM + rotM - Id
    end do

    ! Then add the zero-sum constraint? You should convince yourself that this really is the
    ! matrix representation of adding all forces together.
    rotM = 0.0_r8
    do i = 1, uc%na
        do j = 1, 3
            rotM((i - 1)*3 + j, j) = 1.08
        end do
    end do
    invarM = invarM + rotM

    ! rotational invariance? nobody knows

    ! Now we can project out the nullspace of this matrix
    call lo_nullspace_coefficient_matrix(invarM, coeffM, n_irred, tolerance=lo_sqtol)

    if (talk) then
        write (*, *) '... found ', tochar(n_irred), ' internal degrees of freedom'
        dof = uc%na*3
        allocate (dv(dof))
        allocate (dr(3, uc%na))
        do i = 1, n_irred
            write (*, *) 'x'//tochar(i)//': (fractional, Cartesian)'
            dv = 0.0_r8
            dv(i) = 1.0_r8
            dv = matmul(coeffM, dv(1:n_irred))
            do j = 1, uc%na
                dr(:, j) = dv((j - 1)*3 + 1:j*3)
                write (*, "(6(2X,F15.9))") dr(:, j), uc%fractional_to_cartesian(dr(:, j))
            end do
        end do
    end if

    ! now solve for the irreducible components if there is anything to solve
    if (n_irred .gt. 0) then
        fixforces: block
            real(r8), dimension(:, :), allocatable :: mA, mB
            real(r8), dimension(:), allocatable :: vB
            real(r8), dimension(3) :: v0
            integer :: l, n

            n = size(coeffM, 1)
            if (size(forces, 1)*size(forces, 2) .ne. n) then
                write (*, *) 'Wrong number of forces'
                stop
            end if

            ! and some space for a solution
            allocate (forces_irred(n_irred))
            forces_irred = 0.0_r8

            ! normal least squares for forces in fractional coords
            allocate (mA(n, n_irred))
            allocate (mB(n, 1))
            l = 0
            do i = 1, size(forces, 2)
                associate (f_frac => matmul(uc%inv_latticevectors, forces(:, i)))
                    do j = 1, 3
                        l = l + 1
                        mB(l, 1) = f_frac(j)
                    end do
                end associate
            end do
            mA(:, :) = coeffM
            call lo_dgels(mA, mB, info=lo_status)
            if (lo_status .ne. 0) then
                write (*, *) 'Lapack failed for forces'
                stop
            end if
            forces_irred(:) = mB(1:n_irred, 1)

            ! Store the symmetrized forces
            if (.not. allocated(vB)) allocate (vB(n))
            vB(:) = matmul(coeffM, forces_irred)
            do i = 1, uc%na
                v0 = vB((i - 1)*3 + 1:i*3)
                v0 = matmul(uc%latticevectors, v0)
                forces_sym(:, i) = lo_chop(v0, 1E-12_r8)
            end do
        end block fixforces
    else
        ! In case there are no internal degrees of freedom.
        forces_sym = 0.0_r8
    end if
end subroutine

!> expand a 3x3 symmetry operation to a 9x9 operation, to operate on flattened tensors
pure function expandoperation_pair(o) result(bigo)
    !> the original 3x3 operation
    real(r8), dimension(3, 3), intent(in) :: o
    !> the resulting 9x9 operation
    real(r8), dimension(9, 9) :: bigo
    !
    integer :: i, j, ii, jj, k, l
    !
    k = 0
    do i = 1, 3
    do ii = 1, 3
        k = k + 1
        l = 0
        do j = 1, 3
        do jj = 1, 3
            l = l + 1
            bigo(k, l) = o(ii, jj)*o(i, j)
        end do
        end do
    end do
    end do
    bigo = lo_chop(bigo, 1E-11_r8)
end function

end module
