submodule(type_forceconstant_secondorder) type_forceconstant_secondorder_dynamicalmatrix
!!
!! Sort out everything related to dynamical matrices, phonon frequencies and derivatives thereof.
!!
use konstanter, only: lo_groupvel_HartreeBohr_to_ms, lo_exitcode_blaslapack, lo_imag
use gottochblandat, only: lo_negsqrt, lo_complex_gram_schmidt
use lo_sorting, only: lo_qsort
use type_blas_lapack_wrappers, only: lo_zheev, lo_zgeev
implicit none
contains

!> get frequencies, eigenvectors and group velocities
module subroutine frequencies_eigenvectors_groupvelocities(fc, dynamical_matrix, omega, mem, dynamical_matrix_gradient, dynamical_matrix_hessian, eigenvectors, groupvelocities, grouphessian, qpoint)
    !> forceconstant
    class(lo_forceconstant_secondorder), intent(in) :: fc
    !> dynamical matrix
    complex(r8), dimension(:, :), intent(in) :: dynamical_matrix
    !> frequencies
    real(r8), dimension(:), intent(out) :: omega
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> gradient of dynamical matrix
    complex(r8), dimension(:, :, :), intent(in), optional :: dynamical_matrix_gradient
    !> gradient of dynamical matrix
    complex(r8), dimension(:, :, :), intent(in), optional :: dynamical_matrix_hessian
    !> eigenvectors
    complex(r8), dimension(:, :), intent(out), optional :: eigenvectors
    !> group velocities
    real(r8), dimension(:, :), intent(out), optional :: groupvelocities
    !> group hessian
    real(r8), dimension(:, :, :), intent(out), optional :: grouphessian
    !> the actual q-point
    class(lo_qpoint), intent(in), optional :: qpoint

    logical :: calceig, calcvel, calchess, hermitian
    complex(r8), dimension(:, :, :), allocatable :: wHessian
    complex(r8), dimension(:, :), allocatable :: wEigenvector, wVelocities
    real(r8), dimension(:), allocatable :: wEigenval
    integer, dimension(:, :), allocatable :: subspace_mode
    integer, dimension(:), allocatable :: subspace_ctr
    integer :: nb, n_subspace

    ! Start with some heuristics
    init: block
        complex(r8) :: c0
        integer :: i, j

        ! a little shorthand
        nb = fc%na*3
        ! return the eigenvectors?
        if (present(eigenvectors)) then
            calceig = .true.
        else
            calceig = .false.
        end if
        ! return group velocities?
        if (present(groupvelocities)) then
            calcvel = .true.
            if (.not. present(dynamical_matrix_gradient)) call lo_stop_gracefully(['Need the gradient to calculate group velocities'], lo_exitcode_param, __FILE__, __LINE__)
        else
            calcvel = .false.
        end if
        ! return the Hessian w.r.t k?
        if (present(grouphessian)) then
            calchess = .true.
            calcvel = .true.
            if (.not. present(dynamical_matrix_gradient)) call lo_stop_gracefully(['Need the gradient to calculate group Hessian'], lo_exitcode_param, __FILE__, __LINE__)
            if (.not. present(dynamical_matrix_hessian)) call lo_stop_gracefully(['Need the hessian to calculate group Hessian'], lo_exitcode_param, __FILE__, __LINE__)
        else
            calchess = .false.
        end if

        ! quick test to see if the dynamical matrix is Hermitian
        ! should probably report a warning somewhere if this fails.
        ! there are in fact sensible reasons that are not completely
        ! stupid that cause it to become slightly non-hermitian at times.
        ! I don't have the energy to explain now.
        c0 = 0.0_r8
        do i = 1, nb
        do j = i + 1, nb
            c0 = c0 + dynamical_matrix(i, j) - conjg(dynamical_matrix(j, i))
        end do
        end do
        if (abs(c0) .lt. lo_sqtol*nb) then
            hermitian = .true.
        else
            hermitian = .false.
        end if

        ! Make a bit of workspace.
        call mem%allocate(wEigenval, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(wEigenvector, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(wVelocities, [3, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(wHessian, [3, 3, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        wEigenval = 0.0_r8
        wEigenvector = 0.0_r8
        wVelocities = 0.0_r8
        wHessian = 0.0_r8
    end block init

    ! Diagonalize the Hamiltonian
    diagonalize: block
        complex(r8), dimension(:, :), allocatable :: wDynmat
        complex(r8), dimension(:), allocatable :: wV
        real(r8), dimension(:), allocatable :: dr
        integer, dimension(:), allocatable :: ind
        integer :: i, j

        ! Some temporary workspace
        call mem%allocate(wDynmat, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(dr, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ind, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        wDynmat = 0.0_r8
        dr = 0.0_r8
        ind = 0

        ! do the actual diagonalization.
        if (hermitian) then
            ! Hermitian routine
            wEigenvector = dynamical_matrix
            call lo_zheev(wEigenvector, wEigenval, jobz='V', info=i)
            if (i .ne. 0) call lo_stop_gracefully(['zheev exit status '//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)
            ! Orthogonalize the eigenvectors
            !call lo_complex_gram_schmidt(wEigenvector)
        else
            ! Non-hermitian routine
            call mem%allocate(wV, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            wDynmat = dynamical_matrix
            wV = 0.0_r8
            call lo_zgeev(wDynmat, wV, wEigenvector, info=i)
            if (i .ne. 0) call lo_stop_gracefully(['zgeev exit status '//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)
            ! wV contains the (complex) eigenvalues to the dynamical matrix. If they have any
            ! imaginary components, I set them negative, since that is obviously an unphysical result.
            do i = 1, nb
                if (abs(aimag(wV(i))) .gt. lo_sqtol) then
                    wEigenval(i) = -abs(wV(i))
                else
                    wEigenval(i) = real(wV(i))
                end if
            end do
            call mem%deallocate(wV, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        ! Set the three smallest eigenvalues to zero, at gamma
        if (present(qpoint)) then
        if (lo_sqnorm(qpoint%r) .lt. lo_sqtol) then
            dr = abs(wEigenval)     ! dummy copy of eigenvalues
            call lo_qsort(dr, ind) ! get a sorted list of eigenvalues
            do i = 1, 3
                wEigenval(ind(i)) = 0.0_r8
            end do
        end if
        end if

        ! Then sort eigenvalues according to size, and eigenvectors accordingly.
        wDynmat = wEigenvector
        dr = lo_negsqrt(wEigenval)
        call lo_qsort(dr, ind)
        dr = wEigenval
        do i = 1, nb
            j = ind(i)
            wEigenval(i) = dr(j)
            wEigenvector(:, i) = wDynmat(:, j)
        end do

        call mem%deallocate(wDynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dr, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ind, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block diagonalize

    ! Calculate group velocities.
    if (calcvel) then
        ! figure out what we are doing with degeneracies (or near degeneracies):
        call return_degenerate_subspaces(wEigenval, n_subspace, subspace_ctr, subspace_mode, mem)
        if (n_subspace .lt. 0) then
            ! This is the easy case where we can just apply non-degenerate
            ! perturbation theory directly. Should be the most common one
            ! as well.
            groupvel: block
                complex(r8), dimension(:, :), allocatable :: calpha, cbeta
                complex(r8), dimension(:), allocatable :: wV
                real(r8) :: f0
                integer :: i, j, ia, ib, ih, imode, jmode

                ! Set the velocities to nothing
                wVelocities = 0.0_r8
                wHessian = 0.0_r8

                ! A little workspace
                call mem%allocate(wV, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%allocate(calpha, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%allocate(cbeta, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                wV = 0.0_r8
                calpha = 0.0_r8
                cbeta = 0.0_r8

                ! First the group velocities, just straight-up Hellman-Feynmann:
                if (calcvel) then
                    ! dw^2/dk = < e | dD/dq | e >
                    do imode = 1, nb
                    do ia = 1, 3
                        call zgemv('N', nb, nb, (1.0_r8, 0.0_r8), dynamical_matrix_gradient(:, :, ia), nb, wEigenvector(:, imode), 1, (0.0_r8, 0.0_r8), wV, 1)
                        wVelocities(ia, imode) = dot_product(wEigenvector(:, imode), wV)
                    end do
                    end do
                    ! This was dw^2/dq, convert to dw/dq
                    do imode = 1, nb
                        f0 = lo_negsqrt(wEigenval(imode))
                        if (f0 .gt. lo_freqtol) then
                            wVelocities(:, imode) = wVelocities(:, imode)/(2*f0)
                        else
                            wVelocities(:, imode) = 0.0_r8
                        end if
                    end do
                end if

                ! Then the Hessian.
                if (calchess) then
                    wHessian = 0.0_r8
                    prtloop: do ih = 1, 6
                        ! Which components are we dealing with for this perturbation
                        ! the missing ones are filled in with symmetry later.
                        select case (ih)
                        case (1); ia = 1; ib = 1 ! xx
                        case (2); ia = 1; ib = 2 ! xy
                        case (3); ia = 1; ib = 3 ! xz
                        case (4); ia = 2; ib = 2 ! yy
                        case (5); ia = 2; ib = 3 ! yz
                        case (6); ia = 3; ib = 3 ! zz
                        end select

                        ! So, now ih points to the correct hessian, and
                        ! ia and ib to the correct first order perturbation.
                        ! build the matrices that correct the wavefunctions:
                        calpha = 0.0_r8
                        cbeta = 0.0_r8
                        do i = 1, nb
                        do j = 1, nb
                            if (i .eq. j) cycle
                            ! < e_i | h^a | e_j >
                            call zgemv('N', nb, nb, (1.0_r8, 0.0_r8), dynamical_matrix_gradient(:, :, ia), nb, wEigenvector(:, j), 1, (0.0_r8, 0.0_r8), wV, 1)
                            calpha(i, j) = dot_product(wEigenvector(:, i), wV)
                            ! < e_i | h^b | e_j >
                            call zgemv('N', nb, nb, (1.0_r8, 0.0_r8), dynamical_matrix_gradient(:, :, ib), nb, wEigenvector(:, j), 1, (0.0_r8, 0.0_r8), wV, 1)
                            cbeta(i, j) = dot_product(wEigenvector(:, i), wV)
                        end do
                        end do
                        ! Then build the relevant matrix elements
                        do imode = 1, nb
                            do jmode = 1, nb
                                if (imode .eq. jmode) cycle
                                whessian(ia, ib, imode) = whessian(ia, ib, imode) + calpha(imode, jmode)*cbeta(jmode, imode)/(wEigenval(imode) - wEigenval(jmode))
                                whessian(ia, ib, imode) = whessian(ia, ib, imode) + cbeta(imode, jmode)*calpha(jmode, imode)/(wEigenval(imode) - wEigenval(jmode))
                            end do
                            ! 2 < e | h^ab | e >
                            call zgemv('N', nb, nb, (1.0_r8, 0.0_r8), dynamical_matrix_hessian(:, :, ih), nb, wEigenvector(:, imode), 1, (0.0_r8, 0.0_r8), wV, 1)
                            whessian(ia, ib, imode) = whessian(ia, ib, imode) + dot_product(wEigenvector(:, imode), wV)
                        end do
                    end do prtloop

                    ! This was d^2 w^2/dqi dqj, convert to d^2 w/dqi qdj
                    ! note that I need the group velocities here.
                    do imode = 1, nb
                        f0 = lo_negsqrt(wEigenval(imode))
                        if (f0 .gt. lo_freqtol) then
                            do ia = 1, 3
                            do ib = ia, 3
                                wHessian(ia, ib, imode) = wHessian(ia, ib, imode) - 2*wVelocities(ia, imode)*wVelocities(ib, imode)
                                if (ia .ne. ib) wHessian(ib, ia, imode) = wHessian(ia, ib, imode)
                            end do
                            end do
                            wHessian(:, :, imode) = wHessian(:, :, imode)/(2*f0)
                        else
                            wHessian(:, :, imode) = 0.0_r8
                        end if
                    end do
                end if ! calchess

                ! Cleanup
                call mem%deallocate(wV, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%deallocate(calpha, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%deallocate(cbeta, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            end block groupvel
        else
            ! We have to deal with degeneracies. This is not fun. Should be rare, hopefully.
            ! Now with the correct subspace, proceed like normal.
            degengroupvel: block
                complex(r8), dimension(:, :), allocatable :: wPhiEigenvector, wPhiHamiltonian, calpha, cbeta
                complex(r8), dimension(:), allocatable :: wV, wU
                real(r8), dimension(:), allocatable :: wPhiEigenvalue
                real(r8) :: f0
                integer :: i, j, ia, ib, ih, is, imode, jmode, mb

                ! Set the velocities to nothing
                wVelocities = 0.0_r8
                wHessian = 0.0_r8

                ! A little extra space:
                call mem%allocate(wV, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%allocate(wU, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%allocate(calpha, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%allocate(cbeta, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                wV = 0.0_r8
                wU = 0.0_r8
                calpha = 0.0_r8
                cbeta = 0.0_r8

                ! Start with just the group velocities, I think. Loop over the degenerate subspaces:
                spcloop: do is = 1, n_subspace
                    ! There might be a quick exit:
                    if (subspace_ctr(is) .eq. 1) then
                        ! This is a trivial subspace, only one mode. Same as above.
                        imode = subspace_mode(1, is)
                        do ia = 1, 3
                            call zgemv('N', nb, nb, (1.0_r8, 0.0_r8), dynamical_matrix_gradient(:, :, ia), nb, wEigenvector(:, imode), 1, (0.0_r8, 0.0_r8), wV, 1)
                            wVelocities(ia, imode) = dot_product(wEigenvector(:, imode), wV)
                        end do
                        ! and move on to the next!
                        cycle spcloop
                    end if

                    ! No quick exit, have to deal with the subspace the hard way. No worries.
                    ! Shorthand for the dimension of the degenerate subspace:
                    mb = subspace_ctr(is)
                    ! Space for temporary Hamiltonian guy
                    call mem%allocate(wPhiEigenvector, [nb, mb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                    call mem%allocate(wPhiHamiltonian, [mb, mb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                    call mem%allocate(wPhiEigenvalue, mb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                    wPhiEigenvector = 0.0_r8
                    wPhiHamiltonian = 0.0_r8
                    wPhiEigenvalue = 0.0_r8
                    ! Fetch the eigenvectors for this subspace
                    do i = 1, mb
                        imode = subspace_mode(i, is)
                        wPhiEigenvector(:, i) = wEigenvector(:, imode)
                    end do

                    ! Loop over pertubations
                    prtloop: do ia = 1, 3
                        ! Build the Hamiltonian matrix elements for the perturbing Hamiltonian in
                        ! the degenerate subspace, normal degenerate perturbation theory.
                        wPhiHamiltonian = 0.0_r8
                        do i = 1, mb
                        do j = i, mb
                            call zgemv('N', nb, nb, (1.0_r8, 0.0_r8), dynamical_matrix_gradient(:, :, ia), nb, wPhiEigenvector(:, j), 1, (0.0_r8, 0.0_r8), wV, 1)
                            wPhiHamiltonian(i, j) = dot_product(wPhiEigenvector(:, i), wV)
                            if (i .ne. j) wPhiHamiltonian(j, i) = conjg(wPhiHamiltonian(j, i))
                        end do
                        end do
                        ! Diagonalize this?
                        call lo_zheev(wPhiHamiltonian, wPhiEigenvalue, jobz='V', info=i)
                        if (i .ne. 0) call lo_stop_gracefully(['zheev exit status '//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)

                        ! For each of the degenerate modes, store the average of the eigenvalues:
                        do i = 1, mb
                            imode = subspace_mode(i, is)
                            wVelocities(ia, imode) = sum(wPhiEigenvalue)/real(mb, r8)
                        end do
                    end do prtloop
                    ! And cleanup
                    call mem%deallocate(wPhiEigenvector, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                    call mem%deallocate(wPhiHamiltonian, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                    call mem%deallocate(wPhiEigenvalue, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                end do spcloop

                ! This was dw^2/dq, convert to dw/dq
                do imode = 1, nb
                    f0 = lo_negsqrt(wEigenval(imode))
                    if (f0 .gt. lo_freqtol) then
                        wVelocities(:, imode) = wVelocities(:, imode)/(2*f0)
                    else
                        wVelocities(:, imode) = 0.0_r8
                    end if
                end do

                ! Now start dealing with the Hessian.
                if (calchess) then
                    prtloop2: do ih = 1, 6
                        ! Which components are we dealing with for this perturbation
                        ! the missing ones are filled in with symmetry later.
                        select case (ih)
                        case (1); ia = 1; ib = 1 ! xx
                        case (2); ia = 1; ib = 2 ! xy
                        case (3); ia = 1; ib = 3 ! xz
                        case (4); ia = 2; ib = 2 ! yy
                        case (5); ia = 2; ib = 3 ! yz
                        case (6); ia = 3; ib = 3 ! zz
                        end select
                        ! So, now ih points to the correct hessian, and
                        ! ia and ib to the correct first order perturbation.
                        ! build the matrices that correct the wavefunctions:
                        calpha = 0.0_r8
                        cbeta = 0.0_r8
                        do i = 1, nb
                        do j = i, nb
                            ! < e_i | h^a | e_j >
                            call zgemv('N', nb, nb, (1.0_r8, 0.0_r8), dynamical_matrix_gradient(:, :, ia), nb, wEigenvector(:, j), 1, (0.0_r8, 0.0_r8), wV, 1)
                            calpha(i, j) = dot_product(wEigenvector(:, i), wV)
                            ! < e_i | h^b | e_j >
                            call zgemv('N', nb, nb, (1.0_r8, 0.0_r8), dynamical_matrix_gradient(:, :, ib), nb, wEigenvector(:, j), 1, (0.0_r8, 0.0_r8), wV, 1)
                            cbeta(i, j) = dot_product(wEigenvector(:, i), wV)
                            if (i .ne. j) then
                                calpha(j, i) = conjg(calpha(i, j))
                                cbeta(j, i) = conjg(cbeta(i, j))
                            end if
                        end do
                        end do

                        ! Now try to resolve the subspace thingies:
                        do is = 1, n_subspace
                            if (subspace_ctr(is) .eq. 1) then
                                imode = subspace_mode(1, is)
                                ! Add the two first-order corrections
                                do jmode = 1, nb
                                    if (jmode .eq. imode) cycle
                                    whessian(ia, ib, imode) = whessian(ia, ib, imode) + 2*calpha(imode, jmode)*cbeta(jmode, imode)/(wEigenval(imode) - wEigenval(jmode))
                                end do
                                ! Add the normal second-order guy.
                                call zgemv('N', nb, nb, (1.0_r8, 0.0_r8), dynamical_matrix_hessian(:, :, ih), nb, wEigenvector(:, imode), 1, (0.0_r8, 0.0_r8), wV, 1)
                                whessian(ia, ib, imode) = whessian(ia, ib, imode) + dot_product(wEigenvector(:, imode), wV)
                            else
                                ! Do many many things. Hmm.
                            end if
                        end do

                    end do prtloop2

                    ! This was d^2 w^2/dqi dqj, convert to d^2 w/dqi qdj
                    do imode = 1, nb
                        f0 = lo_negsqrt(wEigenval(imode))
                        if (f0 .gt. lo_freqtol) then
                            do ia = 1, 3
                            do ib = ia, 3
                                wHessian(ia, ib, imode) = wHessian(ia, ib, imode) - 2*wVelocities(ia, imode)*wVelocities(ib, imode)
                                if (ia .ne. ib) wHessian(ib, ia, imode) = wHessian(ia, ib, imode)
                            end do
                            end do
                            wHessian(:, :, imode) = wHessian(:, :, imode)/(2*f0)
                        else
                            wHessian(:, :, imode) = 0.0_r8
                        end if
                    end do

                    ! Kill the Hessian at degenerate points. Can not get it to work.
                    !do is=1,n_subspace
                    !    if ( subspace_ctr(is) .eq. 1 ) cycle
                    !    do i=1,subspace_ctr(is)
                    !        imode=subspace_mode(i,is)
                    !        wHessian(:,:,imode)=0.0_r8
                    !    enddo
                    !enddo
                end if

                ! And clean things we are no longer using
                call mem%deallocate(wV, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%deallocate(wU, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%deallocate(calpha, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%deallocate(cbeta, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            end block degengroupvel
        end if
        ! And some cleanup, perhaps?
        if (allocated(subspace_ctr)) call mem%deallocate(subspace_ctr, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        if (allocated(subspace_mode)) call mem%deallocate(subspace_mode, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
    end if ! calcvel

    ! Decide what to return, and cleanup things
    returnvals: block
        ! Sort the frequencies
        ! Return whatever
        omega = lo_negsqrt(wEigenval)
        if (present(eigenvectors)) eigenvectors = wEigenvector
        if (present(groupvelocities)) groupvelocities = lo_chop(real(wVelocities), 1E-9/lo_groupvel_Hartreebohr_to_ms)
        if (present(grouphessian)) grouphessian = real(wHessian)

        ! And some cleanup
        call mem%deallocate(wEigenval, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(wEigenvector, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(wVelocities, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(wHessian, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block returnvals
end subroutine

!> calculate the dynamical matrix
module subroutine dynamicalmatrix(fc, p, qpoint, dynamical_matrix, mem, dynamical_matrix_gradient, dynamical_matrix_hessian, qdirection, skipnonanalytical)
    !> forceconstant
    class(lo_forceconstant_secondorder), intent(inout) :: fc
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-point
    class(lo_qpoint), intent(in) :: qpoint
    !> dynamical matrix
    complex(r8), dimension(:, :), intent(out) :: dynamical_matrix
    !> memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> gradient of dynamical matrix (:,:,x,y,z)
    complex(r8), dimension(:, :, :), intent(out), optional :: dynamical_matrix_gradient
    !> hessian of dynamical matrix (:,:,xx,xy,xz,yy,yz,zz)
    complex(r8), dimension(:, :, :), intent(out), optional :: dynamical_matrix_hessian
    !> direction approaching gamma
    real(r8), dimension(3), intent(in), optional :: qdirection
    !> skip the non-analytical term
    logical, intent(in), optional :: skipnonanalytical

    complex(r8), dimension(:, :, :, :, :), allocatable :: gradDt, hessDt, gradPt, hessPt
    complex(r8), dimension(:, :, :, :), allocatable :: Dt, Pt
    real(r8), dimension(:, :, :, :), allocatable :: Gt
    real(r8), dimension(3) :: qv, qdir
    integer :: solmode, nb, na
    logical :: nonanalytical

    ! Figure out what to do
    init: block
        integer :: i
        ! Some shorthand
        na = fc%na
        nb = fc%na*3
        qv = lo_chop(qpoint%r*lo_twopi, 1E-13_r8)

        ! Figure out which gradients, if any, I want
        i = 0
        if (present(dynamical_matrix_gradient)) i = i + 1
        if (present(dynamical_matrix_hessian)) i = i + 2
        ! Now we have the options 0,1,2,3. 2 Means only the hessian, and not
        ! the gradient. That should probably never happen.
        select case (i)
        case (0)
            ! Just the dynamical matrix
            solmode = 1
        case (1)
            ! Dynamical matrix and gradient
            solmode = 2
        case (3)
            ! Dynamical matrix, gradient and hessian
            solmode = 3
        case default
            call lo_stop_gracefully(['Can not solve for Hessian without gradient present.'], lo_exitcode_param, __FILE__, __LINE__)
        end select

        ! Sanity test:
        select case (fc%loto%correctiontype)
        case (1:2)
            call lo_stop_gracefully(['This polar correction type is deprecated.'], lo_exitcode_param, __FILE__, __LINE__)
        end select

        ! Are we at Gamma and should add the non-analytical term?
        if (lo_sqnorm(qpoint%r) .lt. lo_sqtol) then
            if (present(skipnonanalytical)) then
                if (skipnonanalytical) then
                    nonanalytical = .false.
                else
                    nonanalytical = .true.
                end if
            else
                nonanalytical = .true.
            end if
        else
            nonanalytical = .false.
        end if

        ! Do I have a direction?
        if (present(qdirection)) then
            qdir = qdirection
        else
            ! Just pick one, not sure what is the correct thing to do.
            qdir = [1.0_r8, 0.0_r8, 0.0_r8]
        end if

        ! Start allocating tons of space. We always need the normal dynamical matrix.
        call mem%allocate(Dt, [3, 3, na, na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        Dt = 0.0_r8
        ! And we might need the polar dynmat
        select case (fc%loto%correctiontype)
        case (2:3)
            call mem%allocate(Pt, [3, 3, na, na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            Pt = 0.0_r8
        end select
        ! And we might need the non-analytical term.
        if (nonanalytical) then
            call mem%allocate(Gt, [3, 3, na, na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            Gt = 0.0_r8
        end if

        ! Now, for the gradients and such it gets worse.
        select case (solmode)
        case (1)
            ! Nothing extra.
        case (2)
            ! Just gradient
            call mem%allocate(gradDt, [3, 3, na, na, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            gradDt = 0.0_r8
            select case (fc%loto%correctiontype)
            case (0)
                ! Nothing
            case (3)
                call mem%allocate(gradPt, [3, 3, na, na, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                gradPt = 0.0_r8
            end select
        case (3)
            ! Gradient and Hessian
            call mem%allocate(gradDt, [3, 3, na, na, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(hessDt, [3, 3, na, na, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            gradDt = 0.0_r8
            hessDt = 0.0_r8
            select case (fc%loto%correctiontype)
            case (0)
                ! Nothing
            case (3)
                call mem%allocate(gradPt, [3, 3, na, na, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%allocate(hessPt, [3, 3, na, na, 6], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                gradPt = 0.0_r8
                hessPt = 0.0_r8
            end select
        end select
    end block init

    ! Now calculate the actual dynamical matrix.
    getdynmat: block
        complex(r8) :: exp_iqr
        real(r8), dimension(3) :: Rv
        real(r8) :: qdotr
        integer :: i, a1, a2
        ! Calculate the Dynamical matrices
        select case (fc%loto%correctiontype)
        case (0)
            ! No electrostatic corrections at all.
            do a1 = 1, fc%na
            do i = 1, fc%atom(a1)%n
                a2 = fc%atom(a1)%pair(i)%i2
                Rv = fc%atom(a1)%pair(i)%lv2
                qdotr = dot_product(Rv, qv)
                exp_iqr = cmplx(cos(qdotr), sin(qdotr), r8)
                select case (solmode)
                case (1)
                    Dt(:, :, a1, a2) = Dt(:, :, a1, a2) + exp_iqr*fc%atom(a1)%pair(i)%m
                case (2)
                    Dt(:, :, a1, a2) = Dt(:, :, a1, a2) + exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 1) = gradDt(:, :, a1, a2, 1) + lo_imag*Rv(1)*exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 2) = gradDt(:, :, a1, a2, 2) + lo_imag*Rv(2)*exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 3) = gradDt(:, :, a1, a2, 3) + lo_imag*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                case (3)
                    Dt(:, :, a1, a2) = Dt(:, :, a1, a2) + exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 1) = gradDt(:, :, a1, a2, 1) + lo_imag*Rv(1)*exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 2) = gradDt(:, :, a1, a2, 2) + lo_imag*Rv(2)*exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 3) = gradDt(:, :, a1, a2, 3) + lo_imag*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 1) = hessDt(:, :, a1, a2, 1) - Rv(1)*Rv(1)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 2) = hessDt(:, :, a1, a2, 2) - Rv(1)*Rv(2)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 3) = hessDt(:, :, a1, a2, 3) - Rv(1)*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 4) = hessDt(:, :, a1, a2, 4) - Rv(2)*Rv(2)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 5) = hessDt(:, :, a1, a2, 5) - Rv(2)*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 6) = hessDt(:, :, a1, a2, 6) - Rv(3)*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                end select
            end do
            end do
        case (3)
            ! Incredibly smart correction. First the shortranged-part:
            do a1 = 1, fc%na
            do i = 1, fc%atom(a1)%n
                a2 = fc%atom(a1)%pair(i)%i2
                Rv = fc%atom(a1)%pair(i)%lv2
                qdotr = dot_product(Rv, qv)
                exp_iqr = cmplx(cos(qdotr), sin(qdotr), r8)
                select case (solmode)
                case (1)
                    Dt(:, :, a1, a2) = Dt(:, :, a1, a2) + exp_iqr*fc%atom(a1)%pair(i)%m
                case (2)
                    Dt(:, :, a1, a2) = Dt(:, :, a1, a2) + exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 1) = gradDt(:, :, a1, a2, 1) + lo_imag*Rv(1)*exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 2) = gradDt(:, :, a1, a2, 2) + lo_imag*Rv(2)*exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 3) = gradDt(:, :, a1, a2, 3) + lo_imag*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                case (3)
                    Dt(:, :, a1, a2) = Dt(:, :, a1, a2) + exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 1) = gradDt(:, :, a1, a2, 1) + lo_imag*Rv(1)*exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 2) = gradDt(:, :, a1, a2, 2) + lo_imag*Rv(2)*exp_iqr*fc%atom(a1)%pair(i)%m
                    gradDt(:, :, a1, a2, 3) = gradDt(:, :, a1, a2, 3) + lo_imag*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 1) = hessDt(:, :, a1, a2, 1) - Rv(1)*Rv(1)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 2) = hessDt(:, :, a1, a2, 2) - Rv(1)*Rv(2)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 3) = hessDt(:, :, a1, a2, 3) - Rv(1)*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 4) = hessDt(:, :, a1, a2, 4) - Rv(2)*Rv(2)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 5) = hessDt(:, :, a1, a2, 5) - Rv(2)*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                    hessDt(:, :, a1, a2, 6) = hessDt(:, :, a1, a2, 6) - Rv(3)*Rv(3)*exp_iqr*fc%atom(a1)%pair(i)%m
                end select
            end do
            end do
            ! Then the polar part
            select case (solmode)
            case (1)
                call longrange_dynamical_matrix(fc, Pt, p, qpoint%r)
            case (2)
                call longrange_dynamical_matrix(fc, Pt, p, qpoint%r, gradPt(:, :, :, :, 1), gradPt(:, :, :, :, 2), gradPt(:, :, :, :, 3))
            case (3)
                call lo_stop_gracefully(['Have not implemented Hessian for polar corrections yet.'], lo_exitcode_param, __FILE__, __LINE__)
            end select

            ! And the final non-analytical part
            if (nonanalytical) then
                call nonanalytical_dynamical_matrix(fc, p, qdir, Gt)
            end if

            ! Add it together
            select case (solmode)
            case (1)
                Dt = Dt + Pt
                if (nonanalytical) Dt = Dt + Gt
            case (2)
                Dt = Dt + Pt
                if (nonanalytical) Dt = Dt + Gt
                gradDt = gradDt + gradPt
            case (3)
                Dt = Dt + Pt
                if (nonanalytical) Dt = Dt + Gt
                gradDt = gradDt + gradPt
                hessDt = hessDt + hessPt
            end select
        end select
    end block getdynmat

    ! And fix some final things and cleanup.
    finalfix: block
        integer :: a1, a2, i, j, ii, jj
        ! Multiply in masses
        do a1 = 1, na
        do a2 = 1, na
            select case (solmode)
            case (1)
                Dt(:, :, a1, a2) = Dt(:, :, a1, a2)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
            case (2)
                Dt(:, :, a1, a2) = Dt(:, :, a1, a2)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                gradDt(:, :, a1, a2, 1) = gradDt(:, :, a1, a2, 1)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                gradDt(:, :, a1, a2, 2) = gradDt(:, :, a1, a2, 2)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                gradDt(:, :, a1, a2, 3) = gradDt(:, :, a1, a2, 3)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
            case (3)
                Dt(:, :, a1, a2) = Dt(:, :, a1, a2)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                gradDt(:, :, a1, a2, 1) = gradDt(:, :, a1, a2, 1)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                gradDt(:, :, a1, a2, 2) = gradDt(:, :, a1, a2, 2)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                gradDt(:, :, a1, a2, 3) = gradDt(:, :, a1, a2, 3)*p%invsqrtmass(a1)*p%invsqrtmass(a2)

                hessDt(:, :, a1, a2, 1) = hessDt(:, :, a1, a2, 1)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                hessDt(:, :, a1, a2, 2) = hessDt(:, :, a1, a2, 2)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                hessDt(:, :, a1, a2, 3) = hessDt(:, :, a1, a2, 3)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                hessDt(:, :, a1, a2, 4) = hessDt(:, :, a1, a2, 4)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                hessDt(:, :, a1, a2, 5) = hessDt(:, :, a1, a2, 5)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
                hessDt(:, :, a1, a2, 6) = hessDt(:, :, a1, a2, 6)*p%invsqrtmass(a1)*p%invsqrtmass(a2)
            end select
        end do
        end do
        ! Flatten it
        do a1 = 1, na
        do a2 = 1, na
            do i = 1, 3
            do j = 1, 3
                ii = (a2 - 1)*3 + j
                jj = (a1 - 1)*3 + i
                select case (solmode)
                case (1)
                    dynamical_matrix(ii, jj) = Dt(i, j, a1, a2)
                case (2)
                    dynamical_matrix(ii, jj) = Dt(i, j, a1, a2)
                    dynamical_matrix_gradient(ii, jj, 1) = gradDt(i, j, a1, a2, 1)
                    dynamical_matrix_gradient(ii, jj, 2) = gradDt(i, j, a1, a2, 2)
                    dynamical_matrix_gradient(ii, jj, 3) = gradDt(i, j, a1, a2, 3)
                case (3)
                    dynamical_matrix(ii, jj) = Dt(i, j, a1, a2)
                    dynamical_matrix_gradient(ii, jj, 1) = gradDt(i, j, a1, a2, 1)
                    dynamical_matrix_gradient(ii, jj, 2) = gradDt(i, j, a1, a2, 2)
                    dynamical_matrix_gradient(ii, jj, 3) = gradDt(i, j, a1, a2, 3)

                    dynamical_matrix_hessian(ii, jj, 1) = hessDt(i, j, a1, a2, 1)
                    dynamical_matrix_hessian(ii, jj, 2) = hessDt(i, j, a1, a2, 2)
                    dynamical_matrix_hessian(ii, jj, 3) = hessDt(i, j, a1, a2, 3)
                    dynamical_matrix_hessian(ii, jj, 4) = hessDt(i, j, a1, a2, 4)
                    dynamical_matrix_hessian(ii, jj, 5) = hessDt(i, j, a1, a2, 5)
                    dynamical_matrix_hessian(ii, jj, 6) = hessDt(i, j, a1, a2, 6)
                end select
            end do
            end do
        end do
        end do

        ! And a little cleanup
        if (allocated(gradDt)) call mem%deallocate(gradDt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (allocated(hessDt)) call mem%deallocate(hessDt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (allocated(gradPt)) call mem%deallocate(gradPt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (allocated(hessPt)) call mem%deallocate(hessPt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (allocated(Dt)) call mem%deallocate(Dt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (allocated(Gt)) call mem%deallocate(Gt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (allocated(Pt)) call mem%deallocate(Pt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block finalfix
end subroutine

! !> Non-analytical contribution at Gamma
! module subroutine dynmat_at_gamma(fc,uc,qdir,dynmat)
!     !> forceconstant
!     type(lo_forceconstant_secondorder), intent(in) :: fc
!     !> structure
!     type(lo_crystalstructure), intent(in) :: uc
!     !> direction of approach to Gamma
!     real(r8), dimension(3), intent(in) :: qdir
!     !> non-analytical dynamical matrix
!     real(r8), dimension(:,:,:,:), intent(out) :: dynmat
!
!     real(r8), dimension(3) :: q
!     real(r8) :: f0
!     integer :: a1,a2,i,j
!     q=qdir/norm2(qdir)
!
!     f0=0.0_r8
!     do j=1,3
!     do i=1,3
!         f0=f0+q(i)*fc%loto%eps(i,j)*q(j)
!     enddo
!     enddo
!
!     f0=(1.0_r8/f0)*4.0_r8*lo_pi/uc%volume
!     dynmat=0.0_r8
!     do a1=1,uc%na
!     do a2=1,uc%na
!         do i=1,3
!         do j=1,3
!             dynmat(i,j,a1,a2)=dot_product(q,fc%loto%born_effective_charges(:,i,a1))*dot_product(q,fc%loto%born_effective_charges(:,j,a2))*f0
!         enddo
!         enddo
!     enddo
!     enddo
! end subroutine

! Below are not exposed

! !> The dumbest possible correction for long-range interactions
! subroutine add_dynmat_to_forceconstant(fc,p,qpoint,mem)
!     !> forceconstant
!     type(lo_forceconstant_secondorder), intent(inout) :: fc
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!     !> q-points
!     class(lo_qpoint), intent(in) :: qpoint
!     !> memory tracker
!     type(lo_mem_helper), intent(inout) :: mem
!
!     real(r8), dimension(:,:,:,:), allocatable :: D
!     real(r8), dimension(3) :: qv
!     real(r8) :: f0,f1,qnorm
!     integer :: a1,a2,i,j
!
!     ! Skip for Gamma, do nothing
!     qnorm=lo_sqnorm(qpoint%r)
!     if ( qnorm .lt. lo_sqtol ) then
!         qv=[1.0_r8,0.0_r8,0.0_r8]
!     else
!         qv=qpoint%r
!     endif
!
!     ! The denominator
!     f0=0.0_r8
!     do i=1,3
!     do j=1,3
!         f0=f0+qv(i)*fc%loto%eps(i,j)*qv(j)
!     enddo
!     enddo
!     f0=(4.0_r8*lo_pi/p%volume)/f0
!
!     ! The dynamical matrix at q->0
!     call mem%allocate(D,[3,3,fc%na,fc%na],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
!     D=0.0_r8
!     do a1=1,fc%na
!     do a2=1,fc%na
!         do i=1,3
!         do j=1,3
!             D(i,j,a1,a2)=dot_product(qv,fc%atom(a1)%Z(:,i))*dot_product(qv,fc%atom(a2)%Z(:,j))*f0
!         enddo
!         enddo
!     enddo
!     enddo
!
!     ! fix the normalization
!     D=D/product(fc%paddim)
!     ! Make sure it's completely gone at the zone edge. There might be some aliasing errors. Probably
!     ! mostly cosmetic, but feels bad. This gives me a 0-1 thingy:
!     f0=1.0_r8-p%bz%distance_to_zone_edge(qpoint%r)
!     ! Below a magic scaling factor, I want it to be 0
!     f1=0.15_r8
!     if ( f0 .lt. f1 ) then
!         f0=0.0_r8
!     else
!         f0=(f0-f1)/(1.0_r8-f1)
!     endif
!     ! Now cosine the thing
!     f0=lo_chop((cos(f0*lo_pi)+1.0_r8)*0.5_r8,lo_sqtol)
!     ! And scale the whole thing
!     D=D*f0
!
!     ! Add it to the forceconstant
!     if ( qnorm .gt. lo_sqtol ) then
!         do a1=1,fc%na
!         do i=1,fc%atom(a1)%n
!             a2=fc%atom(a1)%pair(i)%i2
!             fc%atom(a1)%pair(i)%adjusted_m=fc%atom(a1)%pair(i)%m+D(:,:,a1,a2)*fc%atom(a1)%pair(i)%weight
!         enddo
!         enddo
!     else
!         do a1=1,fc%na
!         do i=1,fc%atom(a1)%n
!             a2=fc%atom(a1)%pair(i)%i2
!             fc%atom(a1)%pair(i)%adjusted_m=fc%atom(a1)%pair(i)%m
!         enddo
!         enddo
!     endif
!     call mem%deallocate(D,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
! end subroutine

!> figure out which the degenerate subspaces are.
subroutine return_degenerate_subspaces(omegasquare, n_subspace, subspace_ctr, subspace_mode, mem)
    !> frequencies, sorted by size
    real(r8), dimension(:), intent(in) :: omegasquare
    !> how many degenerate subspaces are there?
    integer, intent(out) :: n_subspace
    !> how many modes per subspace
    integer, dimension(:), allocatable, intent(out) :: subspace_ctr
    !> which modes in each subspace
    integer, dimension(:, :), allocatable, intent(out) :: subspace_mode
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:), allocatable :: x
    integer, dimension(:, :), allocatable :: dj
    integer, dimension(:), allocatable :: di
    integer :: i, j, k, nb

    ! Shorthand for number of modes
    nb = size(omegasquare)
    call mem%allocate(x, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    x = lo_negsqrt(omegasquare)
    ! First a quick test, if there is nothing degenerate we can exit quickly.
    j = 0
    do i = 1, nb - 1
        if (x(i + 1) - x(i) .lt. lo_freqtol) j = j + 1
    end do
    if (j .eq. 0) then
        ! No degeneracies, nothing to worry about!
        n_subspace = -1
        call mem%deallocate(x, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        return
    end if

    ! This is a degenerate point, we have to be a little more careful.
    x = lo_negsqrt(omegasquare)
    call mem%allocate(di, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(dj, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    di = 0
    dj = 0

    n_subspace = 0
    il: do i = 1, nb
        ! Check if mode belongs to any of the existing subspaces?
        do j = 1, n_subspace
        do k = 1, di(j)
            if (abs(x(dj(k, j)) - x(i)) .lt. lo_freqtol) then
                ! we are in this subspace! increment the counter
                ! for this subspace, and store the mode.
                di(j) = di(j) + 1
                dj(di(j), j) = i
                cycle il
            end if
        end do
        end do
        ! if we make it here, we have a new subspace.
        n_subspace = n_subspace + 1
        di(n_subspace) = di(n_subspace) + 1
        dj(1, n_subspace) = i
    end do il

    i = maxval(di)
    call mem%allocate(subspace_ctr, n_subspace, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(subspace_mode, [i, n_subspace], persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
    subspace_ctr = di(1:n_subspace)
    subspace_mode = dj(1:i, 1:n_subspace)
    ! And clean the temporary things
    call mem%deallocate(x, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

end submodule
