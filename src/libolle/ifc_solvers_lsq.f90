#include "precompilerdefinitions"
submodule(ifc_solvers) ifc_solvers_lsq
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_distancetable, only: lo_distancetable
use quadratic_programming, only: lo_solve_quadratic_program
implicit none

!> Very specialized parallel matrix multiplication thingy
type, abstract :: lo_distmatrix
    ! Block size
    integer :: nb = -lo_hugeint
    ! Total number of rows
    integer :: nrow = -lo_hugeint
    ! Total nomber of columns
    integer :: ncol = -lo_hugeint
end type
! This guy holds (j1:j2, : ) parts of a matrix
type, extends(lo_distmatrix) :: lo_distmatrix_row
    ! How many columns
    integer :: nj = -lo_hugeint
    ! Which columns, in blocks of n
    integer, dimension(:), allocatable :: rowind
    ! The actual columns
    real(r8), dimension(:, :, :), allocatable :: rowbuf
end type

contains

!> Normal least squares solution
module subroutine lsq_solve(map, uc, ss, ih, tp, force_stable, usereference, maxdisp, temperature, mw, mem, verbosity)
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> helper with coefficient matrices
    type(lo_ifc_helper), intent(inout) :: ih
    !> helper with distributions
    type(lo_trainpred_helper), intent(in) :: tp
    !> Force solution table?
    logical, intent(in) :: force_stable
    !> Use reference positions for second order?
    logical, intent(in) :: usereference
    !> Largest displacement allowed
    real(r8), intent(in) :: maxdisp
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> Talk a lot
    integer, intent(in) :: verbosity

    real(r8) :: timer, t0, t1
    integer :: div, solrnk, nf

    timer = walltime()
    t0 = timer
    t1 = timer

    ! Just set some basic things
    init: block
        ! Set forces to something to start from
        if (tp%nt .gt. 0) then
            ih%f1 = ih%f0
            ih%f2 = ih%f0
            ih%f3 = ih%f0
            ih%f4 = ih%f0
        end if
        ! number of forces, local
        nf = tp%nt*ss%na*3
        ! Which rank to run the serial solutions on
        solrnk = mw%n - 1
        ! And start the timer
        if (verbosity .gt. 0) then
            write (lo_iou, *) '... least squares solver'
        end if
    end block init

    if (map%have_fc_pair) then
        second: block
            type(lo_scaling_and_product) :: asc
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: wB, wC
            real(r8), dimension(6) :: stat
            integer :: div, nrow

            ! Normal least squares for the full solution
            call asc%generate(ih%cm2, ih%f2, nf, map%xuc%nx_fc_pair, mw, &
                              noscale=.true., usereference=usereference, n_atom=ss%na)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares( &
                    asc%ATA, asc%ATB, map%xuc%x_fc_pair, map%constraints%eq2, map%constraints%d2, map%constraints%neq2, gramified=.true.)
            end if
            call asc%destroy()
            call mw%bcast(map%xuc%x_fc_pair, solrnk, __FILE__, __LINE__)

            if (force_stable) then
                ! Now we might want to force it into a stable solution, complicated thing
                call force_solution_stable_fifth_time( &
                    map, ih, tp, solrnk, uc, ss, stat, temperature, maxdisp, mw, mem, verbosity)
            else
                ! Normal least squares for the full solution
                call asc%generate(ih%cm2, ih%f2, nf, map%xuc%nx_fc_pair, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares( &
                        asc%ATA, asc%ATB, map%xuc%x_fc_pair, map%constraints%eq2, map%constraints%d2, map%constraints%neq2, gramified=.true.)
                end if
                call asc%destroy()
                call mw%bcast(map%xuc%x_fc_pair, solrnk, __FILE__, __LINE__)
            end if

            ! Then the cross-checking guys
            do div = 1, n_division
                call tp%collect(map, 2, div, .true., nrow, mem, ih=ih, coeff=wA, values=wB)
                call asc%generate(wA, wB, nrow, map%xuc%nx_fc_pair, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares( &
                        asc%ATA, asc%ATB, ih%dx2(:, div), map%constraints%eq2, map%constraints%d2, map%constraints%neq2, gramified=.true.)
                end if
                call asc%destroy()
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            end do
            call mw%bcast(ih%dx2, solrnk, __FILE__, __LINE__)

            ! And subtract forces from the solution
            if (tp%nt .gt. 0) then
                call mem%allocate(wC, nf, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                wC = 0.0_r8
                call lo_gemv(ih%cm2, map%xuc%x_fc_pair, wC)
                ih%f1 = ih%f1 - wC
                ih%f2 = wC
                ih%f3 = ih%f3 - wC
                ih%f4 = ih%f4 - wC
                call mem%deallocate(wC, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            end if
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... solved second order (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block second
    else
        ! No second order, set it to nothing.
        ih%f2 = 0.0_r8
    end if

    if (map%have_fc_singlet) then
        first: block
            type(lo_scaling_and_product) :: asc
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: wB, wC
            integer :: nrow

            ! Normal least squares for the full solution
            call asc%generate(ih%cm1, ih%f1, nf, map%xuc%nx_fc_singlet, mw, noscale=.true.)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares( &
                    asc%ATA, asc%ATB, map%xuc%x_fc_singlet, map%constraints%eq1, map%constraints%d1, map%constraints%neq1, gramified=.true.)
            end if
            call mw%bcast(map%xuc%x_fc_singlet, solrnk, __FILE__, __LINE__)
            call asc%destroy()

            ! Then the cross-checking guys
            do div = 1, n_division
                call tp%collect(map, 1, div, .true., nrow, mem, ih=ih, coeff=wA, values=wB)
                call asc%generate(wA, wB, nrow, map%xuc%nx_fc_singlet, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares( &
                        asc%ATA, asc%ATB, ih%dx1(:, div), map%constraints%eq1, map%constraints%d1, map%constraints%neq1, gramified=.true.)
                end if
                call asc%destroy()
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            end do
            call mw%bcast(ih%dx1, solrnk, __FILE__, __LINE__)

            ! And subtract forces from the solution
            if (tp%nt .gt. 0) then
                call mem%allocate(wC, nf, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                wC = 0.0_r8
                call lo_gemv(ih%cm1, map%xuc%x_fc_singlet, wC)
                ih%f1 = wC
                !ih%f2=ih%f2-wC
                ih%f3 = ih%f3 - wC
                ih%f4 = ih%f4 - wC
                call mem%deallocate(wC, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            end if
            if (verbosity .gt. 0) write (*, *) '... solved first order (', tochar(walltime() - timer), 's)'
        end block first
    else
        ! No first order, set it to nothing.
        ih%f1 = 0.0_r8
    end if

    if (map%have_fc_triplet) then
        third: block
            type(lo_scaling_and_product) :: asc
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: wB, wC
            integer :: nrow

            ! Normal least squares for the full solution
            call asc%generate(ih%cm3, ih%f3, nf, map%xuc%nx_fc_triplet, mw, noscale=.true.)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares( &
                    asc%ATA, asc%ATB, map%xuc%x_fc_triplet, map%constraints%eq3, map%constraints%d3, map%constraints%neq3, gramified=.true.)
            end if
            call asc%destroy()
            call mw%bcast(map%xuc%x_fc_triplet, solrnk, __FILE__, __LINE__)

            ! Then the cross-checking guys
            do div = 1, n_division
                call tp%collect(map, 3, div, .true., nrow, mem, ih=ih, coeff=wA, values=wB)
                call asc%generate(wA, wB, nrow, map%xuc%nx_fc_triplet, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares( &
                        asc%ATA, asc%ATB, ih%dx3(:, div), map%constraints%eq3, map%constraints%d3, map%constraints%neq3, gramified=.true.)
                end if
                call asc%destroy()
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            end do
            call mw%bcast(ih%dx3, solrnk, __FILE__, __LINE__)

            ! And subtract forces from the solution
            if (tp%nt .gt. 0) then
                call mem%allocate(wC, nf, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                wC = 0.0_r8
                call lo_gemv(ih%cm3, map%xuc%x_fc_triplet, wC)
                ih%f3 = wC
                ih%f4 = ih%f4 - wC
                call mem%deallocate(wC, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            end if
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... solved third order (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block third
    else
        ! No third order, set it to nothing.
        ih%f3 = 0.0_r8
    end if

    if (map%have_fc_quartet) then
        fourth: block
            type(lo_scaling_and_product) :: asc
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: wB, wC
            integer :: nrow
            timer = walltime()

            ! Normal least squares for the full solution
            call asc%generate(ih%cm4, ih%f4, nf, map%xuc%nx_fc_quartet, mw, noscale=.true.)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares( &
                    asc%ATA, asc%ATB, map%xuc%x_fc_quartet, map%constraints%eq4, map%constraints%d4, map%constraints%neq4, gramified=.true.)
            end if
            call asc%destroy()
            call mw%bcast(map%xuc%x_fc_quartet, solrnk, __FILE__, __LINE__)

            ! Then the cross-checking guys
            do div = 1, n_division
                call tp%collect(map, 4, div, .true., nrow, mem, ih=ih, coeff=wA, values=wB)
                call asc%generate(wA, wB, nrow, map%xuc%nx_fc_quartet, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares( &
                        asc%ATA, asc%ATB, ih%dx4(:, div), map%constraints%eq4, map%constraints%d4, map%constraints%neq4, gramified=.true.)
                end if
                call asc%destroy()
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            end do
            call mw%bcast(ih%dx4, solrnk, __FILE__, __LINE__)

            ! And subtract forces from the solution
            if (tp%nt .gt. 0) then
                call mem%allocate(wC, nf, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
                wC = 0.0_r8
                call lo_gemv(ih%cm4, map%xuc%x_fc_quartet, wC)
                ih%f4 = wC
                call mem%deallocate(wC, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            end if
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... solved fourth order (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block fourth
    else
        ! No fourth order, set it to nothing.
        ih%f4 = 0.0_r8
    end if
end subroutine

! below are not exposed

!> Try another way to get a stable solution
subroutine force_solution_stable_fifth_time( &
    map, ih, tp, solrnk, uc, ss, stat, temperature, max_displacement, mw, mem, verbosity)
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> helper with data
    type(lo_ifc_helper), intent(in) :: ih
    !> helper with distributions
    type(lo_trainpred_helper), intent(in) :: tp
    !> rank for serial solution
    integer, intent(in) :: solrnk
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> How did it go?
    real(r8), dimension(6), intent(out) :: stat
    !> Maybe provide a temperature, otherwise guess
    real(r8), intent(in) :: temperature
    !> Largest allowed displacement
    real(r8), intent(in) :: max_displacement
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> verbosity
    integer, intent(in) :: verbosity

    !type(lo_forceconstant_secondorder) :: fc
    type(lo_distmatrix_row) :: dB
    real(r8), dimension(:, :), allocatable :: dynmat, lrdynmat, eigenvectors
    real(r8), dimension(:), allocatable :: eigenvalues, ifc
    real(r8) :: t0, t1
    integer :: nx

    ! start timer
    t0 = walltime()
    t1 = t0

    ! Set some basic things
    init: block
        ! Set the status to nothing
        stat = -5.0_r8
        ! Number of independent
        nx = map%xuc%nx_fc_pair
        ! Temporary space for forceconstants
        call mem%allocate(ifc, nx, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        ifc = 0.0_r8
    end block init

    ! First thing first, get the normal solution
    normalsol: block
        type(lo_scaling_and_product) :: asc

        ! Normal least squares for the full solution
        call asc%generate(ih%cm2, ih%f2, tp%nt*ss%na*3, map%xuc%nx_fc_pair, mw, noscale=.true.)
        if (mw%r .eq. solrnk) then
            call lo_linear_least_squares(asc%ATA, asc%ATB, ifc, map%constraints%eq2, map%constraints%d2, map%constraints%neq2, gramified=.true.)
        end if
        call asc%destroy()
        call mw%bcast(ifc, solrnk, __FILE__, __LINE__)
    end block normalsol

    ! Prep eigenvalue solver
    prepegv: block
        integer :: a1, a2, i, j, ii, jj
        real(r8) :: f0, maxdisp

        ! Space for dynamical matrices and stuffs
        call mem%allocate(dynmat, [ss%na*3, ss%na*3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(lrdynmat, [ss%na*3, ss%na*3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(eigenvectors, [ss%na*3, ss%na*3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(eigenvalues, ss%na*3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        dynmat = 0.0_r8
        lrdynmat = 0.0_r8
        eigenvectors = 0.0_r8
        eigenvalues = 0.0_r8

        if (map%polar .gt. 0 .and. map%polarcorrectiontype .eq. 3) then
            ! get the long-ranged interactions
            do a1 = 1, ss%na
            do a2 = 1, ss%na
                f0 = ss%invsqrtmass(a1)*ss%invsqrtmass(a2)
                do i = 1, 3
                do j = 1, 3
                    ii = (a1 - 1)*3 + i
                    jj = (a2 - 1)*3 + j
                    lrdynmat(jj, ii) = ih%polar_fc(i, j, a1, a2)*f0
                end do
                end do
            end do
            end do
        else
            lrdynmat = 0.0_r8
        end if
        ! Calculate eigenvalues and check if I need to do something about them
        call lo_supercell_dynmat_from_irreducible_forceconstant(map, uc, ss, ifc, dynmat, mw, mem, shortrange=.true.)
        dynmat = dynmat + lrdynmat
        if (mw%r .eq. solrnk) then
            call lo_real_symmetric_eigenvalues_eigenvectors(dynmat, eigenvalues, eigenvectors, nzeros=3, tolerance=1E-13_r8)
        else
            eigenvalues = -lo_huge
        end if
        call mw%bcast(eigenvalues, solrnk)
        call mw%bcast(eigenvectors, solrnk)
        stat(1) = lo_negsqrt(minval(eigenvalues))
        stat(2) = lo_negsqrt(maxval(eigenvalues))

        ! Calculate the largest displacements
        call lo_thermal_displacement_matrix_commensurate( &
            ss, uc, eigenvectors, lo_negsqrt(eigenvalues), temperature, lo_tiny, maxdisp=maxdisp)

        ! Check if we are fine?
        if (maxdisp .lt. max_displacement) then
            ! yup, we are fine
            if (verbosity .gt. 0) then
                write (lo_iou, *) '... forceconstants already stable'
            end if
            map%xuc%x_fc_pair = ifc
            stat(3) = 0.0_r8
            stat(4) = lo_huge
            do i = 1, size(eigenvalues)
                if (abs(eigenvalues(i)) .gt. lo_tiny) stat(4) = min(stat(4), eigenvalues(i))
            end do
            stat(4) = lo_negsqrt(stat(4))
            call mem%deallocate(ifc, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(lrdynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(eigenvectors, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(eigenvalues, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            ! NOTE early return here, in case everything is already fine.
            return
        else
            ! Get the huge, explicit coefficient matrix, now I need it.
            f0 = lo_negsqrt(minval(eigenvalues))*lo_frequency_Hartree_to_Thz
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (*, *) '... found unstable mode: ', tochar(f0), 'THz (', tochar(t1 - t0), 's)'
                t0 = t1
            end if

            call get_large_B_matrix(map, mw, dB, verbosity)
        end if
    end block prepegv

    ! Start iterating, maybe?
    qprg: block
        real(r8), parameter :: alpha0 = (0.001_r8*lo_frequency_THz_to_Hartree)**2 ! truly magic number.
        integer, parameter :: nmaxbadstep = 4                 ! Largest number of bad steps
        integer, parameter :: maxiter = 100
        type(lo_scaling_and_product) :: sc
        real(r8), dimension(:, :), allocatable :: inequal_constraints
        real(r8), dimension(:), allocatable :: inequal_alpha, equality_b, dv, cvs, xv, xw, xbest
        real(r8) :: alpha, maxdisp, f0, tt0, tt1
        real(r8) :: timer_eig, timer_constr, timer_solve, bestminom
        integer :: maxnconstraints, badstep
        integer :: nics, iter
        integer :: i, j, l, i1, a1, a2

        if (verbosity .gt. 0) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'Fixing instabilities:'
            write (lo_iou, *) ' Iter       min(w)   Nneg     max(u)    Ncon  Nbad    Time/step(s)'
        end if

        tt0 = walltime()
        ! Get the gramian and stuff. Can use the same routine as always.
        call sc%generate(ih%cm2, ih%f2, tp%nt*ss%na*3, nx, mw, noscale=.true.)

        ! Decide on a reasonable number of extra constraints?
        maxnconstraints = maxiter*nx

        call mem%allocate(inequal_constraints, [maxnconstraints, nx], &
                          persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(inequal_alpha, maxnconstraints, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(dv, dB%nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(cvs, nx, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(xv, nx, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(xw, nx, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(xbest, nx, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        inequal_constraints = 0.0_r8
        inequal_alpha = 0.0_r8
        dv = 0.0_r8
        nics = 0
        if (map%constraints%neq2 .gt. 0) then
            call mem%allocate(equality_b, map%constraints%neq2, &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            equality_b = 0.0_r8
        else
            call mem%allocate(equality_b, 1, &
                              persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            equality_b = -lo_huge
        end if
        ! Reset some timers so that I know what to improve
        timer_eig = 0.0_r8
        timer_constr = 0.0_r8
        timer_solve = 0.0_r8

        ! Something to start from
        xw = ifc
        xv = ifc
        xbest = ifc
        bestminom = -lo_huge
        badstep = 0
        ! Current min eigenvalue, alpha
        alpha = alpha0
        iterloop: do iter = 1, maxiter !100
            t0 = walltime()
            ! Get the current frequencies and stuff
            call lo_supercell_dynmat_from_irreducible_forceconstant(map, uc, ss, xv, dynmat, mw, mem, shortrange=.true.)
            dynmat = dynmat + lrdynmat
            if (mw%r .eq. solrnk) then
                call lo_real_symmetric_eigenvalues_eigenvectors(dynmat, eigenvalues, eigenvectors, nzeros=3, tolerance=1E-15_r8)
                if (minval(eigenvalues) .gt. -lo_tiny) then
                    ! Calculate current largest displacement
                    call lo_thermal_displacement_matrix_commensurate(ss, uc, eigenvectors, lo_negsqrt(eigenvalues), &
                                                                     temperature, lo_tiny, maxdisp=maxdisp)
                else
                    ! Or some other larger number
                    maxdisp = 1000.0_r8
                end if
            end if
            call mw%bcast(eigenvalues, solrnk)
            call mw%bcast(eigenvectors, solrnk)
            call mw%bcast(maxdisp, solrnk)

            ! Store an emergency solution, in case all else fails
            if (minval(eigenvalues) .gt. bestminom) then
                bestminom = minval(eigenvalues)
                xbest = xv
            end if

            timer_eig = timer_eig + walltime() - t0

            ! First check that I only got three acoustic modes
            j = count(abs(eigenvalues) < lo_tiny)
            if (j .gt. 3) then
                ! Ok this means I got more than three translational modes. This is just
                ! bad luck, change alpha a little to push them to either side.
                alpha = alpha*2
                xv = xw
                cycle iterloop
            end if

            ! Report how it is going
            if (verbosity .gt. 0) then
                f0 = lo_negsqrt(minval(eigenvalues))*lo_frequency_Hartree_to_THz
                tt1 = walltime()
                write (*, "(1X,I4,1X,F13.4,1X,I4,1X,F13.4,1X,I5,3X,I2,2X,F12.3)") &
                    iter, f0, count(eigenvalues < -lo_tiny), maxdisp, nics, badstep, tt1 - tt0
                tt0 = tt1
            end if

            ! Check that the MSD is not too large
            if (maxdisp .gt. max_displacement) then
                ! Increase alpha
                alpha = alpha*3
            else
                ! This might be converged!
                if (minval(eigenvalues) .gt. -lo_tiny) then
                    ! Yup, converged!
                    if (verbosity .gt. 0) then
                        write (*, *) '... converged to feasible solution'
                    end if
                    stat(3) = 1.0_r8
                    stat(4) = lo_negsqrt(minval(eigenvalues))
                    call mw%check_and_sync(stat, 1E-12_r8, 0, 2, 'stat', __FILE__, __LINE__)
                    exit iterloop
                end if
            end if

            ! Ok, time to create more constraints!
            t0 = walltime()
            do l = 1, ss%na*3 ! Or some other magic number. This seems to be the safe choice.
                ! skip acoustic
                if (abs(eigenvalues(l)) .lt. lo_tiny) cycle
                ! Only add constraints to modes that are close to being unstable
                if (eigenvalues(l) .gt. alpha*100) cycle
                ! Now build the magical constraint. Hmmm.
                cvs = 0.0_r8
                do i1 = 1, dB%nj
                    ! This is the eigenvector outer product.
                    i = dB%rowind(i1)
                    do j = 1, dB%nb
                        a1 = (i - 1)/3 + 1
                        a2 = (j - 1)/3 + 1
                        dv(j) = eigenvectors(i, l)*eigenvectors(j, l)*ss%invsqrtmass(a1)*ss%invsqrtmass(a2)
                    end do
                    call lo_gemv(dB%rowbuf(:, :, i1), dv, cvs, trans='T', beta=1.0_r8)
                end do
                call mw%allreduce('sum', cvs)
                nics = nics + 1

                ! Idiot test 1, in case we have to abandon ship.
                if (nics .gt. maxnconstraints) then
                    stat(3) = -1.0_r8
                    xv = xbest
                    call mw%check_and_sync(stat, 1E-12_r8, 0, 2, 'stat', __FILE__, __LINE__)
                    exit iterloop ! ABANDON SHIP!!
                end if

                inequal_constraints(nics, :) = cvs
                call lo_gemv(lrdynmat, eigenvectors(:, l), dv)
                inequal_alpha(nics) = alpha - dot_product(eigenvectors(:, l), dv)
            end do
            timer_constr = timer_constr + walltime() - t0
            t0 = walltime()

            ! Actually solve for things
            xv = 0.0_r8
            if (mw%r .eq. solrnk) then
                call lo_solve_quadratic_program(sc%ATA, sc%ATB, xv, map%constraints%eq2, inequal_constraints, &
                                                equality_b, inequal_alpha, map%constraints%neq2, nics, verbosity - 10, 1E-13_r8)
            end if
            call mw%bcast(xv, solrnk)
            timer_solve = timer_solve + walltime() - t0

            ! I will slowly increase alpha, think that is reasonable
            alpha = alpha*1.5_r8

            ! Idiot test 2
            if (iter .eq. maxiter) then
                stat(3) = -1.0_r8
                xv = xbest
                call mw%check_and_sync(stat, 1E-12_r8, 0, 2, 'stat', __FILE__, __LINE__)
                exit iterloop
            end if

            ! Idiot test 3
            if (minval(eigenvalues) .lt. 2*bestminom) then
                badstep = badstep + 1
            else
                ! Not sure about this. Maybe a good idea.
                badstep = 0
            end if
            if (badstep .gt. nmaxbadstep) then
                stat(3) = -1.0_r8
                xv = xbest
                call mw%check_and_sync(stat, 1E-12_r8, 0, 2, 'stat', __FILE__, __LINE__)
                exit iterloop
            end if

        end do iterloop
        ! Return solution?
        stat(4) = lo_huge
        do i = 1, size(eigenvalues)
            if (abs(eigenvalues(i)) .gt. lo_tiny) stat(4) = min(stat(4), eigenvalues(i))
        end do
        stat(4) = lo_negsqrt(stat(4))
        call mpi_barrier(mw%comm, mw%error)
        map%xuc%x_fc_pair = xv
        call mw%check_and_sync(map%xuc%x_fc_pair, 1E-14_r8, 0, 2, 'ifc', __FILE__, __LINE__)
        call mw%check_and_sync(stat, 1E-12_r8, 0, 2, 'stat', __FILE__, __LINE__)

        ! Cleanup
        call mem%deallocate(ifc, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(lrdynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(eigenvectors, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(eigenvalues, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(inequal_constraints, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(inequal_alpha, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(cvs, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(xv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(xw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(xbest, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(equality_b, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        !  Report how it went?
        if (verbosity .gt. 0) then
            write (*, *) '... stored solution'
            write (*, *) '   timer eigenvalue problem:', timer_eig
            write (*, *) '    timer build constraints:', timer_constr
            write (*, *) '    timer quadratic program:', timer_solve
            !write(*,*) '                      total:',walltime()-timer
        end if
    end block qprg
end subroutine

!> Get some diagnostics from the thermal displacement matrix
subroutine lo_thermal_displacement_matrix_commensurate(ss, uc, eigenvectors, omega, temperature, thres, sigma, maxdisp)
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> eigenvectors
    real(r8), dimension(:, :), intent(in) :: eigenvectors
    !> frequencies
    real(r8), dimension(:), intent(in) :: omega
    !> temperature
    real(r8), intent(in) :: temperature
    !> threshold to replace imaginary frequencies with
    real(r8), intent(in), optional :: thres
    !> displacement matrix
    real(r8), dimension(:, :, :), intent(out), optional :: sigma
    !> max spherically averaged mean displacement
    real(r8), intent(out), optional :: maxdisp

    type(lo_distancetable) :: dt
    real(r8), dimension(3, 3) :: m0, m1
    real(r8), dimension(3) :: v0
    real(r8) :: imfreqval, f0, om
    integer :: nuc, uca, i, a1

    !@todo add sanity test
    nuc = maxval(ss%info%index_in_unitcell)

    if (present(thres)) then
        imfreqval = thres
    else
        imfreqval = 0.0_r8
    end if

    if (present(sigma)) sigma = 0.0_r8
    if (present(maxdisp)) then
        call dt%generate(uc%r, uc%latticevectors, ss%nearest_neighbour_distance()*3, verbosity=0)
        maxdisp = 0.0_r8
    end if
    do uca = 1, nuc
        m0 = 0.0_r8
        a1 = -1
        do i = 1, ss%na
            if (ss%info%index_in_unitcell(i) .eq. uca) then
                a1 = i
                exit
            end if
        end do
        ! Now sum over all modes
        do i = 1, ss%na*3
            om = omega(i)
            if (om .lt. -lo_tiny) om = imfreqval
            if (abs(om) .lt. lo_tiny) cycle
            f0 = 0.5_r8*(1.0_r8 + 2.0_r8*lo_planck(temperature, om))/om
            v0 = eigenvectors((a1 - 1)*3 + 1:a1*3, i)
            m0 = m0 + lo_outerproduct(v0, v0)*f0
        end do
        ! And the mass
        m0 = m0/ss%mass(a1)
        ! Maybe store
        if (present(sigma)) sigma(:, :, uca) = m0
        if (present(maxdisp)) then
            ! Use the largest mean square displacement and calculate
            ! a spherically averaged mean displacement from this.
            !call lo_symmetric_eigensystem_3x3matrix(m0,v0,m1)
            !f0=sqrt(maxval(v0))*gammafactor
            !maxdisp=max(maxdisp,f0)

            ! Not spherically averaged, look in the relevant directions instead.
            m1 = lo_invert3x3matrix(m0)
            do i = 2, dt%particle(uca)%n
                v0 = dt%particle(uca)%v(:, i)/dt%particle(uca)%d(i)
                f0 = 1.0_r8/sqrt(dot_product(v0, matmul(m1, v0)))
                f0 = f0/dt%particle(uca)%d(i)
                maxdisp = max(maxdisp, f0)
            end do
        end if
    end do
end subroutine

!> Calculate the explicit, huge, coefficient matrix, in distributed form.
subroutine get_large_B_matrix(map, mw, dB, verbosity)
    !> forcemap
    type(lo_forcemap), intent(in) :: map
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> distributed coefficient matrix
    type(lo_distmatrix_row), intent(out) :: dB
    !> talk a lot
    integer, intent(in) :: verbosity

    real(r8) :: timer
    integer, dimension(:, :, :), allocatable :: dmi

    ! Start the timer
    timer = walltime()

    ! Set some prep things
    init: block
        integer :: l1, l2, a1, a2, i1, i2, i, j, k

        ! Figure out indices. Hurts my brain.
        lo_allocate(dmi(9, map%n_atom_ss, map%n_atom_ss))
        dmi = 0
        l1 = 0
        do a1 = 1, map%n_atom_ss
        do i1 = 1, 3
            l1 = l1 + 1
            l2 = 0
            do a2 = 1, map%n_atom_ss
            do i2 = 1, 3
                l2 = l2 + 1
                j = (i1 - 1)*3 + i2
                k = (l1 - 1)*map%n_atom_ss*3 + l2
                dmi(j, a1, a2) = k
                ! This is the flattening I hope to achieve dynmat(l2,l1)=dynmat(i1,i2,a1,a2)
            end do
            end do
        end do
        end do

        ! Make space for the distributed matrix.
        dB%nb = map%n_atom_ss*3
        dB%nrow = (map%n_atom_ss*3)**2
        dB%ncol = map%xuc%nx_fc_pair
        ! Then specifically for this rank.
        dB%nj = 0
        do i = 1, dB%nb
            if (mod(i, mw%n) .ne. mw%r) cycle
            dB%nj = dB%nj + 1
        end do
        lo_allocate(dB%rowind(db%nj))
        lo_allocate(dB%rowbuf(dB%nb, dB%ncol, db%nj))
        dB%rowind = 0
        dB%rowbuf = 0.0_r8
        j = 0
        do i = 1, dB%nb
            if (mod(i, mw%n) .ne. mw%r) cycle
            j = j + 1
            dB%rowind(j) = i
        end do

        if (verbosity .gt. 0) then
            write (*, *) '... sorted out distribution of coefficient matrix'
            call lo_progressbar_init()
        end if
    end block init

    ! Ok. Now I know which parts of the matrix is supposed to end up on this rank.
    ! Go through the entire thing and put things where they belong.
    sortstuff: block
        real(r8), dimension(9, 9) :: dm0
        integer, dimension(9) :: thetaind, rowind
        integer :: a1, a2, j, k, unsh, unop, ntheta, ipair
        integer :: irow, minind, maxind

        do irow = 1, db%nj
            ! Ok. Now I know which parts of the matrix is supposed to end up on this rank.
            ! Lowest and highest row indices for this part
            minind = (db%rowind(irow) - 1)*db%nb + 1
            maxind = db%rowind(irow)*db%nb
            do ipair = 1, map%xss%n_fc_pair
                unsh = map%xss%ind_fc_pair(1, ipair)
                unop = map%xss%ind_fc_pair(2, ipair)
                a1 = map%xss%ind_fc_pair(3, ipair)
                a2 = map%xss%ind_fc_pair(4, ipair)
                ntheta = map%fc_pair_shell(unsh)%nx
                if (ntheta .eq. 0) cycle
                rowind = dmi(:, a1, a2)
                if (maxval(rowind) .ge. minind) then
                if (minval(rowind) .le. maxind) then
                    thetaind(1:ntheta) = map%fc_pair_shell(unsh)%ind_global
                    dm0(:, 1:ntheta) = matmul(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff)
                    ! Now figure out where to put these. A little tricky.
                    do j = 1, 9
                        if (rowind(j) .ge. minind .and. rowind(j) .le. maxind) then
                            ! To local index
                            k = rowind(j) - (db%rowind(irow) - 1)*db%nb
                            dB%rowbuf(k, thetaind(1:ntheta), irow) = dB%rowbuf(k, thetaind(1:ntheta), irow) + dm0(j, 1:ntheta)
                        end if
                    end do
                end if
                end if
                ! Then try to add the self-term
                rowind = dmi(:, a1, a1)
                if (maxval(rowind) .ge. minind) then
                if (minval(rowind) .le. maxind) then
                    thetaind(1:ntheta) = map%fc_pair_shell(unsh)%ind_global
                    dm0(:, 1:ntheta) = matmul(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff)
                    ! Now figure out where to put these. A little tricky.
                    do j = 1, 9
                        if (rowind(j) .ge. minind .and. rowind(j) .le. maxind) then
                            ! To local index
                            k = rowind(j) - (db%rowind(irow) - 1)*db%nb
                            dB%rowbuf(k, thetaind(1:ntheta), irow) = dB%rowbuf(k, thetaind(1:ntheta), irow) - dm0(j, 1:ntheta)
                        end if
                    end do
                end if
                end if
            end do
            if (verbosity .gt. 0) then
                call lo_progressbar(' ... coefficient matrix', irow, db%nj, walltime() - timer)
            end if
        end do
    end block sortstuff
    ! And a little cleanup
    lo_deallocate(dmi)
end subroutine

end submodule
