submodule(ifc_solvers) ifc_solvers_diel
!!
!! Solve for dielectric things
!!
use konstanter, only: lo_pi
use gottochblandat, only: lo_compress_equations
use type_blas_lapack_wrappers, only: lo_dgelss, lo_dgglse
use type_forcemap, only: lo_coeffmatrix_unitcell_Z_singlet
implicit none

contains

!> normal least squares solution for eps and z to many orders
module subroutine lsq_solve_diel(dh, tp, map, uc, ss, mw, mem, verbosity)
    !> coefficient matrices
    type(lo_dielectric_helper), intent(inout) :: dh
    !> distribution information
    type(lo_trainpred_helper), intent(in) :: tp
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> Talk a lot
    integer, intent(in) :: verbosity

    integer :: solrnk
    real(r8) :: timer, t0, t1

    timer = walltime()
    t0 = timer
    t1 = timer

    ! Just set some basic things
    init: block
        ! Set the measured values to something to start from
        if (tp%nt .gt. 0) then
            dh%Z1 = dh%Z0
            dh%Z2 = dh%Z0
            dh%Z3 = dh%Z0
            dh%eps0 = dh%epso
            dh%eps1 = dh%epso
            dh%eps2 = dh%epso
        end if
        ! Which rank to run the serial solutions on
        solrnk = mw%n - 1
        ! And start the timer
        if (verbosity .gt. 0 .and. mw%talk) then
            write (*, *) '... least squares solver for dielectric tensor and Born charges'
        end if
    end block init

    ! First sort out the dielectric tensor. This goes slightly backwards, but it
    ! makes some sense I promise. Should be done in this order:
    !
    ! 1) Fit eps singlet, subtract
    ! 2) Fit eps pair, subtract
    ! 3) Fit eps global.
    !

    ! First get the normal Born-charges
    if (map%polar .gt. 0) then
    if (map%xuc%nx_eps_global .gt. 0) then
        epsglob: block
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: wB
            real(r8), dimension(3, 3) :: m0
            type(lo_scaling_and_product) :: asc
            integer :: div, nrow, i

            ! Normal least squares, for the raw second derivatives
            call asc%generate(dh%cm_eps0, dh%eps0, tp%nt*9, map%xuc%nx_eps_global, mw, noscale=.true.)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares(asc%ATA, asc%ATB, map%xuc%x_eps_global_deriv, gramified=.true.)
            end if
            call asc%destroy()
            call mw%bcast(map%xuc%x_eps_global_deriv, solrnk, __FILE__, __LINE__)

            ! Now, to confuse everyone, I need the irreducible representation for eps and not the
            ! raw derivatives I just calculated. This is quite fast and easy.
            call mem%allocate(wA, [9, map%xuc%nx_eps_global], persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%allocate(wB, 9, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            wA = 0.0_r8
            wB = 0.0_r8
            wB = matmul(map%eps_global_shell%coeff, map%xuc%x_eps_global_deriv)
            m0 = lo_unflatten_2tensor(wB)
            m0 = -m0*4*lo_pi/ss%volume
            do i = 1, 3
                m0(i, i) = 1.0_r8 + m0(i, i)
            end do
            wA = map%eps_global_shell%coeff
            wB = lo_flattentensor(m0)
            call lo_linear_least_squares(wA, wB, map%xuc%x_eps_global)
            call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)

            ! Then the cross-checking guys
            do div = 1, n_division
                call tp%collect(map, 10, div, .true., nrow, mem, dh=dh, coeff=wA, values=wB)
                call asc%generate(wA, wB, nrow, map%xuc%nx_eps_global, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares(asc%ATA, asc%ATB, dh%dx_eps0(:, div), gramified=.true.)
                end if
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call asc%destroy()
                call mw%bcast(dh%dx_eps0(:, div), solrnk, __FILE__, __LINE__)
            end do

            ! Subtract current results from solution
            if (tp%nt .gt. 0) then
                call lo_gemv(dh%cm_eps0, map%xuc%x_eps_global_deriv, dh%eps0)
                dh%eps1 = dh%eps1 - dh%eps0
                dh%eps2 = dh%eps2 - dh%eps0
            end if

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... solved for global dielectric constant (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block epsglob
    else
        call lo_stop_gracefully(['No components of the dielectric tensor, makes no sense'], &
                                lo_exitcode_param, __FILE__, __LINE__, mw%comm)
    end if

    if (map%xuc%nx_Z_singlet .gt. 0) then
        zsing: block
            type(lo_scaling_and_product) :: asc
            real(r8), dimension(:, :), allocatable :: coeff, wA
            real(r8), dimension(:), allocatable :: wB
            real(r8), dimension(3, 3) :: eps
            integer :: div, nrow, i

            ! Normal least squares
            call asc%generate(dh%cm_Z1, dh%Z1, tp%nt*ss%na*9, map%xuc%nx_Z_singlet, mw, noscale=.true.)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares(asc%ATA, asc%ATB, map%xuc%x_Z_singlet, gramified=.true.)
                ! Make sure we are Hermitian.
                eps = lo_unflatten_2tensor(matmul(map%eps_global_shell%coeff, map%xuc%x_eps_global))
                ! Convert eps to actual eps
                eps = -eps*4*lo_pi/ss%volume
                do i = 1, 3
                    eps(i, i) = eps(i, i) + 1.0_r8
                end do
                call mem%allocate(coeff, [9*map%n_atom_uc, map%xuc%nx_Z_singlet], &
                                  persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                coeff = 0.0_r8
                call lo_coeffmatrix_unitcell_Z_singlet(map, coeff)
                call dh%ew%set(uc, eps, 2, 1E-20_r8, verbosity)
                call dh%ew%force_borncharges_Hermitian(map%xuc%x_Z_singlet, coeff, eps, uc, verbosity)
                call mem%deallocate(coeff, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            end if
            call asc%destroy()
            call mw%bcast(map%xuc%x_Z_singlet, solrnk, __FILE__, __LINE__)

            ! Cross-validation things?
            do div = 1, n_division
                call tp%collect(map, 21, div, .true., nrow, mem, dh=dh, coeff=wA, values=wB)
                call asc%generate(wA, wB, nrow, map%xuc%nx_Z_singlet, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares(asc%ATA, asc%ATB, dh%dx_Z1(:, div), gramified=.true.)
                end if
                call asc%destroy()
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            end do
            call mw%bcast(dh%dx_Z1, solrnk, __FILE__, __LINE__)

            ! Subtract current results from solution
            if (tp%nt .gt. 0) then
                call lo_gemv(dh%cm_Z1, map%xuc%x_Z_singlet, dh%Z1)
                dh%Z2 = dh%Z2 - dh%Z1
                dh%Z3 = dh%Z3 - dh%Z1
            end if

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... solved for Born charge singlets (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block zsing
    end if
    end if

    ! Now for things we might not always have.
    if (map%have_Z_pair) then
        zpair: block
            type(lo_scaling_and_product) :: asc
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: wB
            integer :: div, nrow

            ! Normal least squares
            call asc%generate(dh%cm_Z2, dh%Z2, tp%nt*map%n_atom_ss*9, map%xuc%nx_Z_pair, mw, noscale=.true.)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares(asc%ATA, asc%ATB, map%xuc%x_Z_pair, gramified=.true.)
            end if
            call mw%bcast(map%xuc%x_Z_pair, solrnk, __FILE__, __LINE__)

            ! Cross-validation things?
            do div = 1, n_division
                call tp%collect(map, 22, div, .true., nrow, mem, dh=dh, coeff=wA, values=wB)
                call asc%generate(wA, wB, nrow, map%xuc%nx_Z_pair, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares(asc%ATA, asc%ATB, dh%dx_Z2(:, div), gramified=.true.)
                end if
                call asc%destroy()
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            end do
            call mw%bcast(dh%dx_Z2, solrnk, __FILE__, __LINE__)

            ! Subtract current results from solution
            if (tp%nt .gt. 0) then
                call lo_gemv(dh%cm_Z2, map%xuc%x_Z_pair, dh%Z2)
                dh%Z3 = dh%Z3 - dh%Z2
            end if

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... solved for Born charge pairs (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block zpair
    end if

    ! Now for things we might not always have.
    if (map%have_Z_triplet) then
        ztriplet: block
            type(lo_scaling_and_product) :: asc
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: wB
            integer :: div, nrow

            ! Normal least squares
            call asc%generate(dh%cm_Z3, dh%Z3, tp%nt*map%n_atom_ss*9, map%xuc%nx_Z_triplet, mw, noscale=.true.)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares( &
                    asc%ATA, asc%ATB, map%xuc%x_Z_triplet, map%constraints%eqz3, map%constraints%neqz3, gramified=.true.)
            end if
            call mw%bcast(map%xuc%x_Z_triplet, solrnk, __FILE__, __LINE__)

            ! Cross-validation things?
            do div = 1, n_division
                call tp%collect(map, 23, div, .true., nrow, mem, dh=dh, coeff=wA, values=wB)
                call asc%generate(wA, wB, nrow, map%xuc%nx_Z_triplet, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares( &
                        asc%ATA, asc%ATB, dh%dx_Z3(:, div), map%constraints%eqz3, map%constraints%neqz3, gramified=.true.)
                end if
                call asc%destroy()
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            end do
            call mw%bcast(dh%dx_Z2, solrnk, __FILE__, __LINE__)

            ! Subtract current results from solution
            if (tp%nt .gt. 0) then
                call lo_gemv(dh%cm_Z3, map%xuc%x_Z_triplet, dh%Z3)
                !dh%Z3=dh%Z3-dh%Z2
            end if

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... solved for Born charge triplets (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block ztriplet
    end if

    if (map%xuc%nx_eps_singlet .gt. 0) then
        epssing: block
            type(lo_scaling_and_product) :: asc
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: wB
            integer :: div, nrow

            ! Normal least squares
            call asc%generate(dh%cm_eps1, dh%eps1, tp%nt*9, map%xuc%nx_eps_singlet, mw, noscale=.true.)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares(asc%ATA, asc%ATB, map%xuc%x_eps_singlet, gramified=.true.)
            end if
            call asc%destroy()
            call mw%bcast(map%xuc%x_eps_singlet, solrnk, __FILE__, __LINE__)
            ! Cross-validation things?
            do div = 1, n_division
                call tp%collect(map, 11, div, .true., nrow, mem, dh=dh, coeff=wA, values=wB)
                call asc%generate(wA, wB, nrow, map%xuc%nx_eps_singlet, mw, noscale=.true.)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares(asc%ATA, asc%ATB, dh%dx_eps1(:, div), gramified=.true.)
                end if
                call asc%destroy()
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            end do
            call mw%bcast(dh%dx_eps1, solrnk, __FILE__, __LINE__)
            ! Subtract current results from solution
            if (tp%nt .gt. 0) then
                call lo_gemv(dh%cm_eps1, map%xuc%x_eps_singlet, dh%eps1)
                dh%eps2 = dh%eps2 - dh%eps1
            end if

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... solved for dielectric charge singlets (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block epssing
    end if

    if (map%xuc%nx_eps_pair .gt. 0) then
        epspair: block
            type(lo_scaling_and_product) :: asc
            real(r8), dimension(:, :), allocatable :: wA
            real(r8), dimension(:), allocatable :: wB
            integer :: div, nrow

            ! Normal least squares
            call asc%merge(dh%cm_eps2, dh%eps2, tp%nt*9, map%xuc%nx_eps_pair, mw)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares(asc%ATA, asc%ATB, map%xuc%x_eps_pair)
            end if
            call asc%destroy()
            call mw%bcast(map%xuc%x_eps_pair, solrnk, __FILE__, __LINE__)
            ! Cross-validation things?
            do div = 1, n_division
                call tp%collect(map, 12, div, .true., nrow, mem, dh=dh, coeff=wA, values=wB)
                call asc%merge(wA, wB, nrow, map%xuc%nx_eps_pair, mw)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares(asc%ATA, asc%ATB, dh%dx_eps2(:, div))
                end if
                call asc%destroy()
                call mem%deallocate(wA, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
                call mem%deallocate(wB, persistent=.true., scalable=.true., file=__FILE__, line=__LINE__)
            end do
            call mw%bcast(dh%dx_eps2, solrnk, __FILE__, __LINE__)
            ! Subtract current results from solution
            if (tp%nt .gt. 0) then
                call lo_gemv(dh%cm_eps2, map%xuc%x_eps_pair, dh%eps2)
                ! dh%eps2=dh%eps2-dh%eps1
            end if

            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... solved for dielectric charge pairs (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
        end block epspair
    end if
end subroutine

!> get Born effective charges and dielectric tensor, the simple way.
module subroutine lo_solve_for_borncharges(map, p, Z, eps, filename, mw, mem, verbosity)
    !> forcemap
    class(lo_forcemap), intent(inout) :: map
    !> structure
    type(lo_crystalstructure), intent(in) :: p
    !> perhaps the born-charges are input?
    real(r8), dimension(:, :, :), intent(in), optional :: Z
    !> perhaps also the dielectric tensor?
    real(r8), dimension(3, 3), intent(in), optional :: eps
    !> filename where borncharges are stored
    character(len=*), intent(in), optional :: filename
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    real(r8), dimension(:, :, :), allocatable :: Zbuf
    real(r8), dimension(3, 3) :: epsbuf, derivbuf
    integer :: solrnk

    ! Which rank to solve on? This is so fast I do it serially.
    ! and the nonlinear crap is not parallel yet anyway.
    solrnk = mw%n - 1

    if (mw%r .eq. solrnk) then
        ! start by getting Z and eps, either from input or from a file
        readZandeps: block
            integer :: u, i, j

            call mem%allocate(Zbuf, [3, 3, map%n_atom_uc], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            Zbuf = 0.0_r8

            if (present(Z) .and. present(eps)) then
                Zbuf = Z
                epsbuf = eps
            elseif (present(filename)) then
                u = open_file('in', trim(filename))
                do i = 1, 3
                    read (u, *) epsbuf(:, i)
                end do
                do i = 1, map%n_atom_uc
                    do j = 1, 3
                        read (u, *) Zbuf(:, j, i)
                    end do
                end do
                close (u)
            else
                call lo_stop_gracefully(['You have to provide Born effective charges and dielectric tensors somehow.'], &
                                        lo_exitcode_param, __FILE__, __LINE__)
            end if

            ! For consistancy, I need both eps and the raw derivatives.
            ! Used for different things in different places.
            derivbuf = epsbuf
            do i = 1, 3
                derivbuf(i, i) = 1.0_r8 - derivbuf(i, i)
            end do
            derivbuf = derivbuf*p%volume*0.25_r8/lo_pi

        end block readZandeps

        ! Solve for the irreducible dielectric tensor
        solveeps: block
            real(r8), dimension(9, map%xuc%nx_eps_global) :: wA
            real(r8), dimension(9, 1) :: wB
            integer :: i

            ! Solve for the dielectric tensor
            wA = map%eps_global_shell%coeff
            wB(:, 1) = lo_flattentensor(epsbuf)
            call lo_dgelss(wA, wB, info=i)
            if (i .ne. 0) then
                call lo_stop_gracefully(['dgelss exit status '//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)
            end if
            map%xuc%x_eps_global = wB(1:map%xuc%nx_eps_global, 1)
            ! Solve for derivatives, exactly the same thing again
            wA = map%eps_global_shell%coeff
            wB(:, 1) = lo_flattentensor(derivbuf)
            call lo_dgelss(wA, wB, info=i)
            if (i .ne. 0) then
                call lo_stop_gracefully(['dgelss exit status '//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)
            end if
            map%xuc%x_eps_global_deriv = wB(1:map%xuc%nx_eps_global, 1)

            ! store the symmetrized thing in epsbuf
            wB(:, 1) = matmul(map%eps_global_shell%coeff, map%xuc%x_eps_global)
            epsbuf = lo_unflatten_2tensor(wB(:, 1))
            if (verbosity .gt. 0) write (lo_iou, *) '... got dielectric tensor'
        end block solveeps

        ! solve for the irreducible Z
        if (map%xuc%nx_Z_singlet .gt. 0) then
            solveZ: block
                type(lo_ewald_parameters) :: ew
                real(r8), dimension(:, :), allocatable :: zM, coeffM
                real(r8), dimension(:, :), allocatable :: constraintM
                real(r8), dimension(:), allocatable :: c, d
                integer :: a1, i, n_constraint

                call mem%allocate(coeffM, [9*map%n_atom_uc, map%xuc%nx_Z_singlet], &
                                  persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%allocate(zM, [map%n_atom_uc*9, 1], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                coeffM = 0.0_r8
                zM = 0.0_r8

                call lo_coeffmatrix_unitcell_Z_singlet(map, coeffM)
                do a1 = 1, map%n_atom_uc
                    zM((a1 - 1)*9 + 1:a1*9, 1) = lo_flattentensor(zbuf(:, :, a1))
                end do
                call lo_constraintmatrix_unitcell_Z_singlet(map, n_constraint, constraintM, mem)

                ! fkdev: disable rotational for now, they break too much
                call lo_dgelss(coeffM, zM, info=i)
                map%xuc%x_Z_singlet = zM(1:map%xuc%nx_Z_singlet, 1)

                ! if (n_constraint .gt. 0) then
                !     write (*, *) 'HELLO HELLO SOLVING FOR Z'
                !     allocate (c(map%xuc%nx_Z_singlet))
                !     allocate (d(map%xuc%nx_Z_singlet))
                !     c = ZM(:, 1)
                !     d = 0.0_r8
                !     call lo_dgglse(coeffM, constraintM, c, d, map%xuc%x_Z_singlet, info=i)
                ! else
                !     call lo_dgelss(coeffM, zM, info=i)
                !     map%xuc%x_Z_singlet = zM(1:map%xuc%nx_Z_singlet, 1)
                ! end if

                if (i .ne. 0) then
                    call lo_stop_gracefully(['dgelss exit status '//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)
                end if
                ! store results

                if (verbosity .gt. 0) write (lo_iou, *) '... got Born effective charges'

!write(*,*) 'constr err:',matmul(constraintM,matmul(coeffM,map%xuc%x_Z_singlet))

                ! Make sure it's hermitian? Seems like the sensible thing to do.
                call lo_coeffmatrix_unitcell_Z_singlet(map, coeffM)

                call ew%set(p, epsbuf, 2, 1E-20_r8, verbosity)
                call ew%force_borncharges_Hermitian(map%xuc%x_Z_singlet, coeffM, epsbuf, p, verbosity)

                ! cleanup
                call mem%deallocate(coeffM, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%deallocate(Zbuf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
                call mem%deallocate(zM, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            end block solveZ
        end if
    end if
    ! Spread the solution to everyone
    call mw%bcast(map%xuc%x_eps_global, from=solrnk)
    call mw%bcast(map%xuc%x_eps_global_deriv, from=solrnk)
    if (map%xuc%nx_Z_singlet .gt. 0) then
        call mw%bcast(map%xuc%x_Z_singlet, from=solrnk)
    end if
end subroutine

!> build the rotational constraint matrix
subroutine lo_constraintmatrix_unitcell_Z_singlet(map, n_constraint, constraint, mem)
    !> forcemap
    class(lo_forcemap), intent(inout) :: map
    !> number of constraints
    integer, intent(out) :: n_constraint
    !> constraint matrix
    real(r8), dimension(:, :), allocatable, intent(out) :: constraint
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:, :), allocatable :: coeffM, m0, m1, m2
    real(r8), dimension(18, 9) :: tcm
    integer :: i

    tcm = 0.0_r8
    ! tcm(1,2)=2.0_r8
    ! tcm(2,3)=2.0_r8
    ! tcm(3,6)=1.0_r8
    ! tcm(4,3)=1.0_r8
    ! tcm(5,8)=1.0_r8
    ! tcm(6,2)=1.0_r8
    ! tcm(7,6)=1.0_r8
    ! tcm(8,3)=1.0_r8
    ! tcm(9,4)=2.0_r8
    ! tcm(10,6)=2.0_r8
    ! tcm(11,7)=1.0_r8
    ! tcm(12,4)=1.0_r8
    ! tcm(13,8)=1.0_r8
    ! tcm(14,2)=1.0_r8
    ! tcm(15,7)=1.0_r8
    ! tcm(16,4)=1.0_r8
    ! tcm(17,7)=2.0_r8
    ! tcm(18,8)=2.0_r8
    tcm(1, 4) = 2.0_r8
    tcm(2, 7) = 2.0_r8
    tcm(3, 8) = 1.0_r8
    tcm(4, 7) = 1.0_r8
    tcm(5, 6) = 1.0_r8
    tcm(6, 4) = 1.0_r8
    tcm(7, 8) = 1.0_r8
    tcm(8, 7) = 1.0_r8
    tcm(9, 2) = 2.0_r8
    tcm(10, 8) = 2.0_r8
    tcm(11, 3) = 1.0_r8
    tcm(12, 2) = 1.0_r8
    tcm(13, 6) = 1.0_r8
    tcm(14, 4) = 1.0_r8
    tcm(15, 3) = 1.0_r8
    tcm(16, 2) = 1.0_r8
    tcm(17, 3) = 2.0_r8
    tcm(18, 6) = 2.0_r8

    ! Get the constraint matrix for the full thing
    call mem%allocate(m0, [map%n_atom_uc*18, map%n_atom_uc*9], &
                      persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(m1, [map%n_atom_uc*18, map%xuc%nx_Z_singlet], &
                      persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(coeffM, [map%n_atom_uc*9, map%xuc%nx_Z_singlet], &
                      persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    m0 = 0.0_r8
    m1 = 0.0_r8
    coeffM = 0.0_r8

    m0 = 0.0_r8
    do i = 1, map%n_atom_uc
        m0((i - 1)*18 + 1:i*18, (i - 1)*9 + 1:i*9) = tcm
    end do
    ! Get coefficient matrix
    call lo_coeffmatrix_unitcell_Z_singlet(map, coeffM)
    ! Get constraint matrix acting on the irreducible?
    m1 = matmul(m0, coeffM)
    ! Compress to linearly independent constraints?
    call lo_compress_equations(m1, n_constraint, constraint, .false.)

    call mem%deallocate(m0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(m1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(coeffM, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

end submodule
