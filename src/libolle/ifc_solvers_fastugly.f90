submodule(ifc_solvers) ifc_solvers_fastugly
!!
!! Solve for everything but without any checks or statistics.
!!
use type_blas_lapack_wrappers, only: lo_gemm
use type_forcemap, only: lo_coeffmatrix_singlet, lo_coeffmatrix_pair, lo_coeffmatrix_triplet, &
                         lo_coeffmatrix_quartet, lo_coeffmatrix_supercell_Z_singlet, &
                         lo_coeffmatrix_Z_pair, lo_coeffmatrix_Z_triplet, &
                         lo_coeffmatrix_unitcell_Z_singlet, lo_coeffmatrix_eps_singlet, lo_coeffmatrix_eps_pair
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder

implicit none

contains

!> normal least squares solution but with no checks or anything
module subroutine lo_solve_for_irreducible_ifc_fastugly(map, uc, ss, sim, mw, mem, verbosity, fix_secondorder)
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> simulation, preferrably with everything attached
    type(lo_mdsim), intent(in) :: sim
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> keep the second order fix
    logical, intent(in) :: fix_secondorder

    real(r8) :: timer, t0, t1
    integer :: nstep, solrnk

    timer = walltime()
    t0 = timer
    t1 = timer

    init: block
        integer :: i

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'FAST SOLUTION FOR IFCs'
        end if

        !@TODO add sanity checks that everything is included in the simulation

        ! Count steps per rank
        nstep = 0
        do i = 1, sim%nt
            if (mod(i, mw%n) .ne. mw%r) cycle
            nstep = nstep + 1
        end do

        ! Which rank to solve on
        solrnk = mw%n - 1
    end block init

    ! Polar things first?
    if ( fix_secondorder .eqv. .false. ) then
        select case (map%polar)
        case (0)
            ! Not polar, don't have to care
        case (1)
            ! A little polar, deal with it swiftly
            call lo_solve_for_borncharges(map, uc, filename='infile.lotosplitting', mw=mw, mem=mem, verbosity=verbosity)
        case (2)
            ! Very polar, do lots of things
            diel: block

            end block diel
        end select
    endif

    !> do all the force constant solving things
    ifc: block
        type(lo_forceconstant_secondorder) :: fc
        logical, dimension(81) :: rel_quartet_ntheta
        real(r8), dimension(:,:,:,:), allocatable :: polar_fc
        real(r8), dimension(:, :), allocatable :: buf_cm, buf_c, ATA, polar_f
        real(r8), dimension(:), allocatable :: buf_f, ATB
        real(r8), dimension(3) :: v0
        integer :: it, jt, iatom,jatom, ix, i, nx, ii, jj

        ! Build a force buffer
        if (nstep .gt. 0) then
            call mem%allocate(buf_f, nstep*sim%na*3, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
        else
            call mem%allocate(buf_f, 1, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
        end if
        buf_f = 0.0_r8
        jt = 0
        do it = 1, sim%nt
            if (mod(it, mw%n) .ne. mw%r) cycle
            jt = jt + 1
            do iatom = 1, sim%na
            do ix = 1, 3
                i = (jt - 1)*sim%na*3 + (iatom - 1)*3 + ix
                buf_f(i) = sim%f(ix, iatom, it)
            end do
            end do
        end do

        ! If polar, remove polar forces here.
        call map%get_secondorder_forceconstant(uc, fc, mem, verbosity=-1)
        if (fc%polar) then
            allocate (polar_fc(3, 3, ss%na, ss%na))
            allocate(polar_f(3,ss%na))
            polar_fc = 0.0_r8
            polar_f=0.0_r8
            call fc%supercell_longrange_dynamical_matrix_at_gamma(ss, polar_fc, 1E-15_r8)

            jt=0
            do it = 1, sim%nt
                if (mod(it, mw%n) .ne. mw%r) cycle
                ! Calculate the polar forces
                polar_f=0.0_r8
                do iatom=1,ss%na
                    do jatom=1,ss%na
                        polar_f(:,iatom)=polar_f(:,iatom) - matmul(polar_fc(:,:,iatom,jatom),sim%u(:,jatom,it))
                    enddo
                enddo

                ! Remove polar forces
                jt = jt + 1
                do iatom = 1, sim%na
                do ix = 1, 3
                    i = (jt - 1)*sim%na*3 + (iatom - 1)*3 + ix
                    buf_f(i) = buf_f(i) - polar_f(ix,iatom)
                end do
                end do
            end do
        end if

        ! Fix pairs first
        if (map%have_fc_pair) then
            nx = map%xuc%nx_fc_pair
            call mem%allocate(buf_c, [sim%na*3, nx], persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_cm, [nstep*sim%na*3, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(ATA, [nx, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(ATB, nx, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            buf_c = 0.0_r8
            buf_cm = 0.0_r8
            ATA = 0.0_r8
            ATB = 0.0_r8
            ! Build coefficient matrices and solve
            jt = 0
            do it = 1, sim%nt
                if (mod(it, mw%n) .ne. mw%r) cycle
                jt = jt + 1
                call lo_coeffmatrix_pair(sim%u(:, :, it), buf_c, map)
                ii = (jt - 1)*sim%na*3 + 1
                jj = jt*sim%na*3
                buf_cm(ii:jj, :) = buf_c
            end do

            if ( fix_secondorder .eqv. .false. ) then
                call lo_gemm(buf_cm, buf_cm, ATA, transa='T')
                call lo_gemv(buf_cm, buf_f, ATB, alpha=-1.0_r8, trans='T')
                call mw%allreduce('sum', ATA)
                call mw%allreduce('sum', ATB)
                if (mw%r .eq. solrnk) then
                    call lo_linear_least_squares( &
                        ATA, ATB, map%xuc%x_fc_pair, map%constraints%eq2, map%constraints%d2, map%constraints%neq2, gramified=.true.)
                end if
                call mw%bcast(map%xuc%x_fc_pair, solrnk, __FILE__, __LINE__)
            endif

            ! Subtract forces
            call lo_gemv(buf_cm, map%xuc%x_fc_pair, buf_f, alpha=1.0_r8, beta=1.0_r8)
            !buf_f = buf_f + matmul(buf_cm,map%xuc%x_fc_pair)
            call mem%deallocate(buf_c, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ATA, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ATB, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        if (map%have_fc_triplet) then
            nx = map%xuc%nx_fc_triplet
            call mem%allocate(buf_c, [sim%na*3, nx], persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_cm, [nstep*sim%na*3, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(ATA, [nx, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(ATB, nx, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            buf_c = 0.0_r8
            buf_cm = 0.0_r8
            ATA = 0.0_r8
            ATB = 0.0_r8
            ! Build coefficient matrices and solve
            jt = 0
            do it = 1, sim%nt
                if (mod(it, mw%n) .ne. mw%r) cycle
                jt = jt + 1
                call lo_coeffmatrix_triplet(sim%u(:, :, it), buf_c, map)
                ii = (jt - 1)*sim%na*3 + 1
                jj = jt*sim%na*3
                buf_cm(ii:jj, :) = buf_c
            end do
            call lo_gemm(buf_cm, buf_cm, ATA, transa='T')
            call lo_gemv(buf_cm, buf_f, ATB, alpha=-1.0_r8, trans='T')
            call mw%allreduce('sum', ATA)
            call mw%allreduce('sum', ATB)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares( &
                    ATA, ATB, map%xuc%x_fc_triplet, map%constraints%eq3, map%constraints%d3, map%constraints%neq3, gramified=.true.)
            end if
            call mw%bcast(map%xuc%x_fc_triplet, solrnk, __FILE__, __LINE__)
            ! Subtract forces
            call lo_gemv(buf_cm, map%xuc%x_fc_triplet, buf_f, alpha=1.0_r8, beta=1.0_r8)
            call mem%deallocate(buf_c, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ATA, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ATB, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        if (map%have_fc_quartet) then
            nx = map%xuc%nx_fc_quartet
            call mem%allocate(buf_c, [sim%na*3, nx], persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%allocate(buf_cm, [nstep*sim%na*3, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(ATA, [nx, nx], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(ATB, nx, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            buf_c = 0.0_r8
            buf_cm = 0.0_r8
            ATA = 0.0_r8
            ATB = 0.0_r8

            rel_quartet_ntheta = .false.
            do i = 1, map%n_fc_quartet_shell
                if (map%fc_quartet_shell(i)%nx .eq. 0) cycle
                rel_quartet_ntheta(map%fc_quartet_shell(i)%nx) = .true.
            end do

            ! Build coefficient matrices and solve
            jt = 0
            do it = 1, sim%nt
                if (mod(it, mw%n) .ne. mw%r) cycle
                jt = jt + 1
                call lo_coeffmatrix_quartet(sim%u(:, :, it), buf_c, map, rel_quartet_ntheta)
                ii = (jt - 1)*sim%na*3 + 1
                jj = jt*sim%na*3
                buf_cm(ii:jj, :) = buf_c
            end do
            call lo_gemm(buf_cm, buf_cm, ATA, transa='T')
            call lo_gemv(buf_cm, buf_f, ATB, alpha=-1.0_r8, trans='T')
            call mw%allreduce('sum', ATA)
            call mw%allreduce('sum', ATB)
            if (mw%r .eq. solrnk) then
                call lo_linear_least_squares( &
                    ATA, ATB, map%xuc%x_fc_quartet, map%constraints%eq4, map%constraints%d4, map%constraints%neq4, gramified=.true.)
            end if
            call mw%bcast(map%xuc%x_fc_quartet, solrnk, __FILE__, __LINE__)
            ! Subtract forces
            call lo_gemv(buf_cm, map%xuc%x_fc_quartet, buf_f, alpha=1.0_r8, beta=1.0_r8)
            call mem%deallocate(buf_c, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ATA, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ATB, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if

        call mem%deallocate(buf_f, persistent=.false., scalable=.true., file=__FILE__, line=__LINE__)
    end block ifc

end subroutine


end submodule
