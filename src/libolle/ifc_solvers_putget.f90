#include "precompilerdefinitions"
submodule (ifc_solvers) ifc_solvers_putget
    use type_blas_lapack_wrappers, only: lo_gemm,lo_dgelss
    use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
implicit none
contains

!> Get the second order irreducible representation from a supercell dynamical matrix
module subroutine lo_irreducible_forceconstant_from_supercell_dynmat( &
    map, ss, dynmat_M, lrdynmat_M, dynmat_T, lrdynmat_T, enforce_constraints, &
    fullhavemass, longrangehavemass, subset, mw, verbosity)
    !> forcemap, irreducible representation will be return in map%xuc%x_fc_pair
    type(lo_forcemap), intent(inout) :: map
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> dynamical matrix, with or without masses, in matrix form
    real(r8), dimension(:, :), intent(in), optional :: dynmat_M
    !> dynamical matrix, with or without masses, in tensor form
    real(r8), dimension(:, :, :, :), intent(in), optional :: dynmat_T
    !> long-range polar dynamical matrix, with or without masses, in matrix form
    real(r8), dimension(:, :), intent(in), optional :: lrdynmat_M
    !> long-range polar dynamical matrix, with or without masses, in tensor form
    real(r8), dimension(:, :, :, :), intent(in), optional :: lrdynmat_T
    !> Should I enforce constraints? Not sure about doing it always
    logical, intent(in), optional :: enforce_constraints
    !> Does the full dynamical matrix have masses multiplied in?
    logical, intent(in) :: fullhavemass
    !> Does the long-range dynamical matrix have masses multiplied in?
    logical, intent(in) :: longrangehavemass
    !> If on tensor form, do I have just a subset of the full tensor form?
    integer, dimension(:), intent(in), optional :: subset
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:, :, :, :), allocatable :: dmt
    real(r8) :: timer
    integer, dimension(:), allocatable :: ind
    integer :: fcna
    logical :: enforce

    timer = walltime()

    ! Now figure out exactly what to do
    init: block
        real(r8), dimension(:, :, :, :), allocatable :: lrdmt
        real(r8), dimension(3, 3) :: m0
        integer :: i, a1, a2

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'Converting supercell dynamical matrix to irreducible form'
        end if

        ! First some dummy tests
        if (present(dynmat_M) .and. present(dynmat_T)) then
            call lo_stop_gracefully(['Provide dynamical matrix on either matrix or tensor form, not both'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if
        if (present(dynmat_M) .and. present(dynmat_T)) then
            call lo_stop_gracefully(['Provide longrange dynamical matrix on either matrix or tensor form, not both'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if

        ! Which pairs do I have?
        if (present(subset)) then
            if (present(dynmat_m)) then
                call lo_stop_gracefully(['Subset only works with tensor for for now.'], &
                                        lo_exitcode_param, __FILE__, __LINE__)
            end if
            fcna = size(subset)
            if (fcna .lt. map%n_atom_uc) then
                call lo_stop_gracefully(['Subset too small to determine irreducible representation'], &
                                        lo_exitcode_param, __FILE__, __LINE__)
            end if
            lo_allocate(ind(fcna))
            ind = subset
        else
            fcna = ss%na
            lo_allocate(ind(fcna))
            do i = 1, ss%na
                ind(i) = i
            end do
        end if

        ! Should I constrain?
        if (present(enforce_constraints)) then
            enforce = enforce_constraints
        else
            enforce = .false.
        end if

        ! Some space for the matrices
        lo_allocate(dmt(3, 3, ss%na, fcna))
        lo_allocate(lrdmt(3, 3, ss%na, fcna))
        dmt = 0.0_r8
        lrdmt = 0.0_r8
        if (present(dynmat_M)) then
            ! Sort it into tensor form
            do a1 = 1, ss%na
            do a2 = 1, ss%na
                dmt(:, :, a1, a2) = dynmat_M((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3)
            end do
            end do
        elseif (present(dynmat_T)) then
            ! Just copy
            dmt = dynmat_T
        else
            call lo_stop_gracefully(['You have to provide a dynamical matrix.'], lo_exitcode_param, __FILE__, __LINE__)
        end if
        ! remove masses, if necessary
        if (fullhavemass) then
            do i = 1, ss%na
                a2 = ind(i)
                do a1 = 1, ss%na
                    dmt(:, :, a1, i) = dmt(:, :, a1, i)/(ss%invsqrtmass(a1)*ss%invsqrtmass(a2))
                end do
            end do
        end if

        ! Now the polar part
        if (map%polarcorrectiontype .eq. 3) then
            if (present(lrdynmat_M)) then
                ! Sort it into tensor form
                do a1 = 1, ss%na
                do a2 = 1, ss%na
                    lrdmt(:, :, a1, a2) = lrdynmat_M((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3)
                end do
                end do
            elseif (present(lrdynmat_T)) then
                ! Just copy
                lrdmt = lrdynmat_T
            else
                call lo_stop_gracefully(['You have to provide a dynamical matrix.'], lo_exitcode_param, __FILE__, __LINE__)
            end if
        end if
        ! remove masses, if necessary
        if (longrangehavemass) then
            do i = 1, ss%na
                a2 = ind(i)
                do a1 = 1, ss%na
                    lrdmt(:, :, a1, i) = lrdmt(:, :, a1, i)/(ss%invsqrtmass(a1)*ss%invsqrtmass(a2))
                end do
            end do
        end if

        ! Now I can subtract out the longrange
        dmt = dmt - lrdmt
        ! And enforce the translational invariance
        do i = 1, fcna
            a2 = ind(i)
            m0 = 0.0_r8
            do a1 = 1, ss%na
                if (a1 .ne. a2) m0 = m0 + dmt(:, :, a1, i)
            end do
            dmt(:, :, a2, i) = -m0
        end do
        ! Now I have arranged everything to a neat format.
        lo_deallocate(lrdmt)

        if (verbosity .gt. 0) then
            write (*, *) '... cleaned up the input (', tochar(walltime() - timer), 's)'
        end if
    end block init

    solvepart: block
        ! real(r8), dimension(9,9) :: bm0
        ! real(r8), dimension(9) :: fv,fw
        ! real(r8), dimension(3,3) :: m0
        ! integer, dimension(9) :: thetaind
        integer :: i, j, a1, a2, unsh, unop, ntheta, ctr
        !
        ! ! First count number of pairs I need to bother with
        ! ctr=0
        ! do j=1,fcna
        !     a1=ind(j)
        !     do i=1,map%ss(a1)%n_fc_pair
        !         unsh=map%ss(a1)%pairind(1,i)
        !         if ( map%fc_pair_shell(unsh)%nx .eq. 0 ) cycle
        !         ctr=ctr+1
        !     enddo
        ! enddo

        if (verbosity .gt. 0) call lo_progressbar_init()
        map%xuc%x_fc_pair = 0.0_r8
        do unsh = 1, map%n_fc_pair_shell
            ntheta = map%fc_pair_shell(unsh)%nx
            if (ntheta .eq. 0) cycle
            ! Make parallel?
            if (mod(unsh, mw%n) .ne. mw%r) cycle
            write (*, *) 'FIXME IRR FROM FULL'
            stop
        end do

        ! if ( verbosity .gt. 0 ) call lo_progressbar_init()
        ! map%xuc%x_fc_pair=0.0_r8
        ! do unsh=1,map%n_fc_pair_shell
        !     ntheta=map%fc_pair_shell(unsh)%nx
        !     if ( ntheta .eq. 0 ) cycle
        !     ! Make parallel?
        !     if ( mod(unsh,mw%n) .ne. mw%r ) cycle
        !     ! Locate the prototype for this shell
        !     locateloop: do j=1,fcna
        !         a1=ind(j)
        !         do i=1,map%ss(a1)%n_fc_pair
        !             if ( unsh .ne. map%ss(a1)%pairind(1,i) ) cycle
        !             unop=map%ss(a1)%pairind(2,i)
        !             if ( unop .eq. 1 ) then
        !                 a2=map%ss(a1)%pairind(3,i)
        !                 unop=map%ss(a1)%pairind(2,i)
        !                 thetaind(1:ntheta)=map%fc_pair_shell(unsh)%ind_global
        !                 m0=dmt(:,:,a2,j)
        !                 fv=lo_flattentensor( m0 )
        !                 bm0(:,1:ntheta)=matmul(map%op_pair(unop)%sotr,map%fc_pair_shell(unsh)%coeff)
        !                 call lo_linear_least_squares(bm0(:,1:ntheta),fv,fw(1:ntheta))
        !                 map%xuc%x_fc_pair(thetaind(1:ntheta))=fw(1:ntheta)
        !                 exit locateloop
        !             endif
        !         enddo
        !     enddo locateloop
        !     if ( verbosity .gt. 0 .and. unsh .lt. map%n_fc_pair_shell ) then
        !         call lo_progressbar(' ... grabbing irreducible',unsh,map%n_fc_pair_shell,walltime()-timer)
        !     endif
        ! enddo
        ! call mw%allreduce('sum',map%xuc%x_fc_pair)
        ! if ( verbosity .gt. 0 ) then
        !     call lo_progressbar(' ... grabbing irreducible',map%n_fc_pair_shell,map%n_fc_pair_shell,walltime()-timer)
        ! endif
    end block solvepart
end subroutine

!> Calculate the dynamical matrix from the irreducible representation
module subroutine lo_supercell_dynmat_from_irreducible_forceconstant(map, uc, ss, x, dynmat, mw, mem, dmlr, shortrange)
    !> forcemap
    type(lo_forcemap), intent(in) :: map
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> irreducible represenation
    real(r8), dimension(:), intent(in) :: x
    !> dynamical matrix
    real(r8), dimension(:, :), intent(out) :: dynmat
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> pre-calculated longrange polar forceconstant?
    real(r8), dimension(:, :, :, :), intent(in), optional :: dmlr
    !> Calculate only the short-ranged part
    logical, intent(in), optional :: shortrange

    real(r8), dimension(:, :, :, :), allocatable :: dm
    real(r8), dimension(:, :), allocatable :: shell_fc
    real(r8), dimension(3, 3) :: m0
    real(r8), dimension(9) :: fv
    real(r8) :: f0
    integer :: sh, op, a1, a2, i, ctr
    logical :: longrange

#ifdef AGRESSIVE_SANITY
    !@todo Add more sanity tests
    if (map%n_atom_ss .ne. ss%na) then
        call lo_stop_gracefully(['Incompatible supercell and forcemap'], lo_exitcode_param, __FILE__, __LINE__)
    end if
#endif

    if (present(shortrange)) then
        longrange = .not. shortrange
    else
        longrange = .true.
    end if

    ! Start by evaluating the forceconstants for each shell
    lo_allocate(shell_fc(9, map%n_fc_pair_shell))
    lo_allocate(dm(3, 3, ss%na, ss%na))
    dm = 0.0_r8
    shell_fc = 0.0_r8
    do sh = 1, map%n_fc_pair_shell
        if (map%fc_pair_shell(sh)%nx .eq. 0) cycle
        if (mod(sh, mw%n) .ne. mw%r) cycle
        shell_fc(:, sh) = matmul(map%fc_pair_shell(sh)%coeff, x(map%fc_pair_shell(sh)%ind_global))
    end do
    call mw%allreduce('sum', shell_fc)

    ! Then muliply it out everywhere
    dynmat = 0.0_r8
    ctr = 0
    do i = 1, map%xss%n_fc_pair
        if (mod(i, mw%n) .ne. mw%r) cycle
        sh = map%xss%ind_fc_pair(1, i)
        op = map%xss%ind_fc_pair(2, i)
        a1 = map%xss%ind_fc_pair(3, i)
        a2 = map%xss%ind_fc_pair(4, i)
        if (map%fc_pair_shell(sh)%nx .eq. 0) cycle
        fv = matmul(map%op_pair(op)%sotr, shell_fc(:, sh))
        m0 = lo_unflatten_2tensor(fv)
        dm(:, :, a1, a2) = dm(:, :, a1, a2) + m0
    end do
    call mw%allreduce('sum', dm)

    ! Now for the polar things. Have to think about where to do the heuristics. For now I assume it is input.
    if (map%polar .gt. 0) then
    if (longrange) then
        !@TODO think about non-analytical part
        if (present(dmlr)) then
            ! it was already precalculated
            dm = dm + dmlr
        else
            ! I have to calculate it
            longrangefc: block
                real(r8), dimension(:, :, :, :), allocatable :: ddm
                type(lo_forceconstant_secondorder) :: fc

                lo_allocate(ddm(3, 3, ss%na, ss%na))
                ddm = 0.0_r8
                if (mw%r .eq. mw%n - 1) then
                    call map%get_secondorder_forceconstant(uc, fc, mem, -1)
                    call fc%supercell_longrange_dynamical_matrix_at_gamma(ss, ddm, 1E-10_r8)
                end if
                call mw%bcast(ddm, mw%n - 1)
                dm = dm + ddm
                lo_deallocate(ddm)
            end block longrangefc
        end if
    end if
    end if

    ! Fix acoustic sum rule
    do a2 = 1, ss%na
        m0 = 0.0_r8
        do a1 = 1, ss%na
            if (a1 .ne. a2) m0 = m0 + dm(:, :, a1, a2)
        end do
        dm(:, :, a2, a2) = -m0
    end do

    ! And multiply in masses
    dynmat = 0.0_r8
    do a2 = 1, ss%na
        if (mod(a2, mw%n) .ne. mw%r) cycle
        do a1 = 1, ss%na
            f0 = ss%invsqrtmass(a1)*ss%invsqrtmass(a2)
            dynmat((a1 - 1)*3 + 1:a1*3, (a2 - 1)*3 + 1:a2*3) = dm(:, :, a1, a2)*f0
        end do
    end do
    call mw%allreduce('sum', dynmat)

    ! And cleanup
    lo_deallocate(shell_fc)
    lo_deallocate(dm)
end subroutine

!> Get the irreducible representation from dynamical matrices on a q-mesh
module subroutine lo_irreducible_forceconstant_from_qmesh_dynmat( &
    map, uc, qp, nq, dynmat, lrdynmat, fullhavemasses, longrangehavemasses, mw, verbosity, enforce, weights)
    !> forcemap for symmetry
    type(lo_forcemap), intent(inout) :: map
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> q-points
    class(lo_qpoint), dimension(:), intent(in) :: qp
    !> number of q-points on this rank
    integer, intent(in) :: nq
    !> dynamical matrices (full, long+shortrange)
    complex(r8), dimension(:, :, :), intent(in) :: dynmat
    !> long-ranged dynamical matrices (only the longrange)
    complex(r8), dimension(:, :, :), intent(in) :: lrdynmat
    !> does the full dynamical matrix have masses?
    logical, intent(in) :: fullhavemasses
    !> does the longrange dynamical matrix have masses?
    logical, intent(in) :: longrangehavemasses
    !> mpi-helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> enforce symmetries?
    logical, intent(in), optional :: enforce
    !> add weights
    real(r8), dimension(:), intent(in), optional :: weights

    complex(r8), dimension(:, :, :), allocatable :: wdyn
    real(r8) :: timer
    integer :: nx, nd, solrnk
    logical :: fixsym

    init: block
        real(r8), dimension(:, :, :), allocatable :: selfterm
        real(r8) :: f0, f1
        integer :: q, a1, a2, i1, i2, j1, j2, gammaind

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'GETTING IRREDUCIBLE REPRESENTATION FROM DYNAMICAL MATRICES'
            timer = walltime()
        end if

        ! Some shorthand
        nx = map%xuc%nx_fc_pair
        nd = (map%n_atom_uc*3)**2
        solrnk = mw%n - 1
        ! Enforce all symmetries
        if (present(enforce)) then
            fixsym = enforce
        else
            fixsym = .true.
        end if
        ! Make a copy of the dynamical matrices, and massage a little
        lo_allocate(selfterm(3, 3, map%n_atom_uc))
        selfterm = 0.0_r8
        gammaind = -1
        if (nq .gt. 0) then
            ! Remove longrange part, and take away the masses, if necessary
            lo_allocate(wdyn(map%n_atom_uc*3, map%n_atom_uc*3, nq))
            wdyn = 0.0_r8
            do q = 1, nq
                do a1 = 1, map%n_atom_uc
                do a2 = 1, map%n_atom_uc
                    if (fullhavemasses) then
                        f0 = sqrt(uc%mass(a1)*uc%mass(a2))
                    else
                        f0 = 1.0_r8
                    end if
                    if (longrangehavemasses) then
                        f1 = sqrt(uc%mass(a1)*uc%mass(a2))
                    else
                        f1 = 1.0_r8
                    end if
                    i1 = (a1 - 1)*3 + 1
                    i2 = (a1*3)
                    j1 = (a2 - 1)*3 + 1
                    j2 = (a2*3)
                    wdyn(j1:j2, i1:i2, q) = dynmat(j1:j2, i1:i2, q)*f0 - lrdynmat(j1:j2, i1:i2, q)*f1
                    !do i=1,3
                    !do j=1,3
                    !    ii=(a1-1)*3+i
                    !    jj=(a2-1)*3+j
                    !    wdyn(jj,ii,q)=dynmat(jj,ii,q)*f0-lrdynmat(jj,ii,q)*f1
                    !enddo
                    !enddo
                end do
                end do
            end do
            ! Locate gamma, and calculate the ASR-correction
            do q = 1, nq
                if (lo_sqnorm(qp(q)%r) .gt. lo_sqtol) cycle
                ! Now I found Gamma.
                gammaind = q
                !f0=0.0_r8
                do a1 = 1, map%n_atom_uc
                do a2 = 1, map%n_atom_uc
                    i1 = (a1 - 1)*3 + 1
                    i2 = (a1*3)
                    j1 = (a2 - 1)*3 + 1
                    j2 = (a2*3)
                    selfterm(:, :, a1) = selfterm(:, :, a1) + real(wdyn(j1:j2, i1:i2, q))
                    !f0=f0+sum(abs(aimag(wdyn(j1:j2,i1:i2,q))))
                end do
                end do
                ! if f0 is too large, throw an error? Maybe
                exit
            end do
        end if

        ! Sanity check to make sure I actually found gamma
        call mw%allreduce('max', gammaind, i1)

        if (i1 .lt. 0 .and. mw%talk) then
            write (*, *) 'WARNING: could not locate Gamma in the set of q-points. Will continue'
            write (*, *) 'but this is most likely bad input. And the results probably nonsense.'
        end if
        ! Get the self-term everywhere
        call mw%allreduce('sum', selfterm)

        ! Sanity check
        f0 = sum(abs(selfterm))
        ! And adjust with the ASR correction
        if (nq .gt. 0) then
            do q = 1, nq
                do a1 = 1, map%n_atom_uc
                    i1 = (a1 - 1)*3 + 1
                    i2 = (a1*3)
                    wdyn(i1:i2, i1:i2, q) = wdyn(i1:i2, i1:i2, q) - selfterm(:, :, a1)
                end do
            end do
        end if
        lo_deallocate(selfterm)

        ! Get the total number of q-points
        i1 = nq
        i2 = 0
        call mw%allreduce('sum', i1, i2)

        if (verbosity .gt. 0) then
            write (*, *) '        number of q-points: ', tochar(i2)
            write (*, *) '    enforce all symmetries:', fixsym
        end if

    end block init

    solve: block
        complex(r8), dimension(:, :), allocatable :: cfm
        complex(r8), dimension(:), allocatable :: cdm
        real(r8), dimension(:, :), allocatable :: ATA, rcfm
        real(r8), dimension(:), allocatable :: ATB, rdm
        integer :: q, a1, a2, i, j, k, ii, jj

        lo_allocate(ATA(nx, nx))
        lo_allocate(ATB(nx))
        lo_allocate(cfm(nd, nx))
        lo_allocate(rcfm(2*nd, nx))
        lo_allocate(rdm(2*nd))
        lo_allocate(cdm(nd))
        ATA = 0.0_r8
        ATB = 0.0_r8
        cfm = 0.0_r8
        rcfm = 0.0_r8
        rdm = 0.0_r8
        cdm = 0.0_r8

        if (verbosity .gt. 0) call lo_progressbar_init()
        do q = 1, nq
            if (mod(q, mw%n) .ne. mw%r) cycle
            ! grab coefficient matrix
            call lo_dynamical_matrix_coefficient_matrix_for_single_q(map, qp(q)%r, cfm)
            ! grab flattened dynamical matrix. I hate indices.
            k = 0
            do a1 = 1, map%n_atom_uc
            do a2 = 1, map%n_atom_uc
                do i = 1, 3
                do j = 1, 3
                    ii = (a1 - 1)*3 + i
                    jj = (a2 - 1)*3 + j
                    !k=(ii-1)*map%n_atom_uc*3+jj
                    k = k + 1
                    cdm(k) = wdyn(jj, ii, q)
                end do
                end do
            end do
            end do
            if (present(weights)) then
                cdm = cdm*weights(q)
                cfm = cfm*weights(q)
            end if

            ! Weird sanity test if I get bad numbers?
            do i = 1, nd
                if (abs(cdm(i)) .gt. 1E20_r8) then
                    cdm(i) = 0.0_r8
                    cfm(i, :) = 0.0_r8
                end if
            end do

            ! Put it in real matrices
            rcfm(1:nd, :) = real(cfm)
            rcfm(nd + 1:2*nd, :) = aimag(cfm)
            rdm(1:nd) = real(cdm)
            rdm(nd + 1:2*nd) = aimag(cdm)
            ! Accumulate
            call lo_gemm(rcfm, rcfm, ATA, transa='T', beta=1.0_r8)
            call lo_gemv(rcfm, rdm, ATB, trans='T', beta=1.0_r8)
            if (verbosity .gt. 0 .and. q .lt. nq) call lo_progressbar(' ... generating coefficients', q, nq, walltime() - timer)
        end do
        ! Now accumulate across ranks, eventually
        call mw%allreduce('sum', ATA)
        call mw%allreduce('sum', ATB)

        if (verbosity .gt. 0) call lo_progressbar(' ... generating coefficients', nq, nq, walltime() - timer)

        ! And solve
        if (mw%r .eq. solrnk) then
            if (fixsym) then
                call lo_linear_least_squares( &
                    ATA, ATB, map%xuc%x_fc_pair, map%constraints%eq2, map%constraints%d2, map%constraints%neq2, gramified=.true.)
            else
                call lo_linear_least_squares( &
                    ATA, ATB, map%xuc%x_fc_pair, map%constraints%eq2, map%constraints%d2,  0, gramified=.true.)
            end if
        end if
        call mw%bcast(map%xuc%x_fc_pair, solrnk)

        lo_deallocate(ATA)
        lo_deallocate(ATB)
        lo_deallocate(cfm)
        lo_deallocate(rcfm)
        lo_deallocate(rdm)
        lo_deallocate(cdm)

        if (verbosity .gt. 0) write (*, *) 'solved for irreducible forceconstants (', tochar(walltime() - timer), 's)'

    end block solve
    ! And final cleanup
    lo_deallocate(wdyn)
end subroutine

!> Construct the dynamical matrix coefficient matrix for a specific q-point
subroutine lo_dynamical_matrix_coefficient_matrix_for_single_q(map, qv, coefficientmatrix, uc)
    !> forcemap
    type(lo_forcemap), intent(in) :: map
    !> q-vector
    real(r8), dimension(3), intent(in) :: qv
    !> coefficient matrix
    complex(r8), dimension(:, :), intent(out) :: coefficientmatrix
    !> unitcell, in case I want the masses in there
    type(lo_crystalstructure), intent(in), optional :: uc

    complex(r8), dimension(:, :, :, :, :), allocatable :: Ck
    complex(r8), dimension(9, 9) :: C1
    complex(r8) :: expiqr
    real(r8) :: k_dot_r
    integer :: i, j, k, l, ii, jj, sh, o, a1, a2, ipair
    integer :: nx, na, nfc

    ! Size of things
    na = map%n_atom_uc         ! number of atoms
    nx = map%xuc%nx_fc_pair    ! dimensions of irreducible IFC

    if (size(coefficientmatrix, 1) .ne. 3*3*na*na ) then
        write (*, *) 'bad dimensions in dynmatrixcoeffM'
        stop
    end if
    if (size(coefficientmatrix, 2) .ne. nx) then
        write (*, *) 'bad dimensions in dynmatrixcoeffM'
        stop
    end if

    allocate (Ck(3, 3, na, na, nx))
    Ck = 0.0_r8
    do ipair = 1, map%xuc%n_fc_pair
        !write(*,*) 'DANGER INDEX CHECK THIS',__LINE__,__FILE__
        a1 = map%xuc%fc_pair(ipair)%i1
        a2 = map%xuc%fc_pair(ipair)%i2
        sh = map%xuc%fc_pair(ipair)%irreducible_shell
        o = map%xuc%fc_pair(ipair)%operation_from_shell
        nfc = map%fc_pair_shell(sh)%nx
        if (nfc .eq. 0) cycle
        ! Get the Fourier transform thingy
        k_dot_r = dot_product(map%xuc%fc_pair(ipair)%lv, qv)*lo_twopi
        expiqr = cmplx(cos(k_dot_r), sin(k_dot_r), r8)
        C1 = 0.0_r8
        C1(:, 1:nfc) = matmul(map%op_pair(o)%sotr, map%fc_pair_shell(sh)%coeff)
        ! Add to the self-term
        do k = 1, nfc
            l = map%fc_pair_shell(sh)%ind_global(k)
            do i = 1, 3
            do j = 1, 3
                ii = (i - 1)*3 + j
                Ck(i, j, a1, a1, l) = Ck(i, j, a1, a1, l) - C1(ii, k)
            end do
            end do
        end do
        ! not the self-term
        C1 = C1*expiqr
        do k = 1, nfc
            l = map%fc_pair_shell(sh)%ind_global(k)
            do i = 1, 3
            do j = 1, 3
                ii = (i - 1)*3 + j
                Ck(i, j, a1, a2, l) = Ck(i, j, a1, a2, l) + C1(ii, k)
            end do
            end do
        end do
    end do

    ! do a1=1,map%n_atom_uc
    ! do p=1,map%uc(a1)%n_fc_pair
    !     sh=map%uc(a1)%fc_pair(p)%irreducible_shell
    !     o=map%uc(a1)%fc_pair(p)%operation_from_shell
    !     nfc=map%fc_pair_shell(sh)%nx
    !     if ( nfc .eq. 0 ) cycle
    !     ! Get the Fourier transform thingy
    !     a2=map%uc(a1)%fc_pair(p)%i2
    !     k_dot_r=dot_product( map%uc(a1)%fc_pair(p)%lv, qv )*lo_twopi
    !     expiqr=cmplx( cos(k_dot_r) , sin(k_dot_r) ,r8)
    !     C1=0.0_r8
    !     C1(:,1:nfc)=matmul(map%pairop(o)%sotr,map%fc_pair_shell(sh)%coeff)
    !     ! Add to the self-term
    !     do k=1,nfc
    !         l=map%fc_pair_shell(sh)%ind_global(k)
    !         do i=1,3
    !         do j=1,3
    !             ii=(i-1)*3+j
    !             Ck(i,j,a1,a1,l)=Ck(i,j,a1,a1,l)-C1(ii,k)
    !         enddo
    !         enddo
    !     enddo
    !     ! not the self-term
    !     C1=C1*expiqr
    !     do k=1,nfc
    !         l=map%fc_pair_shell(sh)%ind_global(k)
    !         do i=1,3
    !         do j=1,3
    !             ii=(i-1)*3+j
    !             Ck(i,j,a1,a2,l)=Ck(i,j,a1,a2,l)+C1(ii,k)
    !         enddo
    !         enddo
    !     enddo
    ! enddo
    ! enddo

    ! Put this in the right place
    coefficientmatrix = 0.0_r8
    jj = 0
    do a1 = 1, na
    do a2 = 1, na
        do i = 1, 3
        do j = 1, 3
            jj = jj + 1
            if (present(uc)) then
                coefficientmatrix(jj, :) = lo_chop(Ck(i, j, a1, a2, :), lo_sqtol)*uc%invsqrtmass(a1)*uc%invsqrtmass(a2)
            else
                coefficientmatrix(jj, :) = lo_chop(Ck(i, j, a1, a2, :), lo_sqtol)
            end if
        end do
        end do
    end do
    end do

    deallocate (Ck)
end subroutine

! !> Create a completely fake forceconstant instead, from a maximum frequency
! subroutine lo_create_fake_forceconstant(map,uc,ss,mw,verbosity,max_frequency,max_displacement,debye_temperature,longrange_dynmat,fakefc)
!     !> forcemap (with the irreducible representation returned in map%xuc%x_fc_pair)
!     type(lo_forcemap), intent(inout) :: map
!     !> unitcell
!     type(lo_crystalstructure), intent(in) :: uc
!     !> supercell
!     type(lo_crystalstructure), intent(in) :: ss
!     !> mpi helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> talk a lot?
!     integer, intent(in) :: verbosity
!     !> desired max frequency?
!     real(r8), intent(in), optional :: max_frequency
!     !> desired max thermal displacement?
!     real(r8), intent(in), optional :: max_displacement
!     !> desired debye temperature?
!     real(r8), intent(in), optional :: debye_temperature
!     !> pre-calculated longrange forceconstant?
!     real(r8), dimension(:,:,:,:), intent(in), optional :: longrange_dynmat
!     !> Return normal forceconstant
!     type(lo_forceconstant_secondorder), intent(out), optional :: fakefc
!
!     !type(lo_forceconstant_secondorder) :: fc,fcss
!     !integer, parameter :: expalpha=6
!
! !    real(r8), dimension(:,:), allocatable :: dynmat
! !    real(r8), dimension(:), allocatable :: eigenvalues,x
! !    integer, dimension(:,:,:), allocatable :: ssind
! !    integer :: nb
! !    real(r8) :: f0,f1,f2
! !
! !    ! Set some defaults
! !    init: block
! !        integer :: a1,l
! !
! !        if ( verbosity .gt. 0 ) then
! !            write(*,*) '... creating fake forceconstant'
! !        endif
! !
! !        ! Get the mapping thingy from unitcell to supercell
! !        map%xuc%x_fc_pair=0.0_r8
! !        call map%get_secondorder_forceconstant(uc,fc)
! !        l=0
! !        do a1=1,fc%na
! !            l=max(l,fc%atom(a1)%n)
! !        enddo
! !        lo_allocate(ssind(3,l,ss%na))
! !        ssind=0
! !        call fc%remap(uc,ss,fcss,ind=ssind)
! !
! !        ! Temporary space
! !        nb=ss%na*3
! !        lo_allocate(dynmat(nb,nb))
! !        lo_allocate(eigenvalues(nb))
! !        lo_allocate(x(map%nx_fc_pair))
! !        dynmat=0.0_r8
! !        eigenvalues=0.0_r8
! !        x=0.0_r8
! !    end block init
! !
! !    ! Optimize?
! !    opt: block
! !        real(r8) :: alpha,minom,maxom,f0,okalpha
! !        integer :: iter,nok
! !
! !        ! Start with a rough guess
! !        alpha=1000.0_r8
! !        call set_dynmat(dynmat,fc,ss,alpha,expalpha,ssind,dmlr)
! !        call lo_real_symmetric_eigenvalues_eigenvectors(dynmat,eigenvalues,careful=.false.,nzeros=3)
! !        minom=lo_negsqrt(minval(eigenvalues))
! !        maxom=lo_negsqrt(maxval(eigenvalues))
! !        ! Guess new alpha
! !        alpha=alpha*(max_frequency/maxom)**2
! !        ! Now adjust accordingly
! !        nok=0
! !        okalpha=0
! !        do iter=1,1000
! !            call set_dynmat(dynmat,fc,ss,alpha,expalpha,ssind,dmlr)
! !            call lo_real_symmetric_eigenvalues_eigenvectors(dynmat,eigenvalues,careful=.false.,nzeros=3)
! !            minom=lo_negsqrt(minval(eigenvalues))
! !            maxom=lo_negsqrt(maxval(eigenvalues))
! !            if ( minom .gt. -lo_tiny ) then
! !                nok=nok+1
! !                okalpha=alpha
! !                if ( abs(max_frequency-maxom) .lt. lo_freqtol ) exit
! !            endif
! !            ! We might be done here.
! !            if ( minom .lt. -lo_tiny ) then
! !                alpha=alpha*1.1
! !            elseif ( maxom .gt. max_frequency ) then
! !                alpha=alpha*0.95
! !            endif
! !            if ( nok .ge. 4 ) exit
! !            if ( nok .lt. 1 ) cycle
! !        enddo
! !
! !        ! Now we have a decent alpha. Return what we wanted
! !        call set_dynmat(dynmat,fc,ss,okalpha,expalpha,ssind,dmlr)
! !        call lo_irreducible_forceconstant_from_supercell_dynmat(map,ss,dynmat,dmlr,enforce=.false.)
! !        if ( present(irrifc) ) then
! !            irrifc=map%xuc%x_fc_pair
! !        endif
! !        if ( present(fakefc) ) then
! !            call map%get_secondorder_forceconstant(uc,fakefc)
! !        endif
! !    end block opt
! !
! !    contains
! !    ! Set the dynamical matrix
! !    subroutine set_dynmat(dynmat,fc,ss,x,n,ind,dmlr)
! !        real(r8), dimension(:,:), intent(inout) :: dynmat
! !        type(lo_forceconstant_secondorder), intent(in) :: fc
! !        type(lo_crystalstructure), intent(in) :: ss
! !        real(r8), intent(in) :: x
! !        integer, intent(in) :: n
! !        integer, dimension(:,:,:), intent(in) :: ind
! !        real(r8), dimension(:,:,:,:), intent(in), optional :: dmlr
! !
! !        real(r8), dimension(3,3) :: m0,m1
! !        real(r8) :: f0
! !        integer :: a1,a2,i,uca
! !
! !        dynmat=0.0_r8
! !        do a1=1,ss%na
! !            ! First shortrange
! !            uca=ind(1,1,a1)
! !            m1=0.0_r8
! !            do i=1,ind(2,1,a1)
! !                a2=ind(3,i,a1)
! !                f0=ss%invsqrtmass(a1)*ss%invsqrtmass(a2)
! !                m0=m_from_r(fc%atom(uca)%pair(i)%r,x,n)
! !                m1=m1+m0
! !                dynmat( (a1-1)*3+1:a1*3,(a2-1)*3+1:a2*3 )=dynmat((a1-1)*3+1:a1*3,(a2-1)*3+1:a2*3 )+m0*f0
! !            enddo
! !            dynmat( (a1-1)*3+1:a1*3,(a1-1)*3+1:a1*3 )=-m1*ss%invsqrtmass(a1)**2
! !        enddo
! !        if ( present(dmlr) ) then
! !            ! Longrange part
! !            do a1=1,ss%na
! !            do a2=1,ss%na
! !                f0=ss%invsqrtmass(a1)*ss%invsqrtmass(a2)
! !                dynmat( (a1-1)*3+1:a1*3,(a2-1)*3+1:a2*3 )=dynmat((a1-1)*3+1:a1*3,(a2-1)*3+1:a2*3 )+dmlr(:,:,a1,a2)*f0
! !            enddo
! !            enddo
! !        endif
! !    end subroutine
! !    ! Return a forceconstant component
! !    function m_from_r(v,x,n) result(m)
! !        real(r8), intent(in), dimension(3) :: v
! !        real(r8), intent(in) :: x
! !        integer, intent(in) :: n
! !        real(r8), dimension(3,3) :: m
! !
! !        real(r8) :: r,rx,ry,rz,fpp
! !
! !        r=norm2(v)
! !        if ( r .lt. lo_tol ) then
! !            m=0.0_r8
! !        else
! !            rx=v(1)
! !            ry=v(2)
! !            rz=v(3)
! !            fpp=x/(r**n)
! !            m(1,1) = -rx**2*r*fpp
! !            m(2,2) = -ry**2*r*fpp
! !            m(3,3) = -rz**2*r*fpp
! !            m(1,2) = rx*ry*( -r*fpp )
! !            m(1,3) = rx*rz*( -r*fpp )
! !            m(2,3) = ry*rz*( -r*fpp )
! !            m(2,1) = m(1,2)
! !            m(3,1) = m(1,3)
! !            m(3,2) = m(2,3)
! !            m=m/(r**2)
! !            m=lo_chop(m,1E-10_r8)
! !        endif
! !    end function
! end subroutine

end submodule
