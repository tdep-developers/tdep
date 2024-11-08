#include "precompilerdefinitions"
submodule(type_forcemap) type_forcemap_constraints
implicit none
contains

!> generate the linear constraints on the equations
module subroutine forceconstant_constraints(map, uc, rotational, huanginvariances, hermitian, hermitian_rhs, huang_rhs, rotational_rhs, verbosity)
    !> symmetry stuff rearranged into coordination shells
    class(lo_forcemap), intent(inout) :: map
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> build the rotational constraints
    logical, intent(in) :: rotational
    !> build the huang invariances
    logical, intent(in) :: huanginvariances
    !> build the Hermitian constraints
    logical, intent(in) :: hermitian
    !> rhs for the Hermitian constraints
    real(r8), dimension(:), intent(in) :: hermitian_rhs
    !> rhs for Huang invariances
    real(r8), dimension(:), intent(in) :: huang_rhs
    !> rhs for rotational invariances
    real(r8), dimension(:), intent(in) :: rotational_rhs
    !> how much to talk?
    integer, intent(in) :: verbosity

    real(r8) :: timer
    if (verbosity .gt. 0) then
        timer = walltime()
        write (*, *) ''
        write (*, *) 'GENERATING CONSTRAINTS'
    end if

    ! The clean constraints, that only involves a single order by itself
    cleanconstr: block
        real(r8) :: t0
        !if ( map%have_fc_singlet ) then
        !    t0=walltime()
        !    call lo_firstorder_rotational_translational(map,uc,map%constraints%eq1,map%constraints%neq1,rotational)
        !    if ( verbosity .gt. 0 ) then
        !        write(*,*) '... ',tochar(map%constraints%neq1),' firstorder constraints (',tochar(walltime()-t0),'s)'
        !    endif
        !else
        map%constraints%neq1 = 0
        !endif

        if (map%have_fc_pair) then
            t0 = walltime()
            call lo_secondorder_rot_herm_huang(map, uc, map%constraints%eq2, map%constraints%d2, map%constraints%neq2, rotational, huanginvariances, hermitian, hermitian_rhs, huang_rhs, rotational_rhs)
            if (verbosity .gt. 0) then
                write (*, *) '... ', tochar(map%constraints%neq2), ' secondorder constraints (', tochar(walltime() - t0), 's)'
            end if
        else
            map%constraints%neq2 = 0
        end if

        if (map%have_fc_triplet) then
            t0 = walltime()
            call lo_thirdorder_asr(map, map%constraints%eq3, map%constraints%neq3)
            if ( map%constraints%neq3 .gt. 0 ) then
                allocate(map%constraints%d3(map%constraints%neq3))
            else
                allocate(map%constraints%d3(1))
            endif
            map%constraints%d3=0.0_r8
            if (verbosity .gt. 0) then
                write (*, *) '... ', tochar(map%constraints%neq3), ' thirdorder constraints (', tochar(walltime() - t0), 's)'
            end if
        else
            map%constraints%neq3 = 0
        end if

        if (map%have_fc_quartet) then
            t0 = walltime()
            call lo_fourthorder_asr(map, map%constraints%eq4, map%constraints%neq4)
            if ( map%constraints%neq4 .gt. 0 ) then
                allocate(map%constraints%d4(map%constraints%neq4))
            else
                allocate(map%constraints%d4(1))
            endif
            map%constraints%d4=0.0_r8

            if (verbosity .gt. 0) then
                write (*, *) '... ', tochar(map%constraints%neq4), ' fourthorder constraints (', tochar(walltime() - t0), 's)'
            end if
        else
            map%constraints%neq4 = 0
        end if

        if (map%have_Z_triplet) then
            t0 = walltime()
            call lo_Ztriplet_asr(map, map%constraints%eqz3, map%constraints%neqz3)
            if ( map%constraints%neqz3 .gt. 0 ) then
                allocate(map%constraints%dz(map%constraints%neqz3))
            else
                allocate(map%constraints%dz(1))
            endif
            map%constraints%dz=0.0_r8

            if (verbosity .gt. 0) then
                write (*, *) '... ', tochar(map%constraints%neqz3), ' Z triplet constraints (', tochar(walltime() - t0), 's)'
            end if
        else
            map%constraints%neqz3 = 0
        end if

        ! And keep track of the total number of constraints
        map%constraints%nconstr_tot = map%constraints%neq1 + map%constraints%neq2 + map%constraints%neq3 + map%constraints%neq4
        ! Ifort gets angry sometimes if I send around unallocated matrices. So allocate them to 1.
        if (map%constraints%neq1 .eq. 0) then
            lo_allocate(map%constraints%eq1(1, 1))
            map%constraints%eq1 = -lo_huge
        end if
        if (map%constraints%neq2 .eq. 0) then
            lo_allocate(map%constraints%eq2(1, 1))
            map%constraints%eq2 = -lo_huge
        end if
        if (map%constraints%neq3 .eq. 0) then
            lo_allocate(map%constraints%eq3(1, 1))
            map%constraints%eq3 = -lo_huge
        end if
        if (map%constraints%neq4 .eq. 0) then
            lo_allocate(map%constraints%eq4(1, 1))
            map%constraints%eq4 = -lo_huge
        end if
        if (map%constraints%neqz3 .eq. 0) then
            lo_allocate(map%constraints%eqz3(1, 1))
            map%constraints%eqz3 = -lo_huge
        end if
    end block cleanconstr

    ! Then maybe the cross-constraints? Not now.

    if (verbosity .gt. 0) write (*, *) '... created ', tochar(map%constraints%nconstr_tot), ' constraints (', tochar(walltime() - timer), 's)'
end subroutine

!> Ensure that all invariances hold to very high precision. This is done via Lagrange multipliers, adding a small delta to each forceconstant, minimizing the relative deviation while strictly enforcing the constraints.
module subroutine enforce_constraints(map)
    !> forcemap
    class(lo_forcemap), intent(inout) :: map

    ! Enforce them one by one, smaller matrices to gels that way
    if (map%have_fc_singlet .and. map%constraints%neq1 .gt. 0) then
        call lo_enforce_linear_constraints(map%constraints%eq1, map%xuc%x_fc_singlet)
    end if
    if (map%have_fc_pair .and. map%constraints%neq2 .gt. 0) then
        call lo_enforce_linear_constraints(map%constraints%eq2, map%xuc%x_fc_pair)
    end if
    if (map%have_fc_triplet .and. map%constraints%neq3 .gt. 0) then
        call lo_enforce_linear_constraints(map%constraints%eq3, map%xuc%x_fc_triplet)
    end if
    if (map%have_fc_quartet .and. map%constraints%neq4 .gt. 0) then
        call lo_enforce_linear_constraints(map%constraints%eq4, map%xuc%x_fc_quartet)
    end if
end subroutine

!> Rotational, Hermitian and Huang invariances for the second order
module subroutine lo_secondorder_rot_herm_huang(map, uc, eq2, vD, neq2, rotational, huang, hermitian, hermitian_rhs, huang_rhs, rotational_rhs)
    !> symmetry stuff rearranged into coordination shells
    class(lo_forcemap), intent(in) :: map
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> constraints, C in Cx=D
    real(r8), dimension(:, :), allocatable, intent(out) :: eq2
    !> constraints, D in Cx=D
    real(r8), dimension(:), allocatable, intent(out) :: vD
    !> how many equations did I get?
    integer, intent(out) :: neq2
    !> which constraints should I care about?
    logical, intent(in) :: rotational, huang, hermitian
    !> rhs of hermitian constraints
    real(r8), dimension(:), intent(in) :: hermitian_rhs
    real(r8), dimension(:), intent(in) :: huang_rhs
    real(r8), dimension(:), intent(in) :: rotational_rhs

    real(r8), dimension(:, :), allocatable :: dumeq, dumeq_rot, dumeq_huang, dumeq_herm, m0,U,V
    real(r8), dimension(:), allocatable :: dum_rhs, S
    real(r8) :: wm9x1(9, 1), wm9x2(9, 2), wm9x3(9, 3), wm9x4(9, 4), wm9x5(9, 5), wm9x6(9, 6), wm9x7(9, 7), wm9x8(9, 8), wm9x9(9, 9)
    real(r8) :: m9x81(9, 81), m9x27(9, 27), m9x9(9, 9)
    real(r8), dimension(9, 9) :: CM
    real(r8), dimension(3) :: rij
    real(r8) :: f0
    integer, dimension(9) :: xind, xlocind
    integer :: eqctr, a1, a2, i, j
    integer :: unsh, unop, nfc, ipair, nx

    nfc = map%xuc%nx_fc_pair
    eqctr = 0
    allocate (dumeq(nfc, map%n_atom_uc*27 + 81 + 9*map%n_atom_uc))
    allocate (dumeq_rot(nfc, map%n_atom_uc*27))
    allocate (dumeq_huang(nfc, 81))
    allocate (dumeq_herm(nfc, map%n_atom_uc*9))
    allocate (dum_rhs( map%n_atom_uc*27 + 81 + 9*map%n_atom_uc))
    dumeq = 0.0_r8
    dumeq_rot = 0.0_r8
    dumeq_huang = 0.0_r8
    dumeq_herm = 0.0_r8
    dum_rhs=0.0_r8

    atomloop: do a1 = 1, map%n_atom_uc
        dumeq_rot = 0.0_r8
        pairloop: do ipair = 1, map%xuc%n_fc_pair
            if (map%xuc%fc_pair(ipair)%selfterm) cycle pairloop
            if (map%xuc%fc_pair(ipair)%i1 .ne. a1) cycle pairloop
            a2 = map%xuc%fc_pair(ipair)%i2
            unsh = map%xuc%fc_pair(ipair)%irreducible_shell
            unop = map%xuc%fc_pair(ipair)%operation_from_shell
            nx = map%fc_pair_shell(unsh)%nx
            xind(1:nx) = map%fc_pair_shell(unsh)%ind_global
            xlocind(1:nx) = map%fc_pair_shell(unsh)%ind_local
            CM = 0.0_r8
            select case (nx)
            case (1)
                call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x1)
                CM(:, xlocind(1:nx)) = wm9x1
            case (2)
                call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x2)
                CM(:, xlocind(1:nx)) = wm9x2
            case (3)
                call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x3)
                CM(:, xlocind(1:nx)) = wm9x3
            case (4)
                call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x4)
                CM(:, xlocind(1:nx)) = wm9x4
            case (5)
                call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x5)
                CM(:, xlocind(1:nx)) = wm9x5
            case (6)
                call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x6)
                CM(:, xlocind(1:nx)) = wm9x6
            case (7)
                call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x7)
                CM(:, xlocind(1:nx)) = wm9x7
            case (8)
                call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x8)
                CM(:, xlocind(1:nx)) = wm9x8
            case (9)
                call lo_gemm(map%op_pair(unop)%sotr, map%fc_pair_shell(unsh)%coeff, wm9x9)
                CM(:, xlocind(1:nx)) = wm9x9
            end select
            ! the pair vector
            rij = map%xuc%fc_pair(ipair)%flv + uc%r(:, a2) - uc%r(:, a1)
            rij = lo_chop(uc%fractional_to_cartesian(rij), lo_sqtol)
            ! add equations together
            if (rotational) then
                m9x27 = secondorder_rotational_coefficient(rij, CM)
                dumeq_rot(xind(1:nx), (a1-1)*27+1:a1*27) = dumeq_rot(xind(1:nx), (a1-1)*27+1:a1*27) + m9x27(xlocind(1:nx), :)
            end if
            if (huang) then
                m9x81 = secondorder_huang_coefficient(rij, CM)
                dumeq_huang(xind(1:nx), :) = dumeq_huang(xind(1:nx), :) + m9x81(xlocind(1:nx), :)
            end if
            if (hermitian) then
                m9x9 = secondorder_hermitian_character(CM)
                dumeq_herm(xind(1:nx), (a1-1)*9+1:a1*9) = dumeq_herm(xind(1:nx), (a1-1)*9+1:a1*9) + m9x9(xlocind(1:nx), :)
            end if
        end do pairloop
    end do atomloop
    ! Now stack the other equations to this
    if (rotational) then
        do i = 1, map%n_atom_uc*27
            eqctr = eqctr + 1
            dumeq(:, eqctr) = dumeq_rot(:, i)
            dum_rhs(eqctr) = rotational_rhs(i)
        end do
    end if
    if (huang) then
        do i = 1, 81
            eqctr = eqctr + 1
            dumeq(:, eqctr) = dumeq_huang(:, i)
            dum_rhs(eqctr) = huang_rhs(i)
            !write(*,*) i,eqctr,dum_rhs(eqctr),sum(abs(dumeq))
        end do
    end if
    if (hermitian) then
        do i=1,map%n_atom_uc*9
            eqctr = eqctr + 1
            dumeq(:, eqctr) = dumeq_herm(:, i)
            dum_rhs(eqctr) = hermitian_rhs(i)
            !write(*,*) i,eqctr,dum_rhs(eqctr),sum(abs(dumeq))
        enddo
    end if

    ! Now we want to reduce the equations
    allocate(m0(eqctr,nfc))
    m0=transpose(dumeq(:,1:eqctr))
    m0=lo_chop(m0,1E-13_r8)
    ! Some intermediate cleanup
    deallocate(dumeq)
    deallocate(dumeq_herm)
    deallocate(dumeq_rot)
    deallocate(dumeq_huang)

    ! Reduce equations to full rank
    allocate(S(max(eqctr,nfc)))
    allocate(U(eqctr,eqctr))
    allocate(V(nfc,nfc))
    S=0.0_r8
    U=0.0_r8
    V=0.0_r8
    call lo_dgesvd(m0,s,u,v)

    if ( s(1) .gt. 1E-15_r8 ) then
        ! Means we have non-trivial constraints
        j=0
        do i=1,size(s)
            if ( s(i)/s(1) .gt. 1E-10_r8 ) j=j+1
        enddo
        neq2=j

        allocate(vD(neq2))
        vD=0.0_r8

        ! Fix the RHS of the constraints
        call lo_gemv(U(:,1:neq2),dum_rhs(1:eqctr),vD,trans='T')
        ! And the actual coefficient matrix
        allocate(eq2(neq2,nfc))
        eq2=0.0_r8
        do i=1,neq2
            eq2(i,:)=S(i)*V(i,:)
        enddo
        eq2=lo_chop(eq2,1E-12_r8)
    else
        ! No relevant constraints
        neq2=0
        allocate(vD(1))
        vD=0.0_r8
    endif

contains
    pure function secondorder_rotational_coefficient(r, cl) result(m)
        real(r8), dimension(3), intent(in) :: r
        real(r8), dimension(9, 9), intent(in) :: cl
        real(r8), dimension(9, 27) :: m
        !
        real(r8) :: rx, ry, rz
        rx = r(1)
        ry = r(2)
        rz = r(3)
        m(:, 1) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 2) = [-(cl(2, 1)*rx) + cl(1, 1)*ry, -(cl(2, 2)*rx) + cl(1, 2)*ry, -(cl(2, 3)*rx) + cl(1, 3)*ry, -(cl(2, 4)*rx) + cl(1, 4)*ry, -(cl(2, 5)*rx) + cl(1, 5)*ry, -(cl(2, 6)*rx) + cl(1, 6)*ry, -(cl(2, 7)*rx) + cl(1, 7)*ry, -(cl(2, 8)*rx) + cl(1, 8)*ry, -(cl(2, 9)*rx) + cl(1, 9)*ry]
        m(:, 3) = [-(cl(3, 1)*rx) + cl(1, 1)*rz, -(cl(3, 2)*rx) + cl(1, 2)*rz, -(cl(3, 3)*rx) + cl(1, 3)*rz, -(cl(3, 4)*rx) + cl(1, 4)*rz, -(cl(3, 5)*rx) + cl(1, 5)*rz, -(cl(3, 6)*rx) + cl(1, 6)*rz, -(cl(3, 7)*rx) + cl(1, 7)*rz, -(cl(3, 8)*rx) + cl(1, 8)*rz, -(cl(3, 9)*rx) + cl(1, 9)*rz]
        m(:, 4) = [cl(2, 1)*rx - cl(1, 1)*ry, cl(2, 2)*rx - cl(1, 2)*ry, cl(2, 3)*rx - cl(1, 3)*ry, cl(2, 4)*rx - cl(1, 4)*ry, cl(2, 5)*rx - cl(1, 5)*ry, cl(2, 6)*rx - cl(1, 6)*ry, cl(2, 7)*rx - cl(1, 7)*ry, cl(2, 8)*rx - cl(1, 8)*ry, cl(2, 9)*rx - cl(1, 9)*ry]
        m(:, 5) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 6) = [-(cl(3, 1)*ry) + cl(2, 1)*rz, -(cl(3, 2)*ry) + cl(2, 2)*rz, -(cl(3, 3)*ry) + cl(2, 3)*rz, -(cl(3, 4)*ry) + cl(2, 4)*rz, -(cl(3, 5)*ry) + cl(2, 5)*rz, -(cl(3, 6)*ry) + cl(2, 6)*rz, -(cl(3, 7)*ry) + cl(2, 7)*rz, -(cl(3, 8)*ry) + cl(2, 8)*rz, -(cl(3, 9)*ry) + cl(2, 9)*rz]
        m(:, 7) = [cl(3, 1)*rx - cl(1, 1)*rz, cl(3, 2)*rx - cl(1, 2)*rz, cl(3, 3)*rx - cl(1, 3)*rz, cl(3, 4)*rx - cl(1, 4)*rz, cl(3, 5)*rx - cl(1, 5)*rz, cl(3, 6)*rx - cl(1, 6)*rz, cl(3, 7)*rx - cl(1, 7)*rz, cl(3, 8)*rx - cl(1, 8)*rz, cl(3, 9)*rx - cl(1, 9)*rz]
        m(:, 8) = [cl(3, 1)*ry - cl(2, 1)*rz, cl(3, 2)*ry - cl(2, 2)*rz, cl(3, 3)*ry - cl(2, 3)*rz, cl(3, 4)*ry - cl(2, 4)*rz, cl(3, 5)*ry - cl(2, 5)*rz, cl(3, 6)*ry - cl(2, 6)*rz, cl(3, 7)*ry - cl(2, 7)*rz, cl(3, 8)*ry - cl(2, 8)*rz, cl(3, 9)*ry - cl(2, 9)*rz]
        m(:, 9) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 10) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 11) = [-(cl(5, 1)*rx) + cl(4, 1)*ry, -(cl(5, 2)*rx) + cl(4, 2)*ry, -(cl(5, 3)*rx) + cl(4, 3)*ry, -(cl(5, 4)*rx) + cl(4, 4)*ry, -(cl(5, 5)*rx) + cl(4, 5)*ry, -(cl(5, 6)*rx) + cl(4, 6)*ry, -(cl(5, 7)*rx) + cl(4, 7)*ry, -(cl(5, 8)*rx) + cl(4, 8)*ry, -(cl(5, 9)*rx) + cl(4, 9)*ry]
        m(:, 12) = [-(cl(6, 1)*rx) + cl(4, 1)*rz, -(cl(6, 2)*rx) + cl(4, 2)*rz, -(cl(6, 3)*rx) + cl(4, 3)*rz, -(cl(6, 4)*rx) + cl(4, 4)*rz, -(cl(6, 5)*rx) + cl(4, 5)*rz, -(cl(6, 6)*rx) + cl(4, 6)*rz, -(cl(6, 7)*rx) + cl(4, 7)*rz, -(cl(6, 8)*rx) + cl(4, 8)*rz, -(cl(6, 9)*rx) + cl(4, 9)*rz]
        m(:, 13) = [cl(5, 1)*rx - cl(4, 1)*ry, cl(5, 2)*rx - cl(4, 2)*ry, cl(5, 3)*rx - cl(4, 3)*ry, cl(5, 4)*rx - cl(4, 4)*ry, cl(5, 5)*rx - cl(4, 5)*ry, cl(5, 6)*rx - cl(4, 6)*ry, cl(5, 7)*rx - cl(4, 7)*ry, cl(5, 8)*rx - cl(4, 8)*ry, cl(5, 9)*rx - cl(4, 9)*ry]
        m(:, 14) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 15) = [-(cl(6, 1)*ry) + cl(5, 1)*rz, -(cl(6, 2)*ry) + cl(5, 2)*rz, -(cl(6, 3)*ry) + cl(5, 3)*rz, -(cl(6, 4)*ry) + cl(5, 4)*rz, -(cl(6, 5)*ry) + cl(5, 5)*rz, -(cl(6, 6)*ry) + cl(5, 6)*rz, -(cl(6, 7)*ry) + cl(5, 7)*rz, -(cl(6, 8)*ry) + cl(5, 8)*rz, -(cl(6, 9)*ry) + cl(5, 9)*rz]
        m(:, 16) = [cl(6, 1)*rx - cl(4, 1)*rz, cl(6, 2)*rx - cl(4, 2)*rz, cl(6, 3)*rx - cl(4, 3)*rz, cl(6, 4)*rx - cl(4, 4)*rz, cl(6, 5)*rx - cl(4, 5)*rz, cl(6, 6)*rx - cl(4, 6)*rz, cl(6, 7)*rx - cl(4, 7)*rz, cl(6, 8)*rx - cl(4, 8)*rz, cl(6, 9)*rx - cl(4, 9)*rz]
        m(:, 17) = [cl(6, 1)*ry - cl(5, 1)*rz, cl(6, 2)*ry - cl(5, 2)*rz, cl(6, 3)*ry - cl(5, 3)*rz, cl(6, 4)*ry - cl(5, 4)*rz, cl(6, 5)*ry - cl(5, 5)*rz, cl(6, 6)*ry - cl(5, 6)*rz, cl(6, 7)*ry - cl(5, 7)*rz, cl(6, 8)*ry - cl(5, 8)*rz, cl(6, 9)*ry - cl(5, 9)*rz]
        m(:, 18) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 19) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 20) = [-(cl(8, 1)*rx) + cl(7, 1)*ry, -(cl(8, 2)*rx) + cl(7, 2)*ry, -(cl(8, 3)*rx) + cl(7, 3)*ry, -(cl(8, 4)*rx) + cl(7, 4)*ry, -(cl(8, 5)*rx) + cl(7, 5)*ry, -(cl(8, 6)*rx) + cl(7, 6)*ry, -(cl(8, 7)*rx) + cl(7, 7)*ry, -(cl(8, 8)*rx) + cl(7, 8)*ry, -(cl(8, 9)*rx) + cl(7, 9)*ry]
        m(:, 21) = [-(cl(9, 1)*rx) + cl(7, 1)*rz, -(cl(9, 2)*rx) + cl(7, 2)*rz, -(cl(9, 3)*rx) + cl(7, 3)*rz, -(cl(9, 4)*rx) + cl(7, 4)*rz, -(cl(9, 5)*rx) + cl(7, 5)*rz, -(cl(9, 6)*rx) + cl(7, 6)*rz, -(cl(9, 7)*rx) + cl(7, 7)*rz, -(cl(9, 8)*rx) + cl(7, 8)*rz, -(cl(9, 9)*rx) + cl(7, 9)*rz]
        m(:, 22) = [cl(8, 1)*rx - cl(7, 1)*ry, cl(8, 2)*rx - cl(7, 2)*ry, cl(8, 3)*rx - cl(7, 3)*ry, cl(8, 4)*rx - cl(7, 4)*ry, cl(8, 5)*rx - cl(7, 5)*ry, cl(8, 6)*rx - cl(7, 6)*ry, cl(8, 7)*rx - cl(7, 7)*ry, cl(8, 8)*rx - cl(7, 8)*ry, cl(8, 9)*rx - cl(7, 9)*ry]
        m(:, 23) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 24) = [-(cl(9, 1)*ry) + cl(8, 1)*rz, -(cl(9, 2)*ry) + cl(8, 2)*rz, -(cl(9, 3)*ry) + cl(8, 3)*rz, -(cl(9, 4)*ry) + cl(8, 4)*rz, -(cl(9, 5)*ry) + cl(8, 5)*rz, -(cl(9, 6)*ry) + cl(8, 6)*rz, -(cl(9, 7)*ry) + cl(8, 7)*rz, -(cl(9, 8)*ry) + cl(8, 8)*rz, -(cl(9, 9)*ry) + cl(8, 9)*rz]
        m(:, 25) = [cl(9, 1)*rx - cl(7, 1)*rz, cl(9, 2)*rx - cl(7, 2)*rz, cl(9, 3)*rx - cl(7, 3)*rz, cl(9, 4)*rx - cl(7, 4)*rz, cl(9, 5)*rx - cl(7, 5)*rz, cl(9, 6)*rx - cl(7, 6)*rz, cl(9, 7)*rx - cl(7, 7)*rz, cl(9, 8)*rx - cl(7, 8)*rz, cl(9, 9)*rx - cl(7, 9)*rz]
        m(:, 26) = [cl(9, 1)*ry - cl(8, 1)*rz, cl(9, 2)*ry - cl(8, 2)*rz, cl(9, 3)*ry - cl(8, 3)*rz, cl(9, 4)*ry - cl(8, 4)*rz, cl(9, 5)*ry - cl(8, 5)*rz, cl(9, 6)*ry - cl(8, 6)*rz, cl(9, 7)*ry - cl(8, 7)*rz, cl(9, 8)*ry - cl(8, 8)*rz, cl(9, 9)*ry - cl(8, 9)*rz]
        m(:, 27) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
    end function
    pure function secondorder_hermitian_character(cl) result(m)
        real(r8), dimension(9, 9), intent(in) :: cl
        real(r8), dimension(9, 9) :: m
        m(:, 1) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 2) = [cl(2, 1) - cl(4, 1), cl(2, 2) - cl(4, 2), cl(2, 3) - cl(4, 3), &
                   cl(2, 4) - cl(4, 4), cl(2, 5) - cl(4, 5), cl(2, 6) - cl(4, 6), &
                   cl(2, 7) - cl(4, 7), cl(2, 8) - cl(4, 8), cl(2, 9) - cl(4, 9)]
        m(:, 3) = [cl(3, 1) - cl(7, 1), cl(3, 2) - cl(7, 2), cl(3, 3) - cl(7, 3), &
                   cl(3, 4) - cl(7, 4), cl(3, 5) - cl(7, 5), cl(3, 6) - cl(7, 6), &
                   cl(3, 7) - cl(7, 7), cl(3, 8) - cl(7, 8), cl(3, 9) - cl(7, 9)]
        m(:, 4) = [-cl(2, 1) + cl(4, 1), -cl(2, 2) + cl(4, 2), -cl(2, 3) + cl(4, 3), &
                   -cl(2, 4) + cl(4, 4), -cl(2, 5) + cl(4, 5), -cl(2, 6) + cl(4, 6), &
                   -cl(2, 7) + cl(4, 7), -cl(2, 8) + cl(4, 8), -cl(2, 9) + cl(4, 9)]
        m(:, 5) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 6) = [cl(6, 1) - cl(8, 1), cl(6, 2) - cl(8, 2), cl(6, 3) - cl(8, 3), &
                   cl(6, 4) - cl(8, 4), cl(6, 5) - cl(8, 5), cl(6, 6) - cl(8, 6), &
                   cl(6, 7) - cl(8, 7), cl(6, 8) - cl(8, 8), cl(6, 9) - cl(8, 9)]
        m(:, 7) = [-cl(3, 1) + cl(7, 1), -cl(3, 2) + cl(7, 2), -cl(3, 3) + cl(7, 3), &
                   -cl(3, 4) + cl(7, 4), -cl(3, 5) + cl(7, 5), -cl(3, 6) + cl(7, 6), &
                   -cl(3, 7) + cl(7, 7), -cl(3, 8) + cl(7, 8), -cl(3, 9) + cl(7, 9)]
        m(:, 8) = [-cl(6, 1) + cl(8, 1), -cl(6, 2) + cl(8, 2), -cl(6, 3) + cl(8, 3), &
                   -cl(6, 4) + cl(8, 4), -cl(6, 5) + cl(8, 5), -cl(6, 6) + cl(8, 6), &
                   -cl(6, 7) + cl(8, 7), -cl(6, 8) + cl(8, 8), -cl(6, 9) + cl(8, 9)]
        m(:, 9) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8]
    end function
    function secondorder_huang_coefficient(r, cl) result(m)
        real(r8), dimension(3), intent(in) :: r
        real(r8), dimension(9, 9), intent(in) :: cl
        real(r8), dimension(9, 81) :: m
        !
        real(r8) :: rx, ry, rz
        rx = r(1)
        ry = r(2)
        rz = r(3)

        m(:, 1) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                   0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 2) = [-(cl(2, 1)*rx**2) + cl(1, 1)*rx*ry, -(cl(2, 2)*rx**2) + &
                   cl(1, 2)*rx*ry, -(cl(2, 3)*rx**2) + cl(1, 3)*rx*ry, -(cl(2, 4)*rx**2) + &
                   cl(1, 4)*rx*ry, -(cl(2, 5)*rx**2) + cl(1, 5)*rx*ry, -(cl(2, 6)*rx**2) + &
                   cl(1, 6)*rx*ry, -(cl(2, 7)*rx**2) + cl(1, 7)*rx*ry, -(cl(2, 8)*rx**2) + &
                   cl(1, 8)*rx*ry, -(cl(2, 9)*rx**2) + cl(1, 9)*rx*ry]
        m(:, 3) = [-(cl(3, 1)*rx**2) + cl(1, 1)*rx*rz, -(cl(3, 2)*rx**2) + &
                   cl(1, 2)*rx*rz, -(cl(3, 3)*rx**2) + cl(1, 3)*rx*rz, -(cl(3, 4)*rx**2) + &
                   cl(1, 4)*rx*rz, -(cl(3, 5)*rx**2) + cl(1, 5)*rx*rz, -(cl(3, 6)*rx**2) + &
                   cl(1, 6)*rx*rz, -(cl(3, 7)*rx**2) + cl(1, 7)*rx*rz, -(cl(3, 8)*rx**2) + &
                   cl(1, 8)*rx*rz, -(cl(3, 9)*rx**2) + cl(1, 9)*rx*rz]
        m(:, 4) = [-(cl(4, 1)*rx**2) + cl(1, 1)*rx*ry, -(cl(4, 2)*rx**2) + &
                   cl(1, 2)*rx*ry, -(cl(4, 3)*rx**2) + cl(1, 3)*rx*ry, -(cl(4, 4)*rx**2) + &
                   cl(1, 4)*rx*ry, -(cl(4, 5)*rx**2) + cl(1, 5)*rx*ry, -(cl(4, 6)*rx**2) + &
                   cl(1, 6)*rx*ry, -(cl(4, 7)*rx**2) + cl(1, 7)*rx*ry, -(cl(4, 8)*rx**2) + &
                   cl(1, 8)*rx*ry, -(cl(4, 9)*rx**2) + cl(1, 9)*rx*ry]
        m(:, 5) = [-(cl(5, 1)*rx**2) + cl(1, 1)*ry**2, -(cl(5, 2)*rx**2) + &
                   cl(1, 2)*ry**2, -(cl(5, 3)*rx**2) + cl(1, 3)*ry**2, -(cl(5, 4)*rx**2) + &
                   cl(1, 4)*ry**2, -(cl(5, 5)*rx**2) + cl(1, 5)*ry**2, -(cl(5, 6)*rx**2) + &
                   cl(1, 6)*ry**2, -(cl(5, 7)*rx**2) + cl(1, 7)*ry**2, -(cl(5, 8)*rx**2) + &
                   cl(1, 8)*ry**2, -(cl(5, 9)*rx**2) + cl(1, 9)*ry**2]
        m(:, 6) = [-(cl(6, 1)*rx**2) + cl(1, 1)*ry*rz, -(cl(6, 2)*rx**2) + &
                   cl(1, 2)*ry*rz, -(cl(6, 3)*rx**2) + cl(1, 3)*ry*rz, -(cl(6, 4)*rx**2) + &
                   cl(1, 4)*ry*rz, -(cl(6, 5)*rx**2) + cl(1, 5)*ry*rz, -(cl(6, 6)*rx**2) + &
                   cl(1, 6)*ry*rz, -(cl(6, 7)*rx**2) + cl(1, 7)*ry*rz, -(cl(6, 8)*rx**2) + &
                   cl(1, 8)*ry*rz, -(cl(6, 9)*rx**2) + cl(1, 9)*ry*rz]
        m(:, 7) = [-(cl(7, 1)*rx**2) + cl(1, 1)*rx*rz, -(cl(7, 2)*rx**2) + &
                   cl(1, 2)*rx*rz, -(cl(7, 3)*rx**2) + cl(1, 3)*rx*rz, -(cl(7, 4)*rx**2) + &
                   cl(1, 4)*rx*rz, -(cl(7, 5)*rx**2) + cl(1, 5)*rx*rz, -(cl(7, 6)*rx**2) + &
                   cl(1, 6)*rx*rz, -(cl(7, 7)*rx**2) + cl(1, 7)*rx*rz, -(cl(7, 8)*rx**2) + &
                   cl(1, 8)*rx*rz, -(cl(7, 9)*rx**2) + cl(1, 9)*rx*rz]
        m(:, 8) = [-(cl(8, 1)*rx**2) + cl(1, 1)*ry*rz, -(cl(8, 2)*rx**2) + &
                   cl(1, 2)*ry*rz, -(cl(8, 3)*rx**2) + cl(1, 3)*ry*rz, -(cl(8, 4)*rx**2) + &
                   cl(1, 4)*ry*rz, -(cl(8, 5)*rx**2) + cl(1, 5)*ry*rz, -(cl(8, 6)*rx**2) + &
                   cl(1, 6)*ry*rz, -(cl(8, 7)*rx**2) + cl(1, 7)*ry*rz, -(cl(8, 8)*rx**2) + &
                   cl(1, 8)*ry*rz, -(cl(8, 9)*rx**2) + cl(1, 9)*ry*rz]
        m(:, 9) = [-(cl(9, 1)*rx**2) + cl(1, 1)*rz**2, -(cl(9, 2)*rx**2) + &
                   cl(1, 2)*rz**2, -(cl(9, 3)*rx**2) + cl(1, 3)*rz**2, -(cl(9, 4)*rx**2) + &
                   cl(1, 4)*rz**2, -(cl(9, 5)*rx**2) + cl(1, 5)*rz**2, -(cl(9, 6)*rx**2) + &
                   cl(1, 6)*rz**2, -(cl(9, 7)*rx**2) + cl(1, 7)*rz**2, -(cl(9, 8)*rx**2) + &
                   cl(1, 8)*rz**2, -(cl(9, 9)*rx**2) + cl(1, 9)*rz**2]
        m(:, 10) = [cl(2, 1)*rx**2 - cl(1, 1)*rx*ry, cl(2, 2)*rx**2 - cl(1, 2)*rx*ry, &
                    cl(2, 3)*rx**2 - cl(1, 3)*rx*ry, cl(2, 4)*rx**2 - cl(1, 4)*rx*ry, &
                    cl(2, 5)*rx**2 - cl(1, 5)*rx*ry, cl(2, 6)*rx**2 - cl(1, 6)*rx*ry, &
                    cl(2, 7)*rx**2 - cl(1, 7)*rx*ry, cl(2, 8)*rx**2 - cl(1, 8)*rx*ry, &
                    cl(2, 9)*rx**2 - cl(1, 9)*rx*ry]
        m(:, 11) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                    0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 12) = [-(cl(3, 1)*rx*ry) + cl(2, 1)*rx*rz, -(cl(3, 2)*rx*ry) + &
                    cl(2, 2)*rx*rz, -(cl(3, 3)*rx*ry) + cl(2, 3)*rx*rz, -(cl(3, 4)*rx*ry) + &
                    cl(2, 4)*rx*rz, -(cl(3, 5)*rx*ry) + cl(2, 5)*rx*rz, -(cl(3, 6)*rx*ry) + &
                    cl(2, 6)*rx*rz, -(cl(3, 7)*rx*ry) + cl(2, 7)*rx*rz, -(cl(3, 8)*rx*ry) + &
                    cl(2, 8)*rx*rz, -(cl(3, 9)*rx*ry) + cl(2, 9)*rx*rz]
        m(:, 13) = [cl(2, 1)*rx*ry - cl(4, 1)*rx*ry, cl(2, 2)*rx*ry - &
                    cl(4, 2)*rx*ry, cl(2, 3)*rx*ry - cl(4, 3)*rx*ry, cl(2, 4)*rx*ry - &
                    cl(4, 4)*rx*ry, cl(2, 5)*rx*ry - cl(4, 5)*rx*ry, cl(2, 6)*rx*ry - &
                    cl(4, 6)*rx*ry, cl(2, 7)*rx*ry - cl(4, 7)*rx*ry, cl(2, 8)*rx*ry - &
                    cl(4, 8)*rx*ry, cl(2, 9)*rx*ry - cl(4, 9)*rx*ry]
        m(:, 14) = [-(cl(5, 1)*rx*ry) + cl(2, 1)*ry**2, -(cl(5, 2)*rx*ry) + &
                    cl(2, 2)*ry**2, -(cl(5, 3)*rx*ry) + cl(2, 3)*ry**2, -(cl(5, 4)*rx*ry) + &
                    cl(2, 4)*ry**2, -(cl(5, 5)*rx*ry) + cl(2, 5)*ry**2, -(cl(5, 6)*rx*ry) + &
                    cl(2, 6)*ry**2, -(cl(5, 7)*rx*ry) + cl(2, 7)*ry**2, -(cl(5, 8)*rx*ry) + &
                    cl(2, 8)*ry**2, -(cl(5, 9)*rx*ry) + cl(2, 9)*ry**2]
        m(:, 15) = [-(cl(6, 1)*rx*ry) + cl(2, 1)*ry*rz, -(cl(6, 2)*rx*ry) + &
                    cl(2, 2)*ry*rz, -(cl(6, 3)*rx*ry) + cl(2, 3)*ry*rz, -(cl(6, 4)*rx*ry) + &
                    cl(2, 4)*ry*rz, -(cl(6, 5)*rx*ry) + cl(2, 5)*ry*rz, -(cl(6, 6)*rx*ry) + &
                    cl(2, 6)*ry*rz, -(cl(6, 7)*rx*ry) + cl(2, 7)*ry*rz, -(cl(6, 8)*rx*ry) + &
                    cl(2, 8)*ry*rz, -(cl(6, 9)*rx*ry) + cl(2, 9)*ry*rz]
        m(:, 16) = [-(cl(7, 1)*rx*ry) + cl(2, 1)*rx*rz, -(cl(7, 2)*rx*ry) + &
                    cl(2, 2)*rx*rz, -(cl(7, 3)*rx*ry) + cl(2, 3)*rx*rz, -(cl(7, 4)*rx*ry) + &
                    cl(2, 4)*rx*rz, -(cl(7, 5)*rx*ry) + cl(2, 5)*rx*rz, -(cl(7, 6)*rx*ry) + &
                    cl(2, 6)*rx*rz, -(cl(7, 7)*rx*ry) + cl(2, 7)*rx*rz, -(cl(7, 8)*rx*ry) + &
                    cl(2, 8)*rx*rz, -(cl(7, 9)*rx*ry) + cl(2, 9)*rx*rz]
        m(:, 17) = [-(cl(8, 1)*rx*ry) + cl(2, 1)*ry*rz, -(cl(8, 2)*rx*ry) + &
                    cl(2, 2)*ry*rz, -(cl(8, 3)*rx*ry) + cl(2, 3)*ry*rz, -(cl(8, 4)*rx*ry) + &
                    cl(2, 4)*ry*rz, -(cl(8, 5)*rx*ry) + cl(2, 5)*ry*rz, -(cl(8, 6)*rx*ry) + &
                    cl(2, 6)*ry*rz, -(cl(8, 7)*rx*ry) + cl(2, 7)*ry*rz, -(cl(8, 8)*rx*ry) + &
                    cl(2, 8)*ry*rz, -(cl(8, 9)*rx*ry) + cl(2, 9)*ry*rz]
        m(:, 18) = [-(cl(9, 1)*rx*ry) + cl(2, 1)*rz**2, -(cl(9, 2)*rx*ry) + &
                    cl(2, 2)*rz**2, -(cl(9, 3)*rx*ry) + cl(2, 3)*rz**2, -(cl(9, 4)*rx*ry) + &
                    cl(2, 4)*rz**2, -(cl(9, 5)*rx*ry) + cl(2, 5)*rz**2, -(cl(9, 6)*rx*ry) + &
                    cl(2, 6)*rz**2, -(cl(9, 7)*rx*ry) + cl(2, 7)*rz**2, -(cl(9, 8)*rx*ry) + &
                    cl(2, 8)*rz**2, -(cl(9, 9)*rx*ry) + cl(2, 9)*rz**2]
        m(:, 19) = [cl(3, 1)*rx**2 - cl(1, 1)*rx*rz, cl(3, 2)*rx**2 - cl(1, 2)*rx*rz, &
                    cl(3, 3)*rx**2 - cl(1, 3)*rx*rz, cl(3, 4)*rx**2 - cl(1, 4)*rx*rz, &
                    cl(3, 5)*rx**2 - cl(1, 5)*rx*rz, cl(3, 6)*rx**2 - cl(1, 6)*rx*rz, &
                    cl(3, 7)*rx**2 - cl(1, 7)*rx*rz, cl(3, 8)*rx**2 - cl(1, 8)*rx*rz, &
                    cl(3, 9)*rx**2 - cl(1, 9)*rx*rz]
        m(:, 20) = [cl(3, 1)*rx*ry - cl(2, 1)*rx*rz, cl(3, 2)*rx*ry - &
                    cl(2, 2)*rx*rz, cl(3, 3)*rx*ry - cl(2, 3)*rx*rz, cl(3, 4)*rx*ry - &
                    cl(2, 4)*rx*rz, cl(3, 5)*rx*ry - cl(2, 5)*rx*rz, cl(3, 6)*rx*ry - &
                    cl(2, 6)*rx*rz, cl(3, 7)*rx*ry - cl(2, 7)*rx*rz, cl(3, 8)*rx*ry - &
                    cl(2, 8)*rx*rz, cl(3, 9)*rx*ry - cl(2, 9)*rx*rz]
        m(:, 21) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                    0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 22) = [cl(3, 1)*rx*ry - cl(4, 1)*rx*rz, cl(3, 2)*rx*ry - &
                    cl(4, 2)*rx*rz, cl(3, 3)*rx*ry - cl(4, 3)*rx*rz, cl(3, 4)*rx*ry - &
                    cl(4, 4)*rx*rz, cl(3, 5)*rx*ry - cl(4, 5)*rx*rz, cl(3, 6)*rx*ry - &
                    cl(4, 6)*rx*rz, cl(3, 7)*rx*ry - cl(4, 7)*rx*rz, cl(3, 8)*rx*ry - &
                    cl(4, 8)*rx*rz, cl(3, 9)*rx*ry - cl(4, 9)*rx*rz]
        m(:, 23) = [cl(3, 1)*ry**2 - cl(5, 1)*rx*rz, cl(3, 2)*ry**2 - cl(5, 2)*rx*rz, &
                    cl(3, 3)*ry**2 - cl(5, 3)*rx*rz, cl(3, 4)*ry**2 - cl(5, 4)*rx*rz, &
                    cl(3, 5)*ry**2 - cl(5, 5)*rx*rz, cl(3, 6)*ry**2 - cl(5, 6)*rx*rz, &
                    cl(3, 7)*ry**2 - cl(5, 7)*rx*rz, cl(3, 8)*ry**2 - cl(5, 8)*rx*rz, &
                    cl(3, 9)*ry**2 - cl(5, 9)*rx*rz]
        m(:, 24) = [-(cl(6, 1)*rx*rz) + cl(3, 1)*ry*rz, -(cl(6, 2)*rx*rz) + &
                    cl(3, 2)*ry*rz, -(cl(6, 3)*rx*rz) + cl(3, 3)*ry*rz, -(cl(6, 4)*rx*rz) + &
                    cl(3, 4)*ry*rz, -(cl(6, 5)*rx*rz) + cl(3, 5)*ry*rz, -(cl(6, 6)*rx*rz) + &
                    cl(3, 6)*ry*rz, -(cl(6, 7)*rx*rz) + cl(3, 7)*ry*rz, -(cl(6, 8)*rx*rz) + &
                    cl(3, 8)*ry*rz, -(cl(6, 9)*rx*rz) + cl(3, 9)*ry*rz]
        m(:, 25) = [cl(3, 1)*rx*rz - cl(7, 1)*rx*rz, cl(3, 2)*rx*rz - &
                    cl(7, 2)*rx*rz, cl(3, 3)*rx*rz - cl(7, 3)*rx*rz, cl(3, 4)*rx*rz - &
                    cl(7, 4)*rx*rz, cl(3, 5)*rx*rz - cl(7, 5)*rx*rz, cl(3, 6)*rx*rz - &
                    cl(7, 6)*rx*rz, cl(3, 7)*rx*rz - cl(7, 7)*rx*rz, cl(3, 8)*rx*rz - &
                    cl(7, 8)*rx*rz, cl(3, 9)*rx*rz - cl(7, 9)*rx*rz]
        m(:, 26) = [-(cl(8, 1)*rx*rz) + cl(3, 1)*ry*rz, -(cl(8, 2)*rx*rz) + &
                    cl(3, 2)*ry*rz, -(cl(8, 3)*rx*rz) + cl(3, 3)*ry*rz, -(cl(8, 4)*rx*rz) + &
                    cl(3, 4)*ry*rz, -(cl(8, 5)*rx*rz) + cl(3, 5)*ry*rz, -(cl(8, 6)*rx*rz) + &
                    cl(3, 6)*ry*rz, -(cl(8, 7)*rx*rz) + cl(3, 7)*ry*rz, -(cl(8, 8)*rx*rz) + &
                    cl(3, 8)*ry*rz, -(cl(8, 9)*rx*rz) + cl(3, 9)*ry*rz]
        m(:, 27) = [-(cl(9, 1)*rx*rz) + cl(3, 1)*rz**2, -(cl(9, 2)*rx*rz) + &
                    cl(3, 2)*rz**2, -(cl(9, 3)*rx*rz) + cl(3, 3)*rz**2, -(cl(9, 4)*rx*rz) + &
                    cl(3, 4)*rz**2, -(cl(9, 5)*rx*rz) + cl(3, 5)*rz**2, -(cl(9, 6)*rx*rz) + &
                    cl(3, 6)*rz**2, -(cl(9, 7)*rx*rz) + cl(3, 7)*rz**2, -(cl(9, 8)*rx*rz) + &
                    cl(3, 8)*rz**2, -(cl(9, 9)*rx*rz) + cl(3, 9)*rz**2]
        m(:, 28) = [cl(4, 1)*rx**2 - cl(1, 1)*rx*ry, cl(4, 2)*rx**2 - cl(1, 2)*rx*ry, &
                    cl(4, 3)*rx**2 - cl(1, 3)*rx*ry, cl(4, 4)*rx**2 - cl(1, 4)*rx*ry, &
                    cl(4, 5)*rx**2 - cl(1, 5)*rx*ry, cl(4, 6)*rx**2 - cl(1, 6)*rx*ry, &
                    cl(4, 7)*rx**2 - cl(1, 7)*rx*ry, cl(4, 8)*rx**2 - cl(1, 8)*rx*ry, &
                    cl(4, 9)*rx**2 - cl(1, 9)*rx*ry]
        m(:, 29) = [-(cl(2, 1)*rx*ry) + cl(4, 1)*rx*ry, -(cl(2, 2)*rx*ry) + &
                    cl(4, 2)*rx*ry, -(cl(2, 3)*rx*ry) + cl(4, 3)*rx*ry, -(cl(2, 4)*rx*ry) + &
                    cl(4, 4)*rx*ry, -(cl(2, 5)*rx*ry) + cl(4, 5)*rx*ry, -(cl(2, 6)*rx*ry) + &
                    cl(4, 6)*rx*ry, -(cl(2, 7)*rx*ry) + cl(4, 7)*rx*ry, -(cl(2, 8)*rx*ry) + &
                    cl(4, 8)*rx*ry, -(cl(2, 9)*rx*ry) + cl(4, 9)*rx*ry]
        m(:, 30) = [-(cl(3, 1)*rx*ry) + cl(4, 1)*rx*rz, -(cl(3, 2)*rx*ry) + &
                    cl(4, 2)*rx*rz, -(cl(3, 3)*rx*ry) + cl(4, 3)*rx*rz, -(cl(3, 4)*rx*ry) + &
                    cl(4, 4)*rx*rz, -(cl(3, 5)*rx*ry) + cl(4, 5)*rx*rz, -(cl(3, 6)*rx*ry) + &
                    cl(4, 6)*rx*rz, -(cl(3, 7)*rx*ry) + cl(4, 7)*rx*rz, -(cl(3, 8)*rx*ry) + &
                    cl(4, 8)*rx*rz, -(cl(3, 9)*rx*ry) + cl(4, 9)*rx*rz]
        m(:, 31) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                    0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 32) = [-(cl(5, 1)*rx*ry) + cl(4, 1)*ry**2, -(cl(5, 2)*rx*ry) + &
                    cl(4, 2)*ry**2, -(cl(5, 3)*rx*ry) + cl(4, 3)*ry**2, -(cl(5, 4)*rx*ry) + &
                    cl(4, 4)*ry**2, -(cl(5, 5)*rx*ry) + cl(4, 5)*ry**2, -(cl(5, 6)*rx*ry) + &
                    cl(4, 6)*ry**2, -(cl(5, 7)*rx*ry) + cl(4, 7)*ry**2, -(cl(5, 8)*rx*ry) + &
                    cl(4, 8)*ry**2, -(cl(5, 9)*rx*ry) + cl(4, 9)*ry**2]
        m(:, 33) = [-(cl(6, 1)*rx*ry) + cl(4, 1)*ry*rz, -(cl(6, 2)*rx*ry) + &
                    cl(4, 2)*ry*rz, -(cl(6, 3)*rx*ry) + cl(4, 3)*ry*rz, -(cl(6, 4)*rx*ry) + &
                    cl(4, 4)*ry*rz, -(cl(6, 5)*rx*ry) + cl(4, 5)*ry*rz, -(cl(6, 6)*rx*ry) + &
                    cl(4, 6)*ry*rz, -(cl(6, 7)*rx*ry) + cl(4, 7)*ry*rz, -(cl(6, 8)*rx*ry) + &
                    cl(4, 8)*ry*rz, -(cl(6, 9)*rx*ry) + cl(4, 9)*ry*rz]
        m(:, 34) = [-(cl(7, 1)*rx*ry) + cl(4, 1)*rx*rz, -(cl(7, 2)*rx*ry) + &
                    cl(4, 2)*rx*rz, -(cl(7, 3)*rx*ry) + cl(4, 3)*rx*rz, -(cl(7, 4)*rx*ry) + &
                    cl(4, 4)*rx*rz, -(cl(7, 5)*rx*ry) + cl(4, 5)*rx*rz, -(cl(7, 6)*rx*ry) + &
                    cl(4, 6)*rx*rz, -(cl(7, 7)*rx*ry) + cl(4, 7)*rx*rz, -(cl(7, 8)*rx*ry) + &
                    cl(4, 8)*rx*rz, -(cl(7, 9)*rx*ry) + cl(4, 9)*rx*rz]
        m(:, 35) = [-(cl(8, 1)*rx*ry) + cl(4, 1)*ry*rz, -(cl(8, 2)*rx*ry) + &
                    cl(4, 2)*ry*rz, -(cl(8, 3)*rx*ry) + cl(4, 3)*ry*rz, -(cl(8, 4)*rx*ry) + &
                    cl(4, 4)*ry*rz, -(cl(8, 5)*rx*ry) + cl(4, 5)*ry*rz, -(cl(8, 6)*rx*ry) + &
                    cl(4, 6)*ry*rz, -(cl(8, 7)*rx*ry) + cl(4, 7)*ry*rz, -(cl(8, 8)*rx*ry) + &
                    cl(4, 8)*ry*rz, -(cl(8, 9)*rx*ry) + cl(4, 9)*ry*rz]
        m(:, 36) = [-(cl(9, 1)*rx*ry) + cl(4, 1)*rz**2, -(cl(9, 2)*rx*ry) + &
                    cl(4, 2)*rz**2, -(cl(9, 3)*rx*ry) + cl(4, 3)*rz**2, -(cl(9, 4)*rx*ry) + &
                    cl(4, 4)*rz**2, -(cl(9, 5)*rx*ry) + cl(4, 5)*rz**2, -(cl(9, 6)*rx*ry) + &
                    cl(4, 6)*rz**2, -(cl(9, 7)*rx*ry) + cl(4, 7)*rz**2, -(cl(9, 8)*rx*ry) + &
                    cl(4, 8)*rz**2, -(cl(9, 9)*rx*ry) + cl(4, 9)*rz**2]
        m(:, 37) = [cl(5, 1)*rx**2 - cl(1, 1)*ry**2, cl(5, 2)*rx**2 - cl(1, 2)*ry**2, &
                    cl(5, 3)*rx**2 - cl(1, 3)*ry**2, cl(5, 4)*rx**2 - cl(1, 4)*ry**2, &
                    cl(5, 5)*rx**2 - cl(1, 5)*ry**2, cl(5, 6)*rx**2 - cl(1, 6)*ry**2, &
                    cl(5, 7)*rx**2 - cl(1, 7)*ry**2, cl(5, 8)*rx**2 - cl(1, 8)*ry**2, &
                    cl(5, 9)*rx**2 - cl(1, 9)*ry**2]
        m(:, 38) = [cl(5, 1)*rx*ry - cl(2, 1)*ry**2, cl(5, 2)*rx*ry - cl(2, 2)*ry**2, &
                    cl(5, 3)*rx*ry - cl(2, 3)*ry**2, cl(5, 4)*rx*ry - cl(2, 4)*ry**2, &
                    cl(5, 5)*rx*ry - cl(2, 5)*ry**2, cl(5, 6)*rx*ry - cl(2, 6)*ry**2, &
                    cl(5, 7)*rx*ry - cl(2, 7)*ry**2, cl(5, 8)*rx*ry - cl(2, 8)*ry**2, &
                    cl(5, 9)*rx*ry - cl(2, 9)*ry**2]
        m(:, 39) = [-(cl(3, 1)*ry**2) + cl(5, 1)*rx*rz, -(cl(3, 2)*ry**2) + &
                    cl(5, 2)*rx*rz, -(cl(3, 3)*ry**2) + cl(5, 3)*rx*rz, -(cl(3, 4)*ry**2) + &
                    cl(5, 4)*rx*rz, -(cl(3, 5)*ry**2) + cl(5, 5)*rx*rz, -(cl(3, 6)*ry**2) + &
                    cl(5, 6)*rx*rz, -(cl(3, 7)*ry**2) + cl(5, 7)*rx*rz, -(cl(3, 8)*ry**2) + &
                    cl(5, 8)*rx*rz, -(cl(3, 9)*ry**2) + cl(5, 9)*rx*rz]
        m(:, 40) = [cl(5, 1)*rx*ry - cl(4, 1)*ry**2, cl(5, 2)*rx*ry - cl(4, 2)*ry**2, &
                    cl(5, 3)*rx*ry - cl(4, 3)*ry**2, cl(5, 4)*rx*ry - cl(4, 4)*ry**2, &
                    cl(5, 5)*rx*ry - cl(4, 5)*ry**2, cl(5, 6)*rx*ry - cl(4, 6)*ry**2, &
                    cl(5, 7)*rx*ry - cl(4, 7)*ry**2, cl(5, 8)*rx*ry - cl(4, 8)*ry**2, &
                    cl(5, 9)*rx*ry - cl(4, 9)*ry**2]
        m(:, 41) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                    0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 42) = [-(cl(6, 1)*ry**2) + cl(5, 1)*ry*rz, -(cl(6, 2)*ry**2) + &
                    cl(5, 2)*ry*rz, -(cl(6, 3)*ry**2) + cl(5, 3)*ry*rz, -(cl(6, 4)*ry**2) + &
                    cl(5, 4)*ry*rz, -(cl(6, 5)*ry**2) + cl(5, 5)*ry*rz, -(cl(6, 6)*ry**2) + &
                    cl(5, 6)*ry*rz, -(cl(6, 7)*ry**2) + cl(5, 7)*ry*rz, -(cl(6, 8)*ry**2) + &
                    cl(5, 8)*ry*rz, -(cl(6, 9)*ry**2) + cl(5, 9)*ry*rz]
        m(:, 43) = [-(cl(7, 1)*ry**2) + cl(5, 1)*rx*rz, -(cl(7, 2)*ry**2) + &
                    cl(5, 2)*rx*rz, -(cl(7, 3)*ry**2) + cl(5, 3)*rx*rz, -(cl(7, 4)*ry**2) + &
                    cl(5, 4)*rx*rz, -(cl(7, 5)*ry**2) + cl(5, 5)*rx*rz, -(cl(7, 6)*ry**2) + &
                    cl(5, 6)*rx*rz, -(cl(7, 7)*ry**2) + cl(5, 7)*rx*rz, -(cl(7, 8)*ry**2) + &
                    cl(5, 8)*rx*rz, -(cl(7, 9)*ry**2) + cl(5, 9)*rx*rz]
        m(:, 44) = [-(cl(8, 1)*ry**2) + cl(5, 1)*ry*rz, -(cl(8, 2)*ry**2) + &
                    cl(5, 2)*ry*rz, -(cl(8, 3)*ry**2) + cl(5, 3)*ry*rz, -(cl(8, 4)*ry**2) + &
                    cl(5, 4)*ry*rz, -(cl(8, 5)*ry**2) + cl(5, 5)*ry*rz, -(cl(8, 6)*ry**2) + &
                    cl(5, 6)*ry*rz, -(cl(8, 7)*ry**2) + cl(5, 7)*ry*rz, -(cl(8, 8)*ry**2) + &
                    cl(5, 8)*ry*rz, -(cl(8, 9)*ry**2) + cl(5, 9)*ry*rz]
        m(:, 45) = [-(cl(9, 1)*ry**2) + cl(5, 1)*rz**2, -(cl(9, 2)*ry**2) + &
                    cl(5, 2)*rz**2, -(cl(9, 3)*ry**2) + cl(5, 3)*rz**2, -(cl(9, 4)*ry**2) + &
                    cl(5, 4)*rz**2, -(cl(9, 5)*ry**2) + cl(5, 5)*rz**2, -(cl(9, 6)*ry**2) + &
                    cl(5, 6)*rz**2, -(cl(9, 7)*ry**2) + cl(5, 7)*rz**2, -(cl(9, 8)*ry**2) + &
                    cl(5, 8)*rz**2, -(cl(9, 9)*ry**2) + cl(5, 9)*rz**2]
        m(:, 46) = [cl(6, 1)*rx**2 - cl(1, 1)*ry*rz, cl(6, 2)*rx**2 - cl(1, 2)*ry*rz, &
                    cl(6, 3)*rx**2 - cl(1, 3)*ry*rz, cl(6, 4)*rx**2 - cl(1, 4)*ry*rz, &
                    cl(6, 5)*rx**2 - cl(1, 5)*ry*rz, cl(6, 6)*rx**2 - cl(1, 6)*ry*rz, &
                    cl(6, 7)*rx**2 - cl(1, 7)*ry*rz, cl(6, 8)*rx**2 - cl(1, 8)*ry*rz, &
                    cl(6, 9)*rx**2 - cl(1, 9)*ry*rz]
        m(:, 47) = [cl(6, 1)*rx*ry - cl(2, 1)*ry*rz, cl(6, 2)*rx*ry - &
                    cl(2, 2)*ry*rz, cl(6, 3)*rx*ry - cl(2, 3)*ry*rz, cl(6, 4)*rx*ry - &
                    cl(2, 4)*ry*rz, cl(6, 5)*rx*ry - cl(2, 5)*ry*rz, cl(6, 6)*rx*ry - &
                    cl(2, 6)*ry*rz, cl(6, 7)*rx*ry - cl(2, 7)*ry*rz, cl(6, 8)*rx*ry - &
                    cl(2, 8)*ry*rz, cl(6, 9)*rx*ry - cl(2, 9)*ry*rz]
        m(:, 48) = [cl(6, 1)*rx*rz - cl(3, 1)*ry*rz, cl(6, 2)*rx*rz - &
                    cl(3, 2)*ry*rz, cl(6, 3)*rx*rz - cl(3, 3)*ry*rz, cl(6, 4)*rx*rz - &
                    cl(3, 4)*ry*rz, cl(6, 5)*rx*rz - cl(3, 5)*ry*rz, cl(6, 6)*rx*rz - &
                    cl(3, 6)*ry*rz, cl(6, 7)*rx*rz - cl(3, 7)*ry*rz, cl(6, 8)*rx*rz - &
                    cl(3, 8)*ry*rz, cl(6, 9)*rx*rz - cl(3, 9)*ry*rz]
        m(:, 49) = [cl(6, 1)*rx*ry - cl(4, 1)*ry*rz, cl(6, 2)*rx*ry - &
                    cl(4, 2)*ry*rz, cl(6, 3)*rx*ry - cl(4, 3)*ry*rz, cl(6, 4)*rx*ry - &
                    cl(4, 4)*ry*rz, cl(6, 5)*rx*ry - cl(4, 5)*ry*rz, cl(6, 6)*rx*ry - &
                    cl(4, 6)*ry*rz, cl(6, 7)*rx*ry - cl(4, 7)*ry*rz, cl(6, 8)*rx*ry - &
                    cl(4, 8)*ry*rz, cl(6, 9)*rx*ry - cl(4, 9)*ry*rz]
        m(:, 50) = [cl(6, 1)*ry**2 - cl(5, 1)*ry*rz, cl(6, 2)*ry**2 - cl(5, 2)*ry*rz, &
                    cl(6, 3)*ry**2 - cl(5, 3)*ry*rz, cl(6, 4)*ry**2 - cl(5, 4)*ry*rz, &
                    cl(6, 5)*ry**2 - cl(5, 5)*ry*rz, cl(6, 6)*ry**2 - cl(5, 6)*ry*rz, &
                    cl(6, 7)*ry**2 - cl(5, 7)*ry*rz, cl(6, 8)*ry**2 - cl(5, 8)*ry*rz, &
                    cl(6, 9)*ry**2 - cl(5, 9)*ry*rz]
        m(:, 51) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                    0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 52) = [cl(6, 1)*rx*rz - cl(7, 1)*ry*rz, cl(6, 2)*rx*rz - &
                    cl(7, 2)*ry*rz, cl(6, 3)*rx*rz - cl(7, 3)*ry*rz, cl(6, 4)*rx*rz - &
                    cl(7, 4)*ry*rz, cl(6, 5)*rx*rz - cl(7, 5)*ry*rz, cl(6, 6)*rx*rz - &
                    cl(7, 6)*ry*rz, cl(6, 7)*rx*rz - cl(7, 7)*ry*rz, cl(6, 8)*rx*rz - &
                    cl(7, 8)*ry*rz, cl(6, 9)*rx*rz - cl(7, 9)*ry*rz]
        m(:, 53) = [cl(6, 1)*ry*rz - cl(8, 1)*ry*rz, cl(6, 2)*ry*rz - &
                    cl(8, 2)*ry*rz, cl(6, 3)*ry*rz - cl(8, 3)*ry*rz, cl(6, 4)*ry*rz - &
                    cl(8, 4)*ry*rz, cl(6, 5)*ry*rz - cl(8, 5)*ry*rz, cl(6, 6)*ry*rz - &
                    cl(8, 6)*ry*rz, cl(6, 7)*ry*rz - cl(8, 7)*ry*rz, cl(6, 8)*ry*rz - &
                    cl(8, 8)*ry*rz, cl(6, 9)*ry*rz - cl(8, 9)*ry*rz]
        m(:, 54) = [-(cl(9, 1)*ry*rz) + cl(6, 1)*rz**2, -(cl(9, 2)*ry*rz) + &
                    cl(6, 2)*rz**2, -(cl(9, 3)*ry*rz) + cl(6, 3)*rz**2, -(cl(9, 4)*ry*rz) + &
                    cl(6, 4)*rz**2, -(cl(9, 5)*ry*rz) + cl(6, 5)*rz**2, -(cl(9, 6)*ry*rz) + &
                    cl(6, 6)*rz**2, -(cl(9, 7)*ry*rz) + cl(6, 7)*rz**2, -(cl(9, 8)*ry*rz) + &
                    cl(6, 8)*rz**2, -(cl(9, 9)*ry*rz) + cl(6, 9)*rz**2]
        m(:, 55) = [cl(7, 1)*rx**2 - cl(1, 1)*rx*rz, cl(7, 2)*rx**2 - cl(1, 2)*rx*rz, &
                    cl(7, 3)*rx**2 - cl(1, 3)*rx*rz, cl(7, 4)*rx**2 - cl(1, 4)*rx*rz, &
                    cl(7, 5)*rx**2 - cl(1, 5)*rx*rz, cl(7, 6)*rx**2 - cl(1, 6)*rx*rz, &
                    cl(7, 7)*rx**2 - cl(1, 7)*rx*rz, cl(7, 8)*rx**2 - cl(1, 8)*rx*rz, &
                    cl(7, 9)*rx**2 - cl(1, 9)*rx*rz]
        m(:, 56) = [cl(7, 1)*rx*ry - cl(2, 1)*rx*rz, cl(7, 2)*rx*ry - &
                    cl(2, 2)*rx*rz, cl(7, 3)*rx*ry - cl(2, 3)*rx*rz, cl(7, 4)*rx*ry - &
                    cl(2, 4)*rx*rz, cl(7, 5)*rx*ry - cl(2, 5)*rx*rz, cl(7, 6)*rx*ry - &
                    cl(2, 6)*rx*rz, cl(7, 7)*rx*ry - cl(2, 7)*rx*rz, cl(7, 8)*rx*ry - &
                    cl(2, 8)*rx*rz, cl(7, 9)*rx*ry - cl(2, 9)*rx*rz]
        m(:, 57) = [-(cl(3, 1)*rx*rz) + cl(7, 1)*rx*rz, -(cl(3, 2)*rx*rz) + &
                    cl(7, 2)*rx*rz, -(cl(3, 3)*rx*rz) + cl(7, 3)*rx*rz, -(cl(3, 4)*rx*rz) + &
                    cl(7, 4)*rx*rz, -(cl(3, 5)*rx*rz) + cl(7, 5)*rx*rz, -(cl(3, 6)*rx*rz) + &
                    cl(7, 6)*rx*rz, -(cl(3, 7)*rx*rz) + cl(7, 7)*rx*rz, -(cl(3, 8)*rx*rz) + &
                    cl(7, 8)*rx*rz, -(cl(3, 9)*rx*rz) + cl(7, 9)*rx*rz]
        m(:, 58) = [cl(7, 1)*rx*ry - cl(4, 1)*rx*rz, cl(7, 2)*rx*ry - &
                    cl(4, 2)*rx*rz, cl(7, 3)*rx*ry - cl(4, 3)*rx*rz, cl(7, 4)*rx*ry - &
                    cl(4, 4)*rx*rz, cl(7, 5)*rx*ry - cl(4, 5)*rx*rz, cl(7, 6)*rx*ry - &
                    cl(4, 6)*rx*rz, cl(7, 7)*rx*ry - cl(4, 7)*rx*rz, cl(7, 8)*rx*ry - &
                    cl(4, 8)*rx*rz, cl(7, 9)*rx*ry - cl(4, 9)*rx*rz]
        m(:, 59) = [cl(7, 1)*ry**2 - cl(5, 1)*rx*rz, cl(7, 2)*ry**2 - cl(5, 2)*rx*rz, &
                    cl(7, 3)*ry**2 - cl(5, 3)*rx*rz, cl(7, 4)*ry**2 - cl(5, 4)*rx*rz, &
                    cl(7, 5)*ry**2 - cl(5, 5)*rx*rz, cl(7, 6)*ry**2 - cl(5, 6)*rx*rz, &
                    cl(7, 7)*ry**2 - cl(5, 7)*rx*rz, cl(7, 8)*ry**2 - cl(5, 8)*rx*rz, &
                    cl(7, 9)*ry**2 - cl(5, 9)*rx*rz]
        m(:, 60) = [-(cl(6, 1)*rx*rz) + cl(7, 1)*ry*rz, -(cl(6, 2)*rx*rz) + &
                    cl(7, 2)*ry*rz, -(cl(6, 3)*rx*rz) + cl(7, 3)*ry*rz, -(cl(6, 4)*rx*rz) + &
                    cl(7, 4)*ry*rz, -(cl(6, 5)*rx*rz) + cl(7, 5)*ry*rz, -(cl(6, 6)*rx*rz) + &
                    cl(7, 6)*ry*rz, -(cl(6, 7)*rx*rz) + cl(7, 7)*ry*rz, -(cl(6, 8)*rx*rz) + &
                    cl(7, 8)*ry*rz, -(cl(6, 9)*rx*rz) + cl(7, 9)*ry*rz]
        m(:, 61) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                    0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 62) = [-(cl(8, 1)*rx*rz) + cl(7, 1)*ry*rz, -(cl(8, 2)*rx*rz) + &
                    cl(7, 2)*ry*rz, -(cl(8, 3)*rx*rz) + cl(7, 3)*ry*rz, -(cl(8, 4)*rx*rz) + &
                    cl(7, 4)*ry*rz, -(cl(8, 5)*rx*rz) + cl(7, 5)*ry*rz, -(cl(8, 6)*rx*rz) + &
                    cl(7, 6)*ry*rz, -(cl(8, 7)*rx*rz) + cl(7, 7)*ry*rz, -(cl(8, 8)*rx*rz) + &
                    cl(7, 8)*ry*rz, -(cl(8, 9)*rx*rz) + cl(7, 9)*ry*rz]
        m(:, 63) = [-(cl(9, 1)*rx*rz) + cl(7, 1)*rz**2, -(cl(9, 2)*rx*rz) + &
                    cl(7, 2)*rz**2, -(cl(9, 3)*rx*rz) + cl(7, 3)*rz**2, -(cl(9, 4)*rx*rz) + &
                    cl(7, 4)*rz**2, -(cl(9, 5)*rx*rz) + cl(7, 5)*rz**2, -(cl(9, 6)*rx*rz) + &
                    cl(7, 6)*rz**2, -(cl(9, 7)*rx*rz) + cl(7, 7)*rz**2, -(cl(9, 8)*rx*rz) + &
                    cl(7, 8)*rz**2, -(cl(9, 9)*rx*rz) + cl(7, 9)*rz**2]
        m(:, 64) = [cl(8, 1)*rx**2 - cl(1, 1)*ry*rz, cl(8, 2)*rx**2 - cl(1, 2)*ry*rz, &
                    cl(8, 3)*rx**2 - cl(1, 3)*ry*rz, cl(8, 4)*rx**2 - cl(1, 4)*ry*rz, &
                    cl(8, 5)*rx**2 - cl(1, 5)*ry*rz, cl(8, 6)*rx**2 - cl(1, 6)*ry*rz, &
                    cl(8, 7)*rx**2 - cl(1, 7)*ry*rz, cl(8, 8)*rx**2 - cl(1, 8)*ry*rz, &
                    cl(8, 9)*rx**2 - cl(1, 9)*ry*rz]
        m(:, 65) = [cl(8, 1)*rx*ry - cl(2, 1)*ry*rz, cl(8, 2)*rx*ry - &
                    cl(2, 2)*ry*rz, cl(8, 3)*rx*ry - cl(2, 3)*ry*rz, cl(8, 4)*rx*ry - &
                    cl(2, 4)*ry*rz, cl(8, 5)*rx*ry - cl(2, 5)*ry*rz, cl(8, 6)*rx*ry - &
                    cl(2, 6)*ry*rz, cl(8, 7)*rx*ry - cl(2, 7)*ry*rz, cl(8, 8)*rx*ry - &
                    cl(2, 8)*ry*rz, cl(8, 9)*rx*ry - cl(2, 9)*ry*rz]
        m(:, 66) = [cl(8, 1)*rx*rz - cl(3, 1)*ry*rz, cl(8, 2)*rx*rz - &
                    cl(3, 2)*ry*rz, cl(8, 3)*rx*rz - cl(3, 3)*ry*rz, cl(8, 4)*rx*rz - &
                    cl(3, 4)*ry*rz, cl(8, 5)*rx*rz - cl(3, 5)*ry*rz, cl(8, 6)*rx*rz - &
                    cl(3, 6)*ry*rz, cl(8, 7)*rx*rz - cl(3, 7)*ry*rz, cl(8, 8)*rx*rz - &
                    cl(3, 8)*ry*rz, cl(8, 9)*rx*rz - cl(3, 9)*ry*rz]
        m(:, 67) = [cl(8, 1)*rx*ry - cl(4, 1)*ry*rz, cl(8, 2)*rx*ry - &
                    cl(4, 2)*ry*rz, cl(8, 3)*rx*ry - cl(4, 3)*ry*rz, cl(8, 4)*rx*ry - &
                    cl(4, 4)*ry*rz, cl(8, 5)*rx*ry - cl(4, 5)*ry*rz, cl(8, 6)*rx*ry - &
                    cl(4, 6)*ry*rz, cl(8, 7)*rx*ry - cl(4, 7)*ry*rz, cl(8, 8)*rx*ry - &
                    cl(4, 8)*ry*rz, cl(8, 9)*rx*ry - cl(4, 9)*ry*rz]
        m(:, 68) = [cl(8, 1)*ry**2 - cl(5, 1)*ry*rz, cl(8, 2)*ry**2 - cl(5, 2)*ry*rz, &
                    cl(8, 3)*ry**2 - cl(5, 3)*ry*rz, cl(8, 4)*ry**2 - cl(5, 4)*ry*rz, &
                    cl(8, 5)*ry**2 - cl(5, 5)*ry*rz, cl(8, 6)*ry**2 - cl(5, 6)*ry*rz, &
                    cl(8, 7)*ry**2 - cl(5, 7)*ry*rz, cl(8, 8)*ry**2 - cl(5, 8)*ry*rz, &
                    cl(8, 9)*ry**2 - cl(5, 9)*ry*rz]
        m(:, 69) = [-(cl(6, 1)*ry*rz) + cl(8, 1)*ry*rz, -(cl(6, 2)*ry*rz) + &
                    cl(8, 2)*ry*rz, -(cl(6, 3)*ry*rz) + cl(8, 3)*ry*rz, -(cl(6, 4)*ry*rz) + &
                    cl(8, 4)*ry*rz, -(cl(6, 5)*ry*rz) + cl(8, 5)*ry*rz, -(cl(6, 6)*ry*rz) + &
                    cl(8, 6)*ry*rz, -(cl(6, 7)*ry*rz) + cl(8, 7)*ry*rz, -(cl(6, 8)*ry*rz) + &
                    cl(8, 8)*ry*rz, -(cl(6, 9)*ry*rz) + cl(8, 9)*ry*rz]
        m(:, 70) = [cl(8, 1)*rx*rz - cl(7, 1)*ry*rz, cl(8, 2)*rx*rz - &
                    cl(7, 2)*ry*rz, cl(8, 3)*rx*rz - cl(7, 3)*ry*rz, cl(8, 4)*rx*rz - &
                    cl(7, 4)*ry*rz, cl(8, 5)*rx*rz - cl(7, 5)*ry*rz, cl(8, 6)*rx*rz - &
                    cl(7, 6)*ry*rz, cl(8, 7)*rx*rz - cl(7, 7)*ry*rz, cl(8, 8)*rx*rz - &
                    cl(7, 8)*ry*rz, cl(8, 9)*rx*rz - cl(7, 9)*ry*rz]
        m(:, 71) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                    0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 72) = [-(cl(9, 1)*ry*rz) + cl(8, 1)*rz**2, -(cl(9, 2)*ry*rz) + &
                    cl(8, 2)*rz**2, -(cl(9, 3)*ry*rz) + cl(8, 3)*rz**2, -(cl(9, 4)*ry*rz) + &
                    cl(8, 4)*rz**2, -(cl(9, 5)*ry*rz) + cl(8, 5)*rz**2, -(cl(9, 6)*ry*rz) + &
                    cl(8, 6)*rz**2, -(cl(9, 7)*ry*rz) + cl(8, 7)*rz**2, -(cl(9, 8)*ry*rz) + &
                    cl(8, 8)*rz**2, -(cl(9, 9)*ry*rz) + cl(8, 9)*rz**2]
        m(:, 73) = [cl(9, 1)*rx**2 - cl(1, 1)*rz**2, cl(9, 2)*rx**2 - cl(1, 2)*rz**2, &
                    cl(9, 3)*rx**2 - cl(1, 3)*rz**2, cl(9, 4)*rx**2 - cl(1, 4)*rz**2, &
                    cl(9, 5)*rx**2 - cl(1, 5)*rz**2, cl(9, 6)*rx**2 - cl(1, 6)*rz**2, &
                    cl(9, 7)*rx**2 - cl(1, 7)*rz**2, cl(9, 8)*rx**2 - cl(1, 8)*rz**2, &
                    cl(9, 9)*rx**2 - cl(1, 9)*rz**2]
        m(:, 74) = [cl(9, 1)*rx*ry - cl(2, 1)*rz**2, cl(9, 2)*rx*ry - cl(2, 2)*rz**2, &
                    cl(9, 3)*rx*ry - cl(2, 3)*rz**2, cl(9, 4)*rx*ry - cl(2, 4)*rz**2, &
                    cl(9, 5)*rx*ry - cl(2, 5)*rz**2, cl(9, 6)*rx*ry - cl(2, 6)*rz**2, &
                    cl(9, 7)*rx*ry - cl(2, 7)*rz**2, cl(9, 8)*rx*ry - cl(2, 8)*rz**2, &
                    cl(9, 9)*rx*ry - cl(2, 9)*rz**2]
        m(:, 75) = [cl(9, 1)*rx*rz - cl(3, 1)*rz**2, cl(9, 2)*rx*rz - cl(3, 2)*rz**2, &
                    cl(9, 3)*rx*rz - cl(3, 3)*rz**2, cl(9, 4)*rx*rz - cl(3, 4)*rz**2, &
                    cl(9, 5)*rx*rz - cl(3, 5)*rz**2, cl(9, 6)*rx*rz - cl(3, 6)*rz**2, &
                    cl(9, 7)*rx*rz - cl(3, 7)*rz**2, cl(9, 8)*rx*rz - cl(3, 8)*rz**2, &
                    cl(9, 9)*rx*rz - cl(3, 9)*rz**2]
        m(:, 76) = [cl(9, 1)*rx*ry - cl(4, 1)*rz**2, cl(9, 2)*rx*ry - cl(4, 2)*rz**2, &
                    cl(9, 3)*rx*ry - cl(4, 3)*rz**2, cl(9, 4)*rx*ry - cl(4, 4)*rz**2, &
                    cl(9, 5)*rx*ry - cl(4, 5)*rz**2, cl(9, 6)*rx*ry - cl(4, 6)*rz**2, &
                    cl(9, 7)*rx*ry - cl(4, 7)*rz**2, cl(9, 8)*rx*ry - cl(4, 8)*rz**2, &
                    cl(9, 9)*rx*ry - cl(4, 9)*rz**2]
        m(:, 77) = [cl(9, 1)*ry**2 - cl(5, 1)*rz**2, cl(9, 2)*ry**2 - cl(5, 2)*rz**2, &
                    cl(9, 3)*ry**2 - cl(5, 3)*rz**2, cl(9, 4)*ry**2 - cl(5, 4)*rz**2, &
                    cl(9, 5)*ry**2 - cl(5, 5)*rz**2, cl(9, 6)*ry**2 - cl(5, 6)*rz**2, &
                    cl(9, 7)*ry**2 - cl(5, 7)*rz**2, cl(9, 8)*ry**2 - cl(5, 8)*rz**2, &
                    cl(9, 9)*ry**2 - cl(5, 9)*rz**2]
        m(:, 78) = [cl(9, 1)*ry*rz - cl(6, 1)*rz**2, cl(9, 2)*ry*rz - cl(6, 2)*rz**2, &
                    cl(9, 3)*ry*rz - cl(6, 3)*rz**2, cl(9, 4)*ry*rz - cl(6, 4)*rz**2, &
                    cl(9, 5)*ry*rz - cl(6, 5)*rz**2, cl(9, 6)*ry*rz - cl(6, 6)*rz**2, &
                    cl(9, 7)*ry*rz - cl(6, 7)*rz**2, cl(9, 8)*ry*rz - cl(6, 8)*rz**2, &
                    cl(9, 9)*ry*rz - cl(6, 9)*rz**2]
        m(:, 79) = [cl(9, 1)*rx*rz - cl(7, 1)*rz**2, cl(9, 2)*rx*rz - cl(7, 2)*rz**2, &
                    cl(9, 3)*rx*rz - cl(7, 3)*rz**2, cl(9, 4)*rx*rz - cl(7, 4)*rz**2, &
                    cl(9, 5)*rx*rz - cl(7, 5)*rz**2, cl(9, 6)*rx*rz - cl(7, 6)*rz**2, &
                    cl(9, 7)*rx*rz - cl(7, 7)*rz**2, cl(9, 8)*rx*rz - cl(7, 8)*rz**2, &
                    cl(9, 9)*rx*rz - cl(7, 9)*rz**2]
        m(:, 80) = [cl(9, 1)*ry*rz - cl(8, 1)*rz**2, cl(9, 2)*ry*rz - cl(8, 2)*rz**2, &
                    cl(9, 3)*ry*rz - cl(8, 3)*rz**2, cl(9, 4)*ry*rz - cl(8, 4)*rz**2, &
                    cl(9, 5)*ry*rz - cl(8, 5)*rz**2, cl(9, 6)*ry*rz - cl(8, 6)*rz**2, &
                    cl(9, 7)*ry*rz - cl(8, 7)*rz**2, cl(9, 8)*ry*rz - cl(8, 8)*rz**2, &
                    cl(9, 9)*ry*rz - cl(8, 9)*rz**2]
        m(:, 81) = [0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, 0.0_r8, &
                    0.0_r8, 0.0_r8, 0.0_r8]
    end function
end subroutine

! Below not exposed

!> Translational invariance for the fourth order
subroutine lo_fourthorder_asr(map, eq4, neq4)
    !> symmetry stuff rearranged into coordination shells
    class(lo_forcemap), intent(inout) :: map
    !> fourthorder asr constraints
    real(r8), dimension(:, :), allocatable, intent(out) :: eq4
    !> how many linear constraints
    integer, intent(out) :: neq4

    logical, dimension(81) :: relntheta
    real(r8), dimension(:, :), allocatable :: dumeq1, dumeq2, dumeq3
    ! had to make these allocatable since I ran into lots of stack issues
    real(r8), dimension(:, :), allocatable :: &
        wm81x1, wm81x2, wm81x3, wm81x4, wm81x5, wm81x6, wm81x7, wm81x8, &
        wm81x9, wm81x10, wm81x11, wm81x12, wm81x13, wm81x14, wm81x15, &
        wm81x16, wm81x17, wm81x18, wm81x19, wm81x20, wm81x21, wm81x22, &
        wm81x23, wm81x24, wm81x25, wm81x26, wm81x27, wm81x28, wm81x29, &
        wm81x30, wm81x31, wm81x32, wm81x33, wm81x34, wm81x35, wm81x36, &
        wm81x37, wm81x38, wm81x39, wm81x40, wm81x41, wm81x42, wm81x43, &
        wm81x44, wm81x45, wm81x46, wm81x47, wm81x48, wm81x49, wm81x50, &
        wm81x51, wm81x52, wm81x53, wm81x54, wm81x55, wm81x56, wm81x57, &
        wm81x58, wm81x59, wm81x60, wm81x61, wm81x62, wm81x63, wm81x64, &
        wm81x65, wm81x66, wm81x67, wm81x68, wm81x69, wm81x70, wm81x71, &
        wm81x72, wm81x73, wm81x74, wm81x75, wm81x76, wm81x77, wm81x78, &
        wm81x79, wm81x80, wm81x81
    integer, dimension(81) :: thetaind
    integer, dimension(:), allocatable :: ngroup_per_atom
    integer, dimension(:, :), allocatable :: group_per_atom
    integer :: ngroup_tot, nfc, group, unsh, unop, ntheta, maxntheta
    integer :: a1, i, j, eqctr, ii
    ! which ntheta can happen?
    maxntheta = 0
    relntheta = .false.
    do i = 1, map%n_fc_quartet_shell
        j = map%fc_quartet_shell(i)%nx
        if (j .gt. 0) relntheta(j) = .true.
        maxntheta = max(maxntheta, j)
    end do
    ! total number of groups?
    ngroup_tot = size(map%xuc%fc_quartet_group)
    nfc = map%xuc%nx_fc_quartet
    eqctr = 0

    allocate (ngroup_per_atom(map%n_atom_uc))
    ngroup_per_atom = 0
    do group = 1, size(map%xuc%fc_quartet_group)
        j = map%xuc%fc_quartet_group(group)%ind(1)
        a1 = map%xuc%fc_quartet(j)%i1
        ngroup_per_atom(a1) = ngroup_per_atom(a1) + 1
    end do
    allocate (group_per_atom(maxval(ngroup_per_atom), map%n_atom_uc))
    group_per_atom = -1
    ngroup_per_atom = 0
    do group = 1, size(map%xuc%fc_quartet_group)
        j = map%xuc%fc_quartet_group(group)%ind(1)
        a1 = map%xuc%fc_quartet(j)%i1
        ngroup_per_atom(a1) = ngroup_per_atom(a1) + 1
        group_per_atom(ngroup_per_atom(a1), a1) = group
    end do

    ! Only allocate those that will actually be used, usually quite few
    if (nfc .gt. 0) then
        if (relntheta(1)) allocate (wm81x1(81, 1))
        if (relntheta(2)) allocate (wm81x2(81, 2))
        if (relntheta(3)) allocate (wm81x3(81, 3))
        if (relntheta(4)) allocate (wm81x4(81, 4))
        if (relntheta(5)) allocate (wm81x5(81, 5))
        if (relntheta(6)) allocate (wm81x6(81, 6))
        if (relntheta(7)) allocate (wm81x7(81, 7))
        if (relntheta(8)) allocate (wm81x8(81, 8))
        if (relntheta(9)) allocate (wm81x9(81, 9))
        if (relntheta(10)) allocate (wm81x10(81, 10))
        if (relntheta(11)) allocate (wm81x11(81, 11))
        if (relntheta(12)) allocate (wm81x12(81, 12))
        if (relntheta(13)) allocate (wm81x13(81, 13))
        if (relntheta(14)) allocate (wm81x14(81, 14))
        if (relntheta(15)) allocate (wm81x15(81, 15))
        if (relntheta(16)) allocate (wm81x16(81, 16))
        if (relntheta(17)) allocate (wm81x17(81, 17))
        if (relntheta(18)) allocate (wm81x18(81, 18))
        if (relntheta(19)) allocate (wm81x19(81, 19))
        if (relntheta(20)) allocate (wm81x20(81, 20))
        if (relntheta(21)) allocate (wm81x21(81, 21))
        if (relntheta(22)) allocate (wm81x22(81, 22))
        if (relntheta(23)) allocate (wm81x23(81, 23))
        if (relntheta(24)) allocate (wm81x24(81, 24))
        if (relntheta(25)) allocate (wm81x25(81, 25))
        if (relntheta(26)) allocate (wm81x26(81, 26))
        if (relntheta(27)) allocate (wm81x27(81, 27))
        if (relntheta(28)) allocate (wm81x28(81, 28))
        if (relntheta(29)) allocate (wm81x29(81, 29))
        if (relntheta(30)) allocate (wm81x30(81, 30))
        if (relntheta(31)) allocate (wm81x31(81, 31))
        if (relntheta(32)) allocate (wm81x32(81, 32))
        if (relntheta(33)) allocate (wm81x33(81, 33))
        if (relntheta(34)) allocate (wm81x34(81, 34))
        if (relntheta(35)) allocate (wm81x35(81, 35))
        if (relntheta(36)) allocate (wm81x36(81, 36))
        if (relntheta(37)) allocate (wm81x37(81, 37))
        if (relntheta(38)) allocate (wm81x38(81, 38))
        if (relntheta(39)) allocate (wm81x39(81, 39))
        if (relntheta(40)) allocate (wm81x40(81, 40))
        if (relntheta(41)) allocate (wm81x41(81, 41))
        if (relntheta(42)) allocate (wm81x42(81, 42))
        if (relntheta(43)) allocate (wm81x43(81, 43))
        if (relntheta(44)) allocate (wm81x44(81, 44))
        if (relntheta(45)) allocate (wm81x45(81, 45))
        if (relntheta(46)) allocate (wm81x46(81, 46))
        if (relntheta(47)) allocate (wm81x47(81, 47))
        if (relntheta(48)) allocate (wm81x48(81, 48))
        if (relntheta(49)) allocate (wm81x49(81, 49))
        if (relntheta(50)) allocate (wm81x50(81, 50))
        if (relntheta(51)) allocate (wm81x51(81, 51))
        if (relntheta(52)) allocate (wm81x52(81, 52))
        if (relntheta(53)) allocate (wm81x53(81, 53))
        if (relntheta(54)) allocate (wm81x54(81, 54))
        if (relntheta(55)) allocate (wm81x55(81, 55))
        if (relntheta(56)) allocate (wm81x56(81, 56))
        if (relntheta(57)) allocate (wm81x57(81, 57))
        if (relntheta(58)) allocate (wm81x58(81, 58))
        if (relntheta(59)) allocate (wm81x59(81, 59))
        if (relntheta(60)) allocate (wm81x60(81, 60))
        if (relntheta(61)) allocate (wm81x61(81, 61))
        if (relntheta(62)) allocate (wm81x62(81, 62))
        if (relntheta(63)) allocate (wm81x63(81, 63))
        if (relntheta(64)) allocate (wm81x64(81, 64))
        if (relntheta(65)) allocate (wm81x65(81, 65))
        if (relntheta(66)) allocate (wm81x66(81, 66))
        if (relntheta(67)) allocate (wm81x67(81, 67))
        if (relntheta(68)) allocate (wm81x68(81, 68))
        if (relntheta(69)) allocate (wm81x69(81, 69))
        if (relntheta(70)) allocate (wm81x70(81, 70))
        if (relntheta(71)) allocate (wm81x71(81, 71))
        if (relntheta(72)) allocate (wm81x72(81, 72))
        if (relntheta(73)) allocate (wm81x73(81, 73))
        if (relntheta(74)) allocate (wm81x74(81, 74))
        if (relntheta(75)) allocate (wm81x75(81, 75))
        if (relntheta(76)) allocate (wm81x76(81, 76))
        if (relntheta(77)) allocate (wm81x77(81, 77))
        if (relntheta(78)) allocate (wm81x78(81, 78))
        if (relntheta(79)) allocate (wm81x79(81, 79))
        if (relntheta(80)) allocate (wm81x80(81, 80))
        if (relntheta(81)) allocate (wm81x81(81, 81))
    end if

    lo_allocate(dumeq1(nfc, 81))
    lo_allocate(dumeq2(nfc, 81))
    lo_allocate(dumeq3(nfc, 81*maxntheta))
    dumeq1 = 0.0_r8
    dumeq2 = 0.0_r8
    dumeq3 = 0.0_r8
    atloop: do a1 = 1, map%n_atom_uc
        dumeq2 = 0.0_r8
        grloop: do ii = 1, ngroup_per_atom(a1)
            group = group_per_atom(ii, a1)
            dumeq1 = 0.0_r8
            trloop: do j = 1, map%xuc%fc_quartet_group(group)%n
                i = map%xuc%fc_quartet_group(group)%ind(j)
                unsh = map%xuc%fc_quartet(i)%irreducible_shell
                unop = map%xuc%fc_quartet(i)%operation_from_shell
                ntheta = map%fc_quartet_shell(unsh)%nx
                if (ntheta .eq. 0) cycle trloop
                thetaind(1:ntheta) = map%fc_quartet_shell(unsh)%ind_global
                select case (ntheta)
                case (1)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x1)
                    dumeq1(thetaind(1:1), :) = dumeq1(thetaind(1:1), :) + transpose(wm81x1)
                    dumeq2(thetaind(1:1), :) = dumeq2(thetaind(1:1), :) - transpose(wm81x1)
                case (2)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x2)
                    dumeq1(thetaind(1:2), :) = dumeq1(thetaind(1:2), :) + transpose(wm81x2)
                    dumeq2(thetaind(1:2), :) = dumeq2(thetaind(1:2), :) - transpose(wm81x2)
                case (3)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x3)
                    dumeq1(thetaind(1:3), :) = dumeq1(thetaind(1:3), :) + transpose(wm81x3)
                    dumeq2(thetaind(1:3), :) = dumeq2(thetaind(1:3), :) - transpose(wm81x3)
                case (4)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x4)
                    dumeq1(thetaind(1:4), :) = dumeq1(thetaind(1:4), :) + transpose(wm81x4)
                    dumeq2(thetaind(1:4), :) = dumeq2(thetaind(1:4), :) - transpose(wm81x4)
                case (5)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x5)
                    dumeq1(thetaind(1:5), :) = dumeq1(thetaind(1:5), :) + transpose(wm81x5)
                    dumeq2(thetaind(1:5), :) = dumeq2(thetaind(1:5), :) - transpose(wm81x5)
                case (6)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x6)
                    dumeq1(thetaind(1:6), :) = dumeq1(thetaind(1:6), :) + transpose(wm81x6)
                    dumeq2(thetaind(1:6), :) = dumeq2(thetaind(1:6), :) - transpose(wm81x6)
                case (7)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x7)
                    dumeq1(thetaind(1:7), :) = dumeq1(thetaind(1:7), :) + transpose(wm81x7)
                    dumeq2(thetaind(1:7), :) = dumeq2(thetaind(1:7), :) - transpose(wm81x7)
                case (8)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x8)
                    dumeq1(thetaind(1:8), :) = dumeq1(thetaind(1:8), :) + transpose(wm81x8)
                    dumeq2(thetaind(1:8), :) = dumeq2(thetaind(1:8), :) - transpose(wm81x8)
                case (9)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x9)
                    dumeq1(thetaind(1:9), :) = dumeq1(thetaind(1:9), :) + transpose(wm81x9)
                    dumeq2(thetaind(1:9), :) = dumeq2(thetaind(1:9), :) - transpose(wm81x9)
                case (10)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x10)
                    dumeq1(thetaind(1:10), :) = dumeq1(thetaind(1:10), :) + transpose(wm81x10)
                    dumeq2(thetaind(1:10), :) = dumeq2(thetaind(1:10), :) - transpose(wm81x10)
                case (11)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x11)
                    dumeq1(thetaind(1:11), :) = dumeq1(thetaind(1:11), :) + transpose(wm81x11)
                    dumeq2(thetaind(1:11), :) = dumeq2(thetaind(1:11), :) - transpose(wm81x11)
                case (12)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x12)
                    dumeq1(thetaind(1:12), :) = dumeq1(thetaind(1:12), :) + transpose(wm81x12)
                    dumeq2(thetaind(1:12), :) = dumeq2(thetaind(1:12), :) - transpose(wm81x12)
                case (13)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x13)
                    dumeq1(thetaind(1:13), :) = dumeq1(thetaind(1:13), :) + transpose(wm81x13)
                    dumeq2(thetaind(1:13), :) = dumeq2(thetaind(1:13), :) - transpose(wm81x13)
                case (14)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x14)
                    dumeq1(thetaind(1:14), :) = dumeq1(thetaind(1:14), :) + transpose(wm81x14)
                    dumeq2(thetaind(1:14), :) = dumeq2(thetaind(1:14), :) - transpose(wm81x14)
                case (15)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x15)
                    dumeq1(thetaind(1:15), :) = dumeq1(thetaind(1:15), :) + transpose(wm81x15)
                    dumeq2(thetaind(1:15), :) = dumeq2(thetaind(1:15), :) - transpose(wm81x15)
                case (16)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x16)
                    dumeq1(thetaind(1:16), :) = dumeq1(thetaind(1:16), :) + transpose(wm81x16)
                    dumeq2(thetaind(1:16), :) = dumeq2(thetaind(1:16), :) - transpose(wm81x16)
                case (17)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x17)
                    dumeq1(thetaind(1:17), :) = dumeq1(thetaind(1:17), :) + transpose(wm81x17)
                    dumeq2(thetaind(1:17), :) = dumeq2(thetaind(1:17), :) - transpose(wm81x17)
                case (18)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x18)
                    dumeq1(thetaind(1:18), :) = dumeq1(thetaind(1:18), :) + transpose(wm81x18)
                    dumeq2(thetaind(1:18), :) = dumeq2(thetaind(1:18), :) - transpose(wm81x18)
                case (19)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x19)
                    dumeq1(thetaind(1:19), :) = dumeq1(thetaind(1:19), :) + transpose(wm81x19)
                    dumeq2(thetaind(1:19), :) = dumeq2(thetaind(1:19), :) - transpose(wm81x19)
                case (20)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x20)
                    dumeq1(thetaind(1:20), :) = dumeq1(thetaind(1:20), :) + transpose(wm81x20)
                    dumeq2(thetaind(1:20), :) = dumeq2(thetaind(1:20), :) - transpose(wm81x20)
                case (21)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x21)
                    dumeq1(thetaind(1:21), :) = dumeq1(thetaind(1:21), :) + transpose(wm81x21)
                    dumeq2(thetaind(1:21), :) = dumeq2(thetaind(1:21), :) - transpose(wm81x21)
                case (22)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x22)
                    dumeq1(thetaind(1:22), :) = dumeq1(thetaind(1:22), :) + transpose(wm81x22)
                    dumeq2(thetaind(1:22), :) = dumeq2(thetaind(1:22), :) - transpose(wm81x22)
                case (23)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x23)
                    dumeq1(thetaind(1:23), :) = dumeq1(thetaind(1:23), :) + transpose(wm81x23)
                    dumeq2(thetaind(1:23), :) = dumeq2(thetaind(1:23), :) - transpose(wm81x23)
                case (24)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x24)
                    dumeq1(thetaind(1:24), :) = dumeq1(thetaind(1:24), :) + transpose(wm81x24)
                    dumeq2(thetaind(1:24), :) = dumeq2(thetaind(1:24), :) - transpose(wm81x24)
                case (25)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x25)
                    dumeq1(thetaind(1:25), :) = dumeq1(thetaind(1:25), :) + transpose(wm81x25)
                    dumeq2(thetaind(1:25), :) = dumeq2(thetaind(1:25), :) - transpose(wm81x25)
                case (26)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x26)
                    dumeq1(thetaind(1:26), :) = dumeq1(thetaind(1:26), :) + transpose(wm81x26)
                    dumeq2(thetaind(1:26), :) = dumeq2(thetaind(1:26), :) - transpose(wm81x26)
                case (27)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x27)
                    dumeq1(thetaind(1:27), :) = dumeq1(thetaind(1:27), :) + transpose(wm81x27)
                    dumeq2(thetaind(1:27), :) = dumeq2(thetaind(1:27), :) - transpose(wm81x27)
                case (28)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x28)
                    dumeq1(thetaind(1:28), :) = dumeq1(thetaind(1:28), :) + transpose(wm81x28)
                    dumeq2(thetaind(1:28), :) = dumeq2(thetaind(1:28), :) - transpose(wm81x28)
                case (29)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x29)
                    dumeq1(thetaind(1:29), :) = dumeq1(thetaind(1:29), :) + transpose(wm81x29)
                    dumeq2(thetaind(1:29), :) = dumeq2(thetaind(1:29), :) - transpose(wm81x29)
                case (30)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x30)
                    dumeq1(thetaind(1:30), :) = dumeq1(thetaind(1:30), :) + transpose(wm81x30)
                    dumeq2(thetaind(1:30), :) = dumeq2(thetaind(1:30), :) - transpose(wm81x30)
                case (31)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x31)
                    dumeq1(thetaind(1:31), :) = dumeq1(thetaind(1:31), :) + transpose(wm81x31)
                    dumeq2(thetaind(1:31), :) = dumeq2(thetaind(1:31), :) - transpose(wm81x31)
                case (32)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x32)
                    dumeq1(thetaind(1:32), :) = dumeq1(thetaind(1:32), :) + transpose(wm81x32)
                    dumeq2(thetaind(1:32), :) = dumeq2(thetaind(1:32), :) - transpose(wm81x32)
                case (33)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x33)
                    dumeq1(thetaind(1:33), :) = dumeq1(thetaind(1:33), :) + transpose(wm81x33)
                    dumeq2(thetaind(1:33), :) = dumeq2(thetaind(1:33), :) - transpose(wm81x33)
                case (34)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x34)
                    dumeq1(thetaind(1:34), :) = dumeq1(thetaind(1:34), :) + transpose(wm81x34)
                    dumeq2(thetaind(1:34), :) = dumeq2(thetaind(1:34), :) - transpose(wm81x34)
                case (35)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x35)
                    dumeq1(thetaind(1:35), :) = dumeq1(thetaind(1:35), :) + transpose(wm81x35)
                    dumeq2(thetaind(1:35), :) = dumeq2(thetaind(1:35), :) - transpose(wm81x35)
                case (36)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x36)
                    dumeq1(thetaind(1:36), :) = dumeq1(thetaind(1:36), :) + transpose(wm81x36)
                    dumeq2(thetaind(1:36), :) = dumeq2(thetaind(1:36), :) - transpose(wm81x36)
                case (37)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x37)
                    dumeq1(thetaind(1:37), :) = dumeq1(thetaind(1:37), :) + transpose(wm81x37)
                    dumeq2(thetaind(1:37), :) = dumeq2(thetaind(1:37), :) - transpose(wm81x37)
                case (38)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x38)
                    dumeq1(thetaind(1:38), :) = dumeq1(thetaind(1:38), :) + transpose(wm81x38)
                    dumeq2(thetaind(1:38), :) = dumeq2(thetaind(1:38), :) - transpose(wm81x38)
                case (39)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x39)
                    dumeq1(thetaind(1:39), :) = dumeq1(thetaind(1:39), :) + transpose(wm81x39)
                    dumeq2(thetaind(1:39), :) = dumeq2(thetaind(1:39), :) - transpose(wm81x39)
                case (40)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x40)
                    dumeq1(thetaind(1:40), :) = dumeq1(thetaind(1:40), :) + transpose(wm81x40)
                    dumeq2(thetaind(1:40), :) = dumeq2(thetaind(1:40), :) - transpose(wm81x40)
                case (41)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x41)
                    dumeq1(thetaind(1:41), :) = dumeq1(thetaind(1:41), :) + transpose(wm81x41)
                    dumeq2(thetaind(1:41), :) = dumeq2(thetaind(1:41), :) - transpose(wm81x41)
                case (42)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x42)
                    dumeq1(thetaind(1:42), :) = dumeq1(thetaind(1:42), :) + transpose(wm81x42)
                    dumeq2(thetaind(1:42), :) = dumeq2(thetaind(1:42), :) - transpose(wm81x42)
                case (43)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x43)
                    dumeq1(thetaind(1:43), :) = dumeq1(thetaind(1:43), :) + transpose(wm81x43)
                    dumeq2(thetaind(1:43), :) = dumeq2(thetaind(1:43), :) - transpose(wm81x43)
                case (44)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x44)
                    dumeq1(thetaind(1:44), :) = dumeq1(thetaind(1:44), :) + transpose(wm81x44)
                    dumeq2(thetaind(1:44), :) = dumeq2(thetaind(1:44), :) - transpose(wm81x44)
                case (45)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x45)
                    dumeq1(thetaind(1:45), :) = dumeq1(thetaind(1:45), :) + transpose(wm81x45)
                    dumeq2(thetaind(1:45), :) = dumeq2(thetaind(1:45), :) - transpose(wm81x45)
                case (46)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x46)
                    dumeq1(thetaind(1:46), :) = dumeq1(thetaind(1:46), :) + transpose(wm81x46)
                    dumeq2(thetaind(1:46), :) = dumeq2(thetaind(1:46), :) - transpose(wm81x46)
                case (47)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x47)
                    dumeq1(thetaind(1:47), :) = dumeq1(thetaind(1:47), :) + transpose(wm81x47)
                    dumeq2(thetaind(1:47), :) = dumeq2(thetaind(1:47), :) - transpose(wm81x47)
                case (48)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x48)
                    dumeq1(thetaind(1:48), :) = dumeq1(thetaind(1:48), :) + transpose(wm81x48)
                    dumeq2(thetaind(1:48), :) = dumeq2(thetaind(1:48), :) - transpose(wm81x48)
                case (49)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x49)
                    dumeq1(thetaind(1:49), :) = dumeq1(thetaind(1:49), :) + transpose(wm81x49)
                    dumeq2(thetaind(1:49), :) = dumeq2(thetaind(1:49), :) - transpose(wm81x49)
                case (50)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x50)
                    dumeq1(thetaind(1:50), :) = dumeq1(thetaind(1:50), :) + transpose(wm81x50)
                    dumeq2(thetaind(1:50), :) = dumeq2(thetaind(1:50), :) - transpose(wm81x50)
                case (51)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x51)
                    dumeq1(thetaind(1:51), :) = dumeq1(thetaind(1:51), :) + transpose(wm81x51)
                    dumeq2(thetaind(1:51), :) = dumeq2(thetaind(1:51), :) - transpose(wm81x51)
                case (52)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x52)
                    dumeq1(thetaind(1:52), :) = dumeq1(thetaind(1:52), :) + transpose(wm81x52)
                    dumeq2(thetaind(1:52), :) = dumeq2(thetaind(1:52), :) - transpose(wm81x52)
                case (53)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x53)
                    dumeq1(thetaind(1:53), :) = dumeq1(thetaind(1:53), :) + transpose(wm81x53)
                    dumeq2(thetaind(1:53), :) = dumeq2(thetaind(1:53), :) - transpose(wm81x53)
                case (54)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x54)
                    dumeq1(thetaind(1:54), :) = dumeq1(thetaind(1:54), :) + transpose(wm81x54)
                    dumeq2(thetaind(1:54), :) = dumeq2(thetaind(1:54), :) - transpose(wm81x54)
                case (55)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x55)
                    dumeq1(thetaind(1:55), :) = dumeq1(thetaind(1:55), :) + transpose(wm81x55)
                    dumeq2(thetaind(1:55), :) = dumeq2(thetaind(1:55), :) - transpose(wm81x55)
                case (56)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x56)
                    dumeq1(thetaind(1:56), :) = dumeq1(thetaind(1:56), :) + transpose(wm81x56)
                    dumeq2(thetaind(1:56), :) = dumeq2(thetaind(1:56), :) - transpose(wm81x56)
                case (57)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x57)
                    dumeq1(thetaind(1:57), :) = dumeq1(thetaind(1:57), :) + transpose(wm81x57)
                    dumeq2(thetaind(1:57), :) = dumeq2(thetaind(1:57), :) - transpose(wm81x57)
                case (58)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x58)
                    dumeq1(thetaind(1:58), :) = dumeq1(thetaind(1:58), :) + transpose(wm81x58)
                    dumeq2(thetaind(1:58), :) = dumeq2(thetaind(1:58), :) - transpose(wm81x58)
                case (59)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x59)
                    dumeq1(thetaind(1:59), :) = dumeq1(thetaind(1:59), :) + transpose(wm81x59)
                    dumeq2(thetaind(1:59), :) = dumeq2(thetaind(1:59), :) - transpose(wm81x59)
                case (60)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x60)
                    dumeq1(thetaind(1:60), :) = dumeq1(thetaind(1:60), :) + transpose(wm81x60)
                    dumeq2(thetaind(1:60), :) = dumeq2(thetaind(1:60), :) - transpose(wm81x60)
                case (61)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x61)
                    dumeq1(thetaind(1:61), :) = dumeq1(thetaind(1:61), :) + transpose(wm81x61)
                    dumeq2(thetaind(1:61), :) = dumeq2(thetaind(1:61), :) - transpose(wm81x61)
                case (62)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x62)
                    dumeq1(thetaind(1:62), :) = dumeq1(thetaind(1:62), :) + transpose(wm81x62)
                    dumeq2(thetaind(1:62), :) = dumeq2(thetaind(1:62), :) - transpose(wm81x62)
                case (63)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x63)
                    dumeq1(thetaind(1:63), :) = dumeq1(thetaind(1:63), :) + transpose(wm81x63)
                    dumeq2(thetaind(1:63), :) = dumeq2(thetaind(1:63), :) - transpose(wm81x63)
                case (64)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x64)
                    dumeq1(thetaind(1:64), :) = dumeq1(thetaind(1:64), :) + transpose(wm81x64)
                    dumeq2(thetaind(1:64), :) = dumeq2(thetaind(1:64), :) - transpose(wm81x64)
                case (65)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x65)
                    dumeq1(thetaind(1:65), :) = dumeq1(thetaind(1:65), :) + transpose(wm81x65)
                    dumeq2(thetaind(1:65), :) = dumeq2(thetaind(1:65), :) - transpose(wm81x65)
                case (66)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x66)
                    dumeq1(thetaind(1:66), :) = dumeq1(thetaind(1:66), :) + transpose(wm81x66)
                    dumeq2(thetaind(1:66), :) = dumeq2(thetaind(1:66), :) - transpose(wm81x66)
                case (67)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x67)
                    dumeq1(thetaind(1:67), :) = dumeq1(thetaind(1:67), :) + transpose(wm81x67)
                    dumeq2(thetaind(1:67), :) = dumeq2(thetaind(1:67), :) - transpose(wm81x67)
                case (68)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x68)
                    dumeq1(thetaind(1:68), :) = dumeq1(thetaind(1:68), :) + transpose(wm81x68)
                    dumeq2(thetaind(1:68), :) = dumeq2(thetaind(1:68), :) - transpose(wm81x68)
                case (69)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x69)
                    dumeq1(thetaind(1:69), :) = dumeq1(thetaind(1:69), :) + transpose(wm81x69)
                    dumeq2(thetaind(1:69), :) = dumeq2(thetaind(1:69), :) - transpose(wm81x69)
                case (70)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x70)
                    dumeq1(thetaind(1:70), :) = dumeq1(thetaind(1:70), :) + transpose(wm81x70)
                    dumeq2(thetaind(1:70), :) = dumeq2(thetaind(1:70), :) - transpose(wm81x70)
                case (71)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x71)
                    dumeq1(thetaind(1:71), :) = dumeq1(thetaind(1:71), :) + transpose(wm81x71)
                    dumeq2(thetaind(1:71), :) = dumeq2(thetaind(1:71), :) - transpose(wm81x71)
                case (72)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x72)
                    dumeq1(thetaind(1:72), :) = dumeq1(thetaind(1:72), :) + transpose(wm81x72)
                    dumeq2(thetaind(1:72), :) = dumeq2(thetaind(1:72), :) - transpose(wm81x72)
                case (73)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x73)
                    dumeq1(thetaind(1:73), :) = dumeq1(thetaind(1:73), :) + transpose(wm81x73)
                    dumeq2(thetaind(1:73), :) = dumeq2(thetaind(1:73), :) - transpose(wm81x73)
                case (74)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x74)
                    dumeq1(thetaind(1:74), :) = dumeq1(thetaind(1:74), :) + transpose(wm81x74)
                    dumeq2(thetaind(1:74), :) = dumeq2(thetaind(1:74), :) - transpose(wm81x74)
                case (75)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x75)
                    dumeq1(thetaind(1:75), :) = dumeq1(thetaind(1:75), :) + transpose(wm81x75)
                    dumeq2(thetaind(1:75), :) = dumeq2(thetaind(1:75), :) - transpose(wm81x75)
                case (76)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x76)
                    dumeq1(thetaind(1:76), :) = dumeq1(thetaind(1:76), :) + transpose(wm81x76)
                    dumeq2(thetaind(1:76), :) = dumeq2(thetaind(1:76), :) - transpose(wm81x76)
                case (77)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x77)
                    dumeq1(thetaind(1:77), :) = dumeq1(thetaind(1:77), :) + transpose(wm81x77)
                    dumeq2(thetaind(1:77), :) = dumeq2(thetaind(1:77), :) - transpose(wm81x77)
                case (78)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x78)
                    dumeq1(thetaind(1:78), :) = dumeq1(thetaind(1:78), :) + transpose(wm81x78)
                    dumeq2(thetaind(1:78), :) = dumeq2(thetaind(1:78), :) - transpose(wm81x78)
                case (79)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x79)
                    dumeq1(thetaind(1:79), :) = dumeq1(thetaind(1:79), :) + transpose(wm81x79)
                    dumeq2(thetaind(1:79), :) = dumeq2(thetaind(1:79), :) - transpose(wm81x79)
                case (80)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x80)
                    dumeq1(thetaind(1:80), :) = dumeq1(thetaind(1:80), :) + transpose(wm81x80)
                    dumeq2(thetaind(1:80), :) = dumeq2(thetaind(1:80), :) - transpose(wm81x80)
                case (81)
                    call lo_gemm(map%op_quartet(unop)%sotr, map%fc_quartet_shell(unsh)%coeff, wm81x81)
                    dumeq1(thetaind(1:81), :) = dumeq1(thetaind(1:81), :) + transpose(wm81x81)
                    dumeq2(thetaind(1:81), :) = dumeq2(thetaind(1:81), :) - transpose(wm81x81)
                end select
            end do trloop
            ! What to do with these equations?
            if (map%xuc%fc_quartet_group(group)%contains_selfterm) then
                ! if this was the group with the self-term, accumulate
                dumeq2 = dumeq2 + dumeq1
            else
                ! this was not including the self-term, store these
                eql_asr4_1: do i = 1, 81
                    if (sum(abs(dumeq1(:, i))) .lt. lo_tol) cycle eql_asr4_1
                    do j = 1, eqctr
                        if (sum(abs(dumeq3(:, j) - dumeq1(:, i))) .lt. lo_tol) cycle eql_asr4_1
                    end do
                    eqctr = eqctr + 1
                    dumeq3(:, eqctr) = dumeq1(:, i)
                end do eql_asr4_1
            end if
        end do grloop
        ! store the self-term group
        eql_asr4_2: do i = 1, 81
            if (sum(abs(dumeq2(:, i))) .lt. lo_tol) cycle eql_asr4_2
            do j = 1, eqctr
                if (sum(abs(dumeq3(:, j) - dumeq2(:, i))) .lt. lo_tol) cycle eql_asr4_2
            end do
            eqctr = eqctr + 1
            dumeq3(:, eqctr) = dumeq2(:, i)
        end do eql_asr4_2
    end do atloop

    ! reduce and return equations
    call lo_reduce_equations(dumeq3(:, 1:eqctr), eq4, neq4)
end subroutine

!> Translational symmetry for the third order
subroutine lo_thirdorder_asr(map, eq3, neq3)
    !> symmetry stuff rearranged into coordination shells
    class(lo_forcemap), intent(inout) :: map
    !> thirdorder asr constraints
    real(r8), dimension(:, :), allocatable, intent(out) :: eq3
    !> how many linear constraints
    integer, intent(out) :: neq3

    real(r8), dimension(:, :), allocatable :: dumeq1, dumeq2, dumeq3
    real(r8) :: wm1x27(1, 27), wm2x27(2, 27), wm3x27(3, 27), wm4x27(4, 27), wm5x27(5, 27), wm6x27(6, 27), wm7x27(7, 27), &
                wm8x27(8, 27), wm9x27(9, 27), wm10x27(10, 27), wm11x27(11, 27), wm12x27(12, 27), wm13x27(13, 27), &
                wm14x27(14, 27), wm15x27(15, 27), wm16x27(16, 27), wm17x27(17, 27), wm18x27(18, 27), wm19x27(19, 27), &
                wm20x27(20, 27), wm21x27(21, 27), wm22x27(22, 27), wm23x27(23, 27), wm24x27(24, 27), wm25x27(25, 27), &
                wm26x27(26, 27), wm27x27(27, 27)
    integer, dimension(27) :: xind
    integer :: ngroup_tot, nfc, igroup, itrip, unsh, unop, nx
    integer :: a1, i, j, eqctr

    ! Get total number of groups
    ngroup_tot = size(map%xuc%fc_triplet_group)
    nfc = map%xuc%nx_fc_triplet
    eqctr = 0
    allocate (dumeq1(nfc, 27))
    allocate (dumeq2(nfc, 27))
    allocate (dumeq3(nfc, 27*ngroup_tot)) ! definitely safe. Could probably reduce this a little, but ok for now.
    dumeq1 = 0.0_r8
    dumeq2 = 0.0_r8
    dumeq3 = 0.0_r8
    atomloop: do a1 = 1, map%n_atom_uc
        dumeq2 = 0.0_r8
        grouploop: do igroup = 1, ngroup_tot
            dumeq1 = 0.0_r8
            triploop: do i = 1, map%xuc%fc_triplet_group(igroup)%n
                itrip = map%xuc%fc_triplet_group(igroup)%ind(i)
                if (map%xuc%fc_triplet(itrip)%i1 .ne. a1) cycle triploop ! maybe cycle grouploop
                unsh = map%xuc%fc_triplet(itrip)%irreducible_shell
                unop = map%xuc%fc_triplet(itrip)%operation_from_shell
                nx = map%fc_triplet_shell(unsh)%nx
                if (nx .eq. 0) cycle
                xind(1:nx) = map%fc_triplet_shell(unsh)%ind_global
                select case (nx)
                case (1)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm1x27, transa='T', transb='T')
                    dumeq1(xind(1:1), :) = dumeq1(xind(1:1), :) + wm1x27
                    dumeq2(xind(1:1), :) = dumeq2(xind(1:1), :) - wm1x27
                case (2)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm2x27, transa='T', transb='T')
                    dumeq1(xind(1:2), :) = dumeq1(xind(1:2), :) + wm2x27
                    dumeq2(xind(1:2), :) = dumeq2(xind(1:2), :) - wm2x27
                case (3)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm3x27, transa='T', transb='T')
                    dumeq1(xind(1:3), :) = dumeq1(xind(1:3), :) + wm3x27
                    dumeq2(xind(1:3), :) = dumeq2(xind(1:3), :) - wm3x27
                case (4)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm4x27, transa='T', transb='T')
                    dumeq1(xind(1:4), :) = dumeq1(xind(1:4), :) + wm4x27
                    dumeq2(xind(1:4), :) = dumeq2(xind(1:4), :) - wm4x27
                case (5)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm5x27, transa='T', transb='T')
                    dumeq1(xind(1:5), :) = dumeq1(xind(1:5), :) + wm5x27
                    dumeq2(xind(1:5), :) = dumeq2(xind(1:5), :) - wm5x27
                case (6)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm6x27, transa='T', transb='T')
                    dumeq1(xind(1:6), :) = dumeq1(xind(1:6), :) + wm6x27
                    dumeq2(xind(1:6), :) = dumeq2(xind(1:6), :) - wm6x27
                case (7)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm7x27, transa='T', transb='T')
                    dumeq1(xind(1:7), :) = dumeq1(xind(1:7), :) + wm7x27
                    dumeq2(xind(1:7), :) = dumeq2(xind(1:7), :) - wm7x27
                case (8)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm8x27, transa='T', transb='T')
                    dumeq1(xind(1:8), :) = dumeq1(xind(1:8), :) + wm8x27
                    dumeq2(xind(1:8), :) = dumeq2(xind(1:8), :) - wm8x27
                case (9)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm9x27, transa='T', transb='T')
                    dumeq1(xind(1:9), :) = dumeq1(xind(1:9), :) + wm9x27
                    dumeq2(xind(1:9), :) = dumeq2(xind(1:9), :) - wm9x27
                case (10)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm10x27, transa='T', transb='T')
                    dumeq1(xind(1:10), :) = dumeq1(xind(1:10), :) + wm10x27
                    dumeq2(xind(1:10), :) = dumeq2(xind(1:10), :) - wm10x27
                case (11)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm11x27, transa='T', transb='T')
                    dumeq1(xind(1:11), :) = dumeq1(xind(1:11), :) + wm11x27
                    dumeq2(xind(1:11), :) = dumeq2(xind(1:11), :) - wm11x27
                case (12)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm12x27, transa='T', transb='T')
                    dumeq1(xind(1:12), :) = dumeq1(xind(1:12), :) + wm12x27
                    dumeq2(xind(1:12), :) = dumeq2(xind(1:12), :) - wm12x27
                case (13)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm13x27, transa='T', transb='T')
                    dumeq1(xind(1:13), :) = dumeq1(xind(1:13), :) + wm13x27
                    dumeq2(xind(1:13), :) = dumeq2(xind(1:13), :) - wm13x27
                case (14)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm14x27, transa='T', transb='T')
                    dumeq1(xind(1:14), :) = dumeq1(xind(1:14), :) + wm14x27
                    dumeq2(xind(1:14), :) = dumeq2(xind(1:14), :) - wm14x27
                case (15)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm15x27, transa='T', transb='T')
                    dumeq1(xind(1:15), :) = dumeq1(xind(1:15), :) + wm15x27
                    dumeq2(xind(1:15), :) = dumeq2(xind(1:15), :) - wm15x27
                case (16)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm16x27, transa='T', transb='T')
                    dumeq1(xind(1:16), :) = dumeq1(xind(1:16), :) + wm16x27
                    dumeq2(xind(1:16), :) = dumeq2(xind(1:16), :) - wm16x27
                case (17)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm17x27, transa='T', transb='T')
                    dumeq1(xind(1:17), :) = dumeq1(xind(1:17), :) + wm17x27
                    dumeq2(xind(1:17), :) = dumeq2(xind(1:17), :) - wm17x27
                case (18)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm18x27, transa='T', transb='T')
                    dumeq1(xind(1:18), :) = dumeq1(xind(1:18), :) + wm18x27
                    dumeq2(xind(1:18), :) = dumeq2(xind(1:18), :) - wm18x27
                case (19)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm19x27, transa='T', transb='T')
                    dumeq1(xind(1:19), :) = dumeq1(xind(1:19), :) + wm19x27
                    dumeq2(xind(1:19), :) = dumeq2(xind(1:19), :) - wm19x27
                case (20)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm20x27, transa='T', transb='T')
                    dumeq1(xind(1:20), :) = dumeq1(xind(1:20), :) + wm20x27
                    dumeq2(xind(1:20), :) = dumeq2(xind(1:20), :) - wm20x27
                case (21)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm21x27, transa='T', transb='T')
                    dumeq1(xind(1:21), :) = dumeq1(xind(1:21), :) + wm21x27
                    dumeq2(xind(1:21), :) = dumeq2(xind(1:21), :) - wm21x27
                case (22)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm22x27, transa='T', transb='T')
                    dumeq1(xind(1:22), :) = dumeq1(xind(1:22), :) + wm22x27
                    dumeq2(xind(1:22), :) = dumeq2(xind(1:22), :) - wm22x27
                case (23)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm23x27, transa='T', transb='T')
                    dumeq1(xind(1:23), :) = dumeq1(xind(1:23), :) + wm23x27
                    dumeq2(xind(1:23), :) = dumeq2(xind(1:23), :) - wm23x27
                case (24)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm24x27, transa='T', transb='T')
                    dumeq1(xind(1:24), :) = dumeq1(xind(1:24), :) + wm24x27
                    dumeq2(xind(1:24), :) = dumeq2(xind(1:24), :) - wm24x27
                case (25)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm25x27, transa='T', transb='T')
                    dumeq1(xind(1:25), :) = dumeq1(xind(1:25), :) + wm25x27
                    dumeq2(xind(1:25), :) = dumeq2(xind(1:25), :) - wm25x27
                case (26)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm26x27, transa='T', transb='T')
                    dumeq1(xind(1:26), :) = dumeq1(xind(1:26), :) + wm26x27
                    dumeq2(xind(1:26), :) = dumeq2(xind(1:26), :) - wm26x27
                case (27)
                    call lo_gemm(map%fc_triplet_shell(unsh)%coeff, map%op_triplet(unop)%sotr, wm27x27, transa='T', transb='T')
                    dumeq1(xind(1:27), :) = dumeq1(xind(1:27), :) + wm27x27
                    dumeq2(xind(1:27), :) = dumeq2(xind(1:27), :) - wm27x27
                end select
            end do triploop
            ! What to do with these equations?
            if (map%xuc%fc_triplet_group(igroup)%contains_selfterm) then
                ! if this was the group with the self-term, accumulate
                dumeq2 = dumeq2 + dumeq1
            else
                ! this was not including the self-term, store these
                eql_asr3_1: do i = 1, 27
                    if (sum(abs(dumeq1(:, i))) .lt. lo_tol) cycle eql_asr3_1
                    do j = 1, eqctr
                        if (sum(abs(dumeq3(:, j) - dumeq1(:, i))) .lt. lo_tol) cycle eql_asr3_1
                    end do
                    eqctr = eqctr + 1
                    dumeq3(:, eqctr) = dumeq1(:, i)
                end do eql_asr3_1
            end if
        end do grouploop
        eql_asr3_2: do i = 1, 27
            if (sum(abs(dumeq2(:, i))) .lt. lo_tol) cycle eql_asr3_2
            do j = 1, eqctr
                if (sum(abs(dumeq3(:, j) - dumeq2(:, i))) .lt. lo_tol) cycle eql_asr3_2
            end do
            eqctr = eqctr + 1
            dumeq3(:, eqctr) = dumeq2(:, i)
        end do eql_asr3_2
    end do atomloop
    ! And clean the equations
    call lo_reduce_equations(dumeq3(:, 1:eqctr), eq3, neq3)
end subroutine

!> Translational symmetry for the third order
subroutine lo_Ztriplet_asr(map, eq3, neq3)
    !> symmetry stuff rearranged into coordination shells
    class(lo_forcemap), intent(inout) :: map
    !> thirdorder asr constraints
    real(r8), dimension(:, :), allocatable, intent(out) :: eq3
    !> how many linear constraints
    integer, intent(out) :: neq3

    real(r8), dimension(:, :), allocatable :: dumeq1, dumeq2, dumeq3
    real(r8), dimension(81, 81) :: B, C
    integer, dimension(81) :: xind
    integer :: ngroup_tot, nfc, igroup, itrip, unsh, unop, nx
    integer :: a1, i, j, eqctr

    ! Get total number of groups
    ngroup_tot = size(map%xuc%Z_triplet_group)
    nfc = map%xuc%nx_Z_triplet
    eqctr = 0

    allocate (dumeq1(nfc, 81))
    allocate (dumeq2(nfc, 81))
    allocate (dumeq3(nfc, 81*ngroup_tot)) ! definitely safe. Could probably reduce this a little, but ok for now.
    dumeq1 = 0.0_r8
    dumeq2 = 0.0_r8
    dumeq3 = 0.0_r8
    atomloop: do a1 = 1, map%n_atom_uc
        dumeq2 = 0.0_r8
        grouploop: do igroup = 1, ngroup_tot
            dumeq1 = 0.0_r8
            triploop: do i = 1, map%xuc%Z_triplet_group(igroup)%n
                itrip = map%xuc%Z_triplet_group(igroup)%ind(i)
                if (map%xuc%Z_triplet(itrip)%i1 .ne. a1) cycle triploop ! maybe cycle grouploop
                unsh = map%xuc%Z_triplet(itrip)%irreducible_shell
                unop = map%xuc%Z_triplet(itrip)%operation_from_shell
                nx = map%Z_triplet_shell(unsh)%nx
                if (nx .eq. 0) cycle
                xind(1:nx) = map%Z_triplet_shell(unsh)%ind_global
                call lo_gemm(map%op_quartet(unop)%sotr, map%Z_triplet_shell(unsh)%coeff, B(:, 1:nx))
                C(1:nx, :) = transpose(B(:, 1:nx))
                dumeq1(xind(1:nx), :) = dumeq1(xind(1:nx), :) + C(1:nx, :) !wm1x27
                dumeq2(xind(1:nx), :) = dumeq2(xind(1:nx), :) - C(1:nx, :) !wm1x27
            end do triploop
            ! What to do with these equations?
            if (map%xuc%Z_triplet_group(igroup)%contains_selfterm) then
                ! if this was the group with the self-term, accumulate
                dumeq2 = dumeq2 + dumeq1
            else
                ! this was not including the self-term, store these
                eql_asr3_1: do i = 1, 81
                    if (sum(abs(dumeq1(:, i))) .lt. lo_tol) cycle eql_asr3_1
                    do j = 1, eqctr
                        if (sum(abs(dumeq3(:, j) - dumeq1(:, i))) .lt. lo_tol) cycle eql_asr3_1
                    end do
                    eqctr = eqctr + 1
                    dumeq3(:, eqctr) = dumeq1(:, i)
                end do eql_asr3_1
            end if
        end do grouploop
        eql_asr3_2: do i = 1, 81
            if (sum(abs(dumeq2(:, i))) .lt. lo_tol) cycle eql_asr3_2
            do j = 1, eqctr
                if (sum(abs(dumeq3(:, j) - dumeq2(:, i))) .lt. lo_tol) cycle eql_asr3_2
            end do
            eqctr = eqctr + 1
            dumeq3(:, eqctr) = dumeq2(:, i)
        end do eql_asr3_2
    end do atomloop
    ! And clean the equations
    call lo_reduce_equations(dumeq3(:, 1:eqctr), eq3, neq3)
end subroutine

!> Rotational constraints for the first order forceconstants, as well as constraints to make sure it adds up to zero
subroutine lo_firstorder_rotational_translational(map, uc, firstorder_rot, neq, rotational)
    !> symmetry stuff rearranged into coordination shells
    class(lo_forcemap), intent(inout) :: map
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> the constraints
    real(r8), dimension(:, :), allocatable, intent(out) :: firstorder_rot
    !> how many equations did I get?
    integer, intent(out) :: neq
    !> rotational constraints?
    logical, intent(in) :: rotational

    real(r8), dimension(:, :), allocatable :: dumeq
    real(r8) :: wm3x1(3, 1), wm3x2(3, 2), wm3x3(3, 3)
    real(r8), dimension(3, 9) :: m3x9
    real(r8), dimension(3, 3) :: CM, CMT
    integer, dimension(3) :: thetaind, thetalocind
    integer :: a1, nfc, unsh, unop, ntheta

    nfc = map%xuc%nx_fc_singlet
    lo_allocate(dumeq(nfc, 9 + 3)) ! max n of equations
    dumeq = 0.0_r8
    do a1 = 1, map%n_atom_uc
        ! Grab the indices where to do stuff, and the
        unsh = map%xuc%fc_singlet(a1)%irreducible_shell
        unop = map%xuc%fc_singlet(a1)%operation_from_shell
        ntheta = map%fc_singlet_shell(unsh)%nx
        if (ntheta .eq. 0) cycle
        thetaind(1:ntheta) = map%fc_singlet_shell(unsh)%ind_global
        thetalocind(1:ntheta) = map%fc_singlet_shell(unsh)%ind_local
        CM = 0.0_r8
        select case (ntheta)
        case (1)
            wm3x1 = matmul(map%op_singlet(unop)%m3, map%fc_singlet_shell(unsh)%coeff)
            CM(:, thetalocind(1:ntheta)) = wm3x1
        case (2)
            wm3x2 = matmul(map%op_singlet(unop)%m3, map%fc_singlet_shell(unsh)%coeff)
            CM(:, thetalocind(1:ntheta)) = wm3x2
        case (3)
            wm3x3 = matmul(map%op_singlet(unop)%m3, map%fc_singlet_shell(unsh)%coeff)
            CM(:, thetalocind(1:ntheta)) = wm3x3
        case default
            call lo_stop_gracefully(['many constraint things for singlets'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end select
        ! Rotations first
        if (rotational) then
            m3x9 = firstorder_rotational_coefficient(uc%rcart(:, a1), CM)
            dumeq(thetaind(1:ntheta), 1:9) = dumeq(thetaind(1:ntheta), 1:9) + m3x9(thetalocind(1:ntheta), :)
        end if
        ! And acoustic sum rule
        CMT = transpose(CM)
        dumeq(thetaind(1:ntheta), 10:12) = dumeq(thetaind(1:ntheta), 10:12) + CMT(thetalocind(1:ntheta), :)
    end do
    ! And reduce and return these
    call lo_reduce_equations(dumeq, firstorder_rot, neq)

contains
    pure function firstorder_rotational_coefficient(r, ck) result(m)
        real(r8), dimension(3), intent(in) :: r
        real(r8), dimension(3, 3), intent(in) :: ck
        real(r8), dimension(3, 9) :: m
        !
        real(r8) :: rx, ry, rz
        rx = r(1)
        ry = r(2)
        rz = r(3)
        m(:, 1) = [0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 2) = [-(ck(2, 1)*rx) + ck(1, 1)*ry, -(ck(2, 2)*rx) + ck(1, 2)*ry, -(ck(2, 3)*rx) + ck(1, 3)*ry]
        m(:, 3) = [-(ck(3, 1)*rx) + ck(1, 1)*rz, -(ck(3, 2)*rx) + ck(1, 2)*rz, -(ck(3, 3)*rx) + ck(1, 3)*rz]
        m(:, 4) = [ck(2, 1)*rx - ck(1, 1)*ry, ck(2, 2)*rx - ck(1, 2)*ry, ck(2, 3)*rx - ck(1, 3)*ry]
        m(:, 5) = [0.0_r8, 0.0_r8, 0.0_r8]
        m(:, 6) = [-(ck(3, 1)*ry) + ck(2, 1)*rz, -(ck(3, 2)*ry) + ck(2, 2)*rz, -(ck(3, 3)*ry) + ck(2, 3)*rz]
        m(:, 7) = [ck(3, 1)*rx - ck(1, 1)*rz, ck(3, 2)*rx - ck(1, 2)*rz, ck(3, 3)*rx - ck(1, 3)*rz]
        m(:, 8) = [ck(3, 1)*ry - ck(2, 1)*rz, ck(3, 2)*ry - ck(2, 2)*rz, ck(3, 3)*ry - ck(2, 3)*rz]
        m(:, 9) = [0.0_r8, 0.0_r8, 0.0_r8]
    end function
end subroutine

!> take a bunch of equations and SVD them, to get the irreducible amount
subroutine lo_reduce_equations(allequations, redeq, nredeq)
    !> all equations
    real(r8), dimension(:, :), intent(in) :: allequations
    !> the reduced set
    real(r8), dimension(:, :), allocatable, intent(out) :: redeq
    !> how many in the reduced set
    integer, intent(out) :: nredeq

    real(r8), dimension(:, :), allocatable :: m, u, v
    real(r8), dimension(:), allocatable :: s
    integer :: i, l, nu, ne, ns
    !
    nu = size(allequations, 1)
    ne = size(allequations, 2)
    ns = min(nu, ne)
    allocate (m(nu, ne))
    m = allequations !transpose(allequations)
    ! remove tiny tiny numbers
    m = lo_chop(m, lo_sqtol)
    ! SVD this
    allocate(s(ns))
    allocate(u(nu,nu))
    allocate(v(ne,ne))
    call lo_dgesvd(m,s,u,v)
    l=0
    do i=1,ns
        if ( s(i) .gt. lo_tol ) then
            l=l+1
        endif
    enddo
    nredeq=l
    if ( l .gt. 0 ) then
        ! return the equations with non-zero singular vectors
        allocate (redeq(l, nu))
        redeq = transpose(u(:, 1:l))
        ! remove tiny numbers
        redeq = lo_chop(redeq, lo_sqtol)
    end if
    ! cleanup
    deallocate (m, s, u, v)
end subroutine

end submodule
