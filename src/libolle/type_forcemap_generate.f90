#include "precompilerdefinitions"
submodule(type_forcemap) type_forcemap_generate
use gottochblandat, only: lo_clean_fractional_coordinates
use lo_sorting, only: lo_return_unique_indices
implicit none
contains

!> get the map of how forces are related to specific IFCs
module subroutine generate(map, uc, ss, polarcorrectiontype, st, mw, mem, verbosity, devmode)
    !> forcemap
    class(lo_forcemap), intent(out) :: map
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> what kind of polar correction
    integer, intent(in) :: polarcorrectiontype
    !> symmetry list thing
    type(lo_interaction_tensors), intent(inout) :: st
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity
    !> developer mode?
    logical, intent(in), optional :: devmode

    real(r8) :: t0, t1

    call mem%tick()
    t0 = walltime()
    t1 = t0

    ! fkdev: dev mode?
    if (present(devmode)) then
        map%devmode = devmode
    else
        map%devmode = .false.
    end if

    init: block
        if (verbosity .gt. 0) then
            t0 = walltime()
            write (lo_iou, *) ''
            write (lo_iou, *) 'GENERATING FORCEMAP'
        end if

        ! In general, I will want to know how many atoms there are.
        map%n_atom_uc = uc%na
        map%n_atom_ss = ss%na

        ! Will also likely need symmetry operations in the odd tensorial thing
        map%op_singlet = st%singletop
        map%op_pair = st%pairop
        map%op_triplet = st%tripletop
        map%op_quartet = st%quartetop

        ! Check if it's polar, and if so how much.
        map%polar = 0
        if (st%have_eps_global .and. st%have_Z_singlet) then
            map%polar = 1
        end if
        if (st%have_Z_pair) map%polar = 2
        if (st%have_Z_triplet) map%polar = 2
        !if ( st%have_eps_singlet ) map%polar=2
        if (st%have_eps_pair) map%polar = 2
    end block init

    ! Get the forceconstant singlets
    if (st%nx_fc_singlet .gt. 0) then
        fcsinglet: block
            integer :: ish, i, isinglet, iop, a1

            ! make a note that we have pairs, and store some basic information
            map%have_fc_singlet = .true.                 ! we have fc singlets
            map%n_fc_singlet_shell = st%n_fc_singlet_shell  ! number of shells
            map%xuc%nx_fc_singlet = st%nx_fc_singlet       ! number of irreducible
            allocate (map%xuc%x_fc_singlet(map%xuc%nx_fc_singlet))
            map%xuc%x_fc_singlet = 0.0_r8                 ! initialize irreducible to nothing

            ! Populate the shells, can just do a straight copy f2008 style
            call copyshells(st%fc_singlet_shell, map%fc_singlet_shell)

            ! Count singlets?
            map%xuc%n_fc_singlet = 0
            do ish = 1, map%n_fc_singlet_shell
                map%xuc%n_fc_singlet = map%xuc%n_fc_singlet + map%fc_singlet_shell(ish)%n_unfold
            end do
            ! Space for unitcell singlets
            allocate (map%xuc%fc_singlet(map%xuc%n_fc_singlet))
            ! Store unitcell singlets
            do ish = 1, map%n_fc_singlet_shell
            do i = 1, map%fc_singlet_shell(ish)%n_unfold
                isinglet = st%fc_singlet_shell(ish)%unfold_index(i)  ! singlet index
                iop = st%fc_singlet_shell(ish)%unfold_operation(i)   ! Operation from shell
                map%xuc%fc_singlet(isinglet)%irreducible_shell = ish
                map%xuc%fc_singlet(isinglet)%operation_from_shell = iop
            end do
            end do

            ! Build the supercell things. Much easier in this case.
            map%xss%n_fc_singlet = ss%na
            allocate (map%xss%ind_fc_singlet(2, ss%na))
            do i = 1, map%xss%n_fc_singlet
                a1 = ss%info%index_in_unitcell(i)
                map%xss%ind_fc_singlet(1, i) = map%xuc%fc_singlet(a1)%irreducible_shell
                map%xss%ind_fc_singlet(2, i) = map%xuc%fc_singlet(a1)%operation_from_shell
            end do
        end block fcsinglet
    else
        ! no singlets
        map%have_fc_singlet = .false.
    end if

    ! Start with forceconstant pairs
    if (st%nx_fc_pair .gt. 0) then
        fcpair: block
            real(r8), dimension(3) :: r0, v0, l0
            integer, dimension(:, :), allocatable :: dj
            integer, dimension(:), allocatable :: di, dk
            integer :: ish, i, iop, ipair, a1, a2, l

            ! make a note that we have pairs, and store some basic information
            map%have_fc_pair = .true.                 ! we have fc pairs
            map%n_fc_pair_shell = st%n_fc_pair_shell  ! number of shells
            map%xuc%nx_fc_pair = st%nx_fc_pair        ! number of irreducible
            allocate (map%xuc%x_fc_pair(map%xuc%nx_fc_pair))
            map%xuc%x_fc_pair = 0.0_r8                ! initialize irreducible to nothing

            ! Populate the shells, can just do a straight copy f2008 style
            ! Actually it seems I can not, that made ifort very angry and confused.
            call copyshells(st%fc_pair_shell, map%fc_pair_shell)

            ! Count pairs?
            map%xuc%n_fc_pair = 0
            do ish = 1, map%n_fc_pair_shell
                map%xuc%n_fc_pair = map%xuc%n_fc_pair + map%fc_pair_shell(ish)%n_unfold
            end do

            ! Space for unitcell pairs
            allocate (map%xuc%fc_pair(map%xuc%n_fc_pair))
            ! Store unitcell pairs
            do ish = 1, map%n_fc_pair_shell
            do i = 1, map%fc_pair_shell(ish)%n_unfold
                ipair = st%fc_pair_shell(ish)%unfold_index(i)     ! pair index
                iop = st%fc_pair_shell(ish)%unfold_operation(i)   ! Operation from shell
                a1 = st%uc_pair(ipair)%i1                         ! index in unit cell to atom 1
                a2 = st%uc_pair(ipair)%i2                         ! index in unit cell to atom 2
                r0 = st%uc_pair(ipair)%v                          ! Cartesian vector
                v0 = r0 - uc%rcart(:, a2) + uc%rcart(:, a1)             ! Cartesian lattice vector
                l0 = anint(matmul(uc%inv_latticevectors, v0))      ! Fractional lattice vector
                v0 = matmul(uc%latticevectors, l0)                 ! Cartesian lattice vector, clean

                map%xuc%fc_pair(ipair)%irreducible_shell = ish
                map%xuc%fc_pair(ipair)%operation_from_shell = iop
                map%xuc%fc_pair(ipair)%i1 = a1
                map%xuc%fc_pair(ipair)%i2 = a2
                map%xuc%fc_pair(ipair)%r = r0
                map%xuc%fc_pair(ipair)%lv = v0
                map%xuc%fc_pair(ipair)%flv = l0
                if (lo_sqnorm(r0) .lt. lo_sqtol) then
                    map%xuc%fc_pair(ipair)%selfterm = .true.
                else
                    map%xuc%fc_pair(ipair)%selfterm = .false.
                end if
            end do
            end do

            ! Sort things out for the supercell. I have made this a little confusing, but that
            ! is to make sure it's completely general. Think it works alright. Technically I have
            ! the proper list of supercell pairs already, but I want to generate it procedurally
            ! since then I can generalize to all manner of tuplets and interactions with minimal
            ! change in the code.
            call mem%allocate(di, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            do ipair = 1, map%xuc%n_fc_pair
                a1 = map%xuc%fc_pair(ipair)%i1
                di(a1) = di(a1) + 1
            end do
            call mem%allocate(dj, [maxval(di), uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dk, ss%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            dj = 0
            dk = 0
            do ipair = 1, map%xuc%n_fc_pair
                a1 = map%xuc%fc_pair(ipair)%i1
                di(a1) = di(a1) + 1
                dj(di(a1), a1) = ipair
            end do
            ! So, now di holds the number of pairs per unitcell atom, and dj holds
            ! the list of pairs per unitcell atom, with that I can easily construct
            ! the list of supercell pairs
            map%xss%n_fc_pair = 0
            do a1 = 1, map%n_atom_ss
                a2 = ss%info%index_in_unitcell(a1)
                dk(a1) = map%xss%n_fc_pair
                map%xss%n_fc_pair = map%xss%n_fc_pair + di(a2)
            end do
            allocate (map%xss%ind_fc_pair(4, map%xss%n_fc_pair))
            map%xss%ind_fc_pair = 0
            ! Fetch irreducible shell, operation and first atom.
            do a1 = 1, map%n_atom_ss
                if (mod(a1, mw%n) .ne. mw%r) cycle
                a2 = ss%info%index_in_unitcell(a1)
                do i = 1, di(a2)
                    ipair = dj(i, a2)
                    l = dk(a1) + i
                    map%xss%ind_fc_pair(1, l) = map%xuc%fc_pair(ipair)%irreducible_shell
                    map%xss%ind_fc_pair(2, l) = map%xuc%fc_pair(ipair)%operation_from_shell
                    map%xss%ind_fc_pair(3, l) = a1
                    map%xss%ind_fc_pair(4, l) = -ipair
                end do
            end do
            call mw%allreduce('sum', map%xss%ind_fc_pair)
            ! Get the final index for each supercell pair
            do l = 1, map%xss%n_fc_pair
                if (mod(l, mw%n) .ne. mw%r) cycle
                a1 = map%xss%ind_fc_pair(3, l)
                ipair = -map%xss%ind_fc_pair(4, l)
                v0 = ss%rcart(:, a1) + map%xuc%fc_pair(ipair)%r
                v0 = matmul(ss%inv_latticevectors, v0)
                v0 = lo_clean_fractional_coordinates(v0)
                a2 = 0
                do i = 1, ss%na
                    r0 = lo_clean_fractional_coordinates(ss%r(:, i) - v0 + 0.5_r8) - 0.5_r8
                    if (lo_sqnorm(r0) .lt. lo_sqtol) then
                        a2 = i
                        exit
                    end if
                end do
                if (a2 .gt. 0) then
                    map%xss%ind_fc_pair(4, l) = a2
                else
                    call lo_stop_gracefully(['Could not locate pair'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                end if
            end do
            call mw%allreduce('max', map%xss%ind_fc_pair)
            ! And that should be it for the pairs!
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... built pair forcemap (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
            ! Cleanup
            call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dk, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block fcpair
    end if

    if (st%nx_fc_triplet .gt. 0) then
        fctriplet: block
            real(r8), dimension(3) :: r2, r3, v2, v3, l2, l3
            integer, dimension(:, :), allocatable :: dj
            integer, dimension(:), allocatable :: di, dk
            integer :: ish, i, j, iop, itriplet, a1, a2, a3, l

            ! make a note that we have triplets, and store some basic information
            map%have_fc_triplet = .true.                 ! we have fc triplets
            map%n_fc_triplet_shell = st%n_fc_triplet_shell  ! number of shells
            map%xuc%nx_fc_triplet = st%nx_fc_triplet       ! number of irreducible
            allocate (map%xuc%x_fc_triplet(map%xuc%nx_fc_triplet))
            map%xuc%x_fc_triplet = 0.0_r8                 ! initialize irreducible to nothing

            ! Populate the shells, can just do a straight copy f2008 style
            call copyshells(st%fc_triplet_shell, map%fc_triplet_shell)

            ! Count triplets?
            map%xuc%n_fc_triplet = 0
            do ish = 1, map%n_fc_triplet_shell
                map%xuc%n_fc_triplet = map%xuc%n_fc_triplet + map%fc_triplet_shell(ish)%n_unfold
            end do
            ! Space for unitcell triplets
            allocate (map%xuc%fc_triplet(map%xuc%n_fc_triplet))
            ! Store unitcell triplets
            do ish = 1, map%n_fc_triplet_shell
            do i = 1, map%fc_triplet_shell(ish)%n_unfold
                itriplet = st%fc_triplet_shell(ish)%unfold_index(i)  ! triplet index
                iop = st%fc_triplet_shell(ish)%unfold_operation(i)   ! Operation from shell
                a1 = st%uc_triplet(itriplet)%i1                      ! index in unit cell to atom 1
                a2 = st%uc_triplet(itriplet)%i2                      ! index in unit cell to atom 2
                a3 = st%uc_triplet(itriplet)%i3                      ! index in unit cell to atom 2
                !r1=st%uc_triplet(itriplet)%v1                      ! Cartesian vector
                r2 = st%uc_triplet(itriplet)%v2                      ! Cartesian vector
                r3 = st%uc_triplet(itriplet)%v3                      ! Cartesian vector
                !v1=r1-uc%rcart(:,a1)+uc%rcart(:,a1)                ! Cartesian lattice vector
                v2 = r2 - uc%rcart(:, a2) + uc%rcart(:, a1)                ! Cartesian lattice vector
                v3 = r3 - uc%rcart(:, a3) + uc%rcart(:, a1)                ! Cartesian lattice vector
                !l0=anint(matmul(uc%inv_latticevectors,v1))         ! Fractional lattice vector
                l2 = anint(matmul(uc%inv_latticevectors, v2))         ! Fractional lattice vector
                l3 = anint(matmul(uc%inv_latticevectors, v3))         ! Fractional lattice vector
                !v1=matmul(uc%latticevectors,l1)                    ! Cartesian lattice vector, clean
                v2 = matmul(uc%latticevectors, l2)                    ! Cartesian lattice vector, clean
                v3 = matmul(uc%latticevectors, l3)                    ! Cartesian lattice vector, clean

                map%xuc%fc_triplet(itriplet)%irreducible_shell = ish
                map%xuc%fc_triplet(itriplet)%operation_from_shell = iop
                map%xuc%fc_triplet(itriplet)%i1 = a1
                map%xuc%fc_triplet(itriplet)%i2 = a2
                map%xuc%fc_triplet(itriplet)%i3 = a3
                map%xuc%fc_triplet(itriplet)%v2 = r2
                map%xuc%fc_triplet(itriplet)%v3 = r3
                map%xuc%fc_triplet(itriplet)%lv2 = v2
                map%xuc%fc_triplet(itriplet)%lv3 = v3
                map%xuc%fc_triplet(itriplet)%flv2 = l2
                map%xuc%fc_triplet(itriplet)%flv3 = l3
                if (lo_sqnorm(r2) + lo_sqnorm(r3) .lt. lo_sqtol) then
                    map%xuc%fc_triplet(itriplet)%selfterm = .true.
                else
                    map%xuc%fc_triplet(itriplet)%selfterm = .false.
                end if
            end do
            end do

            ! Sort things out for the supercell. I have made this a little confusing, but that
            ! is to make sure it's completely general. Think it works alright. Technically I have
            ! the proper list of supercell triplets already, but I want to generate it procedurally
            ! since then I can generalize to all manner of tuplets and interactions with minimal
            ! change in the code.
            call mem%allocate(di, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            do itriplet = 1, map%xuc%n_fc_triplet
                a1 = map%xuc%fc_triplet(itriplet)%i1
                di(a1) = di(a1) + 1
            end do
            call mem%allocate(dj, [maxval(di), uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dk, ss%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            dj = 0
            dk = 0
            do itriplet = 1, map%xuc%n_fc_triplet
                a1 = map%xuc%fc_triplet(itriplet)%i1
                di(a1) = di(a1) + 1
                dj(di(a1), a1) = itriplet
            end do
            ! So, now di holds the number of triplets per unitcell atom, and dj holds
            ! the list of triplets per unitcell atom, with that I can easily construct
            ! the list of supercell triplets
            map%xss%n_fc_triplet = 0
            do a1 = 1, map%n_atom_ss
                a2 = ss%info%index_in_unitcell(a1)
                dk(a1) = map%xss%n_fc_triplet
                map%xss%n_fc_triplet = map%xss%n_fc_triplet + di(a2)
            end do
            allocate (map%xss%ind_fc_triplet(5, map%xss%n_fc_triplet))
            map%xss%ind_fc_triplet = 0
            ! Fetch irreducible shell, operation and first atom.
            do a1 = 1, map%n_atom_ss
                if (mod(a1, mw%n) .ne. mw%r) cycle
                a2 = ss%info%index_in_unitcell(a1)
                do i = 1, di(a2)
                    itriplet = dj(i, a2)
                    l = dk(a1) + i
                    map%xss%ind_fc_triplet(1, l) = map%xuc%fc_triplet(itriplet)%irreducible_shell
                    map%xss%ind_fc_triplet(2, l) = map%xuc%fc_triplet(itriplet)%operation_from_shell
                    map%xss%ind_fc_triplet(3, l) = a1
                    map%xss%ind_fc_triplet(4, l) = -itriplet
                    map%xss%ind_fc_triplet(5, l) = -itriplet
                end do
            end do
            call mw%allreduce('sum', map%xss%ind_fc_triplet)
            ! Get the final index for each supercell triplet
            do l = 1, map%xss%n_fc_triplet
                if (mod(l, mw%n) .ne. mw%r) cycle
                a1 = map%xss%ind_fc_triplet(3, l)
                itriplet = -map%xss%ind_fc_triplet(4, l)
                v2 = ss%rcart(:, a1) + map%xuc%fc_triplet(itriplet)%v2
                v3 = ss%rcart(:, a1) + map%xuc%fc_triplet(itriplet)%v3
                v2 = matmul(ss%inv_latticevectors, v2)
                v3 = matmul(ss%inv_latticevectors, v3)
                v2 = lo_clean_fractional_coordinates(v2)
                v3 = lo_clean_fractional_coordinates(v3)
                a2 = 0
                a3 = 0
                j = 0
                do i = 1, ss%na
                    r2 = lo_clean_fractional_coordinates(ss%r(:, i) - v2 + 0.5_r8) - 0.5_r8
                    r3 = lo_clean_fractional_coordinates(ss%r(:, i) - v3 + 0.5_r8) - 0.5_r8
                    if (lo_sqnorm(r2) .lt. lo_sqtol) then
                        a2 = i
                        j = j + 1
                    end if
                    if (lo_sqnorm(r3) .lt. lo_sqtol) then
                        a3 = i
                        j = j + 1
                    end if
                    if (j .eq. 2) exit
                end do
                if (min(a2, a3) .gt. 0) then
                    map%xss%ind_fc_triplet(4, l) = a2
                    map%xss%ind_fc_triplet(5, l) = a3
                else
                    call lo_stop_gracefully(['Could not locate triplet'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                end if
            end do
            call mw%allreduce('max', map%xss%ind_fc_triplet)
            ! And sort out the ASR tuplets
            call thirdorder_asr_tuplets(map, mem)
            ! And that should be it for the triplets!
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... built triplet forcemap (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
            ! Cleanup
            call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dk, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block fctriplet
    end if

    if (st%nx_fc_quartet .gt. 0) then
        fcquartet: block
            real(r8), dimension(3) :: r2, r3, r4, v2, v3, v4, l2, l3, l4
            integer, dimension(:, :), allocatable :: dj
            integer, dimension(:), allocatable :: di, dk
            integer :: ish, i, j, iop, iquartet, a1, a2, a3, a4, l

            ! make a note that we have quartets, and store some basic information
            map%have_fc_quartet = .true.                 ! we have fc quartets
            map%n_fc_quartet_shell = st%n_fc_quartet_shell  ! number of shells
            map%xuc%nx_fc_quartet = st%nx_fc_quartet       ! number of irreducible
            allocate (map%xuc%x_fc_quartet(map%xuc%nx_fc_quartet))
            map%xuc%x_fc_quartet = 0.0_r8                 ! initialize irreducible to nothing

            ! Populate the shells, can just do a straight copy f2008 style
            call copyshells(st%fc_quartet_shell, map%fc_quartet_shell)

            ! Count quartets?
            map%xuc%n_fc_quartet = 0
            do ish = 1, map%n_fc_quartet_shell
                map%xuc%n_fc_quartet = map%xuc%n_fc_quartet + map%fc_quartet_shell(ish)%n_unfold
            end do
            ! Space for unitcell quartets
            allocate (map%xuc%fc_quartet(map%xuc%n_fc_quartet))
            ! Store unitcell quartets
            do ish = 1, map%n_fc_quartet_shell
            do i = 1, map%fc_quartet_shell(ish)%n_unfold
                iquartet = st%fc_quartet_shell(ish)%unfold_index(i)  ! quartet index
                iop = st%fc_quartet_shell(ish)%unfold_operation(i)   ! Operation from shell
                a1 = st%uc_quartet(iquartet)%i1                      ! index in unit cell to atom 1
                a2 = st%uc_quartet(iquartet)%i2                      ! index in unit cell to atom 2
                a3 = st%uc_quartet(iquartet)%i3                      ! index in unit cell to atom 3
                a4 = st%uc_quartet(iquartet)%i4                      ! index in unit cell to atom 4
                !r1=st%uc_quartet(iquartet)%v1                      ! Cartesian vector
                r2 = st%uc_quartet(iquartet)%v2                      ! Cartesian vector
                r3 = st%uc_quartet(iquartet)%v3                      ! Cartesian vector
                r4 = st%uc_quartet(iquartet)%v4                      ! Cartesian vector
                !v1=r1-uc%rcart(:,a1)+uc%rcart(:,a1)                ! Cartesian lattice vector
                v2 = r2 - uc%rcart(:, a2) + uc%rcart(:, a1)                ! Cartesian lattice vector
                v3 = r3 - uc%rcart(:, a3) + uc%rcart(:, a1)                ! Cartesian lattice vector
                v4 = r4 - uc%rcart(:, a4) + uc%rcart(:, a1)                ! Cartesian lattice vector
                !l0=anint(matmul(uc%inv_latticevectors,v1))         ! Fractional lattice vector
                l2 = anint(matmul(uc%inv_latticevectors, v2))         ! Fractional lattice vector
                l3 = anint(matmul(uc%inv_latticevectors, v3))         ! Fractional lattice vector
                l4 = anint(matmul(uc%inv_latticevectors, v4))         ! Fractional lattice vector
                !v1=matmul(uc%latticevectors,l1)                    ! Cartesian lattice vector, clean
                v2 = matmul(uc%latticevectors, l2)                    ! Cartesian lattice vector, clean
                v3 = matmul(uc%latticevectors, l3)                    ! Cartesian lattice vector, clean
                v4 = matmul(uc%latticevectors, l4)                    ! Cartesian lattice vector, clean

                map%xuc%fc_quartet(iquartet)%irreducible_shell = ish
                map%xuc%fc_quartet(iquartet)%operation_from_shell = iop
                map%xuc%fc_quartet(iquartet)%i1 = a1
                map%xuc%fc_quartet(iquartet)%i2 = a2
                map%xuc%fc_quartet(iquartet)%i3 = a3
                map%xuc%fc_quartet(iquartet)%i4 = a4
                map%xuc%fc_quartet(iquartet)%v2 = r2
                map%xuc%fc_quartet(iquartet)%v3 = r3
                map%xuc%fc_quartet(iquartet)%v4 = r4
                map%xuc%fc_quartet(iquartet)%lv2 = v2
                map%xuc%fc_quartet(iquartet)%lv3 = v3
                map%xuc%fc_quartet(iquartet)%lv4 = v4
                map%xuc%fc_quartet(iquartet)%flv2 = l2
                map%xuc%fc_quartet(iquartet)%flv3 = l3
                map%xuc%fc_quartet(iquartet)%flv4 = l4
                if (lo_sqnorm(r2) + lo_sqnorm(r3) + lo_sqnorm(r4) .lt. lo_sqtol) then
                    map%xuc%fc_quartet(iquartet)%selfterm = .true.
                else
                    map%xuc%fc_quartet(iquartet)%selfterm = .false.
                end if
            end do
            end do

            ! Sort things out for the supercell. I have made this a little confusing, but that
            ! is to make sure it's completely general. Think it works alright. Technically I have
            ! the proper list of supercell quartets already, but I want to generate it procedurally
            ! since then I can generalize to all manner of tuplets and interactions with minimal
            ! change in the code.
            call mem%allocate(di, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            do iquartet = 1, map%xuc%n_fc_quartet
                a1 = map%xuc%fc_quartet(iquartet)%i1
                di(a1) = di(a1) + 1
            end do
            call mem%allocate(dj, [maxval(di), uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dk, ss%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            dj = 0
            dk = 0
            do iquartet = 1, map%xuc%n_fc_quartet
                a1 = map%xuc%fc_quartet(iquartet)%i1
                di(a1) = di(a1) + 1
                dj(di(a1), a1) = iquartet
            end do
            ! So, now di holds the number of quartets per unitcell atom, and dj holds
            ! the list of quartets per unitcell atom, with that I can easily construct
            ! the list of supercell quartets
            map%xss%n_fc_quartet = 0
            do a1 = 1, map%n_atom_ss
                a2 = ss%info%index_in_unitcell(a1)
                dk(a1) = map%xss%n_fc_quartet
                map%xss%n_fc_quartet = map%xss%n_fc_quartet + di(a2)
            end do
            allocate (map%xss%ind_fc_quartet(6, map%xss%n_fc_quartet))
            map%xss%ind_fc_quartet = 0
            ! Fetch irreducible shell, operation and first atom.
            do a1 = 1, map%n_atom_ss
                if (mod(a1, mw%n) .ne. mw%r) cycle
                a2 = ss%info%index_in_unitcell(a1)
                do i = 1, di(a2)
                    iquartet = dj(i, a2)
                    l = dk(a1) + i
                    map%xss%ind_fc_quartet(1, l) = map%xuc%fc_quartet(iquartet)%irreducible_shell
                    map%xss%ind_fc_quartet(2, l) = map%xuc%fc_quartet(iquartet)%operation_from_shell
                    map%xss%ind_fc_quartet(3, l) = a1
                    map%xss%ind_fc_quartet(4, l) = -iquartet
                    map%xss%ind_fc_quartet(5, l) = -iquartet
                    map%xss%ind_fc_quartet(6, l) = -iquartet
                end do
            end do
            call mw%allreduce('sum', map%xss%ind_fc_quartet)
            ! Get the final index for each supercell quartet
            do l = 1, map%xss%n_fc_quartet
                if (mod(l, mw%n) .ne. mw%r) cycle
                a1 = map%xss%ind_fc_quartet(3, l)
                iquartet = -map%xss%ind_fc_quartet(4, l)
                v2 = ss%rcart(:, a1) + map%xuc%fc_quartet(iquartet)%v2
                v3 = ss%rcart(:, a1) + map%xuc%fc_quartet(iquartet)%v3
                v4 = ss%rcart(:, a1) + map%xuc%fc_quartet(iquartet)%v4
                v2 = matmul(ss%inv_latticevectors, v2)
                v3 = matmul(ss%inv_latticevectors, v3)
                v4 = matmul(ss%inv_latticevectors, v4)
                v2 = lo_clean_fractional_coordinates(v2)
                v3 = lo_clean_fractional_coordinates(v3)
                v4 = lo_clean_fractional_coordinates(v4)
                a2 = 0
                a3 = 0
                a4 = 0
                j = 0
                do i = 1, ss%na
                    r2 = lo_clean_fractional_coordinates(ss%r(:, i) - v2 + 0.5_r8) - 0.5_r8
                    r3 = lo_clean_fractional_coordinates(ss%r(:, i) - v3 + 0.5_r8) - 0.5_r8
                    r4 = lo_clean_fractional_coordinates(ss%r(:, i) - v4 + 0.5_r8) - 0.5_r8
                    if (lo_sqnorm(r2) .lt. lo_sqtol) then
                        a2 = i
                        j = j + 1
                    end if
                    if (lo_sqnorm(r3) .lt. lo_sqtol) then
                        a3 = i
                        j = j + 1
                    end if
                    if (lo_sqnorm(r4) .lt. lo_sqtol) then
                        a4 = i
                        j = j + 1
                    end if
                    if (j .eq. 3) exit
                end do
                if (min(a2, a3, a4) .gt. 0) then
                    map%xss%ind_fc_quartet(4, l) = a2
                    map%xss%ind_fc_quartet(5, l) = a3
                    map%xss%ind_fc_quartet(6, l) = a4
                else
                    call lo_stop_gracefully(['Could not locate quartet'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                end if
            end do
            call mw%allreduce('max', map%xss%ind_fc_quartet)
            ! And sort out the ASR tuplets
            call fourthorder_asr_tuplets(map, mem)
            ! And that should be it for the quartets!
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... built quartet forcemap (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
            ! Cleanup
            call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dk, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block fcquartet
    end if

    ! Dielectric constant
    if (st%nx_eps_global .gt. 0) then
        epsglobal: block
            ! space for solution
            map%xuc%nx_eps_global = st%nx_eps_global
            allocate (map%xuc%x_eps_global(map%xuc%nx_eps_global))
            allocate (map%xuc%x_eps_global_deriv(map%xuc%nx_eps_global))
            map%xuc%x_eps_global = 0.0_r8
            map%xuc%x_eps_global_deriv = 0.0_r8
            ! copy of shell
            map%eps_global_shell%nx = st%eps_global_shell%nx
            if (map%eps_global_shell%nx .gt. 0) then
                allocate (map%eps_global_shell%coeff(size(st%eps_global_shell%coeff, 1), map%eps_global_shell%nx))
                map%eps_global_shell%coeff = st%eps_global_shell%coeff
            end if
            map%eps_global_shell%n_unfold = -1

            ! how to correct
            map%polarcorrectiontype = polarcorrectiontype
        end block epsglobal
    end if

    ! Get the borncharge singlets
    if (st%nx_eps_singlet .gt. 0) then
        epssinglet: block
            integer :: ish, i, isinglet, iop, a1

            ! make a note that we have singlets, and store some basic information
            map%have_eps_singlet = .true.
            map%n_eps_singlet_shell = st%n_eps_singlet_shell  ! number of shells
            map%xuc%nx_eps_singlet = st%nx_eps_singlet       ! number of irreducible
            allocate (map%xuc%x_eps_singlet(map%xuc%nx_eps_singlet))
            map%xuc%x_eps_singlet = 0.0_r8                 ! initialiepse irreducible to nothing

            ! Populate the shells, can just do a straight copy f2008 style
            call copyshells(st%eps_singlet_shell, map%eps_singlet_shell)

            ! Count singlets?
            map%xuc%n_eps_singlet = 0
            do ish = 1, map%n_eps_singlet_shell
                map%xuc%n_eps_singlet = map%xuc%n_eps_singlet + map%eps_singlet_shell(ish)%n_unfold
            end do
            ! Space for unitcell singlets
            allocate (map%xuc%eps_singlet(map%xuc%n_eps_singlet))
            ! Store unitcell singlets
            do ish = 1, map%n_eps_singlet_shell
            do i = 1, map%eps_singlet_shell(ish)%n_unfold
                isinglet = st%eps_singlet_shell(ish)%unfold_index(i)  ! singlet index
                iop = st%eps_singlet_shell(ish)%unfold_operation(i)   ! Operation from shell
                map%xuc%eps_singlet(isinglet)%irreducible_shell = ish
                map%xuc%eps_singlet(isinglet)%operation_from_shell = iop
            end do
            end do

            ! Build the supercell things. Much easier in this case.
            map%xss%n_eps_singlet = ss%na
            allocate (map%xss%ind_eps_singlet(2, ss%na))
            map%xss%ind_eps_singlet = 0
            do i = 1, map%xss%n_eps_singlet
                a1 = ss%info%index_in_unitcell(i)
                map%xss%ind_eps_singlet(1, i) = map%xuc%eps_singlet(a1)%irreducible_shell
                map%xss%ind_eps_singlet(2, i) = map%xuc%eps_singlet(a1)%operation_from_shell
            end do
        end block epssinglet
    end if

    ! Dielectric pairs
    if (st%nx_eps_pair .gt. 0) then
        epspair: block
            real(r8), dimension(3) :: r0, v0, l0
            integer, dimension(:, :), allocatable :: dj
            integer, dimension(:), allocatable :: di, dk, ind_to_quartet_op
            integer :: ish, iop, jop, ipair, jpair, a1, a2, i, j, k, l, ctr

            ! make a note that we have pairs, and store some basic information
            map%have_eps_pair = .true.              ! we have eps pairs
            map%n_eps_pair_shell = st%n_eps_pair_shell ! number of shells
            map%xuc%nx_eps_pair = st%nx_eps_pair      ! number of irreducible
            allocate (map%xuc%x_eps_pair(map%xuc%nx_eps_pair))
            map%xuc%x_eps_pair = 0.0_r8              ! initialiepse irreducible to nothing

            ! I need to know how to translate from pair operations to triplet operations. On entry, all the
            ! indices to operations point to pair operations, I need to redirect those indices to triplet
            ! operations. This is confusing, I know.
            call mem%allocate(ind_to_quartet_op, st%n_pair_operation, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            ind_to_quartet_op = 0
            do i = 1, st%n_pair_operation
                k = 0
                ctr = 0
                do j = 1, st%n_quartet_operation
                    ! have to be the same rotation
                    if (st%quartetop(j)%opind .ne. st%pairop(i)%opind) cycle
                    ! can not permute the electric field indices
                    if (st%quartetop(j)%perm(1) .ne. 1) cycle
                    if (st%quartetop(j)%perm(2) .ne. 2) cycle
                    ! next one is annoying.
                    if (sum(abs(st%pairop(i)%perm - [1, 2])) .eq. 0) then
                        if (sum(abs(st%quartetop(j)%perm(3:4) - [3, 4])) .ne. 0) cycle
                    else
                        if (sum(abs(st%quartetop(j)%perm(3:4) - [4, 3])) .ne. 0) cycle
                    end if
                    ctr = ctr + 1
                    k = j
                end do
                if (ctr .gt. 1) call lo_stop_gracefully(['More than one operation matches, impossible'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                if (k .eq. 0) then
                    call lo_stop_gracefully(['Could not find operation, impossible'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                else
                    ind_to_quartet_op(i) = k
                end if
            end do

            ! Populate the shells, can just do a straight copy f2008 style
            call copyshells(st%eps_pair_shell, map%eps_pair_shell)

            ! Count pairs?
            map%xuc%n_eps_pair = 0
            do ish = 1, map%n_eps_pair_shell
                map%xuc%n_eps_pair = map%xuc%n_eps_pair + map%eps_pair_shell(ish)%n_unfold
            end do
            ! Space for unitcell pairs
            allocate (map%xuc%eps_pair(map%xuc%n_eps_pair))
            ! Store unitcell pairs
            jpair = 0
            do ish = 1, map%n_eps_pair_shell
            do i = 1, map%eps_pair_shell(ish)%n_unfold
                jpair = jpair + 1
                ipair = st%eps_pair_shell(ish)%unfold_index(i)    ! pair index in all pairs
                iop = st%eps_pair_shell(ish)%unfold_operation(i)  ! Operation from shell, points to pair operation
                jop = ind_to_quartet_op(iop)                      ! Operation from shell, points to quartet operation
                a1 = st%uc_pair(ipair)%i1                         ! index in unit cell to atom 1
                a2 = st%uc_pair(ipair)%i2                         ! index in unit cell to atom 2
                r0 = st%uc_pair(ipair)%v                          ! Cartesian vector
                v0 = r0 - uc%rcart(:, a2) + uc%rcart(:, a1)             ! Cartesian lattice vector
                l0 = matmul(uc%inv_latticevectors, v0)             ! Fractional lattice vector
                if (sum(abs(l0 - anint(l0))) .gt. lo_sqtol) then
                    call lo_stop_gracefully(['Clearly I do not undestand vectors'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                else
                    l0 = anint(l0)
                end if
                v0 = matmul(uc%latticevectors, l0)                 ! Cartesian lattice vector, clean

                map%xuc%eps_pair(jpair)%irreducible_shell = ish
                map%xuc%eps_pair(jpair)%operation_from_shell = jop ! Note, points to triplet operation
                map%xuc%eps_pair(jpair)%i1 = a1
                map%xuc%eps_pair(jpair)%i2 = a2
                map%xuc%eps_pair(jpair)%r = r0
                map%xuc%eps_pair(jpair)%lv = v0
                map%xuc%eps_pair(jpair)%flv = l0
                if (lo_sqnorm(r0) .lt. lo_sqtol) then
                    map%xuc%eps_pair(jpair)%selfterm = .true.
                else
                    map%xuc%eps_pair(jpair)%selfterm = .false.
                end if
                ! and finally, make sure that the information in the shell makes sense
                ! and point to the right kind of pairs.
                map%eps_pair_shell(ish)%unfold_index(i) = jpair
                ! Don't know if useful to keep this, for now I set it to nothing so I don't use it the wrong way.
                map%eps_pair_shell(ish)%unfold_operation(i) = -1
            end do
            end do

            ! Sort things out for the supercell. I have made this a little confusing, but that
            ! is to make sure it's completely general. Think it works alright. Technically I have
            ! the proper list of supercell pairs already, but I want to generate it procedurally
            ! since then I can generaliepse to all manner of tuplets and interactions with minimal
            ! change in the code.
            call mem%allocate(di, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            do ipair = 1, map%xuc%n_eps_pair
                a1 = map%xuc%eps_pair(ipair)%i1
                di(a1) = di(a1) + 1
            end do
            call mem%allocate(dj, [maxval(di), uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dk, ss%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            dj = 0
            dk = 0
            do ipair = 1, map%xuc%n_eps_pair
                a1 = map%xuc%eps_pair(ipair)%i1
                di(a1) = di(a1) + 1
                dj(di(a1), a1) = ipair
            end do
            ! So, now di holds the number of pairs per unitcell atom, and dj holds
            ! the list of pairs per unitcell atom, with that I can easily construct
            ! the list of supercell pairs
            map%xss%n_eps_pair = 0
            do a1 = 1, map%n_atom_ss
                a2 = ss%info%index_in_unitcell(a1)
                dk(a1) = map%xss%n_eps_pair
                map%xss%n_eps_pair = map%xss%n_eps_pair + di(a2)
            end do
            allocate (map%xss%ind_eps_pair(4, map%xss%n_eps_pair))
            map%xss%ind_eps_pair = 0
            ! Fetch irreducible shell, operation and first atom.
            do a1 = 1, map%n_atom_ss
                if (mod(a1, mw%n) .ne. mw%r) cycle
                a2 = ss%info%index_in_unitcell(a1)
                do i = 1, di(a2)
                    ipair = dj(i, a2)
                    l = dk(a1) + i
                    map%xss%ind_eps_pair(1, l) = map%xuc%eps_pair(ipair)%irreducible_shell
                    map%xss%ind_eps_pair(2, l) = map%xuc%eps_pair(ipair)%operation_from_shell
                    map%xss%ind_eps_pair(3, l) = a1
                    map%xss%ind_eps_pair(4, l) = -ipair
                end do
            end do
            call mw%allreduce('sum', map%xss%ind_eps_pair)
            ! Get the final index for each supercell pair
            do l = 1, map%xss%n_eps_pair
                if (mod(l, mw%n) .ne. mw%r) cycle
                a1 = map%xss%ind_eps_pair(3, l)
                ipair = -map%xss%ind_eps_pair(4, l)
                v0 = ss%rcart(:, a1) + map%xuc%eps_pair(ipair)%r
                v0 = matmul(ss%inv_latticevectors, v0)
                v0 = lo_clean_fractional_coordinates(v0)
                a2 = 0
                do i = 1, ss%na
                    r0 = lo_clean_fractional_coordinates(ss%r(:, i) - v0 + 0.5_r8) - 0.5_r8
                    if (lo_sqnorm(r0) .lt. lo_sqtol) then
                        a2 = i
                        exit
                    end if
                end do
                if (a2 .gt. 0) then
                    map%xss%ind_eps_pair(4, l) = a2
                else
                    call lo_stop_gracefully(['Could not locate pair'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                end if
            end do
            call mw%allreduce('max', map%xss%ind_eps_pair)
            ! And that should be it for the pairs!
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... built pair forcemap (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
            ! Cleanup
            call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dk, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ind_to_quartet_op, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block epspair
    end if

    ! Get the borncharge singlets
    if (st%nx_Z_singlet .gt. 0) then
        Zsinglet: block
            integer :: ish, i, isinglet, iop, a1

            ! make a note that we have pairs, and store some basic information
            map%n_Z_singlet_shell = st%n_Z_singlet_shell  ! number of shells
            map%xuc%nx_Z_singlet = st%nx_Z_singlet       ! number of irreducible
            allocate (map%xuc%x_Z_singlet(map%xuc%nx_Z_singlet))
            map%xuc%x_Z_singlet = 0.0_r8                 ! initialize irreducible to nothing

            ! Populate the shells, can just do a straight copy f2008 style
            call copyshells(st%Z_singlet_shell, map%Z_singlet_shell)

            ! Count singlets?
            map%xuc%n_Z_singlet = 0
            do ish = 1, map%n_Z_singlet_shell
                map%xuc%n_Z_singlet = map%xuc%n_Z_singlet + map%Z_singlet_shell(ish)%n_unfold
            end do
            ! Space for unitcell singlets
            allocate (map%xuc%Z_singlet(map%xuc%n_Z_singlet))
            ! Store unitcell singlets
            do ish = 1, map%n_Z_singlet_shell
            do i = 1, map%Z_singlet_shell(ish)%n_unfold
                isinglet = st%Z_singlet_shell(ish)%unfold_index(i)  ! singlet index
                iop = st%Z_singlet_shell(ish)%unfold_operation(i)   ! Operation from shell
                map%xuc%Z_singlet(isinglet)%irreducible_shell = ish
                map%xuc%Z_singlet(isinglet)%operation_from_shell = iop
            end do
            end do

            ! Build the supercell things. Much easier in this case.
            map%xss%n_Z_singlet = ss%na
            allocate (map%xss%ind_Z_singlet(2, ss%na))
            do i = 1, map%xss%n_Z_singlet
                a1 = ss%info%index_in_unitcell(i)
                map%xss%ind_Z_singlet(1, i) = map%xuc%Z_singlet(a1)%irreducible_shell
                map%xss%ind_Z_singlet(2, i) = map%xuc%Z_singlet(a1)%operation_from_shell
            end do
        end block Zsinglet
    end if

    ! Borncharge pairs
    if (st%nx_Z_pair .gt. 0) then
        Zpair: block
            real(r8), dimension(3) :: r0, v0, l0
            integer, dimension(:, :), allocatable :: dj
            integer, dimension(:), allocatable :: di, dk, ind_to_triplet_op
            integer :: ish, iop, jop, ipair, jpair, a1, a2, i, j, k, l, ctr

            ! make a note that we have pairs, and store some basic information
            map%have_Z_pair = .true.             ! we have Z pairs
            map%n_Z_pair_shell = st%n_Z_pair_shell  ! number of shells
            map%xuc%nx_Z_pair = st%nx_Z_pair       ! number of irreducible
            allocate (map%xuc%x_Z_pair(map%xuc%nx_Z_pair))
            map%xuc%x_Z_pair = 0.0_r8              ! initialize irreducible to nothing

            ! I need to know how to translate from pair operations to triplet operations. On entry, all the
            ! indices to operations point to pair operations, I need to redirect those indices to triplet
            ! operations. This is confusing, I know.
            call mem%allocate(ind_to_triplet_op, st%n_pair_operation, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            ind_to_triplet_op = 0
            do i = 1, st%n_pair_operation
                k = 0
                ctr = 0
                do j = 1, st%n_triplet_operation
                    ! have to be the same rotation
                    if (st%tripletop(j)%opind .ne. st%pairop(i)%opind) cycle
                    ! can not permute the electric field indices
                    if (st%tripletop(j)%perm(1) .ne. 1) cycle
                    ! next one is annoying.
                    if (sum(abs(st%pairop(i)%perm - [1, 2])) .eq. 0) then
                        if (sum(abs(st%tripletop(j)%perm(2:3) - [2, 3])) .ne. 0) cycle
                    else
                        if (sum(abs(st%tripletop(j)%perm(2:3) - [2, 3])) .ne. 0) cycle
                    end if
                    ctr = ctr + 1
                    k = j
                end do
                if (ctr .gt. 1) call lo_stop_gracefully(['More than one operation matches, impossible'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                if (k .eq. 0) then
                    call lo_stop_gracefully(['Could not find operation, impossible'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                else
                    ind_to_triplet_op(i) = k
                end if
            end do

            ! Populate the shells, can just do a straight copy f2008 style
            call copyshells(st%Z_pair_shell, map%Z_pair_shell)

            ! Count pairs?
            map%xuc%n_Z_pair = 0
            do ish = 1, map%n_Z_pair_shell
                map%xuc%n_Z_pair = map%xuc%n_Z_pair + map%Z_pair_shell(ish)%n_unfold
            end do
            ! Space for unitcell pairs
            allocate (map%xuc%Z_pair(map%xuc%n_Z_pair))
            ! Store unitcell pairs
            jpair = 0
            do ish = 1, map%n_Z_pair_shell
            do i = 1, map%Z_pair_shell(ish)%n_unfold
                jpair = jpair + 1
                ipair = st%Z_pair_shell(ish)%unfold_index(i)      ! pair index in all pairs
                iop = st%Z_pair_shell(ish)%unfold_operation(i)    ! Operation from shell, points to pair operation
                jop = ind_to_triplet_op(iop)                      ! Operation from shell, points to triplet operation
                a1 = st%uc_pair(ipair)%i1                         ! index in unit cell to atom 1
                a2 = st%uc_pair(ipair)%i2                         ! index in unit cell to atom 2
                r0 = st%uc_pair(ipair)%v                          ! Cartesian vector
                v0 = r0 - uc%rcart(:, a2) + uc%rcart(:, a1)             ! Cartesian lattice vector
                l0 = anint(matmul(uc%inv_latticevectors, v0))      ! Fractional lattice vector
                v0 = matmul(uc%latticevectors, l0)                 ! Cartesian lattice vector, clean

                map%xuc%Z_pair(jpair)%irreducible_shell = ish
                map%xuc%Z_pair(jpair)%operation_from_shell = jop ! Note, points to triplet operation
                map%xuc%Z_pair(jpair)%i1 = a1
                map%xuc%Z_pair(jpair)%i2 = a2
                map%xuc%Z_pair(jpair)%r = r0
                map%xuc%Z_pair(jpair)%lv = v0
                map%xuc%Z_pair(jpair)%flv = l0
                if (lo_sqnorm(r0) .lt. lo_sqtol) then
                    map%xuc%Z_pair(jpair)%selfterm = .true.
                else
                    map%xuc%Z_pair(jpair)%selfterm = .false.
                end if
                ! and finally, make sure that the information in the shell makes sense
                ! and point to the right kind of pairs.
                map%Z_pair_shell(ish)%unfold_index(i) = jpair
                ! Don't know if useful to keep this, for now I set it to nothing so I don't use it the wrong way.
                map%Z_pair_shell(ish)%unfold_operation(i) = -1
            end do
            end do

            ! Sort things out for the supercell. I have made this a little confusing, but that
            ! is to make sure it's completely general. Think it works alright. Technically I have
            ! the proper list of supercell pairs already, but I want to generate it procedurally
            ! since then I can generalize to all manner of tuplets and interactions with minimal
            ! change in the code.
            call mem%allocate(di, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            do ipair = 1, map%xuc%n_Z_pair
                a1 = map%xuc%Z_pair(ipair)%i1
                di(a1) = di(a1) + 1
            end do
            call mem%allocate(dj, [maxval(di), uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dk, ss%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            dj = 0
            dk = 0
            do ipair = 1, map%xuc%n_Z_pair
                a1 = map%xuc%Z_pair(ipair)%i1
                di(a1) = di(a1) + 1
                dj(di(a1), a1) = ipair
            end do
            ! So, now di holds the number of pairs per unitcell atom, and dj holds
            ! the list of pairs per unitcell atom, with that I can easily construct
            ! the list of supercell pairs
            map%xss%n_Z_pair = 0
            do a1 = 1, map%n_atom_ss
                a2 = ss%info%index_in_unitcell(a1)
                dk(a1) = map%xss%n_Z_pair
                map%xss%n_Z_pair = map%xss%n_Z_pair + di(a2)
            end do
            allocate (map%xss%ind_Z_pair(4, map%xss%n_Z_pair))
            map%xss%ind_Z_pair = 0
            ! Fetch irreducible shell, operation and first atom.
            do a1 = 1, map%n_atom_ss
                if (mod(a1, mw%n) .ne. mw%r) cycle
                a2 = ss%info%index_in_unitcell(a1)
                do i = 1, di(a2)
                    ipair = dj(i, a2)
                    l = dk(a1) + i
                    map%xss%ind_Z_pair(1, l) = map%xuc%Z_pair(ipair)%irreducible_shell
                    map%xss%ind_Z_pair(2, l) = map%xuc%Z_pair(ipair)%operation_from_shell
                    map%xss%ind_Z_pair(3, l) = a1
                    map%xss%ind_Z_pair(4, l) = -ipair
                end do
            end do
            call mw%allreduce('sum', map%xss%ind_Z_pair)
            ! Get the final index for each supercell pair
            do l = 1, map%xss%n_Z_pair
                if (mod(l, mw%n) .ne. mw%r) cycle
                a1 = map%xss%ind_Z_pair(3, l)
                ipair = -map%xss%ind_Z_pair(4, l)
                v0 = ss%rcart(:, a1) + map%xuc%Z_pair(ipair)%r
                v0 = matmul(ss%inv_latticevectors, v0)
                v0 = lo_clean_fractional_coordinates(v0)
                a2 = 0
                do i = 1, ss%na
                    r0 = lo_clean_fractional_coordinates(ss%r(:, i) - v0 + 0.5_r8) - 0.5_r8
                    if (lo_sqnorm(r0) .lt. lo_sqtol) then
                        a2 = i
                        exit
                    end if
                end do
                if (a2 .gt. 0) then
                    map%xss%ind_Z_pair(4, l) = a2
                else
                    call lo_stop_gracefully(['Could not locate pair'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                end if
            end do
            call mw%allreduce('max', map%xss%ind_Z_pair)
            ! And that should be it for the pairs!
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... built pair forcemap (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
            ! Cleanup
            call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dk, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ind_to_triplet_op, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block Zpair
    end if

    ! Borncharge pairs
    if (st%nx_Z_triplet .gt. 0) then
        Ztriplet: block
            real(r8), dimension(3) :: r2, r3, v2, v3, l2, l3
            integer, dimension(:, :), allocatable :: dj
            integer, dimension(:), allocatable :: di, dk, ind_to_quartet_op
            integer :: ish, iop, jop, itriplet, jtriplet, a1, a2, a3, i, j, k, l, ctr

            ! make a note that we have triplets, and store some basic information
            map%have_Z_triplet = .true.              ! we have Z triplets
            map%n_Z_triplet_shell = st%n_Z_triplet_shell  ! number of shells
            map%xuc%nx_Z_triplet = st%nx_Z_triplet       ! number of irreducible
            allocate (map%xuc%x_Z_triplet(map%xuc%nx_Z_triplet))
            map%xuc%x_Z_triplet = 0.0_r8              ! initialize irreducible to nothing

            ! I need to know how to translate from triplet operations to triplet operations. On entry, all the
            ! indices to operations point to triplet operations, I need to redirect those indices to triplet
            ! operations. This is confusing, I know.
            call mem%allocate(ind_to_quartet_op, st%n_quartet_operation, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            ind_to_quartet_op = 0
            do i = 1, st%n_triplet_operation
                k = 0
                ctr = 0
                do j = 1, st%n_quartet_operation
                    ! have to be the same rotation
                    if (st%quartetop(j)%opind .ne. st%tripletop(i)%opind) cycle
                    ! can not permute the electric field index
                    if (st%quartetop(j)%perm(1) .ne. 1) cycle
                    ! now it gets annoying. If the triplet permutation is
                    ! [1,2,3], then the quartet operation must be [2,3,4]. I think.
                    if (sum(abs(st%tripletop(i)%perm - (st%quartetop(j)%perm(2:4) - 1))) .eq. 0) then
                        ctr = ctr + 1
                        k = j
                    end if
                end do

                if (ctr .gt. 1) call lo_stop_gracefully(['More than one operation matches, impossible'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                if (k .eq. 0) then
                    call lo_stop_gracefully(['Could not find operation, impossible'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                else
                    ind_to_quartet_op(i) = k
                end if
            end do

            ! Populate the shells, can just do a straight copy f2008 style
            call copyshells(st%Z_triplet_shell, map%Z_triplet_shell)

            ! Count triplets?
            map%xuc%n_Z_triplet = 0
            do ish = 1, map%n_Z_triplet_shell
                map%xuc%n_Z_triplet = map%xuc%n_Z_triplet + map%Z_triplet_shell(ish)%n_unfold
            end do
            ! Space for unitcell triplets
            allocate (map%xuc%Z_triplet(map%xuc%n_Z_triplet))
            ! Store unitcell triplets
            jtriplet = 0
            do ish = 1, map%n_Z_triplet_shell
            do i = 1, map%Z_triplet_shell(ish)%n_unfold
                jtriplet = jtriplet + 1
                itriplet = st%Z_triplet_shell(ish)%unfold_index(i) ! triplet index in all triplets
                iop = st%Z_triplet_shell(ish)%unfold_operation(i)  ! Operation from shell, points to triplet operation
                jop = ind_to_quartet_op(iop)                       ! Operation from shell, points to quartet operation
                a1 = st%uc_triplet(itriplet)%i1                    ! index in unit cell to atom 1
                a2 = st%uc_triplet(itriplet)%i2                    ! index in unit cell to atom 2
                a3 = st%uc_triplet(itriplet)%i3                    ! index in unit cell to atom 3

                r2 = st%uc_triplet(itriplet)%v2                   ! Cartesian vector
                r3 = st%uc_triplet(itriplet)%v3                   ! Cartesian vector
                v2 = r2 - uc%rcart(:, a2) + uc%rcart(:, a1)             ! Cartesian lattice vector
                v3 = r3 - uc%rcart(:, a3) + uc%rcart(:, a1)             ! Cartesian lattice vector
                l2 = anint(matmul(uc%inv_latticevectors, v2))      ! Fractional lattice vector
                l3 = anint(matmul(uc%inv_latticevectors, v3))      ! Fractional lattice vector
                v2 = matmul(uc%latticevectors, l2)                 ! Cartesian lattice vector, clean
                v3 = matmul(uc%latticevectors, l3)                 ! Cartesian lattice vector, clean

                map%xuc%Z_triplet(jtriplet)%irreducible_shell = ish
                map%xuc%Z_triplet(jtriplet)%operation_from_shell = jop ! Note, points to quartet operation
                map%xuc%Z_triplet(jtriplet)%i1 = a1
                map%xuc%Z_triplet(jtriplet)%i2 = a2
                map%xuc%Z_triplet(jtriplet)%i3 = a3
                map%xuc%Z_triplet(jtriplet)%v2 = r2
                map%xuc%Z_triplet(jtriplet)%v3 = r3
                map%xuc%Z_triplet(jtriplet)%lv2 = v2
                map%xuc%Z_triplet(jtriplet)%lv3 = v3
                map%xuc%Z_triplet(jtriplet)%flv2 = l2
                map%xuc%Z_triplet(jtriplet)%flv3 = l3
                if (lo_sqnorm(r2) + lo_sqnorm(r3) .lt. lo_sqtol) then
                    map%xuc%Z_triplet(jtriplet)%selfterm = .true.
                else
                    map%xuc%Z_triplet(jtriplet)%selfterm = .false.
                end if
                ! and finally, make sure that the information in the shell makes sense
                ! and point to the right kind of triplets.
                map%Z_triplet_shell(ish)%unfold_index(i) = jtriplet
                ! Don't know if useful to keep this, for now I set it to nothing so I don't use it the wrong way.
                map%Z_triplet_shell(ish)%unfold_operation(i) = -1
            end do
            end do

            ! Sort things out for the supercell. I have made this a little confusing, but that
            ! is to make sure it's completely general. Think it works alright. Technically I have
            ! the proper list of supercell triplets already, but I want to generate it procedurally
            ! since then I can generalize to all manner of tuplets and interactions with minimal
            ! change in the code.
            call mem%allocate(di, uc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            do itriplet = 1, map%xuc%n_Z_triplet
                a1 = map%xuc%Z_triplet(itriplet)%i1
                di(a1) = di(a1) + 1
            end do
            call mem%allocate(dj, [maxval(di), uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dk, ss%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            di = 0
            dj = 0
            dk = 0
            do itriplet = 1, map%xuc%n_Z_triplet
                a1 = map%xuc%Z_triplet(itriplet)%i1
                di(a1) = di(a1) + 1
                dj(di(a1), a1) = itriplet
            end do
            ! So, now di holds the number of triplets per unitcell atom, and dj holds
            ! the list of triplets per unitcell atom, with that I can easily construct
            ! the list of supercell triplets
            map%xss%n_Z_triplet = 0
            do a1 = 1, map%n_atom_ss
                a2 = ss%info%index_in_unitcell(a1)
                dk(a1) = map%xss%n_Z_triplet
                map%xss%n_Z_triplet = map%xss%n_Z_triplet + di(a2)
            end do
            allocate (map%xss%ind_Z_triplet(5, map%xss%n_Z_triplet))
            map%xss%ind_Z_triplet = 0
            ! Fetch irreducible shell, operation and first atom.
            do a1 = 1, map%n_atom_ss
                if (mod(a1, mw%n) .ne. mw%r) cycle
                a2 = ss%info%index_in_unitcell(a1)
                do i = 1, di(a2)
                    itriplet = dj(i, a2)
                    l = dk(a1) + i
                    map%xss%ind_Z_triplet(1, l) = map%xuc%Z_triplet(itriplet)%irreducible_shell
                    map%xss%ind_Z_triplet(2, l) = map%xuc%Z_triplet(itriplet)%operation_from_shell
                    map%xss%ind_Z_triplet(3, l) = a1
                    map%xss%ind_Z_triplet(4, l) = -itriplet
                    map%xss%ind_Z_triplet(5, l) = -itriplet
                end do
            end do
            call mw%allreduce('sum', map%xss%ind_Z_triplet)
            ! Get the final index for each supercell triplet
            do l = 1, map%xss%n_Z_triplet
                if (mod(l, mw%n) .ne. mw%r) cycle
                a1 = map%xss%ind_Z_triplet(3, l)
                itriplet = -map%xss%ind_Z_triplet(4, l)
                v2 = ss%rcart(:, a1) + map%xuc%Z_triplet(itriplet)%v2
                v3 = ss%rcart(:, a1) + map%xuc%Z_triplet(itriplet)%v3
                v2 = matmul(ss%inv_latticevectors, v2)
                v3 = matmul(ss%inv_latticevectors, v3)
                v2 = lo_clean_fractional_coordinates(v2)
                v3 = lo_clean_fractional_coordinates(v3)
                a2 = 0
                a3 = 0
                j = 0
                do i = 1, ss%na
                    r2 = lo_clean_fractional_coordinates(ss%r(:, i) - v2 + 0.5_r8) - 0.5_r8
                    r3 = lo_clean_fractional_coordinates(ss%r(:, i) - v3 + 0.5_r8) - 0.5_r8
                    if (lo_sqnorm(r2) .lt. lo_sqtol) then
                        a2 = i
                        j = j + 1
                    end if
                    if (lo_sqnorm(r3) .lt. lo_sqtol) then
                        a3 = i
                        j = j + 1
                    end if
                    if (j .eq. 2) exit
                end do
                if (min(a2, a3) .gt. 0) then
                    map%xss%ind_Z_triplet(4, l) = a2
                    map%xss%ind_Z_triplet(5, l) = a3
                else
                    call lo_stop_gracefully(['Could not locate triplet'], lo_exitcode_symmetry, __FILE__, __LINE__, mw%comm)
                end if
            end do
            call mw%allreduce('max', map%xss%ind_Z_triplet)
            ! And sort out the ASR tuplets
            call Ztriplet_asr_tuplets(map, mem)
            ! And that should be it for the triplets!
            if (verbosity .gt. 0) then
                t1 = walltime()
                write (lo_iou, *) '... built triplet forcemap (', tochar(t1 - t0), 's)'
                t0 = t1
            end if
            ! Cleanup
            call mem%deallocate(di, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dk, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ind_to_quartet_op, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block Ztriplet
    end if

    ! and we are done!
    if (verbosity .gt. 0) write (*, *) '... finished creating forcemap (', tochar(walltime() - t0), 's)'
end subroutine

! Not exposed below

!> get the thirdorder tuplets that should be summed to zero!
subroutine thirdorder_asr_tuplets(map, mem)
    !> symmetry stuff rearranged into coordination shells
    class(lo_forcemap), intent(inout) :: map
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:, :), allocatable :: dumv2
    integer, dimension(:), allocatable :: dj, dk, dl
    integer :: ngroup, a1, itrip, i, j, ntrip

    ngroup = 0
    ! First just count, I think.
    atomloop1: do a1 = 1, map%n_atom_uc
        ! Count triplets
        ntrip = 0
        do itrip = 1, map%xuc%n_fc_triplet
            if (map%xuc%fc_triplet(itrip)%i1 .ne. a1) cycle
            ntrip = ntrip + 1
        end do

        if (ntrip .gt. 0) then
            call mem%allocate(dj, ntrip, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dl, ntrip, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dumv2, [3, ntrip], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dj = 0
            dumv2 = 0.0_r8
        else
            call lo_stop_gracefully(['No triplets found'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if

        ! Count triplets
        i = 0
        do itrip = 1, map%xuc%n_fc_triplet
            if (map%xuc%fc_triplet(itrip)%i1 .ne. a1) cycle
            i = i + 1
            dj(i) = itrip
            dumv2(:, i) = map%xuc%fc_triplet(itrip)%v2
        end do
        ! Get the unique
        call lo_return_unique_indices(dumv2, dk, mem, redind=dl, tol=lo_tol)
        ! Count number of groups
        ngroup = ngroup + size(dk)
        ! And cleanup for the next atom
        call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dk, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dl, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dumv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end do atomloop1
    ! Make space for groups
    allocate (map%xuc%fc_triplet_group(ngroup))
    do i = 1, ngroup
        map%xuc%fc_triplet_group(i)%n = 0
        map%xuc%fc_triplet_group(i)%contains_selfterm = .false.
    end do
    ! Now do the same thing again and store the information.
    ngroup = 0
    atomloop2: do a1 = 1, map%n_atom_uc
        ! Count triplets
        ntrip = 0
        do itrip = 1, map%xuc%n_fc_triplet
            if (map%xuc%fc_triplet(itrip)%i1 .ne. a1) cycle
            ntrip = ntrip + 1
        end do
        if (ntrip .gt. 0) then
            call mem%allocate(dj, ntrip, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dl, ntrip, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dumv2, [3, ntrip], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dj = 0
            dumv2 = 0.0_r8
        else
            call lo_stop_gracefully(['No triplets found'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
        i = 0
        do itrip = 1, map%xuc%n_fc_triplet
            if (map%xuc%fc_triplet(itrip)%i1 .ne. a1) cycle
            i = i + 1
            dj(i) = itrip
            dumv2(:, i) = map%xuc%fc_triplet(itrip)%v2
        end do
        ! Get the unique
        call lo_return_unique_indices(dumv2, dk, mem, redind=dl, tol=lo_tol)
        ! Count members
        do i = 1, ntrip
            j = dl(i) + ngroup
            map%xuc%fc_triplet_group(j)%n = map%xuc%fc_triplet_group(j)%n + 1
        end do
        do i = 1, size(dk)
            j = ngroup + i
            allocate (map%xuc%fc_triplet_group(j)%ind(map%xuc%fc_triplet_group(j)%n))
            map%xuc%fc_triplet_group(j)%ind = 0
            map%xuc%fc_triplet_group(j)%n = 0
        end do
        ! Store members
        do i = 1, ntrip
            j = dl(i) + ngroup
            map%xuc%fc_triplet_group(j)%n = map%xuc%fc_triplet_group(j)%n + 1
            map%xuc%fc_triplet_group(j)%ind(map%xuc%fc_triplet_group(j)%n) = dj(i)
        end do
        ! Check if any of the groups contains a selfterm
        do i = 1, size(dk)
            if (lo_sqnorm(dumv2(:, dk(i))) .lt. lo_sqtol) then
                j = ngroup + i
                map%xuc%fc_triplet_group(j)%contains_selfterm = .true.
            end if
        end do
        ! And note how many groups we are at.
        ngroup = ngroup + size(dk)
        ! And cleanup for the next atom
        call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dk, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dl, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dumv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end do atomloop2
end subroutine

!> get the fourthorder tuplets that should be summed to zero!
subroutine fourthorder_asr_tuplets(map, mem)
    !> symmetry stuff rearranged into coordination shells
    class(lo_forcemap), intent(inout) :: map
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:, :), allocatable :: dumv2
    integer, dimension(:), allocatable :: dj, dk, dl
    integer :: ngroup, a1, iquart, i, j, k, l, nquart, igroup

    ngroup = 0
    ! First just count, I think.
    atomloop1: do a1 = 1, map%n_atom_uc
        ! Count quartets
        nquart = 0
        do iquart = 1, map%xuc%n_fc_quartet
            if (map%xuc%fc_quartet(iquart)%i1 .ne. a1) cycle
            nquart = nquart + 1
        end do

        if (nquart .gt. 0) then
            call mem%allocate(dj, nquart, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dl, nquart, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dumv2, [6, nquart], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dj = 0
            dumv2 = 0.0_r8
        else
            call lo_stop_gracefully(['No quartets found'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if

        ! Count quartets
        i = 0
        do iquart = 1, map%xuc%n_fc_quartet
            if (map%xuc%fc_quartet(iquart)%i1 .ne. a1) cycle
            i = i + 1
            dj(i) = iquart
            dumv2(1:3, i) = map%xuc%fc_quartet(iquart)%v2
            dumv2(4:6, i) = map%xuc%fc_quartet(iquart)%v3
        end do
        ! Get the unique
        call lo_return_unique_indices(dumv2, dk, mem, redind=dl, tol=lo_tol)
        ! Count number of groups
        ngroup = ngroup + size(dk)
        ! And cleanup for the next atom
        call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dk, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dl, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dumv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end do atomloop1
    ! Make space for groups
    allocate (map%xuc%fc_quartet_group(ngroup))
    do i = 1, ngroup
        map%xuc%fc_quartet_group(i)%n = 0
        map%xuc%fc_quartet_group(i)%contains_selfterm = .false.
    end do
    ! Now do the same thing again and store the information.
    ngroup = 0
    atomloop2: do a1 = 1, map%n_atom_uc
        ! Count quartets
        nquart = 0
        do iquart = 1, map%xuc%n_fc_quartet
            if (map%xuc%fc_quartet(iquart)%i1 .ne. a1) cycle
            nquart = nquart + 1
        end do
        if (nquart .gt. 0) then
            call mem%allocate(dj, nquart, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dl, nquart, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dumv2, [6, nquart], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dj = 0
            dl = 0
            dumv2 = 0.0_r8
        else
            call lo_stop_gracefully(['No quartets found'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
        i = 0
        do iquart = 1, map%xuc%n_fc_quartet
            if (map%xuc%fc_quartet(iquart)%i1 .ne. a1) cycle
            i = i + 1
            dj(i) = iquart
            dumv2(1:3, i) = map%xuc%fc_quartet(iquart)%v2
            dumv2(4:6, i) = map%xuc%fc_quartet(iquart)%v3
        end do
        ! Get the unique
        call lo_return_unique_indices(dumv2, dk, mem, redind=dl, tol=lo_tol)
        ! Count members
        do i = 1, nquart
            j = dl(i) + ngroup
            map%xuc%fc_quartet_group(j)%n = map%xuc%fc_quartet_group(j)%n + 1
        end do
        do i = 1, size(dk)
            j = ngroup + i
            allocate (map%xuc%fc_quartet_group(j)%ind(map%xuc%fc_quartet_group(j)%n))
            map%xuc%fc_quartet_group(j)%ind = 0
            map%xuc%fc_quartet_group(j)%n = 0
        end do
        ! Store members
        do i = 1, nquart
            j = dl(i) + ngroup
            map%xuc%fc_quartet_group(j)%n = map%xuc%fc_quartet_group(j)%n + 1
            map%xuc%fc_quartet_group(j)%ind(map%xuc%fc_quartet_group(j)%n) = dj(i)
        end do
        ! And note how many groups we are at.
        ngroup = ngroup + size(dk)
        ! And cleanup for the next atom
        call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dk, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dl, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dumv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end do atomloop2

    ! Check for selfterms
    j = 0
    do igroup = 1, ngroup
        j = max(j, map%xuc%fc_quartet_group(igroup)%n)
    end do
    call mem%allocate(dj, j, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    dj = 0
    grloop: do igroup = 1, ngroup
        k = -1
        l = 0
        dj = -1
        do i = 1, map%xuc%fc_quartet_group(igroup)%n
            iquart = map%xuc%fc_quartet_group(igroup)%ind(i)
            if (map%xuc%fc_quartet(iquart)%selfterm) then
                k = iquart
            else
                l = l + 1
                dj(l) = iquart
            end if
        end do
        if (k .gt. 0) then
            ! Found the group with the self-term
            map%xuc%fc_quartet_group(igroup)%n = l
            map%xuc%fc_quartet_group(igroup)%ind = -1
            map%xuc%fc_quartet_group(igroup)%ind(1:l) = dj(1:l)
            map%xuc%fc_quartet_group(igroup)%contains_selfterm = .true.
        end if
    end do grloop

    call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> get the thirdorder tuplets that should be summed to zero!
subroutine Ztriplet_asr_tuplets(map, mem)
    !> symmetry stuff rearranged into coordination shells
    class(lo_forcemap), intent(inout) :: map
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:, :), allocatable :: dumv2
    integer, dimension(:), allocatable :: dj, dk, dl
    integer :: ngroup, a1, itrip, i, j, ntrip

    ngroup = 0
    ! First just count, I think.
    atomloop1: do a1 = 1, map%n_atom_uc
        ! Count triplets
        ntrip = 0
        do itrip = 1, map%xuc%n_Z_triplet
            if (map%xuc%Z_triplet(itrip)%i1 .ne. a1) cycle
            ntrip = ntrip + 1
        end do

        if (ntrip .gt. 0) then
            call mem%allocate(dj, ntrip, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dl, ntrip, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dumv2, [3, ntrip], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dj = 0
            dumv2 = 0.0_r8
        else
            call lo_stop_gracefully(['No triplets found'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if

        ! Count triplets
        i = 0
        do itrip = 1, map%xuc%n_Z_triplet
            if (map%xuc%Z_triplet(itrip)%i1 .ne. a1) cycle
            i = i + 1
            dj(i) = itrip
            dumv2(:, i) = map%xuc%Z_triplet(itrip)%v2
        end do
        ! Get the unique
        call lo_return_unique_indices(dumv2, dk, mem, redind=dl, tol=lo_tol)
        ! Count number of groups
        ngroup = ngroup + size(dk)
        ! And cleanup for the next atom
        call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dk, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dl, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dumv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end do atomloop1
    ! Make space for groups
    allocate (map%xuc%Z_triplet_group(ngroup))
    do i = 1, ngroup
        map%xuc%Z_triplet_group(i)%n = 0
        map%xuc%Z_triplet_group(i)%contains_selfterm = .false.
    end do
    ! Now do the same thing again and store the information.
    ngroup = 0
    atomloop2: do a1 = 1, map%n_atom_uc
        ! Count triplets
        ntrip = 0
        do itrip = 1, map%xuc%n_Z_triplet
            if (map%xuc%Z_triplet(itrip)%i1 .ne. a1) cycle
            ntrip = ntrip + 1
        end do
        if (ntrip .gt. 0) then
            call mem%allocate(dj, ntrip, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dl, ntrip, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dumv2, [3, ntrip], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            dj = 0
            dumv2 = 0.0_r8
        else
            call lo_stop_gracefully(['No triplets found'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
        i = 0
        do itrip = 1, map%xuc%n_Z_triplet
            if (map%xuc%Z_triplet(itrip)%i1 .ne. a1) cycle
            i = i + 1
            dj(i) = itrip
            dumv2(:, i) = map%xuc%Z_triplet(itrip)%v2
        end do
        ! Get the unique
        call lo_return_unique_indices(dumv2, dk, mem, redind=dl, tol=lo_tol)
        ! Count members
        do i = 1, ntrip
            j = dl(i) + ngroup
            map%xuc%Z_triplet_group(j)%n = map%xuc%Z_triplet_group(j)%n + 1
        end do
        do i = 1, size(dk)
            j = ngroup + i
            allocate (map%xuc%Z_triplet_group(j)%ind(map%xuc%Z_triplet_group(j)%n))
            map%xuc%Z_triplet_group(j)%ind = 0
            map%xuc%Z_triplet_group(j)%n = 0
        end do
        ! Store members
        do i = 1, ntrip
            j = dl(i) + ngroup
            map%xuc%Z_triplet_group(j)%n = map%xuc%Z_triplet_group(j)%n + 1
            map%xuc%Z_triplet_group(j)%ind(map%xuc%Z_triplet_group(j)%n) = dj(i)
        end do
        ! Check if any of the groups contains a selfterm
        do i = 1, size(dk)
            if (lo_sqnorm(dumv2(:, dk(i))) .lt. lo_sqtol) then
                j = ngroup + i
                map%xuc%Z_triplet_group(j)%contains_selfterm = .true.
            end if
        end do
        ! And note how many groups we are at.
        ngroup = ngroup + size(dk)
        ! And cleanup for the next atom
        call mem%deallocate(dj, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dk, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dl, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dumv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end do atomloop2
end subroutine

! ifort was not clever enough to copy properly so I have to do it manually
subroutine copyshells(shellin, shellout)
    !> original shells
    type(lo_tensor_shell), dimension(:), intent(in) :: shellin
    !> copy of shells
    type(lo_tensor_shell), dimension(:), allocatable, intent(out) :: shellout

    integer :: ish

    allocate (shellout(size(shellin)))
    do ish = 1, size(shellin)
        shellout(ish)%nx = shellin(ish)%nx
        shellout(ish)%n_unfold = shellin(ish)%n_unfold
        if (shellout(ish)%nx .gt. 0) then
            allocate (shellout(ish)%ind_local(shellout(ish)%nx))
            allocate (shellout(ish)%ind_global(shellout(ish)%nx))
            allocate (shellout(ish)%coeff(size(shellin(ish)%coeff, 1), shellout(ish)%nx))
            shellout(ish)%ind_local = shellin(ish)%ind_local
            shellout(ish)%ind_global = shellin(ish)%ind_global
            shellout(ish)%coeff = shellin(ish)%coeff
        end if
        if (shellout(ish)%n_unfold .gt. 0) then
            allocate (shellout(ish)%unfold_index(shellout(ish)%n_unfold))
            allocate (shellout(ish)%unfold_operation(shellout(ish)%n_unfold))
            shellout(ish)%unfold_index = shellin(ish)%unfold_index
            shellout(ish)%unfold_operation = shellin(ish)%unfold_operation
        end if
    end do
end subroutine

end submodule
