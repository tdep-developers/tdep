module refine
!! Routines to clean up a crystal structure
use konstanter, only: r8, lo_sqtol, lo_tol, lo_huge, lo_exitcode_symmetry, lo_iou, lo_hugeint
use gottochblandat, only: lo_determ, lo_fetch_tolerance, lo_clean_fractional_coordinates, lo_chop, lo_return_unique, &
                          qsort, lo_invert3x3matrix, lo_sqnorm, lo_stop_gracefully, lo_find_rotation_that_makes_strain_diagonal
use type_crystalstructure, only: lo_crystalstructure
use type_symmetryoperation, only: lo_symset, lo_operate_on_vector
use type_distancetable, only: lo_distancetable
use geometryfunctions, only: lo_bounding_sphere_of_box
use lo_memtracker, only: lo_mem_helper
use lo_spacegroup, only: lo_symmetry_group
use lo_sorting, only: lo_qsort, lo_return_unique_indices
implicit none

private
! public :: refine_unitcell
public :: refine_supercell
! public :: refine_lattice
! public :: find_true_unitcell
public :: refine_one_cell

type lo_speciesbox_box
    !> how many of this species
    integer :: n = -lo_hugeint
    !> the points
    real(r8), dimension(:, :), allocatable :: r
    !> the original index of these points
    integer, dimension(:), allocatable :: ind
end type

type lo_speciesbox
    !> how many species
    integer :: n_species = -lo_hugeint
    !> how many atoms, in total
    integer :: n_atom = -lo_hugeint
    !> one box per species
    type(lo_speciesbox_box), dimension(:), allocatable :: species
    !> workspace for transformed positions
    real(r8), dimension(:, :), allocatable, private :: tf
contains
    !> sort atoms into boxes
    procedure :: generate => sort_atoms_into_boxes
    !> mismatch after operations
    procedure :: operation_error => box_operation_error
end type

contains

!> find mismatch and possibly rectify it?
subroutine refine_one_cell(uc, p, tolerance_lattice, tolerance_internal, mem, verbosity)
    !> input structure, possibly ugly
    type(lo_crystalstructure), intent(in) :: uc
    !> output structure
    type(lo_crystalstructure), intent(out) :: p
    !> tolerance
    real(r8), intent(in) :: tolerance_lattice
    real(r8), intent(in) :: tolerance_internal
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_symmetry_group) :: sg
    real(r8), parameter :: tighttolerance = 1E-12_r8
    real(r8), dimension(:, :, :), allocatable :: op
    real(r8), dimension(:, :), allocatable :: tr, positions
    real(r8), dimension(3, 3) :: refinedbasis

    ! Fetch some operations to work with?
    init: block
        if (verbosity .gt. 0) then
            write (lo_iou, *) ''
            write (lo_iou, *) 'REFINING STRUCTURE'
        end if

    end block init

    ! Start by checking the lattice perhaps
    lattice: block
        real(r8), dimension(3, 3, 24) :: hexops
        real(r8), dimension(3, 3, 48) :: octops
        real(r8), dimension(3, 3, 48) :: valops
        real(r8), dimension(:, :, :), allocatable :: bufop
        real(r8), dimension(:, :), allocatable :: cluster, rotcluster, symcluster
        real(r8), dimension(:), allocatable :: symerr
        real(r8), dimension(3, 3) :: m0, im0, m1, m2, m3
        real(r8) :: f0, f1, f2
        integer, dimension(:), allocatable :: di1, di2
        integer, dimension(3) :: lvind
        integer :: i, j, k, l, iter, ctr, nc, iop, initer, errcounter

        ! Get some operations? This is for the initial stuff?
        call sg%generate_pool_of_point_operations(uc%latticevectors, hexops, octops, f0, f1, mem, verbosity)
        if (f0 .lt. f1) then
            allocate (op(3, 3, 24))
            op = hexops
        else
            allocate (op(3, 3, 48))
            op = octops
        end if

        ! Get a cluster thing? Think that makes reasonable sense.
        call cluster_from_lattice_vectors(uc%latticevectors, cluster, lvind)

        ! Calculate mismatch?
        allocate (symerr(size(op, 3)))
        allocate (di1(size(op, 3)))
        di1 = 0
        symerr = 0.0_r8
        do iop = 1, size(op, 3)
            symerr(iop) = operation_transform_error(op(:, :, iop), cluster)
        end do
        call lo_qsort(symerr, di1)

        ! Report mismatch
        ctr = 0
        l = 0
        do i = 1, size(op, 3)
            iop = di1(i)
            if (symerr(i) .lt. tolerance_lattice) then
                ctr = ctr + 1
                write (*, *) 'op', iop, 'err:', symerr(i)
            else
                l = l + 1
                if (l .eq. 1) then
                    write (*, *) '----- tol:', tolerance_lattice
                end if
                write (*, *) 'op', iop, 'err:', symerr(i)
                if (l .ge. 10) exit
            end if
        end do

        ! Try to refine this a little bit?
        errcounter = 0
        m0 = uc%latticevectors
        iterloop: do iter = 1, 10
            ! Cleanup, just in case
            if (allocated(op)) deallocate (op)
            if (allocated(cluster)) deallocate (cluster)
            if (allocated(rotcluster)) deallocate (rotcluster)
            if (allocated(symcluster)) deallocate (symcluster)
            if (allocated(di1)) deallocate (di1)
            if (allocated(di2)) deallocate (di2)
            if (allocated(symerr)) deallocate (symerr)

            write (*, *) ''
            write (*, *) 'iteration:', iter

            ! Re-generate operations?
            call sg%generate_pool_of_point_operations(uc%latticevectors, hexops, octops, f0, f1, mem, verbosity - 1)
            if (f0 .lt. f1) then
                allocate (op(3, 3, 24))
                op = hexops
                allocate (symerr(24))
                allocate (di1(24))
            else
                allocate (op(3, 3, 48))
                op = octops
                allocate (symerr(48))
                allocate (di1(48))
            end if

            ! Re-generate cluster
            call cluster_from_lattice_vectors(m0, cluster, lvind)
            nc = size(cluster, 2)
            allocate (rotcluster(3, nc))
            allocate (symcluster(3, nc))
            allocate (di2(nc))

            ! Find valid subset?
            valops = 0.0_r8
            ctr = 0
            do iop = 1, size(op, 3)
                f0 = operation_transform_error(op(:, :, iop), cluster)
                if (f0 .lt. tolerance_lattice) then
                    ctr = ctr + 1
                    valops(:, :, ctr) = op(:, :, iop)
                end if
            end do
            call expand_valid_operations(valops, ctr, tolerance_lattice)

            f0 = cluster_transform_error(valops(:, :, 1:ctr), cluster)
            if (f0 .lt. tighttolerance) then
                ! We are done!
                exit iterloop
            end if

            ! Now we have a set of valid operations. Make sure
            ! the cluster satisfies these operations?
            di2 = 0
            symcluster = 0.0_r8
            do iop = 1, ctr
                ! Rotate the cluster
                rotcluster = matmul(valops(:, :, iop), cluster)
                do i = 1, nc
                    k = -1
                    f1 = lo_huge
                    do j = 1, nc
                        f0 = norm2(cluster(:, i) - rotcluster(:, j))
                        if (f0 .lt. f1) then
                            f1 = f0
                            k = j
                        end if
                    end do
                    symcluster(:, i) = symcluster(:, i) + rotcluster(:, k)
                    di2(i) = di2(i) + 1
                end do
            end do
            do i = 1, nc
                cluster(:, i) = symcluster(:, i)/real(di2(i), r8)
            end do

            ! Now update the lattice vectors?
            do i = 1, 3
                j = lvind(i)
                m0(:, i) = cluster(:, j)
            end do
            m0 = lo_chop(m0, 1E-13_r8)
        end do iterloop

        ! Store the symmetry-refined lattice vectors.
        refinedbasis = m0

        ! Now we are happy, report the error?
        di1 = 0
        symerr = 0.0_r8
        do iop = 1, size(op, 3)
            symerr(iop) = operation_transform_error(op(:, :, iop), cluster)
        end do
        call lo_qsort(symerr, di1)

        ! Report mismatch
        ctr = 0
        l = 0
        do i = 1, size(op, 3)
            iop = di1(i)
            if (symerr(i) .lt. tolerance_lattice) then
                ctr = ctr + 1
                write (*, *) 'op', iop, 'err:', symerr(i)
            else
                l = l + 1
                if (l .eq. 1) then
                    write (*, *) '-----'
                end if
                write (*, *) 'op', iop, 'err:', symerr(i)
                if (l .ge. 10) exit
            end if
        end do

        write (*, *) 'old basis:'
        do i = 1, 3
            write (*, *) uc%latticevectors(:, i)
        end do
        write (*, *) 'refined basis:'
        do i = 1, 3
            write (*, *) refinedbasis(:, i)
        end do

        ! Make sure the point operations are with respect to
        ! the refined basis, and only the valid operations?
        call sg%generate_pool_of_point_operations(uc%latticevectors, hexops, octops, f0, f1, mem, verbosity)
        deallocate (op)
        ctr = 0
        do iop = 1, 24
            f0 = operation_transform_error(hexops(:, :, iop), cluster)
            if (f0 .lt. tolerance_lattice) ctr = ctr + 1
        end do
        do iop = 1, 48
            f0 = operation_transform_error(octops(:, :, iop), cluster)
            if (f0 .lt. tolerance_lattice) ctr = ctr + 1
        end do
        allocate (op(3, 3, ctr))
        ctr = 0
        do iop = 1, 24
            f0 = operation_transform_error(hexops(:, :, iop), cluster)
            if (f0 .lt. tolerance_lattice) then
                ctr = ctr + 1
                op(:, :, ctr) = hexops(:, :, iop)
            end if
        end do
        do iop = 1, 48
            f0 = operation_transform_error(octops(:, :, iop), cluster)
            if (f0 .lt. tolerance_lattice) then
                ctr = ctr + 1
                op(:, :, ctr) = octops(:, :, iop)
            end if
        end do
        ! Now double check that this forms a group!
        call sg%ensure_point_operations_are_a_group(op, tolerance_lattice)
    end block lattice

    ! Now that we have a neat refined basis, perhaps nudge it to a form
    ! that is more well-known?
    ! finaltouch: block
    ! end block finaltouch

    ! Fix the internal degrees of freedom.
    internal: block
        type(lo_speciesbox) :: box
        integer, parameter :: nmax = 4 ! max fractional translation?
        real(r8), dimension(:, :, :), allocatable :: fop, valop
        real(r8), dimension(:, :), allocatable :: r0, r1, r2, valtr
        real(r8), dimension(:), allocatable :: symerr, rv0
        real(r8), dimension(3, 3) :: inversebasis, m0
        real(r8) :: f0, f1
        integer, dimension(:, :), allocatable :: di1
        integer, dimension(:), allocatable :: di2
        integer :: na, iop, i, j, k, l, nop, iter

        write (*, *) ''
        write (*, *) 'Refining internal positions'

        na = uc%na
        allocate (r0(3, na))
        allocate (r1(3, na))
        allocate (r2(3, na))
        r0 = 0.0_r8
        r1 = 0.0_r8
        r2 = 0.0_r8
        inversebasis = lo_invert3x3matrix(refinedbasis)

        ! First, convert operations to fractional coordinates?
        allocate (fop(3, 3, size(op, 3)))
        fop = 0.0_r8
        do iop = 1, size(op, 3)
            m0 = op(:, :, iop)
            m0 = matmul(inversebasis, m0)
            m0 = matmul(m0, refinedbasis)
            if (norm2(m0 - anint(m0)) .lt. tighttolerance) then
                fop(:, :, iop) = anint(m0)
            else
                call lo_stop_gracefully(['failed generating point operations'], lo_exitcode_symmetry, __FILE__, __LINE__)
            end if
        end do

        ! Get some translations to work with?
        call sg%generate_bruteforce_pool_of_translations(tr, 4)

        ! Pair translations and point operations
        nop = size(tr, 2)*size(op, 3)
        allocate (di1(2, nop))
        l = 0
        do i = 1, size(op, 3)
        do j = 1, size(tr, 2)
            l = l + 1
            di1(:, l) = [i, j]
        end do
        end do
        allocate (symerr(nop))
        allocate (rv0(nop))

        ! Sort the atoms into chunks per species
        call box%generate(uc%r, uc%species, mem)

        iterloopinternal: do iter = 1, 50

            write (*, *) 'Refining internal positions iteration', iter

            ! Measure error
            l = 0
            f1 = 0.0_r8
            do iop = 1, nop
                i = di1(1, iop)
                j = di1(2, iop)
                f0 = box%operation_error(fop(:, :, i), tr(:, j))
                symerr(iop) = f0
                if (f0 .lt. tolerance_internal) then
                    l = l + 1
                    f1 = f1 + f0
                end if
            end do
            if (allocated(valop)) deallocate (valop)
            if (allocated(valtr)) deallocate (valtr)
            allocate (valop(3, 3, l))
            allocate (valtr(3, l))

            rv0 = symerr
            call lo_qsort(rv0)

            l = 0
            do iop = 1, nop
                if (rv0(iop) .gt. tolerance_internal) l = l + 1
                if (l .eq. 1) then
                    write (*, *) '======'
                end if
                write (*, *) 'op', iop, 'err', rv0(iop)
                if (l .ge. 15) exit
            end do

            ! Store valid operataions
            l = 0
            do iop = 1, nop
                if (symerr(iop) .gt. tolerance_internal) cycle
                l = l + 1
                i = di1(1, iop)
                j = di1(2, iop)
                valop(:, :, l) = fop(:, :, i)
                valtr(:, l) = tr(:, j)
            end do

            ! Now, we might be converged already!
            if (f1 .lt. tighttolerance) then
                exit
            end if

            ! Enforce symmetry
            call box_refine_positions(box, valop, valtr)

            ! Perhaps report a little
            if (iter .eq. 1) then
                call lo_qsort(symerr)
                l = 0
                do iop = 1, nop
                    if (symerr(iop) .lt. tolerance_internal) then
                        write (*, *) iop, symerr(iop)
                    else

                        l = l + 1
                        if (l .eq. 1) write (*, *) '----- tolerance:', tolerance_internal
                        write (*, *) iop, symerr(iop)
                        if (l .gt. 15) exit
                    end if
                end do
            end if
        end do iterloopinternal

        ! Store the refined positions?
        deallocate (r0)
        allocate (r0(3, uc%na))
        do i = 1, box%n_species
        do j = 1, box%species(i)%n
            l = box%species(i)%ind(j)
            r0(:, l) = box%species(i)%r(:, j)
        end do
        end do
        call lo_return_unique(r0, positions, tolerance_internal)
        allocate (di2(size(positions, 2)))
        di2 = 0
        do j = 1, size(di2)
            do i = 1, uc%na
                if (norm2(r0(:, i) - positions(:, j)) .lt. tighttolerance) then
                    di2(j) = uc%atomic_number(i)
                end if
            end do
        end do

        ! Update the structure
        call p%generate(refinedbasis, positions, di2, 2)

    end block internal

end subroutine

!> Refine the supercell
subroutine refine_supercell(uc, ss)
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss

    type(lo_crystalstructure) :: p
    real(r8), dimension(3, 3) :: m0
    real(r8), dimension(3) :: v0
    real(r8) :: f0
    integer, dimension(3, 3) :: ssdim
    integer :: i, j, l

    ! Get the supercell matrix
    m0 = matmul(uc%inv_latticevectors, ss%latticevectors)
    ssdim = int(anint(m0))
    ! Build proper supercell from the refined guy
    call uc%build_supercell(p, nondiagdimensions=ssdim)

    ! Replace the latticevectors
    ss%latticevectors = p%latticevectors
    ss%inv_latticevectors = p%inv_latticevectors
    ! Now return the positions, in the same order as they where in the original cell
    do i = 1, ss%na
        f0 = lo_huge
        l = -1
        do j = 1, p%na
            v0 = lo_clean_fractional_coordinates(p%r(:, j) - ss%r(:, i) + 0.5_r8) - 0.5_r8
            if (lo_sqnorm(v0) .lt. f0) then
                f0 = lo_sqnorm(v0)
                l = j
            end if
        end do
        if (l .lt. 0) then
            call lo_stop_gracefully(['Did thinking wrong when building supercell'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
        ss%r(:, i) = p%r(:, l)
    end do
    ! Make sure it actually works?
    call ss%classify('supercell', uc)
end subroutine

!> Clean unitcell
subroutine refine_unitcell(uc, sym)
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> symmetry operations
    type(lo_symset), intent(in) :: sym

    integer, parameter :: maxit = 100
    real(r8), parameter :: symtol = 1E-15_r8
    real(r8), dimension(:, :), allocatable :: r0, r1
    real(r8) :: f0
    integer :: i, j, o, iter

    allocate (r0(3, uc%na))
    allocate (r1(3, uc%na))
    r0 = uc%r
    r1 = 0.0_r8
    iter = 0
    ! Quick test to see if optimization is necessary
    f0 = symdist_internal(sym, uc%r)
    if (f0 .lt. symtol) then
        ! Nothing else to do here!
        write (*, "(1X,A,ES11.4)") '... internal symmetry satisfied to high precision: ', f0
        return
    end if

    write (*, *) '... cleaning up internal degrees of freedom'
    write (*, *) ' iteration      violation'

    do iter = 1, maxit
        f0 = symdist_internal(sym, r0)
        write (*, "(2X,I4,7X,ES15.8)") iter - 1, f0
        if (f0 .lt. symtol) then
            write (*, "(1X,A,ES11.4)") '... internal symmetry satisfied to high precision: ', f0
            uc%r = lo_chop(r0, 1E-12_r8)
            exit
        end if
        do o = 1, sym%n
            do i = 1, uc%na
                j = sym%op(o)%fmap(i)
                r1(:, j) = lo_operate_on_vector(sym%op(o), r0(:, i), fractional=.true.)
            end do
            r1 = lo_clean_fractional_coordinates(r1 - r0 + 0.5_r8, 1E-13_r8) - 0.5_r8
            r0 = lo_clean_fractional_coordinates(r0 + 0.25_r8*r1, 1E-13_r8)
        end do
    end do

end subroutine

!> distance from proper symmetry, internal degrees of freedom
pure function symdist_internal(sym, r) result(dst)
    type(lo_symset), intent(in) :: sym
    real(r8), dimension(:, :), intent(in) :: r
    real(r8) :: dst

    integer :: i, j, o, na
    real(r8), dimension(3) :: v0, v1, v2, v3

    na = size(r, 2)
    dst = 0.0_r8
    do o = 1, sym%n
    do i = 1, na
        v0 = r(:, i)
        j = sym%op(o)%fmap(i)
        v1 = lo_operate_on_vector(sym%op(o), v0, fractional=.true.)
        v2 = r(:, j)
        v3 = lo_clean_fractional_coordinates(v1 - v2 + 0.5_r8, 1E-13_r8) - 0.5_r8
        dst = dst + lo_sqnorm(v3)
    end do
    end do
    dst = sqrt(dst/na)
end function

!!> triple matrix product, m = s * n * s^T for 3x3 matrices
!pure function matmul_SNST(n,s) result(m)
!    !> matrix n
!    real(r8), dimension(3,3), intent(in) :: n
!    !> operation S
!    real(r8), dimension(3,3), intent(in) :: s
!    !> transformed matrix
!    real(r8), dimension(3,3) :: m
!
!    m(1,1)=n(1,1)*s(1,1)**2 + n(1,2)*s(1,1)*s(1,2) + n(2,1)*s(1,1)*s(1,2) + n(2,2)*s(1,2)**2 + n(1,3)*s(1,1)*s(1,3) + &
!           n(3,1)*s(1,1)*s(1,3) + n(2,3)*s(1,2)*s(1,3) + n(3,2)*s(1,2)*s(1,3) + n(3,3)*s(1,3)**2
!    m(1,2)=n(1,1)*s(1,1)*s(2,1) + n(2,1)*s(1,2)*s(2,1) + n(3,1)*s(1,3)*s(2,1) + n(1,2)*s(1,1)*s(2,2) + n(2,2)*s(1,2)*s(2,2) + &
!           n(3,2)*s(1,3)*s(2,2) + n(1,3)*s(1,1)*s(2,3) + n(2,3)*s(1,2)*s(2,3) + n(3,3)*s(1,3)*s(2,3)
!    m(1,3)=n(1,1)*s(1,1)*s(3,1) + n(2,1)*s(1,2)*s(3,1) + n(3,1)*s(1,3)*s(3,1) + n(1,2)*s(1,1)*s(3,2) + n(2,2)*s(1,2)*s(3,2) + &
!           n(3,2)*s(1,3)*s(3,2) + n(1,3)*s(1,1)*s(3,3) + n(2,3)*s(1,2)*s(3,3) + n(3,3)*s(1,3)*s(3,3)
!    m(2,1)=n(1,1)*s(1,1)*s(2,1) + n(1,2)*s(1,2)*s(2,1) + n(1,3)*s(1,3)*s(2,1) + n(2,1)*s(1,1)*s(2,2) + n(2,2)*s(1,2)*s(2,2) + &
!           n(2,3)*s(1,3)*s(2,2) + n(3,1)*s(1,1)*s(2,3) + n(3,2)*s(1,2)*s(2,3) + n(3,3)*s(1,3)*s(2,3)
!    m(2,2)=n(1,1)*s(2,1)**2 + n(1,2)*s(2,1)*s(2,2) + n(2,1)*s(2,1)*s(2,2) + n(2,2)*s(2,2)**2 + n(1,3)*s(2,1)*s(2,3) + &
!           n(3,1)*s(2,1)*s(2,3) + n(2,3)*s(2,2)*s(2,3) + n(3,2)*s(2,2)*s(2,3) + n(3,3)*s(2,3)**2
!    m(2,3)=n(1,1)*s(2,1)*s(3,1) + n(2,1)*s(2,2)*s(3,1) + n(3,1)*s(2,3)*s(3,1) + n(1,2)*s(2,1)*s(3,2) + n(2,2)*s(2,2)*s(3,2) + &
!           n(3,2)*s(2,3)*s(3,2) + n(1,3)*s(2,1)*s(3,3) + n(2,3)*s(2,2)*s(3,3) + n(3,3)*s(2,3)*s(3,3)
!    m(3,1)=n(1,1)*s(1,1)*s(3,1) + n(1,2)*s(1,2)*s(3,1) + n(1,3)*s(1,3)*s(3,1) + n(2,1)*s(1,1)*s(3,2) + n(2,2)*s(1,2)*s(3,2) + &
!           n(2,3)*s(1,3)*s(3,2) + n(3,1)*s(1,1)*s(3,3) + n(3,2)*s(1,2)*s(3,3) + n(3,3)*s(1,3)*s(3,3)
!    m(3,2)=n(1,1)*s(2,1)*s(3,1) + n(1,2)*s(2,2)*s(3,1) + n(1,3)*s(2,3)*s(3,1) + n(2,1)*s(2,1)*s(3,2) + n(2,2)*s(2,2)*s(3,2) + &
!           n(2,3)*s(2,3)*s(3,2) + n(3,1)*s(2,1)*s(3,3) + n(3,2)*s(2,2)*s(3,3) + n(3,3)*s(2,3)*s(3,3)
!    m(3,3)=n(1,1)*s(3,1)**2 + n(1,2)*s(3,1)*s(3,2) + n(2,1)*s(3,1)*s(3,2) + n(2,2)*s(3,2)**2 + n(1,3)*s(3,1)*s(3,3) + &
!           n(3,1)*s(3,1)*s(3,3) + n(2,3)*s(3,2)*s(3,3) + n(3,2)*s(3,2)*s(3,3) + n(3,3)*s(3,3)**2
!end function

!!> triple matrix product m = s^T * n * s for 3x3 matrices
!pure function matmul_STNS(n,s) result(m)
!    !> matrix n
!    real(r8), dimension(3,3), intent(in) :: n
!    !> operation S
!    real(r8), dimension(3,3), intent(in) :: s
!    !> transformed matrix
!    real(r8), dimension(3,3) :: m
!
!    m(1,1)=n(1,1)*s(1,1)**2 + n(1,2)*s(1,1)*s(2,1) + n(2,1)*s(1,1)*s(2,1) + n(2,2)*s(2,1)**2 + n(1,3)*s(1,1)*s(3,1) + &
!           n(3,1)*s(1,1)*s(3,1) + n(2,3)*s(2,1)*s(3,1) + n(3,2)*s(2,1)*s(3,1) + n(3,3)*s(3,1)**2
!    m(1,2)=n(1,1)*s(1,1)*s(1,2) + n(2,1)*s(1,2)*s(2,1) + n(1,2)*s(1,1)*s(2,2) + n(2,2)*s(2,1)*s(2,2) + n(3,1)*s(1,2)*s(3,1) + &
!           n(3,2)*s(2,2)*s(3,1) + n(1,3)*s(1,1)*s(3,2) + n(2,3)*s(2,1)*s(3,2) + n(3,3)*s(3,1)*s(3,2)
!    m(1,3)=n(1,1)*s(1,1)*s(1,3) + n(2,1)*s(1,3)*s(2,1) + n(1,2)*s(1,1)*s(2,3) + n(2,2)*s(2,1)*s(2,3) + n(3,1)*s(1,3)*s(3,1) + &
!           n(3,2)*s(2,3)*s(3,1) + n(1,3)*s(1,1)*s(3,3) + n(2,3)*s(2,1)*s(3,3) + n(3,3)*s(3,1)*s(3,3)
!    m(2,1)=n(1,1)*s(1,1)*s(1,2) + n(1,2)*s(1,2)*s(2,1) + n(2,1)*s(1,1)*s(2,2) + n(2,2)*s(2,1)*s(2,2) + n(1,3)*s(1,2)*s(3,1) + &
!           n(2,3)*s(2,2)*s(3,1) + n(3,1)*s(1,1)*s(3,2) + n(3,2)*s(2,1)*s(3,2) + n(3,3)*s(3,1)*s(3,2)
!    m(2,2)=n(1,1)*s(1,2)**2 + n(1,2)*s(1,2)*s(2,2) + n(2,1)*s(1,2)*s(2,2) + n(2,2)*s(2,2)**2 + n(1,3)*s(1,2)*s(3,2) + &
!           n(3,1)*s(1,2)*s(3,2) + n(2,3)*s(2,2)*s(3,2) + n(3,2)*s(2,2)*s(3,2) + n(3,3)*s(3,2)**2
!    m(2,3)=n(1,1)*s(1,2)*s(1,3) + n(2,1)*s(1,3)*s(2,2) + n(1,2)*s(1,2)*s(2,3) + n(2,2)*s(2,2)*s(2,3) + n(3,1)*s(1,3)*s(3,2) +&
!           n(3,2)*s(2,3)*s(3,2) + n(1,3)*s(1,2)*s(3,3) + n(2,3)*s(2,2)*s(3,3) + n(3,3)*s(3,2)*s(3,3)
!    m(3,1)=n(1,1)*s(1,1)*s(1,3) + n(1,2)*s(1,3)*s(2,1) + n(2,1)*s(1,1)*s(2,3) + n(2,2)*s(2,1)*s(2,3) + n(1,3)*s(1,3)*s(3,1) +&
!           n(2,3)*s(2,3)*s(3,1) + n(3,1)*s(1,1)*s(3,3) + n(3,2)*s(2,1)*s(3,3) + n(3,3)*s(3,1)*s(3,3)
!    m(3,2)=n(1,1)*s(1,2)*s(1,3) + n(1,2)*s(1,3)*s(2,2) + n(2,1)*s(1,2)*s(2,3) + n(2,2)*s(2,2)*s(2,3) + n(1,3)*s(1,3)*s(3,2) +&
!           n(2,3)*s(2,3)*s(3,2) + n(3,1)*s(1,2)*s(3,3) + n(3,2)*s(2,2)*s(3,3) + n(3,3)*s(3,2)*s(3,3)
!    m(3,3)=n(1,1)*s(1,3)**2 + n(1,2)*s(1,3)*s(2,3) + n(2,1)*s(1,3)*s(2,3) + n(2,2)*s(2,3)**2 + n(1,3)*s(1,3)*s(3,3) + &
!           n(3,1)*s(1,3)*s(3,3) + n(2,3)*s(2,3)*s(3,3) + n(3,2)*s(2,3)*s(3,3) + n(3,3)*s(3,3)**2
!end function
!
!!> return a smaller cell
!subroutine return_primitive_cell(p)
!    !> original structure
!    type(lo_crystalstructure), intent(inout) :: p
!
!    real(r8), dimension(3,3) :: primbasis
!    real(r8), dimension(:,:), allocatable :: r
!    integer, dimension(:), allocatable :: species
!    integer, dimension(p%na) :: flavor
!    integer :: i,j,k,l
!    flavor=0
!
!    call find_primitive_cell(p%latticevectors,p%r,p%species,flavor,primbasis,r,species,1E-3_r8)
!
!write(*,*) 'primcell'
!do i=1,3
!    write(*,*) lo_chop(primbasis(:,i),lo_sqtol)
!enddo
!write(*,*) 'positions'
!do i=1,size(r,2)
!    write(*,*) lo_chop(r(:,i),lo_sqtol),species(i)
!enddo
!
!end subroutine

!> Find the proper unit cell, maybe
subroutine find_true_unitcell(p)
    !> proposed unit cell
    type(lo_crystalstructure), intent(in) :: p

    real(r8), dimension(3, 3) :: newbasis
    real(r8), dimension(:, :), allocatable :: newr
    integer, dimension(:), allocatable :: newspecies
    integer, dimension(:), allocatable :: oldflavor !,newflavor
    integer :: na, i

    na = p%na
    allocate (oldflavor(na))
    oldflavor = 0

    write (*, *) ''
    do i = 1, 3
        write (*, *) p%latticevectors(:, i)
    end do

    call find_primitive_cell(p%latticevectors, p%r, p%atomic_number, oldflavor, newbasis, newr, newspecies, 1E-5_r8)

    write (*, *) ''
    do i = 1, 3
        write (*, *) newbasis(:, i)
    end do
    do i = 1, size(newr, 2)
        write (*, *) i, newr(:, i)
    end do

end subroutine

!> Find the primitive cell
subroutine find_primitive_cell(basis, r, species, flavor, primbasis, primr, primspecies, tolerance)
    !> cell to be reduced
    real(r8), dimension(3, 3), intent(in) :: basis
    !> atoms in the cell
    real(r8), dimension(:, :), intent(in) :: r
    !> species classification
    integer, dimension(:), intent(in) :: species
    !> flavor classification
    integer, dimension(:), intent(in) :: flavor
    !> primitive lattice vectors
    real(r8), dimension(3, 3), intent(out) :: primbasis
    !> positions in the new coordinates
    real(r8), dimension(:, :), allocatable, intent(out) :: primr
    !> species in new coordinates
    integer, dimension(:), allocatable, intent(out) :: primspecies
    !> tolerance
    real(r8), intent(in) :: tolerance
    !
    real(r8), dimension(:, :), allocatable :: r0, r1
    real(r8), dimension(3, 3) :: basis0, basis1
    integer, dimension(:), allocatable :: species0, flavor0, species1, flavor1
    integer :: i, na
    !
    na = size(r, 2)
    allocate (r0(3, na))
    allocate (species0(na))
    allocate (flavor0(na))
    basis0 = basis
    r0 = r
    species0 = species
    flavor0 = flavor

    ! Recursively find a smaller cell
    redloop: do
        ! Try to find a smaller cell
        call reduce_cell(basis0, r0, species0, flavor0, basis1, r1, species1, flavor1, tolerance)
        if (abs(abs(lo_determ(basis0)) - abs(lo_determ(basis1))) .lt. tolerance) then
            ! this means no new cell, we are converged.
            exit redloop
        else
            ! this means a new cell, update and try again!
            if (allocated(r0)) deallocate (r0)
            if (allocated(species0)) deallocate (species0)
            if (allocated(flavor0)) deallocate (flavor0)
            i = size(r1, 2)
            allocate (r0(3, i))
            allocate (species0(i))
            allocate (flavor0(i))
            basis0 = basis1
            r0 = r1
            species0 = species1
            flavor0 = flavor1
            if (allocated(r1)) deallocate (r1)
            if (allocated(species1)) deallocate (species1)
            if (allocated(flavor1)) deallocate (flavor1)
        end if
    end do redloop

    ! Here we should be done.
    primbasis = basis0
    na = size(r0, 2)
    allocate (primr(3, na))
    allocate (primspecies(na))
    primr = r0
    primspecies = species0

    if (allocated(r0)) deallocate (r0)
    if (allocated(r1)) deallocate (r1)
    if (allocated(species0)) deallocate (species0)
    if (allocated(species1)) deallocate (species1)
    if (allocated(flavor0)) deallocate (flavor0)
    if (allocated(flavor1)) deallocate (flavor1)

contains

    !> single iteration of the cell reduction
    subroutine reduce_cell(latticevectors, r, species, flavor, newlatticevectors, newr, newspecies, newflavor, tolerance)
        !> original lattice vectors
        real(r8), dimension(3, 3), intent(in) :: latticevectors
        !> original positions
        real(r8), dimension(:, :), intent(in) :: r
        !> original species
        integer, dimension(:), intent(in) :: species
        !> additional flavor for symmetry thing
        integer, dimension(:), intent(in) :: flavor
        !> tolerance
        real(r8), intent(in) :: tolerance

        real(r8), dimension(3, 3), intent(out) :: newlatticevectors
        real(r8), dimension(:, :), allocatable, intent(out) :: newr
        integer, dimension(:), allocatable, intent(out) :: newspecies
        integer, dimension(:), allocatable, intent(out) :: newflavor

        real(r8), dimension(:, :), allocatable :: vecs, dum
        real(r8), dimension(3, 3) :: m0, m1
        real(r8) :: vol, f0, tol, reltol, frtol
        integer :: i, j, k, l, nvecs, na
        logical :: cellvalid

        ! set the tolerance
        call lo_fetch_tolerance(tolerance, latticevectors, realspace_cart_tol=tol, realspace_fractional_tol=frtol, relative_tol=reltol)

        na = size(r, 2)
        ! Get a list of possible lattice vectors: All the vectors in the cell + old latticevectors
        l = na + 3
        allocate (dum(3, na + 3))
        dum(:, 1:na) = lo_clean_fractional_coordinates(r, frtol**2)
        do i = 1, na
            dum(:, i) = matmul(latticevectors, dum(:, i))
        end do
        do i = 1, 3
            dum(:, i + na) = latticevectors(:, i)
        end do
        ! Reduce these to the unique. Should not be necessary, but never hurts.
        call lo_return_unique(dum, vecs, frtol)
        nvecs = size(vecs, 2)

        ! volume of original cell
        vol = abs(lo_determ(latticevectors))
        ! Loop over all vectors to try to find a new cell.
        vl1: do i = 1, nvecs
        vl2: do j = i + 1, nvecs
        vl3: do k = j + 1, nvecs
            ! possible new set of lattice vectors
            m0(:, 1) = vecs(:, i)
            m0(:, 2) = vecs(:, j)
            m0(:, 3) = vecs(:, k)
            ! some quick tests to see if there is any point in continuing:
            f0 = abs(lo_determ(m0))
            if (f0 .lt. tol) cycle vl3                           ! stupid cell if volume is zero
            if (abs(mod(vol/f0, 1.0_r8)) .gt. reltol) cycle vl3 ! has to be an integer division
            if (int(anint(vol/f0)) .eq. 1) cycle vl3             ! has to make the cell smaller, not the same
            ! I can get really retarded cells from this, do the Delaunay thing on it
            call reduce_basis_to_something_decent(m0, m1)
            ! Now do an actual test
            call test_cell(latticevectors, r, species, flavor, m1, cellvalid, newr, newspecies, newflavor, tolerance)
            if (cellvalid) then
                ! If I find a new cell, we will exit this routine and start over. Makes things way way faster.
                newlatticevectors = m1
                deallocate (vecs)
                return
            end if
        end do vl3
        end do vl2
        end do vl1
        ! If we made it here, the cell is already the smallest. Return the same basis again.
        newlatticevectors = latticevectors
    end subroutine

    !> test if a new cell is a smaller unit cell
    subroutine test_cell(oldcell, oldr, oldspecies, oldflavor, cell, valid, newr, newspecies, newflavor, tolerance)
        !> original cell
        real(r8), dimension(3, 3), intent(in) :: oldcell
        !> original positions
        real(r8), dimension(:, :), intent(in) :: oldr
        !> original species
        integer, dimension(:), intent(in) :: oldspecies
        !> original flavor
        integer, dimension(:), intent(in) :: oldflavor
        !> new cell
        real(r8), dimension(3, 3), intent(in) :: cell
        !> is it a new cell?
        logical, intent(out) :: valid
        !> new positions in case of a valid cell
        real(r8), dimension(:, :), allocatable, intent(out) :: newr
        !> new species in case of a valid cell
        integer, dimension(:), allocatable, intent(out) :: newspecies
        !> new flavor in case of a valid cell
        integer, dimension(:), allocatable, intent(out) :: newflavor
        !> tolerance
        real(r8), intent(in) :: tolerance
        !
        real(r8), dimension(3, 3) :: icell, m0
        real(r8), dimension(:, :), allocatable :: dumr1, dumr2
        real(r8) :: tol, sqtol
        integer :: i, j, l, mult, oldna, newna

        ! set the tolerance
        call lo_fetch_tolerance(tolerance, oldcell, realspace_fractional_tol=tol)
        sqtol = tol**2

        valid = .false.
        ! How many atoms should I get of each?
        mult = int(anint(abs(lo_determ(oldcell)/lo_determ(cell))))

        ! m0 is a matrix that converts to fractional coordinates in the new cell
        icell = lo_invert3x3matrix(cell)
        m0 = matmul(icell, oldcell)

        ! convert coordinates to fractional with respect to the new cell. Also
        ! I decorate the positions with the species and flavor, to keep track of that.
        oldna = size(oldr, 2)
        allocate (dumr1(5, oldna))
        do i = 1, oldna
            dumr1(1:3, i) = lo_clean_fractional_coordinates(matmul(m0, oldr(:, i)), sqtol)
            dumr1(4:5, i) = [oldspecies(i), oldflavor(i)]
        end do
        ! Reduce this to the union of the list, the unique
        call lo_return_unique(dumr1, dumr2, tol)

        ! Did it become the correct number of atoms?
        if (abs(size(dumr2, 2)*mult - oldna) .ne. 0) then
            valid = .false.
            return
        end if

        ! Slightly better check: Each of the atoms in the reduced cell must appear
        ! exactly mult times in the old cell.
        do i = 1, size(dumr2, 2)
            l = 0
            do j = 1, oldna
                if (sum(abs(dumr1(:, j) - dumr2(:, i))) .lt. tol) l = l + 1
            end do
            ! kill if not true
            if (l .ne. mult) then
                valid = .false.
                return
            end if
        end do

        ! ok, now I say this is valid, return the new cell
        newna = size(dumr2, 2)
        allocate (newr(3, newna))
        allocate (newspecies(newna))
        allocate (newflavor(newna))
        newr = dumr2(1:3, :)
        newspecies = int(anint(dumr2(4, :)))
        newflavor = int(anint(dumr2(5, :)))
        valid = .true.
    end subroutine
end subroutine

!> reduce a set of lattice vectors to canonical form
subroutine reduce_basis_to_something_decent(basis, canonical_basis)
    !> original lattice vectors
    real(r8), dimension(3, 3), intent(in) :: basis
    !> canonical lattice vectors
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

    ! Look through all variants to get shortest vectors
    dum = 0.0_r8
    dum(:, 1:4) = ext_basis
    dum(:, 5) = ext_basis(:, 1) + ext_basis(:, 2)
    dum(:, 6) = ext_basis(:, 2) + ext_basis(:, 3)
    dum(:, 7) = ext_basis(:, 1) + ext_basis(:, 3)
    do i = 1, 7
        d(i) = norm2(dum(:, i))
    end do
    call qsort(d, ind)
    !
    do i = 3, 7
        canonical_basis(:, 1) = dum(:, ind(1))
        canonical_basis(:, 2) = dum(:, ind(2))
        canonical_basis(:, 3) = dum(:, ind(i))
        if (abs(lo_determ(canonical_basis)) .gt. lo_tol) exit
    end do
    if (lo_determ(canonical_basis) .lt. 0.0_r8) canonical_basis = -canonical_basis
end subroutine

!> get a cluster of atoms
subroutine cluster_from_lattice_vectors(basis, cluster, lvind)
    !> lattice vectors
    real(r8), dimension(3, 3), intent(in) :: basis
    !> cluster
    real(r8), dimension(:, :), allocatable, intent(out) :: cluster
    !> indices to lattice vectors?
    integer, dimension(3), intent(out) :: lvind

    type(lo_distancetable) :: dt
    real(r8), dimension(3, 1) :: r0
    real(r8) :: cutoff
    integer :: i, j, k

    ! Build a cluster?
    cutoff = lo_bounding_sphere_of_box(basis)
    cutoff = max(cutoff, norm2(basis(:, 1)))
    cutoff = max(cutoff, norm2(basis(:, 2)))
    cutoff = max(cutoff, norm2(basis(:, 3)))
    cutoff = cutoff*1.5_r8
    r0 = 0.0_r8
    call dt%generate(r0, basis, cutoff, -1)

    allocate (cluster(3, dt%particle(1)%n))
    cluster = dt%particle(1)%v

    do j = 1, 3
        k = 0
        do i = 1, size(cluster, 2)
            if (norm2(cluster(:, i) - basis(:, j)) .lt. 1E-10_r8) then
                k = i
                exit
            end if
        end do
        if (k .gt. 0) then
            lvind(j) = k
        else
            call lo_stop_gracefully(['Could not locate lattice vector'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
    end do
end subroutine

function cluster_transform_error(op, r) result(err)
    !> operations
    real(r8), dimension(:, :, :), intent(in) :: op
    !> points to transform
    real(r8), dimension(:, :), intent(in) :: r
    !> error
    real(r8) :: err

    real(r8), dimension(3, size(r, 2)) :: u
    real(r8) :: f0, f1
    integer :: i, j, n, iop

    err = 0.0_r8
    n = size(r, 2)
    do iop = 1, size(op, 3)
        ! Rotate cluster
        u = matmul(op(:, :, iop), r)
        ! Measure error
        do i = 1, n
            f1 = lo_huge
            do j = 1, n
                f0 = norm2(r(:, i) - u(:, j))
                f1 = min(f1, f0)
            end do
            err = err + f1
        end do
    end do
    err = err/real(size(r, 2)*size(op, 3), r8)
end function

function operation_cell_transform_error(op, r) result(err)
    !> operation
    real(r8), dimension(3, 3), intent(in) :: op
    !> (fractional) points to transform
    real(r8), dimension(:, :), intent(in) :: r
    !> error
    real(r8) :: err

    real(r8), dimension(3, size(r, 2)) :: u
    real(r8) :: f0, f1
    integer :: i, j, n, iop

    err = 0.0_r8
    n = size(r, 2)
    ! Rotate cluster
    u = matmul(op, r)
    ! Measure error
    do i = 1, n
        f1 = lo_huge
        do j = 1, n
            f0 = norm2(r(:, i) - u(:, j))
            f1 = min(f1, f0)
        end do
        err = err + f1
    end do
    err = err/real(size(r, 2), r8)
end function

function operation_transform_error(op, r) result(err)
    !> operations
    real(r8), dimension(3, 3), intent(in) :: op
    !> points to transform
    real(r8), dimension(:, :), intent(in) :: r
    !> error
    real(r8) :: err

    real(r8), dimension(3, size(r, 2)) :: u
    real(r8) :: f0, f1
    integer :: i, j, n, iop

    err = 0.0_r8
    n = size(r, 2)
    ! Rotate cluster
    u = matmul(op, r)
    ! Measure error
    do i = 1, n
        f1 = lo_huge
        do j = 1, n
            f0 = norm2(r(:, i) - u(:, j))
            f1 = min(f1, f0)
        end do
        err = err + f1
    end do
    err = err/real(size(r, 2), r8)
end function

!> make sure a subset of operations actually form a group
subroutine expand_valid_operations(op, n, tolerance)
    !> operations
    real(r8), dimension(3, 3, 48), intent(inout) :: op
    !> number of valid operations
    integer, intent(inout) :: n
    !> tolerance
    real(r8), intent(in) :: tolerance

    real(r8), dimension(3, 3) :: m0
    integer :: i, j, k, l, iter

    iterl: do iter = 1, 100
        l = n
        do i = 1, l
        vl2: do j = 1, l
            m0 = matmul(op(:, :, i), op(:, :, j))
            do k = 1, l
                if (norm2(m0 - op(:, :, k)) .lt. tolerance) then
                    cycle vl2
                end if
            end do
            n = n + 1
            op(:, :, n) = m0
            ! it was new!
            cycle iterl
        end do vl2
        end do
        if (n .eq. l) then
            ! If we make it here it's done
            exit iterl
        end if
    end do iterl
end subroutine

!> sort atoms into per-species boxes
subroutine sort_atoms_into_boxes(box, r, species, mem)
    !> box sorted per species
    class(lo_speciesbox), intent(out) :: box
    !> coordinates of atoms (fractional)
    real(r8), dimension(:, :), intent(in) :: r
    !> species of atoms
    integer, dimension(:), intent(in) :: species
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    integer, dimension(:), allocatable :: di1, di2
    integer :: i, j
    ! Get the unique species
    call mem%allocate(di2, size(species), persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call lo_return_unique_indices(species, di1, mem, di2)

    box%n_species = size(di1)
    allocate (box%species(box%n_species))
    do i = 1, box%n_species
        box%species(i)%n = 0
    end do
    do i = 1, size(species)
        j = di2(i)
        box%species(j)%n = box%species(j)%n + 1
    end do
    do i = 1, box%n_species
        allocate (box%species(i)%r(3, box%species(i)%n))
        allocate (box%species(i)%ind(box%species(i)%n))
        box%species(i)%r = 0.0_r8
        box%species(i)%ind = 0
        box%species(i)%n = 0
    end do
    do i = 1, size(species)
        j = di2(i)
        box%species(j)%n = box%species(j)%n + 1
        box%species(j)%r(:, box%species(j)%n) = r(:, i)
        box%species(j)%ind(box%species(j)%n) = i
    end do
    call mem%deallocate(di1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(di2, persistent=.true., scalable=.false., file=__FILE__, line=__LINE__)

    ! Temporary space for transformed positions?
    j = 0
    do i = 1, box%n_species
        j = max(j, box%species(i)%n)
    end do
    allocate (box%tf(3, j))
    box%tf = 0.0_r8

    box%n_atom = size(species)
end subroutine

!> measure the mismatch after a symmetry operation
function box_operation_error(box, op, tr) result(err)
    !> box sorted per species
    class(lo_speciesbox), intent(inout) :: box
    !> rotation
    real(r8), dimension(3, 3), intent(in) :: op
    !> translation
    real(r8), dimension(3), intent(in) :: tr
    !> error
    real(r8) :: err

    real(r8), parameter :: tighttolerance = 1E-13_r8
    real(r8), dimension(3) :: v0
    real(r8) :: f0, f1
    integer :: isp, i, j

    err = 0.0_r8
    do isp = 1, box%n_species
        ! Transform positions
        do i = 1, box%species(isp)%n
            v0 = matmul(op, box%species(isp)%r(:, i)) + tr
            v0 = lo_chop(v0, tighttolerance)
            v0 = lo_clean_fractional_coordinates(v0, tighttolerance)
            v0 = lo_chop(v0, tighttolerance)
            v0 = lo_clean_fractional_coordinates(v0, tighttolerance)
            box%tf(:, i) = v0
        end do
        ! Measure error?
        do i = 1, box%species(isp)%n
            f1 = lo_huge
            do j = 1, box%species(isp)%n
                v0 = box%tf(:, j) - box%species(isp)%r(:, i)
                v0 = lo_clean_fractional_coordinates(v0 + 0.5_r8) - 0.5_r8
                f0 = norm2(v0)
                f1 = min(f0, f1)
            end do
            err = err + f1
        end do
    end do
    err = err/real(box%n_atom, r8)
end function

!> refine positions according to a set of symmetry operations
subroutine box_refine_positions(box, op, tr)
    !> box sorted per species
    class(lo_speciesbox), intent(inout) :: box
    !> rotation
    real(r8), dimension(:, :, :), intent(in) :: op
    !> translation
    real(r8), dimension(:, :), intent(in) :: tr

    real(r8), dimension(:, :), allocatable :: delta
    integer, dimension(:), allocatable :: ctr
    real(r8), parameter :: tighttolerance = 1E-13_r8
    real(r8), dimension(3) :: v0, v1
    real(r8) :: f0, f1
    integer :: iop, isp, i, j

    allocate (delta(3, box%n_atom))
    allocate (ctr(box%n_atom))
    delta = 0.0_r8
    ctr = 0.0

    do iop = 1, size(op, 3)
    do isp = 1, box%n_species
        ! Transform positions
        do i = 1, box%species(isp)%n
            v0 = matmul(op(:, :, iop), box%species(isp)%r(:, i)) + tr(:, iop)
            v0 = lo_chop(v0, tighttolerance)
            v0 = lo_clean_fractional_coordinates(v0, tighttolerance)
            v0 = lo_chop(v0, tighttolerance)
            v0 = lo_clean_fractional_coordinates(v0, tighttolerance)
            box%tf(:, i) = v0
        end do
        ! Measure error?
        do i = 1, box%species(isp)%n
            f1 = lo_huge
            do j = 1, box%species(isp)%n
                v0 = box%tf(:, j) - box%species(isp)%r(:, i)
                v0 = lo_clean_fractional_coordinates(v0 + 0.5_r8) - 0.5_r8
                f0 = norm2(v0)
                if (f0 .lt. f1) then
                    v1 = v0
                    f1 = f0
                end if
            end do
            j = box%species(isp)%ind(i)
            ! Accumulate average
            delta(:, j) = delta(:, j) + v1
            ctr(j) = ctr(j) + 1
        end do
    end do
    end do

    ! Sort out the average
    do i = 1, box%n_atom
        delta(:, i) = delta(:, i)/real(ctr(i), r8)
    end do
    do isp = 1, box%n_species
    do i = 1, box%species(isp)%n
        j = box%species(isp)%ind(i)
        v0 = box%species(isp)%r(:, i) + delta(:, j)
        v0 = lo_chop(v0, tighttolerance)
        v0 = lo_clean_fractional_coordinates(v0, tighttolerance)
        v0 = lo_chop(v0, tighttolerance)
        v0 = lo_clean_fractional_coordinates(v0, tighttolerance)
        box%species(isp)%r(:, i) = v0
    end do
    end do

    deallocate (delta)
    deallocate (ctr)
end subroutine

end module
