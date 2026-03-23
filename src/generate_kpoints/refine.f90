#include "precompilerdefinitions"
module refine
!! Refine wedge-based meshes
use konstanter
use gottochblandat
use geometryfunctions
use type_crystalstructure
use type_qpointmesh
use type_symmetryoperation
use lo_sorting, only: lo_qsort
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
implicit none

private
public :: massage_mesh
public :: get_histograms
public :: fake_integration

type lo_fireminimizer
    real(r8) :: P
    real(r8) :: ts
    real(r8) :: timestepmax
    real(r8) :: timestepmin
    real(r8) :: alpha0
    real(r8) :: alpha
    real(r8) :: ftsinc
    real(r8) :: ftsdec
    real(r8) :: falpha
    real(r8) :: nndist
    real(r8) :: tolerance
    integer :: nmin
    integer :: counter
    integer :: niter
    integer :: exitstatus
end type

contains

!> try to make wedge-based meshes slightly more uniform
subroutine massage_mesh(qp, p, idealside, idealvolume, nstep, forcetolerance, sizetolerance, splitedge, mw, mem, verbosity)
    !> wedge mesh
    type(lo_wedge_mesh), intent(inout) :: qp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> ideal side
    real(r8), intent(in) :: idealside
    !> ideal volume
    real(r8), intent(in) :: idealvolume
    !> how many steps to optimize for, maximum
    integer, intent(in) :: nstep
    !> when to stop the relaxation
    real(r8), intent(in) :: forcetolerance
    !> size of tetrahedrons to split
    real(r8), intent(in) :: sizetolerance
    !> split edges
    logical, intent(in) :: splitedge
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> Talk a lot?
    integer, intent(in) :: verbosity

    ! Place to stuff the points
    type(lo_fireminimizer) :: fire
    real(r8), dimension(:, :), allocatable :: r0, r, f, v, linevecs
    real(r8), dimension(:), allocatable :: maxdist_per_point
    integer, dimension(:), allocatable :: plctr, planeind
    real(r8) :: timer, radius, avgvol, mindist, maxdist, maxdr, orig_mindist
    real(r8) :: tolerance

    timer = walltime()
    if (verbosity .gt. 0) then
        write (*, *) ''
        write (*, *) 'OPTIMIZING MESH'
    end if

    ! What is a sensible size tolerance?
    tolerance = idealside*1E-3_r8

    ! First make sure we have a sensible amount of points along each edge?
    if (splitedge) then
        fixedges: block
            real(r8), dimension(:, :), allocatable :: ctrpts, dr0, dr1
            real(r8), dimension(3) :: v0, v1, v2
            real(r8) :: f0, f1
            integer, dimension(:, :), allocatable :: di, dj
            integer :: nctr, nedge
            integer :: ikp, i, j, k, l, ii, jj

            ! First, stow away the points that are not along lines
            allocate (ctrpts(3, qp%n_irr_point))
            ctrpts = lo_huge
            nctr = 0
            do ikp = 1, qp%n_irr_point
                v0 = qp%ip(ikp)%r
                j = 0
                do i = 1, p%irrw%nfaces
                    f0 = p%irrw%face(i)%plane%distance_to_point(v0)
                    if (abs(f0) .lt. tolerance) j = j + 1
                end do
                if (j .lt. 2) then
                    nctr = nctr + 1
                    ctrpts(:, nctr) = v0
                end if
            end do
            if (mw%talk) write (*, *) '... got ', tochar(nctr), ' points not on edges (out of ', tochar(qp%n_irr_point), ')'

            ! Seems sensible. Now get a list of the possible edges?
            j = 0
            do i = 1, p%irrw%nfaces
                j = j + p%irrw%face(i)%n
            end do
            allocate (di(2, j))
            di = 0
            l = 0
            do i = 1, p%irrw%nfaces
                do j = 1, p%irrw%face(i)%n
                    k = lo_index_in_periodic_array(j + 1, p%irrw%face(i)%n)
                    ii = p%irrw%face(i)%ind(j)
                    jj = p%irrw%face(i)%ind(k)
                    l = l + 1
                    if (ii .gt. jj) then
                        di(:, l) = [jj, ii]
                    else
                        di(:, l) = [ii, jj]
                    end if
                end do
            end do
            ! Get the unique edges
            call lo_return_unique(di, dj)
            nedge = size(dj, 2)
            if (mw%talk) write (*, *) '... found ', tochar(nedge), ' edges to split'

            ! First count the number of points?
            l = 0
            do i = 1, nedge
                v0 = p%irrw%r(:, dj(1, i))
                v1 = p%irrw%r(:, dj(2, i))
                v2 = v1 - v0
                f0 = norm2(v2)
                j = 1
                do
                    if (f0/j .lt. idealside*0.75_r8) then
                        exit
                    else
                        j = j + 1
                    end if
                end do
                j = max(j, 3)
                l = l + j
            end do
            allocate (dr0(3, l))
            dr0 = 0.0_r8
            ! Store the points
            l = 0
            do i = 1, nedge
                v0 = p%irrw%r(:, dj(1, i))
                v1 = p%irrw%r(:, dj(2, i))
                v2 = v1 - v0
                f0 = norm2(v2)
                j = 1
                do
                    if (f0/j .lt. idealside*0.75_r8) then
                        exit
                    else
                        j = j + 1
                    end if
                end do
                j = max(j, 3)
                do k = 1, j
                    f1 = real(k - 1, r8)/real(j - 1, r8)
                    l = l + 1
                    dr0(:, l) = v0 + f1*v2
                end do
            end do
            ! And get the unique points
            call lo_return_unique(dr0, dr1, tolerance)
            ! stuff the points together and re-tesselate the mesh
            deallocate (dr0)
            allocate (dr0(3, nctr + size(dr1, 2)))
            dr0 = 0.0_r8
            dr0(:, 1:nctr) = ctrpts(:, 1:nctr)
            dr0(:, nctr + 1:nctr + size(dr1, 2)) = dr1
            call qp%tesselate_wedge_mesh(p, dr0, idealside*0.25_r8, -1.0_r8, mw, mem, verbosity)
        end block fixedges
    end if

    ! If there are some really large tetrahedrons, it makes sense to split them!
    if (sizetolerance .gt. 0.0_r8) then
        fixlarge: block
            integer, parameter :: ntetpts = 31
            real(r8), dimension(4, ntetpts) :: bc
            real(r8), dimension(:, :), allocatable :: dr0, dr1
            real(r8), dimension(3, 4) :: dumt, tmtx
            real(r8), dimension(3) :: v0, v1, v2
            real(r8) :: f0, f1
            integer, dimension(:), allocatable :: badtets
            integer :: itet, icrn, ikp, ctr, i, j, k, l

            ! Sensible test-points to inject
            bc(:, 1) = [0.000000_r8, 0.000000_r8, 0.250000_r8, 0.750000_r8]
            bc(:, 2) = [0.000000_r8, 0.000000_r8, 0.500000_r8, 0.500000_r8]
            bc(:, 3) = [0.000000_r8, 0.000000_r8, 0.750000_r8, 0.250000_r8]
            bc(:, 4) = [0.000000_r8, 0.250000_r8, 0.000000_r8, 0.750000_r8]
            bc(:, 5) = [0.000000_r8, 0.250000_r8, 0.250000_r8, 0.500000_r8]
            bc(:, 6) = [0.000000_r8, 0.250000_r8, 0.500000_r8, 0.250000_r8]
            bc(:, 7) = [0.000000_r8, 0.250000_r8, 0.750000_r8, 0.000000_r8]
            bc(:, 8) = [0.000000_r8, 0.500000_r8, 0.000000_r8, 0.500000_r8]
            bc(:, 9) = [0.000000_r8, 0.500000_r8, 0.250000_r8, 0.250000_r8]
            bc(:, 10) = [0.000000_r8, 0.500000_r8, 0.500000_r8, 0.000000_r8]
            bc(:, 11) = [0.000000_r8, 0.750000_r8, 0.000000_r8, 0.250000_r8]
            bc(:, 12) = [0.000000_r8, 0.750000_r8, 0.250000_r8, 0.000000_r8]
            bc(:, 13) = [0.250000_r8, 0.000000_r8, 0.000000_r8, 0.750000_r8]
            bc(:, 14) = [0.250000_r8, 0.000000_r8, 0.250000_r8, 0.500000_r8]
            bc(:, 15) = [0.250000_r8, 0.000000_r8, 0.500000_r8, 0.250000_r8]
            bc(:, 16) = [0.250000_r8, 0.000000_r8, 0.750000_r8, 0.000000_r8]
            bc(:, 17) = [0.250000_r8, 0.250000_r8, 0.000000_r8, 0.500000_r8]
            bc(:, 18) = [0.250000_r8, 0.250000_r8, 0.250000_r8, 0.250000_r8]
            bc(:, 19) = [0.250000_r8, 0.250000_r8, 0.500000_r8, 0.000000_r8]
            bc(:, 20) = [0.250000_r8, 0.500000_r8, 0.000000_r8, 0.250000_r8]
            bc(:, 21) = [0.250000_r8, 0.500000_r8, 0.250000_r8, 0.000000_r8]
            bc(:, 22) = [0.250000_r8, 0.750000_r8, 0.000000_r8, 0.000000_r8]
            bc(:, 23) = [0.500000_r8, 0.000000_r8, 0.000000_r8, 0.500000_r8]
            bc(:, 24) = [0.500000_r8, 0.000000_r8, 0.250000_r8, 0.250000_r8]
            bc(:, 25) = [0.500000_r8, 0.000000_r8, 0.500000_r8, 0.000000_r8]
            bc(:, 26) = [0.500000_r8, 0.250000_r8, 0.000000_r8, 0.250000_r8]
            bc(:, 27) = [0.500000_r8, 0.250000_r8, 0.250000_r8, 0.000000_r8]
            bc(:, 28) = [0.500000_r8, 0.500000_r8, 0.000000_r8, 0.000000_r8]
            bc(:, 29) = [0.750000_r8, 0.000000_r8, 0.000000_r8, 0.250000_r8]
            bc(:, 30) = [0.750000_r8, 0.000000_r8, 0.250000_r8, 0.000000_r8]
            bc(:, 31) = [0.750000_r8, 0.250000_r8, 0.000000_r8, 0.000000_r8]

            ! Count tetrahedrons that are too large
            ctr = 0
            do itet = 1, qp%n_irr_tet
                ! fetch tetrahedron
                do icrn = 1, 4
                    dumt(:, icrn) = qp%ip(qp%it(itet)%irreducible_index(icrn))%r
                end do
                f0 = lo_unsigned_tetrahedron_volume(dumt)/idealvolume
                if (f0 .gt. sizetolerance) then
                    ctr = ctr + 1
                end if
            end do

            if (ctr .gt. 0) then
                if (mw%talk) write (*, *) '... found ', tochar(ctr), ' points to add'
                ! keep track of which the bad tetrahedrons are
                allocate (badtets(ctr))
                badtets = 0
                ctr = 0
                do itet = 1, qp%n_irr_tet
                    ! fetch tetrahedron
                    do icrn = 1, 4
                        dumt(:, icrn) = qp%ip(qp%it(itet)%irreducible_index(icrn))%r
                    end do
                    f0 = lo_unsigned_tetrahedron_volume(dumt)/idealvolume
                    if (f0 .gt. sizetolerance) then
                        ctr = ctr + 1
                        badtets(ctr) = itet
                    end if
                end do

                ! buffer to store points in
                allocate (dr0(3, qp%n_irr_point + ctr))
                dr0 = 0.0_r8
                do ikp = 1, qp%n_irr_point
                    dr0(:, ikp) = qp%ip(ikp)%r
                end do
                ! Start adding points, as far away as possible from any existing points
                l = qp%n_irr_point
                do i = 1, ctr
                    itet = badtets(i)
                    ! fetch tetrahedron
                    do icrn = 1, 4
                        dumt(:, icrn) = qp%ip(qp%it(itet)%irreducible_index(icrn))%r
                    end do
                    ! Get the conversion from barycentric coordinates to Cartesian
                    tmtx(1, 1) = dumt(1, 1) - dumt(1, 4)
                    tmtx(1, 2) = dumt(1, 2) - dumt(1, 4)
                    tmtx(1, 3) = dumt(1, 3) - dumt(1, 4)
                    tmtx(2, 1) = dumt(2, 1) - dumt(2, 4)
                    tmtx(2, 2) = dumt(2, 2) - dumt(2, 4)
                    tmtx(2, 3) = dumt(2, 3) - dumt(2, 4)
                    tmtx(3, 1) = dumt(3, 1) - dumt(3, 4)
                    tmtx(3, 2) = dumt(3, 2) - dumt(3, 4)
                    tmtx(3, 3) = dumt(3, 3) - dumt(3, 4)
                    tmtx(:, 4) = 0.0_r8
                    ! Reference corner in the tetrahedron
                    v0 = dumt(:, 4)
                    f0 = 0.0_r8
                    v2 = 0.0_r8
                    do j = 1, ntetpts
                        v1 = v0 + matmul(tmtx, bc(:, j))
                        f1 = lo_huge
                        do k = 1, l
                            f1 = min(f1, lo_sqnorm(dr0(:, k) - v1))
                        end do
                        if (f1 .gt. f0) then
                            f0 = f1
                            v2 = v1
                        end if
                    end do
                    ! add this point
                    l = l + 1
                    dr0(:, l) = v2
                end do
                ! retesselate the mesh
                call qp%tesselate_wedge_mesh(p, dr0, -1.0_r8, -1.0_r8, mw, mem, verbosity)
            end if
        end block fixlarge
    end if

    ! set up things for mesh optimization
    init: block
        type(lo_linesegment) :: linesegment
        type(lo_verletbox) :: vb
        real(r8), dimension(3, 2) :: plts
        real(r8) :: f0, f1, ptpldist
        integer :: i, j, k, l, bi, bj, bk, ii, jj, kk

        timer = walltime()

        allocate (r(3, qp%n_irr_point))
        allocate (r0(3, qp%n_irr_point))
        allocate (f(3, qp%n_irr_point))
        allocate (v(3, qp%n_irr_point))
        allocate (linevecs(3, qp%n_irr_point))
        allocate (plctr(qp%n_irr_point))
        allocate (planeind(qp%n_irr_point))
        allocate (maxdist_per_point(qp%n_irr_point))
        r0 = 0.0_r8
        r = 0.0_r8
        f = 0.0_r8
        v = 0.0_r8
        linevecs = 0.0_r8
        maxdist_per_point = 0.0_r8
        plctr = 0
        planeind = 0
        ptpldist = lo_huge
        ! Constrain some points?
        do i = 1, qp%n_irr_point
            if (mod(i, mw%n) .ne. mw%r) cycle
            r(:, i) = qp%ip(i)%r
            do j = 1, p%irrw%nfaces
                f0 = p%irrw%face(j)%plane%distance_to_point(r(:, i))
                if (abs(f0) .lt. lo_sqtol) then
                    plctr(i) = plctr(i) + 1
                else
                    ptpldist = min(abs(f0), ptpldist)
                end if
            end do
            ! Three or more planes is a point.
            plctr(i) = min(plctr(i), 3)

            select case (plctr(i))
            case (1)
                ! on a plane
                do j = 1, p%irrw%nfaces
                    f0 = p%irrw%face(j)%plane%distance_to_point(r(:, i))
                    if (abs(f0) .lt. lo_sqtol) planeind(i) = j
                end do
            case (2)
                ! on a line
                l = 0
                do j = 1, p%irrw%nfaces
                    f0 = p%irrw%face(j)%plane%distance_to_point(r(:, i))
                    if (abs(f0) .lt. lo_sqtol) then
                        l = l + 1
                        plts(:, l) = p%irrw%face(j)%plane%normal
                    end if
                end do
                linevecs(:, i) = lo_cross(plts(:, 1), plts(:, 2))
                linevecs(:, i) = linevecs(:, i)/norm2(linevecs(:, i))
            end select
        end do

        ! Sync things
        call mw%allreduce('sum', ptpldist)
        call mw%allreduce('sum', r)
        call mw%allreduce('sum', plctr)
        call mw%allreduce('sum', planeind)
        call mw%allreduce('sum', linevecs)
        r0 = r

        ! Get the nearest neighbour distance per point?
        ! I will measure the distance to neighbouring points, neighbouring planes
        ! and neighbouring lines and points. Then I will make sure that we are not
        ! allowed to move more than half that distance, seems like a safe choice?

        maxdist_per_point = 0.0_r8
        call vb%generate(r0, [15, 15, 15])
        do i = 1, qp%n_irr_point
            if (mod(i, mw%n) .ne. mw%r) cycle

            f0 = lo_huge
            call vb%boxind(r(:, i), bi, bj, bk)
            do ii = -1, 1
            do jj = -1, 1
            do kk = -1, 1
                do l = 1, vb%box(bi, bj, bk)%n
                    j = vb%box(bi, bj, bk)%ind(l)
                    if (j .eq. i) cycle
                    f1 = lo_sqnorm(r(:, i) - r(:, j))
                    f0 = min(f0, f1)
                end do
            end do
            end do
            end do
            ! Could not find neighbour with boxes, do it the old-fashined way
            if (f0 .gt. 1E10_r8) then
                do j = 1, qp%n_irr_point
                    if (j .eq. i) cycle
                    f1 = lo_sqnorm(r(:, i) - r(:, j))
                    f0 = min(f0, f1)
                end do
            end if
            ! Make sure it's an actual distance
            f0 = sqrt(f0)
            ! Now check with planes and lines and things
            select case (plctr(i))
            case (0)
                ! In the middle, check distance to planes
                do j = 1, p%irrw%nfaces
                    f1 = abs(p%irrw%face(j)%plane%distance_to_point(r(:, i)))
                    f0 = min(f1, f0)
                end do
            case (1)
                ! On a face, check distance to edges!
                do j = 1, p%irrw%nfaces
                do k = 1, p%irrw%face(j)%n
                    ii = p%irrw%face(j)%ind(k)
                    jj = lo_index_in_periodic_array(ii + 1, p%irrw%face(j)%n)
                    call linesegment%generate(p%irrw%r(:, ii), p%irrw%r(:, jj))
                    f1 = linesegment%distance_to_point(r(:, i))
                    f0 = min(f1, f0)
                end do
                end do
            end select
            ! Now we have the closest distance from each point to something it's not allowed to touch!
            maxdist_per_point(i) = f0*0.15_r8
        end do
        call mw%allreduce('sum', maxdist_per_point)

        ! Smallest pair distance
        orig_mindist = minval(maxdist_per_point)

        ! What is the largest distance
        avgvol = lo_determ(p%reciprocal_latticevectors)/p%sym%n/qp%n_irr_tet
        radius = (avgvol*3.0_r8/4.0_r8/lo_pi)**(1.0_r8/3.0_r8)
        avgvol = avgvol*1.1_r8

        ! Get the initial forces
        call forces_on_vertices(qp, p, f, r, r0, plctr, planeind, linevecs, idealside, maxdist_per_point, mindist, maxdist, maxdr, mw)

        ! Starting timestep for the FIRE relaxation to make it not go crazy when starting
        f0 = 0.0_r8
        do i = 1, qp%n_irr_point
            f0 = max(f0, norm2(f(:, i)))
        end do
        fire%ts = sqrt(radius*1E-1_r8/f0)

        ! Set parameters for the FIRE relaxation
        fire%exitstatus = 0
        fire%niter = nstep
        fire%P = 0.0_r8
        fire%timestepmax = 1E5*fire%ts             ! max timestep
        fire%timestepmin = 1E-20_r8*fire%ts      ! fire_ts*1E-4_r8
        fire%alpha0 = 0.5_r8                     ! initial mixing parameter
        fire%alpha = fire%alpha0                   ! current mixing parameter
        fire%nmin = 3                              ! minimum number of steps to take
        fire%ftsinc = 1.02_r8                    ! factor with which to increase the timestep
        fire%ftsdec = 0.75_r8                    ! factor with which to derease the timestep
        fire%falpha = 0.90_r8                    ! factor with which to derease the mixing
        fire%counter = 0                           ! just a counter to know how long ago it was we froze the system
        fire%nndist = radius                       ! measure of nearest neighbour distance
        fire%tolerance = forcetolerance            ! criteria for convergence
    end block init

    ! Do the actual relaxation
    optmesh: block
        real(r8), parameter :: maxdeltar = 0.1_r8
        real(r8), dimension(3) :: v0, v1, dv, dr
        real(r8) :: fx, fy, f0, f1
        integer :: i, j, k, iter, rsctr

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'FIRE MINIMIZER'
            write (*, *) '    rad:', fire%nndist
            write (*, *) '     ts:', fire%ts
            write (*, *) '  max f:', maxval(f)
            write (*, *) '  min d:', orig_mindist/idealside
        end if

        rsctr = 0
        relaxloop: do iter = 1, fire%niter

            ! First a normal Verlet step
            do i = 1, qp%n_irr_point
                ! v(:,i)=v(:,i)+f(:,i)*0.5_r8*fire%ts
                ! ! Make sure the velocities also are constrained properly
                ! select case(plctr(i))
                ! case(1) ! On a plane
                !     j=planeind(i)
                !     fx=dot_product(p%irrw%face(j)%plane%v1,v(:,i))
                !     fy=dot_product(p%irrw%face(j)%plane%v2,v(:,i))
                !     v(:,i)=fx*p%irrw%face(j)%plane%v1+fy*p%irrw%face(j)%plane%v2
                ! case(2) ! On a line
                !     v(:,i)=dot_product(linevecs(:,i),v(:,i))*linevecs(:,i)
                !     !v(:,i)=0.0_r8
                ! case(3) ! An endpoint
                !     v(:,i)=0.0_r8
                ! end select
                ! r(:,i)=r(:,i)+fire%ts*v(:,i)

                ! Make sure we don't do large steps.
                dv = f(:, i)*0.5_r8*fire%ts
                select case (plctr(i))
                case (1) ! On a plane
                    j = planeind(i)
                    fx = dot_product(p%irrw%face(j)%plane%v1, dv)
                    fy = dot_product(p%irrw%face(j)%plane%v2, dv)
                    dv = fx*p%irrw%face(j)%plane%v1 + fy*p%irrw%face(j)%plane%v2
                case (2) ! On a line
                    dv = dot_product(linevecs(:, i), dv)*linevecs(:, i)
                case (3) ! An endpoint
                    dv = 0.0_r8
                end select
                dr = fire%ts*dv
                if (norm2(dr) .gt. maxdist_per_point(i)*maxdeltar) then
                    dv = dv*(maxdist_per_point(i)*maxdeltar/fire%ts)/norm2(dv)
                end if
                v(:, i) = v(:, i) + dv

                ! Make sure the velocities also are constrained properly
                select case (plctr(i))
                case (1) ! On a plane
                    j = planeind(i)
                    fx = dot_product(p%irrw%face(j)%plane%v1, v(:, i))
                    fy = dot_product(p%irrw%face(j)%plane%v2, v(:, i))
                    v(:, i) = fx*p%irrw%face(j)%plane%v1 + fy*p%irrw%face(j)%plane%v2
                case (2) ! On a line
                    v(:, i) = dot_product(linevecs(:, i), v(:, i))*linevecs(:, i)
                    !v(:,i)=0.0_r8
                case (3) ! An endpoint
                    v(:, i) = 0.0_r8
                end select
                r(:, i) = r(:, i) + fire%ts*v(:, i)
            end do

            ! Get some kind of forces!
            call forces_on_vertices(qp, p, f, r, r0, plctr, planeind, linevecs, idealside, maxdist_per_point, mindist, maxdist, maxdr, mw)

            ! This is the second Verlet step. Also here, I do not want to add too much.
            ! hmm. I think the sensible thing to do is to scale down the forces instead?
            do i = 1, qp%n_irr_point
                v(:, i) = v(:, i) + f(:, i)*0.5_r8*fire%ts
            end do

            ! make sure velocities are neat.
            do i = 1, qp%n_irr_point
                f0 = norm2(f(:, i))
                if (f0 .gt. 1E-9_r8) then
                    v(:, i) = (1.0_r8 - fire%alpha)*v(:, i) + fire%alpha*norm2(v(:, i))*f(:, i)/f0
                else
                    v(:, i) = 0.0_r8
                end if
                select case (plctr(i))
                case (1)
                    ! on a plane
                    j = planeind(i)
                    fx = dot_product(p%irrw%face(j)%plane%v1, v(:, i))
                    fy = dot_product(p%irrw%face(j)%plane%v2, v(:, i))
                    v(:, i) = fx*p%irrw%face(j)%plane%v1 + fy*p%irrw%face(j)%plane%v2
                case (2)
                    v(:, i) = dot_product(linevecs(:, i), v(:, i))*linevecs(:, i)
                    !v(:,i)=0.0_r8
                case (3)
                    v(:, i) = 0.0_r8
                end select
            end do

            ! Fire thingy
            fire%P = 0.0_r8
            do i = 1, qp%n_irr_point
                fire%P = fire%P + dot_product(f(:, i), v(:, i))
            end do

            ! Perhaps an emergency freeze?
            if (maxval(abs(v))*fire%ts/idealside .gt. 0.15_r8) fire%P = -1.0_r8

            ! What to do about this?
            if (fire%P .gt. 0) then
                ! another step where P>0
                fire%counter = fire%counter + 1
                if (fire%counter .gt. fire%nmin) then
                    ! it was sufficiently long since we did a reset, increase the timestep
                    fire%ts = min(fire%ts*fire%ftsinc, fire%timestepmax)
                    ! decrease mixing
                    fire%alpha = fire%alpha*fire%falpha
                end if
            else
                ! do a reset, first the counter
                fire%counter = 0
                ! freeze atoms
                v = 0.0_r8
                ! decrease timestep
                fire%ts = max(fire%ts*fire%ftsdec, fire%timestepmin)
                ! reset alpha
                fire%alpha = fire%alpha0
                ! count number of resets
                rsctr = rsctr + 1
            end if

            f0 = sum(abs(f))/(3*qp%n_irr_point)
            f1 = maxval(abs(f))
            if (verbosity .gt. 0) then
            if (lo_trueNtimes(iter, 200, fire%niter) .or. iter .eq. 1) then
                write (*, '(1X,I4,2(1X,E11.3),3(1X,F10.5),1X,I4,1X,E11.3)') iter, f0, f1, mindist/idealside, maxdist/idealside, maxdr/idealside, rsctr, fire%ts
            end if
            end if

            ! check if converged
            if (f0 .lt. fire%tolerance/10.0_r8 .and. f1 .lt. fire%tolerance) then
                exit relaxloop
            end if
        end do relaxloop
    end block optmesh

    ! Retesselate to be on the safe side.
    call qp%tesselate_wedge_mesh(p, r, idealside*0.25_r8, -1.0_r8, verbosity=verbosity, mw=mw, mem=mem)
    ! And update the weights
    call qp%update_integration_weight(p, mw, mem, tolerance, verbosity)

    if (verbosity .gt. 0) then
        write (*, *) 'done optimizing (', tochar(walltime() - timer), ')'
    end if
end subroutine

!> Calculate the forces on the vertices
subroutine forces_on_vertices(qp, p, f, r, r0, plctr, planeind, linevecs, idealside, maxdist_per_point, mindist, maxdist, maxdr, mw)
    !> mesh
    type(lo_wedge_mesh), intent(in) :: qp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> forces
    real(r8), dimension(:, :), intent(inout) :: f
    !> current positions
    real(r8), dimension(:, :), intent(in) :: r
    !> original positions
    real(r8), dimension(:, :), intent(in) :: r0
    !> helper array to protect features
    integer, dimension(:), intent(in) :: plctr
    !> helper array with face indices
    integer, dimension(:), intent(in) :: planeind
    !> helper array with line vectors
    real(r8), dimension(:, :), intent(in) :: linevecs
    !> side of tetrahedron if the tesselation is ideal
    real(r8), intent(in) :: idealside
    !> max distance each point is allowed to move
    real(r8), dimension(:), intent(in) :: maxdist_per_point
    !> smallest distance between two points
    real(r8), intent(out) :: mindist
    !> longest edge in tetrahedron
    real(r8), intent(out) :: maxdist
    !> largest distance travelled from where it started
    real(r8), intent(out) :: maxdr
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    real(r8), parameter :: confinerad = 0.05_r8
    integer, parameter :: rexp = 2
    type(lo_verletbox) :: pbox
    real(r8), dimension(:, :), allocatable :: rall, rirr
    real(r8), dimension(3, 4) :: corner
    real(r8), dimension(3, 3) :: m0
    real(r8), dimension(3) :: v0, v1, v2
    real(r8) :: cutoff, cutoffsq, cospref, pfpref, cfpref, rtol, sidepref
    real(r8) :: f0, f1, fx, fy, rsq, ft1, ft2
    integer :: i, j, k, l, ii, jj, kk, bi, bj, bk, ci, cj, ck, t, nbox

    ! Cutoff distance (pretty short, but should be ok).
    f0 = lo_inscribed_sphere_in_box(p%reciprocal_latticevectors)*0.99_r8
    cutoff = min(f0, 1.2_r8*idealside)
    cutoffsq = cutoff**2
    rtol = lo_tol*(idealside**2)
    cospref = lo_pi/cutoffsq
    pfpref = idealside**2*1E-1_r8
    cfpref = (idealside*sqrt(6.0_r8)*0.25_r8)**2

    ! First step, rotate out all the points
    lo_allocate(rall(3, qp%n_full_point))
    lo_allocate(rirr(3, qp%n_irr_point))
    rall = 0.0_r8
    rirr = 0.0_r8
    do i = 1, qp%n_full_point
        if (mod(i, mw%n) .ne. mw%r) cycle
        k = qp%ap(i)%irreducible_index
        l = qp%ap(i)%operation_from_irreducible
        if (l .gt. 0) then
            rall(:, i) = matmul(p%sym%op(l)%m, r(:, k))
        else
            rall(:, i) = -matmul(p%sym%op(-l)%m, r(:, k))
        end if
        rall(:, i) = matmul(p%inv_reciprocal_latticevectors, rall(:, i))
        rall(:, i) = lo_clean_fractional_coordinates(rall(:, i))
    end do

    ! Get the irreducible guys in fractional coordinates
    do i = 1, qp%n_irr_point
        if (mod(i, mw%n) .ne. mw%r) cycle
        rirr(:, i) = matmul(p%inv_reciprocal_latticevectors, r(:, i))
    end do

    call mw%allreduce('sum', rirr)
    call mw%allreduce('sum', rall)

    ! Get the longest and shortest tetrahedron edge
    mindist = lo_huge
    maxdist = -lo_huge
    do i = 1, qp%n_irr_tet
        if (mod(i, mw%n) .ne. mw%r) cycle
        f0 = norm2(r(:, qp%it(i)%irreducible_index(2)) - r(:, qp%it(i)%irreducible_index(1)))
        mindist = min(f0, mindist)
        maxdist = max(f0, maxdist)
        f0 = norm2(r(:, qp%it(i)%irreducible_index(3)) - r(:, qp%it(i)%irreducible_index(1)))
        mindist = min(f0, mindist)
        maxdist = max(f0, maxdist)
        f0 = norm2(r(:, qp%it(i)%irreducible_index(4)) - r(:, qp%it(i)%irreducible_index(1)))
        mindist = min(f0, mindist)
        maxdist = max(f0, maxdist)
        f0 = norm2(r(:, qp%it(i)%irreducible_index(3)) - r(:, qp%it(i)%irreducible_index(2)))
        mindist = min(f0, mindist)
        maxdist = max(f0, maxdist)
        f0 = norm2(r(:, qp%it(i)%irreducible_index(4)) - r(:, qp%it(i)%irreducible_index(2)))
        mindist = min(f0, mindist)
        maxdist = max(f0, maxdist)
        f0 = norm2(r(:, qp%it(i)%irreducible_index(4)) - r(:, qp%it(i)%irreducible_index(3)))
        mindist = min(f0, mindist)
        maxdist = max(f0, maxdist)
    end do
    call mw%allreduce('min', mindist)
    call mw%allreduce('max', maxdist)

    ! Check if I can subdivide the reciprocal lattice so that Verlet boxes become useful
    nbox = 0
    do i = 1, 100
        m0 = p%reciprocal_latticevectors/i
        f0 = lo_inscribed_sphere_in_box(m0)
        if (f0 .gt. cutoff) then
            nbox = i
        else
            exit
        end if
    end do

    ! Reset forces
    f = 0.0_r8
    ! Loop over all pairs, hopefully using Verlet boxes
    if (nbox .ge. 4) then
        ! Now it's worth using Verlet boxes?
        call pbox%generate(rall, [nbox, nbox, nbox])
        ! Get all the forces
        do i = 1, qp%n_irr_point
            if (mod(i, mw%n) .ne. mw%r) cycle
            ! Reference vector
            v1 = rirr(:, i)
            v2 = lo_clean_fractional_coordinates(v1)
            ! First points, which box am I in?
            call pbox%boxind(v2, bi, bj, bk)
            do ii = -1, 1
            do jj = -1, 1
            do kk = -1, 1
                ci = lo_index_in_periodic_array(bi + ii, nbox)
                cj = lo_index_in_periodic_array(bj + jj, nbox)
                ck = lo_index_in_periodic_array(bk + kk, nbox)
                do l = 1, pbox%box(ci, cj, ck)%n
                    j = pbox%box(ci, cj, ck)%ind(l)
                    ! Now proceed with force calculation
                    v0 = mod(rall(:, j) - v1 + 0.5_r8, 1.0_r8) - 0.5_r8
                    v0 = matmul(p%reciprocal_latticevectors, v0)
                    rsq = lo_sqnorm(v0)
                    ! check the cutoffs
                    if (rsq .lt. rtol) cycle
                    if (rsq .gt. cutoffsq) cycle
                    ! Calculate an actual force, first a Cos term that smoothly cuts it off
                    ft1 = (cos(rsq*cospref) + 1.0_r8)*0.5_r8
                    ! Then the 1/r^n thing
                    ft2 = (pfpref/rsq)**rexp
                    ! Add this to the force
                    f(:, i) = f(:, i) - ft1*ft2*v0
                    ! Also grab the shortest distance
                end do
            end do
            end do
            end do
        end do
    else
        ! No boxes, just loop
        do i = 1, qp%n_irr_point
            if (mod(i, mw%n) .ne. mw%r) cycle
            ! Reference vector
            v1 = rirr(:, i)
            ! First points
            do j = 1, qp%n_full_point
                ! Now proceed with force calculation
                v0 = mod(rall(:, j) - v1 + 0.5_r8, 1.0_r8) - 0.5_r8
                v0 = matmul(p%reciprocal_latticevectors, v0)
                rsq = lo_sqnorm(v0)
                ! check the cutoffs
                if (rsq .lt. rtol) cycle
                if (rsq .gt. cutoffsq) cycle
                ! Calculate an actual force, first a Cos term that smoothly cuts it off
                ft1 = (cos(rsq*cospref) + 1.0_r8)*0.5_r8
                ! Then the 1/r^n thing
                ft2 = (pfpref/rsq)**rexp
                ! Add this to the force
                f(:, i) = f(:, i) - ft1*ft2*v0
            end do
        end do
    end if

    ! This is a punishment term for edges that deviate from the ideal.
    sidepref = -1.0_r8
    do t = 1, qp%n_irr_tet
        if (mod(t, mw%n) .ne. mw%r) cycle
        ! what is the ideal side for this tetrahedron?
        ! corner=0.0_r8
        ! do i=1,4
        !     corner(:,i)=qp%it(t)%irreducible_index(i)
        ! enddo
        ! f1=lo_unsigned_tetrahedron_volume(corner)
        ! f1=(f1*6.0_r8*sqrt(2.0_r8))**(1.0_r8/3.0_r8)

        ii = 2; jj = 1; i = qp%it(t)%irreducible_index(ii); j = qp%it(t)%irreducible_index(jj)
        v0 = r(:, i)
        v1 = r(:, j)
        v2 = (v0 + v1)*0.5_r8
        f0 = 2*(norm2(v1 - v0) - idealside)/idealside
        f0 = f0**3
        !f0=max(lo_sqnorm(v0-v1)-pfpref,0.0_r8)
        f(:, i) = f(:, i) + (v0 - v2)*f0*sidepref
        f(:, j) = f(:, j) + (v1 - v2)*f0*sidepref

        ii = 3; jj = 1; i = qp%it(t)%irreducible_index(ii); j = qp%it(t)%irreducible_index(jj)
        v0 = r(:, i)
        v1 = r(:, j)
        v2 = (v0 + v1)*0.5_r8
        f0 = 2*(norm2(v1 - v0) - idealside)/idealside
        f0 = f0**3
        !f0=max(lo_sqnorm(v0-v1)-pfpref,0.0_r8)
        f(:, i) = f(:, i) + (v0 - v2)*f0*sidepref
        f(:, j) = f(:, j) + (v1 - v2)*f0*sidepref

        ii = 4; jj = 1; i = qp%it(t)%irreducible_index(ii); j = qp%it(t)%irreducible_index(jj)
        v0 = r(:, i)
        v1 = r(:, j)
        v2 = (v0 + v1)*0.5_r8
        f0 = 2*(norm2(v1 - v0) - idealside)/idealside
        f0 = f0**3
        !f0=max(lo_sqnorm(v0-v1)-pfpref,0.0_r8)
        f(:, i) = f(:, i) + (v0 - v2)*f0*sidepref
        f(:, j) = f(:, j) + (v1 - v2)*f0*sidepref

        ii = 3; jj = 2; i = qp%it(t)%irreducible_index(ii); j = qp%it(t)%irreducible_index(jj)
        v0 = r(:, i)
        v1 = r(:, j)
        v2 = (v0 + v1)*0.5_r8
        f0 = 2*(norm2(v1 - v0) - idealside)/idealside
        f0 = f0**3
        !f0=max(lo_sqnorm(v0-v1)-pfpref,0.0_r8)
        f(:, i) = f(:, i) + (v0 - v2)*f0*sidepref
        f(:, j) = f(:, j) + (v1 - v2)*f0*sidepref

        ii = 4; jj = 2; i = qp%it(t)%irreducible_index(ii); j = qp%it(t)%irreducible_index(jj)
        v0 = r(:, i)
        v1 = r(:, j)
        v2 = (v0 + v1)*0.5_r8
        f0 = 2*(norm2(v1 - v0) - idealside)/idealside
        f0 = f0**3
        !f0=max(lo_sqnorm(v0-v1)-pfpref,0.0_r8)
        f(:, i) = f(:, i) + (v0 - v2)*f0*sidepref
        f(:, j) = f(:, j) + (v1 - v2)*f0*sidepref

        ii = 4; jj = 3; i = qp%it(t)%irreducible_index(ii); j = qp%it(t)%irreducible_index(jj)
        v0 = r(:, i)
        v1 = r(:, j)
        v2 = (v0 + v1)*0.5_r8
        f0 = 2*(norm2(v1 - v0) - idealside)/idealside
        f0 = f0**3
        !f0=max(lo_sqnorm(v0-v1)-pfpref,0.0_r8)
        f(:, i) = f(:, i) + (v0 - v2)*f0*sidepref
        f(:, j) = f(:, j) + (v1 - v2)*f0*sidepref
    end do

    ! Add a very agressive force that confines the points so that they don't wander off too far.
    maxdr = 0.0_r8
    do i = 1, qp%n_irr_point
        if (mod(i, mw%n) .ne. mw%r) cycle
        v0 = r(:, i) - r0(:, i)
        f0 = norm2(v0)
        maxdr = max(maxdr, f0)

        if (f0 .gt. maxdist_per_point(i)) then
            f1 = f0 - maxdist_per_point(i)
            f(:, i) = f(:, i) - 1E5_r8*v0*f1**2
            f(:, i) = f(:, i) - 1E8_r8*v0*f1**4
        end if
    end do

    ! Decent time to sync the forces
    call mw%allreduce('sum', f)

    ! then make sure the points on the sharp features are preserved
    do i = 1, qp%n_irr_point
        select case (plctr(i))
        case (1)
            ! on a plane
            j = planeind(i)
            fx = dot_product(p%irrw%face(j)%plane%v1, f(:, i))
            fy = dot_product(p%irrw%face(j)%plane%v2, f(:, i))
            f(:, i) = (fx*p%irrw%face(j)%plane%v1 + fy*p%irrw%face(j)%plane%v2)
        case (2)
            !f(:,i)=dot_product(linevecs(:,i),f(:,i))*linevecs(:,i)
            f(:, i) = 0.0_r8
        case (3)
            f(:, i) = 0.0_r8
        end select
    end do
end subroutine

!> dump histograms over tetrahedron volumes
subroutine get_histograms(qp, vx, vy, sx, sy, idealvolume, idealside, mw)
    !> q-point mesh
    type(lo_wedge_mesh), intent(in) :: qp
    !> x and y-axis for histogram
    real(r8), intent(inout), dimension(:) :: vx, vy
    real(r8), intent(inout), dimension(:) :: sx, sy
    !> limits to the histogram
    real(r8), intent(in) :: idealvolume, idealside
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    real(r8), parameter :: hfac = 4.0_r8
    real(r8), dimension(3, 4) :: dumtet
    real(r8), dimension(:), allocatable :: dv, dw
    real(r8) :: f0, f1, invf_vol, invf_side, vmax, smax
    integer :: n, i, j, k, ii, jj, kk, ll
    ! First set the range. Hmmm.
    n = size(vx) + 1 ! number of x-values
    vmax = idealvolume*hfac
    smax = idealside*hfac
    vmax = hfac
    smax = hfac
    ! Set the values for the center of the bins
    do i = 1, n - 1
        vx(i) = vmax*(i - 0.5_r8)/real(n - 1, r8)
        sx(i) = smax*(i - 0.5_r8)/real(n - 1, r8)
    end do
    invf_vol = real(n, r8)/hfac
    invf_side = real(n, r8)/hfac

    allocate (dv(qp%n_irr_tet))
    allocate (dw(qp%n_irr_point))
    dv = 0.0_r8
    dw = 1E10_r8

    vy = 0.0_r8
    sy = 0.0_r8
    do i = 1, qp%n_irr_tet
        if (mod(i, mw%n) .ne. mw%r) cycle
        do j = 1, 4
            dumtet(:, j) = qp%ip(qp%it(i)%irreducible_index(j))%r
        end do
        f0 = lo_unsigned_tetrahedron_volume(dumtet)/idealvolume
        dv(i) = f0
        j = floor(f0*invf_vol) + 1
        j = max(j, 1)
        j = min(j, n - 1)
        vy(j) = vy(j) + 1
        do ii = 1, 4
            do jj = ii + 1, 4
                f0 = norm2(dumtet(:, ii) - dumtet(:, jj))/idealside
                j = floor(f0*invf_side) + 1
                j = max(j, 1)
                j = min(j, n - 1)
                sy(j) = sy(j) + 1

                kk = qp%it(i)%irreducible_index(ii)
                ll = qp%it(i)%irreducible_index(jj)
                dw(kk) = min(dw(kk), f0)
                dw(ll) = min(dw(ll), f0)
            end do
        end do
    end do
    call mw%allreduce('sum', vy)
    call mw%allreduce('sum', sy)
    call mw%allreduce('sum', dv)
    call mw%allreduce('min', dw)
    vy = vy/lo_trapezoid_integration(vx, vy)
    sy = sy/lo_trapezoid_integration(sx, sy)

    call lo_qsort(dv)
    call lo_qsort(dw)
    ! List the smallest and largest tetrahedrons?
    if (mw%talk) then
        do i = 1, 10
            write (*, *) i, dv(i), dw(i)
        end do
        do i = 1, 10
            write (*, *) i, dv(qp%n_irr_tet - 10 + i), dw(qp%n_irr_point - 10 + i)
        end do
    end if
end subroutine

!> make a fake integration
subroutine fake_integration(uc, qp, m0)
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> integral thingy
    real(r8), dimension(3, 3), intent(out) :: m0

    real(r8), dimension(:, :), allocatable :: fval
    real(r8), dimension(3) :: v0
    real(r8) :: rsq, f0, f1, prefactor
    integer :: i

    allocate (fval(3, qp%n_full_point))
    fval = 0.0_r8

    prefactor = 2*(lo_pi**2 - 6)*(uc%bz%rmin**3)/(9*lo_pi)/abs(lo_determ(uc%reciprocal_latticevectors))
    prefactor = 1.0_r8/prefactor

    rsq = uc%bz%rmin
    do i = 1, qp%n_full_point
        v0 = qp%ap(i)%r
        f1 = norm2(v0)
        if (f1 .le. rsq) then
            f0 = cos(lo_pi*0.5_r8*norm2(v0)/rsq)
        else
            f0 = 0.0_r8
        end if

        if (f1 .gt. lo_sqtol) then
            fval(:, i) = f0*v0/f1
        else
            fval(:, i) = 0.0_r8
        end if
    end do

    m0 = 0.0_r8
    do i = 1, qp%n_full_point
        v0 = fval(:, i)
        m0 = m0 + lo_outerproduct(v0, v0)*qp%ap(i)%integration_weight
    end do
    m0 = m0*prefactor
end subroutine

end module
