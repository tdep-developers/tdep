#include "precompilerdefinitions"
module type_distancetable
!! Type to handle calculation of distance tables, with or without periodic boundary conditions for arbitrary sets of points
use konstanter, only: flyt,lo_huge,lo_hugeint,lo_tol,lo_sqtol,&
                      lo_status,lo_exitcode_param,lo_exitcode_symmetry
use gottochblandat, only: tochar,walltime,lo_clean_fractional_coordinates,lo_chop,qsort,lo_stop_gracefully,&
                          lo_sqnorm,lo_cross,lo_invert3x3matrix
use geometryfunctions, only: lo_inscribed_sphere_in_box,lo_bounding_sphere_of_box,lo_plane,lo_increment_dimensions
use mpi_wrappers, only: lo_mpi_helper
use type_blas_lapack_wrappers, only: lo_gemm
implicit none
private
public :: lo_distancetable
public :: lo_distancetable_particle

!> particle in a distance table
type lo_distancetable_particle
    !> how many neighbours are there?
    integer :: n=-lo_hugeint
    !> What is the vector pointing towards these particles?
    real(flyt), dimension(:,:), allocatable :: v
    !> Vector that wraps the periodic boundary things
    real(flyt), dimension(:,:), allocatable :: lv
    !> What is the index in the original cell of the neighbouring particle?
    integer, dimension(:), allocatable :: ind
    !> What is the distance
    real(flyt), dimension(:), allocatable :: d
    !> What is the weight of this point, in the Wigner-Seitz scheme
    real(flyt), dimension(:), allocatable :: weight
    !> How many nearest neighbours jumps is it here?
    integer, dimension(:), allocatable :: jumps
end type

!> A distance table for a set of points. Handles periodic and non-periodic boundary conditions, and should be reasonably fast. Not classical-MD-fast, but fast enough for most uses. With sufficiently many points and short cutoff, it automatically switches to O(N) Verlet lists. Or at least it will, once I need it.
type lo_distancetable
    !> Number of particles
    integer :: np=-lo_hugeint
    !> The cutoff
    real(flyt) :: cutoff=-lo_huge
    !> The list originating from each particle
    type(lo_distancetable_particle), dimension(:), allocatable :: particle
    contains
        !> create the distance table
        procedure :: generate
        !> largest number of neighbours
        procedure :: max_number_of_neighbours
        !> smallest number of neighbours
        procedure :: min_number_of_neighbours
        !> count neighbours within cutoff
        procedure :: count_neighbours_within_cutoff
        !> remove neighours outside some cutoff
        procedure :: prune
end type

contains

!> generate distance table
subroutine generate(dt,particles,basis,cutoff,verbosity,mw,tolerance)
    !> The distance table
    class(lo_distancetable), intent(out) :: dt
    !> The list of particles in fractional coordinates
    real(flyt), dimension(:,:), intent(in) :: particles
    !> The basis vectors for the particles
    real(flyt), dimension(3,3), intent(in) :: basis
    !> Maximum distance to consider
    real(flyt), intent(in) :: cutoff
    !> How much to talk
    integer, intent(in) :: verbosity
    !> MPI helper
    type(lo_mpi_helper), intent(inout), optional :: mw
    !> Tolerance
    real(flyt), intent(in), optional :: tolerance

    integer, parameter :: verletlist_crossover=800 ! with fever particles than this it's not worth using a verlet list.
    integer, parameter :: maxnbox=50 ! largest number of Verlet boxes. Should be plenty for anything reasonable.
    real(flyt) :: cutoffbuf  ! it is bad if the cutoff is exactly on a shell, so buffer it a little.

    integer :: np,algo
    integer, dimension(3) :: nrep,nbox
    real(flyt), dimension(:,:), allocatable :: cleanpositions,rcart
    real(flyt) :: timer,rc,sqrc,tol,sqtol
    logical :: verletlist,mpi

    ! set some basic stuff
    init: block
        real(flyt), dimension(3,3) :: m0
        integer :: i,j

        ! Start timer if talking
        if ( verbosity .gt. 0 ) then
            timer=walltime()
        endif

        ! Tolerance?
        if ( present(tolerance) ) then
            tol=tolerance
            sqtol=tol**2
        else
            tol=lo_tol
            sqtol=lo_sqtol
        endif

        ! Should it be done in parallel?
        if ( present(mw) ) then
            mpi=.true.
        else
            mpi=.false.
        endif

        ! number of particles
        np=size(particles,2)
        if ( verbosity .gt. 0 ) then
            write(*,*) '... building distance table for '//tochar(np)//' points'
        endif
        ! Add buffer to cutoff
        cutoffbuf=100*tol
        rc=cutoff+cutoffbuf
        sqrc=(cutoff+cutoffbuf)**2

        ! clean the positions on input to make sure they are [0-1[.
        lo_allocate(cleanpositions(3,np))
        cleanpositions=lo_clean_fractional_coordinates(particles)
        if ( verbosity .gt. 0 ) write(*,*) '... cleaned the positions'

#ifdef AGRESSIVE_SANITY
        ! make sure I don't have stupid cutoff
        if ( cutoff .gt. lo_inscribed_sphere_in_box(basis)*1E4_flyt ) then
            call lo_stop_gracefully(['Trying to calculate distances with very long cutoff, probably not intended.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif

        ! make sure the number of particles is not too large?
        if ( np .gt. 2000 ) then
            call lo_stop_gracefully(['Trying to calculate distance table for more than 2000 particles. Nag on me to fix.'],&
                                    lo_exitcode_param,__FILE__,__LINE__)
        endif
#endif

        ! figure out how many times things need to be repeated.
        nrep=0
        if ( lo_inscribed_sphere_in_box(basis) .lt. rc ) then
            do j=1,100
                do i=1,3
                    m0(:,i)=basis(:,i)*(2*nrep(i)+1)
                enddo
                if ( lo_inscribed_sphere_in_box(m0)-lo_bounding_sphere_of_box(basis) .gt. cutoff+cutoffbuf ) exit
                nrep=lo_increment_dimensions(nrep,basis)
            enddo
        endif

        ! Now, if no repetitions are needed, it might be worth trying a Verlet-list type calculation of distances.
        ! For that to be any kind of speedup, I must be able to divide the box into at least 4 sub-boxes, where
        ! each sub-box must be large enough to contain the cutoff. Also, the number of particles needs to be large
        ! enough to validate the extra overhead of switching to an O(N) algorithm.
        if ( np .gt. verletlist_crossover ) then
            nbox=0
            do i=1,3
                do j=maxnbox,1,-1
                    m0=basis
                    m0(:,i)=m0(:,i)/j
                    if ( lo_inscribed_sphere_in_box(m0) .gt. rc ) then
                        nbox(i)=j
                        exit
                    endif
                enddo
            enddo
#ifdef AGRESSIVE_SANITY
            if ( maxval(nbox) .eq. maxnbox ) then
                call lo_stop_gracefully(['Hit the maximum number of Verlet boxes. This is not reasonable for anything I intended.'],&
                                         lo_exitcode_param,__FILE__,__LINE__)
            endif
#endif
            ! Now I can check if the divisions are large enough to use.
            if ( minval(nbox) .gt. 3 ) then
                verletlist=.true.
            else
                verletlist=.false.
            endif
        else
            ! no point in even trying with this few particles.
            verletlist=.false.
        endif

        ! Ok, now I have enough information to determine wich algorithm to use.
        algo=0
        if ( sum(nrep) .eq. 0 ) then
            if ( mpi ) then
                if ( verletlist ) then
                    !algo=-1
                else
                    algo=3
                endif
            else
                if ( verletlist ) then
                    !algo=-1
                else
                    algo=4
                endif
            endif
        else
            if ( mpi ) then
                algo=1
            else
                algo=2
            endif
        endif

        if ( algo .eq. 0 ) then
            call lo_stop_gracefully(['Thinking is hard when it comes to distance tables.'],&
                                     lo_exitcode_param,__FILE__,__LINE__)
        endif

        ! For some algorithms, I pre-convert to Cartesian coordinates
        if ( algo .le. 2 ) then
            lo_allocate(rcart(3,np))
            call lo_gemm(basis,cleanpositions,rcart)
        endif

        if ( verbosity .gt. 0 ) then
            write(*,*) '... using',tochar(nrep),' repetitions to build the distance table'
        endif
    end block init

    ! Now for the individual algorithms
    select case(algo)
    case(1)
    repalgo: block
        ! Use filtered repetitions
        real(flyt), dimension(:,:,:), allocatable :: dumr,dumlv
        real(flyt), dimension(:,:), allocatable :: dumd,latvec,llv
        real(flyt), dimension(:), allocatable :: dr
        real(flyt), dimension(3) :: v0,v1,v2
        real(flyt) :: buf
        integer, dimension(:,:), allocatable :: dumind,dumjnd
        integer, dimension(:), allocatable :: dumctr,atind,lai,rnkctr,offset
        integer :: maxnn,ncells,ctr,ii,jj,kk,a1,a2,l,i,j

        ! Parallel things
        lo_allocate(rnkctr(mw%n))
        lo_allocate(offset(mw%n))
        rnkctr=0
        offset=0

        ! First step is to generate a list of latticevectors and atoms, local
        ! to each rank
        ncells=0
        ctr=0
        do a1=1,np
            v0=rcart(:,a1)
            do ii=-nrep(1),nrep(1)
            do jj=-nrep(2),nrep(2)
            do kk=-nrep(3),nrep(3)
                ctr=ctr+1
                if ( mod(ctr,mw%n) .ne. mw%r ) cycle
                v1=matmul(basis,[ii,jj,kk]*1.0_flyt)
                ! First the fast test
                if ( shortest_distance_cell_to_point(v1,v0,basis) .gt. rc ) cycle
                ! Then angrier test
                ! if ( distance_point_box( v1,v0,basis ) .gt. rc ) cycle
                ncells=ncells+1
            enddo
            enddo
            enddo
        enddo

        lo_allocate(llv(3,ncells))
        lo_allocate(lai(ncells))
        llv=0.0_flyt
        lai=0
        ncells=0
        ctr=0
        do a1=1,np
            v0=rcart(:,a1)
            do ii=-nrep(1),nrep(1)
            do jj=-nrep(2),nrep(2)
            do kk=-nrep(3),nrep(3)
                ctr=ctr+1
                if ( mod(ctr,mw%n) .ne. mw%r ) cycle
                v1=matmul(basis,[ii,jj,kk]*1.0_flyt)
                if ( shortest_distance_cell_to_point(v1,v0,basis) .gt. rc ) cycle
                ! if ( distance_point_box( v1,v0,basis ) .gt. rc ) cycle
                ncells=ncells+1
                lai(ncells)=a1
                llv(:,ncells)=v1
            enddo
            enddo
            enddo
        enddo

        ! Figure out how many cells I have on each rank, and the offset
        rnkctr=0
        offset=0
        rnkctr(mw%r+1)=ncells
        call mw%allreduce('sum',rnkctr)
        do i=1,mw%n-1
            offset(i+1)=sum(rnkctr(1:i))
        enddo

        ! Make large arrays to store all cells on all ranks?
        lo_allocate(atind(sum(rnkctr)))
        lo_allocate(latvec(3,sum(rnkctr)))
        atind=0
        latvec=0.0_flyt
        ! Store the latticevectors and atom indices in the right place
        do i=1,ncells
            j=offset(mw%r+1)+i
            atind(j)=lai(i)
            latvec(:,j)=llv(:,i)
        enddo
        call mw%allreduce('sum',atind)
        call mw%allreduce('sum',latvec)
        lo_deallocate(lai)
        lo_deallocate(llv)
        ! And make sure total number of cells is known
        ncells=sum(rnkctr)

        ! That was the first level of parallelism. Now count the number of neighbours? Not sure what
        ! the best way is, perhaps a big flat array thing?
        lo_allocate(dumctr(np))
        dumctr=0
        do i=1,ncells
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            a1=atind(i)
            v0=rcart(:,a1)
            v1=latvec(:,i)
            do a2=1,np
                v2=rcart(:,a2)-v0+v1
                if ( lo_sqnorm(v2) .gt. sqrc ) cycle
                dumctr(a1)=dumctr(a1)+1
            enddo
        enddo

        ! This would be the size of my buffer thingies.
        ! Pretty sure it's ok, a little waste of memory
        ! but should be fine.
        maxnn=maxval(dumctr)
        call mw%allreduce('max',maxnn)
        do i=1,mw%n
            offset(i)=(i-1)*maxnn
        enddo
        maxnn=maxnn*mw%n
        lo_allocate(dumr(3,maxnn,np))
        lo_allocate(dumlv(3,maxnn,np))
        lo_allocate(dumd(maxnn,np))
        lo_allocate(dumind(maxnn,np))
        lo_allocate(dumjnd(maxnn,np))
        lo_allocate(dr(maxnn))
        dumr=0.0_flyt
        dumlv=0.0_flyt
        dumd=1E10_flyt*rc
        dumind=0
        dumjnd=0
        dumctr=0

        ! Now store vectors and stuff
        do i=1,ncells
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            a1=atind(i)
            v0=rcart(:,a1)
            v1=latvec(:,i)
            do a2=1,np
                v2=rcart(:,a2)-v0+v1
                if ( lo_sqnorm(v2) .gt. sqrc ) cycle
                dumctr(a1)=dumctr(a1)+1
                j=dumctr(a1)+offset(mw%r+1)
                dumr(:,j,a1)=v2
                dumlv(:,j,a1)=v1
                dumd(j,a1)=norm2( v2 )
                dumind(j,a1)=a2
            enddo
        enddo
        ! Communicate everywhere
        call mw%allreduce('sum',dumr)
        call mw%allreduce('sum',dumlv)
        call mw%allreduce('min',dumd)
        call mw%allreduce('sum',dumind)
        call mw%allreduce('sum',dumctr)

        ! Sort things by distance.
        dumjnd=0
        do a1=1,np
            if ( mod(a1,mw%n) .ne. mw%r ) cycle
            dr=dumd(:,a1)
            call qsort(dr,dumjnd(:,a1))
        enddo
        call mw%allreduce('sum',dumjnd)

        ! Now make sure the cutoff makes sense, or if I have to increase it a little
        ! to get off of peaks of the rdf. Think it's ok to do this serially for now.
        buf=0.0_flyt
        do j=1,100
            l=0
            do a1=1,np
            do i=1,dumctr(a1)
                if ( abs( cutoff+buf-dumd(dumjnd(i,a1),a1) ) .lt. 2*tol ) then
                    l=l+1
                endif
            enddo
            enddo
            if ( l .eq. 0 ) then
                exit
            else
                buf=buf+tol
            endif
        enddo
        if ( verbosity .gt. 0 .and. buf .gt. sqtol ) then
            write(*,*) '... increased cutoff a little to accomodate full shells'
        endif

        ! now store everything
        dt%np=np
        dt%cutoff=0.0_flyt
        lo_allocate(dt%particle(dt%np))
        do a1=1,np
            l=0
            do i=1,dumctr(a1)
                j=dumjnd(i,a1)
                if ( dumd(j,a1) .lt. cutoff+buf ) l=l+1
            enddo
            dt%particle(a1)%n=l
            lo_allocate(dt%particle(a1)%d(l))
            lo_allocate(dt%particle(a1)%v(3,l))
            lo_allocate(dt%particle(a1)%lv(3,l))
            lo_allocate(dt%particle(a1)%ind(l))
            l=0
            do i=1,dumctr(a1)
                j=dumjnd(i,a1)
                if ( dumd(j,a1) .lt. cutoff+buf ) then
                    l=l+1
                    dt%particle(a1)%d(l)=lo_chop( dumd(j,a1), sqtol )
                    dt%particle(a1)%v(:,l)=lo_chop( dumr(:,j,a1), sqtol )
                    dt%particle(a1)%lv(:,l)=lo_chop( dumlv(:,j,a1), sqtol )
                    dt%particle(a1)%ind(l)=dumind(j,a1)
                endif
            enddo
            dt%cutoff=max(dt%cutoff,dt%particle(a1)%d(l))
        enddo
        ! And some cleanup
        lo_deallocate(dumd)
        lo_deallocate(dumr)
        lo_deallocate(dumlv)
        lo_deallocate(dumind)
        lo_deallocate(dumjnd)
        lo_deallocate(dumctr)
        lo_deallocate(rnkctr)
        lo_deallocate(offset)
        lo_deallocate(dr)
    end block repalgo
    case(2)
    repalgoser: block
        ! Serial version of filtered repetitions
        real(flyt), dimension(:,:,:), allocatable :: dumr
        real(flyt), dimension(:,:,:), allocatable :: dumlv
        real(flyt), dimension(:,:), allocatable :: dumd,latvec
        real(flyt), dimension(:), allocatable :: dr
        real(flyt), dimension(3) :: v0,v1,v2
        real(flyt) :: buf
        integer, dimension(:,:), allocatable :: dumind,dumjnd
        integer, dimension(:), allocatable :: dumctr,atind
        integer :: maxnn,ncells,a1,a2,i,j,ii,jj,kk,l

        ! First step is to generate a list of latticevectors and atoms, local
        ! to each rank
        ncells=0
        do a1=1,np
            v0=rcart(:,a1)
            do ii=-nrep(1),nrep(1)
            do jj=-nrep(2),nrep(2)
            do kk=-nrep(3),nrep(3)
                v1=matmul(basis,[ii,jj,kk]*1.0_flyt)
                ! First the fast test
                if ( shortest_distance_cell_to_point(v1,v0,basis) .gt. rc ) cycle
                ! Then angrier test
                ! if ( distance_point_box( v1,v0,basis ) .gt. rc ) cycle
                ncells=ncells+1
            enddo
            enddo
            enddo
        enddo

        lo_allocate(latvec(3,ncells))
        lo_allocate(atind(ncells))
        latvec=0.0_flyt
        atind=0
        ncells=0
        do a1=1,np
            v0=rcart(:,a1)
            do ii=-nrep(1),nrep(1)
            do jj=-nrep(2),nrep(2)
            do kk=-nrep(3),nrep(3)
                v1=matmul(basis,[ii,jj,kk]*1.0_flyt)
                if ( shortest_distance_cell_to_point(v1,v0,basis) .gt. rc ) cycle
                ! if ( distance_point_box( v1,v0,basis ) .gt. rc ) cycle
                ncells=ncells+1
                atind(ncells)=a1
                latvec(:,ncells)=v1
            enddo
            enddo
            enddo
        enddo

        ! That was the first level of parallelism. Now count the number of neighbours? Not sure what
        ! the best way is, perhaps a big flat array thing?
        lo_allocate(dumctr(np))
        dumctr=0
        do i=1,ncells
            a1=atind(i)
            v0=rcart(:,a1)
            v1=latvec(:,i)
            do a2=1,np
                v2=rcart(:,a2)-v0+v1
                if ( lo_sqnorm(v2) .gt. sqrc ) cycle
                dumctr(a1)=dumctr(a1)+1
            enddo
        enddo

        ! This would be the size of my buffer thingies.
        ! Pretty sure it's ok, a little waste of memory
        ! but should be fine.
        maxnn=maxval(dumctr)
        lo_allocate(dumr(3,maxnn,np))
        lo_allocate(dumlv(3,maxnn,np))
        lo_allocate(dumd(maxnn,np))
        lo_allocate(dumind(maxnn,np))
        lo_allocate(dumjnd(maxnn,np))
        lo_allocate(dr(maxnn))
        dumr=0.0_flyt
        dumlv=0.0_flyt
        dumd=1E10_flyt*rc
        dumind=0
        dumjnd=0
        dumctr=0

        ! Now store vectors and stuff
        do i=1,ncells
            a1=atind(i)
            v0=rcart(:,a1)
            v1=latvec(:,i)
            do a2=1,np
                v2=rcart(:,a2)-v0+v1
                if ( lo_sqnorm(v2) .gt. sqrc ) cycle
                dumctr(a1)=dumctr(a1)+1
                j=dumctr(a1)
                dumr(:,j,a1)=v2
                dumlv(:,j,a1)=v1
                dumd(j,a1)=norm2( v2 )
                dumind(j,a1)=a2
            enddo
        enddo
        ! Sort things by distance
        dumjnd=0
        do a1=1,np
            dr=dumd(:,a1)
            call qsort(dr,dumjnd(:,a1))
        enddo

        ! Now make sure the cutoff makes sense, or if I have to increase it a little
        ! to get off of peaks of the rdf. Think it's ok to do this serially for now.
        buf=0.0_flyt
        do j=1,100
            l=0
            do a1=1,np
            do i=1,dumctr(a1)
                if ( abs( cutoff+buf-dumd(dumjnd(i,a1),a1) ) .lt. 2*tol ) then
                    l=l+1
                endif
            enddo
            enddo
            if ( l .eq. 0 ) then
                exit
            else
                buf=buf+lo_tol
            endif
        enddo
        if ( verbosity .gt. 0 .and. buf .gt. sqtol ) then
            write(*,*) '... increased cutoff a little to accomodate full shells'
        endif

        ! now store everything
        dt%np=np
        dt%cutoff=0.0_flyt
        lo_allocate(dt%particle(dt%np))
        do a1=1,np
            l=0
            do i=1,dumctr(a1)
                j=dumjnd(i,a1)
                if ( dumd(j,a1) .lt. cutoff+buf ) l=l+1
            enddo
            dt%particle(a1)%n=l
            lo_allocate(dt%particle(a1)%d(l))
            lo_allocate(dt%particle(a1)%v(3,l))
            lo_allocate(dt%particle(a1)%lv(3,l))
            lo_allocate(dt%particle(a1)%ind(l))
            l=0
            do i=1,dumctr(a1)
                j=dumjnd(i,a1)
                if ( dumd(j,a1) .lt. cutoff+buf ) then
                    l=l+1
                    dt%particle(a1)%d(l)=lo_chop( dumd(j,a1), sqtol )
                    dt%particle(a1)%v(:,l)=lo_chop( dumr(:,j,a1), sqtol )
                    dt%particle(a1)%lv(:,l)=lo_chop( dumlv(:,j,a1), sqtol )
                    dt%particle(a1)%ind(l)=dumind(j,a1)
                endif
            enddo
            dt%cutoff=max(dt%cutoff,dt%particle(a1)%d(l))
        enddo
        ! And some cleanup
        lo_deallocate(dumd)
        lo_deallocate(dumr)
        lo_deallocate(dumlv)
        lo_deallocate(dumind)
        lo_deallocate(dumjnd)
        lo_deallocate(dumctr)
        lo_deallocate(dr)
    end block repalgoser
    case(3)
    dirnn: block
        ! Direct calculation, in parallel
        real(flyt), dimension(:,:,:), allocatable :: dumr
        real(flyt), dimension(:,:,:), allocatable :: dumlv
        real(flyt), dimension(:,:), allocatable :: dumd
        real(flyt), dimension(3) :: v0,v1,v2
        real(flyt) :: buf
        integer, dimension(:,:), allocatable :: dumind
        integer, dimension(:), allocatable :: sortind,dumctr,offset
        integer :: maxnn,a1,a2,i,j,l

        ! Count to get an upper bound on the number of neighbours
        lo_allocate(dumctr(np))
        lo_allocate(offset(mw%n))
        dumctr=0
        offset=0
        do a1=1,np
            if ( mod(a1,mw%n) .ne. mw%r ) cycle
            l=0
            do a2=1,np
                ! pairvector without PBC-check
                v0=cleanpositions(:,a2)-cleanpositions(:,a1)
                ! pairvector with pbc-check
                v1=lo_clean_fractional_coordinates(v0+0.5_flyt)-0.5_flyt
                ! Convert to Cartesian
                v1=matmul(basis,v1)
                if ( lo_sqnorm(v1) .gt. sqrc ) cycle
                l=l+1
            enddo
            dumctr(a1)=l
        enddo

        ! Make some space for buffers
        maxnn=maxval(dumctr)
        call mw%allreduce('max',maxnn)
        do i=1,mw%n
            offset(i)=(i-1)*maxnn
        enddo
        maxnn=maxnn*mw%n
        lo_allocate(dumr(3,maxnn,np))
        lo_allocate(dumlv(3,maxnn,np))
        lo_allocate(dumd(maxnn,np))
        lo_allocate(dumind(maxnn,np))
        lo_allocate(sortind(maxnn))
        dumr=0.0_flyt
        dumlv=0.0_flyt
        dumd=1E10_flyt*rc
        dumind=0
        dumctr=0
        sortind=0

        do a1=1,np
            if ( mod(a1,mw%n) .ne. mw%r ) cycle
            l=0
            do a2=1,np
                ! pairvector without PBC-check
                v0=cleanpositions(:,a2)-cleanpositions(:,a1)
                ! pairvector with pbc-check
                v1=lo_clean_fractional_coordinates(v0+0.5_flyt)-0.5_flyt
                ! Convert to Cartesian
                v1=matmul(basis,v1)
                if ( lo_sqnorm(v1) .gt. sqrc ) cycle
                ! lattice vector
                v2=v1-v0
                v2=matmul(basis,anint(v2))
                l=l+1
                dumr(:,l,a1)=v1
                dumlv(:,l,a1)=v2
                dumd(l,a1)=norm2( v1 )
                dumind(l,a1)=a2
            enddo
            dumctr(a1)=l
            ! Sort by distance
            call qsort(dumd(1:l,a1),sortind(1:l))
            dumr(:,1:l,a1)=dumr(:,sortind(1:l),a1)
            dumlv(:,1:l,a1)=dumlv(:,sortind(1:l),a1)
            dumind(1:l,a1)=dumind(sortind(1:l),a1)
        enddo
        call mw%allreduce('sum',dumctr)
        call mw%allreduce('sum',dumr)
        call mw%allreduce('min',dumd)
        call mw%allreduce('sum',dumlv)
        call mw%allreduce('sum',dumind)

        ! Now make sure the cutoff makes sense, or if I have to increase it a little
        ! to get off of peaks of the rdf. Think it's ok to do this serially for now.
        buf=0.0_flyt
        do j=1,100
            l=0
            do a1=1,np
            do i=1,dumctr(a1)
                if ( abs( cutoff+buf-dumd(i,a1) ) .lt. 2*lo_tol ) then
                    l=l+1
                endif
            enddo
            enddo
            if ( l .eq. 0 ) then
                exit
            else
                buf=buf+lo_tol
            endif
        enddo

        ! now store everything
        dt%np=np
        dt%cutoff=0.0_flyt
        lo_allocate(dt%particle(dt%np))
        do a1=1,np
            l=0
            do i=1,dumctr(a1)
                if ( dumd(i,a1) .lt. cutoff+buf ) l=l+1
            enddo
            dt%particle(a1)%n=l
            lo_allocate(dt%particle(a1)%d(l))
            lo_allocate(dt%particle(a1)%v(3,l))
            lo_allocate(dt%particle(a1)%lv(3,l))
            lo_allocate(dt%particle(a1)%ind(l))
            dt%particle(a1)%d=lo_chop( dumd(1:l,a1),sqtol )
            dt%particle(a1)%v=lo_chop( dumr(:,1:l,a1),sqtol )
            dt%particle(a1)%lv=lo_chop( dumlv(:,1:l,a1),sqtol )
            dt%particle(a1)%ind=dumind(1:l,a1)
            dt%cutoff=max(dt%cutoff,dt%particle(a1)%d(l))
        enddo
        ! And some cleanup
        lo_deallocate(dumd)
        lo_deallocate(dumr)
        lo_deallocate(dumlv)
        lo_deallocate(dumind)
        lo_deallocate(dumctr)
    end block dirnn
    case(4)
    dirnnser: block
        ! Direct calculation, no repetition, serial
        integer, parameter :: countcrossover=100
        real(flyt), dimension(:,:,:), allocatable :: dumr,dumlv
        real(flyt), dimension(:,:), allocatable :: dumd
        real(flyt), dimension(3) :: v0,v1,v2
        real(flyt) :: buf
        integer, dimension(:,:), allocatable :: dumind
        integer, dimension(:), allocatable :: sortind,dumctr
        integer :: a1,a2,l,i,j,maxnn

        ! Upper bound for number of neighbours
        if ( np .gt. countcrossover ) then
            maxnn=0
            do a1=1,np
                l=0
                do a2=1,np
                    v0=cleanpositions(:,a2)-cleanpositions(:,a1)
                    v1=lo_clean_fractional_coordinates(v0+0.5_flyt)-0.5_flyt
                    v1=matmul(basis,v1)
                    if ( lo_sqnorm(v1) .gt. sqrc ) cycle
                    l=l+1
                enddo
                maxnn=max(maxnn,l)
            enddo
        else
            maxnn=np
        endif

        lo_allocate(dumctr(np))
        lo_allocate(dumr(3,maxnn,np))
        lo_allocate(dumlv(3,maxnn,np))
        lo_allocate(dumd(maxnn,np))
        lo_allocate(dumind(maxnn,np))
        lo_allocate(sortind(maxnn))
        dumr=0.0_flyt
        dumlv=0.0_flyt
        dumd=1E10_flyt*rc
        dumind=0
        dumctr=0
        sortind=0
        do a1=1,np
            l=0
            do a2=1,np
                ! pairvector without PBC-check
                v0=cleanpositions(:,a2)-cleanpositions(:,a1)
                ! pairvector with pbc-check
                v1=lo_clean_fractional_coordinates(v0+0.5_flyt)-0.5_flyt
                ! Convert to Cartesian
                v1=matmul(basis,v1)
                if ( lo_sqnorm(v1) .gt. sqrc ) cycle
                ! lattice vector
                v2=v1-v0
                v2=matmul(basis,anint(v2))
                l=l+1
                dumr(:,l,a1)=v1
                dumlv(:,l,a1)=v2
                dumd(l,a1)=norm2(v1)
                dumind(l,a1)=a2
            enddo
            dumctr(a1)=l
            ! Sort by distance
            call qsort(dumd(1:l,a1),sortind(1:l))
            dumr(:,1:l,a1)=dumr(:,sortind(1:l),a1)
            dumlv(:,1:l,a1)=dumlv(:,sortind(1:l),a1)
            dumind(1:l,a1)=dumind(sortind(1:l),a1)
        enddo

        ! Now make sure the cutoff makes sense, or if I have to increase it a little
        ! to get off of peaks of the rdf. Think it's ok to do this serially for now.
        buf=0.0_flyt
        do j=1,100
            l=0
            do a1=1,np
            do i=1,dumctr(a1)
                if ( abs( cutoff+buf-dumd(i,a1) ) .lt. 2*lo_tol ) then
                    l=l+1
                endif
            enddo
            enddo
            if ( l .eq. 0 ) then
                exit
            else
                buf=buf+lo_tol
            endif
        enddo
        ! now store everything
        dt%np=np
        dt%cutoff=0.0_flyt
        lo_allocate(dt%particle(dt%np))
        do a1=1,np
            l=0
            do i=1,dumctr(a1)
                if ( dumd(i,a1) .lt. cutoff+buf ) l=l+1
            enddo
            dt%particle(a1)%n=l
            lo_allocate(dt%particle(a1)%d(l))
            lo_allocate(dt%particle(a1)%v(3,l))
            lo_allocate(dt%particle(a1)%lv(3,l))
            lo_allocate(dt%particle(a1)%ind(l))
            dt%particle(a1)%d=lo_chop( dumd(1:l,a1), lo_sqtol )
            dt%particle(a1)%v=lo_chop( dumr(:,1:l,a1), lo_sqtol )
            dt%particle(a1)%lv=lo_chop( dumlv(:,1:l,a1), lo_sqtol )
            dt%particle(a1)%ind=dumind(1:l,a1)
            dt%cutoff=max(dt%cutoff,dt%particle(a1)%d(l))
        enddo
        ! And some cleanup
        lo_deallocate(dumd)
        lo_deallocate(dumr)
        lo_deallocate(dumlv)
        lo_deallocate(dumind)
        lo_deallocate(dumctr)
    end block dirnnser
    end select

    ! And some final cleanup
    lo_deallocate(cleanpositions)
    if ( verbosity .gt. 0 ) write(*,*) '... done generating distance table (',tochar(walltime()-timer),'s)'

end subroutine

!function distance_point_box(latticevector,point,basis) result(distance)
!    !> coordinates of origin of the box
!    real(flyt), dimension(3), intent(in) :: latticevector
!    !> the point
!    real(flyt), dimension(3), intent(in) :: point
!    !> actual box
!    real(flyt), dimension(3,3), intent(in) :: basis
!    !> distance to the point
!    real(flyt) :: distance
!
!    real(flyt), dimension(3) :: spt
!    real(flyt), dimension(3) :: v0,v1,v2,v3,n1,n2,n3,n4,n5,n6
!    real(flyt), dimension(3) :: r1,r2
!    real(flyt) :: d1,d2,d3,d4,d5,d6
!    integer :: i,j,k,l,ii,jj,kk
!
!    ! Shift the point by the origin of the box.
!    spt=point-latticevector
!
!    ! vectors that build up planes and stuff
!    v1=basis(:,1)
!    v2=basis(:,2)
!    v3=basis(:,3)
!    ! Get the normals of the planes
!    n1=lo_cross(v2,v3)
!    n2=lo_cross(v1,v3)
!    n3=lo_cross(v1,v2)
!    n1=n1/norm2(n1)
!    n2=n2/norm2(n2)
!    n3=n3/norm2(n3)
!    ! Make sure my normals point the way I think, i.e. normals pointing out
!    v0=matmul(basis,[0.5_flyt,0.5_flyt,0.5_flyt])
!    if ( dot_product(n1,v0) .gt. 0.0_flyt ) n1=-n1
!    if ( dot_product(n2,v0) .gt. 0.0_flyt ) n2=-n2
!    if ( dot_product(n3,v0) .gt. 0.0_flyt ) n3=-n3
!    ! And the three opposing planes
!    n4=-n1
!    n5=-n2
!    n6=-n3
!    ! Then the distances from the point to the planes
!    d1=dot_product(n1,spt)
!    d2=dot_product(n2,spt)
!    d3=dot_product(n3,spt)
!    d4=dot_product(n4,spt-v1)
!    d5=dot_product(n5,spt-v2)
!    d6=dot_product(n6,spt-v3)
!
!    ! Use the plane distances to work out which of the 27 quadrants I am in.
!    if ( d1 .gt. lo_sqtol ) then
!        ii=1 ! x < 0
!    elseif ( d4 .gt. lo_sqtol  ) then
!        ii=3 ! x < 1
!    else
!        ii=2 ! 0 < x < 1
!    endif
!    if ( d2 .gt. lo_sqtol ) then
!        jj=1 ! y < 0
!    elseif ( d5 .gt. lo_sqtol ) then
!        jj=3 ! y < 1
!    else
!        jj=2 ! 0 < y < 1
!    endif
!    if ( d3 .gt. lo_sqtol ) then
!        kk=1 ! z < 0
!    elseif ( d6 .gt. lo_sqtol ) then
!        kk=3 ! z < 1
!    else
!        kk=2 ! 0 < z < 1
!    endif
!
!    distance=-lo_huge
!    ! Now there are 27 cases. A little confusing, but ok.
!    if ( ii .eq. 2 .and. jj .eq. 2 .and. kk .eq. 2 ) then
!        ! Inside the box
!        distance=0.0_flyt
!    ! These are when the corners are the closest
!    elseif( ii .eq. 1 .and. jj .eq. 1 .and. kk .eq. 1 ) then
!        distance=norm2( spt )
!        distance=-lo_huge ! norm2( spt )
!    elseif( ii .eq. 1 .and. jj .eq. 1 .and. kk .eq. 3 ) then
!        distance=norm2( spt-v3 )
!write(*,*) tochar([d1,d2,d3,d4,d5,d6])
!        distance=-lo_huge ! norm2( spt )
!    elseif( ii .eq. 1 .and. jj .eq. 3 .and. kk .eq. 1 ) then
!        distance=norm2( spt-v2 )
!        distance=-lo_huge ! norm2( spt )
!    elseif( ii .eq. 3 .and. jj .eq. 1 .and. kk .eq. 1 ) then
!        distance=norm2( spt-v1 )
!        distance=-lo_huge ! norm2( spt )
!    elseif( ii .eq. 1 .and. jj .eq. 3 .and. kk .eq. 3 ) then
!        distance=norm2( spt-v2-v3 )
!        distance=-lo_huge ! norm2( spt )
!    elseif( ii .eq. 3 .and. jj .eq. 1 .and. kk .eq. 3 ) then
!        distance=norm2( spt-v1-v3 )
!        distance=-lo_huge ! norm2( spt )
!    elseif( ii .eq. 3 .and. jj .eq. 3 .and. kk .eq. 1 ) then
!        distance=norm2( spt-v1-v2 )
!        distance=-lo_huge ! norm2( spt )
!    elseif( ii .eq. 3 .and. jj .eq. 3 .and. kk .eq. 3 ) then
!        distance=norm2( spt-v1-v2-v3 )
!        distance=-lo_huge ! norm2( spt )
!
!    ! These are when an edge is the closest
!    elseif ( ii .eq. 2 .and. jj .eq. 1 .and. kk .eq. 1 ) then
!        r1=0.0_flyt
!        r2=v1
!        v0=r1+dot_product( spt-r1,r2-r1 )/dot_product(r2-r1,r2-r1)*(r2-r1)
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!    elseif ( ii .eq. 2 .and. jj .eq. 1 .and. kk .eq. 3 ) then
!        r1=v3
!        r2=v3+v1
!        v0=r1+dot_product( spt-r1,r2-r1 )/dot_product(r2-r1,r2-r1)*(r2-r1)
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!    elseif ( ii .eq. 2 .and. jj .eq. 3 .and. kk .eq. 1 ) then
!        r1=v2
!        r2=v2+v1
!        v0=r1+dot_product( spt-r1,r2-r1 )/dot_product(r2-r1,r2-r1)*(r2-r1)
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!    elseif ( ii .eq. 2 .and. jj .eq. 3 .and. kk .eq. 3 ) then
!        r1=v2+v3
!        r2=v2+v3+v1
!        v0=r1+dot_product( spt-r1,r2-r1 )/dot_product(r2-r1,r2-r1)*(r2-r1)
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!
!
!
!    elseif ( ii .eq. 1 .and. jj .eq. 2 .and. kk .eq. 1 ) then
!        r1=0.0_flyt
!        r2=v2
!        v0=r1+dot_product( spt-r1,v2 )/dot_product(v2,v2)*v2
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!    elseif ( ii .eq. 1 .and. jj .eq. 2 .and. kk .eq. 3 ) then
!        r1=v3
!        r2=v3+v2
!        v0=r1+dot_product( spt-r1,v2 )/dot_product(v2,v2)*v2
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!    elseif ( ii .eq. 3 .and. jj .eq. 2 .and. kk .eq. 1 ) then
!        r1=v1
!        r2=v1+v2
!        v0=r1+dot_product( spt-r1,v2 )/dot_product(v2,v2)*v2
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!    elseif ( ii .eq. 3 .and. jj .eq. 2 .and. kk .eq. 3 ) then
!        r1=v1+v3
!        r2=v1+v3+v2
!        v0=r1+dot_product( spt-r1,v2 )/dot_product(v2,v2)*v2
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!
!
!
!
!    elseif ( ii .eq. 1 .and. jj .eq. 1 .and. kk .eq. 2 ) then
!        r1=0.0_flyt
!        r2=v3
!        v0=r1+dot_product( spt-r1,v3 )/dot_product(v3,v3)*v3
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!    elseif ( ii .eq. 1 .and. jj .eq. 3 .and. kk .eq. 2 ) then
!        r1=v2
!        r2=v2+v3
!        v0=r1+dot_product( spt-r1,v3 )/dot_product(v3,v3)*v3
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!    elseif ( ii .eq. 3 .and. jj .eq. 1 .and. kk .eq. 2 ) then
!        r1=v1
!        r2=v1+v3
!        v0=r1+dot_product( spt-r1,v3 )/dot_product(v3,v3)*v3
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!    elseif ( ii .eq. 3 .and. jj .eq. 3 .and. kk .eq. 2 ) then
!        r1=v1+v2
!        r2=v1+v2+v3
!        v0=r1+dot_product( spt-r1,v3 )/dot_product(v3,v3)*v3
!        distance=norm2(spt-v0)
!        distance=-lo_huge ! norm2( spt )
!
!    ! And finally when a face is closest
!    elseif ( ii .eq. 1 .and. jj .eq. 2 .and. kk .eq. 2 ) then
!        distance=d1
!    elseif ( ii .eq. 3 .and. jj .eq. 2 .and. kk .eq. 2 ) then
!        distance=d4
!    elseif ( ii .eq. 2 .and. jj .eq. 1 .and. kk .eq. 2 ) then
!        distance=d2
!    elseif ( ii .eq. 2 .and. jj .eq. 3 .and. kk .eq. 2 ) then
!        distance=d5
!    elseif ( ii .eq. 2 .and. jj .eq. 2 .and. kk .eq. 1 ) then
!        distance=d3
!    elseif ( ii .eq. 2 .and. jj .eq. 2 .and. kk .eq. 3 ) then
!        distance=d6
!    endif
!
!end function

#ifdef AGRESSIVE_SANITY
function shortest_distance_cell_to_point(latticevector,point,basis) result(distance)
#else
pure function shortest_distance_cell_to_point(latticevector,point,basis) result(distance)
#endif
    real(flyt), dimension(3), intent(in) :: latticevector
    real(flyt), dimension(3), intent(in) :: point
    real(flyt), dimension(3,3), intent(in) :: basis
    real(flyt) :: distance

    real(flyt), dimension(3) :: v0
    real(flyt) :: r0

    ! Bounding sphere of box
    r0=lo_bounding_sphere_of_box(basis)
    ! Get the center of the box
    v0=latticevector+matmul(basis,[0.5_flyt,0.5_flyt,0.5_flyt])-point
    distance=norm2(v0)-r0
end function

!> Get the largest number of neighours
pure function max_number_of_neighbours(dt) result(n)
    !> distance table
    class(lo_distancetable), intent(in) :: dt
    !> number
    integer :: n

    integer :: i

    n=-1
    do i=1,dt%np
        n=max(dt%particle(i)%n,n)
    enddo
end function

!> Get the smallest number of neighours
pure function min_number_of_neighbours(dt) result(n)
    !> distance table
    class(lo_distancetable), intent(in) :: dt
    !> number
    integer :: n

    integer :: i

    n=huge(n)
    do i=1,dt%np
        n=min(dt%particle(i)%n,n)
    enddo
end function

!> remove unwanted neighbours from the distance table
subroutine prune(dt,cutoff_per_particle)
    !> distance table
    class(lo_distancetable), intent(inout) :: dt
    !> max cutoff, per particle
    real(flyt), dimension(dt%np), intent(in) :: cutoff_per_particle

    integer :: i,j,n
    real(flyt), dimension(:,:), allocatable :: v
    real(flyt), dimension(:,:), allocatable :: lv
    real(flyt), dimension(:), allocatable :: d
    real(flyt), dimension(:), allocatable :: weight
    integer, dimension(:), allocatable :: ind
    !
    i=dt%max_number_of_neighbours()
    lo_allocate(v(3,i))
    lo_allocate(lv(3,i))
    lo_allocate(d(i))
    lo_allocate(weight(i))
    lo_allocate(ind(i))

    do i=1,dt%np
        ! see how many I will keep. I assume that the distance table is sorted properly.
        n=-1
        do j=1,dt%particle(i)%n
            if ( dt%particle(i)%d(j)-lo_tol .gt. cutoff_per_particle(i) ) then
                n=j
                exit
            endif
        enddo
        ! maybe no pruning needed?
        if ( n .eq. dt%particle(i)%n ) cycle
        ! make a copy and reallocate
        if( allocated(dt%particle(i)%v) ) then
            v(:,1:n)=dt%particle(i)%v(:,1:n)
            deallocate(dt%particle(i)%v)
            allocate(dt%particle(i)%v,source=v(:,1:n))
        endif
        if( allocated(dt%particle(i)%lv) ) then
            lv(:,1:n)=dt%particle(i)%lv(:,1:n)
            deallocate(dt%particle(i)%lv)
            allocate(dt%particle(i)%lv,source=lv(:,1:n))
        endif
        if( allocated(dt%particle(i)%ind) ) then
            ind(1:n)=dt%particle(i)%ind(1:n)
            deallocate(dt%particle(i)%ind)
            allocate(dt%particle(i)%ind,source=ind(1:n))
        endif
        if( allocated(dt%particle(i)%weight) ) then
            weight(1:n)=dt%particle(i)%weight(1:n)
            deallocate(dt%particle(i)%weight)
            allocate(dt%particle(i)%weight,source=weight(1:n))
        endif
        ! and store the new number of neighbours
        dt%particle(i)%n=n
    enddo
end subroutine

!> Count the number of neighbours per atom within a certain cutoff. Not including itself.
subroutine count_neighbours_within_cutoff(dt,cutoff,ctr)
    !> distance table
    class(lo_distancetable), intent(in) :: dt
    !> cutoff
    real(flyt), intent(in) :: cutoff
    !> how many neighbours for each atom?
    integer, dimension(:), intent(out) :: ctr

    integer :: i,j,k

    ctr=0
    ! it might be really simple
    if ( cutoff .gt. dt%cutoff-lo_tol ) then
        do i=1,dt%np
            ctr(i)=dt%particle(i)%n-1
        enddo
        return
    endif
    do i=1,dt%np
        k=0
        do j=2,dt%particle(i)%n
            if ( dt%particle(i)%d(j) .lt. cutoff+lo_tol ) then
                k=k+1
            else
                exit
            endif
        enddo
        ctr(i)=k
    enddo
end subroutine

end module
