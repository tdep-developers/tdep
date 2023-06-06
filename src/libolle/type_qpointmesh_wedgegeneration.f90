submodule (type_qpointmesh) type_qpointmesh_wedgegeneration
!!
!! Generate wedge-based meshes
!!
use cgal_wrappers
use lo_sorting, only: lo_return_unique
use lo_verletboxes, only: lo_verletbox
implicit none
contains

#ifdef usecgal
!> generate fancy new kind of mesh
module subroutine lo_get_wedge_mesh(qp,uc,griddensity,mw,mem,verbosity)
    !> q-mesh
    type(lo_wedge_mesh), intent(inout) :: qp
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: uc
    !> normal grid dimensions
    integer, dimension(3), intent(in) :: griddensity
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot
    integer, intent(in) :: verbosity

    real(r8), dimension(7) :: meshpars
    real(r8) :: timer
    integer :: solrnk

    ! Start timer
    timer=walltime()

    call mem%tick()

    ! Decide on some things first, like how to generate the mesh.
    init: block
        real(r8) :: f0,f1,ideal_side
        integer :: i

        ! First a simple sanity check
        if ( uc%info%havewedge .eqv. .false. ) then
            call lo_stop_gracefully(['Really need the wedge here'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
        endif

        ! Report what we are doing
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) '... building wedge mesh'
        endif

        ! Decide which rank does the solving
        solrnk=mw%n-1

        ! some weird guess for the number of points:
        i=4*product(griddensity)/uc%sym%n
        ! The volume of the wedge
        f0=abs(lo_determ(uc%reciprocal_latticevectors))/uc%sym%n
        ! Volume per tetrahedron, sort of
        f1=f0/i
        ! Ideal side of a tetrahedron with that volume
        ideal_side=(f1*6.0_r8*sqrt(2.0_r8))**(1.0_r8/3.0_r8)

        ! Here I convert these measures to a parameters that CGAL understands
        ! In the end, I just use the ideal side of each tetrahedron as the
        ! relevant measure. Seems to be the most robust.
        meshpars(1)=20_r8       ! max_dihedral_angle, for feature detection
        meshpars(2)=ideal_side  ! edge_size
        meshpars(3)=20_r8       ! facet_angle
        meshpars(4)=ideal_side  ! facet_size
        meshpars(5)=ideal_side  ! facet_distance
        meshpars(6)=2.0_r8      ! cell_radius_edge_ratio
        meshpars(7)=ideal_side  ! cell_size

        !meshpars(1)=20           ! 120  max_dihedral_angle
        !meshpars(2)=f0           ! 0.025  edge_size
        !meshpars(3)=20.0_r8      ! 20     facet_angle
        !meshpars(4)=f0           ! 0.03   facet_size
        !meshpars(5)=10000.0_r8   !  0.05   facet_distance
        !meshpars(6)=3            ! 3      cell_radius_edge_ratio
        !meshpars(7)=f0           ! 0.05   cell_size
    end block init

    ! Do the CGAL things and tesselate the irreducible wedge
    getrawwedge: block
        real(r8), dimension(:,:), allocatable :: wedgevertices,dr0,dr1
        real(r8), dimension(3) :: v0,v1
        real(r8) :: f0 !,f1,f2
        integer, dimension(:,:), allocatable :: wedgetets
        integer :: iface,inode,ctr,ii,jj

        ! Get the proper points for the convex hull. CGAL is far too exact, and gets really confused
        ! if there are collinear points, even with double precision tolerance it gets angry. So I
        ! filter things out with a sensible tolerance to give CGAL a clean polyhedron to work with.
        ctr=0
        do iface=1,uc%irrw%nfaces
            ctr=ctr+uc%irrw%face(iface)%n
        enddo
        call mem%allocate(dr0,[3,ctr],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        dr0=0.0_r8
        ! Remove collinear points
        ctr=0
        do iface=1,uc%irrw%nfaces
            do inode=1,uc%irrw%face(iface)%n
                ! previous point
                ii=lo_index_in_periodic_array(inode-1,uc%irrw%face(iface)%n)
                ! next point
                jj=lo_index_in_periodic_array(inode+1,uc%irrw%face(iface)%n)
                ! Does it form an angle that is not 180?
                v0=uc%irrw%r( :,uc%irrw%face(iface)%ind(ii) )-uc%irrw%r( :,uc%irrw%face(iface)%ind(inode) )
                v1=uc%irrw%r( :,uc%irrw%face(iface)%ind(jj) )-uc%irrw%r( :,uc%irrw%face(iface)%ind(inode) )
                f0=lo_angle_between_vectors(v0,v1)
                if ( abs(f0-lo_pi) .gt. lo_radiantol ) then
                    ctr=ctr+1
                    dr0(:,ctr)=uc%irrw%r(:,uc%irrw%face(iface)%ind(inode))
                endif
            enddo
        enddo
        ! Now keep only the unique from this
        call lo_return_unique(dr0,dr1,mem,lo_tol)

        if ( verbosity .gt. 0 ) write(lo_iou,*) '... built convex hull of wedge, got ',tochar(size(dr1,2)),' points'

        ! Not sure if CGAL is safe to run across MPI, just do it once to be on the safe side
        if ( verbosity .gt. 0 ) write(lo_iou,*) '... calling CGAL to get raw tesselation'

        if ( mw%r .eq. solrnk ) then
            ! This calls CGAL
            call lo_tesselate_polyhedron(dr1,meshpars,wedgevertices,wedgetets)
            ! Throw away the tetrahedrons, will redo that later.
            deallocate(wedgetets)
            ! Make a note of the number of points generated.
            ctr=size(wedgevertices,2)
        else
            ctr=0
        endif

        ! Distribute the points, first the number of points
        call mw%bcast(ctr,from=solrnk)
        ! Then make space to recieve it
        if ( mw%r .ne. solrnk ) then
            allocate(wedgevertices(3,ctr))
            wedgevertices=0.0_r8
        endif
        ! And distribute the actual points
        call mw%bcast(wedgevertices,from=solrnk)

        if ( verbosity .gt. 0 ) write(lo_iou,*) '... obtained raw set of points (',tochar(walltime()-timer),'s)'

        call mem%deallocate(dr0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(dr1,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)

        ! Tesselate and sort out the symmetry
        call qp%tesselate_wedge_mesh(uc,wedgevertices,-1.0_r8,-1.0_r8,mw,mem,verbosity)

        ! Update the integration weights
        !call qp%update_integration_weight(uc,mw,mem,verbosity)

        ! And we might be done! Just a little cleanup left to do.
        deallocate(wedgevertices)
    end block getrawwedge

    call mem%tock(__FILE__,__LINE__,mw%comm)
end subroutine

!> Given a bunch of points in the irreducible wedge, create the tesselation
module subroutine tesselate_wedge_mesh(qp,uc,wedgevertices,prunetol,splittol,mw,mem,verbosity)
    !> q-mesh
    class(lo_wedge_mesh), intent(inout) :: qp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> vertices that are supposed to be turned into q-points
    real(r8), dimension(:,:), allocatable, intent(inout) :: wedgevertices
    !> perhaps remove some points if they are really close together.
    real(r8), intent(in) :: prunetol
    !> perhaps split some tetrahedrons if they are very large with respect to the mean
    real(r8), intent(in) :: splittol
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8) :: t0,t1
    integer, dimension(:,:), allocatable :: wedgetets
    logical :: flipstuff
    integer :: nsym

    ! Start timers
    t0=walltime()
    t1=t0

    ! Safety check for memory consumption
    call mem%tick()

    ! Set some basic stuff
    init: block
        type(lo_verletbox) :: vb
        real(r8), dimension(:,:), allocatable :: dr
        logical, dimension(:), allocatable :: fixedpoint
        integer, dimension(:), allocatable :: ind
        integer, dimension(3) :: boxdim
        integer :: i,j,k,l,npts,bi,bj,bk,solrnk

        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) '... tesselating wedge mesh'
        endif

        ! What rank to solve on?
        solrnk=mw%n-1

        ! Make sure everything is nothing
        qp%n_irr_point=0
        qp%n_full_point=0
        qp%n_irr_tet=0
        qp%n_full_tet=0
        qp%scaledrecbasis=0.0_r8
        !qp%timereversal=uc%sym%timereversal
        if ( allocated(qp%ip) ) deallocate(qp%ip)
        if ( allocated(qp%ap) ) deallocate(qp%ap)
        if ( allocated(qp%it) ) deallocate(qp%it)
        if ( allocated(qp%at) ) deallocate(qp%at)
        if ( allocated(qp%operationok) ) deallocate(qp%operationok)

        ! Now, perhaps remove some points? If not, just keep the list as-is
        if ( prunetol .gt. 0.0_r8 ) then
            ! Make a little workspace
            npts=size(wedgevertices,2)
            call mem%allocate(dr,[3,npts],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(fixedpoint,npts,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%allocate(ind,npts,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            ind=1
            dr=wedgevertices
            fixedpoint=.false.
            ! Choose starting size for Verlet boxes
            boxdim=vb%boxdim(dr,6000,1E-8_r8)

            ! First make sure that I don't remove any point that is on the edge of the wedge
            ptloop1: do i=1,npts
                ! Make it parallel
                if ( mod(i,mw%n) .ne. mw%r ) cycle
                !do j=1,uc%irrw%nfaces
                !    if ( abs(uc%irrw%face(j)%plane%distance_to_point(dr(:,i))) .lt. lo_sqtol ) then
                !        fixedpoint(i)=.true.
                !        cycle ptloop1
                !    endif
                !enddo
                do j=1,uc%irrw%nnodes
                    if ( lo_sqnorm(uc%irrw%r(:,j)-dr(:,i)) .lt. lo_sqtol ) then
                        fixedpoint(i)=.true.
                        cycle ptloop1
                    endif
                enddo
            enddo ptloop1
            ! Sync the fixed points
            call mw%allreduce('or',fixedpoint)

            ! Jam them into boxes
            call vb%generate(dr,boxdim,mem)
            do i=1,npts
                if ( mod(i,mw%n) .ne. mw%r ) cycle
                call vb%boxind(dr(:,i),bi,bj,bk)
                do j=1,vb%box(bi,bj,bk)%n
                    k=vb%box(bi,bj,bk)%ind(j)
                    if ( k .le. i ) cycle
                    if ( fixedpoint(k) ) cycle
                    if ( lo_sqnorm(dr(:,i)-dr(:,k)) .lt. prunetol**2 ) ind(k)=0
                enddo
            enddo
            call vb%destroy(mem)

            ! And again, with a different division
            boxdim=boxdim+1
            call vb%generate(dr,boxdim,mem)
            do i=1,npts
                if ( mod(i,mw%n) .ne. mw%r ) cycle
                call vb%boxind(dr(:,i),bi,bj,bk)
                do j=1,vb%box(bi,bj,bk)%n
                    k=vb%box(bi,bj,bk)%ind(j)
                    if ( k .le. i ) cycle
                    if ( fixedpoint(k) ) cycle
                    if ( lo_sqnorm(dr(:,i)-dr(:,k)) .lt. prunetol**2 ) ind(k)=0
                enddo
            enddo
            call vb%destroy(mem)
            ! Sync the diagnostics
            call mw%allreduce('min',ind)
            ! If the number changed, update the list of points
            if ( npts .ne. sum(ind) ) then
                deallocate(wedgevertices)
                allocate(wedgevertices(3,sum(ind)))
                l=0
                do i=1,npts
                    if ( ind(i) .eq. 0 ) cycle
                    l=l+1
                    wedgevertices(:,l)=dr(:,i)
                enddo
                if ( verbosity .gt. 0 ) then
                    t1=walltime()
                    write(*,*) '... removed ',tochar(npts-sum(ind)),' points that were too close (',tochar(t1-t0),'s)'
                    t0=t1
                endif
            endif
            call mem%deallocate(dr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(fixedpoint,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            call mem%deallocate(ind,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        endif

        ! Do the proper triangulation
        if ( verbosity .gt. 0 ) write(lo_iou,*) '... triangulating'

        ! Call CGAL on one rank for the delaunay triangulation.
        if ( mw%r .eq. solrnk ) then
            call lo_deltri_3d(wedgevertices,wedgetets)
            j=size(wedgetets,2)
        else
            j=0
        endif

        ! Distribute the tetrahedrons everywhere
        call mw%bcast(j,from=solrnk)
        if ( mw%r .ne. solrnk ) then
            allocate(wedgetets(4,j))
            wedgetets=0
        endif
        call mw%bcast(wedgevertices,from=solrnk)
        call mw%bcast(wedgetets,from=solrnk)
        ! Perhaps a sanity test
        if ( verbosity .gt. 0 ) then
            j=size(wedgevertices,2)
            k=size(wedgetets,2)
            t1=walltime()
            write(lo_iou,*) '... got ',tochar(j),' points and ',tochar(k),' tetrahedrons (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block init

    ! The irreducible wedge is more or less ok now, store that.
    storewedge: block
        integer :: i

        qp%n_irr_point=size(wedgevertices,2)
        allocate(qp%ip(qp%n_irr_point))
        do i=1,qp%n_irr_point
            qp%ip(i)%r=wedgevertices(:,i)
            call lo_get_small_group_of_qpoint(qp%ip(i),uc)
            ! some things are not fixed yet
            qp%ip(i)%integration_weight=-1.0_r8
            qp%ip(i)%n_full_point=-1
            qp%ip(i)%full_index=-1
        enddo

        ! And the irreducible tetrahedrons
        qp%n_irr_tet=size(wedgetets,2)
        allocate(qp%it(qp%n_irr_tet))
        do i=1,qp%n_irr_tet
            qp%it(i)%irreducible_index=wedgetets(:,i)
            ! and some things are not fixed yet
            qp%it(i)%full_index=-1
            qp%it(i)%integration_weight=-1.0_r8
        enddo

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... stored initial irreducible wedge (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block storewedge

    ! Use symmetry operations to get all points
    allpoints: block
        type(lo_verletbox) :: vb
        real(r8), dimension(:,:), allocatable :: dr0
        real(r8), dimension(3) :: v0
        real(r8) :: f0
        integer, dimension(:), allocatable :: di1,di2,di3
        integer, dimension(3) :: boxdim
        integer :: i,j,l,o,bi,bj,bk,bl,ii,jj,kk,ll,np

        ! Stupid check, see the fraction of the BZ the wedge occupies
        ! to decide what to do about timereversal symmetry
        f0=0.0_r8
        do i=1,size(wedgetets,2)
            f0=f0+lo_unsigned_tetrahedron_volume(wedgevertices(:,wedgetets(:,i)))
        enddo
        f0=1.0_r8/(uc%volume*f0)
        if ( abs(f0-uc%sym%n) .lt. 0.01_r8 ) then
            if ( verbosity .gt. 0 ) write(lo_iou,*) '... no need to consider timereversal'
            flipstuff=.false.
        elseif ( abs(f0-uc%sym%n*2) .lt. 0.01_r8 .and. qp%timereversal ) then
            if ( verbosity .gt. 0 ) write(lo_iou,*) '... adding time-reversal things'
            flipstuff=.true.
        else
            write(*,*) 'FIXME PROPER RETURN CODE'
            write(*,*) 'sum(wedge)/vol',f0,1.0_r8/f0
            write(*,*) 'too complicated relation stuff with timereversal symmetry. I will fix this at some point. Maybe.'
            stop
        endif

        ! Double the number operations, if I should add the flippy stuff
        if ( flipstuff ) then
            nsym=uc%sym%n*2
        else
            nsym=uc%sym%n
        endif
        ! Raw number of points
        np=nsym*qp%n_irr_point

        ! Temporary space for symmetry reduction
        call mem%allocate(dr0,[3,np],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(di1,np,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(di2,np,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(di3,np,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        dr0=0.0_r8
        di1=0
        di2=0
        di3=0
        ! Build all the points.
        l=0
        do o=1,uc%sym%n
        do i=1,qp%n_irr_point
            l=l+1
            if ( mod(l,mw%n) .ne. mw%r ) cycle
            ! rotate the vector
            v0=lo_operate_on_vector(uc%sym%op(o),qp%ip(i)%r,reciprocal=.true.)
            dr0(:,l)=v0 ! store point
            di1(l)=i   ! store irreducible index
            di2(l)=o   ! store operation
        enddo
        enddo
        if ( flipstuff ) then
            do o=1,uc%sym%n
            do i=1,qp%n_irr_point
                l=l+1
                if ( mod(l,mw%n) .ne. mw%r ) cycle
                ! rotate the vector, and reverse
                v0=-lo_operate_on_vector(uc%sym%op(o),qp%ip(i)%r,reciprocal=.true.)
                dr0(:,l)=v0 ! store point
                di1(l)=i   ! store irreducible index
                di2(l)=-o  ! store operation
            enddo
            enddo
        endif

        call mw%allreduce('sum',dr0)
        call mw%allreduce('sum',di1)
        call mw%allreduce('sum',di2)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... constructed raw full set of ',tochar(np),' points (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! stuff particles in boxes
        boxdim=vb%boxdim(dr0,24**3,1E-8_r8)
        call vb%generate(dr0,boxdim,mem)

        ! Now I should remove redundant points
        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        di3=1
        do i=1,np
            ! Make it a little parallel
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            ! it might already be redundant.
            if ( di3(i) .eq. 0 ) cycle
            ! In which box does this point live?
            call vb%boxind(dr0(:,i),bi,bj,bk)
            ! Check in this, and in neighbouring boxes if we have redundant points.
            do ii=max(1,bi-1),min(vb%nx,bi+1)
            do jj=max(1,bj-1),min(vb%ny,bj+1)
            do kk=max(1,bk-1),min(vb%nz,bk+1)
            do bl=1,vb%box(ii,jj,kk)%n
                ll=vb%box(ii,jj,kk)%ind(bl)
                if ( ll .le. i ) then
                    cycle
                else
                    if ( lo_sqnorm(dr0(:,i)-dr0(:,ll)) .lt. lo_sqtol .and. ll .gt. i ) then
                        di3(ll)=0
                    endif
                endif
            enddo
            enddo
            enddo
            enddo
            ! Report?
            if ( verbosity .gt. 0 .and. lo_trueNtimes(i,25,np) ) then
                call lo_progressbar(' ... removing redundant points',i,np,walltime()-t0)
            endif
        enddo

        ! Now we don't need the boxes anymore
        call vb%destroy(mem)

        ! And sync it up
        call mw%allreduce('min',di3)
        if ( verbosity .gt. 0 ) call lo_progressbar(' ... removing redundant points',np,np,walltime()-t0)
        t0=walltime()

        ! Now store all the points.
        qp%n_full_point=sum(di3)
        allocate(qp%ap(qp%n_full_point))
        i=0
        do l=1,np
            if ( di3(l) .eq. 0 ) cycle
            i=i+1
            ! store the point
            qp%ap(i)%r=dr0(:,l)
            call lo_get_small_group_of_qpoint(qp%ap(i),uc)
            qp%ap(i)%irreducible_index=di1(l)
            qp%ap(i)%operation_from_irreducible=di2(l)
            ! And some things we fix later
            qp%ap(i)%integration_weight=-1.0_r8
        enddo

        ! Sanity-test this, to be on the safe side
        do i=1,qp%n_full_point
            j=qp%ap(i)%irreducible_index
            o=qp%ap(i)%operation_from_irreducible
            if ( o .gt. 0 ) then
                v0=lo_operate_on_vector(uc%sym%op(o),qp%ap(i)%r,reciprocal=.true.,inverse=.true.)-qp%ip(j)%r
            else
                v0=-lo_operate_on_vector(uc%sym%op(abs(o)),qp%ap(i)%r,reciprocal=.true.,inverse=.true.)-qp%ip(j)%r
            endif
            if ( lo_sqnorm(v0) .gt. lo_sqtol ) then
                write(*,*) 'bad operation'
                write(*,*) i,j,o,lo_sqnorm(v0)
                stop
            endif
        enddo

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... got set of full points (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Intermediate cleanup
        call mem%deallocate(dr0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(di1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(di2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(di3,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) write(*,*) '... reduced from ',tochar(np),' to ',tochar(qp%n_full_point),' points ( ~',tochar(real(qp%n_full_point,r8)**(1.0_r8/3.0_r8),1),'^3  )'
    end block allpoints

    ! Figure out the indices on the grid for all the points.
    gridmap: block
        type(lo_verletbox) :: vb
        real(r8), dimension(:,:), allocatable :: r
        real(r8), dimension(3) :: v0,v1
        integer, dimension(:,:), allocatable :: dj,dk
        integer, dimension(:), allocatable :: di
        integer, dimension(3) :: boxdim
        integer :: i,j,k,l,o,ii

        ! Get all the points, and Verlet a Verlet box.
        call mem%allocate(r,[3,qp%n_full_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        r=0.0_r8
        do i=1,qp%n_full_point
            r(:,i)=qp%ap(i)%r
        enddo
        boxdim=vb%boxdim(r,30**3,1E-8_r8)
        call vb%generate(r,boxdim,mem)

        ! Ok, now all the points are boxed up. Start by locating the irreducible on the grid
        call mem%allocate(di,qp%n_irr_point,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        di=0
        do i=1,qp%n_irr_point
            ! make it parallel
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            ! Find this point
            ii=vb%locate(r,qp%ip(i)%r)
            ! If all went well, bp should be the grid index!
            if ( ii .gt. 0 ) then
                di(i)=ii
            else
                ! it failed!!
                call lo_stop_gracefully(['Could not locate irreducible point on the full grid.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
        enddo
        call mw%allreduce('sum',di)
        do i=1,qp%n_irr_point
            qp%ip(i)%full_index=di(i)
        enddo
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... located irreducible in full mesh (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Now build all the tetrahedrons!
        qp%n_full_tet=qp%n_irr_tet*nsym
        call mem%allocate(dj,[4,qp%n_full_tet],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(dk,[4,qp%n_full_tet],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        dj=0
        dk=0
        l=0

        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        do o=1,uc%sym%n
        do i=1,qp%n_irr_tet
            l=l+1
            ! Make it parallel
            if ( mod(l,mw%n) .ne. mw%r ) cycle

            do j=1,4
                ! irreducible index
                k=qp%it(i)%irreducible_index(j)
                dj(j,l)=k

                v0=qp%ip(k)%r   ! irreducible vector
                v1=lo_operate_on_vector(uc%sym%op(o),v0,reciprocal=.true.) ! rotated vector
                ii=vb%locate(r,v1)
                ! If all went well, ii should be the grid index!
                if ( ii .gt. 0 ) then
                    dk(j,l)=ii
                else
                    ! it failed!!
                    call lo_stop_gracefully(['Could not locate rotated point.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                endif
            enddo
            if ( verbosity .gt. 0 .and. lo_trueNtimes(l,87,qp%n_full_tet) ) then
                call lo_progressbar(' ... mapping grid indices',l,qp%n_full_tet,walltime()-t0)
            endif
        enddo
        enddo
        if ( flipstuff ) then
            do o=1,uc%sym%n
            do i=1,qp%n_irr_tet
                l=l+1
                ! Make it parallel
                if ( mod(l,mw%n) .ne. mw%r ) cycle

                do j=1,4
                    ! irreducible index
                    k=qp%it(i)%irreducible_index(j)
                    dj(j,l)=k
                    v0=qp%ip(k)%r   ! irreducible vector
                    v1=-lo_operate_on_vector(uc%sym%op(o),v0,reciprocal=.true.) ! rotated vector
                    ii=vb%locate(r,v1)
                    ! If all went well, ii should be the grid index!
                    if ( ii .gt. 0 ) then
                        dk(j,l)=ii
                    else
                        ! it failed!!
                        call lo_stop_gracefully(['Could not locate rotated point.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                    endif
                enddo
                if ( verbosity .gt. 0 .and. lo_trueNtimes(l,87,qp%n_full_tet) ) then
                    call lo_progressbar(' ... mapping grid indices',l,qp%n_full_tet,walltime()-t0)
                endif
            enddo
            enddo
        endif

        ! Sync
        call mw%allreduce('sum',dj)
        call mw%allreduce('sum',dk)
        ! Store
        allocate(qp%at(qp%n_full_tet))
        do i=1,qp%n_full_tet
            qp%at(i)%irreducible_index=dj(:,i)
            qp%at(i)%full_index=dk(:,i)
        enddo
        call mem%deallocate(dj,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(dk,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            call lo_progressbar(' ... mapping grid indices',qp%n_full_tet,qp%n_full_tet,t1-t0)
            t0=t1
        endif

        ! Also need the grid-indices for the irreducible tetrahedrons
        call mem%allocate(dj,[4,qp%n_irr_tet],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        dj=0
        do i=1,qp%n_irr_tet
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            do j=1,4
                ! irreducible index
                k=qp%it(i)%irreducible_index(j)
                ii=vb%locate(r,qp%ip(k)%r)
                ! If all went well, ii should be the grid index!
                if ( ii .gt. 0 ) then
                    dj(j,i)=ii
                else
                    call lo_stop_gracefully(['Could not locate rotated point.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                endif
            enddo
        enddo
        call mw%allreduce('sum',dj)
        do i=1,qp%n_irr_tet
            qp%it(i)%full_index=dj(:,i)
        enddo
        call mem%deallocate(dj,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        ! Now I can cleanup some more
        call mem%deallocate(r,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call vb%destroy(mem)

        ! I also want a list per irreducible q-point. This list should contain which points on the full grid
        ! it transforms to, and with which operations. First reset the counters.
        do i=1,qp%n_irr_point
            qp%ip(i)%n_full_point=0
        enddo
        ! Count the number of transformations
        do i=1,qp%n_full_point
            ii=qp%ap(i)%irreducible_index
            qp%ip(ii)%n_full_point=qp%ip(ii)%n_full_point+1
        enddo
        ! Space for operations and indices
        do i=1,qp%n_irr_point
            j=qp%ip(i)%n_full_point
            allocate(qp%ip(i)%operation_full_point(j))
            allocate(qp%ip(i)%index_full_point(j))
            qp%ip(i)%operation_full_point=0
            qp%ip(i)%index_full_point=0
            qp%ip(i)%n_full_point=0
        enddo
        ! Store them
        do i=1,qp%n_full_point
            ii=qp%ap(i)%irreducible_index
            qp%ip(ii)%n_full_point=qp%ip(ii)%n_full_point+1
            j=qp%ip(ii)%n_full_point
            qp%ip(ii)%index_full_point(j)=i
            qp%ip(ii)%operation_full_point(j)=qp%ap(i)%operation_from_irreducible
        enddo

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... stored all symmetry operations (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block gridmap

    ! Set all weights properly
    weights: block
        real(r8), dimension(3,4) :: dumtet
        real(r8), dimension(3) :: v0
        real(r8) :: f0,f1
        integer :: i,j,ii

        ! Set them all to zero
        do i=1,qp%n_irr_point
            qp%ip(i)%integration_weight=0.0_r8
        enddo
        do i=1,qp%n_full_point
            qp%ap(i)%integration_weight=0.0_r8
        enddo
        do i=1,qp%n_irr_tet
            qp%it(i)%integration_weight=0.0_r8
        enddo
        do i=1,qp%n_full_tet
            qp%at(i)%integration_weight=0.0_r8
        enddo

        ! Start with the irreducible
        f0=0.0_r8
        do i=1,qp%n_irr_tet
            dumtet=0.0_r8
            do j=1,4
                dumtet(:,j)=qp%ip( qp%it(i)%irreducible_index(j) )%r
            enddo
            qp%it(i)%integration_weight=lo_unsigned_tetrahedron_volume( dumtet )*nsym*uc%volume
            f0=f0+qp%it(i)%integration_weight
        enddo
        do i=1,qp%n_irr_tet
            qp%it(i)%integration_weight=qp%it(i)%integration_weight/f0
        enddo

        ! And all the tetrahedrons
        f0=0.0_r8
        do i=1,qp%n_full_tet
            dumtet=0.0_r8
            do j=1,4
                dumtet(:,j)=qp%ap( qp%at(i)%full_index(j) )%r
            enddo
            qp%at(i)%integration_weight=lo_unsigned_tetrahedron_volume( dumtet )*uc%volume
            f0=f0+qp%at(i)%integration_weight
        enddo
        do i=1,qp%n_full_tet
            qp%at(i)%integration_weight=qp%at(i)%integration_weight/f0
        enddo

        ! And flat weights for the points
        do i=1,qp%n_irr_tet
            f0=qp%it(i)%integration_weight*0.25_r8
            do j=1,4
                ii=qp%it(i)%irreducible_index(j)
                qp%ip(ii)%integration_weight=qp%ip(ii)%integration_weight+f0
            enddo
        enddo
        ! And flat weights for the points
        do i=1,qp%n_full_tet
            f0=qp%at(i)%integration_weight*0.25_r8
            do j=1,4
                ii=qp%at(i)%full_index(j)
                qp%ap(ii)%integration_weight=qp%ap(ii)%integration_weight+f0
                f1=f1+f0
            enddo
        enddo

        ! Get the scaled basis for adaptive gaussian:
        f0=(1.0_r8*qp%n_full_point)**(1.0_r8/3.0_r8) ! points per distance, sort of
        v0(1)=norm2(uc%reciprocal_latticevectors(:,1))
        v0(2)=norm2(uc%reciprocal_latticevectors(:,2))
        v0(3)=norm2(uc%reciprocal_latticevectors(:,3))
        v0=v0*f0/sum(v0) ! get it per axis or somthing
        ! and to normal units
        qp%scaledrecbasis=uc%reciprocal_latticevectors*lo_twopi
        qp%scaledrecbasis(:,1)=qp%scaledrecbasis(:,1)/v0(1)
        qp%scaledrecbasis(:,2)=qp%scaledrecbasis(:,2)/v0(2)
        qp%scaledrecbasis(:,3)=qp%scaledrecbasis(:,3)/v0(3)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... set all integration weights (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block weights

    ! And at the end, check that I have not forgot to allocate
    call mem%tock(__FILE__,__LINE__,mw%comm)
end subroutine

#else
    ! In case CGAL is not linked it, provide some stubs
    module subroutine lo_get_wedge_mesh(qp,uc,griddensity,mw,mem,verbosity)
        !> q-mesh
        type(lo_wedge_mesh), intent(inout) :: qp
        !> crystalstructure
        type(lo_crystalstructure), intent(in) :: uc
        !> normal grid dimensions
        integer, dimension(3), intent(in) :: griddensity
        !> MPI helper
        type(lo_mpi_helper), intent(inout) :: mw
        !> memory tracker
        type(lo_mem_helper), intent(inout) :: mem
        !> talk a lot
        integer, intent(in) :: verbosity

        integer :: i
        ! just stop
        !qp%nq_tot=0
        i=qp%n_irr_point
        i=uc%na
        i=griddensity(1)
        i=mw%n
        i=int(mem%n_allocate)
        i=verbosity
        if ( mw%talk ) write(*,*) 'This needs CGAL, install and recompile to enable'
        call mw%destroy()
        stop
    end subroutine
    module subroutine tesselate_wedge_mesh(qp,uc,wedgevertices,prunetol,splittol,mw,mem,verbosity)
        !> q-mesh
        class(lo_wedge_mesh), intent(inout) :: qp
        !> crystal structure
        type(lo_crystalstructure), intent(in) :: uc
        !> vertices that are supposed to be turned into q-points
        real(r8), dimension(:,:), allocatable, intent(inout) :: wedgevertices
        !> perhaps remove some points if they are really close together.
        real(r8), intent(in) :: prunetol
        !> maybe add some points
        real(r8), intent(in) :: splittol
        !> MPI helper
        type(lo_mpi_helper), intent(inout) :: mw
        !> memory tracker
        type(lo_mem_helper), intent(inout) :: mem
        !> talk a lot
        integer, intent(in) :: verbosity

        integer :: i
        i=size(wedgevertices,1)
        i=uc%na
        i=int(anint(prunetol))
        i=int(anint(splittol))
        !i=qp%nq_tot
        i=verbosity
        i=int(mem%n_allocate)
        i=qp%n_irr_point

        ! just stop
        if ( mw%talk ) write(*,*) 'This needs CGAL, install and recompile to enable'
        call mw%destroy()
        stop
    end subroutine
#endif
end submodule
