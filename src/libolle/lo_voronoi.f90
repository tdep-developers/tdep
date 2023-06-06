module lo_voronoi
!!
!! Generate Voronoi diagrams for periodic systems.
!!
use konstanter, only: r8,lo_iou,lo_status,lo_huge,lo_hugeint,lo_exitcode_symmetry,lo_exitcode_mpi
use gottochblandat, only: tochar,lo_sqnorm,lo_return_unique
use geometryfunctions, only: lo_plane,lo_bounding_sphere_of_box
use lo_sorting, only: lo_return_unique_indices
use type_distancetable, only: lo_distancetable,lo_distancetable_particle
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully,MPI_INTEGER,MPI_DOUBLE_PRECISION,MPI_CHARACTER

implicit none
private
public :: lo_voronoi_diagram
public :: lo_voronoi_cell

!> Face of a cell in the Voronoi diagram
type lo_voronoi_diagram_face
    !> how many points on the face
    integer :: n_node=-lo_hugeint
    !> which nodes make up the face?
    integer, dimension(:), allocatable :: ind
    !> what particle is on the other side of this face?
    integer :: neighbour=-lo_hugeint
    !> what is the vector that points to that neighbour?
    real(r8), dimension(3) :: neighbourvector=lo_huge
    !> what plane defines this face
    type(lo_plane) :: plane
end type

!> Cell in a voronoi diagram
type lo_voronoi_cell
    !> How many nodes?
    integer :: n_node=-lo_hugeint
    !> How many faces?
    integer :: n_face=-lo_hugeint
    !> Nodes, Cartesian coordinates, where zero is the point the cell is defined for.
    real(r8), dimension(:,:), allocatable :: node
    !> Faces of the polyhedron
    type(lo_voronoi_diagram_face), dimension(:), allocatable :: face
    !> smallest sphere inside cell
    real(r8) :: rmin=-lo_huge
    !> bounding sphere
    real(r8) :: rmax=-lo_huge
    !> volume of cell
    real(r8) :: volume=-lo_huge
    contains
        !> create a single cell in the voronoi diagram
        procedure :: generate=>generate_onecell
        !> check if a point is inside this polyhedron
        procedure :: is_point_inside
end type

!> Voroni diagram for a set of points with periodic boundary conditions.
type lo_voronoi_diagram
    !> How point in the diagram
    integer :: n_point=-lo_hugeint
    !> The cells
    type(lo_voronoi_cell), dimension(:), allocatable :: cell
    ! contains
    !     !> generate the diagram
    !     procedure, pass(voro) :: generate=>generate_diagram
end type

! Some parameters for the tesselation, safe defaults and such.
! Will likely have to be adjusted based on corner cases, but a
! problem for the future. Will increase incrementally.
real(r8), parameter :: safe_large_number=1E50_r8
integer, parameter :: max_n_node_per_face=50
integer, parameter :: max_n_face=150
integer, parameter :: max_n_node=120
integer, parameter :: safe_n_node_per_face=30

! Helper type to store the intermediate polygons
type dummypolygon
    !> How many points on this polygon?
    integer :: n=-lo_hugeint
    !> points
    real(r8), dimension(3,max_n_node_per_face) :: r
    !> centroid
    real(r8), dimension(3) :: centroid=-lo_huge
    !> Radius of bounding sphere
    real(r8) :: rmax=lo_huge
    !> has this polygon recently been modified?
    logical :: modified=.false.
    !> The plane the points are confined to
    type(lo_plane) :: plane
end type

contains

!> Generate the Voronoi cell for a single point
subroutine generate_onecell(cell,dtp,cutoff,tol,mem)
    !> Voronoi polyhedron
    class(lo_voronoi_cell), intent(out) :: cell
    !> Distance table for the specific particle
    type(lo_distancetable_particle), intent(in) :: dtp
    !> safe outer cutoff, should contain the cell with a wide margin
    real(r8), intent(in) :: cutoff
    !> tolerance
    real(r8), intent(in) :: tol
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ! This should be enough information to generate the diagram?
    type(dummypolygon), dimension(max_n_face) :: face
    type(lo_plane) :: plane
    real(r8), dimension(:,:), allocatable :: dr0
    real(r8), dimension(3) :: normal,midpoint,v0
    real(r8) :: rmax,f0
    integer, dimension(:), allocatable :: ind,redind
    integer :: i,j,k,l,n_face,n_changed,imax

    ! Initialize polyhedron to something safe
    do i=1,size(face)
        face(i)%n=0
        face(i)%modified=.false.
    enddo
    n_face=0
    n_changed=0
    ! Add a large cube as the starting point
    !@todo Should probably make it a huge tetrahedron instead, would be faster.
    call build_cube(n_face,face,cutoff)
    ! And a very large bounding sphere
    rmax=lo_huge

    ! Start slicing with all the planes defined by half distances to neighbours, not including
    ! itself (that's why I start from i=2). They are sorted by distance, so that once the distance
    ! is larger than the bounding sphere, I can stop!
    imax=0
    do i=2,dtp%n
        ! point half-way to this neighbour
        midpoint=dtp%v(:,i)*0.5_r8
        normal=midpoint
        ! check if we are done?
        imax=i
        if ( norm2(midpoint) .gt. rmax+tol ) exit
        ! if not, construct this plane
        call plane%generate(normal=normal,point=midpoint)
        ! slice with this plane, and in the process update the radius of the bounding sphere
        call slice_polyhedron(n_face,face,plane,rmax,tol,n_changed)
    enddo

    ! Ok, now I have created a neat polyhedron, always useful. But I want a Voronoi cell, so
    ! it makes more sense to rearrange it a little. First count number of nodes
    l=0
    do i=1,size(face)
        if ( face(i)%n .eq. 0 ) cycle
        l=l+face(i)%n
    enddo

    call mem%allocate(dr0,[3,l],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(redind,l,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dr0=0.0_r8
    redind=0

    ! Get the unique nodes
    l=0
    do i=1,size(face)
        if ( face(i)%n .eq. 0 ) cycle
        do j=1,face(i)%n
            l=l+1
            dr0(:,l)=face(i)%r(:,j)
        enddo
    enddo
    ! Get the unique nodes
    call lo_return_unique_indices(dr0,ind,mem,redind,tol)
    cell%n_node=size(ind)
    allocate(cell%node(3,cell%n_node))
    cell%node=0.0_r8
    do i=1,cell%n_node
        cell%node(:,i)=dr0(:,ind(i))
    enddo

    ! Good, now store the faces.
    cell%n_face=n_face
    allocate(cell%face(cell%n_face))
    l=0
    j=0
    do i=1,size(face)
        if ( face(i)%n .eq. 0 ) cycle
        j=j+1
        ! space for indices
        cell%face(j)%n_node=face(i)%n
        allocate(cell%face(j)%ind( cell%face(j)%n_node ))
        cell%face(j)%ind=0
        ! Store plane
        cell%face(j)%plane=face(i)%plane
        ! Store indices
        do k=1,cell%face(j)%n_node
            l=l+1
            cell%face(j)%ind(k)=redind(l)
        enddo
    enddo

    call mem%deallocate(dr0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(ind,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(redind,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! Get some more information that could come handy eventually.
    ! First radius of bounding sphere.
    f0=0.0_r8
    do i=1,cell%n_node
        f0=max(f0,lo_sqnorm(cell%node(:,i)))
    enddo
    cell%rmax=sqrt(f0)+tol
    ! Then the radius of the inscribed sphere
    f0=lo_huge
    do i=1,cell%n_face
        f0=min(f0,abs(cell%face(i)%plane%p))
    enddo
    cell%rmin=f0-tol

    ! Then the volume
    j=0
    do i=1,cell%n_face
        j=max(j,cell%face(i)%n_node)
    enddo
    call mem%allocate(dr0,[2,j+1],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dr0=0.0_r8


    cell%volume=0.0_r8
    do i=1,cell%n_face
        ! Get the centroid
        v0=0.0_r8
        do j=1,cell%face(i)%n_node
            k=cell%face(i)%ind(j)
            v0=v0+cell%node(:,k)/real(cell%face(i)%n_node,r8)
        enddo
        ! Make note that 2*centroid points to another atom?
        cell%face(i)%neighbourvector=2*v0
        ! slice into triangles
        dr0=0.0_r8
        do j=1,cell%face(i)%n_node
            k=cell%face(i)%ind(j)
            dr0(1,j)=dot_product(cell%face(i)%plane%v1,cell%node(:,k)-v0)
            dr0(2,j)=dot_product(cell%face(i)%plane%v2,cell%node(:,k)-v0)
        enddo
        dr0(:,cell%face(i)%n_node+1)=dr0(:,1)
        ! area of this face
        f0=0.0_r8
        do j=1,cell%face(i)%n_node
            f0=f0+dr0(1,j)*dr0(2,j+1)-dr0(1,j+1)*dr0(2,j)
        enddo
        f0=f0*0.5_r8
        ! accumulate volume
        cell%volume=cell%volume+abs(f0*cell%face(i)%plane%p)/3.0_r8
    enddo
    call mem%deallocate(dr0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    ! And we are done for now
end subroutine

! !> Generate a Voronoi diagram for a set of points with periodic boundary conditions. The algorithm is geometrical: given a table of distances (taking periodic boundary conditions into account), I slice a polyhedron with the planes that lie at the midpoint between two points. I stop when the midpoint is outside the bounding sphere of the polyhedron. This works great for atoms and other things that are more or less uniformly distributed, but will obviously fail for non-periodic cases.
! subroutine generate_diagram(voro,points,basis,verbosity,cutoff,mw)
!     !> diagram
!     class(lo_voronoi_diagram), intent(out) :: voro
!     !> points in fractional 0-1 coordinates
!     real(r8), dimension(:,:), intent(in) :: points
!     !> basis that transform the points to Cartesian coordinates
!     real(r8), dimension(3,3), intent(in) :: basis
!     !> how much to talk
!     integer, intent(in) :: verbosity
!     !> max cutoff distance
!     real(r8), intent(in), optional :: cutoff
!     !> MPI helper, maybe?
!     type(lo_mpi_helper), intent(inout), optional :: mw
!     !> weight per point, in case I don't want the normal diagram
!     !real(r8), dimension(:), intent(in), optional :: weight
!     !> use a different metric than the usual?
!     !integer, intent(in), optional :: metric
!
!     type(lo_distancetable) :: dt
!     real(r8) :: rc
!     real(r8) :: timer,t0,t1
!     integer, dimension(:), allocatable :: rnkctr,rnkpts
!     integer :: i,j,l,npts
!
!     ! First some simple heuristics, check dimensions and so on.
!     if ( verbosity .gt. 0 ) then
!         timer=mpi_wtime()
!         t0=timer
!         t1=0.0_r8
!         write(lo_iou,*)
!         write(lo_iou,*) 'Generating Voronoi diagram ( ',tochar(size(points,2)),' points)'
!     endif
!
!     ! how many points are there?
!     voro%n_point=size(points,2)
!     ! make space for the diagram
!     allocate(voro%cell(voro%n_point))
!     do i=1,voro%n_point
!         voro%cell(i)%n_node=-1
!         voro%cell(i)%n_face=-1
!     enddo
!
!     ! Decide on a cutoff, or guesstimate
!     if ( present(cutoff) ) then
!         rc=cutoff
!     else
!         ! this should be a safe upper bound for sure
!         rc=lo_bounding_sphere_of_box(basis)*2.0_r8
!     endif
!
!     ! Use a fix cutoff, or guess it.
!     if ( verbosity .gt. 0 ) then
!         write(lo_iou,*) '... building distance table with cutoff=',tochar(rc,4)
!     endif
!
!     ! Build the distance table
!     if ( present(mw) ) then
!         call dt%generate(points,basis,rc,verbosity=-1,mw=mw)
!     else
!         call dt%generate(points,basis,rc,verbosity=-1)
!     endif
!
!     if ( verbosity .gt. 0 ) then
!         t1=mpi_wtime()
!         write(lo_iou,*) '... got distance table in ',tochar(t1-t0),'s'
!         t0=t1
!     endif
!
!     ! Here I can insert a metric change, I will do that by updating the
!     ! distance variable in the type_distancetable and sort it again, then
!     ! the rest of the routine can stay the same.
!     if ( present(mw) ) then
!         ! When we do it in parallel, I have to think a little about distributing
!         ! the workload. Count number of points per rank.
!         allocate(rnkctr(mw%n))
!         rnkctr=0
!         do i=1,voro%n_point
!             j=mod(i,mw%n)+1
!             rnkctr(j)=rnkctr(j)+1
!         enddo
!         npts=rnkctr(mw%r+1)
!         if ( npts .gt. 0 ) then
!             allocate(rnkpts(npts))
!             rnkpts=0
!         endif
!         l=0
!         do i=1,mw%n
!         do j=1,rnkctr(i)
!             l=l+1
!             if ( mw%r+1 .eq. i ) then
!                 rnkpts(j)=l
!             endif
!         enddo
!         enddo
!
!         do i=1,npts
!             j=rnkpts(i) ! Make sure we pick the right point
!             call polyhedron_for_one_point(voro%cell(j),dt,j,tol=1E-6_r8)
!         enddo
!         ! Spread it out everywhere.
!         call communicate_diagram(voro,mw)
!     else
!         ! Solve it serially, much less to bother about.
!         do i=1,voro%n_point
!             call polyhedron_for_one_point(voro%cell(i),dt,i,tol=1E-6_r8)
!         enddo
!     endif
!
!     if ( verbosity .gt. 0 ) then
!         write(lo_iou,*) '... got tesselation in ',tochar(mpi_wtime()-timer),'s'
!     endif
! end subroutine

! !> Generate the Voronoi cell for a single point
! subroutine polyhedron_for_one_point(cell,dt,point,tol)
!     !> the Voronoi polyhedron
!     type(lo_voronoi_diagram_cell), intent(out) :: cell
!     !> the Distance table
!     type(lo_distancetable), intent(in) :: dt
!     !> which point in the distance table to calculate the diagram for
!     integer, intent(in) :: point
!     !> tolerance
!     real(r8), intent(in) :: tol
!
!     ! This should be enough information to generate the diagram?
!     type(dummypolygon), dimension(max_n_face) :: face
!     type(lo_plane) :: plane
!     real(r8), dimension(:,:), allocatable :: dr0
!     real(r8), dimension(3) :: normal,midpoint,v0
!     real(r8) :: rmax,f0
!     integer :: i,j,k,l,ii,n_face,n_changed,imax
!
!     ! Initialize polyhedron to something safe
!     do i=1,size(face)
!         face(i)%n=0
!         face(i)%modified=.false.
!     enddo
!     ! Add a large cube as the starting point
!     call build_cube(n_face,face,3*dt%cutoff)
!     ! And a very large bounding sphere
!     rmax=lo_huge
!
!     ! Start slicing with all the planes defined by half distances to neighbours, not including
!     ! itself (that's why I start from i=2). They are sorted by distance, so that once the distance
!     ! is larger than the bounding sphere, I can stop!
!     imax=0
!     do i=2,dt%particle(point)%n
!         ! point half-way to this neighbour
!         midpoint=dt%particle(point)%v(:,i)*0.5_r8
!         normal=midpoint
!         ! check if we are done?
!         imax=i
!         if ( norm2(midpoint) .gt. rmax+tol ) exit
!         ! if not, construct this plane
!         call plane%generate(normal=normal,point=midpoint)
!         ! slice with this plane, and in the process update the radius of the bounding sphere
!         call slice_polyhedron(n_face,face,plane,rmax,tol,n_changed)
!     enddo
!
!     ! Ok, now I have created a neat polyhedron, always useful. But I want a Voronoi cell, so
!     ! it makes more sense to rearrange it a little. First count number of nodes
!     l=0
!     do i=1,size(face)
!         if ( face(i)%n .eq. 0 ) cycle
!         l=l+face(i)%n
!     enddo
!     allocate(dr0(3,l))
!     ! Get the unique nodes
!     l=0
!     do i=1,size(face)
!         if ( face(i)%n .eq. 0 ) cycle
!         do j=1,face(i)%n
!             l=l+1
!             dr0(:,l)=face(i)%r(:,j)
!         enddo
!     enddo
!     ! Now we have the unique nodes!
!     call lo_return_unique(dr0,cell%node)
!     cell%n_node=size(cell%node,2)
!
!     ! Good, now store the faces.
!     ! Space for the faces?
!     cell%n_face=n_face
!     allocate(cell%face(cell%n_face))
!     l=0
!     do i=1,size(face)
!         if ( face(i)%n .eq. 0 ) cycle
!         l=l+1
!         cell%face(l)%n_node=face(i)%n
!         allocate(cell%face(l)%ind( cell%face(l)%n_node ))
!         ! Match the points with the nodes
!         do j=1,face(i)%n
!             ii=0
!             v0=face(i)%r(:,j)
!             do k=1,cell%n_node
!                 if ( lo_sqnorm(cell%node(:,k)-v0) .lt. tol**2 ) then
!                     ii=k
!                     exit
!                 endif
!             enddo
!             if ( ii .gt. 0 ) then
!                 cell%face(l)%ind(j)=ii
!             else
!                 call lo_stop_gracefully(['Could not locate point on face.'],lo_exitcode_symmetry)
!             endif
!         enddo
!
!         ! Now match the face with a neighbour. Not sure how to do it fast, but I
!         ! don't think this is a bottleneck anyway.
!         ii=0
!         do k=2,imax
!             if ( abs(dt%particle(point)%d(k)*0.5_r8+face(i)%plane%p) .lt. tol ) then
!                 ! same distance from origin
!                 v0=dt%particle(point)%v(:,k)
!                 v0=v0/norm2(v0)
!                 if ( lo_sqnorm(v0-face(i)%plane%normal) .lt. tol**2 ) then
!                     ! Same normal, it's a match!
!                     ii=k
!                     exit
!                 endif
!             endif
!         enddo
!         if ( ii .gt. 0 ) then
!             cell%face(l)%plane=face(i)%plane
!             cell%face(l)%neighbourvector=dt%particle(point)%v(:,ii)
!             cell%face(l)%neighbour=dt%particle(point)%ind(ii)
!         else
!             call lo_stop_gracefully(['Could not locate neighbour for face.'],lo_exitcode_symmetry)
!         endif
!     enddo
!
!     ! Get some more information that could come handy eventually.
!     ! First radius of bounding sphere.
!     f0=0.0_r8
!     do i=1,cell%n_node
!         f0=max(f0,lo_sqnorm(cell%node(:,i)))
!     enddo
!     cell%rmax=sqrt(f0)+tol
!     ! Then the radius of the inscribed sphere
!     f0=lo_huge
!     do i=1,cell%n_face
!         f0=min(f0,abs(cell%face(i)%plane%p))
!     enddo
!     cell%rmin=f0-tol
!     ! And we are done for now
! end subroutine

!> slice a polyhedron with a plane
subroutine slice_polyhedron(n_face,face,plane,rmax,tol,n_changed)
    !> how many faces do we currently have
    integer, intent(inout) :: n_face
    !> list of faces
    type(dummypolygon), dimension(:), intent(inout) :: face
    !> plane to slice with
    type(lo_plane), intent(in) :: plane
    !> radius of bounding sphere
    real(r8), intent(inout) :: rmax
    !> not the default tolerance?
    real(r8), intent(in) :: tol
    !> how many faces got changed?
    integer, intent(out) :: n_changed

    real(r8), dimension(3,max_n_node_per_face*2) :: df1,df2
    real(r8), dimension(3) :: rj,ri,v0
    real(r8) :: f0,f1
    integer :: i_face,i_node
    integer :: ctr_new_face,ctr_curr_face,ctr_sliced
    integer :: prev_cond,curr_cond
    integer :: i,j,l
    logical :: sliced

    n_changed=0
    ctr_new_face=0
    ctr_sliced=0
    df1=safe_large_number
    df2=safe_large_number

    ! So, this is a little tricky to understand. On entry, we have a convex polyhedron
    ! defined as a list of faces, i.e. convex polygons. Below takes a plane and slices
    ! each face. I keep what is below the plane, and throw away what is above.
    faceloop: do i_face=1,size(face)
        ! if the face has been dropped, skip it.
        if ( face(i_face)%n .eq. 0 ) cycle

        ! Check distance to centroid
        f0=plane%distance_to_point(face(i_face)%centroid)

        if ( f0 .lt. -(face(i_face)%rmax+10*tol) ) then
            ! way below the plane, we don't have to care or adjust
            cycle faceloop
        endif

        if ( f0 .gt. (face(i_face)%rmax+10*tol) ) then
            ! way above, drop this face
            face(i_face)%n=0
            cycle faceloop
        endif

        ! Ok, have to slice this face. That means I will slice each line segment.
        ! Each line segment is from r(:,i_node-1) to r(:,i_node).
        sliced=.false.
        ctr_curr_face=0
        rj=face(i_face)%r(:,face(i_face)%n)
        prev_cond=planecondition(plane,rj,tol)
        do i_node=1,face(i_face)%n
            ri=face(i_face)%r(:,i_node)
            curr_cond=planecondition(plane,ri,tol)
            if ( prev_cond*curr_cond .eq. -1 ) then
                ! this linesegment is sliced by the plane, we get a new point.
                ! figure out where this point is.
                v0=quick_line_intersection(ri,rj,plane)

                ! store in the appropriate places, first in the updated current polygon
                ctr_curr_face=ctr_curr_face+1
                df1(:,ctr_curr_face)=v0
                ! and store in the counter for the new face
                ctr_new_face=ctr_new_face+1
                df2(:,ctr_new_face)=v0
                ! make to note that we have sliced
                sliced=.true.
            endif

            ! Make note to keep the points that are on or below the plane for the current face
            if ( curr_cond .le. 0 ) then
                ctr_curr_face=ctr_curr_face+1
                df1(:,ctr_curr_face)=ri
            endif

            ! If the node is on the plane, it should be added to the list of new points.
            if ( curr_cond .eq. 0 ) then
                ctr_new_face=ctr_new_face+1
                df2(:,ctr_new_face)=ri
            endif
            ! Now update what the previous point is
            prev_cond=curr_cond
            rj=ri
        enddo

        ! Now some thing could have happened. Figure out what.
        if ( ctr_curr_face .le. 2 ) then
            ! this face got dropped. It happens. I drop it when
            ! it has two or less points, since then it is irrelevant.
            face(i_face)%n=0
            cycle faceloop
        endif

        ! The current face could be different, if so, update it!
        if ( sliced .or. ctr_curr_face .ne. face(i_face)%n ) then
            ! Make note of the number or slices
            ctr_sliced=ctr_sliced+1
            ! Yep, current face has changed. Make sure to only
            ! store the unique points.
            face(i_face)%r=safe_large_number
            l=0
            il0: do i=1,ctr_curr_face
                do j=1,l
                    v0=df1(:,i)-face(i_face)%r(:,j)
                    if ( v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3) .lt. tol**2 ) then
                        cycle il0
                    endif
                enddo
                l=l+1
                face(i_face)%r(:,l)=df1(:,i)
            enddo il0
            ! If we have two or fewer points left, drop the face.
            if ( l .le. 2 ) then
                face(i_face)%n=0
                cycle faceloop
            else
                face(i_face)%n=l
                face(i_face)%modified=.true.
            endif
        endif
    enddo faceloop

    ! Now, do we have something new? Collect the unique new points.
    df1=safe_large_number
    l=0
    il1: do i=1,ctr_new_face
        ! Sanity test, remove later
        ! if ( abs(plane%distance_to_point(df2(:,i))) .gt. tol ) then
        !     call lo_stop_gracefully(['New points not on plane'],lo_exitcode_symmetry,__FILE__,__LINE__)
        ! endif

        do j=1,l
            v0=df2(:,i)-df1(:,j)
            if ( v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3) .lt. tol**2 ) then
                cycle il1
            endif
        enddo
        l=l+1
        df1(:,l)=df2(:,i)
    enddo il1

    ! This test is fast, can keep.
    if ( ctr_sliced .gt. 0 .and. l .le. 2 ) then
        call lo_stop_gracefully(['Should have points.'],lo_exitcode_symmetry,__FILE__,__LINE__)
    endif

    ! Two or fewer points is not a polygon, so it might as well be zero.
    if ( l .le. 2 ) then
        l=0
    endif

    ! In case there are some new points, add those as a new face.
    if ( l .gt. 0 ) then
        ! Store the new face in the first empty spot
        i_face=0
        do i=1,size(face)
            if ( face(i)%n .eq. 0 ) then
                i_face=i
                exit
            endif
        enddo
        if ( i_face .eq. 0 ) then
            call lo_stop_gracefully(['Could not find empty spot for new face'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        face(i_face)%n=l
        face(i_face)%r(:,1:l)=df1(:,1:l)
        face(i_face)%plane=plane
        face(i_face)%modified=.true.
    else
        ! This means there was nothing new, and nothing happened.
        n_changed=0
    endif

    ! Now go through the polyhedron and check if I have to update
    ! some metadata on any faces.
    n_face=0
    do i=1,size(face)
        if ( face(i)%n .eq. 0 ) cycle
        n_face=n_face+1
        if ( face(i)%modified ) then
            n_changed=n_changed+1
            ! Seems this face needs to be taken care of.
            ! first make sure the points are sorted clockwise.
            call face(i)%plane%anglesort(face(i)%r(:,1:face(i)%n))
            ! Get the centroid
            v0=0.0_r8
            do j=1,face(i)%n
                v0=v0+face(i)%r(:,j)
            enddo
            face(i)%centroid=v0/real(face(i)%n,r8)
            ! and the bounding sphere
            f1=0.0_r8
            do j=1,face(i)%n
                v0=face(i)%r(:,j)-face(i)%centroid
                f0=v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3)
                f1=max(f0,f1)
            enddo
            face(i)%rmax=sqrt(f1)
            face(i)%modified=.false.
        endif
    enddo

    ! While I'm at it, calculate the radius of the bounding sphere of the polyhedron.
    if ( n_changed .gt. 0 ) then
        rmax=0.0_r8
        do i=1,size(face)
        do j=1,face(i)%n
            v0=face(i)%r(:,j)
            f0=v0(1)*v0(1)+v0(2)*v0(2)+v0(3)*v0(3)
            rmax=max(rmax,f0)
        enddo
        enddo
        rmax=sqrt(rmax)
    endif
end subroutine

!> linesegment plane intersection with all safety-checks removed
pure function quick_line_intersection(r1,r2,plane) result(point)
    real(r8), dimension(3), intent(in) :: r1,r2
    type(lo_plane), intent(in) :: plane
    real(r8), dimension(3) :: point
    !
    real(r8), dimension(3) :: dv
    real(r8) :: f0,f1
    dv=r2-r1
    f0=dot_product(dv,plane%normal)
    f1=-(plane%p+dot_product(r1,plane%normal))
    point=r1+(f1/f0)*dv
end function

!> tells where a point is with respect to a plane, 0 on the plane, -1 below and 1 above.
pure integer function planecondition(plane,point,tol)
    type(lo_plane), intent(in) :: plane
    real(r8), dimension(3), intent(in) :: point
    real(r8), intent(in) :: tol
    !
    real(r8) :: f0
    f0=plane%distance_to_point(point)
    if ( abs(f0) .lt. tol ) then
        planecondition=0
    elseif ( f0 .gt. 0.0_r8 ) then
        planecondition=1
    else
        planecondition=-1
    endif
end function

! !> communicate the partial diagram to everyone
! subroutine communicate_diagram(voro,mw)
!     !> voronoi diagram
!     type(lo_voronoi_diagram), intent(inout) :: voro
!     !> MPI helper
!     type(lo_mpi_helper), intent(inout) :: mw
!
!     integer :: i,j
!     integer, dimension(:), allocatable :: rnksize,rnkoff
!     integer :: bytes,totbytes,pos
!     character, dimension(:), allocatable :: sbuf,rbuf
!
!     ! Count the number of bytes on each rank
!     bytes=0
!     do i=1,voro%n_point
!         if ( voro%cell(i)%n_node .lt. 0 ) cycle
!         bytes=bytes+storage_size(i)
!         bytes=bytes+storage_size(voro%cell(i)%n_node)
!         bytes=bytes+storage_size(voro%cell(i)%n_face)
!         bytes=bytes+storage_size(voro%cell(i)%node)*size(voro%cell(i)%node)
!         bytes=bytes+storage_size(voro%cell(i)%rmin)
!         bytes=bytes+storage_size(voro%cell(i)%rmax)
!         do j=1,voro%cell(i)%n_face
!             bytes=bytes+storage_size(voro%cell(i)%face(j)%n_node)
!             bytes=bytes+storage_size(voro%cell(i)%face(j)%ind)*size(voro%cell(i)%face(j)%ind)
!             bytes=bytes+storage_size(voro%cell(i)%face(j)%neighbour)
!             bytes=bytes+storage_size(voro%cell(i)%face(j)%neighbourvector)*size(voro%cell(i)%face(j)%neighbourvector)
!             bytes=bytes+storage_size(voro%cell(i)%face(j)%plane%normal)*size(voro%cell(i)%face(j)%plane%normal)
!             bytes=bytes+storage_size(voro%cell(i)%face(j)%plane%v1)*size(voro%cell(i)%face(j)%plane%v1)
!             bytes=bytes+storage_size(voro%cell(i)%face(j)%plane%v2)*size(voro%cell(i)%face(j)%plane%v2)
!             bytes=bytes+storage_size(voro%cell(i)%face(j)%plane%p)
!         enddo
!     enddo
!     bytes=bytes/8
!
!     ! Figure out the size per rank and the offsets
!     allocate(rnksize(mw%n))
!     allocate(rnkoff(mw%n))
!     rnkoff=0
!     rnksize=0
!     rnksize(mw%r+1)=bytes
!     call mw%allreduce('sum',rnksize)
!     totbytes=sum(rnksize)
!     j=0
!     do i=1,mw%n
!         rnkoff(i)=j
!         j=j+rnksize(i)
!     enddo
!     ! Make space for the recieve buffer
!     allocate(rbuf(sum(rnksize)))
!     ! space for the send buffer
!     if ( bytes .gt. 0 ) then
!         allocate(sbuf(bytes))
!     else
!         ! dummy allocation in case rank is empty
!         allocate(sbuf(1))
!     endif
!     ! Now, pack the send buffer
!     pos=0
!     lo_status=0
!     do i=1,voro%n_point
!         if ( voro%cell(i)%n_node .lt. 0 ) cycle
!         call mpi_pack(i,1,MPI_INTEGER,sbuf,bytes,pos,mw%comm,mw%error)
!         call mpi_pack(voro%cell(i)%n_node,1,MPI_INTEGER,sbuf,bytes,pos,mw%comm,mw%error)
!         call mpi_pack(voro%cell(i)%n_face,1,MPI_INTEGER,sbuf,bytes,pos,mw%comm,mw%error)
!         call mpi_pack(voro%cell(i)%node,size(voro%cell(i)%node),MPI_DOUBLE_PRECISION,sbuf,bytes,pos,mw%comm,mw%error)
!         call mpi_pack(voro%cell(i)%rmin,1,MPI_DOUBLE_PRECISION,sbuf,bytes,pos,mw%comm,mw%error)
!         call mpi_pack(voro%cell(i)%rmax,1,MPI_DOUBLE_PRECISION,sbuf,bytes,pos,mw%comm,mw%error)
!         do j=1,voro%cell(i)%n_face
!             call mpi_pack(voro%cell(i)%face(j)%n_node,1,MPI_INTEGER,sbuf,bytes,pos,mw%comm,mw%error)
!             call mpi_pack(voro%cell(i)%face(j)%ind,size(voro%cell(i)%face(j)%ind),MPI_INTEGER,sbuf,bytes,pos,mw%comm,mw%error)
!             call mpi_pack(voro%cell(i)%face(j)%neighbour,1,MPI_INTEGER,sbuf,bytes,pos,mw%comm,mw%error)
!             call mpi_pack(voro%cell(i)%face(j)%neighbourvector,3,MPI_DOUBLE_PRECISION,sbuf,bytes,pos,mw%comm,mw%error)
!             call mpi_pack(voro%cell(i)%face(j)%plane%normal,3,MPI_DOUBLE_PRECISION,sbuf,bytes,pos,mw%comm,mw%error)
!             call mpi_pack(voro%cell(i)%face(j)%plane%v1,3,MPI_DOUBLE_PRECISION,sbuf,bytes,pos,mw%comm,mw%error)
!             call mpi_pack(voro%cell(i)%face(j)%plane%v2,3,MPI_DOUBLE_PRECISION,sbuf,bytes,pos,mw%comm,mw%error)
!             call mpi_pack(voro%cell(i)%face(j)%plane%p,1,MPI_DOUBLE_PRECISION,sbuf,bytes,pos,mw%comm,mw%error)
!         enddo
!     enddo
!
!     ! Double-check that I have packed correctly
!     if ( pos .ne. bytes ) then
!         call lo_stop_gracefully(['I could not pack properly.'],lo_exitcode_mpi,mw%comm)
!     endif
!
!     ! Spread it around
!     call mpi_allgatherv(sbuf,bytes,MPI_CHARACTER,rbuf,rnksize,rnkoff,MPI_CHARACTER,mw%comm,mw%error)
!     if ( mw%error .ne. 0 ) then
!         call lo_stop_gracefully(['mpi_allgatherv exit code: '//tochar(mw%error)],lo_exitcode_mpi,mw%comm)
!     endif
!     ! We do a little cleanup now!
!     deallocate(sbuf)
!     deallocate(rnksize)
!     deallocate(rnkoff)
!
!     ! Now, unpack everything
!     pos=0
!     do i=1,voro%n_point
!         call mpi_unpack(rbuf,totbytes,pos,j,1,MPI_INTEGER,mw%comm,mw%error)
!         if ( i .ne. j ) then
!             call lo_stop_gracefully(['I could not communicate.'],lo_exitcode_mpi,mw%comm)
!         endif
!         call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%n_node,1,MPI_INTEGER,mw%comm,mw%error)
!         call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%n_face,1,MPI_INTEGER,mw%comm,mw%error)
!
!         if ( .not.allocated(voro%cell(i)%node) ) then
!             allocate(voro%cell(i)%node( 3,voro%cell(i)%n_node ))
!         endif
!         if ( .not.allocated(voro%cell(i)%face) ) then
!             allocate(voro%cell(i)%face( voro%cell(i)%n_face ))
!         endif
!         call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%node,voro%cell(i)%n_node*3,MPI_DOUBLE_PRECISION,mw%comm,mw%error)
!         call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%rmin,1,MPI_DOUBLE_PRECISION,mw%comm,mw%error)
!         call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%rmax,1,MPI_DOUBLE_PRECISION,mw%comm,mw%error)
!         do j=1,voro%cell(i)%n_face
!             call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%face(j)%n_node,1,MPI_INTEGER,mw%comm,mw%error)
!             if ( .not.allocated(voro%cell(i)%face(j)%ind) ) then
!                 allocate(voro%cell(i)%face(j)%ind( voro%cell(i)%face(j)%n_node ))
!             endif
!             call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%face(j)%ind,voro%cell(i)%face(j)%n_node,MPI_INTEGER,mw%comm,mw%error)
!             call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%face(j)%neighbour,1,MPI_INTEGER,mw%comm,mw%error)
!             call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%face(j)%neighbourvector,3,MPI_DOUBLE_PRECISION,mw%comm,mw%error)
!             call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%face(j)%plane%normal,3,MPI_DOUBLE_PRECISION,mw%comm,mw%error)
!             call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%face(j)%plane%v1,3,MPI_DOUBLE_PRECISION,mw%comm,mw%error)
!             call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%face(j)%plane%v2,3,MPI_DOUBLE_PRECISION,mw%comm,mw%error)
!             call mpi_unpack(rbuf,totbytes,pos,voro%cell(i)%face(j)%plane%p,1,MPI_DOUBLE_PRECISION,mw%comm,mw%error)
!         enddo
!     enddo
!     ! And clean the buffer
!     deallocate(rbuf)
! end subroutine

!> initialize the polyhedron to a cube
subroutine build_cube(n_face,face,side)
    !> polyhedron
    type(dummypolygon), dimension(max_n_face), intent(inout) :: face
    !> number of faces
    integer, intent(out) :: n_face
    !> Side of cube
    real(r8), intent(in) :: side

    real(r8), dimension(3,8) :: corners
    real(r8), dimension(3,6) :: normals
    real(r8), dimension(3) :: v0
    integer, dimension(4,6) :: faces
    integer :: i,j,k,l

    ! a cube has 6 faces
    n_face=6
    ! Get the corners
    l=0
    do i=0,1
    do j=0,1
    do k=0,1
        l=l+1
        corners(:,l)=((/i,j,k/)*1.0_r8-0.5_r8)*side
    enddo
    enddo
    enddo

    ! I know from this what the faces are
    faces(:,1)=[1,2,3,4]
    faces(:,2)=[1,2,5,6]
    faces(:,3)=[1,3,5,7]
    faces(:,4)=[5,6,7,8]
    faces(:,5)=[3,4,7,8]
    faces(:,6)=[2,4,6,8]
    normals(:,1)=[-1.0_r8, 0.0_r8, 0.0_r8]
    normals(:,2)=[ 0.0_r8,-1.0_r8, 0.0_r8]
    normals(:,3)=[ 0.0_r8, 0.0_r8,-1.0_r8]
    normals(:,4)=[ 1.0_r8, 0.0_r8, 0.0_r8]
    normals(:,5)=[ 0.0_r8, 1.0_r8, 0.0_r8]
    normals(:,6)=[ 0.0_r8, 0.0_r8, 1.0_r8]

    do i=1,n_face
        ! store the points
        face(i)%n=4
        !lo_allocate(polyhedron%face(i)%r(3,4))
        face(i)%r(:,1:4)=corners(:,faces(:,i))
        ! get the plane that the planes are on
        call face(i)%plane%generate(normal=normals(:,i),point=corners(:,faces(1,i)))
        call face(i)%plane%anglesort(face(i)%r(:,1:4))
        ! the centroid
        v0=0.0_r8
        do j=1,4
            v0=v0+face(i)%r(:,j)*0.25_r8
        enddo
        face(i)%centroid=v0
        ! the radius of bounding sphere
        face(i)%rmax=0
        do j=1,4
            face(i)%rmax=max(face(i)%rmax,norm2(face(i)%centroid-face(i)%r(:,j)))
        enddo
    enddo
    ! And it's done
end subroutine

!> test if a point is inside the voronoi polyhedron
function is_point_inside(cell,point,tol) result(inside)
    !> voronoi cell
    class(lo_voronoi_cell), intent(in) :: cell
    !> point, in Cartesian coordinates
    real(r8), dimension(3), intent(in) :: point
    !> tolerance
    real(r8), intent(in) :: tol
    !> is the point inside?
    logical :: inside

    integer :: i

    inside=.true.
    do i=1,cell%n_face
        ! if it's outside any face, it's false
        if ( cell%face(i)%plane%distance_to_point(point) .gt. tol ) then
            inside=.false.
            return
        endif
    enddo
end function

end module
