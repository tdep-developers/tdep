#include "precompilerdefinitions"
module type_voronoi
!! Generate Voronoi diagrams for periodic systems.
use konstanter, only: flyt,lo_huge,lo_tol,lo_sqtol,lo_status,lo_hugeint,lo_exitcode_param
use gottochblandat, only: tochar,walltime,lo_chop,lo_sqnorm,lo_index_in_periodic_array,lo_mean,lo_progressbar,&
                   lo_progressbar_init,lo_return_unique,qsort,lo_stop_gracefully
use geometryfunctions, only: lo_geometryobject,lo_plane,lo_bounding_sphere_of_box,lo_polygon,lo_polyhedron,lo_linesegment
use type_distancetable, only: lo_distancetable
use type_linkedlist, only: lo_linked_list

implicit none
private
public :: lo_voronoi_diagram
public :: lo_voronoi_diagram_cell

!> Node in the diagram
type lo_voronoi_diagram_node
    !> in how many faces share this node?
    integer :: nsharedfaces=-lo_hugeint
    !> in which faces is it a part of?
    integer, dimension(:), allocatable :: faceind
    !> how many cells share this node
    integer :: nsharedcells=-lo_hugeint
    !> which particles share this node, except for myself
    integer, dimension(:), allocatable :: neighbourind
    !> what's the vector that points to the connected cell?
    real(flyt), dimension(:,:), allocatable :: neighbourvector
#ifdef AGRESSIVE_SANITY
    contains
        procedure :: size_in_mem=>node_size_in_mem
#endif
end type

!> Edge in the diagram
type lo_voronoi_diagram_edge
    !> first point on edge
    integer :: i1=-lo_hugeint
    !> second point on edge
    integer :: i2=-lo_hugeint
    !> how many cells share this edge
    integer :: nshared=-lo_hugeint
end type

!> Face in the diagram
type lo_voronoi_diagram_face
    !> how many points on the face
    integer :: n=-lo_hugeint
    !> which nodes make up the face?
    integer, dimension(:), allocatable :: ind
    !> what particle is on the other side of this face?
    integer :: neighbour=-lo_hugeint
    !> what is the vector that points to that neighbour?
    real(flyt), dimension(3) :: neighbourvector=lo_huge
    !> what plane defines this face
    type(lo_plane) :: plane
#ifdef AGRESSIVE_SANITY
    contains
        procedure :: size_in_mem=>face_size_in_mem
#endif
end type

!> Cell in a voronoi diagram
type lo_voronoi_diagram_cell
    !> The polyhedron passed by the generating routines
    type(lo_polyhedron) :: polyhedron
    !> how many nodes?
    integer :: nnodes=-lo_hugeint
    !> how many faces
    integer :: nfaces=-lo_hugeint
    !> how many edges
    integer :: nedges=-lo_hugeint
    !> the actual faces
    type(lo_voronoi_diagram_face), dimension(:), allocatable :: face
    !> the nodes
    type(lo_voronoi_diagram_node), dimension(:), allocatable :: node
    !> the edges
    type(lo_voronoi_diagram_edge), dimension(:), allocatable :: edge
    !> the coordinates for the nodes
    real(flyt), dimension(:,:), allocatable :: r
    !> smallest sphere inside cell
    real(flyt) :: rmin=lo_huge
    !> bounding sphere
    real(flyt) :: rmax=lo_huge
#ifdef AGRESSIVE_SANITY
    contains
        procedure :: size_in_mem=>cell_size_in_mem
#endif
end type

!> Voroni diagram for a set of points with periodic boundary conditions.
type lo_voronoi_diagram
    !> How many particles
    integer :: n=-lo_hugeint
    !> The cells
    type(lo_voronoi_diagram_cell), dimension(:), allocatable :: cell
    contains
        !> generate the diagram
        procedure :: generate
#ifdef AGRESSIVE_SANITY
        procedure :: size_in_mem=>diagram_size_in_mem
#endif
end type

contains

!> Generate a Voronoi diagram for a set of points with periodic boundary conditions. The algorithm is geometrical: given a table of distances (taking periodic boundary conditions into account), I slice a polyhedron with the planes that lie at the midpoint between two points. I stop when the midpoint is outside the bounding sphere of the polyhedron. This works great for atoms and other things that are more or less uniformly distributed, but will obviously fail for non-periodic cases.
subroutine generate(voro,points,basis,cutoff,verbosity,dt,minjump)
    !> diagram
    class(lo_voronoi_diagram), intent(out) :: voro
    !> points in fractional 0-1 coordinates
    real(flyt), dimension(:,:), intent(in) :: points
    !> basis that transform the points to Cartesian coordinates
    real(flyt), dimension(3,3), intent(in) :: basis
    !> max cutoff distance
    real(flyt), intent(in), optional :: cutoff
    !> how much to talk
    integer, intent(in), optional :: verbosity
    !> send in a distance table instead of generating it
    type(lo_distancetable), intent(inout), optional :: dt
    !> smallest jump distance
    integer, intent(in), optional :: minjump
    !
    type(lo_distancetable) :: ddt
    real(flyt) :: rc,t0
    integer :: i,verb,neighbourthreshold

    ! How many neighbours should minimally be represented
    neighbourthreshold=60
    ! How much should I talk?
    if ( present(verbosity) ) then
        verb=verbosity
    else
        verb=0
    endif

    ! how many points are there?
    voro%n=size(points,2)
    lo_allocate(voro%cell(voro%n))
    if ( verb .gt. 0 ) then
        write(*,*) '... generating Voronoi diagram for ',tochar(voro%n),' points'
        t0=walltime()
    endif

    ! If no distance table is provided, calculate one
    if ( present(dt) .eqv. .false. ) then
        ! Use a fix cutoff, or guess it.
        if ( present(cutoff) ) then
            rc=cutoff
            if ( verb .gt. 1 ) write(*,*) '... building distance table with cutoff=',tochar(rc,4)
            call ddt%generate(points,basis,rc,verbosity=verb-1)
        else
            ! Don't ask where I got this from
            rc=lo_bounding_sphere_of_box(basis)*2.0_flyt
            if ( verb .gt. 1 ) write(*,*) '... building distance table with cutoff=',tochar(rc)
            call ddt%generate(points,basis,rc,verbosity=verb-1)
        endif

        if ( verb .gt. 0 ) write(*,*) '... got distance table in ',tochar(walltime()-t0),'s'
    else
        if ( verb .gt. 0 ) write(*,*) '... using provided distancetable'
    endif

    ! Build the polyhedrons. Parallel if it's over 100 points or something
    ! Actually no. Some compilers get really angry when it's parallel.
    ! I also had to switch off optimization for this. Not used for anything performance-critical
    ! anyway, so it does not matter. Future Olle will have to deal with this some other time.
    t0=walltime()
    if( verb .gt. 0 ) call lo_progressbar_init()
    do i=1,voro%n
        if ( present(dt) ) then
            if ( present(minjump) ) then
                call polyhedron_for_one_point(voro%cell(i)%polyhedron,dt,i,0,minjump)
            else
                call polyhedron_for_one_point(voro%cell(i)%polyhedron,dt,i,0)
            endif
            call classify_voronoi_cell(voro%cell(i),dt,i)
        else
            call polyhedron_for_one_point(voro%cell(i)%polyhedron,ddt,i,0)
            call classify_voronoi_cell(voro%cell(i),ddt,i)
        endif
        if ( verb .gt. 0 ) call lo_progressbar(' ... creating Voronoi tesselation',i,voro%n)
    enddo

    ! Extra, heavy classification
    if ( present(dt) ) then
        call check_diagram_connectivity(voro,dt)
    else
        call check_diagram_connectivity(voro,ddt)
    endif
    if ( verb .gt. 0 ) then
        write(*,*) '... got tesselation in ',tochar(walltime()-t0),'s'
    endif
end subroutine

!> Generate the Voronoi cell for a single point
subroutine polyhedron_for_one_point(polyhedron,dt,point,verbosity,minjump)
    !> the Voronoi polyhedron
    type(lo_polyhedron), intent(out) :: polyhedron
    !> the Distance table
    type(lo_distancetable), intent(in) :: dt
    !> which point in the distance table to calculate the diagram for
    integer, intent(in) :: point
    !> how much to talk
    integer, intent(in) :: verbosity
    !> jump to consider
    integer, intent(in), optional :: minjump

    ! Work variables
    ! after a slice, there will be something above, below and on the slicing plane.
    class(lo_geometryobject), allocatable :: ob_below,ob_above,ob_on
    ! linked list to keep track of all the faces
    type(lo_linked_list) :: list
    ! linked list to keep track of the linesegments that result from slices, and how to join them into polygons.
    type(lo_linked_list) :: segmentlist
    ! temporary polygon
    type(lo_polygon) :: dumpolygon
    ! the slicing plane
    type(lo_plane) :: plane

    ! Dummy variables
    real(flyt), dimension(:,:), allocatable :: dum1,dum2
    real(flyt), dimension(3) :: v0,v1
    real(flyt) :: f0,rmax,t0
    integer :: ctr_dropped,ctr_kept,ctr_modified,ctr_new,ctr_linesegment
    integer :: i,j,l,slicedecision,checklevel
    logical :: modifications

    ! initialize a polyhedron
    call polyhedron%createcube(4*dt%cutoff)
    ! Add the faces of the polyhedron to the linked list
    do i=1,polyhedron%n
        call list%add(polyhedron%face(i))
    enddo
    ! set the size of the bounding sphere
    rmax=dt%cutoff*10.0_flyt !lo_huge

    checklevel=0
    if ( verbosity .gt. 0 ) then
        t0=walltime()
        write(*,*) 'Built polyhedron as a box, slicing'
        checklevel=3 ! do a lot more sanity checks
    endif

    ! Start slicing with all the planes defined by half distances to neighbours, not including
    ! itself (that's why I start from i=2). They are sorted by distance, so that once the distance
    ! is larger than the bounding sphere, I can stop!
    do i=2,dt%particle(point)%n
        if ( present(minjump) ) then
            ! exclude some neighbours from the diagram generation
            if ( dt%particle(point)%jumps(i) .lt. minjump ) cycle
        endif


        ! Should I stop?
        list%current=>list%first
        if ( mod(i+1,5) .eq. 0 ) rmax=radius_of_bounding_sphere(list)
        f0=norm2(dt%particle(point)%v(:,i))*0.5_flyt
        if ( rmax-f0 .lt. -lo_tol ) then
            ! The plane is outside the bounding sphere of the polyhedron
            ! I have found the cell!
            exit
        endif

        ! Perhaps some info
        if ( verbosity .gt. 1 ) then
            l=list%countitems()
            write(*,*) 'Slicing with neighbour '//tochar(i)//&
                       ' out of '//tochar(dt%particle(point)%n)
            write(*,*) 'currently '//tochar(l)//' faces'
            !call listcontents(list)
            write(*,"(1X,A,F10.6,A,F10.6)") 'rmax: ',rmax,' distance to plane: ',f0
        endif
        ! Reset counters, mostly for debugging.
        ctr_dropped=0
        ctr_kept=0
        ctr_modified=0
        ctr_new=0
        ctr_linesegment=0
        ! construct the slicing plane
        v0=dt%particle(point)%v(:,i)
        v1=dt%particle(point)%v(:,i)*0.5_flyt
        call plane%generate(normal=v0,point=v1)
        if ( verbosity .gt. 1 ) write(*,*) '... got plane to slice with'

        ! Start slicing, first reset the list
        list%current=>list%first
        sliceeachface: do
            ! just make sure I have a clean slate
            if ( allocated(ob_below) ) deallocate(ob_below)
            if ( allocated(ob_on)    ) deallocate(ob_on)
            if ( allocated(ob_above) ) deallocate(ob_above)
            if ( verbosity .gt. 1 ) write(*,*) '... starting slicing polygon with plane'
            ! Ensure that the face is a polygon.
            select type(a=>list%current%obj)
            type is(lo_polygon)
                if ( verbosity .gt. 1 ) write(*,*) '... the face is a polygon, good'
                ! Make a decision what to do:
                ! slicedecision = -1 means not decided yet
                ! slicedecision = 0 means do nothing, iterate
                ! slicedecision = 1 means drop the face, remove it completely
                ! slicedecision = 2 means a face should be updated and points stored

                ! Test if I can make a quick decision:
                slicedecision=-1
                f0=plane%distance_to_point(a%centroid)
                if ( f0 .gt. a%rmax+lo_tol ) then
                    ! the face is far above the plane, drop it
                    slicedecision=1
                elseif ( f0 .lt. -a%rmax-lo_tol ) then
                    ! the face is way below the plane, no need to slice
                    slicedecision=0
                else
                    ! No fast decision, we have to do the slice
                    call a%slice(plane,ob_above,ob_on,ob_below,modifications,verbosity=0)
                    ! if something happened, do something about it.
                    if ( modifications ) then
                        select type(ob_below)
                        type is(lo_polygon)
                            ! something changed in this polygon
                            slicedecision=2
                        class default
                            ! This is not a polygon, then I remove it. I don't bother
                            ! with points and lines, they are redundant.
                            slicedecision=1
                        end select
                    else
                        ! no modifications means do nothing.
                        slicedecision=0
                    endif
                    ! Check if the object on the plane is a line segment, in that case I keep it.
                    select type(ob_on)
                    type is(lo_linesegment)
                        call segmentlist%add(ob_on)
                        ctr_linesegment=ctr_linesegment+1
                        ! maybe a sanity check
                        !if ( checklevel .gt. 0 ) call sanity_check_segment_on_plane(ob_on,plane)
                    end select
                endif

                ! Now I know what to do. Then do that.
                select case(slicedecision)
                    case(0) ! do nothing and iterate
                        ctr_kept=ctr_kept+1
                        if ( associated(list%current%next) ) then
                            list%current=>list%current%next
                        else
                            exit sliceeachface
                        endif
                    case(1) ! drop the face!, no iteration
                        ctr_dropped=ctr_dropped+1
                        call list%remove(j)
                    case(2) ! update the current face and iterate
                        ctr_modified=ctr_modified+1
                        call list%replace(ob_below)
                        ! Move forward in list
                        if ( associated(list%current%next) ) then
                            list%current=>list%current%next
                        else
                            exit sliceeachface
                        endif
                    case default
                        call lo_stop_gracefully(['No decision could be made. Should be impossible.'],lo_exitcode_param,__FILE__,__LINE__)
                end select
            end select

            if ( verbosity .gt. 1 ) then
                write(*,*) '... done slicing this face with neighbour ',tochar(i)
                write(*,*) ' '
            endif

        enddo sliceeachface

        ! Construct the new face: Collect all points
        if ( ctr_linesegment .gt. 1 ) then
            if ( verbosity .gt. 1 ) write(*,*) 'There should be a new face',ctr_linesegment
            ! make some space
            allocate(dum1(3,ctr_linesegment*2))
            dum1=0.0_flyt
            ! reset the segmentlist to be on the safe side. It should be empty though.
            segmentlist%current=>segmentlist%first
            ! extract all the points
            l=0
            do
                select type(a=>segmentlist%current%obj)
                    type is(lo_linesegment)
                        l=l+1
                        dum1(:,l)=a%r1
                        l=l+1
                        dum1(:,l)=a%r2
                end select
                ! iterate
                if ( associated(segmentlist%current%next) ) then
                    segmentlist%current=>segmentlist%current%next
                else
                    exit
                endif
            enddo
            ! Remove redundant points
            call lo_return_unique(dum1,dum2,lo_tol)
            ! Add it to the list of faces if it is a new face
            if ( size(dum2,2) .gt. 2 ) then
                ! create the polygon
                call dumpolygon%generate(dum2,plane)
                ! maybe check it
                !if ( checklevel .gt. 0 ) call sanity_check_polygon(dumpolygon)
                ! add it to the list
                call list%add(dumpolygon)
                ! clean up in the arrays
                if ( allocated(dumpolygon%r) ) deallocate(dumpolygon%r)
                ctr_new=ctr_new+1
            endif
        endif

        ! Clean up
        if ( allocated(dum1) ) deallocate(dum1)
        if ( allocated(dum2) ) deallocate(dum2)
        call segmentlist%destroy()
        if ( verbosity .gt. 1 ) write(*,*) 'done with face ',tochar(i)
    enddo

    ! Convert the list to a polyhedron
    lo_deallocate(polyhedron%face)
    polyhedron%n=list%countitems()
    lo_allocate(polyhedron%face(polyhedron%n))
    list%current=>list%first
    l=0
    do
        ! store this face
        l=l+1
        select type(a=>list%current%obj)
            type is(lo_polygon)
                polyhedron%face(l)%n=a%n
                lo_allocate(polyhedron%face(l)%r(3,a%n))
                polyhedron%face(l)%r=a%r
                do i=1,3
                    polyhedron%face(l)%centroid(i)=lo_mean(a%r(:,i))
                enddo
                polyhedron%face(l)%plane=a%plane
        end select
        ! iterate
        if ( associated(list%current%next) ) then
            list%current=>list%current%next
        else
            exit
        endif
    enddo

    ! Cleanup
    call list%destroy()

    ! Some helper functions
    contains
    !> Returns the current radius of the bounding sphere
    function radius_of_bounding_sphere(list) result(rmax)
        !> list of polygons
        type(lo_linked_list), intent(inout) :: list
        !> the largest radius
        real(flyt) :: rmax
        !
        integer :: l,i
        rmax=0.0_flyt
        if ( associated(list%first) ) then
            list%current=>list%first
            ! iterate the list, show what I have
            l=0
            do
                l=l+1
                !
                select type(a=>list%current%obj)
                type is(lo_polygon)
                    do i=1,a%n
                        rmax=max(rmax,norm2(a%r(:,i)))
                    enddo
                end select
                ! iterate
                if ( associated(list%current%next) ) then
                    list%current=>list%current%next
                else
                    exit
                endif
            enddo
        else
            ! the list might be empty
            rmax=lo_huge
        endif
    end function
end subroutine

!> Sort out the metadata for a voronoi cell
subroutine classify_voronoi_cell(cell,dt,particle)
    !> the cell
    type(lo_voronoi_diagram_cell), intent(inout) :: cell
    !> distance table
    type(lo_distancetable), intent(in) :: dt
    !> which particle in the distance table is this?
    integer, intent(in) :: particle
    !
    integer :: i,j,k,l
    integer, dimension(:,:), allocatable :: dumi,dumj
    real(flyt), dimension(:,:), allocatable :: dum
    real(flyt), dimension(3) :: v0
    real(flyt) :: f0

    ! Get the list of nodes

    ! first count total number of points in polyhedron
    l=0
    do i=1,cell%polyhedron%n
        do j=1,cell%polyhedron%face(i)%n
            l=l+1
        enddo
    enddo
    ! make a nice array of them
    lo_allocate(dum(3,l))
    l=0
    do i=1,cell%polyhedron%n
        do j=1,cell%polyhedron%face(i)%n
            l=l+1
            dum(:,l)=cell%polyhedron%face(i)%r(:,j)
        enddo
    enddo
    dum=lo_chop(dum,1E-12_flyt)
    ! keep the unique
    call lo_return_unique(dum,cell%r,lo_tol)

    lo_deallocate(dum)
    ! Now I have the nodes
    cell%nnodes=size(cell%r,2)
    allocate(cell%node(cell%nnodes))

    ! Build the faces
    cell%nfaces=cell%polyhedron%n
    lo_allocate(cell%face(cell%nfaces))
    do i=1,cell%nfaces
        ! how many points
        cell%face(i)%n=cell%polyhedron%face(i)%n
        lo_allocate(cell%face(i)%ind( cell%face(i)%n ))
        ! locate the unique nodes
        do j=1,cell%face(i)%n
            do k=1,cell%nnodes
                v0=cell%polyhedron%face(i)%r(:,j)-cell%r(:,k)
                if ( lo_sqnorm( v0 ) .lt. lo_sqtol ) then
                    ! this is the node it corresponds to
                    cell%face(i)%ind(j)=k
                    exit
                endif
            enddo
        enddo
    enddo

    ! For each face, figure out what particle is on the other side
    do i=1,cell%polyhedron%n
        ! go through the distance list and find a neighbour whose half-distance lies on the plane
        ! defining the polygon.
        do j=1,dt%particle(particle)%n
            ! probably new record in number of % on one line
            v0=dt%particle(particle)%v(:,j)*0.5_flyt
            f0=cell%polyhedron%face(i)%plane%distance_to_point(v0)
            if ( abs(f0) .lt. lo_tol ) then
                ! got it!
                cell%face(i)%neighbour=dt%particle(particle)%ind(j)
                cell%face(i)%neighbourvector=dt%particle(particle)%v(:,j)
                exit
            endif
        enddo
    enddo

    ! Get some more info on the nodes
    do i=1,cell%nnodes
        cell%node(i)%nsharedfaces=0
    enddo

    ! Count how many faces are joined at each node
    do i=1,cell%nfaces
        do j=1,cell%face(i)%n
            k=cell%face(i)%ind(j)
            cell%node(k)%nsharedfaces=cell%node(k)%nsharedfaces+1
        enddo
    enddo
    ! make space
    do i=1,cell%nnodes
        lo_allocate(cell%node(i)%faceind( cell%node(i)%nsharedfaces ))
        cell%node(i)%nsharedfaces=0
    enddo
    ! Store the info about which node is in which face
    do i=1,cell%nfaces
        do j=1,cell%face(i)%n
            k=cell%face(i)%ind(j)
            cell%node(k)%nsharedfaces=cell%node(k)%nsharedfaces+1
            cell%node(k)%faceind( cell%node(k)%nsharedfaces ) = i
        enddo
    enddo

    ! Get a list of all the edges
    ! count edges
    l=0
    do i=1,cell%nfaces
        l=l+cell%face(i)%n
    enddo
    ! make space
    allocate(dumi(2,l))
    ! store all edges
    l=0
    do i=1,cell%nfaces
        do j=1,cell%face(i)%n
            l=l+1
            dumi(1,l)=cell%face(i)%ind(j)
            k=lo_index_in_periodic_array(j+1,cell%face(i)%n)
            dumi(2,l)=cell%face(i)%ind(k)
            ! sort it
            call qsort(dumi(:,l))
        enddo
    enddo
    ! get the unique edges
    call lo_return_unique(dumi,dumj)
    cell%nedges=size(dumj,2)
    lo_allocate(cell%edge(cell%nedges))
    do i=1,cell%nedges
        cell%edge(i)%i1=dumj(1,i)
        cell%edge(i)%i2=dumj(2,i)
    enddo

    ! Store info about the planes, and get the smallest distance
    cell%rmin=lo_huge
    do i=1,cell%nfaces
        cell%face(i)%plane=cell%polyhedron%face(i)%plane
        cell%rmin=min(cell%rmin,abs(cell%face(i)%plane%p))
    enddo
    cell%rmin=cell%rmin-2.0_flyt*lo_tol
    ! And largest distance
    cell%rmax=-lo_huge
    do i=1,cell%nnodes
        cell%rmax=max(cell%rmax,norm2(cell%r(:,i)))
    enddo
    cell%rmax=cell%rmax+2.0_flyt*lo_tol
    !
end subroutine

!> Figure out how all cells are connected, sort of
subroutine check_diagram_connectivity(voro,dt)
    !> the cell
    type(lo_voronoi_diagram), intent(inout) :: voro
    !> distance table
    type(lo_distancetable), intent(in) :: dt
    !
    integer :: i,j,k,l,cell,node,ctr
    real(flyt), dimension(3) :: v0,v1
    integer, dimension(:), allocatable :: dum

    ! Figure out how nodes are shared
    do cell=1,voro%n
        lo_allocate(dum(dt%particle(cell)%n))
        do node=1,voro%cell(cell)%nnodes
            ! For each node, check neighbouring cells if they also have this node.
            voro%cell(cell)%node(node)%nsharedcells=0
            ! go to each neighbour
            ctr=0
            v0=voro%cell(cell)%r(:,node)
            do j=2,dt%particle(cell)%n ! skip itself
                ! vector to neighbouring cell minus reference
                v1=dt%particle(cell)%v(:,j)-v0
                ! index to neighbouring cell
                k=dt%particle(cell)%ind(j)
                ! can this neighbouring cell be relevant?
                if ( lo_sqnorm(v1) .lt. voro%cell(k)%rmax**2 ) then
                    do l=1,voro%cell(k)%nnodes
                        ! coordinate to the node in the frame of reference of the original cell
                        ctr=ctr+1
                        dum(ctr)=j
                        exit
                    enddo
                endif
            enddo
            ! Store the information
            voro%cell(cell)%node(node)%nsharedcells=ctr
            lo_allocate(voro%cell(cell)%node(node)%neighbourind( ctr ))
            lo_allocate(voro%cell(cell)%node(node)%neighbourvector( 3,ctr ))
            do i=1,ctr
                j=dum(i)
                voro%cell(cell)%node(node)%neighbourind(i)=dt%particle(cell)%ind(j)
                voro%cell(cell)%node(node)%neighbourvector(:,i)=dt%particle(cell)%v(:,j)
            enddo
        enddo
        lo_deallocate(dum)
    enddo

    ! With this info, I can figure out how edges are shared as well
    do cell=1,voro%n
        do i=1,voro%cell(cell)%nedges
            j=voro%cell(cell)%edge(i)%i1
            k=voro%cell(cell)%edge(i)%i2
            ! the ones that exist in both is what matters!
            voro%cell(cell)%edge(i)%nshared=count_vector_intersection( voro%cell(cell)%node(j)%neighbourvector, voro%cell(cell)%node(k)%neighbourvector )
        enddo
    enddo

    contains
    !> Return the number of intersections between two lists of vectors
    pure function count_vector_intersection(l1,l2) result(n)
        real(flyt), dimension(:,:), intent(in) :: l1,l2
        integer :: n
        !
        integer :: i,j
        real(flyt), dimension(3) :: v
        n=0
        do i=1,size(l1,2)
        do j=1,size(l2,2)
            v=l1(:,i)-l2(:,j)
            if ( lo_sqnorm(v) .lt. lo_sqtol ) then
                n=n+1
            endif
        enddo
        enddo
    end function
end subroutine


#ifdef AGRESSIVE_SANITY
!> Size in memory, in bytes
function node_size_in_mem(n) result(mem)
    !> node
    class(lo_voronoi_diagram_node), intent(in) :: n
    !> size in bytes
    integer :: mem

    mem=0
    mem=mem+storage_size(n)
    if ( allocated(n%faceind) ) mem=mem+storage_size(n%faceind)*size(n%faceind)
    if ( allocated(n%neighbourind) ) mem=mem+storage_size(n%neighbourind)*size(n%neighbourind)
    if ( allocated(n%neighbourvector) ) mem=mem+storage_size(n%neighbourvector)*size(n%neighbourvector)
    mem=mem/8
end function
!> Size in memory, in bytes
function face_size_in_mem(f) result(mem)
    !> node
    class(lo_voronoi_diagram_face), intent(in) :: f
    !> size in bytes
    integer :: mem

    mem=0
    mem=mem+storage_size(f)
    if ( allocated(f%ind) ) mem=mem+storage_size(f%ind)*size(f%ind)
    mem=mem/8
end function
!> Size in memory, in bytes
function cell_size_in_mem(c) result(mem)
    !> node
    class(lo_voronoi_diagram_cell), intent(in) :: c
    !> size in bytes
    integer :: mem,i

    mem=0
    mem=mem+storage_size(c)
    if ( allocated(c%r) ) mem=mem+storage_size(c%r)*size(c%r)
    mem=mem/8

    if ( allocated(c%face) ) then
        do i=1,size(c%face)
            mem=mem+c%face(i)%size_in_mem()
        enddo
    endif
    if ( allocated(c%node) ) then
        do i=1,size(c%node)
            mem=mem+c%node(i)%size_in_mem()
        enddo
    endif
    if ( allocated(c%edge) ) then
        mem=mem+storage_size(c%edge)*size(c%edge)/8
    endif
end function
!> Size in memory, in bytes
function diagram_size_in_mem(d) result(mem)
    !> diagram
    class(lo_voronoi_diagram), intent(in) :: d
    !> size in bytes
    integer :: mem,i

    mem=0
    mem=mem+storage_size(d)
    mem=mem/8
    if ( allocated(d%cell) ) then
        do i=1,size(d%cell)
            mem=mem+d%cell(i)%size_in_mem()
        enddo
    endif
end function
#endif


! subroutine dump_polyhedron_for_debugging(list,n)
!    type(lo_linked_list), intent(inout) :: list
!    integer, intent(in) :: n
!    !
!    integer :: l,u,ll,i
!    !
!    u=open_file('out','polygon_debug_'//trim(int2char(n))//'.m')
!    if ( associated(list%first) ) then
!        list%current=>list%first
!        ! iterate the list, show what I have
!        l=0
!        ll=0
!        do
!            l=l+1
!            select type(a=>list%current%obj)
!            type is(lo_polygon)
!                do i=1,a%n
!                    ll=ll+1
!                    write(u,*) 'node(',ll,',:)=[' !];
!                    write(u,*) a%r(:,i)
!                    write(u,*) '];'
!                    write(u,*) 'face{',l,'}(',i,')=[',ll,'];'
!                enddo
!            end select
!            ! iterate
!            if ( associated(list%current%next) ) then
!                list%current=>list%current%next
!            else
!                exit
!            endif
!        enddo
!    else
!        write(*,*) 'list is empty'
!    endif
!    close(u)
! end subroutine

! subroutine listcontents(list)
!     type(lo_linked_list), intent(inout) :: list
!     !
!     integer :: l
!     !
!     if ( associated(list%first) ) then
!         list%current=>list%first
!         ! iterate the list, show what I have
!         l=0
!         do
!             l=l+1
!             select type(a=>list%current%obj)
!             type is(lo_polygon)
!                 write(*,*) 'element '//tochar(l)//' is a polygon'
!             class default
!                 write(*,*) 'element '//tochar(l)//' is not a polygon'
!             end select
!             ! iterate
!             if ( associated(list%current%next) ) then
!                 list%current=>list%current%next
!             else
!                 exit
!             endif
!         enddo
!     else
!         write(*,*) 'list is empty'
!     endif
! end subroutine

! subroutine sanity_check_segment_on_plane(segment,plane)
!     type(lo_plane), intent(in) :: plane
!     type(lo_linesegment), intent(in) :: segment
!     ! is it longer than 0?
!     if ( norm2(segment%r1-segment%r2) .lt. lo_tol ) then
!         write(*,*) 'linesegment is a point'
!         stop
!     endif
!     ! is it on the plane?
!     if ( abs(plane%distance_to_point(segment%r1)) .gt. lo_tol .or. &
!          abs(plane%distance_to_point(segment%r2)) .gt. lo_tol ) then
!         write(*,*) 'linesegment is not on the plane'
!         stop
!     endif
! end subroutine

! subroutine sanity_check_polygon(polygon)
!     type(lo_polygon), intent(in) :: polygon
!     !
!     integer :: i
!     ! is it enough points?
!     if ( polygon%n .lt. 3 ) then
!         write(*,*) 'polygon has only '//tochar(polygon%n)//' points'
!         stop
!     endif
!     ! are the points on the plane?
!     do i=1,polygon%n
!         if ( abs(polygon%plane%distance_to_point(polygon%r(:,i))) .gt. lo_tol ) then
!             write(*,*) 'polygon is not on the plane it is supposed to be on.'
!             stop
!         endif
!     enddo
! end subroutine

! subroutine listcontents(list)
!     type(lo_linked_list), intent(inout) :: list
!     !
!     integer :: l
!     !
!     if ( associated(list%first) ) then
!         list%current=>list%first
!         ! iterate the list, show what I have
!         l=0
!         do
!             l=l+1
!             select type(a=>list%current%obj)
!             type is(lo_polygon)
!                 write(*,*) 'element '//tochar(l)//' is a polygon'
!             class default
!                 write(*,*) 'element '//tochar(l)//' is not a polygon'
!             end select
!             ! iterate
!             if ( associated(list%current%next) ) then
!                 list%current=>list%current%next
!             else
!                 exit
!             endif
!         enddo
!     else
!         write(*,*) 'list is empty'
!     endif
! end subroutine

end module
