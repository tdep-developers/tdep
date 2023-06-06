#include "precompilerdefinitions"
submodule (geometryfunctions) geometryfunctions_slicingroutines
implicit none
contains

!> Divide an object with a plane. Returns three new objects, the part above the plane, the part on the plane and the part below the plane. If there are points on the plane, they will be in all three objects.
module subroutine slice_geometrical_object(g,plane,aboveplane,onplane,belowplane,modifications,verbosity)
    !> the object to be sliced
    class(lo_geometryobject), intent(in) :: g
    !> the plane to slice with
    type(lo_plane), intent(in) :: plane
    !> The object above the plane
    class(lo_geometryobject), intent(out), allocatable :: aboveplane
    !> Object below the plane
    class(lo_geometryobject), intent(out), allocatable :: belowplane
    !> Object on the plane
    class(lo_geometryobject), intent(out), allocatable :: onplane
    !> Did the object actually get sliced?
    logical, intent(out) :: modifications
    !> How much should be written?
    integer, intent(in), optional :: verbosity
    !
    integer :: verb
    if ( present(verbosity) ) then
        verb=verbosity
    else
        verb=0
    endif
    !> @todo Obivously, add routines for more things than polygons. But it should be much easier.
    select type(g)
        type is (lo_polygon)
            call slice_polygon_with_plane(g,plane,aboveplane,onplane,belowplane,modifications,verb)
        class default
            call lo_stop_gracefully(['have not implemented slicing for all objects yet'],lo_exitcode_param,__FILE__,__LINE__)
    end select
end subroutine

!> Divide a polygon with a plane. Will return three objects, the part above, below and on the plane. The points that are on the plane will be included in all three objects.
subroutine slice_polygon_with_plane(polygon,plane,aboveplane,onplane,belowplane,modifications,verbosity)
    !> The polygon
    type(lo_polygon), intent(in) :: polygon
    !> The plane
    type(lo_plane),  intent(in) :: plane
    !> The resulting objects above, on and below the plane
    class(lo_geometryobject), allocatable, intent(out) :: aboveplane,onplane,belowplane
    !> Shorthand to keep track if something changed or not.
    logical, intent(out) :: modifications
    !> How much to talk
    integer, intent(in) :: verbosity

    integer :: i,j,l,n
    integer :: prev,curr,nmax
    integer, dimension(:), allocatable :: above,below,on,slice,ind
    real(flyt) :: t0
    real(flyt), dimension(3) :: v0
    real(flyt), dimension(3,2) :: slicedpoints
    real(flyt), dimension(:,:), allocatable :: dum
    type(lo_line) :: line
    !
    t0=0.0_flyt
    if ( verbosity .gt. 0 ) then
        t0=walltime()
        write(*,*) 'slice_polygon_with_plane'
        if ( verbosity .gt. 1 ) then
            write(*,*) '... the polygon to be sliced:'
            write(*,*) '    edges: ',tochar(polygon%n)
            do i=1,polygon%n
                write(*,*) '   point',i,real(polygon%r(:,i))
            enddo
            write(*,*) '... with this plane:'
            write(*,*) 'plane'
        endif
    endif
    ! angle sort the polygon
    !
    nmax=polygon%n+2
    ! make some space for the classifiers
    lo_allocate(above(nmax))
    lo_allocate(below(nmax))
    lo_allocate(on(nmax))
    lo_allocate(slice(nmax))
    ! sort the points into above, on and below.
    on=0
    above=0
    below=0
    slice=0
    ! the first point on the first line
    v0=polygon%r(:,polygon%n)
    prev=planecondition(plane,v0)
    l=0
    do i=1,polygon%n
        ! where is this point in relation to the plane? -1 is below, 0 on, and 1 above
        curr=planecondition(plane,polygon%r(:,i))
        ! should it be sliced?
        if ( prev*curr .eq. -1 ) then
            slice(i)=1
            l=l+1
            ! small sanity check, will get caught later if not true
            if ( l .le. 2 ) then
                j=lo_index_in_periodic_array(i-1,polygon%n)
                call line%generate(polygon%r(:,i),polygon%r(:,j))
                slicedpoints(:,l)=quick_intersection(plane,line)
            endif
        endif
        ! locate it with respect to the plane
        select case(curr)
            case(-1)
                below(i)=1
            case(0)
                on(i)=1
            case(1)
                above(i)=1
        end select
        !
        prev=curr
    enddo

    ! What can come of this?
    if ( verbosity .gt. 1 ) then
        write(*,*) '... classified points of polygon:'
        write(*,*) '       above:',sum(above)
        write(*,*) '       below:',sum(below)
        write(*,*) '          on:',sum(on)
        write(*,*) '       slice:',sum(slice)
        write(*,*) ' slicedlines:',l
        if ( l .gt. 2 ) then
            write(*,*) 'POLYGON SLICED IN MORE THAN TWO PLACES, BAD!'
            write(*,*) 'dumping some stuff to debug.'
            write(*,*) 'plane'
            write(*,*) plane%normal
            write(*,*) plane%p
            write(*,*) 'polygon'
            do i=1,polygon%n
                write(*,"(1X,3(1X,F12.8))") polygon%r(:,i)
            enddo
            write(*,*) 'polygon distance to plane'
            do i=1,polygon%n
                write(*,*) i,plane%distance_to_point(polygon%r(:,i))
            enddo
            write(*,*) 'slice points:'
            do i=1,l
                write(*,"(1X,3(1X,F12.8))") slicedpoints(:,l)
            enddo
            stop
        endif
    endif
    ! Everything above
    if ( sum(above) .eq. polygon%n ) then
        lo_allocate(lo_null::belowplane)
        lo_allocate(lo_null::onplane)
        allocate(aboveplane,source=polygon,stat=lo_status)
        if ( lo_status .ne. 0 ) then
            write(*,*) 'failed allocation in slicing routine'
            stop
        endif
        modifications=.true.
        return
    endif
    ! Nothing above
    if ( sum(above) .eq. 0 ) then
        allocate(belowplane,source=polygon,stat=lo_status)
        if ( lo_status .ne. 0 ) then
            write(*,*) 'failed allocation in slicing routine'
            stop
        endif
        lo_allocate(lo_null::onplane)
        lo_allocate(lo_null::aboveplane)
        modifications=.false.
        return
    endif

    ! Everything below
    if ( sum(below) .eq. polygon%n ) then
        allocate(belowplane,source=polygon,stat=lo_status)
        if ( lo_status .ne. 0 ) then
            write(*,*) 'failed allocation in slicing routine'
            stop
        endif
        lo_allocate(lo_null::onplane)
        lo_allocate(lo_null::aboveplane)
        modifications=.false.
        return
    endif
    ! Everything on the plane
    if ( sum(on) .eq. polygon%n ) then
        lo_allocate(lo_null::belowplane)
        allocate(onplane,source=polygon,stat=lo_status)
        if ( lo_status .ne. 0 ) then
            write(*,*) 'falied allocation in slicing routine'
            stop
        endif
        lo_allocate(lo_null::aboveplane)
        modifications=.false.
        return
    endif

    ! Get a list of all the points
    n=polygon%n+sum(slice)
    lo_allocate(dum(3,n))
    lo_allocate(ind(n))
    dum(:,1:polygon%n)=polygon%r
    if ( sum(slice) .gt. 0 ) then
        do i=1,sum(slice)
            dum(:,polygon%n+i)=slicedpoints(:,i)
            ! the sliced points should be in all polygons
            above(polygon%n+i)=1
            on(polygon%n+i)=1
            below(polygon%n+i)=1
        enddo
    endif

    ! Get points below:
    l=0
    ind=0
    do i=1,n
        if ( below(i) .eq. 1 .or. on(i) .eq. 1 ) then
            l=l+1
            ind(l)=i
        endif
    enddo
    if ( l .gt. 0 ) then
        call object_from_planar_points(belowplane,dum(:,ind(1:l)),polygon%plane)
    else
        lo_allocate(lo_null::belowplane)
    endif

    ! Get points on:
    l=0
    ind=0
    do i=1,n
        if ( on(i) .eq. 1 ) then
            l=l+1
            ind(l)=i
        endif
    enddo
    if ( l .gt. 0 ) then
        call object_from_planar_points(onplane,dum(:,ind(1:l)),plane)
    else
        lo_allocate(lo_null::onplane)
    endif

    ! Get points above:
    l=0
    ind=0
    do i=1,n
        if ( above(i) .eq. 1 ) then
            l=l+1
            ind(l)=i
        endif
    enddo
    if ( l .gt. 0 ) then
        call object_from_planar_points(aboveplane,dum(:,ind(1:l)),polygon%plane)
    else
        lo_allocate(lo_null::aboveplane)
    endif

    ! If I reached this point, something has changed
    modifications=.true.

    if ( verbosity .gt. 0 ) then
        write(*,*) 'done slice_polygon_with_plane in',walltime()-t0
    endif
    ! And we are done!
end subroutine

! Test where a point is in reference to a plane.
! 0 is on the plane, -1 below and +1 above.
pure function planecondition(plane,point) result(cnd)
    type(lo_plane), intent(in) :: plane
    real(flyt), dimension(3), intent(in) :: point
    integer :: cnd
    !
    real(flyt) :: f0
    f0=plane%distance_to_point(point)
    if ( abs(f0) .lt. lo_tol ) then
        cnd=0
    elseif ( f0 .gt. 0.0_flyt ) then
        cnd=1
    else
        cnd=-1
    endif
end function

! This is a special case of line-plane intersection, since I have already taken care
! of the edge cases once this function is called.
pure function quick_intersection(plane,line) result(point)
    type(lo_plane), intent(in) :: plane
    type(lo_line), intent(in) :: line
    real(flyt), dimension(3) :: point
    !
    real(flyt), dimension(3) :: dv
    real(flyt) :: f0,f1
    dv=line%r2-line%r1
    f0=dot_product(dv,plane%normal)
    f1=-(plane%p+dot_product(line%r1,plane%normal))
    point=line%r1+(f1/f0)*dv
end function

!> Create point,linesegment or polygon from a set of points in the same plane.
subroutine object_from_planar_points(g,r,plane)
    class(lo_geometryobject), intent(out), allocatable :: g
    real(flyt), dimension(:,:), intent(in) :: r
    type(lo_plane), intent(in) :: plane

    integer :: n
#ifdef AGRESSIVE_SANITY
    real(flyt), dimension(:,:), allocatable :: dum
    integer :: i
    do i=1,size(r,2)
        if ( abs(plane%distance_to_point(r(:,i))) .gt. lo_tol ) then
            call lo_stop_gracefully(['Trying to create an object from planar points, but they are not on the same plane.'],lo_exitcode_param)
        endif
    enddo
    n=size(r,2)
    call lo_return_unique(r,dum)
    if ( size(dum,2) .ne. n ) then
        call lo_stop_gracefully(['Trying to create an object from planar points, but there are redundant points.'],lo_exitcode_param)
    endif
#endif

    n=size(r,2)
    if ( n .lt. 1 ) then
        ! it's nothing
        allocate(lo_null::g)
    elseif ( n .eq. 1 ) then
        ! it's a point
        allocate(lo_point::g)
        select type(g)
            type is (lo_point)
            g%r=r(:,1)
        end select
    elseif ( n .eq. 2 ) then
        ! a line segment
        allocate(lo_linesegment::g)
        select type(g)
            type is (lo_linesegment)
            call g%generate(r(:,1),r(:,2))
        end select
    else
        ! a polygon
        allocate(lo_polygon::g)
        select type(g)
        type is(lo_polygon)
            call g%generate(r,plane)
        end select
    endif
end subroutine

end submodule
