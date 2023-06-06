#include "precompilerdefinitions"
submodule (geometryfunctions) geometryfunctions_objectprocedures
implicit none
contains

!> Distance between line segment and a point
#ifdef AGRESSIVE_SANITY
module function distance_linesegment_point(linesegment,point) result(d)
#else
module pure function distance_linesegment_point(linesegment,point) result(d)
#endif
    !> the line segment
    class(lo_linesegment), intent(in) :: linesegment
    !> the point
    real(flyt), dimension(3), intent(in) :: point
    !> the distance
    real(flyt) :: d
    !
    real(flyt), dimension(3) :: u,v,p
    real(flyt) :: a
    !
    u=linesegment%r2-linesegment%r1
    v=point-linesegment%r1
#ifdef AGRESSIVE_SANITY
    if ( dot_product(u,u) .lt. lo_sqtol ) then
        call lo_stop_gracefully(['Trying to calculate a linesegment point distance, but the linesegment is a point. Probably not what you want.'],lo_exitcode_param)
    endif
#endif
    a=dot_product(u,v)/dot_product(u,u)
    if ( a .lt. 0.0_flyt ) then
        d=norm2(point-linesegment%r1)
    elseif ( a .gt. 1.0_flyt ) then
        d=norm2(point-linesegment%r2)
    else
        p=linesegment%r1+a*u
        d=norm2(point-p)
    endif
end function

!> Translate an object with a vector
module recursive subroutine translate_geometrical_object(g,translation)
    !> Object to be moved
    class(lo_geometryobject), intent(inout) :: g
    !> The translation vector
    real(flyt), dimension(3), intent(in) :: translation

    integer :: i
    ! What kind of object is it?
    select type(g)
        type is (lo_null)
            ! do nothing
        type is (lo_point)
            ! move the point
            g%r=g%r+translation
        type is (lo_line)
            ! translate the points
            g%r1=g%r1+translation
            g%r2=g%r2+translation
        type is (lo_linesegment)
            ! translate the points
            g%r1=g%r1+translation
            g%r2=g%r2+translation
        type is (lo_plane)
            ! translate the plane
            g%p=-dot_product(g%normal,g%normal+g%p+translation)
        type is (lo_polygon)
            ! translate the polygons
            do i=1,g%n
                g%r(:,i)=g%r(:,i)+translation
            enddo
            g%centroid=g%centroid+translation
            g%plane%p=-dot_product(g%plane%normal,g%plane%normal+g%plane%p+translation)
        type is (lo_polyhedron)
            ! translate the polyhedron
            do i=1,g%n
                ! apparently, one can do this in a cool recursive way
                call g%face(i)%translate(translation)
            enddo
        class default
            call lo_stop_gracefully(["no translation defined for this class"],lo_exitcode_param)
    end select
end subroutine

!> Signed distance from plane to a point.
pure module function lo_plane_point_distance(plane,point) result(r)
    !> The plane
    class(lo_plane), intent(in) :: plane
    !> The point
    real(flyt), dimension(3), intent(in) :: point
    !> The (signed) distance
    real(flyt) :: r
    !
    r=dot_product(plane%normal,point)+plane%p
end function

!> Sort points clockwise with reference to this plane
module subroutine anglesort_with_plane(plane,points,order)
    !> The plane the points will get projected to
    class(lo_plane), intent(in) :: plane
    !> Points to sort
    real(flyt), dimension(:,:), intent(inout) :: points
    !> The order of the points
    integer, dimension(:), intent(out), optional :: order
    !
    integer :: i,n
    real(flyt), dimension(3) :: ctr
    real(flyt), dimension(:,:), allocatable :: dum
    real(flyt), dimension(:), allocatable :: angle
    integer, dimension(:), allocatable :: dumind
    !
    n=size(points,2)
    lo_allocate(dum(2,n))
    lo_allocate(angle(n))
    lo_allocate(dumind(n))
    ! get the central point, I don't trust the centriod.
    do i=1,3
        ctr(i)=lo_mean(points(i,:))
    enddo
    ! project onto plane
    do i=1,n
        dum(1,i)=dot_product(plane%v1,points(:,i)-ctr)
        dum(2,i)=dot_product(plane%v2,points(:,i)-ctr)
        angle(i)=atan2(dum(1,i),dum(2,i))
    enddo
    ! Get a signed angle
    call qsort(angle,dumind)
    points=points(:,dumind)
    if ( present(order) ) order=dumind
    ! Cleanup
    lo_deallocate(dum)
    lo_deallocate(angle)
    lo_deallocate(dumind)
end subroutine

!> Reflection matrix for a plane. Assumes the plane passes through the origin.
pure module function reflection_matrix(plane) result(m)
    !> The plane with which to reflect
    class(lo_plane), intent(in) :: plane
    !> the reflection matrix
    real(flyt), dimension(3,3) :: m

    real(flyt) :: a,b,c
    a=plane%normal(1)
    b=plane%normal(2)
    c=plane%normal(3)
    m(1,1)=1.0_flyt-2.0_flyt*a**2; m(1,2)=-2.0_flyt*a*b; m(1,3)=-2.0_flyt*a*c
    m(2,1)=-2.0_flyt*a*b; m(2,2)=1.0_flyt-2.0_flyt*b**2; m(2,3)=-2.0_flyt*b*c
    m(3,1)=-2.0_flyt*a*c; m(3,2)=-2.0_flyt*b*c; m(3,3)=1.0_flyt-2.0_flyt*c**2
end function

!> Unsigned area of a convex polygon.
module pure function area_of_polygon(polygon) result(A)
    !> The polygon
    class(lo_polygon), intent(in) :: polygon
    !> The area
    real(flyt) :: A

    real(flyt), dimension(:,:), allocatable :: dum
    integer :: i

    allocate(dum(2,polygon%n+1))
    ! project onto plane
    do i=1,polygon%n
        dum(1,i)=dot_product(polygon%plane%v1,polygon%r(:,i)-polygon%centroid)
        dum(2,i)=dot_product(polygon%plane%v2,polygon%r(:,i)-polygon%centroid)
    enddo
    dum(:,polygon%n+1)=dum(:,1)
    ! get the area
    A=0.0_flyt
    do i=1,polygon%n
        A=A+dum(1,i)*dum(2,i+1)-dum(1,i+1)*dum(2,i)
    enddo
    A=abs(A)*0.5_flyt
    deallocate(dum)
end function

!> Volume of convex polyhedron
pure module function volume_of_polyhedron(polyhedron) result(V)
    !> The polyhedron
    class(lo_polyhedron), intent(in) :: polyhedron
    !> The volume
    real(flyt) :: V
    !
    integer :: i
    !
    V=0.0_flyt
    do i=1,polyhedron%n
        V=V+abs(polyhedron%face(i)%plane%p)*polyhedron%face(i)%area()/3.0_flyt
    enddo
end function

!> Check if a point is inside a polyhedron. Any point within tol from the edge is considered inside
pure module function is_point_inside_polyhedron(polyhedron,point,tol) result(inside)
    !> The polyhedron
    class(lo_polyhedron), intent(in) :: polyhedron
    !> the point
    real(flyt), dimension(3), intent(in) :: point
    !> The tolerance
    real(flyt), intent(in) :: tol
    !> Is the point inside or not
    logical :: inside

    integer :: i

    inside=.true.
    do i=1,polyhedron%n
        ! if it's outside any face, it's false
        if ( polyhedron%face(i)%plane%distance_to_point(point) .gt. tol ) then
            inside=.false.
            return
        endif
    enddo
end function

end submodule
