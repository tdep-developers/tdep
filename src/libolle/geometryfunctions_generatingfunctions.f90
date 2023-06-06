#include "precompilerdefinitions"
submodule (geometryfunctions) geometryfunctions_generatingfunctions
implicit none
contains

!> Initialize a polyhedron as a cube
module subroutine createcube(polyhedron,side)
    !> The polyhedron
    class(lo_polyhedron), intent(inout) :: polyhedron
    !> Side of cube
    real(flyt), intent(in) :: side
    !
    real(flyt), dimension(3,8) :: corners
    real(flyt), dimension(3,6) :: normals
    integer, dimension(4,6) :: faces
    integer :: i,j,k,l

    ! a cube has 6 faces
    polyhedron%n=6
    lo_allocate(polyhedron%face( polyhedron%n ) )

    ! Get the corners
    l=0
    do i=0,1
    do j=0,1
    do k=0,1
        l=l+1
        corners(:,l)=((/i,j,k/)*1.0_flyt-0.5_flyt)*side
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
    normals(:,1)=[-1.0_flyt, 0.0_flyt, 0.0_flyt]
    normals(:,2)=[ 0.0_flyt,-1.0_flyt, 0.0_flyt]
    normals(:,3)=[ 0.0_flyt, 0.0_flyt,-1.0_flyt]
    normals(:,4)=[ 1.0_flyt, 0.0_flyt, 0.0_flyt]
    normals(:,5)=[ 0.0_flyt, 1.0_flyt, 0.0_flyt]
    normals(:,6)=[ 0.0_flyt, 0.0_flyt, 1.0_flyt]

    do i=1,polyhedron%n
        ! store the points
        polyhedron%face(i)%n=4
        lo_allocate(polyhedron%face(i)%r(3,4))
        polyhedron%face(i)%r=corners(:,faces(:,i))
        ! get the plane that the planes are on
        call polyhedron%face(i)%plane%generate(normal=normals(:,i),point=corners(:,faces(1,i)))
        call polyhedron%face(i)%plane%anglesort(polyhedron%face(i)%r)
        ! the centroid
        do j=1,3
            polyhedron%face(i)%centroid(j)=lo_mean(polyhedron%face(i)%r(j,:))
        enddo
        ! the radius of bounding sphere
        polyhedron%face(i)%rmax=lo_huge
    enddo
    ! And it's done
end subroutine

!> Creates a plane, either from three points or a point and a normal.
module subroutine create_plane_from_points(plane,points,point,normal)
    !> The plane
    class(lo_plane), intent(out) :: plane
    !> Three points defining a plane
    real(flyt), dimension(3,3), intent(in), optional :: points
    !> Point in the plane
    real(flyt), dimension(3), intent(in), optional :: point
    !> Normal of the plane
    real(flyt), dimension(3), intent(in), optional :: normal
    !
    real(flyt), dimension(3) :: v0,v1
    real(flyt), dimension(3,6), parameter :: stupidpoints=reshape([&
        1.0_flyt,0.0_flyt,0.0_flyt, &
        0.0_flyt,1.0_flyt,0.0_flyt, &
        0.0_flyt,0.0_flyt,1.0_flyt, &
        0.0_flyt,0.7071067811865475_flyt,0.7071067811865475_flyt,&
        0.7071067811865475_flyt,0.0_flyt,0.7071067811865475_flyt,&
        0.7071067811865475_flyt,0.7071067811865475_flyt,0.0_flyt],[3,6])
    integer :: i

    ! If it's defined from three points
    if ( present(points) ) then
        v0=points(:,1)-points(:,3)
        v1=points(:,2)-points(:,3)
        plane%normal=lo_cross(v0,v1)
#ifdef AGRESSIVE_SANITY
    if ( norm2(plane%normal) .lt. lo_tol ) then
        call lo_stop_gracefully(['Trying to create a plane from points that badly define a plane'],lo_exitcode_param)
    endif
#endif
        plane%normal=plane%normal/norm2(plane%normal)
        plane%p=-dot_product(plane%normal,points(:,3))
    endif

    ! Define it by a normal and a point. I don't assume that the
    ! normal is of norm 1.
    if ( present(point) .and. present(normal) ) then
#ifdef AGRESSIVE_SANITY
    if ( norm2(normal) .lt. lo_tol ) then
        call lo_stop_gracefully(['Trying to create a plane from points that badly define a plane'],lo_exitcode_param)
    endif
#endif
        plane%normal=normal/norm2(normal)
        plane%p=-dot_product(plane%normal,point)
    endif

    ! Create the vectors that span the plane. First get a vector orthogonal
    ! To the normal. I want to keep the routine pure, so I generate some strange numbers.
    ! It's impossible that they are all parallel to the normal.
    do i=1,6
        v1=lo_cross(plane%normal,stupidpoints(:,i))
        if ( lo_sqnorm(v1) .gt. lo_tol ) exit
    enddo
    v1=v1/norm2(v1)
    ! And then a third orthogonal vector
    v0=lo_cross(v1,plane%normal)
    v0=v0/norm2(v0)
    ! I think this becoms a right handed system.
    plane%v1=v0
    plane%v2=v1
end subroutine

!> Generate a line or line segment from two points
module subroutine create_line_from_points(line,point1,point2)
    !> the resulting object
    class(lo_line), intent(out) :: line
    !> first point
    real(flyt), dimension(3), intent(in) :: point1
    !> second point
    real(flyt), dimension(3), intent(in) :: point2
    !
#ifdef AGRESSIVE_SANITY
    if ( norm2(point1-point2) .lt. lo_tol ) then
        call lo_stop_gracefully(['Trying to create a line from points that are too close to each other'],lo_exitcode_param)
    endif
#endif
    line%r1=point1
    line%r2=point2
end subroutine

!> Generate a line or line segment from two points
module subroutine create_linesegment_from_points(line,point1,point2)
    !> the resulting object
    class(lo_linesegment), intent(out) :: line
    !> first point
    real(flyt), dimension(3), intent(in) :: point1
    !> second point
    real(flyt), dimension(3), intent(in) :: point2
    !
#ifdef AGRESSIVE_SANITY
    if ( norm2(point1-point2) .lt. lo_tol ) then
        call lo_stop_gracefully(['Trying to create a linesegment from points that are too close to each other'],lo_exitcode_param)
    endif
#endif
    line%r1=point1
    line%r2=point2
    line%length=norm2(point1-point2)
end subroutine

!> Create a polygon from a set of points in the same plane
module subroutine polygon_from_points(polygon,points,plane)
    !> The resulting polygon
    class(lo_polygon), intent(out) :: polygon
    !> The points on the plane. No check that they are on the plane
    real(flyt), dimension(:,:), intent(in) :: points
    !> The plane that the points lie on.
    type(lo_plane), intent(in) :: plane
    !
    integer :: i

#ifdef AGRESSIVE_SANITY
    do i=1,size(points,2)
        if ( abs(plane%distance_to_point(points(:,i))) .gt. lo_tol ) then
            call lo_stop_gracefully(['Trying to create a polygon from planar points, but they are not on the same plane.'],lo_exitcode_param)
        endif
    enddo
    if ( size(points,2) .lt. 3 ) then
        call lo_stop_gracefully(['Trying to create a polygon from planar points, but only '//tochar(size(points,2))//' points provided'],lo_exitcode_param)
    endif
#endif

    polygon%n=size(points,2)
    allocate(polygon%r(3,polygon%n))
    polygon%r=points
    call plane%anglesort(polygon%r)
    polygon%plane=plane
    do i=1,3
        polygon%centroid(i)=lo_mean(polygon%r(i,:))
    enddo
    polygon%rmax=0.0_flyt
    do i=1,polygon%n
        polygon%rmax=max(polygon%rmax,norm2(polygon%centroid-polygon%r(:,i)))
    enddo
end subroutine

end submodule
