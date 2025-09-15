!#include "precompilerdefinitions"
module geometryfunctions
!! Basic geometric functions and derived types for lines, planes, polygons and convex polyhedrons. Routines to calculate distances, areas, volumes, as well as more advanced ones such as slicing polyhedrons with planes.
use konstanter, only: r8, lo_huge, lo_hugeint, lo_status, lo_tol, lo_sqtol, lo_pi, lo_exitcode_param
use gottochblandat, only: walltime,tochar,qsort,lo_cross,lo_mean,lo_sqnorm,lo_index_in_periodic_array, lo_stop_gracefully, lo_return_unique
implicit none
private
! the public derived types
public :: lo_geometryobject
public :: lo_null
public :: lo_point
public :: lo_line
public :: lo_linesegment
public :: lo_plane
public :: lo_polygon
public :: lo_polyhedron
! the public functions
public :: lo_angle_between_vectors
public :: lo_rotation_matrix_from_axis_and_angle
public :: lo_improper_rotation_matrix_from_axis_and_angle
public :: lo_rotation_matrix_from_vector_a_to_b
public :: lo_inscribed_sphere_in_box
public :: lo_bounding_sphere_of_box
public :: lo_convex_hull_2d
public :: lo_increment_dimensions

!> General type for geometrical objects in R3.
type, abstract :: lo_geometryobject
    contains
        !> Move the object
        procedure :: translate => translate_geometrical_object
        !> Slice the object with a plane
        procedure :: slice => slice_geometrical_object
end type lo_geometryobject

!> The null geometrical object, nothing
type, extends(lo_geometryobject) :: lo_null
end type lo_null

!> A point
type, extends(lo_geometryobject) :: lo_point
    !> Point, believe it or not.
    real(r8), dimension(3) :: r=lo_huge
end type lo_point

!> An infinite line
type, extends(lo_geometryobject) :: lo_line
    !> First point
    real(r8), dimension(3) :: r1=lo_huge
    !> Second point
    real(r8), dimension(3) :: r2=lo_huge
    contains
        !> generate line from a set of points
        procedure :: generate => create_line_from_points
end type lo_line

!> A line segment. Exactly the same as a line, only differs when intersections and similar are calculated.
type, extends(lo_geometryobject) :: lo_linesegment
    !> First point
    real(r8), dimension(3) :: r1=lo_huge
    !> Second point
    real(r8), dimension(3) :: r2=lo_huge
    !> the length of the segment
    real(r8) :: length=lo_huge
    contains
        !> The distance to a point
        procedure :: distance_to_point => distance_linesegment_point
        !> generate linesegment from a set of points
        procedure :: generate => create_linesegment_from_points
end type

!> A plane on Hessian normal form, that is $$ \hat{\mathbf{n}} \cdot \mathbf{x} + p = 0 $$ Also contains contains an orthornmal coordinate system aligned with the plane, useful for projections.
type, extends(lo_geometryobject) :: lo_plane
    !> The plane normal
    real(r8), dimension(3) :: normal=lo_huge
    !> The distance from the origin
    real(r8) :: p=lo_huge
    !> first vector that spans the plane. Orthogonal to the normal and v2.
    real(r8), dimension(3) :: v1=lo_huge
    !> second vector that spans the plane. Orthogonal to the normal and v1.
    real(r8), dimension(3) :: v2=lo_huge
    contains
        !> create the plane from a set of points.
        procedure :: generate => create_plane_from_points
        !> distance to a point
        procedure :: distance_to_point => lo_plane_point_distance
        !> sort a set of points clockwise on this plane. Or counterclockwise, not sure. In order at least.
        procedure :: anglesort => anglesort_with_plane
        !> get the reflection matrix
        procedure :: reflection_matrix
end type

!> A convex polygon in R3
type, extends(lo_geometryobject) :: lo_polygon
    !> How many points on this polygon?
    integer :: n=-lo_hugeint
    !> The points
    real(r8), dimension(:,:), allocatable :: r
    !> The centroid
    real(r8), dimension(3) :: centroid=lo_huge
    !> Radius of bounding sphere
    real(r8) :: rmax=lo_huge
    !> The plane the points are confined to
    type(lo_plane) :: plane
    contains
        !> Create polygon from a set of points
        procedure :: generate => polygon_from_points
        !> Calculate the area
        procedure :: area => area_of_polygon
end type

!> A convex polyhedron in R3, defined as a list of polygons.
type, extends(lo_geometryobject) :: lo_polyhedron
    !> How many faces on this polygon?
    integer :: n=-lo_hugeint
    !> The faces
    type(lo_polygon), dimension(:), allocatable :: face
    contains
        !> create a simple big box
        procedure :: createcube
        !> Calculate the volume
        procedure :: volume => volume_of_polyhedron
        !> check if a point is inside
        procedure :: is_point_inside => is_point_inside_polyhedron
end type

! Interface to slicing routines
interface
    module subroutine slice_geometrical_object(g,plane,aboveplane,onplane,belowplane,modifications,verbosity)
        class(lo_geometryobject), intent(in) :: g
        type(lo_plane), intent(in) :: plane
        class(lo_geometryobject), intent(out), allocatable :: aboveplane
        class(lo_geometryobject), intent(out), allocatable :: belowplane
        class(lo_geometryobject), intent(out), allocatable :: onplane
        logical, intent(out) :: modifications
        integer, intent(in), optional :: verbosity
    end subroutine
end interface

! Interface to generating routines
interface
    module subroutine createcube(polyhedron,side)
        class(lo_polyhedron), intent(inout) :: polyhedron
        real(r8), intent(in) :: side
    end subroutine
    module subroutine create_plane_from_points(plane,points,point,normal)
        class(lo_plane), intent(out) :: plane
        real(r8), dimension(3,3), intent(in), optional :: points
        real(r8), dimension(3), intent(in), optional :: point
        real(r8), dimension(3), intent(in), optional :: normal
    end subroutine
    module subroutine create_line_from_points(line,point1,point2)
        class(lo_line), intent(out) :: line
        real(r8), dimension(3), intent(in) :: point1
        real(r8), dimension(3), intent(in) :: point2
    end subroutine
    module subroutine create_linesegment_from_points(line,point1,point2)
        class(lo_linesegment), intent(out) :: line
        real(r8), dimension(3), intent(in) :: point1
        real(r8), dimension(3), intent(in) :: point2
    end subroutine
    module subroutine polygon_from_points(polygon,points,plane)
        class(lo_polygon), intent(out) :: polygon
        real(r8), dimension(:,:), intent(in) :: points
        type(lo_plane), intent(in) :: plane
    end subroutine
end interface

! Interfaces to type-bound things
interface
    module pure function distance_linesegment_point(linesegment,point) result(d)
        class(lo_linesegment), intent(in) :: linesegment
        real(r8), dimension(3), intent(in) :: point
        real(r8) :: d
    end function
    module recursive subroutine translate_geometrical_object(g,translation)
        class(lo_geometryobject), intent(inout) :: g
        real(r8), dimension(3), intent(in) :: translation
    end subroutine
    module pure function lo_plane_point_distance(plane,point) result(r)
        class(lo_plane), intent(in) :: plane
        real(r8), dimension(3), intent(in) :: point
        real(r8) :: r
    end function
    module subroutine anglesort_with_plane(plane,points,order)
        class(lo_plane), intent(in) :: plane
        real(r8), dimension(:,:), intent(inout) :: points
        integer, dimension(:), intent(out), optional :: order
    end subroutine
    module pure function reflection_matrix(plane) result(m)
        class(lo_plane), intent(in) :: plane
        real(r8), dimension(3,3) :: m
    end function
    module pure function area_of_polygon(polygon) result(A)
        class(lo_polygon), intent(in) :: polygon
        real(r8) :: A
    end function
    module pure function volume_of_polyhedron(polyhedron) result(V)
        class(lo_polyhedron), intent(in) :: polyhedron
        real(r8) :: V
    end function
    module pure function is_point_inside_polyhedron(polyhedron,point,tol) result(inside)
        class(lo_polyhedron), intent(in) :: polyhedron
        real(r8), dimension(3), intent(in) :: point
        real(r8), intent(in) :: tol
        logical :: inside
    end function
end interface

contains

!> Angle between two vectors. Returns the (unsigned) angle between the vectors in radians.
pure function lo_angle_between_vectors(a,b) result(angle)
    !> first vector
    real(r8), dimension(3), intent(in) :: a
    !> second vector
    real(r8), dimension(3), intent(in) :: b
    !> the angle
    real(r8) :: angle
    angle=acos(dot_product(a,b)/(norm2(a)*norm2(b)))
end function

!> Returns a rotation matrix. Defined as rotation about an axis u with angle alpha. Returns 3x3 rotation matrix. Use as `matmul(rotationmatrix,vector)`.
pure function lo_rotation_matrix_from_axis_and_angle(u,alpha) result(m)
    !> the axis
    real(r8), dimension(3), intent(in) :: u
    !> the angle, in radians
    real(r8), intent(in) :: alpha
    !> the rotation matrix
    real(r8), dimension(3,3) :: m
    !
    real(r8) :: st,ct,a,b,c,invnrm

    invnrm=norm2(u)
    ! If the vector is really tiny, I can get strange stuff.  makes little sense if the vector is supersmall
    ! not sure what to do there.
    if ( invnrm .gt. lo_tol ) then
        invnrm=1.0_r8/invnrm
        a=u(1)*invnrm
        b=u(2)*invnrm
        c=u(3)*invnrm
    else
        a=0.0_r8
        b=0.0_r8
        c=0.0_r8
    endif
    !
    st=sin(alpha)
    ct=cos(alpha)
    ! Apparently it's called Rodriguez formula.
    m(1,1)=ct+a*a*(1.0_r8-ct)
    m(1,2)=a*b*(1.0_r8-ct)-c*st
    m(1,3)=a*c*(1.0_r8-ct)+b*st
    !
    m(2,1)=b*a*(1.0_r8-ct)+c*st
    m(2,2)=ct+b*b*(1.0_r8-ct)
    m(2,3)=b*c*(1.0_r8-ct)-a*st
    !
    m(3,1)=a*c*(1.0_r8-ct)-b*st
    m(3,2)=c*b*(1.0_r8-ct)+a*st
    m(3,3)=ct+c*c*(1.0_r8-ct)
end function

!> Returns an improper rotation matrix. Defined as rotation about an axis u with angle alpha, followed by reflection through the origin. Use as `matmul(rotationmatrix,vector)`.
pure function lo_improper_rotation_matrix_from_axis_and_angle(u,alpha) result(m)
    !> the axis
    real(r8), dimension(3), intent(in) :: u
    !> the angle, in radians
    real(r8), intent(in) :: alpha
    !> the rotation matrix
    real(r8), dimension(3,3) :: m
    !
    real(r8) :: st,ct,a,b,c,invnrm

    invnrm=norm2(u)
    ! If the vector is really tiny, I can get strange stuff.  makes little sense if the vector is supersmall
    ! not sure what to do there. Now I return NaN perhaps?
    if ( invnrm .gt. lo_tol ) then
        invnrm=1.0_r8/invnrm
        a=u(1)*invnrm
        b=u(2)*invnrm
        c=u(3)*invnrm
    else
        a=0.0_r8
        b=0.0_r8
        c=0.0_r8
    endif
    !
    st=sin(alpha)
    ct=cos(alpha)
    ! quarternion thing from wikipedia
    m(1,1)=ct-a*a*(1.0_r8+ct)
    m(1,2)=-a*b*(1.0_r8+ct)-c*st
    m(1,3)=-a*c*(1.0_r8+ct)+b*st
    !
    m(2,1)=-b*a*(1.0_r8+ct)+c*st
    m(2,2)=ct-b*b*(1.0_r8+ct)
    m(2,3)=-b*c*(1.0_r8+ct)-a*st
    !
    m(3,1)=-a*c*(1.0_r8+ct)-b*st
    m(3,2)=-c*b*(1.0_r8+ct)+a*st
    m(3,3)=ct-c*c*(1.0_r8+ct)
end function

!> Rotation matrix that takes vector a to vector b. Returns a 3x3 rotation matrix, use as `matmul(rotmat,vector)`. The rotation matrix satisifes rotmat*a-b=0 in the special case when |a|=|b|, otherwise rotmat*a || b.
pure function lo_rotation_matrix_from_vector_a_to_b(a,b) result(m)
    !> vector a
    real(r8), dimension(3), intent(in) :: a
    !> vector b
    real(r8), dimension(3), intent(in) :: b
    !> the rotation matrix
    real(r8), dimension(3,3) :: m
    !
    real(r8), dimension(3) :: axis,v0,v1
    real(r8) :: angle
    integer :: i

    ! Get the rotation axis
    v0=a/norm2(a)
    v1=b/norm2(b)
    axis=lo_cross(v0,v1)
    ! sanity check
    if ( norm2(axis) .lt. lo_sqtol ) then
        ! the rotation is too small, just return the identity matrix
        m=0.0_r8
        do i=1,3
            m(i,i)=1.0_r8
        enddo
    else ! get the actual rotation matrix
        axis=axis/norm2(axis)
        ! Get the angle
        angle=lo_angle_between_vectors(a,b)
        ! And the rotation matrix
        m=lo_rotation_matrix_from_axis_and_angle(axis,angle)
    endif
end function

!> given a parallelepiped, what is the radius of the largest sphere that fits inside?
pure function lo_inscribed_sphere_in_box(box) result(r)
    !> vectors that define the box
    real(r8), dimension(3,3), intent(in) :: box
    !> largest sphere that fit inside box
    real(r8) :: r

    real(r8), dimension(3) :: na,nb,nc
    ! the normals of the faces of the box
    na=lo_cross(box(:,2),box(:,3))
    nb=lo_cross(box(:,3),box(:,1))
    nc=lo_cross(box(:,1),box(:,2))

    na=na/norm2(na)
    nb=nb/norm2(nb)
    nc=nc/norm2(nc)
    ! distances between opposing planes
    r=lo_huge
    r=min(r,abs(dot_product(na,box(:,1))))
    r=min(r,abs(dot_product(nb,box(:,2))))
    r=min(r,abs(dot_product(nc,box(:,3))))
    r=r*0.5_r8
end function

!> given a parallelepiped, what is the radius of the smallest sphere that fits the box?
pure function lo_bounding_sphere_of_box(box) result(r)
    !> vectors that define the box
    real(r8), dimension(3,3), intent(in) :: box
    !> radius of smallest sphere that contains the box
    real(r8) :: r

    real(r8), dimension(3) :: a,b,c,v
    real(r8), dimension(4) :: diags

    a=box(:,1)
    b=box(:,2)
    c=box(:,3)

    ! had to do this in an awkward way to avoid segfaults on intel compilers. Odd.
    v=a-b-c
    diags(1)=dot_product(v,v)
    v=a+b-c
    diags(2)=dot_product(v,v)
    v=a-b+c
    diags(3)=dot_product(v,v)
    v=a+b+c
    diags(4)=dot_product(v,v)
    r=maxval(diags)
    r=sqrt(r)*0.5_r8
end function

!> Calculates the convex hull for a set of 2d points. The algorithm is a Graham scan @cite Graham1972a and checks for counterclockwise turns within a tolerance. The output is an index array to the points on the hull, and optionally an additional array with the points inside the hull.
subroutine lo_convex_hull_2d(points,tol,points_on_hull,points_in_hull)
    !> the points
    real(r8), dimension(:,:), intent(in) :: points
    !> the tolerance
    real(r8), intent(in) :: tol
    !> the points on the hull
    integer, dimension(:), allocatable, intent(out) :: points_on_hull
    !> points inside the hull
    integer, dimension(:), allocatable, intent(out), optional :: points_in_hull

    real(r8), dimension(:,:), allocatable :: dum
    real(r8), dimension(:), allocatable :: d1
    real(r8), dimension(2) :: p0
    integer, dimension(:), allocatable :: ind_y,ind_angle,ind_tot,ind_on
    integer :: i,k,l,np,ii

    np=size(points,2)
    ! maybe we can stop early
    if ( np .le. 3 ) then
        ! The hull is trivial if there are three or less points.
        allocate(points_on_hull(np))
        do i=1,np
            points_on_hull(i)=i
        enddo
        return
    endif
    ! Sort points according to y-coordinate
    allocate(d1(np))
    allocate(ind_y(np))
    d1=points(2,:)
    call qsort(d1,ind_y)
    ! calculate angle to first point
    allocate(dum(2,np))
    dum=points(:,ind_y)
    p0=dum(:,1)
    do i=1,np
        d1(i)=atan2(dum(2,i)-p0(2),dum(1,i)-p0(1))
    enddo

    ! sort according do this angle.
    allocate(ind_angle(np))
    call qsort(d1,ind_angle)
    ! get the total sorting-index-mapping-thing
    allocate(ind_tot(np))
    ind_tot=ind_y(ind_angle)
    ! clean a bit
    deallocate(ind_y)
    deallocate(ind_angle)

    ! start the actual algorithm
    allocate(ind_on(np))
    dum=0.0_r8
    dum(:,1)=points(:,ind_tot(1))
    ind_on(1)=ind_tot(1)
    l=1
    ii=1
    do i=2,np
        ii=ii+1
        p0=points(:,ind_tot(i))
        do
            ! exit if only one point on the hull
            if ( l .lt. 2 ) exit
            ! check if previous point should be removed
            if ( cross2d(dum(:,l-1),dum(:,l),p0,tol) .le. 0 ) then
                l=l-1
                cycle
            else
                exit
            endif
        enddo
        ! this point might be on the hull
        l=l+1
        dum(:,l)=p0
        ind_on(l)=ind_tot(i)
    enddo

    ! Return the hull
    allocate(points_on_hull(l))
    points_on_hull=ind_on(1:l)
    if ( present(points_in_hull) ) then
        ! Returning the points inside the hull.
        k=np-l
        allocate(points_in_hull(k))
        ! reusing this dummy array to create a flat list of indices
        do i=1,np
            ind_tot(i)=i
        enddo
        ! set the ones that are on the hull to -1
        do i=1,l
            ind_tot(ind_on(i))=-1
        enddo
        ! return those that are not -1, that's the points inside the hull.
        l=0
        do i=1,np
            if ( ind_tot(i) .gt. 0 ) then
                l=l+1
                points_in_hull(l)=ind_tot(i)
            endif
        enddo
    endif
    ! cleanup
    deallocate(ind_on)
    deallocate(ind_tot)
    deallocate(dum)

    contains
    !> Check if three points make a clockwise or counterclockwise turn
    function cross2d(p1,p2,p3,tolerance) result(r)
        !> first point
        real(r8), dimension(2), intent(in) :: p1
        !> second point
        real(r8), dimension(2), intent(in) :: p2
        !> third point
        real(r8), dimension(2), intent(in) :: p3
        !> tolerance for what is 0
        real(r8), intent(in) :: tolerance
        !> A value < 0 means it's counterclockwise, > 0 clockwise and 0 that the points are on the same line.
        integer :: r

        real(r8), dimension(2) :: v1,v2
        real(r8) :: alpha

        v1=p2-p1
        v2=p2-p3
        ! signed angle between vectors
        alpha=( atan2(v2(2),v2(1))-atan2(v1(2),v1(1)) )*180.0_r8/lo_pi
        alpha=mod(alpha+360.0_r8,360.0_r8)
        if ( abs(alpha-180.0_r8) .lt. tolerance .or. abs(alpha) .lt. tolerance ) then
            r=0
        elseif ( alpha .gt. 180_r8 ) then
            r=1
        else
            r=-1
        endif
    end function
end subroutine

function lo_increment_dimensions(dimin, box) result(dimut)
    integer, dimension(3), intent(in) :: dimin
    real(r8), dimension(3, 3), intent(in) :: box
    integer, dimension(3) :: dimut
    !
    real(r8), dimension(3, 3) :: m0
    integer, dimension(3) :: di
    integer :: i, j, k
    real(r8) :: f0, f1, ff0

    ! try with the sphere thing. First get a baseline
    do j = 1, 3
        m0(:, j) = box(:, j)*(2*dimin(j) + 1)
    end do
    ff0 = lo_inscribed_sphere_in_box(m0)
    ! Increment the dimension that gives the biggest increase in radii
    f0 = 0.0_r8
    dimut = 0
    do i = 1, 3
        di = dimin
        di(i) = di(i) + 1
        do j = 1, 3
            m0(:, j) = box(:, j)*(2*di(j) + 1)
        end do
        f1 = lo_inscribed_sphere_in_box(m0)
        if (f1 .gt. f0 .and. abs(f1 - ff0) .gt. lo_tol) then
            dimut = di
            f0 = f1
        end if
    end do

    ! if nothing helped, increment the lowest number
    if (dimut(1) .eq. 0) then
        j = lo_hugeint
        do i = 1, 3
            if (di(i) .lt. j) then
                j = di(i)
                k = i
            end if
        end do
        dimut = dimin
        dimut(k) = dimut(k) + 1
    end if
end function

end module
