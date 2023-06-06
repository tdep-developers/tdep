#include "precompilerdefinitions"
module cgal_wrappers
!! Wrappers to the c++ library CGAL. These are pretty low-level, higher level routines are elsewhere.
use konstanter, only: flyt,lo_sqtol,lo_status,lo_exitcode_baddim
use gottochblandat, only: lo_stop_gracefully,tochar
use iso_c_binding

implicit none
private

! massive precompiler remove thingy since they Cray machines got really angry otherwise
#ifdef usecgal
public :: lo_chull_2d
public :: lo_chull_3d
public :: lo_chull_3d_intersection
public :: lo_deltri_2d
public :: lo_deltri_3d
public :: lo_tesselate_polyhedron

!> Interfaces to CGAL routines
interface

!> interface to 2D convex hull
subroutine cgal_chull2d(npts,r_in,npts_on_hull,r_hull) bind(C,name="lo_cgal_chull2")
    import ! no idea what this does, but I think it has to be there.
    !> how many points
    integer(c_int), intent(in) :: npts
    !> the points for which to calculate the hull
    real(c_double), dimension(2,npts), intent(in) :: r_in
    !> how many points are on the hull?
    integer(c_int), intent(out) :: npts_on_hull
    !> the points on the hull
    type(c_ptr), intent(out) :: r_hull
end subroutine

!> interface to 3D convex hull
subroutine cgal_chull3d(npts,r_in,npts_on_hull,r_hull) bind(C,name="lo_cgal_chull3")
    import
    !> how many points
    integer(c_int), intent(in) :: npts
    !> the points for which to calculate the hull
    real(c_double), dimension(3,npts), intent(in) :: r_in
    !> how many points are on the hull?
    integer(c_int), intent(out) :: npts_on_hull
    !> the points on the hull
    type(c_ptr), intent(out) :: r_hull
end subroutine

!> interface to intersection of two 3D convex hulls
subroutine cgal_chull3d_intersection(npts1,npts2,r_in1,r_in2,npts_on_hull,r_hull,volume) bind(C,name="lo_cgal_chull3_intersection")
    import
    !> how many points
    integer(c_int), intent(in) :: npts1,npts2
    !> the points that describe each hull
    real(c_double), dimension(3,npts1), intent(in) :: r_in1
    real(c_double), dimension(3,npts2), intent(in) :: r_in2
    !> how many points are intersected hull?
    integer(c_int), intent(out) :: npts_on_hull
    !> the points on the intersected
    type(c_ptr), intent(out) :: r_hull
    !> the intersected volume
    real(c_double), intent(out) :: volume
end subroutine

!> interface to 2D Delaunay Triangulation
subroutine cgal_deltri2d(npts,r_in,ntri,tri,triarea) bind(C,name="lo_cgal_deltri2")
    import
    !> how many points
    integer(c_int), intent(in) :: npts
    !> the points for which to calculate the hull
    real(c_double), dimension(2,npts), intent(in) :: r_in
    !> how many triangles did I get?
    integer(c_int), intent(out) :: ntri
    !> the triangles represented as indices to the original points
    type(c_ptr), intent(out) :: tri
    !> the area of each triangle for filtering
    type(c_ptr), intent(out) :: triarea
end subroutine

!> interface to 3D Delaunay Triangulation
subroutine cgal_deltri3d(npts,r_in,ntet,tet,tetvol) bind(C,name="lo_cgal_deltri3")
    import
    !> how many points
    integer(c_int), intent(in) :: npts
    !> the points for which to calculate the hull
    real(c_double), dimension(3,npts), intent(in) :: r_in
    !> how many triangles did I get?
    integer(c_int), intent(out) :: ntet
    !> the tetrahedrons represented as indices to the original points
    type(c_ptr), intent(out) :: tet
    !> the volume of each tetrahedron
    type(c_ptr), intent(out) :: tetvol
end subroutine

!> interface to tesselate a convex polyhedron
subroutine cgal_tesselate_polyhedron(npts,r_in,nvert,vertices,ntet,tetrahedrons,&
           max_dihedral_angle,edge_size,facet_angle,facet_size,facet_distance,cell_radius_edge_ratio,cell_size) bind(C,name="lo_cgal_tesselate_polyhedron")
    import
    !> how many points
    integer(c_int), intent(in) :: npts
    !> the nodes of the polyhedron
    real(c_double), dimension(3,npts), intent(in) :: r_in
    !> how many vertices did I get?
    integer(c_int), intent(out) :: nvert
    !> the actual vertices
    type(c_ptr), intent(out) :: vertices
    !> how many tetrahedrons did I get?
    integer(c_int), intent(out) :: ntet
    !> the actual tetrahedrons
    type(c_ptr), intent(out) :: tetrahedrons
    !> the parameters for the meshing
    real(c_double), intent(in) :: max_dihedral_angle
    real(c_double), intent(in) :: edge_size
    real(c_double), intent(in) :: facet_angle
    real(c_double), intent(in) :: facet_size
    real(c_double), intent(in) :: facet_distance
    real(c_double), intent(in) :: cell_radius_edge_ratio
    real(c_double), intent(in) :: cell_size
end subroutine

!> deallocate a c++ int pointer.
subroutine cgal_delete_int_pointer(A) bind(C,name="lo_cgal_cleanup_int_pointer")
    import
    !> the pointer
    type(c_ptr), intent(in) :: A
end subroutine

!> deallocate a c++ double pointer.
subroutine cgal_delete_double_pointer(A) bind(C,name="lo_cgal_cleanup_double_pointer")
    import
    !> the pointer
    type(c_ptr), intent(in) :: A
end subroutine

end interface

contains


!> Wrapper to calculate the 2D convex hull
subroutine lo_chull_2d(r,rhull)
    real(flyt), dimension(:,:), intent(in) :: r
    real(flyt), dimension(:,:), allocatable, intent(out) :: rhull
    !
    integer(c_int) :: cnpts,cnpts_hull
    type(c_ptr) :: ptrhull
    real(c_double), dimension(:,:), allocatable :: cr
    real(c_double), dimension(:,:), pointer :: frhull

    ! number of points
    cnpts=size(r,2)
#ifdef AGRESSIVE_SANITY
    if ( cnpts .le. 2 ) then
        call lo_stop_gracefully(['Makes no sense to calculate the convex hull for '//tochar(size(r,2))//' or less points'],&
            lo_exitcode_baddim)
    endif
#endif
    ! the arrays passed to cgal needs to be pre-allocated
    lo_allocate(cr(2,cnpts))
    cr=r
    ! call CGAL routine via a C wrapper of c++
    call cgal_chull2d(cnpts,cr,cnpts_hull,ptrhull)
    ! convert from C to fortran
    call c_f_pointer(ptrhull,frhull,[2,cnpts_hull])
    ! make it a nice allocatable array
    lo_allocate(rhull(2,cnpts_hull))
    rhull=frhull
    !
    ! And some cleanup
    !
    call cgal_delete_double_pointer(ptrhull)
    lo_deallocate(cr)
end subroutine

!> Wrapper to calculate the 3D convex hull
subroutine lo_chull_3d(r,rhull)
    real(flyt), dimension(:,:), intent(in) :: r
    real(flyt), dimension(:,:), allocatable, intent(out), optional :: rhull

    integer(c_int) :: cnpts,cnpts_hull
    type(c_ptr) :: ptrhull
    real(c_double), dimension(:,:), allocatable :: cr
    real(c_double), dimension(:,:), pointer :: frhull

    ! number of points
    cnpts=size(r,2)
#ifdef AGRESSIVE_SANITY
    if ( cnpts .le. 2 ) then
        call lo_stop_gracefully(['Makes no sense to calculate the convex hull for '//tochar(size(r,2))//' or less points'],&
            lo_exitcode_baddim)
    endif
#endif
    ! the arrays passed to cgal needs to be pre-allocated
    lo_allocate(cr(3,cnpts))
    cr=r
    ! call CGAL routine via a C wrapper of c++
    call cgal_chull3d(cnpts,cr,cnpts_hull,ptrhull)
    ! convert from C to fortran
    call c_f_pointer(ptrhull,frhull,[3,cnpts_hull])
    ! make a neat allocatable array
    lo_allocate(rhull(3,cnpts_hull))
    rhull=frhull

    ! And some cleanup
    call cgal_delete_double_pointer(ptrhull)
    lo_deallocate(cr)
end subroutine

!> Wrapper to calculate intersection of two convex hulls
subroutine lo_chull_3d_intersection(r1,r2,rhull,volume)
    real(flyt), dimension(:,:), intent(in) :: r1,r2
    real(flyt), dimension(:,:), allocatable, intent(out), optional :: rhull
    real(flyt), intent(out), optional :: volume
    !
    integer(c_int) :: cnpts1,cnpts2,cnpts_hull
    type(c_ptr) :: ptrhull
    real(c_double), dimension(:,:), allocatable :: cr1,cr2
    real(c_double), dimension(:,:), pointer :: frhull
    real(c_double) :: cvol

    ! number of points
    cnpts1=size(r1,2)
    cnpts2=size(r2,2)
    ! the arrays passed to cgal needs to be pre-allocated
    lo_allocate(cr1(3,cnpts1))
    lo_allocate(cr2(3,cnpts2))
    cr1=r1
    cr2=r2
    ! call CGAL routine via a C wrapper of c++
    call cgal_chull3d_intersection(cnpts1,cnpts2,cr1,cr2,cnpts_hull,ptrhull,cvol)
    ! make sure there is at least one point. Or two. Whatever.
    if ( cnpts_hull .gt. 0 ) then
        ! convert from C to fortran
        call c_f_pointer(ptrhull,frhull,[3,cnpts_hull])
        ! make it a nice allocatable array
        if ( present(rhull) ) then
            lo_allocate(rhull(3,cnpts_hull))
            rhull=frhull
        endif
        !
        if ( present(volume) ) then
            volume=cvol
        endif
        ! And some cleanup
        call cgal_delete_double_pointer(ptrhull)
    else
        ! set the volume to zero
        if ( present(volume) ) then
            volume=0.0_flyt
        endif
        ! not sure what to do about rhull. I think I leave it unallocated.
    endif
    ! a little more cleanup
    lo_deallocate(cr1)
    lo_deallocate(cr2)
end subroutine

!> Wrapper to calculate the 2D delaunay triangulation
subroutine lo_deltri_2d(r,tri,tolerance,area)
    !> the points to triangulate
    real(flyt), dimension(:,:), intent(in) :: r
    !> the triangles as quartets of indicies to the original points
    integer, dimension(:,:), allocatable, intent(out) :: tri
    !> a tolerance to skip slivers
    real(flyt), intent(in), optional :: tolerance
    !> the total area of the triangulation
    real(flyt), intent(out), optional :: area

    ! for the c++ part
    integer(c_int) :: cnpts,cntri
    type(c_ptr) :: ptrtri,ptrtriarea
    real(c_double), dimension(:,:), allocatable :: cr
    integer(c_int), dimension(:,:), pointer :: ctri
    real(c_double), dimension(:), pointer :: ctriarea
    ! for the normal part
    integer :: i,j
    real(flyt) :: f0,tol
    integer, dimension(:), allocatable :: ind

    ! Set the tolerance
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_sqtol
    endif
    ! number of points
    cnpts=size(r,2)
    if ( cnpts .lt. 3 ) then
        write(*,*) 'cannot triangulate with less than 3 points'
        stop
    endif
    ! the arrays passed to cgal needs to be pre-allocated, both in and out.
    lo_allocate(cr(2,cnpts))
    ! according to some dude on the internet, this should be a safe upper bound to the number of triangles.
    cr=r
    ! call CGAL routine via a C wrapper of c++
    call cgal_deltri2d(cnpts,cr,cntri,ptrtri,ptrtriarea)
    ! convert from C to fortran
    call c_f_pointer(ptrtri,ctri,[3,cntri])
    call c_f_pointer(ptrtriarea,ctriarea,[cntri])
    ! maybe store the area of the convex hull
    if ( present(area) ) then
        area=sum(ctriarea)
    endif

    ! Ok, now I have all the triangles. Skip the stupid ones.
    f0=sum(ctriarea)*cnpts
    lo_allocate(ind(cntri))
    ind=1
    do i=1,cntri
        if ( ctriarea(i)/f0 .lt. tol ) ind(i)=0
    enddo

    ! store the good ones
    lo_allocate(tri(3,sum(ind)))
    j=0
    do i=1,cntri
        if ( ind(i) .eq. 1 ) then
            j=j+1
            tri(:,j)=ctri(:,i)
        endif
    enddo

    ! And some cleanup
    call cgal_delete_int_pointer(ptrtri)
    call cgal_delete_double_pointer(ptrtriarea)
    lo_deallocate(cr)
    lo_deallocate(ind)
end subroutine

!> Wrapper to calculate the 3D delaunay triangulation
subroutine lo_deltri_3d(r,tet,tolerance,volume)
    !> points to triangulate
    real(flyt), dimension(:,:), intent(in) :: r
    !> tetrahedrons as quartets of indicies to the original points
    integer, dimension(:,:), allocatable, intent(out) :: tet
    !> a tolerance to skip slivers
    real(flyt), intent(in), optional :: tolerance
    !> the volume of the convex hull
    real(flyt), intent(out), optional :: volume

    ! for the c++ part
    integer(c_int) :: cnpts,cntet
    type(c_ptr) :: ptrtet,ptrtetvol
    real(c_double), dimension(:,:), allocatable :: cr
    integer(c_int), dimension(:,:), pointer :: ctet
    real(c_double), dimension(:), pointer :: ctetvol
    ! for the normal part
    integer :: i,j
    real(flyt) :: f0,tol
    integer, dimension(:), allocatable :: ind

    ! Set the tolerance
    if ( present(tolerance) ) then
        tol=tolerance
    else
        tol=lo_sqtol
    endif

    ! Recast the points into C-form, to make CGAL happy
    cnpts=size(r,2)
    if ( cnpts .lt. 4 ) then
        write(*,*) 'cannot triangulate with less than 4 points'
        stop
    endif
    lo_allocate(cr(3,cnpts))
    cr=r

    ! call CGAL routine via a C wrapper of c++
    call cgal_deltri3d(cnpts,cr,cntet,ptrtet,ptrtetvol)
    ! convert from C to fortran
    call c_f_pointer(ptrtet,ctet,[4,cntet])
    call c_f_pointer(ptrtetvol,ctetvol,[cntet])
    ! maybe store the volume
    if ( present(volume) ) then
        volume=sum(ctetvol)
    endif

    ! Ok, now I have all the tetrahedrons. There might be a bunch of them that are stupid. I want to skip
    ! those. First get the average volume/tetrahedron
    f0=tol*sum(ctetvol)/cntet
    lo_allocate(ind(cntet))
    ind=1
    do i=1,cntet
        if ( ctetvol(i) .lt. f0 ) ind(i)=0
    enddo

    ! store the good ones
    lo_allocate(tet(4,sum(ind)))
    j=0
    do i=1,cntet
        if ( ind(i) .eq. 1 ) then
            j=j+1
            tet(:,j)=ctet(:,i)
        endif
    enddo

    ! And some cleanup
    call cgal_delete_int_pointer(ptrtet)
    call cgal_delete_double_pointer(ptrtetvol)
    lo_deallocate(cr)
    lo_deallocate(ind)
end subroutine

!> create a mesh from convex polyhedron
subroutine lo_tesselate_polyhedron(r,meshpars,vertices,tetrahedrons)
    !> the points that define a convex polyhedron
    real(flyt), dimension(:,:), intent(in) :: r
    !> settings for the tesselation, see the CGAL manual
    real(flyt), dimension(7), intent(in) :: meshpars
    !> resulting points
    real(flyt), dimension(:,:), allocatable, intent(out) :: vertices
    !> resulting tetrahedrons
    integer, dimension(:,:), allocatable, intent(out) :: tetrahedrons

    ! the c-stuff
    integer(c_int) :: cnpts,cnvert,cntet
    real(c_double), dimension(:,:), allocatable :: cr
    type(c_ptr) :: cvertices,ctetrahedrons
    real(c_double) :: max_dihedral_angle,edge_size,facet_angle,facet_size,facet_distance,cell_radius_edge_ratio,cell_size
    ! when it gets pointed back to fortran
    real(c_double), dimension(:,:), pointer :: fvert
    integer(c_int), dimension(:,:), pointer :: ftets
    ! and the actual arrays I want.
    integer :: ntet,npts
    !
    ! I hope that my irreducible wedge is equal to the convex hull of my nodes.
    ! If this breaks, I have to triangulate the surfaces, which should be trivial.
    ! I commented out a snippet of code in the c++ routine that takes a triangulated
    ! convex polyhedron in my format and gets it to cgal format.
    !
    ! But let's hope it works anyway! Start by stuffing points into c format
    cnpts=size(r,2)
    lo_allocate(cr(3,cnpts))
    cr=r
    ! The settings for the tesselation. Should probably read the CGAL manual more carefully.
    max_dihedral_angle=     meshpars(1) ! 120
    edge_size=              meshpars(2) ! 0.025
    facet_angle=            meshpars(3) ! 20
    facet_size=             meshpars(4) ! 0.03
    facet_distance=         meshpars(5) ! 0.05
    cell_radius_edge_ratio= meshpars(6) ! 3
    cell_size=              meshpars(7) ! 0.05

    ! Run the cgal routine
    call cgal_tesselate_polyhedron(cnpts,cr,cnvert,cvertices,cntet,ctetrahedrons,&
                                   max_dihedral_angle,&
                                   edge_size,&
                                   facet_angle,&
                                   facet_size,&
                                   facet_distance,&
                                   cell_radius_edge_ratio,&
                                   cell_size&
                                   )

    ! I don't trust the c_int kind of int, it might mess stuff up.
    npts=cnvert
    ntet=cntet

    ! Repoint the pointer or whatever it's called, fix so that I can access the
    ! c++ pointer from fortran.
    call c_f_pointer(cvertices,fvert,[3,npts])
    call c_f_pointer(ctetrahedrons,ftets,[4,ntet])

    ! I don't like pointers. You should not use them in Fortran. So I put it in an
    ! allocatable array instead.
    lo_allocate(vertices(3,npts))
    lo_allocate(tetrahedrons(4,ntet))
    vertices=fvert
    tetrahedrons=ftets

    ! free the c++ memory (I checked, it seems to work!!! I am so smart! SMRT!)
    call cgal_delete_double_pointer(cvertices)
    call cgal_delete_int_pointer(ctetrahedrons)
    ! and now it should be done!!
end subroutine

#endif

end module
