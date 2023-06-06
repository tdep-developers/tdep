module lo_verletboxes
!!
!! Verlet lists for fast distance lookups.
!! @TODO update for periodic boundary conditions when needed
!! @TODO make dimensions automatic
!! @TODO automatic determination of point location tolerance
!! @TODO automatically determine largest possible number of neighbours, check box + neighbouring boxes
!! @TODO make it MPI distributed?
!!
use konstanter, only: r8,lo_huge,lo_hugeint,lo_sqtol,lo_exitcode_symmetry,lo_exitcode_param
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_stop_gracefully

implicit none
private
public :: lo_verletbox

!> list of points in a verlet box
type lo_verletbox_box
    !> how many points in this box
    integer :: n=-lo_hugeint
    !> indices to points in this box
    integer, dimension(:), allocatable :: ind
end type

!> minimal Verlet-box to generate distancelists and things like that
type lo_verletbox
    !> box divisions
    integer :: nx=-lo_hugeint
    integer :: ny=-lo_hugeint
    integer :: nz=-lo_hugeint
    !> lower bounds
    real(r8), dimension(3), private :: rmin=lo_huge
    !> upper bounds
    real(r8), dimension(3), private :: rmax=lo_huge
    !> scalefactor per dimension
    real(r8), dimension(3), private :: ird=lo_huge
    !> boxes with points in them
    type(lo_verletbox_box), dimension(:,:,:), allocatable :: box
    contains
        !> stuff the particles into boxes
        procedure :: generate=>add_particles_in_boxes
        !> box-indices from a point
        procedure :: boxind=>boxind_from_coordinate
        !> locate index of a point
        procedure :: locate=>locate_index_of_point
        !> guess sensible box dimensions
        procedure, nopass :: boxdim=>sensible_box_dimensions
        !> size in memory
        procedure :: size_in_mem
        !> destroy
        procedure :: destroy
end type

contains

!> decide on sensible box dimensions
function sensible_box_dimensions(r,n,tol) result(boxdim)
    !> set of points to chop into boxes
    real(r8), dimension(:,:), intent(in) :: r
    !> total number of boxes I want
    integer, intent(in) :: n
    !> tolerance for checking zero values
    real(r8), intent(in) :: tol
    !> the dimensions in each direction
    integer, dimension(3) :: boxdim

    real(r8), dimension(3) :: rmin,rmax,v0
    real(r8) :: f0
    integer :: i

    ! Get the ranges
    rmin=lo_huge
    rmax=-lo_huge
    do i=1,size(r,2)
        rmin(1)=min(rmin(1),r(1,i))
        rmin(2)=min(rmin(2),r(2,i))
        rmin(3)=min(rmin(3),r(3,i))
        rmax(1)=max(rmax(1),r(1,i))
        rmax(2)=max(rmax(2),r(2,i))
        rmax(3)=max(rmax(3),r(3,i))
    enddo

    ! This version is just if we want some constant number of boxes.
    v0(1)=rmax(1)-rmin(1)
    v0(2)=rmax(2)-rmin(2)
    v0(3)=rmax(3)-rmin(3)
    if ( minval(abs(v0)) .lt. tol ) then
        call lo_stop_gracefully(['Building Verlet boxes that are flat in one or more dimensions.'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    f0=real(n,r8)
    v0=v0*(f0/product(v0))**(1.0_r8/3.0_r8)
    boxdim=ceiling(v0)
    !@TODO add sensible box size for certain cutoff.
end function

!> get the indices of a box from the coordinate of a point
pure subroutine boxind_from_coordinate(vb,r,bi,bj,bk)
    !> verlet boxes
    class(lo_verletbox), intent(in) :: vb
    !> coordinates of point
    real(r8), dimension(3), intent(in) :: r
    !> box-indices
    integer, intent(out) :: bi,bj,bk

    bi=floor( (r(1)-vb%rmin(1))*vb%ird(1) )+1
    bj=floor( (r(2)-vb%rmin(2))*vb%ird(2) )+1
    bk=floor( (r(3)-vb%rmin(3))*vb%ird(3) )+1
end subroutine

!> Locate the index of a point using Verlet boxes. Index -1 means the point does not exist.
function locate_index_of_point(vb,points,r,singlebox) result(ind)
    !> verlet boxes
    class(lo_verletbox), intent(in) :: vb
    !> all possible points
    real(r8), dimension(:,:), intent(in) :: points
    !> coordinates of point
    real(r8), dimension(3), intent(in) :: r
    !> check only the first box that comes to mind
    logical, intent(in), optional :: singlebox
    !> index of points
    integer :: ind

    real(r8), dimension(3) :: w
    real(r8) :: f0
    integer :: bi,bj,bk,i,j,k,l,ii

    ! Assume I don't find anything
    ind=-1
    ! First some sanity tests
    j=0
    do i=1,3
        if ( r(i) .lt. vb%rmin(i) ) j=j+1
        if ( r(i) .gt. vb%rmax(i) ) j=j+1
    enddo
    if ( j .ne. 0 ) then
        ! This means the point is not in any of the boxes.
        return
    endif

    ! Which box to start looking?
    bi=floor( (r(1)-vb%rmin(1))*vb%ird(1) )+1
    bj=floor( (r(2)-vb%rmin(2))*vb%ird(2) )+1
    bk=floor( (r(3)-vb%rmin(3))*vb%ird(3) )+1
    ! First look in this box, decent odds the point is there
    do l=1,vb%box(bi,bj,bk)%n
        ii=vb%box(bi,bj,bk)%ind(l)
        w=points(:,ii)-r
        f0=w(1)*w(1)+w(2)*w(2)+w(3)*w(3)
        if ( f0 .lt. lo_sqtol ) then
            ! Yup, found it
            ind=ii
            return
        endif
    enddo

    if ( present(singlebox) ) then
        ! Maybe I only want to check one box at a time.
        if ( singlebox ) return
    endif

    ! Then look for it in the adjacent boxes
    do i=max(1,bi-1),min(vb%nx,bi+1)
    do j=max(1,bj-1),min(vb%ny,bj+1)
    do k=max(1,bk-1),min(vb%nz,bk+1)
        do l=1,vb%box(i,j,k)%n
            ii=vb%box(i,j,k)%ind(l)
            w=points(:,ii)-r
            f0=w(1)*w(1)+w(2)*w(2)+w(3)*w(3)
            if ( f0 .lt. lo_sqtol ) then
                ind=ii
                return
            endif
        enddo
    enddo
    enddo
    enddo
end function

!> Put particles in Verlet boxes
subroutine add_particles_in_boxes(vb,r,ndim,mem)
    !> verlet boxes
    class(lo_verletbox), intent(out) :: vb
    !> particles
    real(r8), dimension(:,:), intent(in) :: r
    !> number of boxes
    integer, dimension(3), intent(in) :: ndim
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(3) :: v0
    integer :: n,i,j,k,l

    ! number of particles
    n=size(r,2)

    ! First get the bounds
    vb%rmin=lo_huge
    vb%rmax=-lo_huge
    do i=1,n
    do j=1,3
        vb%rmin(j)=min(vb%rmin(j),r(j,i))
        vb%rmax(j)=max(vb%rmax(j),r(j,i))
    enddo
    enddo
    ! Slight tolerance, should be fine.
    v0=(vb%rmax-vb%rmin)*lo_sqtol+lo_sqtol
    vb%rmin=vb%rmin-v0
    vb%rmax=vb%rmax+v0
    ! Scaling factor
    vb%ird=(ndim*1.0_r8)/(vb%rmax-vb%rmin)
    ! Dimensions
    vb%nx=ndim(1)
    vb%ny=ndim(2)
    vb%nz=ndim(3)

    ! Make space
    allocate(vb%box( vb%nx,vb%ny,vb%nz ))

    ! Reset the counter
    do i=1,vb%nx
    do j=1,vb%ny
    do k=1,vb%nz
        vb%box(i,j,k)%n=0
    enddo
    enddo
    enddo
    ! Count members per box
    do l=1,n
        v0=r(:,l)
        i=floor( (v0(1)-vb%rmin(1))*vb%ird(1) )+1
        j=floor( (v0(2)-vb%rmin(2))*vb%ird(2) )+1
        k=floor( (v0(3)-vb%rmin(3))*vb%ird(3) )+1
        vb%box(i,j,k)%n=vb%box(i,j,k)%n+1
    enddo
    ! Make some space and reset counter
    do i=1,vb%nx
    do j=1,vb%ny
    do k=1,vb%nz
        if ( vb%box(i,j,k)%n .gt. 0 ) then
            call mem%allocate(vb%box(i,j,k)%ind,vb%box(i,j,k)%n,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        endif
        vb%box(i,j,k)%n=0
    enddo
    enddo
    enddo
    ! Store members
    do l=1,n
        v0=r(:,l)
        i=floor( (v0(1)-vb%rmin(1))*vb%ird(1) )+1
        j=floor( (v0(2)-vb%rmin(2))*vb%ird(2) )+1
        k=floor( (v0(3)-vb%rmin(3))*vb%ird(3) )+1
        vb%box(i,j,k)%n=vb%box(i,j,k)%n+1
        vb%box(i,j,k)%ind( vb%box(i,j,k)%n )=l
    enddo
    ! Check that everything got stored
    l=0
    do i=1,vb%nx
    do j=1,vb%ny
    do k=1,vb%nz
        l=l+vb%box(i,j,k)%n
    enddo
    enddo
    enddo
    if ( l .ne. n ) then
        call lo_stop_gracefully(['Failed putting all particles in boxes.'],lo_exitcode_symmetry)
    endif
end subroutine

!> measure size in memory, in bytes
pure function size_in_mem(vb) result(bytes)
    !> verlet boxes
    class(lo_verletbox), intent(in) :: vb
    !> size in bytesory, in bytes
    integer :: bytes

    integer :: i,j,k

    bytes=0
    bytes=bytes+storage_size(vb)
    if ( allocated(vb%box) ) then
        do k=1,size(vb%box,3)
        do j=1,size(vb%box,2)
        do i=1,size(vb%box,1)
            bytes=bytes+storage_size(vb%box(i,j,k))
            if ( allocated(vb%box(i,j,k)%ind) ) bytes=bytes+storage_size(vb%box(i,j,k)%ind)*size(vb%box(i,j,k)%ind)
        enddo
        enddo
        enddo
    endif
    bytes=bytes/8
end function

!> destroy the object
subroutine destroy(vb,mem)
    !> verlet boxes
    class(lo_verletbox), intent(inout) :: vb
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    integer :: i,j,k

    if ( allocated(vb%box) ) then
        do k=1,size(vb%box,3)
        do j=1,size(vb%box,2)
        do i=1,size(vb%box,1)
            if ( allocated(vb%box(i,j,k)%ind) ) then
                call mem%deallocate(vb%box(i,j,k)%ind,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
            endif
        enddo
        enddo
        enddo
        deallocate(vb%box)
    endif
end subroutine

end module
