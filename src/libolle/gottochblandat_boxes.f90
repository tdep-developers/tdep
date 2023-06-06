#include "precompilerdefinitions"
submodule (gottochblandat) gottochblandat_boxes
implicit none
contains

!> box indices
module subroutine boxind_from_coordinate(vb,r,bi,bj,bk)
    !> verlet boxes
    class(lo_verletbox), intent(in) :: vb
    !> coordinates of point
    real(flyt), dimension(3), intent(in) :: r
    !> box-indices
    integer, intent(out) :: bi,bj,bk

    bi=floor( (r(1)-vb%rmin(1))*vb%ird(1) )+1
    bj=floor( (r(2)-vb%rmin(2))*vb%ird(2) )+1
    bk=floor( (r(3)-vb%rmin(3))*vb%ird(3) )+1
end subroutine

!> Locate the index of a point using Verlet boxes. Index -1 means the point does not exist.
module function locate_index_of_point(vb,points,r,singlebox) result(ind)
    !> verlet boxes
    class(lo_verletbox), intent(in) :: vb
    !> all possible points
    real(flyt), dimension(:,:), intent(in) :: points
    !> coordinates of point
    real(flyt), dimension(3), intent(in) :: r
    !> check only the first box that comes to mind
    logical, intent(in), optional :: singlebox
    !> index of points
    integer :: ind


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
        return
    endif

    ! Which box to start looking?
    bi=floor( (r(1)-vb%rmin(1))*vb%ird(1) )+1
    bj=floor( (r(2)-vb%rmin(2))*vb%ird(2) )+1
    bk=floor( (r(3)-vb%rmin(3))*vb%ird(3) )+1
    ! First look in this box, decent odds the point is there
    do l=1,vb%box(bi,bj,bk)%n
        ii=vb%box(bi,bj,bk)%ind(l)
        if ( lo_sqnorm(points(:,ii)-r) .lt. lo_sqtol ) then
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
            if ( lo_sqnorm(points(:,ii)-r) .lt. lo_sqtol ) then
                ind=ii
                return
            endif
        enddo
    enddo
    enddo
    enddo
end function

!> Put particles in Verlet boxes
module subroutine add_particles_in_boxes(vb,r,ndim)
    !> verlet boxes
    class(lo_verletbox), intent(out) :: vb
    !> particles
    real(flyt), dimension(:,:), intent(in) :: r
    !> number of boxes
    integer, dimension(3), intent(in) :: ndim

    real(flyt), dimension(3) :: v0
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
    v0=(vb%rmax-vb%rmin)*lo_sqtol
    vb%rmin=vb%rmin-v0
    vb%rmax=vb%rmax+v0
    ! Scaling factor
    vb%ird=(ndim*1.0_flyt)/(vb%rmax-vb%rmin)
    ! Dimensions
    vb%nx=ndim(1)
    vb%ny=ndim(2)
    vb%nz=ndim(3)
    ! Make space
    lo_allocate(vb%box( vb%nx,vb%ny,vb%nz ))

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
        if ( vb%box(i,j,k)%n .gt. 0 ) allocate(vb%box(i,j,k)%ind( vb%box(i,j,k)%n ))
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
        call lo_stop_gracefully(['Failed putting all particles in boxes.'],lo_exitcode_symmetry,__FILE__,__LINE__)
    endif
end subroutine

end submodule
