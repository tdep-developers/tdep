#include "precompilerdefinitions"
module type_fast_interpolation
!!
!! Tetrahedral interpolation
!!
use konstanter, only: flyt, lo_huge, lo_hugeint, lo_tol
use gottochblandat, only: lo_cross, walltime, lo_sqnorm, lo_signed_tetrahedron_volume
use geometryfunctions, only: lo_plane
use type_crystalstructure, only: lo_crystalstructure
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_qpointmesh, only: lo_qpoint_mesh, lo_wedge_mesh, lo_fft_mesh
use type_symmetryoperation, only: lo_operate_on_vector
implicit none

private
public :: lo_fancy_deltri_box

!> tetrahedron defined as planes
type lo_fancy_deltri_tet
    type(lo_plane), dimension(4) :: face
end type

! a box containing tetrahedrons
type lo_fancy_deltri_onebox
    !> how many
    integer :: n = -lo_hugeint
    !> which tetrahedrons
    integer, dimension(:), allocatable :: ind
    !> the actual tetrahedrons in an odd representation
    type(lo_fancy_deltri_tet), dimension(:), allocatable :: tet
end type

type lo_fancy_deltri_box
    !> physical dimensions of the box
    real(flyt), dimension(3) :: boxdim = -lo_huge, boxmin = -lo_huge, boxmax = -lo_huge
    !> integer dimensions of the box
    integer, dimension(3) :: nbox = -lo_hugeint
    !> the boxes
    type(lo_fancy_deltri_onebox), dimension(:, :, :), allocatable :: b

    !> also, I keep the Fn interpolation here
    !> reference k0  first the reference e0 for each tetrahedron
!    real(flyt), dimension(:,:), allocatable :: k0
    !> baseline energy per tetrahedron and band
!    real(flyt), dimension(:,:), allocatable :: e0
    !> the gradient per tetrahedron and direction
!    real(flyt), dimension(:,:,:), allocatable :: grad
    !> and the number of bands
!    integer :: nb=-lo_hugeint
contains
    !> find the relevant box from the coordinates
    !procedure :: boxind_from_coordinate
    !> find the correct tetrahedron from coordinates
    !procedure :: tetind
    !> interpolate electron energy to arbitrary r
!        procedure :: interpolate_Fn
    !> generation
    !procedure :: generate
end type

contains

!!> Interpolate Fn to an abritrary q-point
!subroutine interpolate_Fn(box,uc,qp,dr,q,mode,Fn)
!    !> interpolator thing
!    class(lo_fancy_deltri_box), intent(in) :: box
!    !> crystal structure
!    type(lo_crystalstructure), intent(in) :: uc
!    !> q-mesh
!    class(lo_qpoint_mesh), intent(in) :: qp
!    !> dispersions
!    type(lo_phonon_dispersions), intent(in) :: dr
!    !> q-point
!    real(flyt), dimension(3), intent(in) :: q
!    !> mode
!    integer, intent(in) :: mode
!    !> resulting Fn
!    real(flyt), dimension(3), intent(out) :: Fn
!
!    real(flyt), dimension(3,4) :: tetv
!    real(flyt), dimension(3,4) :: tet
!    real(flyt), dimension(3) :: r,v0,v1,v2,v3,r1,r2,r3,b
!    real(flyt) :: f0,e0
!
!    integer, dimension(4) :: gi
!    integer :: t,i,j,k
!
!    Fn=0.0_flyt
!    ! First translate q to the proper domain
!    select type(qp)
!    class is(lo_wedge_mesh)
!        ! in this case, move it to the first BZ
!        if ( lo_sqnorm(q) .gt. uc%bz%rmin**2 ) then
!            r=q-uc%bz%gshift(q)
!        else
!            r=q
!        endif
!    class is(lo_fft_mesh)
!        ! or, move it to 0-1 coordinates
!        r=uc%cartesian_to_fractional(q,reciprocal=.true.)
!        r=uc%fractional_to_cartesian(r,reciprocal=.true.)
!    end select
!
!    ! Figure out the relevant tetrahedron
!    t=-1
!    t=box%tetind(r)
!
!    ! Fetch tetrahedron corners, and their respective grid indices
!    select type(qp)
!    type is(lo_wedge_mesh)
!        do j=1,4
!            tet(:,j)=qp%ap( qp%atet(t)%gridind(j) )%w
!            gi(j)=qp%atet(t)%gridind(j)
!        enddo
!    type is(lo_fft_mesh)
!        v0=0.0_flyt
!        do j=1,4
!            gi(j)=qp%tet(t)%gridind(j)
!            tet(:,j)=matmul(uc%inv_reciprocal_latticevectors, qp%ap( qp%tet(t)%gridind(j) )%v )
!            do k=1,3
!                v0(k)=max(v0(k),tet(k,j))
!            enddo
!        enddo
!        ! pbc-adjust
!        do j=1,4
!            do k=1,3
!                if ( tet(k,j)-v0(k) .lt. -0.5_flyt ) tet(k,j)=tet(k,j)+1.0_flyt
!                if ( tet(k,j)-v0(k) .gt. 0.5_flyt ) tet(k,j)=tet(k,j)-1.0_flyt
!            enddo
!            tet(:,j)=matmul(uc%reciprocal_latticevectors, tet(:,j))
!        enddo
!    end select
!
!    ! Get Fn at the corners of the tetrahedron
!    do i=1,4
!        j=qp%ap( gi(i) )%operation
!        k=qp%ap( gi(i) )%irrind
!        if ( j .gt. 0 ) then
!            v0=lo_operate_on_vector(uc%sym%op(j),dr%iq(k)%Fn(:,mode),reciprocal=.true.)
!        else
!            v0=-lo_operate_on_vector(uc%sym%op(abs(j)),dr%iq(k)%Fn(:,mode),reciprocal=.true.)
!        endif
!        tetv(:,i)=v0
!    enddo
!
!    ! Start building the contragradient
!    f0=-lo_signed_tetrahedron_volume(tet)*6
!    v0=tet(:,1)
!    v1=tet(:,2)-v0
!    v2=tet(:,3)-v0
!    v3=tet(:,4)-v0
!    r1=lo_cross(v2,v3)/f0
!    r2=lo_cross(v3,v1)/f0
!    r3=lo_cross(v1,v2)/f0
!
!    do i=1,3 ! cartesian dimension
!        e0=tetv(i,1)
!        b=0.0_flyt
!        b=b+(tetv(i,2)-e0)*r1
!        b=b+(tetv(i,3)-e0)*r2
!        b=b+(tetv(i,4)-e0)*r3
!        Fn(i)=e0+dot_product(r-v0,b)
!    enddo
!end subroutine

! !> build a decently fast point locator for a wedge mesh.
! subroutine generate(box,qp,uc)
!     !> the box with useful stuff
!     class(lo_fancy_deltri_box), intent(out) :: box
!     !> the wedge mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> the crystalstructure
!     type(lo_crystalstructure), intent(in) :: uc
!     !
!     type(lo_plane) :: plane
!     integer, dimension(3,4) :: faceind
!     integer, dimension(3) :: bi
!     integer :: i,j,k,l,ii,jj,kk,nt
!     real(flyt), dimension(:,:), allocatable :: tetctr
!     real(flyt), dimension(:,:,:), allocatable :: tet
!     real(flyt), dimension(3) :: v0,v1,v2,rmin,rmax
!     real(flyt) :: maxtetrad,t0
!
!     ! Fetch the tetrahedron corners
!     t0=walltime()
!     select type(qp)
!     type is(lo_wedge_mesh)
!         ! For a wedge mesh
!         nt=qp%ntet_tot
!         lo_allocate(tet(3,4,nt))
!         do i=1,qp%ntet_tot
!             ! the centers of the tetrahdrons
!             do j=1,4
!                 tet(:,j,i)=qp%ap( qp%atet(i)%gridind(j) )%w
!             enddo
!         enddo
!     type is(lo_fft_mesh)
!         ! For an fft mesh
!         nt=qp%ntet
!         lo_allocate(tet(3,4,nt))
!         do i=1,qp%ntet
!             ! the centers of the tetrahdrons
!             v0=0.0_flyt
!             do j=1,4
!                 tet(:,j,i)=matmul(uc%inv_reciprocal_latticevectors, qp%ap( qp%tet(i)%gridind(j) )%v )
!                 do k=1,3
!                     v0(k)=max(v0(k),tet(k,j,i))
!                 enddo
!             enddo
!             ! pbc-adjust
!             do j=1,4
!                 do k=1,3
!                     if ( tet(k,j,i)-v0(k) .lt. -0.5_flyt ) tet(k,j,i)=tet(k,j,i)+1.0_flyt
!                     if ( tet(k,j,i)-v0(k) .gt. 0.5_flyt ) tet(k,j,i)=tet(k,j,i)-1.0_flyt
!                 enddo
!                 tet(:,j,i)=matmul(uc%reciprocal_latticevectors, tet(:,j,i))
!             enddo
!         enddo
!     class default
!         write(*,*) 'not done yet'
!         stop
!     end select
!
!     ! Get centers of tetrahedrons
!     lo_allocate(tetctr(3,nt))
!     tetctr=0.0_flyt
!     maxtetrad=0.0_flyt
!     do i=1,nt
!         do j=1,4
!             tetctr(:,i)=tetctr(:,i)+tet(:,j,i)*0.25_flyt
!         enddo
!         ! max diameter of bounding sphere for a tetrahdron
!         do j=1,4
!         do k=j+1,4
!             maxtetrad=max(maxtetrad,norm2(tet(:,j,i)-tet(:,k,i)))
!         enddo
!         enddo
!     enddo
!
!     ! get the bounding region
!     rmin=lo_huge
!     rmax=-lo_huge
!     do i=1,nt
!         do k=1,4
!             v0=tet(:,k,i)
!             do j=1,3
!                 rmin(j)=min(rmin(j),v0(j)-maxtetrad*0.01_flyt)
!                 rmax(j)=max(rmax(j),v0(j)+maxtetrad*0.01_flyt)
!             enddo
!         enddo
!     enddo
!     ! Make the box
!     do j=1,3
!         box%boxmin(j)=rmin(j)
!         box%boxmax(j)=rmax(j)
!         box%boxdim(j)=rmax(j)-rmin(j)
!         box%nbox(j)=floor(box%boxdim(j)/maxtetrad)
!     enddo
!     lo_allocate(box%b(box%nbox(1),box%nbox(2),box%nbox(3)))
!     ! stuff tetrahdrons into boxes, first reset the counter
!     do i=1,box%nbox(1)
!     do j=1,box%nbox(2)
!     do k=1,box%nbox(3)
!         box%b(i,j,k)%n=0
!     enddo
!     enddo
!     enddo
!     ! count tetrahdrons per box
!     do i=1,nt
!         bi=box%boxind_from_coordinate(tetctr(:,i))
!         box%b(bi(1),bi(2),bi(3))%n=box%b(bi(1),bi(2),bi(3))%n+1
!     enddo
!     ! make space
!     do i=1,box%nbox(1)
!     do j=1,box%nbox(2)
!     do k=1,box%nbox(3)
!         l=box%b(i,j,k)%n
!         if ( l .gt. 0 ) then
!             lo_allocate(box%b(i,j,k)%ind( l ))
!             lo_allocate(box%b(i,j,k)%tet( l ))
!         endif
!         box%b(i,j,k)%n=0
!     enddo
!     enddo
!     enddo
!     ! indices to faces
!     faceind(:,1)=[1,2,3]
!     faceind(:,2)=[1,2,4]
!     faceind(:,3)=[1,3,4]
!     faceind(:,4)=[2,3,4]
!     ! store tetrahdrons per box
!     do i=1,nt
!         bi=box%boxind_from_coordinate(tetctr(:,i))
!         box%b(bi(1),bi(2),bi(3))%n=box%b(bi(1),bi(2),bi(3))%n+1
!         box%b(bi(1),bi(2),bi(3))%ind( box%b(bi(1),bi(2),bi(3))%n )=i
!         v0=tetctr(:,i)
!         ! one plane per face, where a the normal points away from the center
!         do j=1,4
!             ii=faceind(1,j)
!             jj=faceind(2,j)
!             kk=faceind(3,j)
!             v1=tet(:,ii,i)
!             v2=lo_cross(tet(:,jj,i)-tet(:,ii,i),tet(:,kk,i)-tet(:,ii,i))
!             call plane%generate(normal=v2,point=v1)
!             if ( plane%distance_to_point(v0) .gt. 0.0_flyt ) then
!                 call plane%generate(normal=-v2,point=v1)
!             endif
!             box%b(bi(1),bi(2),bi(3))%tet( box%b(bi(1),bi(2),bi(3))%n )%face(j)=plane
!         enddo
!     enddo
! end subroutine
!
! !> get the index of the relevant box from coordinates
! function boxind_from_coordinate(box,r) result(ind)
!     !> the box
!     class(lo_fancy_deltri_box), intent(in) :: box
!     !> the point
!     real(flyt), dimension(3) :: r
!     !> the indices
!     integer, dimension(3) :: ind
!     !
!     integer :: j
!     ! no sanity checks at all.
!     do j=1,3
!         ind(j)=floor( ( (r(j)-box%boxmin(j))/box%boxdim(j) )*box%nbox(j) )+1
!     enddo
! end function
!
! !> locate the correct tetrahedron for an arbitrary point
! function tetind(box,r) result(ind)
!     !> the box
!     class(lo_fancy_deltri_box), intent(in) :: box
!     !> the point
!     real(flyt), dimension(3) :: r
!     !> the index
!     integer :: ind
!     !
!     type(lo_fancy_deltri_tet) :: tet
!     integer :: i,j,k,l,ii,jj
!     integer, dimension(3) :: bi
!     real(flyt) :: f0
!
!     ! First, get what box the point is in:
!     bi=box%boxind_from_coordinate(r)
!     l=0
!     ! first check just in this box, since the odds are quite large it is there.
!     tetloop1: do ii=1,box%b(bi(1),bi(2),bi(3))%n
!         tet=box%b(bi(1),bi(2),bi(3))%tet(ii)
!         l=0
!         do i=1,4
!             f0=tet%face(i)%distance_to_point(r)
!             if ( f0 .gt. lo_tol ) then
!                 cycle tetloop1
!             else
!                 l=l+1
!             endif
!         enddo
!         if ( l .eq. 4 ) then
!             ! found it
!             ind=box%b(bi(1),bi(2),bi(3))%ind(ii)
!             return
!         endif
!     enddo tetloop1
!     ! if we made it here, we have to check all other boxes
!     do i=max(1,bi(1)-1),min(box%nbox(1),bi(1)+1)
!     do j=max(1,bi(2)-1),min(box%nbox(2),bi(2)+1)
!     do k=max(1,bi(3)-1),min(box%nbox(3),bi(3)+1)
!         ! not the middle box again
!         if ( sum(abs([i,j,k]-bi)) .eq. 0 ) cycle
!         ! check all the tetrahedrons here
!         tetloop2: do ii=1,box%b(i,j,k)%n
!             tet=box%b(i,j,k)%tet(ii)
!             l=0
!             do jj=1,4
!                 f0=tet%face(jj)%distance_to_point(r)
!                 if ( f0 .gt. lo_tol ) then
!                     cycle tetloop2
!                 else
!                     l=l+1
!                 endif
!             enddo
!             if ( l .eq. 4 ) then
!                 ! found it
!                 ind=box%b(i,j,k)%ind(ii)
!                 return
!             endif
!         enddo tetloop2
!     enddo
!     enddo
!     enddo
!
!     ! if I made it here, something is seriously wrong
!     write(*,*) 'found nothing'
!     write(*,*) bi
!     write(*,*) box%nbox
!     ind=0
! end function

end module
