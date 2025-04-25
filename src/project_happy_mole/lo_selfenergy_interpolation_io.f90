module lo_selfenergy_interpolation
!!
!! Realspace representation of phonon self-energies
!!
use konstanter, only: r8, lo_iou, lo_hugeint, lo_huge, lo_sqtol, lo_exitcode_symmetry, lo_twopi, &
                      lo_bohr_to_A, lo_freqtol, lo_exitcode_param, lo_pi, lo_imag, lo_groupvel_ms_to_Hartreebohr
use gottochblandat, only: tochar, walltime, lo_progressbar_init, lo_progressbar, lo_chop, &
                          lo_linear_least_squares, lo_rsquare, lo_linspace, lo_linear_interpolation, lo_trapezoid_integration,&
                          lo_clean_fractional_coordinates, lo_cross
use geometryfunctions, only: lo_inscribed_sphere_in_box, lo_plane
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh, lo_wedge_mesh, lo_qpoint, lo_read_qmesh_from_file
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use lo_symmetry_of_interactions, only: lo_interaction_tensors
use type_symmetryoperation, only: lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_gemm, lo_dgels

use lineshape_helper, only: lo_spectralfunction_helper
implicit none

private
public :: lo_interpolated_selfenergy_grid

!> tetrahedron defined as planes
type lo_fancy_deltri_tet
    type(lo_plane), dimension(4) :: face
end type

!> a box containing tetrahedrons
type lo_fancy_deltri_onebox
    !> how many
    integer :: n = -lo_hugeint
    !> which tetrahedrons (points to full indices)
    integer, dimension(:), allocatable :: ind
    !> the actual tetrahedrons in an odd representation
    type(lo_fancy_deltri_tet), dimension(:), allocatable :: tet
end type

type lo_fancy_deltri_box
    !> physical dimensions of the box
    real(r8), dimension(3) :: boxdim = -lo_huge, boxmin = -lo_huge, boxmax = -lo_huge
    !> integer dimensions of the box
    integer, dimension(3) :: nbox = -lo_hugeint
    !> the boxes
    type(lo_fancy_deltri_onebox), dimension(:, :, :), allocatable :: b

    !> also, I keep the Fn interpolation here
    !> reference k0  first the reference e0 for each tetrahedron
!    real(r8), dimension(:,:), allocatable :: k0
    !> baseline energy per tetrahedron and band
!    real(r8), dimension(:,:), allocatable :: e0
    !> the gradient per tetrahedron and direction
!    real(r8), dimension(:,:,:), allocatable :: grad
    !> and the number of bands
!    integer :: nb=-lo_hugeint
    contains
        !> find the relevant box from the coordinates
        procedure :: boxind_from_coordinate
    !> find the correct tetrahedron from coordinates
    !procedure :: tetind
    !> interpolate electron energy to arbitrary r
!        procedure :: interpolate_Fn
        !> create locator
        procedure :: generate=>generate_triangulation
end type

type lo_interpolated_selfenergy_grid
    !> point locator thingy
    type(lo_fancy_deltri_box) :: box
    !> q-point mesh used to generate the self-energy
    class(lo_qpoint_mesh), allocatable :: qp
    !> number of energies
    integer :: n_energy=-lo_hugeint
    !> energy axis
    real(r8), dimension(:), allocatable :: omega
    !> real part of self-energy
    real(r8), dimension(:,:,:,:), allocatable :: sigma_Re
    !> imaginary part of self-energy
    real(r8), dimension(:,:,:,:), allocatable :: sigma_Im
    contains
        procedure :: generate=>generate_interpolated_selfenergy
        procedure :: read_from_hdf5=>read_interpolated_selfenergy_from_hdf5
end type

contains

!> create the interpolated self-energy object and store it to file
subroutine generate_interpolated_selfenergy(ise,sf,p,qp,dr,filename,mw,mem,verbosity)
    !> self-energy arranged on a q-point mesh, for interpolation
    class(lo_interpolated_selfenergy_grid), intent(out) :: ise
    !> self-energies on grid
    type(lo_spectralfunction_helper), intent(in) :: sf
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> filename
    character(len=*), intent(in) :: filename
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    type(lo_hdf5_helper) :: h5
    integer :: writerank

    init: block
        writerank=mw%n-1

        if ( mw%r .eq. writerank ) then
            call h5%init(__FILE__, __LINE__)
            call h5%open_file('write', trim(filename))
        endif

        ! Dump the q-grid, we need to have it embedded.
        if ( mw%r .eq. writerank ) then
            call h5%open_group('write','qmesh')
            call qp%write_to_file(p,'null',mem,verbosity,input_id=h5%group_id)
            call h5%close_group()
        endif
    end block init

    ! Dump the self-energy
    if ( mw%r .eq. writerank ) then
    dumpse: block
        real(r8), dimension(:,:,:,:), allocatable :: buf
        real(r8), dimension(:), allocatable :: x
        complex(r8), dimension(:,:), allocatable :: cm0,cm1,cm2
        integer :: iq,ie,imode

        allocate(buf(dr%n_mode,dr%n_mode,sf%n_energy,qp%n_irr_point))
        allocate(cm0(dr%n_mode,dr%n_mode))
        allocate(cm1(dr%n_mode,dr%n_mode))
        allocate(cm2(dr%n_mode,dr%n_mode))

        buf=0.0_r8
        do iq=1,qp%n_irr_point
        do ie=1,sf%n_energy
            cm0=0.0_r8
            do imode=1,dr%n_mode
                cm0(imode,imode)=sf%sigma_Re(ie,imode,iq)
                cm1(:,imode)=dr%iq(iq)%egv(:,imode) !*sqrt( dr%iq(iq)%omega(imode) )
            enddo
            cm2=matmul(cm1,cm0)
            cm1=matmul(cm2,conjg(transpose(cm1)))
            buf(:,:,ie,iq)=real(cm1,r8)
        enddo
        enddo
        call h5%store_data(buf,h5%file_id,'sigma_Re')

        if ( verbosity .gt. 0 ) write(*,*) '... stored real self-energy'

        buf=0.0_r8
        do iq=1,qp%n_irr_point
        do ie=1,sf%n_energy
            cm0=0.0_r8
            do imode=1,dr%n_mode
                cm0(imode,imode)=sf%sigma_Im(ie,imode,iq)
                cm1(:,imode)=dr%iq(iq)%egv(:,imode) !*sqrt( dr%iq(iq)%omega(imode) )
            enddo
            cm2=matmul(cm1,cm0)
            cm1=matmul(cm2,conjg(transpose(cm1)))
            buf(:,:,ie,iq)=real(cm1,r8)
        enddo
        enddo
        call h5%store_data(buf,h5%file_id,'sigma_Im')

        deallocate(buf)

        if ( verbosity .gt. 0 ) write(*,*) '... stored imaginary self-energy'

        allocate(x(sf%n_energy))
        call lo_linspace(0.0_r8,sf%max_energy,x)
        call h5%store_data(x,h5%file_id,'omega')
    end block dumpse
    endif

    if ( mw%r .eq. writerank ) then
        call h5%close_file()
        call h5%destroy(__FILE__, __LINE__)
    endif
end subroutine

subroutine read_interpolated_selfenergy_from_hdf5(ise,p,filename,mw,mem,verbosity)
    !> self-energy
    class(lo_interpolated_selfenergy_grid), intent(out) :: ise
    !> crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> filename
    character(len=*), intent(in) :: filename
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory helper
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! First we grab the raw data from file.
    readfile: block
        type(lo_hdf5_helper) :: h5

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'Reading interpolated self-energy from file'
        endif

        call h5%init(__FILE__,__LINE__)
        call h5%open_file('read',trim(filename))

        ! First we read the q-point mesh from file
        call h5%open_group('read','qmesh')
        call lo_read_qmesh_from_file(ise%qp,p,'null',mem,0,h5%group_id)
        call h5%close_group()

        if ( verbosity .gt. 0 ) then
            write(*,*) '... read q-mesh'
        endif

        call h5%read_data(ise%omega,h5%file_id,'omega')
        ise%n_energy = size(ise%omega)

        call h5%read_data(ise%sigma_Re,h5%file_id,'sigma_Re')
        call h5%read_data(ise%sigma_Im,h5%file_id,'sigma_Im')

        if ( verbosity .gt. 0 ) then
            write(*,*) '... read self-energy'
        endif

        call h5%close_file()
        call h5%destroy()

        if ( verbosity .gt. 0 ) then
            write(*,*) 'Done reading self-energy from file'
        endif
    end block readfile

    ! Generate triangulation thingy
    call ise%box%generate(ise%qp,p)
end subroutine

! Consistent index flattening? Impossibru to get consistent.
function flattenind(a1, a2, ix, iy, nb) result(i)
    integer, intent(in) :: a1, a2, ix, iy, nb
    integer :: i

    integer :: ia, ib

    ia = (a1 - 1)*3 + ix
    ib = (a2 - 1)*3 + iy
    i = (ib - 1)*nb + ia
end function

!> build a decently fast point locator for a wedge mesh.
subroutine generate_triangulation(box,qp,uc)
    !> the box with useful stuff
    class(lo_fancy_deltri_box), intent(out) :: box
    !> the wedge mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> the crystalstructure
    type(lo_crystalstructure), intent(in) :: uc
    !
    type(lo_plane) :: plane
    integer, dimension(3,4) :: faceind
    integer, dimension(3) :: bi
    integer :: i,j,k,l,ii,jj,kk
    real(r8), dimension(:,:), allocatable :: tetctr
    real(r8), dimension(:,:,:), allocatable :: tet
    real(r8), dimension(3) :: v0,v1,v2,rmin,rmax
    real(r8) :: maxtetrad,t0

    ! ! Fetch the tetrahedron corners
    ! t0=walltime()
    ! select type(qp)
    ! type is(lo_wedge_mesh)
    !     ! For a wedge mesh
    !     allocate(tet(3,4,qp%n_full_tet))
    !     tet=-ps_huge
    !     do i=1,qp%ntet_tot
    !         ! the centers of the tetrahdrons
    !         do j=1,4
    !             tet(:,j,i)=qp%ap( qp%atet(i)%gridind(j) )%w
    !         enddo
    !     enddo
    ! type is(lo_fft_mesh)
    !     ! For an fft mesh
    !     nt=qp%ntet
    !     allocate(tet(3,4,nt))
    !     tet=-ps_huge
    !     do i=1,qp%ntet
    !         ! the centers of the tetrahdrons
    !         v0=0.0_r8
    !         do j=1,4
    !             tet(:,j,i)=matmul(uc%inv_reciprocal_latticevectors, qp%ap( qp%tet(i)%gridind(j) )%v )
    !             do k=1,3
    !                 v0(k)=max(v0(k),tet(k,j,i))
    !             enddo
    !         enddo
    !         ! pbc-adjust
    !         do j=1,4
    !             do k=1,3
    !                 if ( tet(k,j,i)-v0(k) .lt. -0.5_r8 ) tet(k,j,i)=tet(k,j,i)+1.0_r8
    !                 if ( tet(k,j,i)-v0(k) .gt. 0.5_r8 ) tet(k,j,i)=tet(k,j,i)-1.0_r8
    !             enddo
    !             tet(:,j,i)=matmul(uc%reciprocal_latticevectors, tet(:,j,i))
    !         enddo
    !     enddo
    ! class default
    !     write(*,*) 'not done yet'
    !     stop
    ! end select

    allocate(tet(3,4,qp%n_full_tet))
    tet=-lo_huge
    do i=1,qp%n_full_tet
        ! the centers of the tetrahdrons
        v0=0.0_r8
        do j=1,4
            tet(:,j,i)=matmul(uc%inv_reciprocal_latticevectors, qp%ap( qp%at(i)%full_index(j) )%r )
            do k=1,3
                v0(k)=max(v0(k),tet(k,j,i))
            enddo
        enddo
        ! pbc-adjust
        do j=1,4
            do k=1,3
                if ( tet(k,j,i)-v0(k) .lt. -0.5_r8 ) tet(k,j,i)=tet(k,j,i)+1.0_r8
                if ( tet(k,j,i)-v0(k) .gt. 0.5_r8 ) tet(k,j,i)=tet(k,j,i)-1.0_r8
            enddo
            tet(:,j,i)=matmul(uc%reciprocal_latticevectors, tet(:,j,i))
        enddo
    enddo

    ! Get centers of tetrahedrons
    allocate(tetctr(3,qp%n_full_tet))
    tetctr=0.0_r8
    maxtetrad=0.0_r8
    do i=1,qp%n_full_tet
        do j=1,4
            tetctr(:,i)=tetctr(:,i)+tet(:,j,i)*0.25_r8
        enddo
        ! max diameter of bounding sphere for a tetrahdron
        do j=1,4
        do k=j+1,4
            maxtetrad=max(maxtetrad,norm2(tet(:,j,i)-tet(:,k,i)))
        enddo
        enddo
    enddo

    ! get the bounding region
    rmin=lo_huge
    rmax=-lo_huge
    do i=1,qp%n_full_tet
        do k=1,4
            v0=tet(:,k,i)
            do j=1,3
                rmin(j)=min(rmin(j),v0(j)-maxtetrad*0.01_r8)
                rmax(j)=max(rmax(j),v0(j)+maxtetrad*0.01_r8)
            enddo
        enddo
    enddo
    ! Make the box
    do j=1,3
        box%boxmin(j)=rmin(j)
        box%boxmax(j)=rmax(j)
        box%boxdim(j)=rmax(j)-rmin(j)
        box%nbox(j)=floor(box%boxdim(j)/maxtetrad)
    enddo
    allocate(box%b(box%nbox(1),box%nbox(2),box%nbox(3)))

    ! stuff tetrahdrons into boxes, first reset the counter
    do i=1,box%nbox(1)
    do j=1,box%nbox(2)
    do k=1,box%nbox(3)
        box%b(i,j,k)%n=0
    enddo
    enddo
    enddo
    ! count tetrahdrons per box
    do i=1,qp%n_full_tet
        bi=box%boxind_from_coordinate(tetctr(:,i))
        box%b(bi(1),bi(2),bi(3))%n=box%b(bi(1),bi(2),bi(3))%n+1
    enddo
    ! make space
    do i=1,box%nbox(1)
    do j=1,box%nbox(2)
    do k=1,box%nbox(3)
        l=box%b(i,j,k)%n
        if ( l .gt. 0 ) then
            allocate(box%b(i,j,k)%ind( l ))
            allocate(box%b(i,j,k)%tet( l ))
        endif
        box%b(i,j,k)%n=0
    enddo
    enddo
    enddo
    ! indices to faces
    faceind(:,1)=[1,2,3]
    faceind(:,2)=[1,2,4]
    faceind(:,3)=[1,3,4]
    faceind(:,4)=[2,3,4]
    ! store tetrahdrons per box
    do i=1,qp%n_full_tet
        bi=box%boxind_from_coordinate(tetctr(:,i))
        box%b(bi(1),bi(2),bi(3))%n=box%b(bi(1),bi(2),bi(3))%n+1
        box%b(bi(1),bi(2),bi(3))%ind( box%b(bi(1),bi(2),bi(3))%n )=i
        v0=tetctr(:,i)
        ! one plane per face, where a the normal points away from the center
        do j=1,4
            ii=faceind(1,j)
            jj=faceind(2,j)
            kk=faceind(3,j)
            v1=tet(:,ii,i)
            v2=lo_cross(tet(:,jj,i)-tet(:,ii,i),tet(:,kk,i)-tet(:,ii,i))
            call plane%generate(normal=v2,point=v1)
            if ( plane%distance_to_point(v0) .gt. 0.0_r8 ) then
                call plane%generate(normal=-v2,point=v1)
            endif
            box%b(bi(1),bi(2),bi(3))%tet( box%b(bi(1),bi(2),bi(3))%n )%face(j)=plane
        enddo
    enddo
end subroutine

!> get the index of the relevant box from coordinates
function boxind_from_coordinate(box,r) result(ind)
    !> the box
    class(lo_fancy_deltri_box), intent(in) :: box
    !> the point
    real(r8), dimension(3) :: r
    !> the indices
    integer, dimension(3) :: ind
    !
    integer :: j

    ! no sanity checks at all.
    do j=1,3
        ind(j)=floor( ( (r(j)-box%boxmin(j))/box%boxdim(j) )*box%nbox(j) )+1
    enddo
end function

!> locate the correct tetrahedron for an arbitrary point
function tetind(box,r) result(ind)
    !> the box
    class(lo_fancy_deltri_box), intent(in) :: box
    !> the point
    real(r8), dimension(3) :: r
    !> the index
    integer :: ind
    !
    type(lo_fancy_deltri_tet) :: tet
    integer :: i,j,k,l,ii,jj
    integer, dimension(3) :: bi
    real(r8) :: f0

    ! First, get what box the point is in:
    bi=box%boxind_from_coordinate(r)
    l=0
    ! first check just in this box, since the odds are quite large it is there.
    tetloop1: do ii=1,box%b(bi(1),bi(2),bi(3))%n
        tet=box%b(bi(1),bi(2),bi(3))%tet(ii)
        l=0
        do i=1,4
            f0=tet%face(i)%distance_to_point(r)
            if ( f0 .gt. 1E-6_r8 ) then
                cycle tetloop1
            else
                l=l+1
            endif
        enddo
        if ( l .eq. 4 ) then
            ! found it
            ind=box%b(bi(1),bi(2),bi(3))%ind(ii)
            return
        endif
    enddo tetloop1
    ! if we made it here, we have to check all other boxes
    do i=max(1,bi(1)-1),min(box%nbox(1),bi(1)+1)
    do j=max(1,bi(2)-1),min(box%nbox(2),bi(2)+1)
    do k=max(1,bi(3)-1),min(box%nbox(3),bi(3)+1)
        ! not the middle box again
        if ( sum(abs([i,j,k]-bi)) .eq. 0 ) cycle
        ! check all the tetrahedrons here
        tetloop2: do ii=1,box%b(i,j,k)%n
            tet=box%b(i,j,k)%tet(ii)
            l=0
            do jj=1,4
                f0=tet%face(jj)%distance_to_point(r)
                if ( f0 .gt. 1E-6_r8 ) then
                    cycle tetloop2
                else
                    l=l+1
                endif
            enddo
            if ( l .eq. 4 ) then
                ! found it
                ind=box%b(i,j,k)%ind(ii)
                return
            endif
        enddo tetloop2
    enddo
    enddo
    enddo

    ! if I made it here, something is seriously wrong
    write(*,*) 'found nothing'
    write(*,*) bi
    write(*,*) box%nbox
    ind=0
end function

end module
