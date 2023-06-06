submodule (type_qpointmesh) type_qpointmesh_commensurate
use konstanter, only: lo_tiny
use geometryfunctions, only: lo_inscribed_sphere_in_box,lo_bounding_sphere_of_box
use gottochblandat, only: lo_reciprocal_basis,lo_return_unique_indices

implicit none
contains

!> Generate a q-point mesh from a supercell
module subroutine lo_generate_commensurate_qmesh(qp,uc,supercellmatrix,mw,mem,verbosity)
    !> resulting q-point mesh
    class(lo_commensurate_mesh), intent(out) :: qp
    !> unit cell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercellmatrix
    real(r8), dimension(3,3), intent(in) :: supercellmatrix
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    ! type(lo_verletbox) :: vb
    ! real(r8), dimension(:,:), allocatable :: qvec
    real(r8) :: timer,t0,t1
    integer :: i

    ! Dummy things so that I can compile with warnings.
    timer=walltime()
    t0=timer
    t1=timer
    call mem%tick()
    call mw%barrier()
    i=uc%na
    i=verbosity
    t0=supercellmatrix(1,1)
    qp%n_irr_point=0
    call lo_stop_gracefully(['Implementation of commensurate meshes not quite done.'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)

    ! ! Some basic things
    ! init: block
    !     integer :: i
    !     if ( verbosity .gt. 0 ) then
    !         write(*,*) ''
    !         write(*,*) 'Generating q-point mesh commensurate with supercell'
    !         write(*,*) 'supercellmatrix:'
    !         do i=1,3
    !             write(*,*) supercellmatrix(:,i)
    !         enddo
    !     endif
    ! end block init

    ! ! Generate the set of k-points needed
    ! genpts: block
    !     real(r8), dimension(:,:), allocatable :: dr0,dr1,dr2
    !     real(r8), dimension(3,3) :: reclat,m0
    !     real(r8), dimension(3) :: v0
    !     real(r8) :: f0,f1
    !     integer, dimension(3) :: nrep
    !     integer :: i,j,k,l,ctr1,ctr2
    !
    !     ! Get the basis vectors for the supercell reciprocal lattice
    !     reclat=matmul(uc%latticevectors,supercellmatrix)
    !     reclat=lo_reciprocal_basis(reclat)
    !
    !     ! Get the number of repetitions needed to cover the reciprocal unit cell
    !     f0=lo_bounding_sphere_of_box(uc%reciprocal_latticevectors)*2
    !     nrep=0
    !     do
    !         do i=1,3
    !             m0(:,i)=(2*nrep(i)+1)*reclat(:,i)
    !         enddo
    !         f1=lo_inscribed_sphere_in_box(m0)-1E-5_r8
    !         if ( f1 .gt. f0 ) then
    !             ! got it
    !             exit
    !         else
    !             nrep=increment_dimensions(nrep,reclat)
    !         endif
    !     enddo
    !
    !     ! repeat the supercell to get a bunch of points
    !     m0=matmul(uc%inv_reciprocal_latticevectors,reclat)
    !     ctr1=0
    !     ctr2=0
    !     do i=-nrep(1),nrep(1)
    !     do j=-nrep(2),nrep(2)
    !     do k=-nrep(3),nrep(3)
    !         v0=[i,j,k]
    !         v0=matmul(m0,v0)
    !         if ( minval(v0) .lt. -1E-10_r8 ) cycle
    !         if ( maxval(v0) .gt. 1.0_r8+1E-10_r8 ) cycle
    !         v0=lo_chop(v0,1E-11_r8)
    !         v0=lo_clean_fractional_coordinates(v0)
    !         v0=lo_chop(v0,1E-11_r8)
    !         v0=lo_clean_fractional_coordinates(v0)
    !         if ( minval(v0) .lt. lo_tiny ) then
    !             ctr1=ctr1+1
    !         else
    !             ctr2=ctr2+1
    !         endif
    !     enddo
    !     enddo
    !     enddo
    !
    !     if ( ctr1 .gt. 0 ) then
    !         allocate(dr0(3,ctr1))
    !         dr0=0.0_r8
    !     endif
    !     if ( ctr2 .gt. 0 ) then
    !         allocate(dr1(3,ctr2))
    !         dr1=0.0_r8
    !     endif
    !     ctr1=0
    !     ctr2=0
    !     do i=-nrep(1),nrep(1)
    !     do j=-nrep(2),nrep(2)
    !     do k=-nrep(3),nrep(3)
    !         v0=[i,j,k]
    !         v0=matmul(m0,v0)
    !         if ( minval(v0) .lt. -1E-10_r8 ) cycle
    !         if ( maxval(v0) .gt. 1.0_r8+1E-10_r8 ) cycle
    !
    !         v0=lo_chop(v0,1E-11_r8)
    !         v0=lo_clean_fractional_coordinates(v0)
    !         v0=lo_chop(v0,1E-11_r8)
    !         v0=lo_clean_fractional_coordinates(v0)
    !
    !         if ( minval(v0) .lt. lo_tiny ) then
    !             ctr1=ctr1+1
    !             dr0(:,ctr1)=v0
    !         else
    !             ctr2=ctr2+1
    !             dr1(:,ctr2)=v0
    !         endif
    !     enddo
    !     enddo
    !     enddo
    !
    !     ! Return the unique from the edges
    !     call lo_return_unique(dr0,dr2)
    !
    !     ! Now check that I got the correct number:
    !     l=size(dr2,2)+ctr2
    !     if ( abs(lo_determ(supercellmatrix)-l) .gt. lo_tol ) then
    !
    !         write(*,*) 'ERROR GENERATING GRID'
    !         stop
    !     endif
    !
    !     ! It was fine, store the points
    !     allocate(qvec(3,l))
    !     qvec=0.0_r8
    !     l=0
    !     do i=1,size(dr2,2)
    !         l=l+1
    !         qvec(:,l)=dr2(:,i)
    !     enddo
    !     do i=1,size(dr1,2)
    !         l=l+1
    !         qvec(:,l)=dr1(:,i)
    !     enddo
    !
    !     ! And cleanup
    !     deallocate(dr0)
    !     deallocate(dr1)
    !     deallocate(dr2)
    !
    !     if ( verbosity .gt. 0 ) then
    !         t1=walltime()
    !         write(*,*) '... generated ',tochar(size(qvec,2)),' points (',tochar(t1-t0),'s)'
    !         t0=t1
    !     endif
    ! end block genpts
    !
    !     ! Given the list of k-points, I can figure out the symmetry, I think.
    !     findsym: block
    !         real(r8), dimension(:,:), allocatable :: dr,dn
    !         integer, dimension(:,:), allocatable :: perm,shell_member
    !         integer, dimension(:), allocatable :: shell_ctr,shell_index,prototype_shell,di,dj
    !         integer :: iop,ct1,ct2,ctr_op,i,k
    !
    !         ! At least I know the total number of points
    !         qp%nq_tot = size(qvec,2)
    !
    !         ! Get a verlet box for fast locating
    !         call vb%generate(qvec,[19,19,19])
    !
    !         ! Get some space for the permutation list
    !         allocate(perm(qp%nq_tot,2*uc%sym%n))
    !         allocate(di(qp%nq_tot))
    !         allocate(dj(qp%nq_tot))
    !         allocate(dr(3,qp%nq_tot))
    !         allocate(dn(3,qp%nq_tot))
    !         perm=0
    !         di=0
    !         dj=0
    !         dr=0.0_r8
    !         dn=0.0_r8
    !         ctr_op=0
    !         do iop=1,uc%sym%n
    !             ! Reset counters
    !             ct1=0
    !             ct2=0
    !             di=0
    !             dj=0
    !             ! Rotate the points
    !             call dgemm('N','N',3,qp%nq_tot,3,1.0_r8,uc%sym%op(iop)%rfm,3,qvec,3,0.0_r8,dr,3)
    !
    !             dr=lo_clean_fractional_coordinates(lo_chop(dr,1E-10_r8))
    !             dr=lo_clean_fractional_coordinates(lo_chop(dr,1E-10_r8))
    !             ! Then we have the same operation but with time-reversal added on.
    !             dn=-dr
    !             dn=lo_clean_fractional_coordinates(lo_chop(dn,1E-10_r8))
    !             dn=lo_clean_fractional_coordinates(lo_chop(dn,1E-10_r8))
    !             ! Use the Verlet boxes to check where each point ended up.
    !             do i=1,qp%nq_tot
    !                 ! locate index of the transformed point in the original array
    !                 k=vb%locate(qvec,dr(:,i))
    !                 if ( k .gt. 0 ) then
    !                     ct1=ct1+1
    !                     di(i)=k
    !                 else
    !                     ! Make note that it's not a good operation
    !                     ct1=0
    !                 endif
    !
    !                 k=vb%locate(qvec,dn(:,i))
    !                 if ( k .gt. 0 ) then
    !                     ct2=ct2+1
    !                     dj(i)=k
    !                 else
    !                     ! Make note that it's not a good operation
    !                     ct1=0
    !                 endif
    !             enddo
    ! write(*,*) 'iop',iop,ct1,qp%nq_tot
    !             ! Check if the operation was ok, and in that case store it.
    !             if ( ct1 .eq. qp%nq_tot ) then
    !                 ctr_op=ctr_op+1
    !                 perm(:,ctr_op)=di
    !             endif
    !             if ( ct2 .eq. qp%nq_tot ) then
    !                 ctr_op=ctr_op+1
    !                 perm(:,ctr_op)=dj
    !             endif
    !         enddo
    !
    ! write(*,*) 'ctr_operation',ctr_op
    !
    !         if ( verbosity .gt. 0 ) then
    !             t1=walltime()
    !             write(*,*) '... built permutation list (',tochar(t1-t0),')'
    !             t0=t1
    !         endif
    !
    !         ! With the permutation list I can slice the points into shells and stuff.
    !         call coordination_shells_from_permutation_list(perm(:,1:ctr_op),shell_ctr,shell_member,shell_index,prototype_shell) !,mw=mw)
    !
    !         ! Then store some information
    !         qp%nq_irr = size(prototype_shell)
    !         allocate(qp%ip(qp%nq_irr))
    !         allocate(qp%ap(qp%nq_tot))
    !         do i=1,qp%nq_tot
    !             qp%ap(i)%v=qvec(:,i)
    !             qp%ap(i)%w=0.0_r8
    !             qp%ap(i)%noperations=0
    !             qp%ap(i)%irrind=shell_index(i)
    !             qp%ap(i)%operation=0
    !             qp%ap(i)%weight=0.0_r8
    !         enddo
    !         do i=1,qp%nq_irr
    !             qp%ip(i)%gridind=prototype_shell(i)
    !             qp%ip(i)%v=qvec(:,qp%ip(i)%gridind)
    !             qp%ip(i)%w=0.0_r8
    !             qp%ip(i)%weight=0.0_r8
    !             qp%ip(i)%ngridpoints=shell_ctr(i)
    !             allocate(qp%ip(i)%gridpointind( qp%ip(i)%ngridpoints ))
    !             allocate(qp%ip(i)%gridpointoperation( qp%ip(i)%ngridpoints ))
    !             qp%ip(i)%gridpointind=shell_member(1:shell_ctr(i),i)
    !             qp%ip(i)%gridpointoperation=0
    !         enddo
    !
    !         if ( verbosity .gt. 0 ) then
    !             t1=walltime()
    !             write(*,*) '... found ',tochar(qp%nq_irr),' irreducible points out of ',tochar(qp%nq_tot),' (',tochar(t1-t0),')'
    !             t0=t1
    !         endif
    !
    !
    !     end block findsym
    !
    !
    !     write(*,*) 'DONE HERE'
    !     stop
end subroutine

!> return a sensible supercell matrix
module subroutine get_supercell_matrix(uc,griddensity,supercellmatrix)
    !> unit cell
    type(lo_crystalstructure), intent(in) :: uc
    !> desired grid density
    integer, dimension(3), intent(in) :: griddensity
    !> supercellmatrix
    real(r8), dimension(3,3), intent(out) :: supercellmatrix

    real(r8), parameter :: perfect_fill=0.523598775598299_r8
    real(r8), parameter :: sphvolf=4.0_r8*lo_pi/3.0_r8
    integer, parameter :: nrep=2

    real(r8), dimension(3,3) :: m0,m1,m2
    real(r8) :: det,f0,f1,f2,fillratio,r0,volfactor
    integer, dimension(3,3) :: guessmatrix
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9

    ! Some parameters. First the ideal fill-ratio of a unit sphere in a unit cube.
    ! if I achieve this, I can cancel the loop or something
    volfactor=product(griddensity)

    ! First get a guess to start searching from:
    f0=(uc%volume*product(griddensity))**(1.0_r8/3.0_r8)
    guessmatrix=int(anint(uc%inv_latticevectors*f0))

    ! ! Now that I have an approximate supercell matrix, search a bit
    ! ! around the guess to find a really good one
    fillratio=lo_huge                   ! reset the fillratio
    m2=0.0_r8                         ! reset the guess matrix
    do i1=nrep,-nrep,-1
    do i2=nrep,-nrep,-1
    do i3=nrep,-nrep,-1
    do i4=nrep,-nrep,-1
    do i5=nrep,-nrep,-1
    do i6=nrep,-nrep,-1
    do i7=nrep,-nrep,-1
    do i8=nrep,-nrep,-1
    do i9=nrep,-nrep,-1
        m0(:,1)=[i1,i2,i3]
        m0(:,2)=[i4,i5,i6]
        m0(:,3)=[i7,i8,i9]
        m0=m0+guessmatrix
        m1=matmul(uc%latticevectors,m0)    ! apply the supercell matrix
        det=lo_determ(m1)                  ! get the volume of the new cell
        if ( det .lt. lo_tol ) cycle
        r0=lo_inscribed_sphere_in_box(m1) ! radius of inscribed sphere
        f0=sphvolf*r0**3/det              ! filling ratio
        f0=perfect_fill-f0                     ! zero means perfect filling
        f1=(det/uc%volume-volfactor)/volfactor ! zero means perfect number of atoms
        f1=f1*1E-2_r8
        ! Convert this to one number
        f2=f0**2+f1**2
        if ( f2 .lt. fillratio ) then
            fillratio=f2
            m2=m0
        endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    ! And return the matrix we need!
    supercellmatrix=m2
end subroutine

function increment_dimensions(dimin, box) result(dimout)
    integer, dimension(3), intent(in) :: dimin
    real(r8), dimension(3,3), intent(in) :: box
    integer, dimension(3) :: dimout
    !
    real(r8), dimension(3,3) :: m0
    integer, dimension(3) :: di
    integer :: i,j,k
    real(r8) :: f0,f1,ff0

    ! try with the sphere thing. First get a baseline
    do j=1,3
        m0(:,j)=box(:,j)*(2*dimin(j)+1)
    enddo
    ff0=lo_inscribed_sphere_in_box(m0)
    ! Increment the dimension that gives the biggest increase in radii
    f0=0.0_r8
    dimout=0
    do i=1,3
        di=dimin
        di(i)=di(i)+1
        do j=1,3
            m0(:,j)=box(:,j)*(2*di(j)+1)
        enddo
        f1=lo_inscribed_sphere_in_box(m0)
        if ( f1 .gt. f0 .and. abs(f1-ff0) .gt. lo_tol ) then
            dimout=di
            f0=f1
        endif
    enddo

    ! if nothing helped, increment the lowest number
    if ( dimout(1) .eq. 0 ) then
        j=lo_hugeint
        k=0
        do i=1,3
            if ( di(i) .lt. j ) then
                j=di(i)
                k=i
            endif
        enddo
        dimout=dimin
        dimout(k)=dimout(k)+1
    endif
end function

end submodule
