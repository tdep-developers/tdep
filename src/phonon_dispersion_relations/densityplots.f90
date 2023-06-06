module densityplots
use konstanter, only: r8, lo_sqtol, lo_frequency_Hartree_to_meV, lo_twopi, lo_frequency_THz_to_Hartree
use gottochblandat, only: lo_sqnorm, lo_gauss, walltime, open_file, lo_linspace, lo_progressbar_init, &
                          lo_progressbar, tochar, lo_planck
use options, only: lo_opts
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_integration_weights_for_one_tetrahedron
use type_phonon_dispersions, only: lo_phonon_dispersions, lo_phonon_dispersions_qpoint
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
implicit none

private
!public :: omega_vs_norm_q_density
public :: inelastic_spectra

interface
    module subroutine inelastic_spectra(bs, uc, fc, qp, dr, mw, mem, verbosity)
        !> bandstructure
        type(lo_phonon_bandstructure), intent(inout) :: bs
        !> unitcell
        type(lo_crystalstructure), intent(in) :: uc
        !> forceconstant
        type(lo_forceconstant_secondorder), intent(inout) :: fc
        !> q-mesh
        class(lo_qpoint_mesh), intent(in) :: qp
        !> phonon dispersions on a grid
        type(lo_phonon_dispersions), intent(in) :: dr
        !> MPI helper
        type(lo_mpi_helper), intent(inout) :: mw
        !> memory tracker
        type(lo_mem_helper), intent(inout) :: mem
        !> Talk a lot?
        integer, intent(in) :: verbosity
    end subroutine
end interface

contains

!#include "densityplots_stuntscattering.f90"

! !> Calculate a norm(q) vs omega density plot
! subroutine omega_vs_norm_q_density(qp,dr,p)
!     !> q-points
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> phonon dispersions
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> structure
!     type(lo_crystalstructure), intent(in) :: p
!
!     real(r8), dimension(:,:,:), allocatable :: tetvecs,tetvals
!     real(r8), dimension(:,:), allocatable :: tetctr,tetref,tetqnorm,hist
!
!     fixtet: block
!         real(r8), dimension(3) :: v0,v1,vref,vrefcart,vctr
!         integer :: t,i,j
!
!         ! Get the center of each tetrahedron
!         lo_allocate(tetvecs(3,4,qp%ntet))
!         lo_allocate(tetctr(3,qp%ntet))
!         lo_allocate(tetref(3,qp%ntet))
!         tetvecs=0.0_r8
!         tetctr=0.0_r8
!         tetref=0.0_r8
!         do t=1,qp%ntet
!             ! reference for this tetrahedron
!             vrefcart=qp%ap( qp%tet(t)%gridind(1) )%v
!             vctr=0.0_r8
!             vref=p%cartesian_to_fractional( vrefcart ,reciprocal=.true.)
!             do i=1,4
!                 v1=p%cartesian_to_fractional( qp%ap( qp%tet(t)%gridind(i) )%v , reciprocal=.true. )-vref
!                 v1=p%displacement_fractional_to_cartesian(v1,reciprocal=.true.)
!                 tetvecs(:,i,t)=v1
!                 vctr=vctr+v1*0.25_r8
!             enddo
!             tetctr(:,t)=vctr+vrefcart
!             tetref(:,t)=vrefcart
!             ! And shift it to the first BZ
!             if ( lo_sqnorm(tetctr(:,t)) .gt. p%bz%rmin**2 ) then
!                 v0=p%bz%gshift( tetctr(:,t) )
!                 tetctr(:,t)=tetctr(:,t)-v0
!                 tetref(:,t)=tetref(:,t)-v0
!             endif
!         enddo
!         ! Get frequencies at tetrahedron corners, as well as q-norms at the corners
!         lo_allocate(tetvals(4,qp%ntet,dr%nb))
!         lo_allocate(tetqnorm(4,qp%ntet))
!         tetvals=0.0_r8
!         tetqnorm=0.0_r8
!         do t=1,qp%ntet
!             do i=1,4
!                 do j=1,dr%nb
!                     tetvals(i,t,j)=dr%iq( qp%tet(t)%irrind(i) )%omega(j)
!                 enddo
!                 tetqnorm(i,t)=norm2( qp%ap( qp%tet(t)%gridind(i) )%w )
!             enddo
!         enddo
!     end block fixtet
!
!     ! start slicing, or something
!     slice: block
!         real(r8), dimension(:,:), allocatable :: tetweight_qnorm,hbuf,kernel
!         real(r8), dimension(:), allocatable :: qnorm,omega
!         real(r8), dimension(4) :: tete,wt,hwt
!         real(r8), dimension(3) :: v0
!         real(r8) :: f0,t0,tetsph,qnormmax,sigma
!         real(r8) :: mine,maxe,invf
!         integer, dimension(:), allocatable :: tetind
!         integer :: t,q,i,j,k,u
!         integer :: ctr,nq,nf,ii,jj,nk
!
!         nq=400
!         nf=600
!         sigma=dr%default_smearing()
!
!         ! radius of the bounding sphere of the largest tetrahedron
!         tetsph=0.0_r8
!         do t=1,qp%ntet
!             v0=0.0_r8
!             do i=1,4
!                 v0=v0+tetvecs(:,i,t)*0.25_r8
!             enddo
!             do i=1,4
!                 tetsph=max(tetsph,norm2(tetvecs(:,i,t)))
!             enddo
!         enddo
!         ! longest possible q
!         qnormmax=0.0_r8
!         do t=1,qp%ntet
!             qnormmax=max(qnormmax,norm2(tetctr(:,t)))
!         enddo
!         qnormmax=qnormmax+tetsph
!
!         ! q-norm axis
!         lo_allocate(qnorm(nq))
!         call lo_linspace(0.0_r8,qnormmax,qnorm)
!         ! omega-axis
!         lo_allocate(omega(nf))
!         call lo_linspace(0.0_r8,dr%omega_max+dr%default_smearing(),omega)
!         invf=(1.0_r8*nf)/maxval(omega)
!         ! Space for the q-weights
!         lo_allocate(tetweight_qnorm(4,qp%ntet))
!         lo_allocate(tetind(qp%ntet))
!         ! And the histogram
!         lo_allocate(hist(nf,nq))
!         hist=0.0_r8
!
!         t0=walltime()
!         call lo_progressbar_init()
!         ! First, let's try to calculate the volume of the BZ inside the |q|-sphere
!         do q=2,nq
!             ! Get the q-norm integral weights
!             tetweight_qnorm=0.0_r8
!             tetind=0
!             ctr=0
!             do t=1,qp%ntet
!                 f0=norm2( tetctr(:,t) )
!                 ! First make an easy decision. If the tetrahedron is completely inside the sphere, it's easy!
!                 if ( f0+tetsph .lt. qnorm(q) ) then
!                     ! completely inside
!                     !tetweight_qnorm(:,t)=qp%tet(t)%weight*0.25_r8
!                 elseif ( f0-tetsph .gt. qnorm(q) ) then
!                     ! completely outside
!                     !tetweight_qnorm(:,t)=0.0_r8
!                 else
!                     ! Get the weights
!                     call lo_integration_weights_for_one_tetrahedron(qp%tet(t),tetqnorm(:,t),qnorm(q),wt,&
!                                                                     lo_sqtol,.false.,lo_sqtol,hwt)
!                     !tetweight_qnorm(:,t)=wt
!                     if ( sum(wt) .gt. lo_sqtol ) then
!                         ctr=ctr+1
!                         tetind(ctr)=t
!                         tetweight_qnorm(:,ctr)=wt
!                     endif
!                 endif
!             enddo
!
!             ! Integrate over omega
!             do j=1,dr%nb
!             do k=1,ctr
!                 t=tetind(k)
!                 tete=tetvals(:,t,j)
!                 mine=minval(tete)
!                 maxe=maxval(tete)
!                 ii=max(floor( mine*invf ),2)
!                 jj=min(ceiling( maxe*invf ),nf)
!                 do i=ii,jj
!                     call lo_integration_weights_for_one_tetrahedron(qp%tet(t),tete,omega(i),wt,sigma,&
!                                                                     blochlcorrections=.false.,hweights=hwt)
!                     hist(i,q)=hist(i,q)+sum(wt*tetweight_qnorm(:,k))
!                 enddo
!             enddo
!             enddo
!             call lo_progressbar(' ... building density',q,nq,walltime()-t0)
!         enddo
!
!         ! Try the magic negative indexing in fortran, do a slight smearing.
!         sigma=1.25_r8
!         nk=3
!         lo_allocate(hbuf(-nk:nf+nk,-nk:nq+nk))
!         lo_allocate(kernel(-nk:nk,-nk:nk))
!         do ii=-nk,nk
!         do jj=-nk,nk
!             kernel=lo_gauss( norm2([ii*1.0_r8,jj*1.0_r8]),0.0_r8,sigma )
!         enddo
!         enddo
!         hbuf=0.0_r8
!         hbuf(1:nf,1:nq)=hist
!         hist=0.0_r8
!         call lo_progressbar_init()
!         do i=1,nf
!             do j=1,nq
!                 do ii=-nk,nk
!                 do jj=-nk,nk
!                     hist(i,j)=hist(i,j)+hbuf(i+ii,j+jj)*kernel(ii,jj)
!                 enddo
!                 enddo
!             enddo
!             call lo_progressbar(' ... smearing a little',i,nf)
!         enddo
!         !
!         u=open_file('out','dumdens')
!             do i=1,nq
!             do j=1,nf
!                 write(u,*) qnorm(i),omega(j),hist(j,i)
!             enddo
!             enddo
!         close(u)
!         !
!     end block slice
!
! end subroutine

end module
