#include "precompilerdefinitions"
module unfold_phonons
use konstanter, only: flyt, lo_twopi
use gottochblandat, only: walltime, lo_chop, lo_clean_fractional_coordinates, tochar, &
                          lo_linspace, lo_gauss, lo_trapezoid_integration, &
                          lo_progressbar_init, lo_progressbar, lo_complex_hermitian_eigenvalues_eigenvectors, lo_negsqrt, &
                          lo_signed_tetrahedron_volume
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use lo_randomnumbers, only: lo_mersennetwister
use type_forceconstant_secondorder
use type_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint, lo_fft_mesh, lo_qpoint_mesh, lo_generate_qmesh
use type_phonon_dispersions
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use options, only: lo_opts
use type_fast_interpolation, only: lo_fancy_deltri_box
use type_forcemap, only: lo_forcemap
use type_sqs, only: lo_sqs
implicit none

private
!public :: unfold_bandstructure
!public :: unfold_brutal

contains

! subroutine unfold_brutal(uc,fc,bs,mw,opts,verbosity)
!     !> crystal structure
!     type(lo_crystalstructure), intent(inout) :: uc
!     !> forceconstant
!     type(lo_forceconstant_secondorder), intent(inout) :: fc
!     !> bandstructure
!     type(lo_phonon_bandstructure), intent(inout) :: bs
!     !> mpi helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> Options
!     type(lo_opts), intent(in) :: opts
!     !> Talk a lot?
!     integer, intent(in) :: verbosity
!
!     ! Common things
!     type(lo_fancy_deltri_box) :: box
!     type(lo_mersennetwister) :: tw
!     type(lo_crystalstructure) :: ss
!     type(lo_forceconstant_secondorder) :: fcss
!     class(lo_qpoint_mesh), allocatable :: qgu,qgs
!     type(lo_phonon_dispersions) :: dr
!     real(flyt), dimension(:,:,:), allocatable :: hist
!     real(flyt), dimension(:,:), allocatable :: Gvecs
!     real(flyt) :: timer,maxomega
!     integer, dimension(:,:,:), allocatable :: rmind
!     integer, dimension(:,:), allocatable :: histctr
!     integer :: nGvecs,ne
!     !type(lo_fft_mesh) :: quc,qss
!
!     init: block
!         real(flyt), dimension(3,3) :: I3
!         real(flyt), dimension(3) :: v0
!         integer, dimension(3) :: qdim
!         integer :: i,j,k,l
!
!         ! What is the unitcell q-grid?
!         qdim=opts%ssdim*opts%qgrid
!
!         if ( verbosity .gt. 0 ) then
!             timer=walltime()
!             write(*,*) ''
!             write(*,*) 'Creating unfolded bandstructure'
!             write(*,*) ' supercell dimensions: ',tochar(opts%ssdim)
!             write(*,*) '     supercell q-grid: ',tochar(opts%qgrid)
!             write(*,*) '      unitcell q-grid: ',tochar(qdim)
!         endif
!
!         ! Build supercell
!         call uc%build_supercell(ss,opts%ssdim)
!         call ss%classify('supercell',uc)
!         ! Get supercell forceconstant
!         l=0
!         do i=1,fc%na
!             l=max(l,fc%atom(i)%n)
!         enddo
!         lo_allocate(rmind(3,l,ss%na))
!         rmind=0
!         call fc%remap(uc,ss,fcss,ind=rmind)
!
!         ! Create a fake spacegroup for the supercell
!         I3=0.0_flyt
!         do i=1,3
!             I3(i,i)=1.0_flyt
!         enddo
!         ss%sym%n=1
!         lo_allocate(ss%sym%op(1))
!         ss%sym%op(1)%m=I3
!         ss%sym%op(1)%im=I3
!         ss%sym%op(1)%fm=I3
!         ss%sym%op(1)%rfm=I3
!         ss%sym%op(1)%irfm=I3
!         ss%sym%op(1)%tr=0.0_flyt
!         ss%sym%op(1)%ftr=0.0_flyt
!         lo_allocate(ss%sym%op(1)%fmap(ss%na))
!         do i=1,ss%na
!             ss%sym%op(1)%fmap(i)=i
!         enddo
!         ss%sym%timereversal=.true.
!         ss%info%havespacegroup=.true.
!         ss%info%decidedtimereversal=.true.
!         call ss%classify('supercell',uc)
!         call ss%classify('wedge',timereversal=.true.)
!
!         ! Generate the relevant q-meshes
!         call lo_generate_qmesh(qgu,uc,qdim,'fft',verbosity=opts%verbosity,timereversal=.true.)
!         call lo_generate_qmesh(qgs,ss,opts%qgrid,'fft',verbosity=opts%verbosity,timereversal=.true.)
!
!         ! And the supercell G-vectors
!         lo_allocate(Gvecs(3,product(opts%ssdim)))
!         l=0
!         do i=1,opts%ssdim(1)
!         do j=1,opts%ssdim(2)
!         do k=1,opts%ssdim(3)
!             l=l+1
!             v0=[i,j,k]-1
!             Gvecs(:,l)=lo_chop( matmul(ss%reciprocal_latticevectors,v0), 1E-11_flyt )
!         enddo
!         enddo
!         enddo
!         nGvecs=l
!
!         ! And space for the histograms
!         ne=opts%nenergy
!         lo_allocate(hist(ne,uc%na*3,qgu%nq_irr))
!         hist=0.0_flyt
!         ! Set the energy axis?
!         lo_allocate(bs%faxis(ne))
!         maxomega=0.0_flyt
!         do i=1,bs%nptot
!             maxomega=max(maxomega,maxval(bs%p(i)%omega))
!         enddo
!         call lo_linspace(0.0_flyt,maxomega*2,bs%faxis)
!         ! Counter for the histogram
!         lo_allocate(histctr(uc%na*3,qgu%nq_irr))
!         histctr=0.0_flyt
!         ! Get the empty intensity
!         lo_allocate(bs%intensity(bs%nptot,ne))
!         bs%intensity=0.0_flyt
!
!         ! Generate the unitcell dispersions at least?
!         call dr%generate(qgu,fc,uc,verbosity=opts%verbosity,timereversal=.true.,mpi_communicator=mw%comm)
!
!         if ( verbosity .gt. 0 ) write(*,*) ' ... created grids and unitcell dispersions (',tochar(walltime()-timer),'s)'
!     end block init
!
!     fold1: block
!         complex(flyt), dimension(:,:), allocatable :: transoperator
!         complex(flyt), dimension(:,:), allocatable :: dynmat,eigenvectors
!         complex(flyt), dimension(:), allocatable :: cv0
!         complex(flyt) :: expiqr,c0
!         real(flyt), dimension(:), allocatable :: omega
!         real(flyt), dimension(3) :: v0,v1,v2
!         real(flyt) :: qdotr,invf,foursigma,sigma,f0,f1,t0
!         integer, dimension(3) :: gi
!         integer :: qi,qj,ig,si,irrind,mode
!         integer :: i,j,ii,jj,ctr
!         ! Figure out the folding map somehow
!         t0=walltime()
!         ! Some temporary space
!         lo_allocate(dynmat(ss%na*3,ss%na*3))
!         lo_allocate(eigenvectors(ss%na*3,ss%na*3))
!         lo_allocate(omega(ss%na*3))
!         lo_allocate(transoperator(ss%na*3,uc%na*3))
!         lo_allocate(cv0(ss%na*3))
!
!         ! Some shorthand for meaningful stuff
!         sigma=(bs%faxis(2)-bs%faxis(1))*5
!         foursigma=6*sigma
!         invf=real(ne,flyt)/bs%faxis(ne)
!
!         if ( mw%talk ) call lo_progressbar_init()
!
!         ctr=0
!         do si=1,1
!         ! Insert randomizer here
!         call randomize_masses(ss,tw)
!         do qi=1,qgs%nq_tot
!             ctr=ctr+1
!             ! make parallel
!             if ( mod(qi,mw%n) .ne. mw%r ) cycle
!             ! Solve the supercell problem
!             call fcss%dynamicalmatrix(ss,qgs%ap(qi),dynmat)
!             call lo_complex_hermitian_eigenvalues_eigenvectors(dynmat,omega,eigenvectors,careful=.false.,tolerance=1E-12_flyt)
!             omega=lo_negsqrt(omega)
!             ! Build
!
!             ! Now unfold this point somehow. What small q-points does it unfold to?
!             c0=0.0_flyt
!             do ig=1,nGvecs
!                 v0=qgs%ap(qi)%w+Gvecs(:,ig)
!                 v1=matmul(uc%inv_reciprocal_latticevectors,v0)
!                 select type(qgu); type is(lo_fft_mesh)
!                     gi=qgu%index_from_coordinate(v1)
!                     qj=qgu%gridind2ind(gi(1),gi(2),gi(3))
!                     irrind=qgu%ap(qj)%irrind
!                 end select
!                 ! Ok now I know where to fold it, I suppose. Build the translation operator?
!                 transoperator=0.0_flyt
!                 do mode=1,dr%nb
!                 do i=1,ss%na
!                     j=ss%info%index_in_unitcell(i)
!                     v2=ss%info%cellindex(:,i)-1
!                     v2=matmul(uc%latticevectors,v2)
!                     qdotr=-dot_product(v0,v2)*lo_twopi
!                     expiqr=cmplx(cos(qdotr),sin(qdotr),flyt)
!                     transoperator( (i-1)*3+1:i*3,mode )=dr%aq(qj)%egv( (j-1)*3+1:j*3, mode )*expiqr
!                 enddo
!                 enddo
!                 transoperator=transoperator/sqrt(real(nGvecs,flyt))
!                 ! Bin it in the right place
!                 do mode=1,dr%nb
!                     do i=1,ss%na*3
!                         !cv0=transoperator(:,mode)*eigenvectors(:,i)
!                         !f0=real(dot_product(cv0,cv0)) ! Spectral weight
!                         c0=dot_product(transoperator(:,mode),eigenvectors(:,i))
!                         f0=abs(conjg(c0)*c0)
!                         f1=omega(i)                   ! Energy
!                         ii=max(floor((f1-foursigma)*invf),1)
!                         jj=min(ceiling((f1+foursigma)*invf),ne)
!                         do j=ii,jj
!                             hist(j,mode,irrind)=hist(j,mode,irrind)+lo_gauss( bs%faxis(j),f1,sigma )*f0
!                         enddo
!                     enddo
!                     histctr(mode,irrind)=histctr(mode,irrind)+1
!                 enddo
!             enddo
!
!             if ( mw%talk .and. qi .lt. qgs%nq_tot ) then
!                 call lo_progressbar(' ... frequencies',ctr,qgs%nq_tot*opts%nsample,walltime()-t0)
!             endif
!         enddo
!         enddo
!
!         call mw%allreduce('sum',hist)
!         call mw%allreduce('sum',histctr)
!         do i=1,qgu%nq_irr
!         do j=1,dr%nb
!             hist(:,j,i)=hist(:,j,i)/histctr(j,i)
!         enddo
!         enddo
!         if ( mw%talk ) then
!             call lo_progressbar(' ... frequencies',qgs%nq_tot*opts%nsample,qgs%nq_tot*opts%nsample,walltime()-t0)
!         endif
!
!     end block fold1
!
!     specfun: block
!         real(flyt), dimension(:,:), allocatable :: y1
!         real(flyt), dimension(:), allocatable :: y0
!         real(flyt), dimension(3,4) :: tet
!         real(flyt), dimension(3) :: qv,v0,v1,v2,v3,r,wts
!         real(flyt) :: f0,t0
!         integer, dimension(4) :: gi
!         integer :: q,ti,i,j,k,mode
!
!         t0=walltime()
!         ! First build the locator thingy
!         call box%generate(qgu,uc)
!         ! Some temp space
!         lo_allocate(y0(ne))
!         lo_allocate(y1(3,ne))
!
!         if ( mw%talk ) call lo_progressbar_init()
!
!         do q=1,bs%nptot
!             ! Make MPI parallel
!             if ( mod(q,mw%n) .ne. mw%r ) cycle
!             ! Grab tetrahedron
!             qv=lo_clean_fractional_coordinates( matmul(uc%inv_reciprocal_latticevectors,bs%q(q)%w) )
!             qv=matmul(uc%reciprocal_latticevectors,qv)
!             ti=box%tetind(qv)
!             select type(qgu)
!             type is(lo_fft_mesh)
!                 v0=matmul(uc%inv_reciprocal_latticevectors,qv)
!                 do j=1,4
!                     gi(j)=qgu%tet(ti)%irrind(j)
!                     tet(:,j)=matmul(uc%inv_reciprocal_latticevectors, qgu%ap( qgu%tet(ti)%gridind(j) )%v )
!                 enddo
!                 ! pbc-adjust
!                 do j=1,4
!                     do k=1,3
!                         if ( tet(k,j)-v0(k) .lt. -0.5_flyt ) tet(k,j)=tet(k,j)+1.0_flyt
!                         if ( tet(k,j)-v0(k) .gt. 0.5_flyt ) tet(k,j)=tet(k,j)-1.0_flyt
!                     enddo
!                     tet(:,j)=matmul(uc%reciprocal_latticevectors, tet(:,j))
!                 enddo
!             end select
!             ! Build contragradient
!             f0=-lo_signed_tetrahedron_volume(tet)*6
!             v0=tet(:,1)
!             v1=tet(:,2)-v0
!             v2=tet(:,3)-v0
!             v3=tet(:,4)-v0
!             r=qv-v0
!             wts=ewts(v1,v2,v3,r,f0)
!             ! Start interpolating
!             do mode=1,bs%nb
!                 ! Fetch values
!                 y0=hist(:,mode,gi(1))
!                 do i=1,3
!                     y1(i,:)=hist(:,mode,gi(i+1))
!                 enddo
!                 ! Begin operation actual operation
!                 do i=1,ne
!                     y0(i)=y0(i)+dot_product(wts,y1(:,i))
!                 enddo
!                 ! Store intensiry
!                 bs%intensity(q,:)=bs%intensity(q,:)+y0
!             enddo
!
!             if ( mw%talk .and. q .lt. bs%nptot ) then
!                 call lo_progressbar(' ... spectralfunction',q,bs%nptot,walltime()-t0)
!             endif
!         enddo
!         call mw%allreduce('sum',bs%intensity)
!         if ( mw%talk ) call lo_progressbar(' ... spectralfunction',bs%nptot,bs%nptot,walltime()-t0)
!
!         ! Dump to file?
!         if ( mw%talk ) then
!             call bs%write_intensity('thz','outfile.unfolded_sqe')
!         endif
!     end block specfun
! end subroutine
!
! ! tetrahedral interpolation weights
! function ewts(v1,v2,v3,r,vol) result(wt)
!     real(flyt), dimension(3), intent(in) :: v1,v2,v3,r
!     real(flyt), intent(in) :: vol
!     real(flyt), dimension(3) :: wt
!
!     real(flyt) :: rx,ry,rz
!     real(flyt) :: v1x,v1y,v1z
!     real(flyt) :: v2x,v2y,v2z
!     real(flyt) :: v3x,v3y,v3z
!
!     v1x=v1(1); v1y=v1(2); v1z=v1(3);
!     v2x=v2(1); v2y=v2(2); v2z=v2(3);
!     v3x=v3(1); v3y=v3(2); v3z=v3(3);
!     rx=r(1);   ry=r(2);   rz=r(3);
!
!     wt(1)=(-(rz*v2y*v3x) + ry*v2z*v3x + rz*v2x*v3y - rx*v2z*v3y - ry*v2x*v3z + rx*v2y*v3z)/vol
!     wt(2)=(rz*v1y*v3x - ry*v1z*v3x - rz*v1x*v3y + rx*v1z*v3y + ry*v1x*v3z - rx*v1y*v3z)/vol
!     wt(3)=(-(rz*v1y*v2x) + ry*v1z*v2x + rz*v1x*v2y - rx*v1z*v2y - ry*v1x*v2z + rx*v1y*v2z)/vol
! end function
!
! subroutine unfold_bandstructure(uc,fc,bs,mw,opts,verbosity)
!     !> crystal structure
!     type(lo_crystalstructure), intent(inout) :: uc
!     !> forceconstant
!     type(lo_forceconstant_secondorder), intent(inout) :: fc
!     !> bandstructure
!     type(lo_phonon_bandstructure), intent(inout) :: bs
!     !> mpi helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> Options
!     type(lo_opts), intent(in) :: opts
!     !> Talk a lot?
!     integer, intent(in) :: verbosity
!
!     type(lo_sqs) :: sqs
!     type(lo_mersennetwister) :: tw
!     type(lo_crystalstructure) :: ss
!     type(lo_forceconstant_secondorder) :: fcss
!     integer, dimension(:,:,:), allocatable :: rmind
!     integer, dimension(3) :: ssdim
!     integer :: nqss,ne
!     real(flyt) :: timer,maxomega
!
!     init: block
!         integer :: i,l
!
!         if ( verbosity .gt. 0 ) then
!             timer=walltime()
!             write(*,*) ''
!             write(*,*) 'Creating unfolded bandstructure'
!         endif
!
!         ! Init random numbers
!         call tw%init(iseed=mw%r,rseed=walltime())
!
!         ! Set the supercell size
!         ssdim=opts%ssdim
!         nqss=product(ssdim)
!         ! Build supercell
!         call uc%build_supercell(ss,ssdim)
!         call ss%classify('supercell',uc)
!         ! Get supercell forceconstant
!         l=0
!         do i=1,fc%na
!             l=max(l,fc%atom(i)%n)
!         enddo
!         lo_allocate(rmind(3,l,ss%na))
!         rmind=0
!         call fc%remap(uc,ss,fcss,ind=rmind)
!
!         ! Set some options
!         ne=opts%nenergy
!         ! Set the energy axis?
!         lo_allocate(bs%faxis(ne))
!         maxomega=0.0_flyt
!         do i=1,bs%nptot
!             maxomega=max(maxomega,maxval(bs%p(i)%omega))
!         enddo
!         call lo_linspace(0.0_flyt,maxomega*2,bs%faxis)
!         ! Get the empty intensity
!         lo_allocate(bs%intensity(bs%nptot,ne))
!         bs%intensity=0.0_flyt
!     end block init
!
!     ! Generate neat SQS structures!
!     buildsqs: block
!         type(lo_symlist) :: sl
!         type(lo_forcemap) :: map
!         real(flyt) :: f0
!
!         f0=13.0_flyt
!         call sl%generate(uc,ss,f0,-1.0_flyt,-1.0_flyt,verbosity=opts%verbosity,&
!                          polar=.false.,wraparound=.true.)
!         call map%generate(sl,uc,ss,polarcorrectiontype=0,verbosity=opts%verbosity)
!
!         call sqs%generate(uc,ss,sl,max(opts%nsample,5),opts%verbosity,mw)
!     end block buildsqs
!
!     ! For every unitcell q-vector, there are a bunch of supercell q-vectors that will unfold to that q-vector, I think
!     findqvecs: block
!         type(lo_qpoint) :: qpoint
!         type(lo_crystalstructure) :: dss
!         complex(flyt), dimension(:,:), allocatable :: ucegv,dynmat,eigenvectors
!         complex(flyt), dimension(:), allocatable :: cv0
!         complex(flyt) :: expiqr,c0
!         real(flyt), dimension(:,:), allocatable :: modehist,hmbas
!         real(flyt), dimension(:), allocatable :: omega
!         real(flyt), dimension(3) :: v0
!         real(flyt) :: qdotr,sigma,foursigma,invf,f0,f1
!         integer :: i,j,ii,jj
!         integer :: q,mode,sample
!
!         ! So how does the transformation work? We have a bunch of vectors:
!         !    k     unitcell k-vector
!         !    K     supercell k-vector
!         !    g     unitcell reciprocal latticevector
!         !    G     supercell reciprocal latticevector
!         ! I think I want to find vectors such that
!         !    K = G + k
!
!         sigma=(bs%faxis(2)-bs%faxis(1))*5
!         foursigma=6*sigma
!         invf=real(ne,flyt)/bs%faxis(ne)
!
!         lo_allocate(ucegv(ss%na*3,bs%nb))
!         lo_allocate(dynmat(ss%na*3,ss%na*3))
!         lo_allocate(eigenvectors(ss%na*3,ss%na*3))
!         lo_allocate(omega(ss%na*3))
!         lo_allocate(modehist(ne,bs%nb))
!         lo_allocate(cv0(ss%na*3))
!         ucegv=0.0_flyt
!         modehist=0.0_flyt
!
!         lo_allocate(hmbas(uc%na*3,uc%na*3))
!         hmbas=0.0_flyt
!         do i=1,uc%na*3
!             hmbas(i,i)=1.0_flyt
!         enddo
!
!         if ( mw%talk ) call lo_progressbar_init()
!
!         do q=1,bs%nptot
!             ! Make it MPI parallel
!             if ( mod(q,mw%n) .ne. mw%r ) cycle
!
!             ! Build the unitcell eigenvectors in the supercell so that I can project properly
!             do mode=1,bs%nb
!                 do i=1,ss%na
!                     j=ss%info%index_in_unitcell(i)
!                     v0=ss%info%cellindex(:,i)-1
!                     v0=matmul(uc%latticevectors,v0)
!                     !qdotr=-dot_product(bs%q(q)%w,ss%rcart(:,i))*lo_twopi
!                     qdotr=-dot_product(bs%q(q)%w,v0)*lo_twopi
!
!                     expiqr=cmplx(cos(qdotr),sin(qdotr),flyt)
!                     ucegv( (i-1)*3+1:i*3,mode )=bs%p(q)%egv( (j-1)*3+1:j*3, mode )*expiqr
!                     ucegv( (i-1)*3+1:i*3,mode )=ucegv( (i-1)*3+1:i*3,mode )*ss%invsqrtmass(i)
!                     !ucegv( (i-1)*3+1:i*3,mode )=hmbas( (j-1)*3+1:j*3, mode )*expiqr
!                 enddo
!             enddo
!             ucegv=ucegv/sqrt(real(nqss,flyt))
!
!             ! Randomize the masses in the supercell?
!             qpoint%v=bs%q(q)%v
!             qpoint%w=bs%q(q)%w
!             qpoint%noperations=0
!             modehist=0.0_flyt
!             do sample=1,opts%nsample
!                 !call randomize_masses(dss,.false.,fc,fcss,rmind,tw)
!                 call sqs%returnstructure(ss,dss,sample,.false.) !,opts%verbosity)
!                 call fcss%dynamicalmatrix(dss,qpoint,dynmat)
!                 call lo_complex_hermitian_eigenvalues_eigenvectors(dynmat,omega,eigenvectors,careful=.false.,tolerance=1E-12_flyt)
!                 omega=lo_negsqrt(omega)
!                 !
!                 do mode=1,ss%na*3
!                 do i=1,ss%na
!                     eigenvectors( (i-1)*3+1:i*3,mode )=eigenvectors( (i-1)*3+1:i*3,mode )*dss%invsqrtmass(i)
!                 enddo
!                 enddo
!
!                 ! Ok now I have randomized guys ... accumulate spectral function?
!                 do mode=1,bs%nb
!                     c0=0.0_flyt
!                     do i=1,ss%na*3
!                         c0=dot_product(ucegv(:,mode),eigenvectors(:,i))
!                         f0=abs(conjg(c0)*c0)
!                         f1=omega(i)
!                         ii=max(floor((f1-foursigma)*invf),1)
!                         jj=min(ceiling((f1+foursigma)*invf),ne)
!                         do j=ii,jj
!                             modehist(j,mode)=modehist(j,mode)+lo_gauss( bs%faxis(j),f1,sigma )*f0
!                         enddo
!                     enddo
!                 enddo
!             enddo
!             ! Add it up?
!             do mode=1,bs%nb
!                 bs%intensity(q,:)=bs%intensity(q,:)+modehist(:,mode)/opts%nsample
!             enddo
!             if ( mw%talk .and. q .lt. bs%nptot ) then
!                 call lo_progressbar(' ... unfolding bandstructure',q,bs%nptot,walltime()-timer)
!             endif
!         enddo
!
!         call mw%allreduce('sum',bs%intensity)
!         if ( mw%talk ) then
!             call lo_progressbar(' ... unfolding bandstructure',bs%nptot,bs%nptot,walltime()-timer)
!         endif
!
!         ! Dump to file?
!         if ( mw%talk ) then
!             call bs%write_intensity('thz','outfile.unfolded_sqe')
!         endif
!
!     end block findqvecs
! end subroutine
!
! !> Assign masses to every atom according to the proper probability
! subroutine randomize_masses(p,tw)
!     !> Structure to assign random masses
!     type(lo_crystalstructure), intent(inout) :: p
!     !> Should each sublattice get exactly the correct amount, or random?
!     !logical, intent(in) :: preserve
!     !> remapping index
!     !integer, dimension(:,:,:), intent(in) :: rmind
!     !> Random number generator
!     type(lo_mersennetwister), intent(inout) :: tw
!
!     integer :: a1,i,j,l,ii,nsp,n1,n2
!     integer, dimension(:,:), allocatable :: di
!     integer, dimension(:), allocatable :: ctr,occ
!     integer, dimension(10) :: cns
!
!     ! Number of species
!     nsp=maxval(p%species)
!     lo_allocate(ctr(nsp))
!     lo_allocate(di(p%na,nsp))
!     lo_allocate(occ(p%na))
!     ctr=0
!     di=0
!     do a1=1,p%na
!         ii=p%species(a1)
!         ctr(ii)=ctr(ii)+1
!         di( ctr(ii), ii ) = a1
!     enddo
!
!     do ii=1,nsp
!         ! Get rough distribution
!         cns=0
!         do i=1,p%alloyspecies(ii)%n
!             cns(i)=max( ceiling(p%alloyspecies(ii)%concentration(i)*ctr(ii)), 1 )
!         enddo
!         ! Make sure it makes sense
!         n1=sum(cns)
!         if ( n1 .gt. ctr(ii) ) then
!             do
!                 n1=0
!                 n2=0
!                 do j=1,p%alloyspecies(ii)%n
!                     n1= max(cns(j),n1)
!                     if ( n1 .ge. n2 ) then
!                         n2=n1
!                         i=j
!                     endif
!                 enddo
!                 cns(i)=cns(i)-1
!                 if ( sum(cns) .eq. ctr(ii) ) exit
!             enddo
!         endif
!         ! Shuffle positions
!         call tw%shuffle_int_array(di( 1:ctr(ii),ii))
!         l=0
!         do i=1,p%alloyspecies(ii)%n
!         do j=1,cns(i)
!             l=l+1
!             a1=di(l,ii)
!             occ(a1)=i
!             p%mass( a1 )=p%alloyspecies(ii)%mass(i)
!             p%invsqrtmass( a1 )=1.0_flyt/sqrt(p%mass(a1))
!         enddo
!         enddo
!     enddo
!
! !    ! Ok now masses are set, fiddle a bit with the forceconstants
! !    do a1=1,fcss%na
! !        ii=rmind(1,1,a1)
! !        do i=1,fcss%atom(a1)%n
! !            a2=fcss%atom(a1)%pair(i)%i2
! !            if ( occ(a1) .eq. occ(a2) ) then
! !                if ( occ(a1) .eq. 1 ) then
! !                fcss%atom(a1)%pair(i)%m = fc%atom(ii)%pair(i)%m*1.2_flyt
! !                else
! !                fcss%atom(a1)%pair(i)%m = fc%atom(ii)%pair(i)%m*0.8_flyt
! !                endif
! !            else
! !                fcss%atom(a1)%pair(i)%m = fc%atom(ii)%pair(i)%m*0.5_flyt
! !            endif
! !        enddo
! !    enddo
! !    call fcss%setsumtozero()
!
! !stop
! !
! !    do a1=1,p%na
! !        ii=p%species(a1)
! !        f0=tw%rnd_real()
! !        do i=1,p%alloyspecies(ii)%n
! !            f1=sum(p%alloyspecies(ii)%concentration(1:i))
! !            if ( f1 .ge. f0 ) exit
! !        enddo
! !        p%mass(a1)=p%alloyspecies(ii)%mass(i)
! !        p%invsqrtmass(a1)=1.0_flyt/sqrt(p%mass(a1))
! !    enddo
!
! end subroutine

! Return the set of supercell K-vectors
!subroutine kvecs_from_k(kv,gvecs)
!end subroutine

end module
