
! !> Get the spectral anharmonic free energy
! subroutine spectral_free_energy(qp,dr,uc,F3,F4,npts,unitfactor,mw)
!     !> qpoint mesh
!     class(lo_qpoint_mesh), intent(inout) :: qp
!     !> phonon dispersions
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> crystal structure
!     type(lo_crystalstructure), intent(in) :: uc
!     !> anharmonic free energy
!     real(r8), dimension(:,:), intent(in) :: F3,F4
!     !> number of points on the x-axis
!     integer, intent(in) :: npts
!     !> output units
!     real(r8), intent(in) :: unitfactor
!     !> MPI helper
!     type(lo_mpi_helper), intent(in) :: mw
!
!     !
!     complex(r8), dimension(3) :: v0
!     real(r8), dimension(:,:,:), allocatable :: tetvals,tetF3,tetF4,buf3,buf4,buf3b,buf4b,tetsiteproj
!     real(r8), dimension(:,:), allocatable :: yp3,yp4
!     real(r8), dimension(:), allocatable :: x,yt3,yt4
!     real(r8), dimension(4) :: tete,tv3,tv4,w1,w2
!     real(r8), dimension(9) :: kernel
!     real(r8), dimension(uc%na) :: tetpr
!     real(r8) :: invf,sigma,deltae,t0,ftot3,ftot4,f0,f1,f2,nf3,nf4
!     integer :: q,t,i,j,k,l,ii,jj,li,u
!     character(len=1000) :: opf
!
!     ! Distribute q-points and tetrahedrons over MPI
!     t0=walltime()
!     !if ( qp%mpi%initialized .eqv. .false. ) call qp%divide_across_mpi(mw%comm)
!     ! Prefetch things for tetrahedrons.
!     lo_allocate(tetsiteproj(uc%na,dr%n_mode,qp%mpi%ntet))
!     lo_allocate(tetvals(4,dr%n_mode,qp%mpi%ntet))
!     lo_allocate(tetF3(4,dr%n_mode,qp%mpi%ntet))
!     lo_allocate(tetF4(4,dr%n_mode,qp%mpi%ntet))
!     tetvals=0.0_r8
!     tetsiteproj=0.0_r8
!     tetF3=0.0_r8
!     tetF4=0.0_r8
!     do li=1,qp%mpi%ntet
!         i=qp%mpi%tet(li)
!         do j=1,4
!             ! energies
!             l=qp%tet(i)%irrind(j)
!             do k=1,dr%n_mode
!                 tetvals(j,k,li)=dr%iq(l)%omega(k)
!                 tetF3(j,k,li)=F3(k,l)
!                 tetF4(j,k,li)=F4(k,l)
!             enddo
!             ! projections
!             do ii=1,uc%na
!             do jj=1,dr%n_mode
!                 v0=dr%iq(l)%egv((ii-1)*3+1:ii*3,jj)
!                 tetsiteproj(ii,jj,li)=tetsiteproj(ii,jj,li)+abs(dot_product(v0,v0))*0.25_r8
!             enddo
!             enddo
!         enddo
!     enddo
!
!     ! Calculate it by energy, mode, site. Make some space
!     lo_allocate(buf3(uc%na,npts,dr%n_mode))
!     lo_allocate(buf4(uc%na,npts,dr%n_mode))
!     buf3=0.0_r8
!     buf4=0.0_r8
!     ! Get the x-axis
!     lo_allocate(x(npts))
!     call lo_linspace(0.0_r8,dr%omega_max*1.1_r8,x)
!     invf=npts/x(npts)
!     sigma=dr%default_smearing()
!     deltae=(x(2)-x(1))*0.25_r8
!
!     ! Now do the actual integration
!     if ( mw%talk ) call lo_progressbar_init()
!     do li=1,qp%mpi%ntet
!         i=qp%mpi%tet(li)
!         do k=1,dr%n_mode
!             ! prefetch some stuff
!             tete=tetvals(:,k,li)
!             tv3=tetF3(:,k,li)
!             tv4=tetF4(:,k,li)
!             tetpr=tetsiteproj(:,k,li)
!             ii=max(floor( (minval(tete)-2*deltae)*invf ),1)
!             jj=min(ceiling( (maxval(tete)+2*deltae)*invf ),npts)
!             ! Add it up
!             do l=ii,jj
!                 call lo_integration_weights_for_one_tetrahedron(qp%tet(i),tete,x(l)-deltae,w1,sigma,blochlcorrections=.false.)
!                 call lo_integration_weights_for_one_tetrahedron(qp%tet(i),tete,x(l)+deltae,w2,sigma,blochlcorrections=.false.)
!                 buf3(:,l,k)=buf3(:,l,k)+tetpr*(sum( (w1+w2)*tv3 ))*0.5_r8
!                 buf4(:,l,k)=buf4(:,l,k)+tetpr*(sum( (w1+w2)*tv4 ))*0.5_r8
!             enddo
!         enddo
!         if ( mw%talk ) then
!             if ( lo_trueNtimes(li,60,qp%mpi%ntet) ) call lo_progressbar(' ... spectral free energy',li,qp%mpi%ntet)
!         endif
!     enddo
!     ! Sum it ip
!     lo_allocate(buf3b(uc%na,npts,dr%n_mode))
!     lo_allocate(buf4b(uc%na,npts,dr%n_mode))
!     buf3b=0.0_r8
!     buf4b=0.0_r8
!     call mpi_allreduce(buf3,buf3b,uc%na*dr%n_mode*npts,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!     call mpi_allreduce(buf4,buf4b,uc%na*dr%n_mode*npts,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!
!     if ( mw%talk ) then
!         call lo_progressbar(' ... spectral free energy',li,qp%mpi%ntet,walltime()-t0)
!     endif
!
!     ! Scale to reasonable units
!     x=x*unitfactor
!     buf3b=buf3b/unitfactor
!     buf4b=buf4b/unitfactor
!     ! Normalize things properly. First calculate the total contribution.
!     ftot3=0.0_r8
!     ftot4=0.0_r8
!     l=0
!     do i=1,qp%n_irr_point
!     do j=1,dr%n_mode
!         l=l+1
!         ftot3=ftot3+qp%ip(i)%weight*F3(j,i)
!         ftot4=ftot4+F4(j,i)*qp%ip(i)%weight/uc%na
!     enddo
!     enddo
!     ! Sum up the spectral to a total.
!     lo_allocate(yt3(npts))
!     lo_allocate(yt4(npts))
!     yt3=0.0_r8
!     yt4=0.0_r8
!     do j=1,dr%n_mode
!     do i=1,uc%na
!         yt3=yt3+buf3b(i,:,j)
!         yt4=yt4+buf4b(i,:,j)
!     enddo
!     enddo
!     nf3=ftot3/lo_trapezoid_integration(x,yt3)
!     nf4=ftot4/lo_trapezoid_integration(x,yt4)
!     yt3=yt3*nf3
!     yt4=yt4*nf4
!     ! Now construct the site-projected
!     lo_allocate(yp3(npts,uc%na))
!     lo_allocate(yp4(npts,uc%na))
!     yp3=0.0_r8
!     yp4=0.0_r8
!     do j=1,dr%n_mode
!     do i=1,uc%na
!         yp3(:,i)=yp3(:,i)+buf3b(i,:,j)
!         yp4(:,i)=yp4(:,i)+buf4b(i,:,j)
!     enddo
!     enddo
!     ! Make sure the site-projected add up to the total
!     do i=1,npts
!         f0=abs(sum(yp3(i,:)))
!         if ( f0 .gt. 1E-30_r8 ) then
!             yp3(i,:)=yp3(i,:)*abs(yt3(i))/f0
!         else
!             yp3(i,:)=0.0_r8
!         endif
!         f0=abs(sum(yp4(i,:)))
!         if ( f0 .gt. 1E-30_r8 ) then
!             yp4(i,:)=yp4(i,:)*yt4(i)/f0
!         else
!             yp4(i,:)=0.0_r8
!         endif
!     enddo
!
!     ! Dump to file
!     if ( mw%talk ) then
!         ! output format
!         u=open_file('out','dumspf')
!             opf="("//tochar(uc%na*2+3)//"(2X,E18.12))"
!             do i=1,npts
!                 write(u,opf) x(i),yt3(i),yt4(i),yp3(i,:),yp4(i,:)
!             enddo
!         close(u)
!     endif
!
! end subroutine

! !> Calculate all scattering rates, MPI-distributed
! subroutine getscatteringrates(sr,qpoint,ompoint,qp,dr,uc,fc,fct,fcf,verbosity,mw)
!     !> scattering rates
!     class(lo_listofscatteringrates), intent(out) :: sr
!     !> q-point in question
!     type(lo_qpoint) :: qpoint
!     !> harmonic properties
!     type(lo_phonon_dispersions_qpoint), intent(in) :: ompoint
!     !> qpoint mesh
!     class(lo_qpoint_mesh), intent(in) :: qp
!     !> phonon dispersions
!     type(lo_phonon_dispersions), intent(in) :: dr
!     !> crystal structure
!     type(lo_crystalstructure), intent(in) :: uc
!     !> second order forceconstant
!     type(lo_forceconstant_secondorder), intent(inout) :: fc
!     !> third order force constant
!     type(lo_forceconstant_thirdorder), intent(in) :: fct
!     !> fourth order force constants
!     type(lo_forceconstant_fourthorder), intent(in) :: fcf
!     !> how much to talk
!     integer, intent(in) :: verbosity
!     !> MPI communicator
!     type(lo_mpi_helper), intent(in) :: mw
!
!     ! Just keep a list of q3, temporarily
!     type(lo_qpoint) :: qpgamma
!     type(lo_phonon_dispersions_qpoint) :: ompgamma
!     real(r8), dimension(:,:), allocatable :: qvec3
!     complex(r8), dimension(:,:,:), allocatable :: egv3
!
!     ! First get the third frequency at all q-points
!     thirdfreq: block
!         type(lo_qpoint) :: dq3
!         type(lo_phonon_dispersions_qpoint) :: ompoint
!         complex(r8), dimension(:,:,:), allocatable :: egv3buf
!         complex(r8), dimension(:,:,:), allocatable :: Dq
!         complex(r8), dimension(:,:), allocatable :: D
!         real(r8), dimension(:,:), allocatable :: om3buf,qv3buf
!         real(r8), dimension(:,:,:), allocatable :: vel3buf
!         real(r8), dimension(3) :: qv1,qv2,qv3
!         real(r8) :: t0
!         integer :: lq,q,b1,b2,l
!
!         ! Gamma-point things
!         qpgamma%r=0.0_r8
!         qpgamma%r=0.0_r8
!         call lo_get_small_group_of_qpoint(qpgamma,uc)
!         lo_allocate(ompgamma%omega(dr%n_mode))
!         lo_allocate(ompgamma%vel(3,dr%n_mode))
!         lo_allocate(ompgamma%egv(dr%n_mode,dr%n_mode))
!         lo_allocate(ompgamma%degeneracy(dr%n_mode))
!         lo_allocate(ompgamma%degenmode(dr%n_mode,dr%n_mode))
! write(*,*) 'fixme generate thingy'
!         !call ompgamma%generate(fc,uc,qpgamma)
!
!         lo_allocate(sr%omega_gamma(dr%n_mode))
!         sr%omega_gamma=ompgamma%omega
!
!         t0=walltime()
!         lo_allocate(om3buf(dr%n_mode,qp%n_full_point))
!         lo_allocate(vel3buf(3,dr%n_mode,qp%n_full_point))
!         lo_allocate(egv3buf(dr%n_mode,dr%n_mode,qp%n_full_point))
!         lo_allocate(qv3buf(3,qp%n_full_point))
!         qv3buf=0.0_r8
!         om3buf=0.0_r8
!         vel3buf=0.0_r8
!         egv3buf=0.0_r8
!         lo_allocate(Dq(dr%n_mode,dr%n_mode,3))
!         lo_allocate(D(dr%n_mode,dr%n_mode))
!         ! Calculate all the frequencies and stuff for the third q-point
! write(*,*) 'FIXME MPI'
!         ! do lq=1,qp%mpi%nqtot
!         !     q=qp%mpi%qtot(lq)
!         !
!         !     ! Get the q-vectors
!         !     qv1=qpoint%r
!         !     qv2=qp%ap(q)%r
!         !     qv3=-qv1-qv2
!         !     ! Frequencies and stuff at the third point
!         !     dq3%r=qv3
!         !     dq3%r=dq3%r-uc%bz%gshift(dq3%r)
!         !     ! store third q-point for later
!         !     qv3buf(:,q)=dq3%r
!         !     ! get phonon frequencies
!         !     call lo_get_small_group_of_qpoint(dq3,uc)
!         !     call fc%dynamicalmatrix(uc,dq3,D,Dq)
!         !     call fc%frequencies_eigenvectors_groupvelocities(uc,D,Dq,om3buf(:,q),egv3buf(:,:,q),vel3buf(:,:,q),dq3)
!         ! enddo
!         ! Communicate this to all ranks
!         lo_allocate(sr%omega3(dr%n_mode,qp%n_full_point))
!         lo_allocate(sr%vel3(3,dr%n_mode,qp%n_full_point))
!         lo_allocate(egv3(dr%n_mode,dr%n_mode,qp%n_full_point))
!         lo_allocate(qvec3(3,qp%n_full_point))
!         sr%omega3=0.0_r8
!         sr%vel3=0.0_r8
!         egv3=0.0_r8
!         qvec3=0.0_r8
!         call mpi_allreduce(om3buf,sr%omega3,dr%n_mode*qp%n_full_point,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!         call mpi_allreduce(vel3buf,sr%vel3,3*dr%n_mode*qp%n_full_point,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!         call mpi_allreduce(egv3buf,egv3,dr%n_mode*dr%n_mode*qp%n_full_point,MPI_DOUBLE_COMPLEX,MPI_SUM,mw%comm,mw%error)
!         call mpi_allreduce(qv3buf,qvec3,3*qp%n_full_point,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!         if ( verbosity .gt. 0 ) then
!             write(*,*) ''
!             write(*,*) "... got q'' frequencies in ",tochar(walltime()-t0),'s'
!         endif
!         ! cleanup
!         lo_deallocate(om3buf)
!         lo_deallocate(vel3buf)
!         lo_deallocate(egv3buf)
!         lo_deallocate(qv3buf)
!         lo_deallocate(Dq)
!         lo_deallocate(D)
!     end block thirdfreq
!
!     !> get the threephonon matrix elements
!     threephsc: block
!         complex(r8), dimension(:,:), allocatable :: egv
!         complex(r8) :: c0
!         real(r8), dimension(:,:,:,:), allocatable :: psibuf,dblbuf
!         real(r8), dimension(3) :: qv1,qv2,qv3,omega
!         real(r8) :: t0,f0
!         integer :: q,b1,b2,b3,l,k,ii,jj,kk,bb1,bb2,bb3,n1,n2,n3
!         logical, dimension(:,:,:), allocatable :: fixdeg
!
!         lo_allocate(psibuf(dr%n_mode,dr%n_mode,dr%n_mode,qp%n_full_point))
!         lo_allocate(dblbuf(dr%n_mode,dr%n_mode,dr%n_mode,qp%n_full_point))
!         lo_allocate(egv(dr%n_mode,3))
!         psibuf=0.0_r8
!         dblbuf=0.0_r8
!         egv=0.0_r8
!
!         t0=walltime()
!         if ( verbosity .gt. 0 ) call lo_progressbar_init()
!         l=0
!         k=0
!
!         do q=1,qp%n_full_point
!             qv1=-qpoint%r*lo_twopi
!             qv2=-qp%ap(q)%r*lo_twopi
!             qv3=-qvec3(:,q)*lo_twopi
!             do b1=1,dr%n_mode
!             do b2=1,dr%n_mode
!             do b3=1,dr%n_mode
!                 ! MPI division
!                 l=l+1
!                 if ( mod(l,mw%n) .ne. mw%r ) cycle
!                 omega(1)=ompoint%omega(b1)
!                 omega(2)=dr%aq(q)%omega(b2)
!                 omega(3)=sr%omega3(b3,q)
!                 ! prefetch eigenvectors
!                 egv(:,1)=ompoint%egv(:,b1)
!                 egv(:,2)=dr%aq(q)%egv(:,b2)
!                 egv(:,3)=egv3(:,b3,q)
!                 ! To finally get the scattering rates
!                 if ( minval(omega) .gt. dr%omega_min*0.5_r8 ) then
!                     c0=fct%scatteringamplitude(omega,egv,qv2,qv3)
!                     psibuf(b1,b2,b3,q)=abs(conjg(c0)*c0)
!                 else
!                     psibuf(b1,b2,b3,q)=0.0_r8
!                 endif
!                 ! And perhaps the double scattering rate thingy
!                 if ( ompgamma%omega(b3) .gt. dr%omega_min*0.5_r8 ) then
!                     ! b3 needs to be an optical mode
!                     omega(3)=ompgamma%omega(b3)
!                     egv(:,3)=ompgamma%egv(:,b3)
!                     ! And the weird double scatteringrates
!                     dblbuf(b1,b2,b3,q)=0.0_r8 !real(double_scatteringamplitude(fct,omega,egv,qv1,qv2))
!                 else
!                     dblbuf(b1,b2,b3,q)=0.0_r8
!                 endif
!             enddo
!             enddo
!             if ( verbosity .gt. 0 ) then
!                 k=k+1
!                 if ( lo_trueNtimes(k,127,qp%n_full_point*dr%n_mode) ) call lo_progressbar(' ... threephonon matrixelements',k,qp%n_full_point*dr%n_mode)
!             endif
!             enddo
!         enddo
!         lo_allocate(sr%psi_3ph_sq(dr%n_mode,dr%n_mode,dr%n_mode,qp%n_full_point))
!         lo_allocate(sr%psi_3ph_db(dr%n_mode,dr%n_mode,dr%n_mode,qp%n_full_point))
!         sr%psi_3ph_sq=0.0_r8
!         sr%psi_3ph_db=0.0_r8
!         ! sum it up
!         call mpi_allreduce(psibuf,sr%psi_3ph_sq,dr%n_mode*dr%n_mode*dr%n_mode*qp%n_full_point,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!         call mpi_allreduce(dblbuf,sr%psi_3ph_db,dr%n_mode*dr%n_mode*dr%n_mode*qp%n_full_point,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!         if ( verbosity .gt. 0 ) call lo_progressbar(' ... threephonon matrixelements',qp%n_full_point*dr%n_mode,qp%n_full_point*dr%n_mode,walltime()-t0)
!     end block threephsc
!
!     ! And the four-phonon scattering elements
!     fourphsc: block
!         complex(r8), dimension(:,:), allocatable :: dumegv1,dumegv2,dumegv3,dumegv4,egv
!         complex(r8) :: c0
!         real(r8), dimension(:,:,:), allocatable :: psibuf
!         real(r8), dimension(:), allocatable :: dumom1,dumom2,dumom3,dumom4
!         real(r8), dimension(4) :: omega
!         real(r8), dimension(3) :: qv1,qv2,qv3,qv4
!         real(r8) :: omegathres,t0
!         integer :: b1,b2,b3,b4,q1,ctr
!
!         t0=walltime()
!         lo_allocate(psibuf(dr%n_mode,dr%n_mode,qp%n_full_point))
!         psibuf=0.0_r8
!         ! For the other q-points
!         lo_allocate(dumom1(dr%n_mode))
!         lo_allocate(dumom2(dr%n_mode))
!         lo_allocate(dumom3(dr%n_mode))
!         lo_allocate(dumom4(dr%n_mode))
!         lo_allocate(dumegv1(dr%n_mode,dr%n_mode))
!         lo_allocate(dumegv2(dr%n_mode,dr%n_mode))
!         lo_allocate(dumegv3(dr%n_mode,dr%n_mode))
!         lo_allocate(dumegv4(dr%n_mode,dr%n_mode))
!         lo_allocate(egv(dr%n_mode,4))
!         dumom1=0.0_r8
!         dumom2=0.0_r8
!         dumom3=0.0_r8
!         dumom4=0.0_r8
!         dumegv1=0.0_r8
!         dumegv2=0.0_r8
!         dumegv3=0.0_r8
!         dumegv4=0.0_r8
!         omegathres=dr%omega_min*0.5_r8
!         ctr=0
!         if ( verbosity .gt. 0 ) call lo_progressbar_init()
!         do q1=1,qp%n_full_point
!             dumom1=ompoint%omega
!             dumom2=ompoint%omega
!             dumom3=dr%aq(q1)%omega
!             dumom4=dr%aq(q1)%omega
!             dumegv1=ompoint%egv
!             dumegv2=conjg(ompoint%egv)
!             dumegv3=dr%aq(q1)%egv
!             dumegv4=conjg(dr%aq(q1)%egv)
!             qv1=qpoint%r*lo_twopi
!             qv2=-qpoint%r*lo_twopi
!             qv3=qp%ap(q1)%r*lo_twopi
!             qv4=-qp%ap(q1)%r*lo_twopi
!             do b1=1,dr%n_mode
!                 b2=b1
!                 ! mpit thingy
!                 ctr=ctr+1
!                 if ( mod(ctr,mw%n) .ne. mw%r ) cycle
!                 do b3=1,dr%n_mode
!                     b4=b3
!                     omega=[dumom1(b1),dumom2(b2),dumom3(b3),dumom4(b4)]
!                     egv(:,1)=dumegv1(:,b1)
!                     egv(:,2)=dumegv2(:,b2)
!                     egv(:,3)=dumegv3(:,b3)
!                     egv(:,4)=dumegv4(:,b4)
!                     if ( minval(omega) .gt. omegathres ) then
!                         c0=fcf%scatteringamplitude(omega,egv,qv2,qv3,qv4)
!                     else
!                         c0=0.0_r8
!                     endif
!                     psibuf(b1,b3,q1)=real(c0)
!                 enddo
!                 !
!                 if ( verbosity .gt. 0 ) then
!                     if ( lo_trueNtimes(ctr,127,qp%n_full_point*dr%n_mode) ) call lo_progressbar(' ... fourphonon matrixelements',ctr,qp%n_full_point*dr%n_mode)
!                 endif
!             enddo
!         enddo
!         lo_allocate(sr%psi_4ph(dr%n_mode,dr%n_mode,qp%n_full_point))
!         sr%psi_4ph=0.0_r8
!         call mpi_allreduce(psibuf,sr%psi_4ph,dr%n_mode*dr%n_mode*qp%n_full_point,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!         if ( verbosity .gt. 0 ) call lo_progressbar(' ... fourphonon matrixelements',qp%n_full_point*dr%n_mode,qp%n_full_point*dr%n_mode,walltime()-t0)
!     end block fourphsc
!
! end subroutine

! !> harmonic properties at specific q-vector
! subroutine harmonic_things_at_single_q(q,uc,fc,loto,ompoint)
!     !> Cartesian coordinates of the q-point in question
!     real(r8), dimension(3), intent(in) :: q
!     !> crystal structure
!     type(lo_crystalstructure), intent(in) :: uc
!     !> second order force constant
!     type(lo_forceconstant_secondorder), intent(in) :: fc
!     !> electrostatic corrections
!     type(lo_loto), intent(in) :: loto
!     !> harmonic properties at this q-point
!     type(lo_phonon_dispersions_qpoint), intent(out) :: ompoint
!
!     type(lo_qpoint) :: qpoint
!     complex(r8), dimension(:,:,:), allocatable :: Dq
!     complex(r8), dimension(:,:), allocatable :: D
!     integer :: nb
!     !
!     nb=uc%na*3
!     qpoint%v=q
!     qpoint%w=qpoint%v-uc%bz%gshift(q+lo_degenvector)
!     call lo_get_small_group_of_qpoint(qpoint,uc)
!     lo_allocate(ompoint%omega(nb))
!     lo_allocate(ompoint%vel(3,nb))
!     lo_allocate(ompoint%egv(nb,nb))
!     lo_allocate(D(nb,nb))
!     lo_allocate(Dq(3,nb,nb))
!     call lo_get_dynamical_matrix(fc,uc,qpoint,loto,D,Dq)
!     call lo_get_omega_and_velocities(D,Dq,uc,ompoint%omega,ompoint%egv,ompoint%vel,qpoint=qpoint)
!     lo_deallocate(D)
!     lo_deallocate(Dq)
! end subroutine

! !> The three-phonon matrix element, but product of two complex thingies
! pure complex(r8) function double_scatteringamplitude(fc,omega,egv,q1,q2)
!     !> the third order force constant
!     class(lo_forceconstant_thirdorder), intent(in) :: fc
!     !> the frequencies in question
!     real(r8), dimension(3), intent(in) :: omega
!     !> the eigenvectors
!     complex(r8), dimension(:,:), intent(in) :: egv
!     !> the two q-vectors that matter
!     real(r8), dimension(3), intent(in) :: q1,q2
!     !
!     integer :: atom1,atom2,atom3,i,j,k,k1,k2,k3,t
!     complex(r8) :: expiqr,c0,c1
!     real(r8), dimension(3) :: rv2,rv3
!     real(r8) :: iqr,omegaprod1,omegaprod2,enhet
!     !
!     complex(r8), dimension(:,:,:,:,:,:), allocatable :: egvprod1,egvprod2
!     !
!     ! Not enough room on the stack for large unitcells, have to allocate space
!     allocate(egvprod1(3,3,3,fc%na,fc%na,fc%na))
!     allocate(egvprod2(3,3,3,fc%na,fc%na,fc%na))
!     ! Eigenvector product
!     do atom3=1,fc%na
!     do atom2=1,fc%na
!     do atom1=1,fc%na
!         do k=1,3
!         do j=1,3
!         do i=1,3
!             k1=(atom1-1)*3+i
!             k2=(atom2-1)*3+j
!             k3=(atom3-1)*3+k
!             egvprod1(i,j,k,atom1,atom2,atom3)=egv(k1,1)*conjg(egv(k2,1))*egv(k3,3)
!             egvprod2(i,j,k,atom1,atom2,atom3)=egv(k1,2)*conjg(egv(k2,2))*conjg(egv(k3,3))
!         enddo
!         enddo
!         enddo
!     enddo
!     enddo
!     enddo
!     ! Frequency product
!     omegaprod1=omega(1)*omega(1)*omega(3)
!     omegaprod2=omega(2)*omega(2)*omega(3)
!     omegaprod1=sqrt(omegaprod1)
!     omegaprod2=sqrt(omegaprod2)
!     ! The actual matrix element
!     c0=0.0_r8
!     c1=0.0_r8
!     do atom1=1,fc%na
!     do t=1,fc%atom(atom1)%n
!         atom2=fc%atom(atom1)%triplet(t)%i2
!         atom3=fc%atom(atom1)%triplet(t)%i3
!         rv2=fc%atom(atom1)%triplet(t)%lv2
!         rv3=fc%atom(atom1)%triplet(t)%lv3
!         !
!         iqr=dot_product(-q1,rv2) !+dot_product(q3,rv3)
!         expiqr=cmplx(cos(iqr),sin(iqr),r8)
!         c0=c0+sum(fc%atom(atom1)%triplet(t)%mwm*egvprod1(:,:,:,atom1,atom2,atom3))*expiqr
!         !
!         iqr=dot_product(-q2,rv2) !+dot_product(-q3,rv3)
!         expiqr=cmplx(cos(iqr),sin(iqr),r8)
!         c1=c1+sum(fc%atom(atom1)%triplet(t)%mwm*egvprod2(:,:,:,atom1,atom2,atom3))*expiqr
!     enddo
!     enddo
!     ! I predid the multiplication
!     enhet=2.367755374064161E51_r8
!     c0=enhet*c0/omegaprod1
!     c1=enhet*c1/omegaprod2
!     !
!     double_scatteringamplitude=c0*c1
! end function
