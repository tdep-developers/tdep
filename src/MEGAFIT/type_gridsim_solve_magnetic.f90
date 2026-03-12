
! !> Brutal non-linear opimization to set the on-site term across the grid.
! subroutine nonlinear_opt_onsite(gs,map,mw,verbosity)
!     !> all the simulations
!     class(lo_gridsim), intent(inout) :: gs
!     !> forcemap
!     type(lo_forcemap), intent(inout) :: map
!     !> MPI helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> how much to talk
!     integer, intent(in) :: verbosity
!
!     init: block
!     end block init
!
! end subroutine

!> Get the magnetic moment thingies, the onsite terms
subroutine solve_magnetic_onsite(gs,map,mw,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> how much to talk
    integer, intent(in) :: verbosity

    integer :: ns,nq

    write(*,*) 'FIXME MAGNETIC ONSITE'
    stop

    ! setup: block
    !     if ( verbosity .gt. 0 ) then
    !         write(*,*) 'SOLVING MAGNETIC ONSITE TERMS'
    !     endif
    !     ns=map%nmagsingletshells
    !     nq=map%ntheta_magpair_qij
    ! end block setup

    ! Do this the stupid way first, fit at each point then interpolate.
    coeff: block
    !     integer, parameter :: niter=3
    !     type(lo_mdsim) :: sim
    !     real(flyt), dimension(:,:), allocatable :: coeffmatrix_qij,tqc,scaledcoord
    !     real(flyt), dimension(:,:), allocatable :: gqij,gmi,gsi,hmi,hsi,hqij
    !     real(flyt), dimension(:,:), allocatable :: sicoeff,magnetic_moments
    !     real(flyt), dimension(:), allocatable :: flat_moments,flat_adjusted_moments,flat_moment_difference
    !     real(flyt), dimension(2) :: histpar
    !     real(flyt) :: t0,timer_coeff,f0,f1,magtref
    !     integer, dimension(:), allocatable :: momentcounter
    !     integer, dimension(gs%ndim) :: maxord
    !     integer :: i,j,k,l,si,u,iter
    !
    !     real(flyt), dimension(gs%ndim) :: d1
    !     character(len=1) :: dum
    !     character(len=5000), dimension(:), allocatable :: filenames
    !
    !     ! start timers
    !     t0=walltime()
    !     timer_coeff=t0
    !
    !     ! grab filenames
    !     lo_allocate(filenames(gs%nsim))
    !     u=open_file('in','infile.simulations')
    !         read(u,*) dum
    !         read(u,*) dum
    !         read(u,*) dum
    !         read(u,*) dum
    !         read(u,*) dum
    !         read(u,*) dum
    !         do i=1,gs%nsim
    !             read(u,*) d1,filenames(i)
    !         enddo
    !     close(u)
    !
    !     ! Scaled coordinates for fit
    !     lo_allocate(scaledcoord(gs%ndim,gs%nsim))
    !     do i=1,gs%nsim
    !         call gs%coordinate_transformation(gs%grid_coordinates(:,i),scaledcoord(:,i))
    !     enddo
    !
    !     ! Space for the solutions
    !     lo_allocate(gqij(gs%nsim,nq))
    !     lo_allocate(gmi(gs%nsim,ns))
    !     lo_allocate(gsi(gs%nsim,ns))
    !     lo_allocate(hqij(gs%nsim,nq))
    !     lo_allocate(hmi(gs%nsim,ns))
    !     lo_allocate(hsi(gs%nsim,ns))
    !     lo_allocate( momentcounter(ns) )
    !     gqij=0.0_flyt
    !     gmi=0.0_flyt
    !     gsi=0.0_flyt
    !     hqij=0.0_flyt
    !     hmi=0.0_flyt
    !     hsi=0.0_flyt
    !     momentcounter=0
    !
    !     ! Guess a reference temperature for the fluctuations
    !     magtref=300.0_flyt
    !
    ! iterloop: do iter=1,niter
    !     ! Solve for each simulation
    !     hmi=0.0_flyt
    !     hsi=0.0_flyt
    !     hqij=0.0_flyt
    !     do si=1,gs%nsim
    !         if ( mod(si,mw%n) .ne. mw%r ) cycle
    !         ! read a simulation
    !         call sim%read_from_hdf5(trim(filenames(si)),verbosity=0)
    !
    !         ! Grab out the magnetic moments in a convenient format. First count the number of moments
    !         ! per kind of thing.
    !         momentcounter=0
    !         do i=1,sim%nt
    !         do j=1,sim%na
    !             k=map%ss(j)%J%irreducible_shell
    !             momentcounter(k)=momentcounter(k)+1
    !         enddo
    !         enddo
    !         ! make some temporary space
    !         lo_allocate( magnetic_moments(maxval(momentcounter),ns) )
    !         lo_allocate( flat_moments(sim%nt*sim%na) )
    !         lo_allocate( flat_adjusted_moments(sim%nt*sim%na) )
    !         lo_allocate( flat_moment_difference(sim%nt*sim%na) )
    !         lo_allocate( coeffmatrix_qij(sim%nt*sim%na,nq) )
    !         lo_allocate( tqc(sim%na,nq) )
    !
    !         ! Lists of the magnetic moments
    !         flat_moments=0.0_flyt
    !         do i=1,sim%nt
    !         do j=1,sim%na
    !             l=(i-1)*sim%na+j
    !             k=map%ss(j)%J%irreducible_shell
    !             flat_moments(l)=norm2(sim%m(:,j,i))
    !             flat_moment_difference(l)=flat_moments(l)-gsi(si,k)
    !         enddo
    !         enddo
    !         ! Coefficient matrix
    !         coeffmatrix_qij=0.0_flyt
    !         tqc=0.0_flyt
    !         do i=1,sim%nt
    !             call lo_coeffmatrix_magnetic_longitudinal(sim%u(:,:,i),tqc,map)
    !             ! Coefficient matrix for the longitudinal guys
    !             coeffmatrix_qij( (i-1)*map%nss+1:i*map%nss,:)=tqc
    !         enddo
    !         ! Get adjusted moments from previous qij
    !         flat_adjusted_moments=flat_moments-matmul(coeffmatrix_qij,gqij(si,:))
    !
    !         magnetic_moments=0.0_flyt
    !         momentcounter=0
    !         do i=1,sim%nt
    !         do j=1,sim%na
    !             k=map%ss(j)%J%irreducible_shell
    !             l=(i-1)*sim%na+j
    !             momentcounter(k)=momentcounter(k)+1
    !             magnetic_moments( momentcounter(k),k )=flat_adjusted_moments(l)
    !         enddo
    !         enddo
    !
    !         ! Fit the onsite term
    !         do i=1,ns
    !             histpar=fit_onsite_model_to_histogram( magnetic_moments(1:momentcounter(i),i) )
    !             if ( iter .eq. 1 ) then
    !                 hmi(si,i)=histpar(1)
    !                 hsi(si,i)=histpar(2)
    !             else
    !                 f0=tempscaler( gs%grid_coordinates( gs%info%dim_temperature,si ),magtref )
    !                 ! Evaluate the old, and mix!
    !                 hmi(si,i)=histpar(1)*0.9_flyt+0.1_flyt*gmi(si,i)
    !                 hsi(si,i)=histpar(2)*0.9_flyt+0.1_flyt*gsi(si,i)/f0
    !             endif
    !         enddo
    !
    !         ! Get the delta-guys
    !         do i=1,sim%nt
    !         do j=1,sim%na
    !             l=(i-1)*sim%na+j
    !             k=map%ss(j)%J%irreducible_shell
    !             flat_moment_difference(l)=flat_moments(l)-gmi(si,k)
    !         enddo
    !         enddo
    !         ! Solve for qij
    !         call linear_least_squares(coeffmatrix_qij,flat_moment_difference,hqij(si,:))
    !
    !         lo_deallocate( magnetic_moments )
    !         lo_deallocate( flat_moments )
    !         lo_deallocate( flat_adjusted_moments )
    !         lo_deallocate( flat_moment_difference )
    !         lo_deallocate( coeffmatrix_qij )
    !         lo_deallocate( tqc )
    !     enddo
    !     ! Now sync things up
    !     call mpi_allreduce(MPI_IN_PLACE,hsi,size(hsi),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    !     call mpi_allreduce(MPI_IN_PLACE,hmi,size(hmi),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    !     call mpi_allreduce(MPI_IN_PLACE,hqij,size(hqij),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    !
    !     ! Optimize reference temperature for the sigma-term
    !     call optimize_tref_for_sigma(gs,hsi,magtref)
    !     ! Scale sigma with reference temperature before interpolation
    !     do i=1,gs%nsim
    !         f0=magtempscaler( gs%grid_coordinates( gs%info%dim_temperature,i ),magtref )
    !         hsi(i,:)=hsi(i,:)*f0
    !     enddo
    !
    !     ! Create interpolations!
    !     maxord=2
    !     call gs%magpair%ipM0%generate( scaledcoord, hmi, 2, 0.15_flyt, maxord )
    !     call gs%magpair%ipMD%generate( scaledcoord, hsi, 2, 0.1_flyt, maxord )
    !     call gs%magpair%ipqij%generate( scaledcoord, hqij, 2, 0.15_flyt, maxord )
    !
    !     ! Evaluate them for next iteration
    !     do i=1,gs%nsim
    !         f0=magtempscaler( gs%grid_coordinates( gs%info%dim_temperature,i ),magtref )
    !         do j=1,ns
    !             gmi(i,j)=gs%magpair%ipM0%eval( j,scaledcoord(:,i),rexp=1 )
    !             gsi(i,j)=gs%magpair%ipMD%eval( j,scaledcoord(:,i),rexp=1 )
    !         enddo
    !         do j=1,nq
    !             gqij(i,j)=gs%magpair%ipqij%eval( j,scaledcoord(:,i),rexp=1 )
    !         enddo
    !     enddo
    !
    !     if ( mw%talk ) then
    !         write(*,*) '... iterating onsite solution:',iter,sum(gqij),sum(gmi),sum(gsi)
    !     endif
    ! enddo iterloop

!if ( mw%talk ) then
!
!    do i=1,gs%nsim
!!        gmi(i,j)=gs%magpair%ipM0%eval( j,scaledcoord(:,i),1 )
!!        gsi(i,j)=gs%magpair%ipMD%eval( j,scaledcoord(:,i),1 )
!        write(*,*) gs%grid_coordinates(:,i),gmi(i,1),hmi(i,1)
!    enddo
!u=open_file('out','dumsig')
!    do i=1,gs%nsim
!        f0=tempscaler( gs%grid_coordinates( gs%info%dim_temperature,i ),magtref )
!!        gmi(i,j)=gs%magpair%ipM0%eval( j,scaledcoord(:,i),1 )
!!        gsi(i,j)=gs%magpair%ipMD%eval( j,scaledcoord(:,i),1 )
!        write(*,*) gs%grid_coordinates(:,i),gsi(i,1),hsi(i,1)
!        write(u,*) gs%grid_coordinates(:,i),gsi(i,1),hsi(i,1)
!    enddo
!close(u)
!!    do i=1,gs%nsim
!!!        gmi(i,j)=gs%magpair%ipM0%eval( j,scaledcoord(:,i),1 )
!!!        gsi(i,j)=gs%magpair%ipMD%eval( j,scaledcoord(:,i),1 )
!!        write(*,*) gs%grid_coordinates(:,i),gsi(i,1),hsi(i,1)
!!    enddo
!
!endif

    end block coeff
end subroutine

!> solve to get equations
subroutine solve_magnetic_pair_gridfit(gs,map,mw,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> how much to talk
    integer, intent(in) :: verbosity

    real(flyt) :: timer
    integer :: njij,ntij

    write(*,*) 'FIXME MAGNETIC PAIR'
    stop

!     init: block
!         timer=walltime()
!         if ( verbosity .gt. 0 ) then
!             write(*,*) ''
!             write(*,*) 'SOLVING MAGNETIC PAIR INTERACTIONS PER POINT'
!         endif
!         ! Some shorthand
!         njij=map%ntheta_magpair_jij
!         ntij=map%ntheta_magpair_tij
!     end block init
!
!     ! Do this the stupid way first, fit at each point then interpolate.
!     jijcoeff: block
!         type(lo_mdsim) :: sim
!         real(flyt), dimension(:,:), allocatable :: jijcoeff,tijcoeff
!         real(flyt), dimension(:,:), allocatable :: raw_jijcoeff,raw_tijcoeff,coeffM,scaledcoord
!         real(flyt), dimension(:), allocatable :: energy_difference,dx
!         real(flyt) :: t0,timer_coeff
!         integer, dimension(:,:), allocatable :: groupind
!         integer, dimension(:), allocatable :: groupcounter
!         integer, dimension(gs%ndim) :: maxord
!         integer :: ngroup,nconf
!         integer :: si,t,i,j,l,ii,jj,ll,u
!
!         real(flyt), dimension(gs%ndim) :: d1
!         character(len=1) :: dum
!         character(len=5000), dimension(:), allocatable :: filenames
!
!         ! start timers
!         t0=walltime()
!         timer_coeff=t0
!
!         ! grab filenames
!         lo_allocate(filenames(gs%nsim))
!         u=open_file('in','infile.simulations')
!             read(u,*) dum
!             read(u,*) dum
!             read(u,*) dum
!             read(u,*) dum
!             read(u,*) dum
!             read(u,*) dum
!             do i=1,gs%nsim
!                 read(u,*) d1,filenames(i)
!             enddo
!         close(u)
!
!         ! Space for the solutions
!         lo_allocate(jijcoeff(gs%nsim,njij))
!         lo_allocate(tijcoeff(gs%nsim,ntij))
!         lo_allocate(dx(njij+ntij))
!         jijcoeff=0.0_flyt
!         tijcoeff=0.0_flyt
!         dx=0.0_flyt
!
!         ! Solve for each simulation
!         do si=1,gs%nsim
!             if ( mod(si,mw%n) .ne. mw%r ) cycle
!             ! read a simulation
!             call sim%read_from_hdf5(trim(filenames(si)),verbosity=0)
!             ! Group the steps
!             call sim%group_timesteps(groupcounter,groupind,ngroup,nconf)
!             ! Throw an error of there are too few timesteps
!             if ( nconf .lt. max(njij,ntij) ) then
!                 jijcoeff(si,:)=0.0_flyt
!                 tijcoeff(si,:)=0.0_flyt
!                 cycle
!             endif
!             ! First get the coefficients, per config.
!             lo_allocate(raw_jijcoeff(njij,sim%nt))
!             lo_allocate(raw_tijcoeff(ntij,sim%nt))
!             raw_jijcoeff=0.0_flyt
!             raw_tijcoeff=0.0_flyt
!             do t=1,sim%nt
!                 call lo_coeffmatrix_magnetic_pair(sim%m(:,:,t),raw_jijcoeff(:,t),map)
!                 call lo_coeffmatrix_magnetic_crossterm_energy(sim%u(:,:,t),sim%m(:,:,t),raw_tijcoeff(:,t),map)
!             enddo
!             ! Build the full coefficient matrix
!             lo_allocate(coeffM(nconf,njij+ntij))
!             lo_allocate(energy_difference(nconf))
!             coeffM=0.0_flyt
!             energy_difference=0.0_flyt
!             l=0
!             do ll=1,ngroup
!                 do i=1,groupcounter(ll)
!                 do j=i+1,groupcounter(ll)
!                     l=l+1
!                     ii=groupind(i,ll)
!                     jj=groupind(j,ll)
!                     coeffM(l,1:njij)=raw_jijcoeff(:,ii)-raw_jijcoeff(:,jj)
!                     coeffM(l,njij+1:njij+ntij)=raw_tijcoeff(:,ii)-raw_tijcoeff(:,jj)
!                     energy_difference(l)=sim%stat%potential_energy(ii)-sim%stat%potential_energy(jj)
!                 enddo
!                 enddo
!             enddo
!             ! Solve it
!             call linear_least_squares(coeffM,energy_difference,dx)
!             jijcoeff(si,:)=dx(1:njij)
!             tijcoeff(si,:)=dx(njij+1:ntij+njij)
!             ! Cleanup
!             lo_deallocate(raw_jijcoeff)
!             lo_deallocate(raw_tijcoeff)
!             lo_deallocate(coeffM)
!             lo_deallocate(energy_difference)
!             ! Talk?
!             if ( verbosity .gt. 0 ) then
!                  if ( walltime()-t0 .gt. 2.0_flyt ) then
!                     call lo_looptimer('... Jij coefficients',timer_coeff,walltime(),si,gs%nsim)
!                     t0=walltime()
!                 endif
!             endif
!          enddo
!
!          ! Ok, we have all the coefficients. Communicate them and create the interpolations
!          call mpi_allreduce(MPI_IN_PLACE,jijcoeff,size(jijcoeff),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!          call mpi_allreduce(MPI_IN_PLACE,tijcoeff,size(tijcoeff),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!          lo_allocate(scaledcoord(gs%ndim,gs%nsim))
!          do i=1,gs%nsim
!              call gs%coordinate_transformation(gs%grid_coordinates(:,i),scaledcoord(:,i))
!          enddo
!          maxord=2
!          call gs%magpair%ipJ%generate( scaledcoord, jijcoeff, 2, 0.75_flyt, maxord )
!          call gs%magpair%ipT%generate( scaledcoord, tijcoeff, 2, 0.75_flyt, maxord )
!      end block jijcoeff
!
! end subroutine
!
! !> solve to get equations
! subroutine solve_magnetic_pair_polyfit(gs,map,mw,verbosity)
!     !> all the simulations
!     class(lo_gridsim), intent(inout) :: gs
!     !> forcemap
!     type(lo_forcemap), intent(inout) :: map
!     !> MPI helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> how much to talk
!     integer, intent(in) :: verbosity
!
!     real(flyt), dimension(:,:), allocatable :: jijCTC,jijCTF
!     real(flyt), dimension(:,:), allocatable :: tijCTC,tijCTF
!     integer, dimension(:,:,:), allocatable :: groupind
!     integer, dimension(:,:), allocatable :: groupctr,confctr
!     integer, dimension(:), allocatable :: ngroup
!     real(flyt) :: tt0
!     integer :: nfc_jij,nfc_tij,nfc_qij,nfc_si,ncoeff
!     integer :: nphi_jij,nphi_tij,nphi_qij,nphi_si
!
!     tt0=walltime()
!     ! figure out some stuff.
!     init: block
!         ! Some shorthand indices that might be useful
!         nfc_jij=map%ntheta_magpair_jij
!         nfc_tij=map%ntheta_magpair_tij
!         nfc_qij=map%ntheta_magpair_qij
!         nfc_si=map%nmagsingletshells
!         ncoeff=gs%poly%ncoeff
!         nphi_jij=nfc_jij*ncoeff
!         nphi_tij=nfc_tij*ncoeff
!         nphi_qij=nfc_qij*ncoeff
!         nphi_si=nfc_si*ncoeff
!         ! I need to know the number of equations on each rank, I think
!         if ( verbosity .gt. 0 ) then
!             write(*,*) ''
!             write(*,*) 'SOLVING FOR MAGNETIC PAIR INTERACTIONS'
!             write(*,*) '      jij:',nfc_jij,nphi_jij
!             write(*,*) '      tij:',nfc_tij,nphi_tij
!             write(*,*) '      qij:',nfc_qij,nphi_qij
!             write(*,*) '       si:',nfc_si,nphi_si
!         endif
!         ! Chop into groups
!         call group_distributed_calculations(gs,groupind,groupctr,confctr,ngroup)
!     end block init
!
! properjijcoeff: block
!     type(lo_sparsematrix) :: sJM,sTM
!     real(flyt), dimension(:,:), allocatable :: jijpartA,jijpartB,jijpartC
!     real(flyt), dimension(:,:), allocatable :: tijpartA,tijpartB,tijpartC
!
!          real(flyt), dimension(:,:,:), allocatable :: cf2
!          real(flyt), dimension(:,:), allocatable :: cf1,cf3
!          real(flyt), dimension(:), allocatable :: ce1
!          real(flyt) :: t0,timer_coeff
!          integer :: i,j,k,l,ii,jj,t
!          integer :: ia,ib,ic,ig,nf,ng,ctr
!
!          t0=walltime()
!          ! Space for the sparse representation of the augmentation matrix
!          sJM%nrow=nfc_jij
!          sJM%ncol=nphi_jij
!          sJM%n=nfc_jij*ncoeff
!          lo_allocate(sJM%rowind(sJM%n))
!          lo_allocate(sJM%colind(sJM%n))
!          lo_allocate(sJM%val(sJM%n))
!          sJM%rowind=0
!          sJM%colind=0
!          sJM%val=0.0_flyt
!
!          sTM%nrow=nfc_tij
!          sTM%ncol=nphi_tij
!          sTM%n=nfc_tij*ncoeff
!          lo_allocate(sTM%rowind(sTM%n))
!          lo_allocate(sTM%colind(sTM%n))
!          lo_allocate(sTM%val(sTM%n))
!          sTM%rowind=0
!          sTM%colind=0
!          sTM%val=0.0_flyt
!
!          ! Space for the thingy to solve
!          lo_allocate(jijCTC(nphi_jij,nphi_jij))
!          lo_allocate(jijCTF(nphi_jij,1))
!          lo_allocate(tijCTC(nphi_tij,nphi_tij))
!          lo_allocate(tijCTF(nphi_tij,1))
!          jijCTC=0.0_flyt
!          jijCTF=0.0_flyt
!          tijCTC=0.0_flyt
!          tijCTF=0.0_flyt
!          lo_allocate(ce1(maxval(groupctr)))
!          lo_allocate(cf1(nfc_jij,maxval(groupctr)))
!          lo_allocate(cf2(map%nss*3,nfc_tij,maxval(groupctr)))
!          lo_allocate(cf3(map%nss*3,maxval(groupctr)))
!          ce1=0.0_flyt
!          cf1=0.0_flyt
!          cf2=0.0_flyt
!          cf3=0.0_flyt
!
!          ctr=0
!          do ii=1,gs%nsim
!          do ig=1,ngroup(ii)
!              nf=confctr(ig,ii)
!              ng=confctr(ig,ii)*map%nss*3
!              lo_allocate(jijpartA(nf,nfc_jij))
!              lo_allocate(jijpartB(nf,nphi_jij))
!              lo_allocate(jijpartC(nf,1))
!              lo_allocate(tijpartA(ng,nfc_tij))
!              lo_allocate(tijpartB(ng,nphi_tij))
!              lo_allocate(tijpartC(ng,1))
!              jijpartA=0.0_flyt
!              jijpartB=0.0_flyt
!              jijpartC=0.0_flyt
!              tijpartA=0.0_flyt
!              tijpartB=0.0_flyt
!              tijpartC=0.0_flyt
!              ! Get raw coefficient matrices and energies
!              ce1=0.0_flyt
!              cf1=0.0_flyt
!              cf2=0.0_flyt
!              cf3=0.0_flyt
!              do i=1,groupctr(ig,ii)
!                  ia=groupind(i,ig,ii)
!                  call lo_coeffmatrix_magnetic_pair(gs%raw%m(:,:,ia),cf1(:,i),map)
!                  call lo_coeffmatrix_magnetic_crossterm_forces(gs%raw%m(:,:,ia),cf2(:,:,i),map)
!                  ce1(i)=gs%raw%e(ia)
!                  do j=1,map%nss
!                  do k=1,3
!                      l=(j-1)*3+k
!                      cf3(l,i)=gs%raw%f(k,j,ia)
!                  enddo
!                  enddo
!              enddo
!              ! Weight by how many there are in the group
!              ce1=ce1/real(groupctr(ig,ii),flyt)
!              cf1=cf1/real(groupctr(ig,ii),flyt)
!              cf2=cf2/real(groupctr(ig,ii),flyt)
!              cf3=cf3/real(groupctr(ig,ii),flyt)
!
!              ! Then the difference coefficient matrix
!              l=0
!              do i=1,groupctr(ig,ii)
!              do j=i+1,groupctr(ig,ii)
!                  l=l+1
!                  jijpartA(l,:)=cf1(:,i)-cf1(:,j)
!                  jijpartC(l,1)=ce1(i)-ce1(j)
!                  !
!                  tijpartA( (l-1)*map%nss*3+1:l*map%nss*3, : )=cf2(:,:,i)-cf2(:,:,j)
!                  tijpartC( (l-1)*map%nss*3+1:l*map%nss*3, 1 )=-( cf3(:,i)-cf3(:,j) )
!              enddo
!              enddo
!
!              ! Augment the coefficient matrix, first construct augmentation matrix
!              sJM%val=lo_huge
!              l=0
!              do i=1,nfc_jij
!              do j=1,ncoeff
!                  k=(i-1)*ncoeff+j
!                  l=l+1
!                  sJM%rowind(l)=i
!                  sJM%colind(l)=k
!                  sJM%val(l)=gs%poly%coeffM(ii,j)
!              enddo
!              enddo
!              sTM%val=lo_huge
!              l=0
!              do i=1,nfc_tij
!              do j=1,ncoeff
!                  k=(i-1)*ncoeff+j
!                  l=l+1
!                  sTM%rowind(l)=i
!                  sTM%colind(l)=k
!                  sTM%val(l)=gs%poly%coeffM(ii,j)
!              enddo
!              enddo
!
!              ! Matrix multiplicataion, but sparse, and manual. Probably fast enough anyway
!              jijpartB=0.0_flyt
!              do l=1,sJM%n
!                  k=sJM%rowind(l)
!                  j=sJM%colind(l)
!                  do i=1,nf
!                      jijpartB(i,j)=jijpartB(i,j)+jijpartA(i,k)*sJM%val(l)
!                  enddo
!              enddo
!              tijpartB=0.0_flyt
!              do l=1,sTM%n
!                  k=sTM%rowind(l)
!                  j=sTM%colind(l)
!                  do i=1,ng
!                      tijpartB(i,j)=tijpartB(i,j)+tijpartA(i,k)*sTM%val(l)
!                  enddo
!              enddo
!              ! Multadd this together!
!              call lo_gemm(jijpartB,jijpartB,jijCTC,transa='T',transb='N',alpha=1.0_flyt,beta=1.0_flyt)
!              call lo_gemm(jijpartB,jijpartC,jijCTF,transa='T',transb='N',alpha=1.0_flyt,beta=1.0_flyt)
!              call lo_gemm(tijpartB,tijpartB,tijCTC,transa='T',transb='N',alpha=1.0_flyt,beta=1.0_flyt)
!              call lo_gemm(tijpartB,tijpartC,tijCTF,transa='T',transb='N',alpha=1.0_flyt,beta=1.0_flyt)
!              ! Report?
!              ctr=ctr+1
!              if ( verbosity .gt. 0 .and. ctr .lt. sum(ngroup) ) then
!                   if ( walltime()-t0 .gt. timereport ) then
!                      call lo_looptimer('... magnetic coefficients',timer_coeff,walltime(),ctr,sum(ngroup))
!                      t0=walltime()
!                  endif
!              endif
!              ! Cleanup
!              lo_deallocate(jijpartA)
!              lo_deallocate(jijpartB)
!              lo_deallocate(jijpartC)
!              lo_deallocate(tijpartA)
!              lo_deallocate(tijpartB)
!              lo_deallocate(tijpartC)
!          enddo
!          enddo
!
!          ! Skip SCALAPACK altogether. I am so smart! I am so smart! SMRT! Just build stuff!
!          t0=walltime()
!          call mpi_allreduce(MPI_IN_PLACE,jijCTC,size(jijCTC),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!          call mpi_allreduce(MPI_IN_PLACE,jijCTF,size(jijCTF),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!          call mpi_allreduce(MPI_IN_PLACE,tijCTC,size(tijCTC),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!          call mpi_allreduce(MPI_IN_PLACE,tijCTF,size(tijCTF),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!
!          if ( verbosity .gt. 0 ) write(*,*) '... communicated matrices (',tochar(walltime()-t0),')'
!     end block properjijcoeff
!
!      ! Store solution in a reasonable form
!      lsq: block
!         real(flyt), dimension(:), allocatable :: jijsolution,tijsolution
!         real(flyt) :: t0
!         integer :: i,j,k,ii
!         t0=walltime()
!
!         ! Normal solver
!         lo_allocate(jijsolution(nphi_jij))
!         lo_allocate(tijsolution(nphi_tij))
!         call lo_dgels(jijCTC,jijCTF,info=lo_status)
!         if ( lo_status .ne. 0 ) then
!             call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
!         endif
!         call lo_dgels(tijCTC,tijCTF,info=lo_status)
!         if ( lo_status .ne. 0 ) then
!             call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
!         endif
!         jijsolution=jijCTF(:,1)
!         tijsolution=tijCTF(:,1)
!
!         if ( verbosity .gt. 0 ) write(*,*) '... solved for magnetic coefficients (',tochar(walltime()-t0),')'
!         ! Store the solution for pair
!         gs%magpair%nfc=nfc_jij
!         lo_allocate(gs%magpair%coeff( gs%poly%ncoeff, gs%magpair%nfc ))
!         gs%magpair%coeff=0.0_flyt
!         k=0
!         do j=1,gs%magpair%nfc
!         do i=1,gs%poly%ncoeff
!             k=k+1
!             gs%magpair%coeff(i,j)=jijsolution(k)
!         enddo
!         enddo
!
!         ! Store the solution for first crossterm
!         gs%magpair%nfc_T=nfc_tij
!         lo_allocate(gs%magpair%coeff_T( gs%poly%ncoeff, gs%magpair%nfc_T ))
!         gs%magpair%coeff_T=0.0_flyt
!         k=0
!         do j=1,gs%magpair%nfc_T
!         do i=1,gs%poly%ncoeff
!             k=k+1
!             gs%magpair%coeff_T(i,j)=tijsolution(k)
!         enddo
!         enddo
!
!      end block lsq
end subroutine

!> remove magnetic forces and energies
subroutine subtract_magnetic_forces(gs,map,mw,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> how much to talk
    integer, intent(in) :: verbosity

write(*,*) 'FIXME MAGNETIC FORCES'
stop

!     type(lo_jij_secondorder) :: jijss
!     real(flyt), dimension(:,:), allocatable :: force_rsq,energy_rsq_jij,energy_rsq_tij,energy_rsq_qij
!     real(flyt), dimension(:,:), allocatable :: f
!     real(flyt), dimension(3) :: v0,v1,v2
!     real(flyt) :: energy,crossenergy,longenergy,tt0,f0,f1,deltaM
!     integer, dimension(:,:,:), allocatable :: groupind
!     integer, dimension(:,:), allocatable :: groupctr,confctr
!     integer, dimension(:), allocatable :: ngroup
!     integer :: gpoint,ii,t,a1,a2,i,j,ig,ia,ib
!
!     tt0=walltime()
!     lo_allocate(force_rsq(2,gs%nsim))
!     lo_allocate(energy_rsq_jij(2,gs%nsim))
!     lo_allocate(energy_rsq_tij(2,gs%nsim))
!     lo_allocate(energy_rsq_qij(2,gs%nsim))
!     lo_allocate(f(3,map%nss))
!     f=0.0_flyt
!     force_rsq=0.0_flyt
!     energy_rsq_jij=0.0_flyt
!     energy_rsq_tij=0.0_flyt
!     energy_rsq_qij=0.0_flyt
!
!     if ( verbosity .gt. 0 ) call lo_progressbar_init()
!     do gpoint=1,gs%raw%nrelevant_gridpoints
!         ii=gs%raw%relevant_gridpoints(gpoint)
!         ! get the Jij
!         forceconst: block
!             type(lo_crystalstructure) :: uc,ss
!             type(lo_jij_secondorder) :: jij
!
!             call uc%generate( gs%ref(ii)%unitcell_latticevectors,gs%ref(ii)%unitcell_positions,gs%ref(ii)%unitcell_atomic_numbers,enhet=2 )
!             call ss%generate( gs%ref(ii)%supercell_latticevectors,gs%ref(ii)%supercell_positions,gs%ref(ii)%supercell_atomic_numbers,enhet=2 )
!             call ss%classify( 'supercell',uc )
!             call gs%eval( map,gs%grid_coordinates(:,ii) )
!             call map%get_secondorder_jij(uc,jij)
!             call jij%remap(uc,ss,jijss)
!         end block forceconst
!
!         ! and forces/energies
!         do t=1,gs%raw%nconf
!             if ( gs%raw%gridind(t) .ne. ii ) cycle
!             f=0.0_flyt
!             energy=0.0
!             crossenergy=0.0
!             longenergy=0.0
!             do a1=1,jijss%na ! map%nss
!                 v1=0.0_flyt
!                 deltaM=0.0_flyt
!                 do i=1,jijss%atom(a1)%n
!                     a2=jijss%atom(a1)%pair(i)%i2
!                     ! jij energy
!                     v1=v1+matmul(jijss%atom(a1)%pair(i)%J,gs%raw%m(:,a2,t))
!                     ! tij thingy
!                     v0=matmul(jijss%atom(a1)%pair(i)%bird,gs%raw%m(:,a2,t))
!                     f(:,a1)=f(:,a1)-v0
!                     ! qij thingy
!                     deltaM=deltaM+dot_product(jijss%atom(a1)%pair(i)%Q,gs%raw%u(:,a2,t))
!                 enddo
!                 f(:,a1)=f(:,a1)-matmul(jijss%atom(a1)%selfterm_bird,gs%raw%m(:,a1,t))
!                 energy=energy+dot_product(gs%raw%m(:,a1,t),v1)*0.5_flyt
!                 crossenergy=crossenergy+dot_product(gs%raw%u(:,a1,t),f(:,a1))*0.5_flyt
!                 !longenergy=longenergy+jijss%longitudinal_energy( norm2(gs%raw%m(:,a1,t)),jijss%atom(a1)%m0+deltaM,jijss%atom(a1)%m0dev)
!             enddo
!             ! sanity check that they add up to zero
!             if ( abs(sum(f)) .gt. lo_tol ) then
!                 call lo_stop_gracefully(['Magnetic crossterm forces do not add up to zero.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
!             endif
! !            ! Calculate R^2 values
! !            do a1=1,map%nss
! !            do i=1,3
! !                ! I say that the average force is defined as zero?
! !                force_rsq(1,ii)=force_rsq(1,ii)+( f(i,a1)-gs%raw%f(i,a1,t) )**2
! !                force_rsq(2,ii)=force_rsq(2,ii)+( gs%raw%f0(i,a1,t) )**2
! !            enddo
! !            enddo
!             ! subtract
!             !gs%raw%f(:,:,t)=gs%raw%f(:,:,t)-f
!             ! store energy
!             gs%raw%e_magnetic(t)=energy
!             gs%raw%e_magnetic_jij(t)=energy
!             gs%raw%e_magnetic_tij(t)=crossenergy
!             gs%raw%e_magnetic_qij(t)=longenergy
!         enddo
!
!         if ( verbosity .gt. 0 ) call lo_progressbar(' ... subtracting magnetic forces',gpoint,gs%raw%nrelevant_gridpoints,walltime()-tt0)
!     enddo
!
!     ! Calculate the energy R^2, first chop into groupd
!     call group_distributed_calculations(gs,groupind,groupctr,confctr,ngroup)
!
!     ! Now try to lump it together somehow. Hmm.
!     do ii=1,gs%nsim
!     do ig=1,ngroup(ii)
!         do i=1,groupctr(ig,ii)
!         do j=i+1,groupctr(ig,ii)
!             ia=groupind(i,ig,ii)
!             ib=groupind(j,ig,ii)
!             f1=gs%raw%e(ia)-gs%raw%e(ib)
!             f0=gs%raw%e_magnetic_jij(ia)-gs%raw%e_magnetic_jij(ib)
!             energy_rsq_jij(1,ii)=energy_rsq_jij(1,ii)+( f0-f1 )**2
!             energy_rsq_jij(2,ii)=energy_rsq_jij(2,ii)+( f1 )**2
!             f0=f0+gs%raw%e_magnetic_tij(ia)-gs%raw%e_magnetic_tij(ib)
!             energy_rsq_tij(1,ii)=energy_rsq_tij(1,ii)+( f0-f1 )**2
!             energy_rsq_tij(2,ii)=energy_rsq_tij(2,ii)+( f1 )**2
!             f0=f0+gs%raw%e_magnetic_qij(ia)-gs%raw%e_magnetic_qij(ib)
!             energy_rsq_qij(1,ii)=energy_rsq_qij(1,ii)+( f0-f1 )**2
!             energy_rsq_qij(2,ii)=energy_rsq_qij(2,ii)+( f1 )**2
!         enddo
!         enddo
!     enddo
!     enddo
!     call mpi_allreduce(MPI_IN_PLACE,energy_rsq_jij,gs%nsim*2,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!     call mpi_allreduce(MPI_IN_PLACE,energy_rsq_tij,gs%nsim*2,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!     call mpi_allreduce(MPI_IN_PLACE,energy_rsq_qij,gs%nsim*2,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!     ! Store the R^2 somewhere
!     lo_allocate(gs%magpair%rsquare_energy_jij(gs%nsim))
!     lo_allocate(gs%magpair%rsquare_energy_tij(gs%nsim))
!     lo_allocate(gs%magpair%rsquare_energy_qij(gs%nsim))
!     lo_allocate(gs%magpair%rsquare_force(gs%nsim))
!     do ii=1,gs%nsim
!         gs%magpair%rsquare_energy_jij(ii)=1.0_flyt-energy_rsq_jij(1,ii)/energy_rsq_jij(2,ii)
!         gs%magpair%rsquare_energy_tij(ii)=1.0_flyt-energy_rsq_tij(1,ii)/energy_rsq_tij(2,ii)
!         gs%magpair%rsquare_energy_qij(ii)=1.0_flyt-energy_rsq_qij(1,ii)/energy_rsq_qij(2,ii)
!         gs%magpair%rsquare_force(ii)=1.0_flyt-force_rsq(1,ii)/force_rsq(2,ii)
!     enddo
end subroutine
