
!> read a lot of simulations from file and distribute over ranks
subroutine read_simulations_from_file(gs,filename,mw,mem,verbosity)
    !> grid of simulations
    class(lo_gridsim), intent(inout) :: gs
    !> filename
    character(len=*), intent(in) :: filename
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    ! temporary list of simulations
    type(lo_mdsim), dimension(:), allocatable :: sim
    ! timer
    real(r8) :: timer,t0,t1
    ! which rank will read all the simulations?
    integer :: readrank

    timer=walltime()
    t0=timer
    t1=timer
    ! which rank will read everything?
    readrank=mw%n-1

    ! read raw data from file
    readsims: block
        integer, parameter :: maxfilenamelength=5000
        real(r8), dimension(:), allocatable :: d1
        character(len=maxfilenamelength), dimension(:), allocatable :: sfn
        character(len=1) :: dum
        integer :: u,i,ndim,nfiles

        ! grab the list of files
        if ( mw%r .eq. readrank ) then
            u=open_file('in',trim(filename))
                read(u,*) ndim
                read(u,*) nfiles
                read(u,*) dum
                read(u,*) dum
                read(u,*) dum
                if ( allocated(gs%eos) ) then
                    read(u,*) dum
                endif
                allocate(d1(ndim))
                allocate(sfn(nfiles))
                do i=1,nfiles
                    read(u,*) d1,sfn(i)
                enddo
            close(u)
            ! keep track of the total number of simulations
            gs%nsim=nfiles
            allocate(sim(gs%nsim))
        endif

        ! Spread the number of simulations
        call mw%bcast(gs%nsim,from=readrank)

        ! Read the simulations on one rank
        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        if ( mw%r .eq. readrank ) then
            do i=1,nfiles
                call sim(i)%read_from_hdf5(trim(sfn(i)),verbosity=0)
                if ( i .lt. nfiles ) call lo_progressbar(' ... reading simulations',i,nfiles,walltime()-t0)
            enddo
            deallocate(d1)
            deallocate(sfn)
            t1=walltime()
            call lo_progressbar(' ... reading simulations',nfiles,nfiles,t1-t0)
            t0=t1
        endif
    end block readsims

    !@TODO: add grouping of timesteps here.

    ! Distribute the simulations all over the place.
    distributesims: block
        character, dimension(:), allocatable :: sbuf,rbuf
        integer, dimension(:), allocatable :: steps_per_rank,offset
        integer :: size_per_step,natom
        integer :: isim,it,ctr,irnk,pos
        logical :: have_mag,have_diel

        call mem%allocate(steps_per_rank,mw%n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(offset,mw%n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        steps_per_rank=0
        offset=0
        ! Get the total number of steps
        if ( mw%r .eq. readrank ) then
            ! Steps per rank?
            ctr=0
            do isim=1,gs%nsim
            do it=1,sim(isim)%nt
                ctr=ctr+1
                irnk=mod(ctr,mw%n)
                steps_per_rank(irnk+1)=steps_per_rank(irnk+1)+1
            enddo
            enddo
            ! Offset
            offset=0
            ctr=0
            do irnk=1,mw%n
                offset(irnk)=ctr
                ctr=ctr+steps_per_rank(irnk)
            enddo

            ! Fields that exist?
            have_mag =sim(1)%have_magnetic_moments
            have_diel=sim(1)%have_dielectric

            ! size of one step?
            size_per_step=0
            ! size of displacements/forces
            size_per_step=size_per_step+storage_size(sim(1)%u(:,:,1))*size(sim(1)%u(:,:,1))/8
            size_per_step=size_per_step+storage_size(sim(1)%f(:,:,1))*size(sim(1)%f(:,:,1))/8
            if ( have_mag ) then
                size_per_step=size_per_step+storage_size(sim(1)%m(:,:,1))*size(sim(1)%m(:,:,1))/8
            endif
            if ( have_diel ) then
                size_per_step=size_per_step+storage_size(sim(1)%eps(:,:,1))*size(sim(1)%eps(:,:,1))/8
                size_per_step=size_per_step+storage_size(sim(1)%Z(:,:,:,1))*size(sim(1)%Z(:,:,:,1))/8
            endif
            ! size of energy
            size_per_step=size_per_step+storage_size(sim(1)%stat%potential_energy(1))/8
            ! size of index the simulation came from
            size_per_step=size_per_step+storage_size(ctr)/8

            ! Make a little space for the send buffer
            allocate(sbuf(sum(steps_per_rank)*size_per_step))
            ! Pack it up
            pos=0
            do isim=1,gs%nsim
            do it=1,sim(isim)%nt
                call mw%pack(sim(isim)%u(:,:,it)    ,sbuf,pos)
                call mw%pack(sim(isim)%f(:,:,it)    ,sbuf,pos)
                if ( have_mag ) then
                    call mw%pack(sim(isim)%m(:,:,it)    ,sbuf,pos)
                endif
                if ( have_diel ) then
                    call mw%pack(sim(isim)%eps(:,:,it)  ,sbuf,pos)
                    call mw%pack(sim(isim)%Z(:,:,:,it)    ,sbuf,pos)
                endif
                call mw%pack(sim(isim)%stat%potential_energy(it)  ,sbuf,pos)
                call mw%pack(isim                                 ,sbuf,pos)
            enddo
            enddo

            ! Also make note of the number of atoms
            natom=sim(1)%na

            t1=walltime()
            write(*,*) '... packed simulations (',tochar(t1-t0),'s)'
            t0=t1
        else
            ! Not the readrank, just allocate a buffer to be on the safe side
            allocate(sbuf(1))
        endif
        ! Send some things around
        call mw%bcast(natom,from=readrank)
        call mw%bcast(have_mag,from=readrank)
        call mw%bcast(have_diel,from=readrank)
        call mw%bcast(steps_per_rank,from=readrank)
        call mw%bcast(offset,from=readrank)
        call mw%bcast(size_per_step,from=readrank)
        ! Space in the recieve buffer
        if ( steps_per_rank(mw%r+1) .gt. 0 ) then
            allocate(rbuf(steps_per_rank(mw%r+1)*size_per_step))
        else
            allocate(rbuf(1))
        endif
        ! Scatter all the timesteps
        call mw%scatterv(sbuf,rbuf,steps_per_rank*size_per_step,offset*size_per_step,readrank,__FILE__,__LINE__)

        ! Now I can start storing the information recieved. First make some space
        if ( steps_per_rank(mw%r+1) .gt. 0 ) then
            gs%raw%nconf=steps_per_rank(mw%r+1)
            allocate(gs%raw%u(3,natom,gs%raw%nconf))
            allocate(gs%raw%f(3,natom,gs%raw%nconf))
            allocate(gs%raw%f0(3,natom,gs%raw%nconf))
            allocate(gs%raw%e(gs%raw%nconf))
            allocate(gs%raw%e_polar(gs%raw%nconf))
            allocate(gs%raw%e_fc_pair(gs%raw%nconf))
            allocate(gs%raw%e_fc_triplet(gs%raw%nconf))
            allocate(gs%raw%e_fc_quartet(gs%raw%nconf))
            allocate(gs%raw%gridind(gs%raw%nconf))
            gs%raw%u=0.0_r8
            gs%raw%f=0.0_r8
            gs%raw%f0=0.0_r8
            gs%raw%e=0.0_r8
            gs%raw%e_polar=0.0_r8
            gs%raw%e_fc_pair=0.0_r8
            gs%raw%e_fc_triplet=0.0_r8
            gs%raw%e_fc_quartet=0.0_r8
            gs%raw%gridind=0.0_r8

            if ( have_mag ) then
                allocate(gs%raw%m(3,natom,gs%raw%nconf))
                gs%raw%m=0.0_r8
            endif
            if ( have_diel ) then
                allocate(gs%raw%Z(3,3,natom,gs%raw%nconf))
                allocate(gs%raw%eps(3,3,gs%raw%nconf))
                gs%raw%Z=0.0_r8
                gs%raw%eps=0.0_r8
            endif

            pos=0
            do it=1,gs%raw%nconf
                call mw%unpack(gs%raw%u(:,:,it)    ,rbuf,pos)
                call mw%unpack(gs%raw%f(:,:,it)    ,rbuf,pos)
                if ( have_mag ) then
                    call mw%unpack(gs%raw%m(:,:,it) ,rbuf,pos)
                endif
                if ( have_diel ) then
                    call mw%unpack(gs%raw%eps(:,:,it) ,rbuf,pos)
                    call mw%unpack(gs%raw%Z(:,:,:,it) ,rbuf,pos)
                endif
                call mw%unpack(gs%raw%e(it),        rbuf,pos)
                call mw%unpack(gs%raw%gridind(it),  rbuf,pos)
            enddo
        else
            ! No steps on this rank.
            gs%raw%nconf=0
            gs%raw%nrelevant_gridpoints=0
        endif

        ! And cleanup
        deallocate(sbuf)
        deallocate(rbuf)
        call mem%deallocate(steps_per_rank,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(offset,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated and unpacked simulations (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block distributesims

    ! Spread out the reference data
    fixref: block
        real(r8), dimension(:), allocatable :: wts
        real(r8) :: f0
        integer :: n_atom_uc,n_atom_ss,isim
        integer :: i,j

        allocate(gs%ref(gs%nsim))
        if ( mw%r .eq. readrank ) then
            n_atom_uc = size(sim(1)%extra%unitcell_atomic_numbers)
            n_atom_ss = size(sim(1)%extra%supercell_atomic_numbers)
        endif
        call mw%bcast(n_atom_uc,from=readrank)
        call mw%bcast(n_atom_ss,from=readrank)

        allocate(wts(gs%nsim))
        wts=0.0_r8
        ! Make a little space and send things around

        do isim=1,gs%nsim
            allocate(gs%ref(isim)%unitcell_positions(3,n_atom_uc))
            allocate(gs%ref(isim)%supercell_positions(3,n_atom_ss))
            allocate(gs%ref(isim)%unitcell_atomic_numbers(n_atom_uc))
            allocate(gs%ref(isim)%supercell_atomic_numbers(n_atom_ss))
            allocate(gs%ref(isim)%born_effective_charges(3,3,n_atom_uc))
            gs%ref(isim)%unitcell_latticevectors=0.0_r8
            gs%ref(isim)%supercell_latticevectors=0.0_r8
            gs%ref(isim)%unitcell_positions=0.0_r8
            gs%ref(isim)%supercell_positions=0.0_r8
            gs%ref(isim)%unitcell_atomic_numbers=0
            gs%ref(isim)%supercell_atomic_numbers=0
            gs%ref(isim)%dielectric_tensor=0.0_r8
            gs%ref(isim)%born_effective_charges=0.0_r8
            gs%ref(isim)%weight=0.0_r8

            if ( mw%r .eq. readrank ) then
                gs%ref(isim)%unitcell_latticevectors  =sim(isim)%extra%unitcell_latticevectors
                gs%ref(isim)%supercell_latticevectors =sim(isim)%extra%supercell_latticevectors
                gs%ref(isim)%unitcell_positions       =sim(isim)%extra%unitcell_positions
                gs%ref(isim)%supercell_positions      =sim(isim)%extra%supercell_positions
                gs%ref(isim)%unitcell_atomic_numbers  =sim(isim)%extra%unitcell_atomic_numbers
                gs%ref(isim)%supercell_atomic_numbers =sim(isim)%extra%supercell_atomic_numbers
                gs%ref(isim)%dielectric_tensor        =sim(isim)%extra%dielectric_tensor
                gs%ref(isim)%born_effective_charges   =sim(isim)%extra%born_effective_charges
            endif

            ! And spread it around
            call mw%bcast(gs%ref(isim)%unitcell_latticevectors ,from=readrank)
            call mw%bcast(gs%ref(isim)%supercell_latticevectors,from=readrank)
            call mw%bcast(gs%ref(isim)%unitcell_positions      ,from=readrank)
            call mw%bcast(gs%ref(isim)%supercell_positions     ,from=readrank)
            call mw%bcast(gs%ref(isim)%unitcell_atomic_numbers ,from=readrank)
            call mw%bcast(gs%ref(isim)%supercell_atomic_numbers,from=readrank)
            call mw%bcast(gs%ref(isim)%dielectric_tensor       ,from=readrank)
            call mw%bcast(gs%ref(isim)%born_effective_charges  ,from=readrank)

            ! Calculate some kind of weight per timestep
            if ( mw%r .eq. readrank ) then
            if ( gs%info%weighted ) then
                f0=0.0_r8
                do i=1,sim(isim)%nt
                do j=1,sim(isim)%na
                    f0=f0+norm2(sim(isim)%f(:,j,i))
                enddo
                enddo
                wts(isim)=f0
            endif
            endif
        enddo

        ! Fix weights
        if ( gs%info%weighted ) then
            if ( mw%r .eq. readrank ) then
                wts=1.0_r8/wts
                wts=gs%nsim*wts/sum(wts)
            endif
            call mw%bcast(wts,from=readrank)
            do isim=1,gs%nsim
                gs%ref(isim)%weight=wts(isim)
                if ( mw%talk ) write(*,*) 'point:',isim,'weight:',gs%ref(isim)%weight
            enddo
        else
            ! No weights
            do isim=1,gs%nsim
                gs%ref(isim)%weight=1.0_r8
            enddo
        endif

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated reference data (',tochar(t1-t0),'s)'
            t0=t1
        endif
    end block fixref

    ! Do the final preparations
    finalize: block
        integer, dimension(:), allocatable :: di
        integer :: i,j

        call mem%allocate(di,gs%nsim,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        di=0

        do j=1,gs%raw%nconf
            i=gs%raw%gridind(j)
            di(i)=1
        enddo
        gs%raw%nrelevant_gridpoints=sum(di)
        if ( gs%raw%nrelevant_gridpoints .gt. 0 ) then
            allocate(gs%raw%relevant_gridpoints(gs%raw%nrelevant_gridpoints))
            gs%raw%relevant_gridpoints=0
            j=0
            do i=1,gs%nsim
                if ( di(i) .eq. 0 ) cycle
                j=j+1
                gs%raw%relevant_gridpoints(j)=i
            enddo
        endif

        ! do i=1,gs%raw%nconf
        !     gs%raw%eps(:,:,i)=( gs%raw%eps(:,:,i)+transpose(gs%raw%eps(:,:,i)) )*0.5_r8
        ! enddo

    end block finalize

    ! And summarize
    if ( verbosity .gt. 0 ) write(*,*) '... read and balanced simulations (',tochar(walltime()-timer),'s)'
end subroutine

! !> divide simluation into groups with constant positions, for magnetic things
! subroutine group_simulation(sim,groupcounter,groupind)
!     !> the simulation
!     type(lo_mdsim), intent(in) :: sim
!     !> how many simulations in each group
!     integer, dimension(:), allocatable, intent(out) :: groupcounter
!     !> which simulations are in each group
!     integer, dimension(:,:), allocatable, intent(out) :: groupind
!     !
!     integer, parameter :: grthres=1 ! how many structure for a group to count
!     integer, dimension(:,:), allocatable :: ddi
!     integer, dimension(:), allocatable :: ctr,di
!     integer :: i,j,ngroup
!     real(r8) :: f0
!
!     lo_allocate(ctr(sim%nt))
!     lo_allocate(di(sim%nt))
!     lo_allocate(ddi(sim%nt,sim%nt))
!     ctr=0
!     di=1
!     ddi=0
!     ! Figure out how many unique positions there are?
!     do i=1,sim%nt
!         if ( di(i) .eq. 0 ) cycle
!         do j=i,sim%nt
!             ! really simple check
!             f0=abs(sim%r(1,1,j)-sim%r(1,1,i))
!             if ( f0 .gt. lo_sqtol ) cycle ! not the same
!             ! slightly better check
!             f0=sum(abs(sim%r(:,:,j)-sim%r(:,:,i)))/sim%na/3
!             if ( f0 .gt. lo_sqtol ) cycle ! not the same
!             ! if I made it here, it's equal
!             di(j)=0
!             ctr(i)=ctr(i)+1
!             ddi(ctr(i),i)=j
!         enddo
!     enddo
!
!     ! Build the actual groups
!     ngroup=0
!     do i=1,sim%nt
!         if ( ctr(i) .ge. grthres ) ngroup=ngroup+1
!     enddo
!     if ( ngroup .gt. 0 ) then
!         lo_allocate(groupcounter(ngroup))
!         lo_allocate(groupind(maxval(ctr),ngroup))
!         groupcounter=0
!         groupind=0
!         j=0
!         do i=1,sim%nt
!             if ( ctr(i) .ge. grthres ) then
!                 j=j+1
!                 groupcounter(j)=ctr(i)
!                 groupind(1:groupcounter(j),j)=ddi(1:ctr(i),i)
!             endif
!         enddo
!     endif
! end subroutine
