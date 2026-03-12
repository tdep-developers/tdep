
!> remove secondorder forces
subroutine subtract_secondorder_forces(gs,map,mw,mem,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    type(lo_forceconstant_secondorder) :: fcss
    real(flyt), dimension(:,:), allocatable :: force_rsq
    real(flyt), dimension(:,:), allocatable :: f
    real(flyt), dimension(3) :: v0
    real(flyt) :: energy,tt0
    integer :: gpoint,ii,t,a1,a2,i

    tt0=walltime()
    lo_allocate(force_rsq(2,gs%nsim))
    lo_allocate(f(3,map%n_atom_ss))
    f=0.0_flyt
    force_rsq=0.0_flyt

    if ( verbosity .gt. 0 ) call lo_progressbar_init()
    do gpoint=1,gs%raw%nrelevant_gridpoints
        ii=gs%raw%relevant_gridpoints(gpoint)
        ! get the forceconstant
        forceconst: block
            type(lo_crystalstructure) :: uc,ss
            type(lo_forceconstant_secondorder) :: fc
            real(flyt), dimension(:,:), allocatable :: constraints
            integer :: nconstr

            call uc%generate( gs%ref(ii)%unitcell_latticevectors,gs%ref(ii)%unitcell_positions,gs%ref(ii)%unitcell_atomic_numbers,enhet=2 )
            call ss%generate( gs%ref(ii)%supercell_latticevectors,gs%ref(ii)%supercell_positions,gs%ref(ii)%supercell_atomic_numbers,enhet=2 )
            call ss%classify( 'supercell',uc )
            !call lo_secondorder_rot_herm_huang( map,uc,constraints,nconstr,.true.,.true.,.true. )
            if ( nconstr .gt. 0 ) then
                call gs%eval( map,gs%grid_coordinates(:,ii),constraints )
            else
                call gs%eval( map,gs%grid_coordinates(:,ii) )
            endif
            call map%get_secondorder_forceconstant(uc,fc,mem,-1)
            call fc%remap(uc,ss,fcss)
        end block forceconst
        ! and forces/energies
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            f=0.0_flyt
            energy=0.0
            do a1=1,fcss%na
                do i=1,fcss%atom(a1)%n
                    a2=fcss%atom(a1)%pair(i)%i2
                    v0=matmul(fcss%atom(a1)%pair(i)%m,gs%raw%u(:,a2,t))
                    f(:,a1)=f(:,a1)-v0
                enddo
                energy=energy-dot_product(gs%raw%u(:,a1,t),f(:,a1))*0.5_flyt
            enddo
            ! sanity check that they add up to zero
            if ( abs(sum(f)) .gt. lo_tol ) then
                call lo_stop_gracefully(['Second order forces do not add up to zero.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
            ! Calculate R^2 values
            do a1=1,map%n_atom_ss
            do i=1,3
                ! I say that the average force is defined as zero?
                force_rsq(1,ii)=force_rsq(1,ii)+( f(i,a1)-gs%raw%f(i,a1,t) )**2
                force_rsq(2,ii)=force_rsq(2,ii)+( gs%raw%f0(i,a1,t) )**2
            enddo
            enddo
            ! subtract
            gs%raw%f(:,:,t)=gs%raw%f(:,:,t)-f
            ! store energy
            gs%raw%e_fc_pair(t)=energy
        enddo
        if ( verbosity .gt. 0 ) call lo_progressbar(' ... subtracting secondorder forces',gpoint,gs%raw%nrelevant_gridpoints,walltime()-tt0)
    enddo
    ! Add the R^2 things up over MPI
    call mpi_allreduce(MPI_IN_PLACE,force_rsq,gs%nsim*2,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    ! Store the R^2 somewhere
    lo_allocate(gs%pair%rsquare(gs%nsim))
    do ii=1,gs%nsim
        gs%pair%rsquare(ii)=1.0_flyt-force_rsq(1,ii)/force_rsq(2,ii)
    enddo
end subroutine

!> get secondorder coefficients, but with gridfit instead.
subroutine solve_secondorder_gridfit(gs,map,mw,mem,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    character(len=5000), dimension(:), allocatable :: filenames

    ! set up things
    init: block
        real(flyt), dimension(gs%ndim) :: d1
        character(len=1) :: dum
        integer :: i,u

        allocate(filenames(gs%nsim))
        ! grab simulations again. Don't have the energy to sort again.
        u=open_file('in','infile.simulations')
            read(u,*) dum
            read(u,*) dum
            read(u,*) dum
            read(u,*) dum
            read(u,*) dum
            if ( allocated(gs%eos) ) then
                read(u,*) dum
            endif
            do i=1,gs%nsim
                read(u,*) d1,filenames(i)
            enddo
        close(u)
    end block init

    getcoeff: block
        type(lo_crystalstructure) :: uc,ss
        type(lo_forceconstant_secondorder) :: fc
        type(lo_mdsim) :: sim
        real(flyt), dimension(:,:,:,:), allocatable :: dipole_forceconstant
        real(flyt), dimension(:,:), allocatable :: CM,wC,theta_pair,f
        real(flyt), dimension(:,:), allocatable :: wA,wB
        real(flyt), dimension(:), allocatable :: vB,vD,vsol
        real(flyt) :: timer_coeff,t0
        integer, dimension(gs%ndim) :: maxord
        integer :: ii,i,j,t,l,nu,nf,nc,nt,a1,a2

        timer_coeff=walltime()
        t0=timer_coeff

        !nu=map%nfc_pair
        nu=map%xuc%nx_fc_pair
        nf=map%n_atom_ss*3
        lo_allocate(CM(nf,nu))
        lo_allocate(vsol(nu))
        lo_allocate(theta_pair(gs%nsim,nu))
        CM=0.0_flyt
        vsol=0.0_flyt
        theta_pair=0.0_flyt
        if ( map%polar .gt. 0 ) then
            lo_allocate(dipole_forceconstant(3,3,map%n_atom_ss,map%n_atom_ss))
            lo_allocate(f(3,map%n_atom_ss))
        endif

        do ii=1,gs%nsim
            if ( mod(ii,mw%n) .ne. mw%r ) cycle
            ! read a simulation
            call sim%read_from_hdf5(trim(filenames(ii)),verbosity=0)
            ! get structure and constraints
            call uc%generate( sim%extra%unitcell_latticevectors, sim%extra%unitcell_positions, sim%extra%unitcell_atomic_numbers, 2 )
            ! fake forceconstant thingy
            !call lo_secondorder_rot_herm_huang(map,uc,wC,nc,.true.,.true.,.true.)

            if ( map%polar .gt. 0 ) then
                ! subtract the polar stuff here instead. Can't figure out how to not do it twice.

                call ss%generate(gs%ref(ii)%supercell_latticevectors,gs%ref(ii)%supercell_positions, &
                                 gs%ref(ii)%supercell_atomic_numbers, enhet=2 )
                call ss%classify('supercell',uc)
                do i=1,gs%polar%nZ
                    map%xuc%x_Z_singlet(i)=gs%polar%ipZ%eval(i,gs%grid_coordinates(:,ii) )
                    !map%theta_Z(i)=gs%polar%ipZ%eval(i,gs%grid_coordinates(:,ii) )
                enddo
                do i=1,gs%polar%neps
                    map%xuc%x_eps_global(i)=gs%polar%ipeps%eval(i,gs%grid_coordinates(:,ii) )
                    !map%theta_eps(i)=gs%polar%ipeps%eval(i,gs%grid_coordinates(:,ii) )
                enddo
                call map%get_secondorder_forceconstant(uc,fc,mem,-1)
                call fc%supercell_longrange_dynamical_matrix_at_gamma(ss,dipole_forceconstant,1E-12_flyt)
                do t=1,sim%nt
                    f=0.0_flyt
                    do a1=1,sim%na
                    do a2=1,sim%na
                        f(:,a1)=f(:,a1)-matmul(dipole_forceconstant(:,:,a1,a2),sim%u(:,a2,t))
                    enddo
                    enddo
                    ! sanity check that they add up to zero
                    if ( abs(sum(f)) .gt. lo_sqtol ) then
                        call lo_stop_gracefully(['dipole-dipole forces do not add up to zero'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                    endif
                    sim%f(:,:,t)=sim%f(:,:,t)-f
                enddo
                ! Also, while I have this, subtract it from the massive grid
            endif

            ! first get the coefficient matrix
            nt=sim%nt*nf
            lo_allocate(wA(nt,nu))
            lo_allocate(vB(nt))
            CM=0.0_flyt
            wA=0.0_flyt
            vB=0.0_flyt
            l=0
            do t=1,sim%nt
                call lo_coeffmatrix_pair(sim%u(:,:,t),CM,map)
                wA( (t-1)*nf+1:t*nf,: )=CM
                do i=1,sim%na
                do j=1,3
                    l=l+1
                    vB(l)=-sim%f(j,i,t)
                enddo
                enddo
            enddo

            ! solve everything?
            if ( nc .gt. 0 ) then
                lo_allocate(vD(nc))
                vD=0.0_flyt
                call lo_dgglse(wA,wC,vB,vD,vsol,info=lo_status)
                if ( lo_status .ne. 0 ) then
                    call lo_stop_gracefully(['dgglse exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                endif
                lo_deallocate(vD)
            else
                lo_allocate(wB(nt,1))
                wB(:,1)=vB
                call lo_dgels(wA,wB,info=lo_status)
                if ( lo_status .ne. 0 ) then
                    call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                else
                    vsol=wB(1:nu,1)
                endif
                lo_deallocate(wB)
            endif

            ! store solution
            theta_pair(ii,:)=vsol

            lo_deallocate(wA)
            lo_deallocate(vB)
            ! Report?
            if ( verbosity .gt. 0 ) then
                 if ( walltime()-t0 .gt. 2.0_flyt ) then
                    call lo_looptimer('... pair coefficients',timer_coeff,walltime(),ii,gs%nsim)
                    t0=walltime()
                endif
            endif
        enddo
        lo_deallocate(CM)
        lo_deallocate(vsol)
        ! Add up all the pair solutions
        call mpi_allreduce(MPI_IN_PLACE,theta_pair,gs%nsim*nu,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)

        ! And initialize the interpolation
        maxord=2
        call gs%gridpair%ip%generate( gs%grid_coordinates, theta_pair, 2, gs%info%distance_scale, maxord )

        ! store some extra things
        gs%gridpair%nfc=nu
    end block getcoeff
end subroutine

!> solve to get equations
subroutine solve_secondorder(gs,map,mw,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> how much to talk
    integer, intent(in) :: verbosity

    real(flyt), dimension(:,:), allocatable :: CTC,CTF,phi_constraints
    real(flyt) :: tt0
    integer :: nfc,nphi,ncoeff,neq,neqtot,nf,nconstr

    tt0=walltime()
    ! figure out some stuff.
    init: block
        ! Some shorthand indices that might be useful
        nfc=map%xuc%nx_fc_pair
        nphi=nfc*gs%poly%ncoeff
        ncoeff=gs%poly%ncoeff
        nf=map%n_atom_ss*3
        neq=gs%raw%nconf*nf
        neqtot=0
        call mpi_allreduce(neq,neqtot,1,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
        ! I need to know the number of equations on each rank, I think
        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'SOLVING FOR PAIR FORCECONSTANTS'
        endif
    end block init

    ! Fix all the constraints: different at every gridpoint
    constr: block
        type(lo_crystalstructure) :: uc
        real(flyt), dimension(:,:,:), allocatable :: buf !rb,sb
        real(flyt), dimension(:,:), allocatable :: A
        integer, dimension(:), allocatable :: cts,ctr
        real(flyt) :: t0
        integer :: ii,jj,kk

        t0=walltime()
        lo_allocate(cts(gs%nsim))
        lo_allocate(ctr(gs%nsim))
        cts=0
        ctr=0
        ! Different constraints for every gridpoint. Can I rely on the idea that all ranks always will have the same
        ! number of constraints? Probably not.
        do ii=1,gs%nsim
            if ( mod(ii,mw%n) .ne. mw%r ) cycle
            ! Generate structures for this point
            call uc%generate( gs%ref(ii)%unitcell_latticevectors, gs%ref(ii)%unitcell_positions, gs%ref(ii)%unitcell_atomic_numbers ,enhet=1 )
            ! Grab the constraints at this point
            
            ! call lo_secondorder_rot_herm_huang(map,uc,A,jj,.true.,.true.,.true.)
            jj=0
            ! Store these away, temporarily
            if ( jj .gt. 0 ) then
                gs%ref(ii)%nconstr_pair=jj
                lo_allocate(gs%ref(ii)%constr_pair(size(A,1),size(A,2)))
                gs%ref(ii)%constr_pair=A
            else
                ! no constraints
                gs%ref(ii)%nconstr_pair=0
            endif
            ! And store the counter
            cts(ii)=jj
        enddo

        call mpi_allreduce(cts,ctr,gs%nsim,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
        ! largest number of constraints
        jj=maxval(ctr)
        if ( jj .gt. 0 ) then
            lo_allocate(buf(jj,nfc,gs%nsim))
            buf=0.0_flyt
            ! grab the constraints I stored away
            do ii=1,gs%nsim
                if ( mod(ii,mw%n) .ne. mw%r ) cycle
                kk=gs%ref(ii)%nconstr_pair
                if ( kk .gt. 0 ) then
                    buf(1:kk,:,ii)=gs%ref(ii)%constr_pair
                    lo_deallocate(gs%ref(ii)%constr_pair)
                endif
            enddo
            ! sum them up
            call mpi_allreduce(MPI_IN_PLACE,buf,gs%nsim*nfc*jj,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
            ! store in appropriate place
            do ii=1,gs%nsim
                kk=ctr(ii)
                if ( kk .gt. 0 ) then
                    gs%ref(ii)%nconstr_pair=kk
                    lo_allocate(gs%ref(ii)%constr_pair(kk,nfc))
                    gs%ref(ii)%constr_pair=buf(1:kk,:,ii)
                else
                    gs%ref(ii)%nconstr_pair=0
                endif
            enddo
            lo_deallocate(buf)
        else
            ! There are no constraints, at all!
            do ii=1,gs%nsim
                gs%ref(ii)%nconstr_pair=0
            enddo
        endif
        ! and the total number of constraints, across all gridpoints
        nconstr=sum(ctr)
        if ( verbosity .gt. 0 ) write(*,*) '... built constraints (',tochar(walltime()-t0),')'
    end block constr

    ! Get the coefficient matrices
    paircoeff: block
        type(lo_sparsematrix) :: sAM
        real(flyt), dimension(:,:), allocatable :: partA,partB,partC,dfA
        real(flyt) :: t0,timer_coeff
        integer, dimension(:), allocatable :: cstartind
        integer :: i,j,k,l,ii,jj,t
#ifdef AGRESSIVE_SANITY
        real(flyt) :: f0,mem,toMiB
        if ( mw%talk ) then
            write(*,*) ''
            write(*,*) 'projected memory use:'
            toMiB=1.0_flyt/8.0_flyt/1024.0_flyt**2
            mem=0.0_flyt

            f0=real(nfc,flyt)*real(ncoeff,flyt)*(storage_size(f0)+storage_size(i)*2)*toMiB
            mem=mem+f0
            write(*,"(1X,A,2X,F12.5)") ' sparse augmentM:',f0
            f0=real(nconstr,flyt)*real(nphi,flyt)*real(gs%nsim,flyt)*storage_size(f0)*toMiB
            mem=mem+f0
            write(*,"(1X,A,2X,F12.5)") '             dfA:',f0
            f0=real(nphi,flyt)*real(nphi,flyt)*storage_size(f0)*toMiB
            mem=mem+f0
            write(*,"(1X,A,2X,F12.5)") '             CTC:',f0
            f0=real(nphi,flyt)*storage_size(f0)*toMiB
            mem=mem+f0
            write(*,"(1X,A,2X,F12.5)") '             CTF:',f0
            f0=real(nf,flyt)*real(nfc,flyt)*storage_size(f0)*toMiB
            mem=mem+f0
            write(*,"(1X,A,2X,F12.5)") '           partA:',f0
            f0=real(nf,flyt)*real(nphi,flyt)*storage_size(f0)*toMiB
            mem=mem+f0
            write(*,"(1X,A,2X,F12.5)") '           partB:',f0
            f0=real(nf,flyt)*storage_size(f0)*toMiB
            mem=mem+f0
            write(*,"(1X,A,2X,F12.5)") '           partC:',f0
            write(*,"(1X,A,2X,F12.5)") '           TOTAL:',mem
        endif
#endif

        t0=walltime()

        ! Space for the sparse representation of the augmentation matrix
        sAM%nrow=nfc
        sAM%ncol=nphi
        sAM%n=nfc*ncoeff
        lo_allocate(sAM%rowind(sAM%n))
        lo_allocate(sAM%colind(sAM%n))
        lo_allocate(sAM%val(sAM%n))
        sAM%rowind=0
        sAM%colind=0
        sAM%val=0.0_flyt

        ! space for the constraints
        if ( nconstr .gt. 0 ) then
            lo_allocate(dfA(nconstr,nphi))
            dfA=0.0_flyt
            lo_allocate(cstartind(gs%nsim))
            cstartind=0
            j=0
            do i=1,gs%nsim-1
                j=j+gs%ref(i)%nconstr_pair
                cstartind(i+1)=j
            enddo
        endif

        ! Space for the thingy to solve
        lo_allocate(CTC(nphi,nphi))
        lo_allocate(CTF(nphi,1))
        CTC=0.0_flyt
        CTF=0.0_flyt
        ! Space for the (local) coefficient matrix
        lo_allocate(partA(nf,nfc))
        lo_allocate(partB(nf,nphi))
        lo_allocate(partC(nf,1))
        partA=0.0_flyt
        partB=0.0_flyt
        partC=0.0_flyt

        timer_coeff=walltime()
        do t=1,gs%raw%nconf
            ! which gridpoint are we on?
            ii=gs%raw%gridind(t)
            ! normal coefficient matrix
            call lo_coeffmatrix_pair(gs%raw%u(:,:,t),partA,map)

            ! Augment the coefficient matrix, first construct augmentation matrix
            sAM%val=lo_huge
            l=0
            do i=1,nfc
            do j=1,ncoeff
                k=(i-1)*ncoeff+j
                l=l+1
                sAM%rowind(l)=i
                sAM%colind(l)=k
                sAM%val(l)=gs%poly%coeffM(ii,j)
            enddo
            enddo
            ! Matrix multiplicataion, but sparse, and manual. Probably fast enough anyway
            partB=0.0_flyt
            do l=1,sAM%n
                k=sAM%rowind(l)
                j=sAM%colind(l)
                do i=1,nf
                    partB(i,j)=partB(i,j)+partA(i,k)*sAM%val(l)
                enddo
            enddo

            ! fetch forces
            k=0
            do i=1,map%n_atom_ss
            do j=1,3
                k=k+1
                partC(k,1)=-gs%raw%f(j,i,t)
            enddo
            enddo

            ! Add weights
            partB=partB*gs%ref(ii)%weight
            partC=partC*gs%ref(ii)%weight
            ! Multadd this together!
            call lo_gemm(partB,partB,CTC,transa='T',transb='N',alpha=1.0_flyt,beta=1.0_flyt)
            call lo_gemm(partB,partC,CTF,transa='T',transb='N',alpha=1.0_flyt,beta=1.0_flyt)
            ! And maybe constraints?
            if ( nconstr .gt. 0 ) then
                jj=cstartind(ii)
                ! Manual matrix multiplication again
                do l=1,sAM%n
                    k=sAM%rowind(l)
                    j=sAM%colind(l)
                    do i=1,gs%ref(ii)%nconstr_pair
                        dfA(jj+i,j)=dfA(jj+i,j)+gs%ref(ii)%constr_pair(i,k)*sAM%val(l)
                    enddo
                enddo
            endif
            ! Report?
            if ( verbosity .gt. 0 ) then
                 if ( walltime()-t0 .gt. timereport ) then
                    call lo_looptimer('... pair coefficients',timer_coeff,walltime(),t,gs%raw%nconf)
                    t0=walltime()
                endif
            endif
        enddo

        ! Skip SCALAPACK altogether. I am so smart! I am so smart! SMRT! Just build stuff!
        t0=walltime()
        call mpi_allreduce(MPI_IN_PLACE,CTC,nphi*nphi,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
        call mpi_allreduce(MPI_IN_PLACE,CTF,nphi,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
        if ( nconstr .gt. 0 ) then
            ! and compress the constraints
            ii=size(dfA,1)*size(dfA,2)
            call mpi_allreduce(MPI_IN_PLACE,dfA,ii,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
            call reduce_equations(dfA,phi_constraints,i)
            lo_deallocate(dfA)
        endif
        lo_deallocate(partA)
        lo_deallocate(partB)
        lo_deallocate(partC)
        if ( verbosity .gt. 0 ) write(*,*) '... communicated matrices (',tochar(walltime()-t0),')'
    end block paircoeff

    ! Store solution in a reasonable form
    lsq: block
        real(flyt), dimension(:,:), allocatable :: wmA,wmB
        real(flyt), dimension(:), allocatable :: solution
        real(flyt) :: t0
        integer :: i,j,k,ii
        ! solve the new system instead. Don't have to bother with anything parallel here, it's small enough. I can also constrain it!!!
        ! Choose solver. Wisely.
        t0=walltime()
        if ( nconstr .gt. 0 ) then
            ! Normal solver
            lo_allocate(solution(nphi))
            ! number of constraining equations
            ii=size(phi_constraints,2)
            ! Build big matrix
            lo_allocate(wmA(nphi+ii,nphi+ii))
            lo_allocate(wmB(nphi+ii,1))
            wmA=0.0_flyt
            wmB=0.0_flyt
            ! build the big matrix
            wmA( 1:nphi,1:nphi )=CTC
            do i=1,nphi
            do j=1,ii
                wmA(nphi+j,i)=phi_constraints(i,j)
                wmA(i,nphi+j)=phi_constraints(i,j)
            enddo
            enddo
            ! and the force thing
            wmB(1:nphi,1)=CTF(:,1)
            ! and solve
            call lo_dgels(wmA,wmB,info=lo_status)
            if ( lo_status .ne. 0 ) then
                call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
            endif
            solution=wmB(1:nphi,1)
        else
            ! Normal solver
            lo_allocate(solution(nphi))
            call lo_dgels(CTC,CTF,info=lo_status)
            if ( lo_status .ne. 0 ) then
                call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
            endif
            solution=CTF(:,1)
        endif
        if ( verbosity .gt. 0 ) write(*,*) '... solved for secondorder coefficients (',tochar(walltime()-t0),')'
        ! Store the solution for pair
        gs%pair%nvar=nfc
        gs%pair%nconstr=0
        lo_allocate(gs%pair%coeff( gs%poly%ncoeff, gs%pair%nvar ))
        gs%pair%coeff=0.0_flyt
        k=0
        do j=1,gs%pair%nvar
        do i=1,gs%poly%ncoeff
            k=k+1
            gs%pair%coeff(i,j)=solution(k)
        enddo
        enddo
    end block lsq
end subroutine

! !> solve to get equations, using scalack. Useful for reference.
! subroutine solve_secondorder(gs,uc,ss,map,mw,verbosity)
!     !> all the simulations
!     class(lo_gridsim), intent(inout) :: gs
!     !> unitcell
!     type(lo_crystalstructure), intent(in) :: uc
!     !> supercell
!     type(lo_crystalstructure), intent(in) :: ss
!     !> forcemap
!     type(lo_forcemap), intent(inout) :: map
!     !> MPI helper
!     type(lo_mpi_helper), intent(inout) :: mw
!     !> how much to talk
!     integer, intent(in) :: verbosity
!
!     type(lo_blacs_helper) :: bw
!     type(lo_scalapack_matrix) :: mA,mB,mC
!     real(flyt), dimension(:), allocatable :: solution
!     real(flyt) :: tt0
!     integer :: Bnrow,Bncol,nfc,nphi,ncoeff,neq,neqtot,nf
!     integer, dimension(:), allocatable :: eqctr,eqoffset
!
!     tt0=walltime()
!     ! figure out some stuff.
!     init: block
!         integer :: i,j
!         ! How big of a matrix are we going to get?
!         ! If it's constrained, how are we going to deal with that? Later problem.
!
!         ! Some shorthand indices that might be useful
!         nfc=map%nfc_pair
!         nphi=gs%info%nphi_pair
!         ncoeff=gs%poly%ncoeff
!         nf=map%n_atom_ss*3
!         neq=gs%raw%nconf*nf
!         neqtot=0
!         call mpi_allreduce(neq,neqtot,1,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
!         ! I need to know the number of equations on each rank, I think
!         lo_allocate(eqctr(mw%n))
!         lo_allocate(eqoffset(mw%n))
!         eqctr=0
!         eqoffset=0
!         eqoffset(mw%r+1)=neq
!         call mpi_allreduce(eqoffset,eqctr,mw%n,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
!         eqoffset=0
!         j=0
!         do i=1,mw%n-1
!             j=j+eqctr(i)
!             eqoffset(i+1)=j
!         enddo
!         ! Start scalapack communication
!         call bw%init()
!         if ( verbosity .gt. 0 ) then
!             write(*,*) ''
!             write(*,*) 'SOLVING FOR PAIR FORCECONSTANTS'
!         endif
!     end block init
!
!     paircoeff: block
!         real(flyt), dimension(:,:,:), allocatable :: augmentM
!         real(flyt), dimension(:,:), allocatable :: CM,partA,partB,FM
!         real(flyt), dimension(:,:), allocatable :: CTC,CTF,CTCbuf,CTFbuf
!         real(flyt) :: t0
!         integer, dimension(:), allocatable :: rowind
!         integer :: i,j,k,l,ii,t
!
!         t0=walltime()
!         if ( verbosity .gt. 0 ) call lo_progressbar_init()
!
!         ! Get the augmentation matrices
!         lo_allocate(augmentM(nfc,nphi,gs%nsim))
!         augmentM=0.0_flyt
!         do ii=1,gs%nsim
!             do i=1,nfc
!             do j=1,ncoeff
!                 k=(i-1)*ncoeff+j
!                 augmentM(i,k,ii)=gs%poly%coeffM(ii,j)
!             enddo
!             enddo
!         enddo
!
!         ! Space for the (local) coefficient matrix
!         lo_allocate(CM(neq,nphi))
!         lo_allocate(partA(nf,nfc))
!         lo_allocate(partB(nf,nphi))
!         lo_allocate(FM(neq,1))
!         lo_allocate(rowind(neq))
!         CM=0.0_flyt
!         partA=0.0_flyt
!         partB=0.0_flyt
!         FM=0.0_flyt
!         l=0
!         do t=1,gs%raw%nconf
!             ! which gridpoint are we on?
!             ii=gs%raw%gridind(t)
!             ! normal coefficient matrix
!             call coeffmatrix_pair(gs%raw%u(:,:,t),partA,map)
!             ! augment it
!             call lo_gemm(partA,augmentM(:,:,ii),partB)
!             ! store it
!             CM( (t-1)*nf+1:t*nf,: )=partB
!             ! And the forces, and the global row index
!             do i=1,map%n_atom_ss
!             do j=1,3
!                 l=l+1
!                 FM(l,1)=-gs%raw%f(j,i,t)
!                 rowind(l)=l+eqoffset(mw%r+1)
!             enddo
!             enddo
!             if ( verbosity .gt. 0 ) call lo_progressbar(' ... building pair coefficient matrix',t,gs%raw%nconf,walltime()-t0)
!         enddo
!
! ! Skip SCALAPACK altogether. I am so smart! I am so smart! SMRT!
! t0=walltime()
! lo_allocate(CTCbuf(nphi,nphi))
! lo_allocate(CTC(nphi,nphi))
! lo_allocate(CTFbuf(nphi,1))
! lo_allocate(CTF(nphi,1))
! CTC=0.0_flyt
! CTCbuf=0.0_flyt
! CTF=0.0_flyt
! CTFbuf=0.0_flyt
! if ( verbosity .gt. 0 ) write(*,*) '... multiplying really big matrices'
! call lo_gemm(CM,CM,CTCbuf,transa='T')
! call lo_gemm(CM,FM,CTFbuf,transa='T')
! ! add them together
! call mpi_allreduce(CTCbuf,CTC,nphi*nphi,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
! call mpi_allreduce(CTFbuf,CTF,nphi,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
!
! ! solve the new system instead. Don't have to bother with anything parallel here, it's small enough. I can also constrain it!!!
! lo_allocate(solution(nphi))
! call lo_dgels(CTC,CTF)
! solution=CTF(:,1)
!
! if ( verbosity .gt. 0 ) write(*,*) '... multiplied really big matrices (',tochar(walltime()-t0),')'
!
!
! !        ! Init the scalapack matrices
! !        t0=walltime()
! !        call mA%init( neqtot,nphi,bw )
! !        call mB%init( neqtot,1   ,bw )
! !        call mC%init(   nphi,nphi,bw )
! !
! !        ! Gather some statistics to give a short report on the progress
! !
! !        ! Distribute the matrices the way ScaLAPACK likes it.
! !        call mA%reshuffle_matrix( CM,rowind,mw,bw,verbosity )
! !        call mB%reshuffle_matrix( FM,rowind,mw,bw,verbosity )
! !        ! Remove unnecessary space.
! !        lo_deallocate(CM)
! !        lo_deallocate(partA)
! !        lo_deallocate(partB)
! !        lo_deallocate(FM)
! !        lo_deallocate(rowind)
! !        if ( verbosity .gt. 0 ) write(*,*) '... rearranged matrices (',tochar(walltime()-t0),'s)'
! !        t0=walltime()
! !        ! Actually solve the system?
! !        lo_allocate(solution(nphi))
! !        call lo_pdgels(mA,mB,solution,bw,mw)
! !        if ( verbosity .gt. 0 ) write(*,*) '... solved linear system (',tochar(walltime()-t0),'s)'
!     end block paircoeff
!
!     ! Store solution in a reasonable form
!     sortsolution: block
!         integer :: i,j,k,l
!         ! Store the solution for pair
!         gs%pair%nvar=nfc
!         gs%pair%nconstr=0 !map%neq_pair
!         lo_allocate(gs%pair%coeff( gs%poly%ncoeff, gs%pair%nvar ))
!         gs%pair%coeff=0.0_flyt
!         k=0
!         do j=1,gs%pair%nvar
!         do i=1,gs%poly%ncoeff
!             k=k+1
!             gs%pair%coeff(i,j)=solution(k)
!         enddo
!         enddo
!     end block sortsolution
! end subroutine
