
!> remove fourthorder forces
subroutine subtract_fourthorder_forces(gs,map,mw,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> how much to talk
    integer, intent(in) :: verbosity

    type(lo_forceconstant_fourthorder) :: fcss
    real(flyt), dimension(:,:), allocatable :: force_rsq
    real(flyt), dimension(:,:), allocatable :: f
    real(flyt), dimension(3,3,3,3) :: m
    real(flyt), dimension(3) :: v0,u2,u3,u4
    real(flyt) :: energy,tt0
    integer :: gpoint,ii,t,a1,a2,a3,a4,i,i1,i2,i3,i4

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
            type(lo_forceconstant_fourthorder) :: fc
            call uc%generate( gs%ref(ii)%unitcell_latticevectors, gs%ref(ii)%unitcell_positions, gs%ref(ii)%unitcell_atomic_numbers ,enhet=1 )
            call ss%generate( gs%ref(ii)%supercell_latticevectors, gs%ref(ii)%supercell_positions, gs%ref(ii)%supercell_atomic_numbers, enhet=1 )
            call ss%classify( 'supercell',uc )
            ! Evaluate the forceconstant
            call gs%eval( map,gs%grid_coordinates(:,ii) )
            call map%get_fourthorder_forceconstant(uc,fc)
            call fc%remap(uc,ss,fcss)
        end block forceconst
        ! and forces/energies
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            f=0.0_flyt
            energy=0.0
            do a1=1,fcss%na
                v0=0.0_flyt
                do i=1,fcss%atom(a1)%n
                    m=fcss%atom(a1)%quartet(i)%m
                    a2=fcss%atom(a1)%quartet(i)%i2
                    a3=fcss%atom(a1)%quartet(i)%i3
                    a4=fcss%atom(a1)%quartet(i)%i4
                    u2=gs%raw%u(:,a2,t)
                    u3=gs%raw%u(:,a3,t)
                    u4=gs%raw%u(:,a4,t)
                    do i1=1,3
                    do i2=1,3
                    do i3=1,3
                    do i4=1,3
                        v0(i1)=v0(i1)-m(i1,i2,i3,i4)*u2(i2)*u3(i3)*u4(i4)
                    enddo
                    enddo
                    enddo
                    enddo
                enddo
                f(:,a1)=v0/6.0_flyt
                energy=energy-dot_product(gs%raw%u(:,a1,t),f(:,a1))/4.0_flyt
            enddo
            ! sanity check that they add up to zero
            if ( abs(sum(f)) .gt. lo_tol ) then
                call lo_stop_gracefully(['Fourth order forces do not add up to zero.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
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
            gs%raw%e_fc_quartet(t)=energy
        enddo

        if ( verbosity .gt. 0 ) call lo_progressbar(' ... subtracting fourthorder forces',gpoint,gs%raw%nrelevant_gridpoints,walltime()-tt0)
    enddo
    ! Add the R^2 things up over MPI
    call mpi_allreduce(MPI_IN_PLACE,force_rsq,gs%nsim*2,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    ! Store the R^2 somewhere
    lo_allocate(gs%quartet%rsquare(gs%nsim))
    do ii=1,gs%nsim
        gs%quartet%rsquare(ii)=1.0_flyt-force_rsq(1,ii)/force_rsq(2,ii)
    enddo
end subroutine

!> solve equations for polynomial coefficients
subroutine solve_fourthorder(gs,map,mw,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> how much to talk
    integer, intent(in) :: verbosity

    real(flyt), dimension(:,:), allocatable :: CTC,CTF,fc_constraints,phi_constraints
    real(flyt) :: tt0
    integer :: nfc,nphi,ncoeff,neq,neqtot,nf,nconstr
    logical, dimension(81) :: relntheta
    logical :: constrain

    tt0=walltime()
    ! figure out some stuff.
    init: block
        real(flyt) :: mem
        integer :: i,j
        ! Some shorthand indices that might be useful
        nfc=map%xuc%nx_fc_quartet
        nphi=nfc*gs%poly%ncoeff
        ncoeff=gs%poly%ncoeff
        nf=map%n_atom_ss*3
        neq=gs%raw%nconf*nf
        neqtot=0
        call mpi_allreduce(neq,neqtot,1,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
        ! I need to know the number of equations on each rank, I think
        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'SOLVING FOR QUARTET FORCECONSTANTS'
        endif

        ! Do the constraint decision, depending on memory usage
        mem=real(map%constraints%neq4,flyt)*real(nphi,flyt)*real(gs%nsim,flyt)*storage_size(mem)/8.0_flyt/1024.0_flyt**2
        ! if small enough, proceed anyway!
        if ( mem .lt. maxmem .and. map%constraints%neq4 .gt. 0 ) then
            constrain=.true.
            nconstr=map%constraints%neq4
            allocate(fc_constraints,source=map%constraints%eq4)
        else
            constrain=.false.
            nconstr=0
        endif

        ! And a helper list with the possible ntheta
        relntheta=.false.
        do i=1,map%n_fc_quartet_shell
            j=map%fc_quartet_shell(i)%nx !quartetshell(i)%ntheta
            if ( j .gt. 0 ) relntheta(j)=.true.
        enddo
    end block init

    ! Get the coefficient matrices
    quartcoeff: block
        type(lo_sparsematrix) :: sAM
        real(flyt), dimension(:,:), allocatable :: partA,partB,partC,dfA
        real(flyt) :: t0,timer_coeff
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
            f0=real(map%constraints%neq3,flyt)*real(nphi,flyt)*real(gs%nsim,flyt)*storage_size(f0)*toMiB
            if ( constrain ) then
                mem=mem+f0
                write(*,"(1X,A,2X,F12.5)") '             dfA:',f0
            else
                write(*,"(1X,A,2X,F12.5,2X,A)") '             dfA:',f0,'(not used)'
            endif
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

        timer_coeff=walltime()

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
        if ( constrain ) then
            lo_allocate(dfA(nconstr*gs%nsim,nphi))
            dfA=0.0_flyt
        endif
        ! space for the buffer
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

        t0=walltime()
        do t=1,gs%raw%nconf
            ! which gridpoint are we on?
            ii=gs%raw%gridind(t)
            ! normal coefficient matrix
            call lo_coeffmatrix_quartet(gs%raw%u(:,:,t),partA,map,relntheta)

            ! Augment the coefficient matrix, first construct augmentation matrix
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
            if ( constrain ) then
                ! Manual matrix multiplication again
                jj=(ii-1)*nconstr
                do l=1,sAM%n
                    k=sAM%rowind(l)
                    j=sAM%colind(l)
                    do i=1,nconstr
                        dfA(jj+i,j)=dfA(jj+i,j)+fc_constraints(i,k)*sAM%val(l)
                    enddo
                enddo
            endif
            ! And report?
            if ( verbosity .gt. 0 ) then
                 if ( walltime()-t0 .gt. timereport ) then
                    call lo_looptimer('... quartet coefficients',timer_coeff,walltime(),t,gs%raw%nconf)
                    t0=walltime()
                endif
            endif
        enddo
        lo_deallocate(partA)
        lo_deallocate(partB)
        lo_deallocate(partC)

        ! Skip SCALAPACK altogether. I am so smart! I am so smart! SMRT! Just build stuff!
        t0=walltime()
        ! add them together
        call mpi_allreduce(MPI_IN_PLACE,CTC,nphi*nphi,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
        call mpi_allreduce(MPI_IN_PLACE,CTF,nphi,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
        if ( nconstr .gt. 0 ) then
            ! and compress the constraints
            call mpi_allreduce(MPI_IN_PLACE,dfA,gs%nsim*nconstr*nphi,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
            call reduce_equations(dfA,phi_constraints,i)
            lo_deallocate(dfA)
        endif
        if ( verbosity .gt. 0 ) write(*,*) '... communicated big matrices (',tochar(walltime()-t0),')'
    end block quartcoeff

    ! Store solution in a reasonable form
    lsq: block
        real(flyt), dimension(:,:), allocatable :: wmA,wmB
        real(flyt), dimension(:), allocatable :: solution
        real(flyt) :: t0
        integer :: i,j,k,ii
        ! solve the new system instead. Don't have to bother with anything parallel here, it's small enough. I can also constrain it!!!
        ! Choose solver. Wisely.
        t0=walltime()
        if ( constrain ) then
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
        if ( verbosity .gt. 0 ) write(*,*) '... solved for fourthorder coefficients (',tochar(walltime()-t0),')'
        ! Store the solution for pair
        gs%quartet%nvar=nfc
        lo_allocate(gs%quartet%coeff( gs%poly%ncoeff, gs%quartet%nvar ))
        gs%quartet%coeff=0.0_flyt
        k=0
        do j=1,gs%quartet%nvar
        do i=1,gs%poly%ncoeff
            k=k+1
            gs%quartet%coeff(i,j)=solution(k)
        enddo
        enddo
        gs%quartet%nconstr=map%constraints%neq4
        if ( gs%quartet%nconstr .gt. 0 ) then
            allocate(gs%quartet%constraints,source=map%constraints%eq4)
        endif
    end block lsq
end subroutine
