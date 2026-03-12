
!> remove the polar forces from the input data.
subroutine subtract_polar_forces(gs,map,mw,mem,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    real(r8), dimension(:,:,:,:), allocatable :: dipole_forceconstant
    real(r8), dimension(:,:), allocatable :: f,force_rsq
    real(r8), dimension(:), allocatable :: phiZ,phieps
    real(r8), dimension(3) :: v0
    real(r8), parameter :: ewaldtol=1E-12_r8
    real(r8) :: tt0,energy
    integer :: gpoint,i,ii,a1,a2,t

    tt0=walltime()

    ! Start by making some space for things
    lo_allocate(dipole_forceconstant(3,3,map%n_atom_ss,map%n_atom_ss))
    lo_allocate(f(3,map%n_atom_ss))
    lo_allocate(force_rsq(2,gs%nsim))
    lo_allocate(phiZ(gs%polar%nZ))
    lo_allocate(phieps(gs%polar%neps))

    dipole_forceconstant=0.0_r8
    f=0.0_r8
    force_rsq=0.0_r8
    phiZ=0.0_r8
    phieps=0.0_r8

    if ( verbosity .gt. 0 ) call lo_progressbar_init()
    ! go through the relevant gridpoint, in order
    do gpoint=1,gs%raw%nrelevant_gridpoints
        ii=gs%raw%relevant_gridpoints(gpoint)
        ! get the dipole forceconstant at this point
        forceconst: block
            type(lo_crystalstructure) :: uc,ss
            type(lo_forceconstant_secondorder) :: fc
            real(r8), dimension(gs%ndim) :: transformed_gridcoord
            real(r8) :: f0,f1

            select case(map%polar)
            case(1)
                do i=1,gs%polar%nZ
                    map%xuc%x_Z_singlet(i)=gs%polar%ipZ%eval(i,gs%grid_coordinates(:,ii) )
                enddo
                do i=1,gs%polar%neps
                    map%xuc%x_eps_global(i)=gs%polar%ipeps%eval(i,gs%grid_coordinates(:,ii) )
                enddo
            case(2)
                call coordinate_transformation(gs,gs%grid_coordinates(:,ii),transformed_gridcoord)
                do i=1,map%xuc%nx_eps_global
                    f0=gs%polar%ipeps%eval(i,gs%grid_coordinates(:,ii) )
                    f1=gs%poly%eval(transformed_gridcoord,gs%eps_global%coeff(:,i))
                    map%xuc%x_eps_global(i)=f0+f1 !gs%poly%eval(transformed_gridcoord,gs%eps_global%coeff(:,i))
                enddo
                do i=1,map%xuc%nx_Z_singlet
                    map%xuc%x_Z_singlet(i)=gs%poly%eval(transformed_gridcoord,gs%Z_singlet%coeff(:,i))
                enddo
            end select

            ! Generate structures for this point
            call uc%generate( gs%ref(ii)%unitcell_latticevectors, gs%ref(ii)%unitcell_positions, &
                              gs%ref(ii)%unitcell_atomic_numbers ,enhet=2 )
            call ss%generate( gs%ref(ii)%supercell_latticevectors, gs%ref(ii)%supercell_positions, &
                              gs%ref(ii)%supercell_atomic_numbers, enhet=2 )
            call ss%classify('supercell',uc)
            ! Generate a fake forceconstant thingy
            call map%get_secondorder_forceconstant(uc,fc,mem,-1)
            ! Get the supercell dynamical matrix
            call fc%supercell_longrange_dynamical_matrix_at_gamma(ss,dipole_forceconstant,ewaldtol)
        end block forceconst

        ! subtract forces
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            f=0.0_r8
            energy=0.0_r8
            do a1=1,map%n_atom_ss
                do a2=1,map%n_atom_ss
                    v0=matmul(dipole_forceconstant(:,:,a1,a2),gs%raw%u(:,a2,t))
                    f(:,a1)=f(:,a1)-v0
                enddo
                energy=energy-dot_product(gs%raw%u(:,a1,t),f(:,a1))*0.5_r8
            enddo
            ! sanity check that they add up to zero
            if ( abs(sum(f)) .gt. lo_sqtol ) then
                call lo_stop_gracefully(['dipole-dipole forces do not add up to zero'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
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
            gs%raw%e_polar(t)=energy
        enddo

        if ( verbosity .gt. 0 .and. gpoint .lt. gs%raw%nrelevant_gridpoints ) then
            call lo_progressbar(' ... subtracting polar forces',gpoint,gs%raw%nrelevant_gridpoints,walltime()-tt0)
        endif
    enddo
    ! Add the R^2 things up over MPI
    call mpi_allreduce(MPI_IN_PLACE,force_rsq,gs%nsim*2,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    ! Store the R^2 somewhere? I do it in Z. Also free memory
    lo_allocate(gs%polar%rsquare(gs%nsim))
    do ii=1,gs%nsim
        gs%polar%rsquare(ii)=1.0_r8-force_rsq(1,ii)/force_rsq(2,ii)
    enddo
    if ( verbosity .gt. 0 ) call lo_progressbar(' ... subtracting polar forces',gs%nsim,gs%nsim,walltime()-tt0)
end subroutine

!> set up interpolation for the polar part.
subroutine create_polar_interpolation(gs,map,mw,mem,verbosity)
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

    type(lo_crystalstructure) :: uc
    type(lo_mpi_helper) :: ml
    real(r8), dimension(:,:), allocatable :: grid_Z,grid_eps
    real(r8) :: t0
    integer, dimension(gs%ndim) :: maxorder
    integer :: ii


    if ( verbosity .gt. 0 .and. mw%talk ) then
        write(*,*) ''
        write(*,*) 'SOLVING FOR BORN EFFECTIVE CHARGES AND DIELECTRIC TENSOR'
        t0=walltime()
    endif

    ! Split the MPI communicator to get non-parallel ones temporarily
    call mw%split(ml,mw%r,__FILE__,__LINE__)

    ! grab the raw data
    if ( map%xuc%nx_Z_singlet .gt. 0 ) then
        allocate(grid_Z(gs%nsim,map%xuc%nx_Z_singlet))
        grid_Z=0.0_r8
    endif
    allocate(grid_eps(gs%nsim,map%xuc%nx_eps_global))
    grid_eps=0.0_r8
    do ii=1,gs%nsim
        if ( mod(ii,mw%n) .ne. mw%r ) cycle
        call uc%generate( gs%ref(ii)%unitcell_latticevectors, gs%ref(ii)%unitcell_positions, &
                          gs%ref(ii)%unitcell_atomic_numbers ,enhet=2 )
        call lo_solve_for_borncharges(map,uc,Z=gs%ref(ii)%born_effective_charges,eps=gs%ref(ii)%dielectric_tensor,mw=ml,mem=mem,verbosity=0)
        ! these should already satisfy hermiticity, since I do that in pack simulation.
        if ( map%xuc%nx_Z_singlet .gt. 0 ) then
            grid_Z(ii,:)=map%xuc%x_Z_singlet
        endif
        grid_eps(ii,:)=map%xuc%x_eps_global
    enddo
    if ( map%xuc%nx_Z_singlet .gt. 0 ) then
        call mw%allreduce('sum',grid_Z)
    endif
    call mw%allreduce('sum',grid_eps)

    gs%polar%nZ=map%xuc%nx_Z_singlet
    gs%polar%neps=map%xuc%nx_eps_global

    ! set the polynomical coefficients:
    maxorder=3
    !if ( gs%info%dim_temperature .gt. 0 ) maxorder(gs%info%dim_temperature)=0
    if ( map%xuc%nx_Z_singlet .gt. 0 ) then
        call gs%polar%ipZ%generate( gs%grid_coordinates,grid_Z,2,0.05_r8,maxorder )
    endif
    call gs%polar%ipeps%generate( gs%grid_coordinates,grid_eps,2,0.05_r8,maxorder )

    ! And destroy the temporary communicators
    call ml%free(__FILE__,__LINE__)

    if ( verbosity .gt. 0 ) write(*,*) '... created polar interpolations (',tochar(walltime()-t0),'s)'
end subroutine

!> get dielectric interactions
subroutine solve_dielectric(gs,map,ss,mw,mem,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    real(r8) :: timer,t0,t1
    logical :: withdiff

    timer=walltime()
    t0=timer
    t1=timer

    init: block
        ! Initialize the coefficients to nothing, for now.
        if ( map%xuc%nx_eps_global .gt. 0 ) then
            allocate(gs%eps_global%coeff(gs%poly%ncoeff,map%xuc%nx_eps_global))
            gs%eps_global%coeff=0.0_r8
        endif
        if ( map%xuc%nx_eps_singlet .gt. 0 ) then
            allocate(gs%eps_singlet%coeff(gs%poly%ncoeff,map%xuc%nx_eps_singlet))
            gs%eps_singlet%coeff=0.0_r8
        endif
        if ( map%xuc%nx_eps_pair .gt. 0 ) then
            allocate(gs%eps_pair%coeff(gs%poly%ncoeff,map%xuc%nx_eps_pair))
            gs%eps_pair%coeff=0.0_r8
        endif
        if ( map%xuc%nx_Z_singlet .gt. 0 ) then
            allocate(gs%Z_singlet%coeff(gs%poly%ncoeff,map%xuc%nx_Z_singlet))
            gs%Z_singlet%coeff=0.0_r8
        endif
        if ( map%xuc%nx_Z_pair .gt. 0 ) then
            allocate(gs%Z_pair%coeff(gs%poly%ncoeff,map%xuc%nx_Z_pair))
            gs%Z_pair%coeff=0.0_r8
        endif
        if ( map%xuc%nx_Z_triplet .gt. 0 ) then
            allocate(gs%Z_triplet%coeff(gs%poly%ncoeff,map%xuc%nx_Z_triplet))
            gs%Z_triplet%coeff=0.0_r8
        endif
    end block init

    ! Chopped this into several parts. It is quite useful with the
    ! standard things as a baseline.
    !call create_polar_interpolation(gs,map,mw,mem,verbosity)
    ! Remove the baseline from the data
    !call subtract_eps_baseline(gs,map,verbosity)

    ! Create the eps fit iteratively.
    iterative: block
        integer, parameter :: niter=0 !-1 !15
        real(r8), parameter :: mixpar=0.5_r8
        type(lo_grid_interpolation) :: ipeps
        real(r8), dimension(:,:,:), allocatable :: buf_eps,buf_eps_baseline
        real(r8), dimension(:,:), allocatable :: buf_mix_sing,buf_mix_pair
        real(r8), dimension(3,3) :: I3,m0
        real(r8) :: volume
        integer, dimension(:), allocatable :: ctr
        integer :: i,j,k,l,gpoint,ii,t
        integer :: iter

        ! Will need temporary space for this
        allocate(buf_mix_sing(gs%poly%ncoeff,map%xuc%nx_eps_singlet))
        allocate(buf_mix_pair(gs%poly%ncoeff,map%xuc%nx_eps_pair))
        allocate(buf_eps(3,3,gs%raw%nconf))
        allocate(buf_eps_baseline(3,3,gs%nsim))
        allocate(ctr(gs%nsim))
        buf_mix_sing=0.0_r8
        buf_mix_pair=0.0_r8
        buf_eps=0.0_r8
        buf_eps_baseline=0.0_r8
        ctr=0

        I3=0.0_r8
        do i=1,3
            I3(i,i)=1.0_r8
        enddo

        ! Convert from eps to raw derivatives
        buf_eps=gs%raw%eps
        do gpoint=1,gs%raw%nrelevant_gridpoints
            ii=gs%raw%relevant_gridpoints(gpoint)
            volume=abs(lo_determ(gs%ref(ii)%supercell_latticevectors))
            do t=1,gs%raw%nconf
                if ( gs%raw%gridind(t) .ne. ii ) cycle
                m0=buf_eps(:,:,t)
                m0=(I3-m0)*volume*0.25_r8/lo_pi
                !m0=I3-m0*4*lo_pi/volume
                buf_eps(:,:,t)=m0
            enddo
        enddo

        ! Need to have polar interpolation for things not to crash. Will be reconstructed eventually.
        call create_polar_interpolation(gs,map,mw,mem,verbosity)
        ! Start with the zeroth iteration, fit to differences. Reset values first.
        gs%raw%eps=buf_eps
        if ( map%have_eps_singlet ) then
            call coeff_eps_singlet_diff(gs,map,mw,mem,verbosity)
            call subtract_eps_singlet(gs,map,mw,mem,verbosity)
        endif
        if ( map%have_eps_pair ) then
            call coeff_eps_pair_diff(gs,map,mw,mem,verbosity)
            call subtract_eps_pair(gs,map,mw,mem,verbosity)
        endif
        ! Calculate the average remaining eps
        ctr=0
        buf_eps_baseline=0.0_r8
        do i=1,gs%raw%nconf
            j=gs%raw%gridind(i)
            ctr(j)=ctr(j)+1
            buf_eps_baseline(:,:,j)=buf_eps_baseline(:,:,j)+gs%raw%eps(:,:,i)
        enddo
        call mw%allreduce('sum',buf_eps_baseline)
        call mw%allreduce('sum',ctr)
        do j=1,gs%nsim
            buf_eps_baseline(:,:,j)=buf_eps_baseline(:,:,j)/real(ctr(j),r8)
        enddo
        call dummy_eps_interpolation(gs,map,buf_eps_baseline,ipeps,mw,mem)

iter=0
if ( mw%talk ) then
    write(*,*) 'did iter',iter
    do i=1,size(gs%eps_singlet%coeff,2)
        write(*,*) 'sing',i,gs%eps_singlet%coeff(:,i)
    enddo
    do i=1,size(gs%eps_pair%coeff,2)
        write(*,*) 'pair',i,gs%eps_pair%coeff(:,i)
    enddo
endif

        ! Now start iterating a few times? Maybe makes sense.
        iterloop: do iter=1,niter
            ! Reset eps values
            gs%raw%eps=buf_eps
            ! Remove current global fit
            call subtract_dummy_eps_baseline(gs,map,ipeps)
            ! New detailed fit
            if ( map%have_eps_singlet ) then
                buf_mix_sing=gs%eps_singlet%coeff
                call coeff_eps_singlet(gs,map,mw,mem,verbosity-1)

                call subtract_eps_singlet(gs,map,mw,mem,verbosity-1)
            endif
            if ( map%have_eps_pair ) then
                buf_mix_pair=gs%eps_pair%coeff
                call coeff_eps_pair(gs,map,mw,mem,verbosity-1)
if ( mw%talk ) write(*,*) 'diff:',norm2(gs%eps_pair%coeff-buf_mix_pair),norm2(gs%eps_pair%coeff-buf_mix_pair)/norm2(gs%eps_pair%coeff)
                gs%eps_pair%coeff=gs%eps_pair%coeff*mixpar+buf_mix_pair*(1.0_r8-mixpar)
            endif

            ! Get new residual, first reset then remove
            gs%raw%eps=buf_eps
            if ( map%have_eps_singlet ) then
                call subtract_eps_singlet(gs,map,mw,mem,verbosity)
            endif
            if ( map%have_eps_pair ) then
                call subtract_eps_pair(gs,map,mw,mem,verbosity)
            endif
            ! Get a new baseline fit.
            ctr=0
            buf_eps_baseline=0.0_r8
            do i=1,gs%raw%nconf
                j=gs%raw%gridind(i)
                ctr(j)=ctr(j)+1
                buf_eps_baseline(:,:,j)=buf_eps_baseline(:,:,j)+gs%raw%eps(:,:,i)
            enddo
            call mw%allreduce('sum',buf_eps_baseline)
            call mw%allreduce('sum',ctr)
            do j=1,gs%nsim
                buf_eps_baseline(:,:,j)=buf_eps_baseline(:,:,j)/real(ctr(j),r8)
            enddo
            call dummy_eps_interpolation(gs,map,buf_eps_baseline,ipeps,mw,mem)

if ( mw%talk ) then
    write(*,*) 'did iter',iter
    do i=1,size(gs%eps_singlet%coeff,2)
        write(*,*) 'sing',i,gs%eps_singlet%coeff(:,i)
    enddo
    do i=1,size(gs%eps_pair%coeff,2)
        write(*,*) 'pair',i,gs%eps_pair%coeff(:,i)
    enddo
endif
        enddo iterloop

        ! Final touch, convert the residual back to eps from derivatives
        do gpoint=1,gs%raw%nrelevant_gridpoints
            ii=gs%raw%relevant_gridpoints(gpoint)
            volume=abs(lo_determ(gs%ref(ii)%supercell_latticevectors))
            do t=1,gs%raw%nconf
                if ( gs%raw%gridind(t) .ne. ii ) cycle
                m0=gs%raw%eps(:,:,t)
                !m0=(I3-m0)*volume*0.25_r8/lo_pi
                m0=I3-m0*4*lo_pi/volume
                gs%raw%eps(:,:,t)=m0
            enddo
        enddo
        ! And a new baseline fit.
        ctr=0
        buf_eps_baseline=0.0_r8
        do i=1,gs%raw%nconf
            j=gs%raw%gridind(i)
            ctr(j)=ctr(j)+1
            buf_eps_baseline(:,:,j)=buf_eps_baseline(:,:,j)+gs%raw%eps(:,:,i)
        enddo
        call mw%allreduce('sum',buf_eps_baseline)
        call mw%allreduce('sum',ctr)
        do j=1,gs%nsim
            gs%ref(j)%dielectric_tensor=buf_eps_baseline(:,:,j)/real(ctr(j),r8)
        enddo
        call create_polar_interpolation(gs,map,mw,mem,verbosity)
    end block iterative

    ! Then we do the Born charge stuff
    if ( map%xuc%nx_Z_singlet .gt. 0 ) then
        call coeff_Z_singlet(gs,map,ss,mw,mem,verbosity)
        call subtract_Z_singlet(gs,map,mw,mem,verbosity)
    endif

    if ( map%have_Z_pair ) then
        call coeff_Z_pair(gs,map,mw,mem,verbosity)
        call subtract_Z_pair(gs,map,mw,mem,verbosity)
    endif

    if ( map%have_Z_triplet ) then
        call coeff_Z_triplet(gs,map,mw,mem,verbosity)
    endif

    ! If we have both Born charges and dielectric tensor I should subtract
    ! the polar forces from the normal forces.
    if ( map%xuc%nx_Z_singlet .gt. 0 .and. map%xuc%nx_eps_global .gt. 0 ) then
        call gs%subtract_polar_forces(map,mw,mem,verbosity)
    endif
end subroutine

!> global epsilon
subroutine coeff_eps_global(gs,map,mw,mem,verbosity)
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

    real(r8), dimension(:,:), allocatable :: CTC,CTF
    real(r8) :: timer,t0,t1
    integer :: nf,nx,nphi,ncoeff,neq

    ! start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! set some simple things
    init: block
        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CREATING DIELECTRIC TENSOR'
        endif

        ! Number of variables to solve for
        nx=map%xuc%nx_eps_global
        ncoeff=gs%poly%ncoeff
        nphi=nx*ncoeff
        nf=9
        neq=gs%raw%nconf*nf
    end block init

    ! build coefficient matrices
    coeff: block
        type(lo_sparsematrix) :: sAM
        real(r8), dimension(:,:), allocatable :: partA,partB,partC
        integer :: i,j,k,l,ii,t

        ! Space for the sparse representation of the augmentation matrix
        call sAM%init(nrow=nx,ncol=nphi)

        ! Space for the thingy to solve
        call mem%allocate(CTC  ,[nphi,nphi],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(CTF  ,[nphi,1   ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partA,[nf,nx    ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partB,[nf,nphi  ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partC,[nf,1     ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        CTC=0.0_r8
        CTF=0.0_r8
        partA=0.0_r8
        partB=0.0_r8
        partC=0.0_r8

        do t=1,gs%raw%nconf
            ! which gridpoint are we on?
            ii=gs%raw%gridind(t)
            ! normal coefficient matrix
            partA=map%eps_global_shell%coeff
            ! Augment the coefficient matrix, first construct augmentation matrix
            sAM%val=lo_huge
            l=0
            do i=1,nx
            do j=1,ncoeff
                k=(i-1)*ncoeff+j
                l=l+1
                sAM%rowind(l)=i
                sAM%colind(l)=k
                sAM%val(l)=gs%poly%coeffM(ii,j)
            enddo
            enddo
            ! Matrix multiplicataion, but sparse, and manual. Probably fast enough anyway
            partB=0.0_r8
            do l=1,sAM%n
                k=sAM%rowind(l)
                j=sAM%colind(l)
                do i=1,nf
                    partB(i,j)=partB(i,j)+partA(i,k)*sAM%val(l)
                enddo
            enddo

            ! fetch epsilon
            do j=1,3
            do i=1,3
                k=(j-1)*3+i
                partC(k,1)=gs%raw%eps(i,j,t)
            enddo
            enddo
            ! Multadd this together!
            call lo_gemm(partB,partB,CTC,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            call lo_gemm(partB,partC,CTF,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            ! Report?
            if ( verbosity .gt. 0 ) then
            if ( walltime()-t0 .gt. timereport ) then
                call lo_looptimer('... global eps coefficients',timer,walltime(),t,gs%raw%nconf)
                t0=walltime()
            endif
            endif
        enddo

        t0=walltime()
        call mw%allreduce('sum',CTC)
        call mw%allreduce('sum',CTF)

        call mem%deallocate(partA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partB,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated matrices (',tochar(t1-t0),')'
            t0=t1
        endif
    end block coeff

    ! solve least squares problem
    slv: block
        real(r8), dimension(:), allocatable :: solution
        integer :: i,j,k

        call mem%allocate(solution,nphi,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        solution=0.0_r8
        if ( mw%r .eq. mw%n-1 ) then
            ! solve serially
            call lo_dgels(CTC,CTF,info=lo_status)
            if ( lo_status .ne. 0 ) then
                call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
            endif
            solution=CTF(:,1)
        endif
        call mw%bcast(solution,from=mw%n-1)

        ! store the solution
        gs%eps_global%nvar=nx
        gs%eps_global%nconstr=0
        gs%eps_global%coeff=0.0_r8
        k=0
        do j=1,nx
        do i=1,ncoeff
            k=k+1
            gs%eps_global%coeff(i,j)=solution(k)
        enddo
        enddo

        ! Cleanup
        call mem%deallocate(solution,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTF,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... solved for global epsilon (',tochar(t1-timer),')'
        endif
    end block slv
end subroutine

!> eps singlets
subroutine coeff_eps_singlet(gs,map,mw,mem,verbosity)
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

    real(r8), dimension(:,:), allocatable :: CTC,CTF
    real(r8) :: timer,t0,t1
    integer :: nf,nx,nphi,ncoeff,neq

    ! start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! set some simple things
    init: block
        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CREATING EPS SINGLET'
        endif

        ! Number of variables to solve for
        nx=map%xuc%nx_eps_singlet
        ncoeff=gs%poly%ncoeff
        nphi=nx*ncoeff
        nf=9
        neq=gs%raw%nconf*nf
    end block init

    ! build coefficient matrices
    coeff: block
        type(lo_sparsematrix) :: sAM
        real(r8), dimension(:,:), allocatable :: partA,partB,partC
        integer :: i,j,k,l,ii,t

        ! Space for the sparse representation of the augmentation matrix
        call sAM%init(nrow=nx,ncol=nphi)

        ! Space for the thingy to solve
        call mem%allocate(CTC  ,[nphi,nphi],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(CTF  ,[nphi,1   ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partA,[nf,nx    ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partB,[nf,nphi  ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partC,[nf,1     ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        CTC=0.0_r8
        CTF=0.0_r8
        partA=0.0_r8
        partB=0.0_r8
        partC=0.0_r8

        do t=1,gs%raw%nconf
            ! which gridpoint are we on?
            ii=gs%raw%gridind(t)
            ! normal coefficient matrix
            call lo_coeffmatrix_eps_singlet(map,gs%raw%u(:,:,t),partA)

            ! Augment the coefficient matrix, first construct augmentation matrix
            sAM%val=lo_huge
            l=0
            do i=1,nx
            do j=1,ncoeff
                k=(i-1)*ncoeff+j
                l=l+1
                sAM%rowind(l)=i
                sAM%colind(l)=k
                sAM%val(l)=gs%poly%coeffM(ii,j)
            enddo
            enddo
            ! Matrix multiplicataion, but sparse, and manual. Probably fast enough anyway
            partB=0.0_r8
            do l=1,sAM%n
                k=sAM%rowind(l)
                j=sAM%colind(l)
                do i=1,nf
                    partB(i,j)=partB(i,j)+partA(i,k)*sAM%val(l)
                enddo
            enddo

            ! fetch
            do j=1,3
            do i=1,3
                k=(j-1)*3+i
                partC(k,1)=gs%raw%eps(i,j,t)
            enddo
            enddo
            ! Multadd this together!
            call lo_gemm(partB,partB,CTC,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            call lo_gemm(partB,partC,CTF,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            ! Report?
            if ( verbosity .gt. 0 ) then
            if ( walltime()-t0 .gt. timereport ) then
                call lo_looptimer('... eps singlet',timer,walltime(),t,gs%raw%nconf)
                t0=walltime()
            endif
            endif
        enddo

        t0=walltime()
        call mw%allreduce('sum',CTC)
        call mw%allreduce('sum',CTF)

        call mem%deallocate(partA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partB,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated matrices (',tochar(t1-t0),')'
            t0=t1
        endif
    end block coeff

    ! solve least squares problem
    slv: block
        real(r8), dimension(:), allocatable :: solution
        integer :: i,j,k

        call mem%allocate(solution,nphi,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        solution=0.0_r8
        if ( mw%r .eq. mw%n-1 ) then
            ! solve serially
            call lo_dgels(CTC,CTF,info=lo_status)
            if ( lo_status .ne. 0 ) then
                call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
            endif
            solution=CTF(:,1)
        endif
        call mw%bcast(solution,from=mw%n-1)

        ! store the solution
        gs%eps_singlet%nvar=nx
        gs%eps_singlet%nconstr=0
        gs%eps_singlet%coeff=0.0_r8
        k=0
        do j=1,nx
        do i=1,ncoeff
            k=k+1
            gs%eps_singlet%coeff(i,j)=solution(k)
        enddo
        enddo

        ! Cleanup
        call mem%deallocate(solution,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTF,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... solved for eps singlet (',tochar(t1-timer),')'
        endif
    end block slv
end subroutine

!> eps singlets from differences
subroutine coeff_eps_singlet_diff(gs,map,mw,mem,verbosity)
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

    real(r8), dimension(:,:), allocatable :: CTC,CTF
    real(r8) :: timer,t0,t1
    integer :: nf,nx,nphi,ncoeff,neq

    ! start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! set some simple things
    init: block

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CREATING EPS SINGLET'
        endif

        ! Number of variables to solve for
        nx=map%xuc%nx_eps_singlet
        ncoeff=gs%poly%ncoeff
        nphi=nx*ncoeff
        nf=9
        neq=gs%raw%nconf*nf
    end block init

    ! build coefficient matrices
    coeff: block
        type(lo_sparsematrix) :: sAM
        real(r8), dimension(:,:,:), allocatable :: bufU,bufeps,bufcoeff
        real(r8), dimension(:,:), allocatable :: bufC
        real(r8), dimension(:,:), allocatable :: partA,partB,partC
        integer, dimension(:), allocatable :: ctr,offset
        integer :: i,j,k,l,ii,t,i1,i2,i3

        ! Space for the sparse representation of the augmentation matrix
        call sAM%init(nrow=nx,ncol=nphi)

        ! Space for the thingy to solve
        call mem%allocate(CTC  ,[nphi,nphi],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(CTF  ,[nphi,1   ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partA,[nf,nx    ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partB,[nf,nphi  ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partC,[nf,1     ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        CTC=0.0_r8
        CTF=0.0_r8
        partA=0.0_r8
        partB=0.0_r8
        partC=0.0_r8

        ! Another version: do the fit on differences?
        allocate(ctr(mw%n))
        allocate(offset(mw%n))

        t0=walltime()
        do ii=1,gs%nsim
            ! collect data from this sim to everyone. Count steps per rank and offset.
            ctr=0
            do t=1,gs%raw%nconf
                if ( gs%raw%gridind(t) .ne. ii ) cycle
                ctr(mw%r+1)=ctr(mw%r+1)+1
            enddo
            call mw%allreduce('sum',ctr)
            offset=0
            j=0
            do i=1,mw%n
                offset(i)=j
                j=j+ctr(i)
            enddo

            ! make some temporary space
            allocate(bufU(3,map%n_atom_ss,sum(ctr)))
            allocate(bufeps(3,3,sum(ctr)))
            allocate(bufcoeff(nf,nx,sum(ctr)))
            allocate(bufC(nf,sum(ctr)))
            bufU=0.0_r8
            bufeps=0.0_r8
            bufcoeff=0.0_r8
            bufC=0.0_r8

            ! collect displacements and eps
            i=offset(mw%r+1)
            do t=1,gs%raw%nconf
                if ( gs%raw%gridind(t) .ne. ii ) cycle
                i=i+1
                bufU(:,:,i)=gs%raw%u(:,:,t)
                bufeps(:,:,i)=gs%raw%eps(:,:,t)
            enddo
            call mw%allreduce('sum',bufU)
            call mw%allreduce('sum',bufeps)

            ! initialize sparse matrix for multiplication
            sAM%val=lo_huge
            l=0
            do i=1,nx
            do j=1,ncoeff
                k=(i-1)*ncoeff+j
                l=l+1
                sAM%rowind(l)=i
                sAM%colind(l)=k
                sAM%val(l)=gs%poly%coeffM(ii,j)
            enddo
            enddo

            ! Get coefficients
            bufcoeff=0.0_r8
            bufC=0.0_r8
            do t=1,sum(ctr)
                if ( mod(t,mw%n) .ne. mw%r ) cycle
                call lo_coeffmatrix_eps_singlet(map,bufU(:,:,t),bufcoeff(:,:,t))
                do j=1,3
                do i=1,3
                    k=(j-1)*3+i
                    bufC(k,t)=bufeps(i,j,t)
                enddo
                enddo
            enddo
            call mw%allreduce('sum',bufcoeff)
            call mw%allreduce('sum',bufC)

            ! Build all pairs and accumulate
            i3=0
            do i1=1,sum(ctr)
            do i2=i1+1,sum(ctr)
                i3=i3+1
                if ( mod(i3,mw%n) .ne. mw%r ) cycle
                ! Get difference
                partA=bufcoeff(:,:,i1)-bufcoeff(:,:,i2)
                partC(:,1)=bufC(:,i1)-bufC(:,i2)
                ! Sparse multiplication
                partB=0.0_r8
                do l=1,sAM%n
                    k=sAM%rowind(l)
                    j=sAM%colind(l)
                    do i=1,nf
                        partB(i,j)=partB(i,j)+partA(i,k)*sAM%val(l)
                    enddo
                enddo
                ! add together
                call lo_gemm(partB,partB,CTC,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
                call lo_gemm(partB,partC,CTF,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            enddo
            enddo

            ! Cleanup for the next round
            deallocate(bufU)
            deallocate(bufeps)
            deallocate(bufcoeff)
            deallocate(bufC)

            ! Report?
            if ( verbosity .gt. 0 ) then
            if ( walltime()-t0 .gt. timereport ) then
                call lo_looptimer('... eps singlet',timer,walltime(),t,gs%raw%nconf)
                t0=walltime()
            endif
            endif
        enddo

        ! Add everything together
        call mw%allreduce('sum',CTC)
        call mw%allreduce('sum',CTF)

        call mem%deallocate(partA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partB,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated matrices (',tochar(t1-t0),')'
            t0=t1
        endif
    end block coeff

    ! solve least squares problem
    slv: block
        real(r8), dimension(:), allocatable :: solution
        integer :: i,j,k

        call mem%allocate(solution,nphi,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        solution=0.0_r8
        if ( mw%r .eq. mw%n-1 ) then
            ! solve serially
            call lo_dgels(CTC,CTF,info=lo_status)
            if ( lo_status .ne. 0 ) then
                call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
            endif
            solution=CTF(:,1)
        endif
        call mw%bcast(solution,from=mw%n-1)

        ! store the solution
        gs%eps_singlet%nvar=nx
        gs%eps_singlet%nconstr=0
        gs%eps_singlet%coeff=0.0_r8
        k=0
        do j=1,nx
        do i=1,ncoeff
            k=k+1
            gs%eps_singlet%coeff(i,j)=solution(k)
        enddo
        enddo

        ! Cleanup
        call mem%deallocate(solution,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTF,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... solved for eps singlet (',tochar(t1-timer),')'
        endif
    end block slv
end subroutine

!> eps pairs
subroutine coeff_eps_pair(gs,map,mw,mem,verbosity)
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

    real(r8), dimension(:,:), allocatable :: CTC,CTF
    real(r8) :: timer,t0,t1
    integer :: nf,nx,nphi,ncoeff,neq

    ! start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! set some simple things
    init: block

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CREATING EPS PAIR'
        endif

        ! Number of variables to solve for
        nx=map%xuc%nx_eps_pair
        ncoeff=gs%poly%ncoeff
        nphi=nx*ncoeff
        nf=9
        neq=gs%raw%nconf*nf
    end block init

    ! build coefficient matrices
    coeff: block
        type(lo_sparsematrix) :: sAM
        real(r8), dimension(:,:), allocatable :: partA,partB,partC
        integer :: i,j,k,l,ii,t

        ! Space for the sparse representation of the augmentation matrix
        call sAM%init(nrow=nx,ncol=nphi)

        ! Space for the thingy to solve
        call mem%allocate(CTC  ,[nphi,nphi],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(CTF  ,[nphi,1   ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partA,[nf,nx    ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partB,[nf,nphi  ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partC,[nf,1     ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        CTC=0.0_r8
        CTF=0.0_r8
        partA=0.0_r8
        partB=0.0_r8
        partC=0.0_r8

        ! Normal version, fit to values only
        do t=1,gs%raw%nconf
            ! which gridpoint are we on?
            ii=gs%raw%gridind(t)
            ! normal coefficient matrix
            call lo_coeffmatrix_eps_pair(map,gs%raw%u(:,:,t),partA)

            ! Augment the coefficient matrix, first construct augmentation matrix
            sAM%val=lo_huge
            l=0
            do i=1,nx
            do j=1,ncoeff
                k=(i-1)*ncoeff+j
                l=l+1
                sAM%rowind(l)=i
                sAM%colind(l)=k
                sAM%val(l)=gs%poly%coeffM(ii,j)
            enddo
            enddo
            ! Matrix multiplicataion, but sparse, and manual. Probably fast enough anyway
            partB=0.0_r8
            do l=1,sAM%n
                k=sAM%rowind(l)
                j=sAM%colind(l)
                do i=1,nf
                    partB(i,j)=partB(i,j)+partA(i,k)*sAM%val(l)
                enddo
            enddo

            ! fetch
            do j=1,3
            do i=1,3
                k=(j-1)*3+i
                partC(k,1)=gs%raw%eps(i,j,t)
            enddo
            enddo

            ! Multadd this together!
            call lo_gemm(partB,partB,CTC,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            call lo_gemm(partB,partC,CTF,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            ! Report?
            if ( verbosity .gt. 0 ) then
            if ( walltime()-t0 .gt. timereport ) then
                call lo_looptimer('... eps pair',timer,walltime(),t,gs%raw%nconf)
                t0=walltime()
            endif
            endif
        enddo

        t0=walltime()
        call mw%allreduce('sum',CTC)
        call mw%allreduce('sum',CTF)

        call mem%deallocate(partA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partB,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated matrices (',tochar(t1-t0),')'
            t0=t1
        endif
    end block coeff

    ! solve least squares problem
    slv: block
        real(r8), dimension(:), allocatable :: solution
        integer :: i,j,k

        call mem%allocate(solution,nphi,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        solution=0.0_r8
        if ( mw%r .eq. mw%n-1 ) then
            ! solve serially
            call lo_dgels(CTC,CTF,info=lo_status)
            if ( lo_status .ne. 0 ) then
                call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
            endif
            solution=CTF(:,1)
        endif
        call mw%bcast(solution,from=mw%n-1)

        ! store the solution
        gs%eps_pair%nvar=nx
        gs%eps_pair%nconstr=0
        gs%eps_pair%coeff=0.0_r8
        k=0
        do j=1,nx
        do i=1,ncoeff
            k=k+1
            gs%eps_pair%coeff(i,j)=solution(k)
        enddo
        enddo

        ! Cleanup
        call mem%deallocate(solution,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTF,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... solved for eps pair (',tochar(t1-timer),')'
        endif
    end block slv
end subroutine

!> eps pairs from differences
subroutine coeff_eps_pair_diff(gs,map,mw,mem,verbosity)
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

    real(r8), dimension(:,:), allocatable :: CTC,CTF
    real(r8) :: timer,t0,t1
    integer :: nf,nx,nphi,ncoeff,neq

    ! start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! set some simple things
    init: block

        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CREATING EPS PAIR'
        endif

        ! Number of variables to solve for
        nx=map%xuc%nx_eps_pair
        ncoeff=gs%poly%ncoeff
        nphi=nx*ncoeff
        nf=9
        neq=gs%raw%nconf*nf
    end block init

    ! build coefficient matrices
    coeff: block
        type(lo_sparsematrix) :: sAM
        real(r8), dimension(:,:,:), allocatable :: bufU,bufeps,bufcoeff
        real(r8), dimension(:,:), allocatable :: bufC
        real(r8), dimension(:,:), allocatable :: partA,partB,partC
        integer, dimension(:), allocatable :: ctr,offset
        integer :: i,j,k,l,ii,t,i1,i2,i3

        ! Space for the sparse representation of the augmentation matrix
        call sAM%init(nrow=nx,ncol=nphi)

        ! Space for the thingy to solve
        call mem%allocate(CTC  ,[nphi,nphi],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(CTF  ,[nphi,1   ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partA,[nf,nx    ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partB,[nf,nphi  ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partC,[nf,1     ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        CTC=0.0_r8
        CTF=0.0_r8
        partA=0.0_r8
        partB=0.0_r8
        partC=0.0_r8

        ! Another version: do the fit on differences?
        allocate(ctr(mw%n))
        allocate(offset(mw%n))

        t0=walltime()
        do ii=1,gs%nsim
            ! collect data from this sim to everyone. Count steps per rank and offset.
            ctr=0
            do t=1,gs%raw%nconf
                if ( gs%raw%gridind(t) .ne. ii ) cycle
                ctr(mw%r+1)=ctr(mw%r+1)+1
            enddo
            call mw%allreduce('sum',ctr)
            offset=0
            j=0
            do i=1,mw%n
                offset(i)=j
                j=j+ctr(i)
            enddo

            ! make some temporary space
            allocate(bufU(3,map%n_atom_ss,sum(ctr)))
            allocate(bufeps(3,3,sum(ctr)))
            allocate(bufcoeff(nf,nx,sum(ctr)))
            allocate(bufC(nf,sum(ctr)))
            bufU=0.0_r8
            bufeps=0.0_r8
            bufcoeff=0.0_r8
            bufC=0.0_r8

            ! collect displacements and eps
            i=offset(mw%r+1)
            do t=1,gs%raw%nconf
                if ( gs%raw%gridind(t) .ne. ii ) cycle
                i=i+1
                bufU(:,:,i)=gs%raw%u(:,:,t)
                bufeps(:,:,i)=gs%raw%eps(:,:,t)
            enddo
            call mw%allreduce('sum',bufU)
            call mw%allreduce('sum',bufeps)

            ! initialize sparse matrix for multiplication
            sAM%val=lo_huge
            l=0
            do i=1,nx
            do j=1,ncoeff
                k=(i-1)*ncoeff+j
                l=l+1
                sAM%rowind(l)=i
                sAM%colind(l)=k
                sAM%val(l)=gs%poly%coeffM(ii,j)
            enddo
            enddo

            ! Get coefficients
            bufcoeff=0.0_r8
            bufC=0.0_r8
            do t=1,sum(ctr)
                if ( mod(t,mw%n) .ne. mw%r ) cycle
                call lo_coeffmatrix_eps_pair(map,bufU(:,:,t),bufcoeff(:,:,t))
                do j=1,3
                do i=1,3
                    k=(j-1)*3+i
                    bufC(k,t)=bufeps(i,j,t)
                enddo
                enddo
            enddo
            call mw%allreduce('sum',bufcoeff)
            call mw%allreduce('sum',bufC)

            ! Build all pairs and accumulate
            i3=0
            do i1=1,sum(ctr)
            do i2=i1+1,sum(ctr)
                i3=i3+1
                if ( mod(i3,mw%n) .ne. mw%r ) cycle
                ! Get difference
                partA=bufcoeff(:,:,i1)-bufcoeff(:,:,i2)
                partC(:,1)=bufC(:,i1)-bufC(:,i2)
                ! Sparse multiplication
                partB=0.0_r8
                do l=1,sAM%n
                    k=sAM%rowind(l)
                    j=sAM%colind(l)
                    do i=1,nf
                        partB(i,j)=partB(i,j)+partA(i,k)*sAM%val(l)
                    enddo
                enddo
                ! add together
                call lo_gemm(partB,partB,CTC,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
                call lo_gemm(partB,partC,CTF,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            enddo
            enddo

            ! Cleanup for the next round
            deallocate(bufU)
            deallocate(bufeps)
            deallocate(bufcoeff)
            deallocate(bufC)

            ! Report?
            if ( verbosity .gt. 0 ) then
            if ( walltime()-t0 .gt. timereport ) then
                call lo_looptimer('... eps pair',timer,walltime(),t,gs%raw%nconf)
                t0=walltime()
            endif
            endif
        enddo

        ! Add everything together
        call mw%allreduce('sum',CTC)
        call mw%allreduce('sum',CTF)

        call mem%deallocate(partA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partB,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated matrices (',tochar(t1-t0),')'
            t0=t1
        endif
    end block coeff

    ! solve least squares problem
    slv: block
        real(r8), dimension(:), allocatable :: solution
        integer :: i,j,k

        call mem%allocate(solution,nphi,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        solution=0.0_r8
        if ( mw%r .eq. mw%n-1 ) then
            ! solve serially
            call lo_dgels(CTC,CTF,info=lo_status)
            if ( lo_status .ne. 0 ) then
                call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
            endif
            solution=CTF(:,1)
        endif
        call mw%bcast(solution,from=mw%n-1)

        ! store the solution
        gs%eps_pair%nvar=nx
        gs%eps_pair%nconstr=0
        gs%eps_pair%coeff=0.0_r8
        k=0
        do j=1,nx
        do i=1,ncoeff
            k=k+1
            gs%eps_pair%coeff(i,j)=solution(k)
        enddo
        enddo

        ! Cleanup
        call mem%deallocate(solution,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTF,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... solved for eps pair (',tochar(t1-timer),')'
        endif
    end block slv
end subroutine

!> given an interpolation for global epsilon, subtract it and also convert from epsilon to derivatives.
subroutine subtract_eps_baseline(gs,map,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> how much to talk
    integer, intent(in) :: verbosity

    real(r8), dimension(3,3) :: m0
    real(r8), dimension(9) :: dr
    real(r8) :: t0
    integer :: gpoint,ii,i,t

    t0=walltime()
    if ( verbosity .gt. 0 ) call lo_progressbar_init()

    do gpoint=1,gs%raw%nrelevant_gridpoints
        ii=gs%raw%relevant_gridpoints(gpoint)
        do i=1,gs%polar%neps
            map%xuc%x_eps_global(i)=gs%polar%ipeps%eval( i,gs%grid_coordinates(:,ii) )
        enddo
        dr=matmul(map%eps_global_shell%coeff,map%xuc%x_eps_global)
        m0=lo_unflatten_2tensor(dr)

        ! Subtract from the data
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            gs%raw%eps(:,:,t)=gs%raw%eps(:,:,t)-m0
        enddo

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... subtracting dielectric baseline',gpoint,gs%raw%nrelevant_gridpoints,walltime()-t0)
        endif
    enddo
end subroutine

!> given an interpolation for global epsilon, subtract it and also convert from epsilon to derivatives.
subroutine subtract_eps_global(gs,map,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> how much to talk
    integer, intent(in) :: verbosity

    real(r8), dimension(3,3) :: m0,m1,I3
    real(r8) :: t0,volume
    integer :: gpoint,ii,t

    t0=walltime()
    if ( verbosity .gt. 0 ) call lo_progressbar_init()

    I3=0.0_r8
    do ii=1,3
        I3(ii,ii)=1.0_r8
    enddo

    do gpoint=1,gs%raw%nrelevant_gridpoints
        ii=gs%raw%relevant_gridpoints(gpoint)

        ! Grab the volume at this gridpoint, volume of the supercell that is
        volume=abs(lo_determ(gs%ref(ii)%supercell_latticevectors))

        ! Interpolate the dielectric constant here.
        geteps: block
            real(r8), dimension(9) :: dr
            real(r8), dimension(gs%ndim) :: transformed_gridcoord
            integer :: i

            call coordinate_transformation(gs,gs%grid_coordinates(:,ii),transformed_gridcoord)
            do i=1,map%xuc%nx_eps_global
                map%xuc%x_eps_global(i)=gs%poly%eval(transformed_gridcoord,gs%eps_global%coeff(:,i))
            enddo
            dr=matmul(map%eps_global_shell%coeff,map%xuc%x_eps_global)
            m0=lo_unflatten_2tensor(dr)
            ! Convert from eps to raw derivatives
            m0=(I3-m0)*volume*0.25_r8/lo_pi
        end block geteps

        ! Subtract from the data, and simultaneously convert from dielectric tensor to derivatives
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            m1=gs%raw%eps(:,:,t)
            ! Convert from eps to raw derivatives
            m1=(I3-m1)*volume*0.25_r8/lo_pi
            ! Store the difference
            gs%raw%eps(:,:,t)=m1-m0
        enddo

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... subtracting dielectric tensor',gpoint,gs%raw%nrelevant_gridpoints,walltime()-t0)
        endif
    enddo
end subroutine

!> given an interpolation for global epsilon, subtract it and also convert from epsilon to derivatives.
subroutine subtract_eps_singlet(gs,map,mw,mem,verbosity)
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

    type(lo_dielectric_tensor) :: diss
    type(lo_mpi_helper) :: ml
    real(r8), dimension(3,3) :: m0
    real(r8) :: t0
    integer :: gpoint,ii,t,j,a1

    t0=walltime()
    if ( verbosity .gt. 0 ) call lo_progressbar_init()

    call mw%split(ml,mw%r,__FILE__,__LINE__)

    do gpoint=1,gs%raw%nrelevant_gridpoints
        ii=gs%raw%relevant_gridpoints(gpoint)

        ! Interpolate the dielectric interactions here
        getZ: block
            type(lo_dielectric_tensor) :: di
            type(lo_crystalstructure) :: uc,ss

            call uc%generate( gs%ref(ii)%unitcell_latticevectors,gs%ref(ii)%unitcell_positions,gs%ref(ii)%unitcell_atomic_numbers,enhet=2 )
            call ss%generate( gs%ref(ii)%supercell_latticevectors,gs%ref(ii)%supercell_positions,gs%ref(ii)%supercell_atomic_numbers,enhet=2 )
            call ss%classify( 'supercell',uc )
            call gs%eval( map,gs%grid_coordinates(:,ii) )
            call map%get_dielectric_tensors(uc,di)
            call di%remap(diss,uc,ss,ml,mem,-1)
        end block getZ

        ! Subtract from the original data
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            m0=0.0_r8
            do a1=1,diss%n_eps_singlet
            do j=1,3
                m0=m0+diss%eps_singlet(a1)%m(:,:,j)*gs%raw%u(j,a1,t)
            enddo
            enddo
            gs%raw%eps(:,:,t)=gs%raw%eps(:,:,t)-m0
        enddo

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... subtracting eps singlet pairs',gpoint,gs%raw%nrelevant_gridpoints,walltime()-t0)
        endif
    enddo

    call ml%free(__FILE__,__LINE__)
end subroutine

!> given an interpolation for global epsilon, subtract it and also convert from epsilon to derivatives.
subroutine subtract_eps_pair(gs,map,mw,mem,verbosity)
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

    type(lo_dielectric_tensor) :: diss
    type(lo_mpi_helper) :: ml
    real(r8), dimension(3,3) :: m0
    real(r8) :: t0
    integer :: gpoint,ii,t,i,i1,i2,a1,a2

    t0=walltime()
    if ( verbosity .gt. 0 ) call lo_progressbar_init()

    call mw%split(ml,mw%r,__FILE__,__LINE__)

    do gpoint=1,gs%raw%nrelevant_gridpoints
        ii=gs%raw%relevant_gridpoints(gpoint)

        ! Interpolate the dielectric interactions here
        getZ: block
            type(lo_dielectric_tensor) :: di
            type(lo_crystalstructure) :: uc,ss

            call uc%generate( gs%ref(ii)%unitcell_latticevectors,gs%ref(ii)%unitcell_positions,gs%ref(ii)%unitcell_atomic_numbers,enhet=2 )
            call ss%generate( gs%ref(ii)%supercell_latticevectors,gs%ref(ii)%supercell_positions,gs%ref(ii)%supercell_atomic_numbers,enhet=2 )
            call ss%classify( 'supercell',uc )
            call gs%eval( map,gs%grid_coordinates(:,ii) )
            call map%get_dielectric_tensors(uc,di)
            call di%remap(diss,uc,ss,ml,mem,-1)
        end block getZ

        ! Subtract from the original data
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            m0=0.0_r8
            do i=1,diss%n_eps_pair
                a1=diss%eps_pair(i)%a1
                a2=diss%eps_pair(i)%a2
                do i1=1,3
                do i2=1,3
                    m0=m0+diss%eps_pair(i)%m(:,:,i1,i2)*gs%raw%u(i1,a1,t)*gs%raw%u(i2,a2,t)
                enddo
                enddo
            enddo
            gs%raw%eps(:,:,t)=gs%raw%eps(:,:,t)-m0
        enddo

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... subtracting eps pairs',gpoint,gs%raw%nrelevant_gridpoints,walltime()-t0)
        endif
    enddo

    call ml%free(__FILE__,__LINE__)
end subroutine

!> born charges
subroutine coeff_Z_singlet(gs,map,ss,mw,mem,verbosity)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity

    real(r8), dimension(:,:), allocatable :: CTC,CTF
    real(r8) :: timer,t0,t1
    integer :: nf,nx,nphi,ncoeff,neq

    ! start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! set some simple things
    init: block
        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CREATING BORN CHARGES'
        endif

        ! Number of variables to solve for
        nx=map%xuc%nx_Z_singlet
        ncoeff=gs%poly%ncoeff
        nphi=nx*ncoeff
        nf=ss%na*9
        neq=gs%raw%nconf*nf
    end block init

    ! build coefficient matrices
    coeff: block
        type(lo_sparsematrix) :: sAM
        real(r8), dimension(:,:), allocatable :: partA,partB,partC
        integer :: i,j,k,l,ii,t

        ! Space for the sparse representation of the augmentation matrix
        call sAM%init(nrow=nx,ncol=nphi)

        ! Space for the thingy to solve
        call mem%allocate(CTC  ,[nphi,nphi],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(CTF  ,[nphi,1   ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partA,[nf,nx    ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partB,[nf,nphi  ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partC,[nf,1     ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        CTC=0.0_r8
        CTF=0.0_r8
        partA=0.0_r8
        partB=0.0_r8
        partC=0.0_r8

        do t=1,gs%raw%nconf
            ! which gridpoint are we on?
            ii=gs%raw%gridind(t)
            ! normal coefficient matrix
            call lo_coeffmatrix_supercell_Z_singlet(map,ss,partA)
            ! Augment the coefficient matrix, first construct augmentation matrix
            sAM%val=lo_huge
            l=0
            do i=1,nx
            do j=1,ncoeff
                k=(i-1)*ncoeff+j
                l=l+1
                sAM%rowind(l)=i
                sAM%colind(l)=k
                sAM%val(l)=gs%poly%coeffM(ii,j)
            enddo
            enddo
            ! Matrix multiplicataion, but sparse, and manual. Probably fast enough anyway
            partB=0.0_r8
            do l=1,sAM%n
                k=sAM%rowind(l)
                j=sAM%colind(l)
                do i=1,nf
                    partB(i,j)=partB(i,j)+partA(i,k)*sAM%val(l)
                enddo
            enddo

            ! fetch borncharge
            do k=1,map%n_atom_ss
            do j=1,3
            do i=1,3
                l=(k-1)*9+(j-1)*3+i
                partC(l,1)=gs%raw%Z(i,j,k,t)
            enddo
            enddo
            enddo
            ! Multadd this together!
            call lo_gemm(partB,partB,CTC,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            call lo_gemm(partB,partC,CTF,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            ! Report?
            if ( verbosity .gt. 0 ) then
            if ( walltime()-t0 .gt. timereport ) then
                call lo_looptimer('... born charge coefficients',timer,walltime(),t,gs%raw%nconf)
                t0=walltime()
            endif
            endif
        enddo

        t0=walltime()
        call mw%allreduce('sum',CTC)
        call mw%allreduce('sum',CTF)

        call mem%deallocate(partA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partB,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated matrices (',tochar(t1-t0),')'
            t0=t1
        endif
    end block coeff

    ! solve least squares problem
    slv: block
        real(r8), dimension(:), allocatable :: solution
        integer :: i,j,k

        call mem%allocate(solution,nphi,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        solution=0.0_r8
        if ( mw%r .eq. mw%n-1 ) then
            ! solve serially
            call lo_dgels(CTC,CTF,info=lo_status)
            if ( lo_status .ne. 0 ) then
                call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
            endif
            solution=CTF(:,1)
        endif
        call mw%bcast(solution,from=mw%n-1)

        ! store the solution
        gs%Z_singlet%nvar=nx
        gs%Z_singlet%nconstr=0
        gs%Z_singlet%coeff=0.0_r8
        k=0
        do j=1,nx
        do i=1,ncoeff
            k=k+1
            gs%Z_singlet%coeff(i,j)=solution(k)
        enddo
        enddo

        ! Cleanup
        call mem%deallocate(solution,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTF,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... solved for Born charges (',tochar(t1-timer),')'
        endif
    end block slv
end subroutine

!> born charge pairs
subroutine coeff_Z_pair(gs,map,mw,mem,verbosity)
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

    real(r8), dimension(:,:), allocatable :: CTC,CTF
    real(r8) :: timer,t0,t1
    integer :: nf,nx,nphi,ncoeff,neq

    ! start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! set some simple things
    init: block
        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CREATING BORN CHARGE PAIRS'
        endif

        ! Number of variables to solve for
        nx=map%xuc%nx_Z_pair
        ncoeff=gs%poly%ncoeff
        nphi=nx*ncoeff
        nf=map%n_atom_ss*9
        neq=gs%raw%nconf*nf
    end block init

    ! build coefficient matrices
    coeff: block
        type(lo_sparsematrix) :: sAM
        real(r8), dimension(:,:), allocatable :: partA,partB,partC
        integer :: i,j,k,l,ii,t

        ! Space for the sparse representation of the augmentation matrix
        call sAM%init(nrow=nx,ncol=nphi)

        ! Space for the thingy to solve
        call mem%allocate(CTC  ,[nphi,nphi],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(CTF  ,[nphi,1   ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partA,[nf,nx    ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partB,[nf,nphi  ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partC,[nf,1     ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        CTC=0.0_r8
        CTF=0.0_r8
        partA=0.0_r8
        partB=0.0_r8
        partC=0.0_r8

        do t=1,gs%raw%nconf
            ! which gridpoint are we on?
            ii=gs%raw%gridind(t)
            ! normal coefficient matrix
            call lo_coeffmatrix_Z_pair(map,gs%raw%u(:,:,t),partA)

            ! Augment the coefficient matrix, first construct augmentation matrix
            sAM%val=lo_huge
            l=0
            do i=1,nx
            do j=1,ncoeff
                k=(i-1)*ncoeff+j
                l=l+1
                sAM%rowind(l)=i
                sAM%colind(l)=k
                sAM%val(l)=gs%poly%coeffM(ii,j)
            enddo
            enddo
            ! Matrix multiplicataion, but sparse, and manual. Probably fast enough anyway
            partB=0.0_r8
            do l=1,sAM%n
                k=sAM%rowind(l)
                j=sAM%colind(l)
                do i=1,nf
                    partB(i,j)=partB(i,j)+partA(i,k)*sAM%val(l)
                enddo
            enddo

            ! fetch borncharge
            do k=1,map%n_atom_ss
            do j=1,3
            do i=1,3
                l=(k-1)*9+(j-1)*3+i
                partC(l,1)=gs%raw%Z(i,j,k,t)
            enddo
            enddo
            enddo
            ! Multadd this together!
            call lo_gemm(partB,partB,CTC,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            call lo_gemm(partB,partC,CTF,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            ! Report?
            if ( verbosity .gt. 0 ) then
            if ( walltime()-t0 .gt. timereport ) then
                call lo_looptimer('... born charge pairs',timer,walltime(),t,gs%raw%nconf)
                t0=walltime()
            endif
            endif
        enddo

        t0=walltime()
        call mw%allreduce('sum',CTC)
        call mw%allreduce('sum',CTF)

        call mem%deallocate(partA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partB,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated matrices (',tochar(t1-t0),')'
            t0=t1
        endif
    end block coeff

    ! solve least squares problem
    slv: block
        real(r8), dimension(:), allocatable :: solution
        integer :: i,j,k

        call mem%allocate(solution,nphi,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        solution=0.0_r8
        if ( mw%r .eq. mw%n-1 ) then
            ! solve serially
            call lo_dgels(CTC,CTF,info=lo_status)
            if ( lo_status .ne. 0 ) then
                call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
            endif
            solution=CTF(:,1)
        endif
        call mw%bcast(solution,from=mw%n-1)

        ! store the solution
        gs%Z_pair%nvar=nx
        gs%Z_pair%nconstr=0
        gs%Z_pair%coeff=0.0_r8
        k=0
        do j=1,nx
        do i=1,ncoeff
            k=k+1
            gs%Z_pair%coeff(i,j)=solution(k)
        enddo
        enddo

        ! Cleanup
        call mem%deallocate(solution,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTF,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... solved for Born charge pairs (',tochar(t1-timer),')'
        endif
    end block slv
end subroutine

!> born charge triplets
subroutine coeff_Z_triplet(gs,map,mw,mem,verbosity)
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

    real(r8), dimension(:,:), allocatable :: CTC,CTF,phi_constraints
    real(r8) :: timer,t0,t1
    integer :: nf,nx,nphi,ncoeff,neq,nconstr

    ! start timers
    timer=walltime()
    t0=timer
    t1=timer

    ! set some simple things
    init: block
        if ( verbosity .gt. 0 ) then
            write(*,*) ''
            write(*,*) 'CREATING BORN CHARGE TRIPLETS'
        endif

        ! Number of variables to solve for
        nx=map%xuc%nx_Z_triplet
        ncoeff=gs%poly%ncoeff
        nphi=nx*ncoeff
        nf=map%n_atom_ss*9
        neq=gs%raw%nconf*nf
        nconstr=map%constraints%neqz3

    end block init

    ! build coefficient matrices
    coeff: block
        type(lo_sparsematrix) :: sAM
        real(r8), dimension(:,:), allocatable :: partA,partB,partC,dfA
        integer :: i,j,k,l,ii,jj,t

        ! Space for the sparse representation of the augmentation matrix
        call sAM%init(nrow=nx,ncol=nphi)

        ! Space for the thingy to solve
        call mem%allocate(CTC  ,[nphi,nphi],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(CTF  ,[nphi,1   ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partA,[nf,nx    ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partB,[nf,nphi  ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(partC,[nf,1     ],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( nconstr .gt. 0 ) then
            call mem%allocate(dfA,[nconstr*gs%nsim,nphi],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            dfA=0.0_r8
        endif
        CTC=0.0_r8
        CTF=0.0_r8
        partA=0.0_r8
        partB=0.0_r8
        partC=0.0_r8

        do t=1,gs%raw%nconf
            ! which gridpoint are we on?
            ii=gs%raw%gridind(t)
            ! normal coefficient matrix
            call lo_coeffmatrix_Z_triplet(map,gs%raw%u(:,:,t),partA)

            ! Augment the coefficient matrix, first construct augmentation matrix
            sAM%val=lo_huge
            l=0
            do i=1,nx
            do j=1,ncoeff
                k=(i-1)*ncoeff+j
                l=l+1
                sAM%rowind(l)=i
                sAM%colind(l)=k
                sAM%val(l)=gs%poly%coeffM(ii,j)
            enddo
            enddo
            ! Matrix multiplicataion, but sparse, and manual. Probably fast enough anyway
            partB=0.0_r8
            do l=1,sAM%n
                k=sAM%rowind(l)
                j=sAM%colind(l)
                do i=1,nf
                    partB(i,j)=partB(i,j)+partA(i,k)*sAM%val(l)
                enddo
            enddo

            ! fetch borncharge
            do k=1,map%n_atom_ss
            do j=1,3
            do i=1,3
                l=(k-1)*9+(j-1)*3+i
                partC(l,1)=gs%raw%Z(i,j,k,t)
            enddo
            enddo
            enddo

            ! Maybe constraints?
            if ( nconstr .gt. 0 ) then
                ! Manual matrix multiplication again
                jj=(ii-1)*nconstr
                do l=1,sAM%n
                    k=sAM%rowind(l)
                    j=sAM%colind(l)
                    do i=1,nconstr
                        dfA(jj+i,j)=dfA(jj+i,j)+map%constraints%eqz3(i,k)*sAM%val(l)
                    enddo
                enddo
            endif

            ! Multadd this together!
            call lo_gemm(partB,partB,CTC,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            call lo_gemm(partB,partC,CTF,transa='T',transb='N',alpha=1.0_r8,beta=1.0_r8)
            ! Report?
            if ( verbosity .gt. 0 ) then
            if ( walltime()-t0 .gt. timereport ) then
                call lo_looptimer('... born charge triplets',timer,walltime(),t,gs%raw%nconf)
                t0=walltime()
            endif
            endif
        enddo

        t0=walltime()
        call mw%allreduce('sum',CTC)
        call mw%allreduce('sum',CTF)
        if ( nconstr .gt. 0 ) then
            call mw%allreduce('sum',dfA)
            call reduce_equations(dfA,phi_constraints,i)
        endif

        call mem%deallocate(partA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partB,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(partC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(dfA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... communicated matrices (',tochar(t1-t0),')'
            t0=t1
        endif
    end block coeff

    ! solve least squares problem
    slv: block
        real(flyt), dimension(:,:), allocatable :: wmA,wmB
        real(r8), dimension(:), allocatable :: solution
        integer :: i,j,k,nc

        call mem%allocate(solution,nphi,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        solution=0.0_r8
        if ( mw%r .eq. mw%n-1 ) then
            ! solve serially
            if ( nconstr .gt. 0 ) then
                ! Constrained solution
                nc=size(phi_constraints,2)
                call mem%allocate(wmA,[nphi+nc,nphi+nc],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                call mem%allocate(wmB,[nphi+nc,1],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                wmA=0.0_r8
                wmB=0.0_r8
                wmA(1:nphi,1:nphi)=CTC
                wmB(1:nphi,1)=CTF(:,1)
                do i=1,nphi
                do j=1,nc
                    wmA(nphi+j,i)=phi_constraints(i,j)
                    wmA(i,nphi+j)=phi_constraints(i,j)
                enddo
                enddo
                call lo_dgels(CTC,CTF,info=lo_status)
                if ( lo_status .ne. 0 ) then
                    call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                endif
                solution=CTF(1:nphi,1)
                call mem%deallocate(wmA,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
                call mem%deallocate(wmB,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            else
                ! Unconstrained solution
                call lo_dgels(CTC,CTF,info=lo_status)
                if ( lo_status .ne. 0 ) then
                    call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                endif
                solution=CTF(:,1)
            endif
        endif
        call mw%bcast(solution,from=mw%n-1)

        ! store the solution
        gs%Z_triplet%nvar=nx
        gs%Z_triplet%nconstr=0
        gs%Z_triplet%coeff=0.0_r8
        k=0
        do j=1,nx
        do i=1,ncoeff
            k=k+1
            gs%Z_triplet%coeff(i,j)=solution(k)
        enddo
        enddo

        ! Cleanup
        call mem%deallocate(solution,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTC,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(CTF,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(*,*) '... solved for Born charge triplets (',tochar(t1-timer),')'
        endif
    end block slv
end subroutine

!> given an interpolation for global epsilon, subtract it and also convert from epsilon to derivatives.
subroutine subtract_Z_singlet(gs,map,mw,mem,verbosity)
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

    type(lo_dielectric_tensor) :: diss
    type(lo_mpi_helper) :: ml
    real(r8) :: t0
    integer :: gpoint,ii,t,i

    call mw%split(ml,mw%r,__FILE__,__LINE__)

    t0=walltime()
    if ( verbosity .gt. 0 ) call lo_progressbar_init()

    do gpoint=1,gs%raw%nrelevant_gridpoints
        ii=gs%raw%relevant_gridpoints(gpoint)

        ! Interpolate the dielectric interactions here
        getZ: block
            type(lo_dielectric_tensor) :: di
            type(lo_crystalstructure) :: uc,ss

            call uc%generate( gs%ref(ii)%unitcell_latticevectors,gs%ref(ii)%unitcell_positions,gs%ref(ii)%unitcell_atomic_numbers,enhet=2 )
            call ss%generate( gs%ref(ii)%supercell_latticevectors,gs%ref(ii)%supercell_positions,gs%ref(ii)%supercell_atomic_numbers,enhet=2 )
            call ss%classify( 'supercell',uc )
            call gs%eval( map,gs%grid_coordinates(:,ii) )
            call map%get_dielectric_tensors(uc,di)
            call di%remap(diss,uc,ss,ml,mem,-1)
        end block getZ

        ! Subtract from the original data
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            do i=1,map%n_atom_ss
                gs%raw%Z(:,:,i,t)=gs%raw%Z(:,:,i,t)-diss%Z_singlet(:,:,i)
            enddo
        enddo

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... subtracting Born charges',gpoint,gs%raw%nrelevant_gridpoints,walltime()-t0)
        endif
    enddo

    call ml%free(__FILE__,__LINE__)
end subroutine

!> given an interpolation for global epsilon, subtract it and also convert from epsilon to derivatives.
subroutine subtract_Z_pair(gs,map,mw,mem,verbosity)
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

    type(lo_dielectric_tensor) :: diss
    type(lo_mpi_helper) :: ml
    real(r8), dimension(:,:,:), allocatable :: dz
    real(r8) :: t0
    integer :: gpoint,ii,t,i,j,a1,a2

    t0=walltime()
    if ( verbosity .gt. 0 ) call lo_progressbar_init()

    call mw%split(ml,mw%r,__FILE__,__LINE__)
    call mem%allocate(dz,[3,3,map%n_atom_ss],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dz=0.0_r8

    do gpoint=1,gs%raw%nrelevant_gridpoints
        ii=gs%raw%relevant_gridpoints(gpoint)

        ! Interpolate the dielectric interactions here
        getZ: block
            type(lo_dielectric_tensor) :: di
            type(lo_crystalstructure) :: uc,ss

            call uc%generate( gs%ref(ii)%unitcell_latticevectors,gs%ref(ii)%unitcell_positions,gs%ref(ii)%unitcell_atomic_numbers,enhet=2 )
            call ss%generate( gs%ref(ii)%supercell_latticevectors,gs%ref(ii)%supercell_positions,gs%ref(ii)%supercell_atomic_numbers,enhet=2 )
            call ss%classify( 'supercell',uc )
            call gs%eval( map,gs%grid_coordinates(:,ii) )
            call map%get_dielectric_tensors(uc,di)
            call di%remap(diss,uc,ss,ml,mem,-1)
        end block getZ

        ! Subtract from the original data
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            dz=0.0_r8
            do i=1,diss%n_Z_pair
                a1=diss%Z_pair(i)%a1
                a2=diss%Z_pair(i)%a2
                do j=1,3
                    dz(:,:,a1)=dz(:,:,a1)+diss%z_pair(i)%m(:,:,j)*gs%raw%u(j,a2,t)
                enddo
            enddo
            gs%raw%Z(:,:,:,t)=gs%raw%Z(:,:,:,t)-dz
        enddo

        if ( verbosity .gt. 0 ) then
            call lo_progressbar(' ... subtracting Born charge pairs',gpoint,gs%raw%nrelevant_gridpoints,walltime()-t0)
        endif
    enddo
    call mem%deallocate(dz,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call ml%free(__FILE__,__LINE__)
end subroutine

!> intermediate interpolation for the dielectric tensor
subroutine dummy_eps_interpolation(gs,map,bufeps,ipeps,mw,mem)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> buffer with eps
    real(r8), dimension(:,:,:), intent(in) :: bufeps
    !> interpolation
    type(lo_grid_interpolation), intent(out) :: ipeps
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:,:), allocatable :: grid_eps
    real(r8), dimension(9,map%xuc%nx_eps_global) :: wA
    real(r8), dimension(9) :: wB
    real(r8), dimension(map%xuc%nx_eps_global) :: wX
    integer, dimension(gs%ndim) :: maxorder
    integer :: ii

    allocate(grid_eps(gs%nsim,map%xuc%nx_eps_global))
    grid_eps=0.0_r8
    do ii=1,gs%nsim
        if ( mod(ii,mw%n) .ne. mw%r ) cycle
        wA=map%eps_global_shell%coeff
        wB=lo_flattentensor(bufeps(:,:,ii))
        call lo_linear_least_squares(wA,wB,wX)
        grid_eps(ii,:)=wX
    enddo
    call mw%allreduce('sum',grid_eps)

    ! set the polynomical coefficients:
    maxorder=3
    call ipeps%generate( gs%grid_coordinates,grid_eps,2,0.05_r8,maxorder )
    deallocate(grid_eps)
end subroutine

!> given an interpolation for global epsilon, subtract it and also convert from epsilon to derivatives.
subroutine subtract_dummy_eps_baseline(gs,map,ipeps)
    !> all the simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(inout) :: map
    !> interpolation
    type(lo_grid_interpolation), intent(inout) :: ipeps

    real(r8), dimension(3,3) :: m0
    real(r8), dimension(9) :: dr

    integer :: gpoint,ii,i,t

    do gpoint=1,gs%raw%nrelevant_gridpoints
        ii=gs%raw%relevant_gridpoints(gpoint)
        do i=1,map%xuc%nx_eps_global
            map%xuc%x_eps_global(i)=ipeps%eval( i,gs%grid_coordinates(:,ii) )
        enddo
        dr=matmul(map%eps_global_shell%coeff,map%xuc%x_eps_global)
        m0=lo_unflatten_2tensor(dr)
        ! Subtract from the data
        do t=1,gs%raw%nconf
            if ( gs%raw%gridind(t) .ne. ii ) cycle
            gs%raw%eps(:,:,t)=gs%raw%eps(:,:,t)-m0
        enddo
    enddo
end subroutine
