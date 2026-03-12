
!> get secondorder coefficients, but with gridfit instead.
subroutine solve_gridfit(gs,map,mw,mem,verbosity)
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
        type(lo_forceconstant_secondorder) :: fc,fcss
        type(lo_forceconstant_thirdorder) :: fct,fctss
        type(lo_forceconstant_fourthorder) :: fcf,fcfss
        type(lo_mdsim) :: sim
        real(flyt), dimension(:,:,:,:), allocatable :: dipole_forceconstant
        real(flyt), dimension(:,:), allocatable :: f,CM,wC,wA,wB
        real(flyt), dimension(:,:), allocatable :: theta_pair,theta_triplet,theta_quartet
        real(flyt), dimension(:), allocatable :: vB,vD,vsol
        real(flyt), dimension(3,3,3,3) :: m4
        real(flyt), dimension(3,3,3) :: m3
        real(flyt), dimension(3) :: v0,u2,u3,u4
        real(flyt) :: timer_coeff,t0
        integer, dimension(gs%ndim) :: maxord
        integer :: ii,i,j,t,l,nu,nf,nc,nt,a1,a2,a3,a4,i1,i2,i3,i4
        logical, dimension(81) :: rel_quartet_ntheta

        timer_coeff=walltime()
        t0=timer_coeff

        lo_allocate(f(3,map%n_atom_ss))
        f=0.0_flyt
        ! space for solution
        if ( map%have_fc_pair    ) allocate(theta_pair(gs%nsim,map%xuc%nx_fc_pair))
        if ( map%have_fc_triplet ) allocate(theta_triplet(gs%nsim,map%xuc%nx_fc_triplet))
        if ( map%have_fc_quartet ) allocate(theta_quartet(gs%nsim,map%xuc%nx_fc_quartet))
        if ( map%have_fc_pair    ) theta_pair=0.0_flyt
        if ( map%have_fc_triplet ) theta_triplet=0.0_flyt
        if ( map%have_fc_quartet ) theta_quartet=0.0_flyt
        if ( map%have_fc_quartet ) then
            rel_quartet_ntheta=.false.
            do i=1,map%n_fc_quartet_shell
                j=map%fc_quartet_shell(i)%nx
                if ( j .gt. 0 ) rel_quartet_ntheta(j)=.true.
            enddo
        endif

        do ii=1,gs%nsim
            if ( mod(ii,mw%n) .ne. mw%r ) cycle

            ! read a simulation
            call sim%read_from_hdf5(trim(filenames(ii)),verbosity=0)
            nf=sim%na*3
            nt=sim%nt*nf
            ! get structures and constraints
            call uc%generate(sim%extra%unitcell_latticevectors, sim%extra%unitcell_positions, &
                             sim%extra%unitcell_atomic_numbers, 2 )
            call ss%generate(sim%extra%supercell_latticevectors,sim%extra%supercell_positions, &
                             sim%extra%supercell_atomic_numbers, enhet=2 )
            call ss%classify('supercell',uc)

            ! Should probably subtract dipole-dipole interactions here
            if ( map%polar .gt. 0 ) then
                do i=1,gs%polar%nZ
                    map%xuc%x_Z_singlet(i)=gs%polar%ipZ%eval(i,gs%grid_coordinates(:,ii) )
                    !map%theta_Z(i)=gs%polar%ipZ%eval(i,gs%grid_coordinates(:,ii) )
                enddo
                do i=1,gs%polar%neps
                    map%xuc%x_eps_global(i)=gs%polar%ipeps%eval(i,gs%grid_coordinates(:,ii) )
                    !map%theta_eps(i)=gs%polar%ipeps%eval(i,gs%grid_coordinates(:,ii) )
                enddo
                call map%get_secondorder_forceconstant(uc,fc,mem,-1)
                lo_allocate(dipole_forceconstant(3,3,map%n_atom_ss,map%n_atom_ss))
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
                        call lo_stop_gracefully(['dipole-dipole forces do not add up to zero'],&
                                                lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
                    endif
                    sim%f(:,:,t)=sim%f(:,:,t)-f
                enddo
                lo_deallocate(dipole_forceconstant)
            endif

            ! Now fit and subtract the second order
            if ( map%have_fc_pair ) then
                ! constraints
                ! call lo_secondorder_rot_herm_huang(map,uc,wC,nc,.true.,.true.,.true.)
                nu=map%xuc%nx_fc_pair
                ! first get the coefficient matrix
                lo_allocate(wA(nt,nu))
                lo_allocate(vB(nt))
                lo_allocate(CM(nf,nu))
                lo_allocate(vsol(nu))
                wA=0.0_flyt
                vB=0.0_flyt
                CM=0.0_flyt
                vsol=0.0_flyt
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
                        call lo_stop_gracefully(['dgglse exit status '//tochar(lo_status)],&
                             lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                    endif
                    lo_deallocate(vD)
                else
                    lo_allocate(wB(nt,1))
                    wB(:,1)=vB
                    call lo_dgels(wA,wB,info=lo_status)
                    if ( lo_status .ne. 0 ) then
                        call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],&
                             lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                    else
                        vsol=wB(1:nu,1)
                    endif
                    lo_deallocate(wB)
                endif
                ! store solution
                theta_pair(ii,:)=vsol
                ! Subtract those forces

                do i=1,gs%polar%nZ
                    map%xuc%x_Z_singlet(i)=gs%polar%ipZ%eval(i,gs%grid_coordinates(:,ii) )
                    !map%theta_Z(i)=gs%polar%ipZ%eval(i,gs%grid_coordinates(:,ii) )
                enddo
                do i=1,gs%polar%neps
                    map%xuc%x_eps_global(i)=gs%polar%ipeps%eval(i,gs%grid_coordinates(:,ii) )
                    !map%theta_eps(i)=gs%polar%ipeps%eval(i,gs%grid_coordinates(:,ii) )
                enddo
                map%xuc%x_fc_pair=vsol
                !map%ifc_pair=vsol
                call map%get_secondorder_forceconstant(uc,fc,mem,-1)
                call fc%remap(uc,ss,fcss)
                do t=1,sim%nt
                    f=0.0_flyt
                    do a1=1,sim%na
                        a2=fcss%atom(a1)%pair(i)%i2
                        v0=matmul(fcss%atom(a1)%pair(i)%m,sim%u(:,a2,t))
                        f(:,a1)=f(:,a1)-v0
                    enddo
                    sim%f(:,:,t)=sim%f(:,:,t)-f
                enddo

                if ( allocated(wA) )    deallocate(wA)
                if ( allocated(vB) )    deallocate(vB)
                if ( allocated(CM) )    deallocate(CM)
                if ( allocated(vsol) )  deallocate(vsol)
                if ( allocated(wC) )    deallocate(wC)

                ! Report?
                if ( verbosity .gt. 0 ) then
                     if ( walltime()-t0 .gt. 2.0_flyt ) then
                        call lo_looptimer('... coefficients',timer_coeff,walltime(),ii,gs%nsim)
                        t0=walltime()
                    endif
                endif
            endif

            ! Now fit and subtract the third order
            if ( map%have_fc_triplet ) then
                nu=map%xuc%nx_fc_triplet
                nc=map%constraints%neq3
                if ( nc .gt. 0 ) then
                    lo_allocate(wC(nc,nu))
                    wC=map%constraints%eq3
                endif
                ! first get the coefficient matrix
                lo_allocate(wA(nt,nu))
                lo_allocate(vB(nt))
                lo_allocate(CM(nf,nu))
                lo_allocate(vsol(nu))
                wA=0.0_flyt
                vB=0.0_flyt
                CM=0.0_flyt
                vsol=0.0_flyt
                l=0
                do t=1,sim%nt
                    call lo_coeffmatrix_triplet(sim%u(:,:,t),CM,map)
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
                        call lo_stop_gracefully(['dgglse exit status '//tochar(lo_status)],&
                             lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                    endif
                    lo_deallocate(vD)
                else
                    lo_allocate(wB(nt,1))
                    wB(:,1)=vB
                    call lo_dgels(wA,wB,info=lo_status)
                    if ( lo_status .ne. 0 ) then
                        call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],&
                             lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                    else
                        vsol=wB(1:nu,1)
                    endif
                    lo_deallocate(wB)
                endif
                ! store solution
                theta_triplet(ii,:)=vsol

                ! Subtract those forces
                map%xuc%x_fc_triplet=vsol
                !map%ifc_triplet=vsol
                call map%get_thirdorder_forceconstant(uc,fct)
                call fct%remap(uc,ss,fctss)
                do t=1,sim%nt
                    do a1=1,fctss%na
                        v0=0.0_flyt
                        do i=1,fctss%atom(a1)%n
                            m3=fctss%atom(a1)%triplet(i)%m
                            a2=fctss%atom(a1)%triplet(i)%i2
                            a3=fctss%atom(a1)%triplet(i)%i3
                            u2=sim%u(:,a2,t)
                            u3=sim%u(:,a3,t)
                            do i1=1,3
                            do i2=1,3
                            do i3=1,3
                                v0(i1)=v0(i1)-m3(i1,i2,i3)*u2(i2)*u3(i3)
                            enddo
                            enddo
                            enddo
                        enddo
                        f(:,a1)=v0*0.5_flyt
                    enddo
                    sim%f(:,:,t)=sim%f(:,:,t)-f
                enddo

                if ( allocated(wA) )    deallocate(wA)
                if ( allocated(vB) )    deallocate(vB)
                if ( allocated(CM) )    deallocate(CM)
                if ( allocated(vsol) )  deallocate(vsol)
                if ( allocated(wC) )    deallocate(wC)

                ! Report?
                if ( verbosity .gt. 0 ) then
                     if ( walltime()-t0 .gt. 2.0_flyt ) then
                        call lo_looptimer('... coefficients',timer_coeff,walltime(),ii,gs%nsim)
                        t0=walltime()
                    endif
                endif
            endif

            ! Now fit and subtract the fourth order
            if ( map%have_fc_quartet ) then
                nu=map%xuc%nx_fc_quartet
                nc=map%constraints%neq4
                if ( nc .gt. 0 ) then
                    lo_allocate(wC(nc,nu))
                    wC=map%constraints%eq4
                endif
                ! first get the coefficient matrix
                lo_allocate(wA(nt,nu))
                lo_allocate(vB(nt))
                lo_allocate(CM(nf,nu))
                lo_allocate(vsol(nu))
                wA=0.0_flyt
                vB=0.0_flyt
                CM=0.0_flyt
                vsol=0.0_flyt
                l=0
                do t=1,sim%nt
                    call lo_coeffmatrix_quartet(sim%u(:,:,t),CM,map,rel_quartet_ntheta)
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
                        call lo_stop_gracefully(['dgglse exit status '//tochar(lo_status)],&
                             lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                    endif
                    lo_deallocate(vD)
                else
                    lo_allocate(wB(nt,1))
                    wB(:,1)=vB
                    call lo_dgels(wA,wB,info=lo_status)
                    if ( lo_status .ne. 0 ) then
                        call lo_stop_gracefully(['dgels exit status '//tochar(lo_status)],&
                             lo_exitcode_blaslapack,__FILE__,__LINE__,mw%comm)
                    else
                        vsol=wB(1:nu,1)
                    endif
                    lo_deallocate(wB)
                endif
                ! store solution
                theta_quartet(ii,:)=vsol
                ! Subtract those forces

                !map%ifc_quartet=vsol
                map%xuc%x_fc_quartet=vsol
                call map%get_fourthorder_forceconstant(uc,fcf)
                call fcf%remap(uc,ss,fcfss)
                do t=1,sim%nt
                    do a1=1,fcfss%na
                        v0=0.0_flyt
                        do i=1,fcfss%atom(a1)%n
                            m4=fcfss%atom(a1)%quartet(i)%m
                            a2=fcfss%atom(a1)%quartet(i)%i2
                            a3=fcfss%atom(a1)%quartet(i)%i3
                            a4=fcfss%atom(a1)%quartet(i)%i4
                            u2=sim%u(:,a2,t)
                            u3=sim%u(:,a3,t)
                            u4=sim%u(:,a4,t)
                            do i1=1,3
                            do i2=1,3
                            do i3=1,3
                            do i4=1,3
                                v0(i1)=v0(i1)-m4(i1,i2,i3,i4)*u2(i2)*u3(i3)*u4(i4)
                            enddo
                            enddo
                            enddo
                            enddo
                        enddo
                        f(:,a1)=v0/6.0_flyt
                    enddo
                    sim%f(:,:,t)=sim%f(:,:,t)-f
                enddo

                if ( allocated(wA) )    deallocate(wA)
                if ( allocated(vB) )    deallocate(vB)
                if ( allocated(CM) )    deallocate(CM)
                if ( allocated(vsol) )  deallocate(vsol)
                if ( allocated(wC) )    deallocate(wC)

                ! Report?
                if ( verbosity .gt. 0 ) then
                     if ( walltime()-t0 .gt. 2.0_flyt ) then
                        call lo_looptimer('... coefficients',timer_coeff,walltime(),ii,gs%nsim)
                        t0=walltime()
                    endif
                endif
            endif
        enddo

        ! Add up all the solutions
        if ( map%have_fc_pair )    call mpi_allreduce(MPI_IN_PLACE,theta_pair,    gs%nsim*map%xuc%nx_fc_pair,    MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
        if ( map%have_fc_triplet ) call mpi_allreduce(MPI_IN_PLACE,theta_triplet, gs%nsim*map%xuc%nx_fc_triplet, MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
        if ( map%have_fc_quartet ) call mpi_allreduce(MPI_IN_PLACE,theta_quartet, gs%nsim*map%xuc%nx_fc_quartet, MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)

        ! And initialize the interpolation
        maxord=2
        if ( map%have_fc_pair )    call gs%gf%ip2%generate( gs%grid_coordinates, theta_pair,    2, gs%info%distance_scale, maxord )
        if ( map%have_fc_triplet ) call gs%gf%ip3%generate( gs%grid_coordinates, theta_triplet, 2, gs%info%distance_scale, maxord )
        if ( map%have_fc_quartet ) call gs%gf%ip4%generate( gs%grid_coordinates, theta_quartet, 2, gs%info%distance_scale, maxord )
    end block getcoeff

end subroutine
