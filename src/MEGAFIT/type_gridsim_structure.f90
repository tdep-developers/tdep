
!> this is fast enough to do serially I believe.
subroutine create_structure_interpolation(gs,map,uc)
    !> grid of simulations
    class(lo_gridsim), intent(inout) :: gs
    !> forcemap
    type(lo_forcemap), intent(in) :: map
    !> prototype structure
    type(lo_crystalstructure), intent(in) :: uc

    ! internal brute force variant
    bfinternal: block
        real(flyt), dimension(:,:), allocatable :: gridvals
        real(flyt), dimension(3) :: v0
        integer, dimension(gs%ndim) :: maxorder
        integer :: ii,i,j

        allocate(gs%structure%atomic_numbers(uc%na))
        gs%structure%atomic_numbers=uc%atomic_number
        ! reference lattice vectors
        gs%structure%lv0=uc%latticevectors
        gs%structure%ilv0=uc%inv_latticevectors

        allocate(gs%structure%r0(3,uc%na))
        gs%structure%r0=uc%r

        gs%structure%ninternal=3*uc%na
        allocate(gs%structure%wV(3*uc%na))
        allocate(gs%structure%wM(3,uc%na))
        gs%structure%wV=0.0_flyt
        gs%structure%wM=0.0_flyt

        allocate(gridvals(gs%nsim,gs%structure%ninternal))
        gridvals=0.0_flyt

        ! Create interpolation:
        do ii=1,gs%nsim
            do i=1,map%n_atom_uc
                v0=gs%ref(ii)%unitcell_positions(:,i)-gs%structure%r0(:,i)
                v0=lo_clean_fractional_coordinates(v0+0.5_flyt)-0.5_flyt
                v0=lo_chop(v0,1E-10_flyt)
                do j=1,3
                    gridvals( ii, (i-1)*3+j ) = v0(j)
                enddo
            enddo
        enddo
        ! Set the orders, and remove any T-dependence.
        maxorder=3
        if ( gs%info%dim_temperature .gt. 0 ) maxorder(gs%info%dim_temperature)=0
        call gs%structure%ipint%generate( gs%grid_coordinates, gridvals, 3, 0.1_flyt, maxorder )
    end block bfinternal

    ! ! first fix the internal degrees of freedom
    ! internal: block
    !     real(flyt), dimension(:,:), allocatable :: wA,wB
    !     real(flyt), dimension(3) :: v0
    !     real(flyt) :: wm3x1(3,1),wm3x2(3,2),wm3x3(3,3)
    !     integer, dimension(gs%ndim) :: maxorder
    !     integer, dimension(3) :: thetaind
    !     integer :: ii,i,j,l,a1,unsh,unop,ntheta
    !
    !     ! store atomic numbers
    !     lo_allocate(gs%structure%atomic_numbers(uc%na))
    !     gs%structure%atomic_numbers=uc%atomic_number
    !     ! reference lattice vectors
    !     gs%structure%lv0=uc%latticevectors
    !     gs%structure%ilv0=uc%inv_latticevectors
    !     ! store reference positions
    !     allocate(gs%structure%r0(3,uc%na))
    !     gs%structure%r0=uc%r
    !     ! number of internal degrees of freedom
    !     gs%structure%ninternal=map%xuc%nx_fc_singlet
    !     ! work arrays for evaluation
    !     allocate(gs%structure%wV(3*uc%na))
    !     allocate(gs%structure%wM(3,uc%na))
    !     gs%structure%wV=0.0_flyt
    !     gs%structure%wM=0.0_flyt
    !     ! space for the values
    !     lo_allocate(gridvals(gs%nsim,gs%structure%ninternal))
    !     gridvals=0.0_flyt
    !     ! number of constraints on the internal degrees
    !     !gs%structure%nconstr=map%constraints%neq1
    !
    !     ! if there are some internal degrees of freedom, create the interpolation
    !     if ( gs%structure%ninternal .gt. 0 ) then
    !         ! coefficient matrix and work vector
    !         allocate(gs%structure%coeffM(3*uc%na,gs%structure%ninternal))
    !         allocate(gs%structure%wU(gs%structure%ninternal))
    !         gs%structure%coeffM=0.0_flyt
    !         do a1=1,map%n_atom_uc
    !             unsh=map%xuc%fc_singlet(a1)%irreducible_shell
    !             unop=map%xuc%fc_singlet(a1)%operation_from_shell
    !             ntheta=map%fc_singlet_shell(unsh)%nx
    !             if ( ntheta .lt. 1 ) cycle
    !             thetaind(1:ntheta)=map%fc_singlet_shell(unsh)%ind_global
    !             !thetaind(1:ntheta)=map%singletshell(unsh)%thetaind
    !             select case(ntheta)
    !             case(1)
    !                 !wm3x1=matmul( map%singletop(unop)%m3,map%singletshell(unsh)%redcoeffM )
    !                 wm3x1=matmul( map%op_singlet(unop)%m3,map%fc_singlet_shell(unsh)%coeff )
    !                 gs%structure%coeffM( (a1-1)*3+1:a1*3,thetaind(1:1) )=gs%structure%coeffM( (a1-1)*3+1:a1*3,thetaind(1:1) )+wm3x1
    !             case(2)
    !                 !wm3x2=matmul( map%singletop(unop)%m3,map%singletshell(unsh)%redcoeffM )
    !                 wm3x2=matmul( map%op_singlet(unop)%m3,map%fc_singlet_shell(unsh)%coeff )
    !                 gs%structure%coeffM( (a1-1)*3+1:a1*3,thetaind(1:2) )=gs%structure%coeffM( (a1-1)*3+1:a1*3,thetaind(1:2) )+wm3x2
    !             case(3)
    !                 !wm3x3=matmul( map%singletop(unop)%m3,map%singletshell(unsh)%redcoeffM )
    !                 wm3x3=matmul( map%op_singlet(unop)%m3,map%fc_singlet_shell(unsh)%coeff )
    !                 gs%structure%coeffM( (a1-1)*3+1:a1*3,thetaind(1:3) )=gs%structure%coeffM( (a1-1)*3+1:a1*3,thetaind(1:3) )+wm3x3
    !             end select
    !         enddo
    !         gs%structure%coeffM=lo_chop(gs%structure%coeffM,lo_sqtol)
    !
    !         ! some work arrays
    !         allocate(wA(3*uc%na,gs%structure%ninternal))
    !         allocate(wB(3*uc%na,1))
    !
    !         do ii=1,gs%nsim
    !             l=0
    !             do i=1,map%n_atom_uc
    !                 v0=gs%ref(ii)%unitcell_positions(:,i)-gs%structure%r0(:,i)
    !                 v0=lo_clean_fractional_coordinates(v0+0.5_flyt)-0.5_flyt
    !                 v0=matmul(uc%latticevectors,v0)
    !                 do j=1,3
    !                     l=l+1
    !                     wB(l,1)=v0(j)
    !                 enddo
    !             enddo
    !             wA=gs%structure%coeffM
    !             call lo_dgels(wA,wB)
    !             ! store the irreducible representation
    !             gridvals(ii,:)=wB(1:gs%structure%ninternal,1)
    !         enddo
    !         ! Set the orders, and remove any T-dependence.
    !         maxorder=3
    !         if ( gs%info%dim_temperature .gt. 0 ) maxorder(gs%info%dim_temperature)=0
    !         call gs%structure%ipint%generate( gs%grid_coordinates, gridvals, 3, 0.1_flyt, maxorder )
    !     endif
    ! end block internal

    ! figure out a neat irreducible representation of the shape of the cell
    lattvecs: block
        real(flyt), dimension(:,:), allocatable :: gridvals
        real(flyt), dimension(3,3) :: m0,m1,m2
        integer, dimension(gs%ndim) :: maxorder
        integer :: ii

        m0=uc%latticevectors
        lo_allocate(gridvals(gs%nsim,9))

        do ii=1,gs%nsim
            m1=gs%ref(ii)%unitcell_latticevectors
            m2=matmul(m1,uc%inv_latticevectors)
            m2=lo_chop(m2/abs(lo_determ(m2)**(1.0_flyt/3.0_flyt)),1E-12_flyt)
            gridvals(ii,:)=lo_flattentensor(m2)
        enddo
        maxorder=3
        call gs%structure%iplv%generate( gs%grid_coordinates, gridvals, 2, 0.1_flyt,maxorder )
    end block lattvecs

end subroutine

!> return a unitcell from depvar
subroutine interpolate_structure_from_depvar(st,depvar,uc,dim_volume)
    !> the raw structure data
    class(lo_gridsim_structure), intent(inout) :: st
    !> where to interpolate it
    real(flyt), dimension(:), intent(in) :: depvar
    !> resulting crystal structure
    type(lo_crystalstructure), intent(out) :: uc
    !> dimension that is volume
    integer, intent(in) :: dim_volume

    ! Actual things I need
    real(flyt) :: volume
    real(flyt), dimension(3,3) :: latticevectors,m0
    real(flyt), dimension(9) :: v9
    real(flyt), dimension(3) :: v0
    ! Helper things
    integer :: i,l,na

    na=size(st%atomic_numbers,1)
    ! first get the lattice vectors
    do i=1,9
        v9(i)=st%iplv%eval( i,depvar )
    enddo
    m0=lo_unflatten_2tensor(v9)
    latticevectors=matmul(m0,st%lv0)
    if ( dim_volume .gt. 0 ) then
        volume=depvar(dim_volume)*na
        latticevectors=latticevectors*( volume/abs(lo_determ(latticevectors)) )**(1.0_flyt/3.0_flyt)
        latticevectors=lo_chop(latticevectors,lo_sqtol)
    else
        !write(*,*) 'FIXME INTERPOLATE VOLUME'
        !stop
    endif

    ! and the internal positions
    if ( st%ninternal .gt. 0 ) then
        do i=1,st%ninternal
            st%wV(i)=st%ipint%eval( i,depvar )
        enddo
        !st%wV=matmul(st%coeffM,st%wU)
        l=0
        do i=1,na
            v0=st%wV( (i-1)*3+1:i*3 )
            !v0=matmul(st%ilv0,v0)
            st%wM(:,i)=st%r0(:,i)+v0
        enddo
    else
        st%wM=st%r0
    endif
    ! Now construct a structure from this!
    call uc%generate(latticevectors,st%wM,st%atomic_numbers,2)
end subroutine
