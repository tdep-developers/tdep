
!> Get an interpolation for the renormalized baseline energy
subroutine create_energy_interpolation(gs,temperature_scale,mw,verbosity)
    !> grid of simulations
    class(lo_gridsim), intent(inout) :: gs
    !> temperature scale
    real(flyt), intent(in) :: temperature_scale
    !> MPI helper
    type(lo_mpi_helper), intent(in) :: mw
    !> talk?
    integer, intent(in) :: verbosity

    ! baseline energies, can really help with the interpolation
    real(flyt) :: tt0

    if ( verbosity .gt. 0 ) then
        tt0=walltime()
        write(*,*) ''
        write(*,*) 'INTERPOLATING INTERNAL ENERGY'
    endif

    ! collect energies
    collect: block
        real(flyt), dimension(:), allocatable :: rbuf
        real(flyt) :: f0,volume,eta,baseline
        integer, dimension(:), allocatable :: ibuf
        integer :: i,j,t,ii

        lo_allocate(rbuf(gs%nsim))
        lo_allocate(ibuf(gs%nsim))
        rbuf=0.0_flyt
        ibuf=0

        do t=1,gs%raw%nconf
            ii=gs%raw%gridind(t)
            ! Calculate U0, sort of
            f0=gs%raw%e(t)-gs%raw%e_polar(t)-gs%raw%e_fc_pair(t)-gs%raw%e_fc_triplet(t)-gs%raw%e_fc_quartet(t)
            ! subtract static energy:
            select type(eos=>gs%eos)
            class is(lo_eos_1d)
                volume=gs%grid_coordinates( gs%info%dim_volume, ii )
                baseline=eos%energy_from_volume( volume )
            class is(lo_eos_2d)
                volume=gs%grid_coordinates( gs%info%dim_volume, ii )
                eta=gs%grid_coordinates( gs%info%dim_eta, ii )
                baseline=eos%energy_from_volume_eta( volume,eta )
            class default
                baseline=0.0_flyt
            end select
            f0=f0/gs%na_ss-baseline
            ibuf(ii)=ibuf(ii)+1
            rbuf(ii)=rbuf(ii)+f0
        enddo
        ! Add it up
        call mpi_allreduce(MPI_IN_PLACE,rbuf,gs%nsim,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
        call mpi_allreduce(MPI_IN_PLACE,ibuf,gs%nsim,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)

        ! Store
        lo_allocate(gs%energy%energy(gs%nsim))
        lo_allocate(gs%energy%xi(gs%ndim,gs%nsim))
        do i=1,gs%nsim
            gs%energy%energy(i)=rbuf(i)/real(ibuf(i),flyt)
        enddo
        lo_allocate(gs%energy%xi_nontr(gs%ndim,gs%nsim))
        gs%energy%xi_nontr=gs%grid_coordinates

        ! Apply the coordinate transformation. Repeated here, since I'm not sure I want the same
        ! transformation for U0 as for the forceconstants.
        gs%energy%temperature_scale=temperature_scale
        gs%energy%dim_temperature=gs%info%dim_temperature
        gs%energy%xi=gs%grid_coordinates
        if ( gs%energy%temperature_scale .gt. 0.0_flyt .and. gs%energy%dim_temperature .gt. 0 ) then
            do i=1,gs%nsim
                j=gs%energy%dim_temperature
                gs%energy%xi( j,i )=tempscaler( gs%grid_coordinates(j,i), gs%energy%temperature_scale )
            enddo
        endif

        ! Scale all the input variables
        lo_allocate(gs%energy%x0(gs%poly%ndim))
        lo_allocate(gs%energy%ixs(gs%poly%ndim))
        lo_allocate(gs%energy%x0_nontr(gs%poly%ndim))
        lo_allocate(gs%energy%ixs_nontr(gs%poly%ndim))
        do i=1,gs%poly%ndim
            gs%energy%x0(i)=minval(gs%energy%xi(i,:))
            gs%energy%ixs(i)=1.0_flyt/(maxval(gs%energy%xi(i,:))-minval(gs%energy%xi(i,:)))
            gs%energy%x0_nontr(i)=minval(gs%energy%xi_nontr(i,:))
            gs%energy%ixs_nontr(i)=1.0_flyt/(maxval(gs%energy%xi_nontr(i,:))-minval(gs%energy%xi_nontr(i,:)))
        enddo
        do i=1,size(gs%energy%xi,2)
            gs%energy%xi(:,i)=(gs%energy%xi(:,i)-gs%energy%x0(:))*gs%energy%ixs(:)
            gs%energy%xi_nontr(:,i)=(gs%energy%xi_nontr(:,i)-gs%energy%x0_nontr(:))*gs%energy%ixs_nontr(:)
        enddo
        f0=lo_huge
        do i=1,size(gs%energy%xi,2)
        do j=i+1,size(gs%energy%xi,2)
            !f0=min(f0,norm2(gs%energy%xi(:,i)-gs%energy%xi(:,j)))
            f0=min(f0,norm2(gs%energy%xi_nontr(:,i)-gs%energy%xi_nontr(:,j)))
        enddo
        enddo
        gs%energy%xdist=f0
        if ( verbosity .gt. 0 ) write(*,*) '... distributed energies'
    end block collect

    ! Initialize the interpolating polynomial
    poly: block
        character(len=2), dimension(:), allocatable :: cnames
        integer :: i,ndim,order

        order=2
        ndim=size(gs%energy%xi,1)
        lo_allocate(cnames(ndim))
        do i=1,ndim
            cnames(i)='x'//tochar(i)
        enddo
        call gs%energy%pl%init(order,ndim,gs%energy%xi,cnames)
        if ( verbosity .gt. 0 ) write(*,*) '... created energy interpolation (',tochar(walltime()-tt0),'s)'
    end block poly
end subroutine

!> return an energy from depvar
subroutine interpolate_energy_from_depvar(en,gridcoord,energy)
    !> the raw structure data
    class(lo_gridsim_energy), intent(in) :: en
    !> where to interpolate it
    real(flyt), dimension(:), intent(in) :: gridcoord
    !> resulting energy
    real(flyt), intent(out) :: energy

    ! Helper things
    real(flyt), dimension(en%pl%ndim) :: dv,dv_nontr
    real(flyt), dimension(:,:), allocatable :: wA,wB
    real(flyt), dimension(:), allocatable :: weights
    real(flyt) :: eta,f0
    integer :: i,npts,ncoeff

    ! First temperate-scale?
    dv=gridcoord
    if ( en%temperature_scale .gt. 0.0_flyt .and. en%dim_temperature .gt. 0 ) then
        i=en%dim_temperature
        dv(i)=tempscaler( gridcoord(i),en%temperature_scale )
    endif

    ! Scale the variable
    dv=(dv-en%x0)*en%ixs
    dv_nontr=(gridcoord-en%x0_nontr)*en%ixs_nontr

    ncoeff=en%pl%ncoeff
    npts=size(en%xi,2)
    ! calculate the weights
    allocate(weights(npts))
    weights=0.0_flyt
    eta=en%xdist*2.0_flyt
    do i=1,npts
        !@todo make informed decision about the metric here
        f0=norm2(dv-en%xi(:,i))+eta
        weights(i)=1.0_flyt/f0
    enddo
    ! Space for solver
    allocate(wA(npts,ncoeff))
    allocate(wB(npts,1))
    wA=0.0_flyt
    wB=0.0_flyt
    ! Weighted coefficient matrix
    do i=1,npts
        wA(i,:)=en%pl%coeffM(i,:)*weights(i)
    enddo
    ! Interpolate the volume first.
    do i=1,npts
        wB(i,1)=en%energy(i)*weights(i)
    enddo
    call lo_dgels(wA,wB)
    energy=en%pl%eval(dv,wB(1:ncoeff,1))
    deallocate(wA)
    deallocate(wB)
    deallocate(weights)
end subroutine
