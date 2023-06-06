submodule (type_qpointmesh) type_qpointmesh_bandstructure
use lo_sorting, only: lo_return_unique
implicit none

contains

!> Interpolate a path given all the start and end-points
module subroutine interpolate_kpoint_path(bs,p,mw,verbosity)
    !> path in kspace
    class(lo_bandstructure), intent(inout) :: bs
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(3) :: v1,v2
    real(r8) :: f0,f1,length
    integer :: i,j,l,ctr

    ! The structure really needs all the symmetry stuff
    if ( p%info%havewedge .eqv. .false. ) then
        call lo_stop_gracefully(['Irreducuble wedge is needed to interpolate a path'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
    endif
    if ( p%info%pointslabelled .eqv. .false. ) then
        call lo_stop_gracefully(['The points needs to be labelled to interpolate a path'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
    endif

    ! Build all q-points
    if ( verbosity .ge. 2 ) write(lo_iou,*) '... interpolating points for path'

    if ( bs%stride .gt. 1 ) then
        ! We might have to adjust the number of points
        ! so that the stride makes sense.
        ctr=0
        do i=1,bs%n_point_per_path
            ctr=ctr+bs%stride
            if ( ctr .ge. bs%n_point_per_path-bs%stride ) exit
        enddo
        ctr=ctr+1
        bs%n_point_per_path=ctr

        ! Now we can figure out which the hard q-points are, and
        ! the rest are supposed to be interpolated at a later point.
        ctr=0
        do i=1,bs%n_point_per_path,bs%stride
            ctr=ctr+1
        enddo
        allocate(bs%stride_q_ind(ctr,bs%n_path))
        bs%stride_q_ind=-1

        do i=1,bs%n_path
            l=0
            do j=1,bs%n_point_per_path,bs%stride
                l=l+1
                bs%stride_q_ind(l,i)=(i-1)*bs%n_point_per_path + j
            enddo
        enddo
    else
        bs%stride=-1
    endif

    ! Make space for the q-points
    bs%n_point=bs%n_path*bs%n_point_per_path
    allocate(bs%q(bs%n_point))
    l=0
    do i=1,bs%n_path
        v1=bs%segment(i)%r1
        v2=bs%segment(i)%r2
        do j=1,bs%n_point_per_path
            f0=(j*1.0_r8-1.0_r8)/(real(bs%n_point_per_path,r8)-1.0_r8)
            l=l+1
            ! absolute coordinate
            bs%q(l)%r_abs=lo_chop( (1.0_r8-f0)*v1+f0*v2,1E-14_r8 )
            ! coordinate in the first BZ
            bs%q(l)%r=lo_chop( bs%q(l)%r_abs - p%bz%gshift( bs%q(l)%r_abs +lo_degenvector),1E-14_r8)
            ! which path is it on
            bs%q(l)%path=i
        enddo
    enddo
    if ( verbosity .ge. 2 ) write(lo_iou,*) '... getting small group for each q-point'

    ! Get the small group of each q-point
    do l=1,bs%n_point
        call lo_get_small_group_of_qpoint(bs%q(l),p)
    enddo

    ! Get the q-axis for plots
    if ( verbosity .ge. 2 ) write(lo_iou,*) '... building axis for plots'

    allocate(bs%q_axis(bs%n_point))
    allocate(bs%q_axis_ticks(bs%n_path+1))
    allocate(bs%q_axis_tick_labels(bs%n_path+1))
    l=0
    f0=0.0_r8
    bs%q_axis_ticks(1)=0.0_r8
    do i=1,bs%n_path
        v1=bs%segment(i)%r1
        v2=bs%segment(i)%r2
        length=norm2(v1-v2)
        ! set the tick
        bs%q_axis_ticks(i+1)=f0+length
        ! get the x-coordinate
        do j=1,bs%n_point_per_path
            l=l+1
            f1=((j-1)*1.0_r8)/(1.0_r8*(bs%n_point_per_path-1))
            bs%q_axis(l)=f0+f1*length
        enddo
        f0=f0+length
        bs%q_axis_tick_labels(i)=bs%symb_q_start(i)
    enddo

    ! And the tick labels
    if ( verbosity .ge. 2 ) write(lo_iou,*) '... setting ticks'
    bs%q_axis_tick_labels(1)=bs%symb_q_start(1)
    do i=1,bs%n_path-1
        if ( trim(bs%symb_q_end(i)) .eq. trim(bs%symb_q_start(i+1)) ) then
            bs%q_axis_tick_labels(i+1)=bs%symb_q_start(i+1)
        else
            bs%q_axis_tick_labels(i+1)=trim(bs%symb_q_end(i))//'|'//trim(bs%symb_q_start(i+1))
        endif
    enddo
    bs%q_axis_tick_labels(bs%n_path+1)=bs%symb_q_end( bs%n_path )

    ! Get proper greek Gamma
    do i=1,bs%n_path+1
        if ( trim(bs%q_axis_tick_labels(i)) .eq. 'GM' ) bs%q_axis_tick_labels(i)='Γ'
    enddo
end subroutine

!> Read a Brillouin zone path from file
module subroutine read_path_from_file(bs,p,mw,verbosity)
    !> the bandstructure
    class(lo_bandstructure), intent(inout) :: bs
    !> the crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(3) :: v0,v1
    integer :: i,u
    character(len=1000) :: dum

    ! Doublecheck that the structure has what is needed
    if ( p%info%havewedge .eqv. .false. ) then
        call lo_stop_gracefully(['Irreducuble wedge is needed to read a path from file'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
    endif
    if ( p%info%pointslabelled .eqv. .false. ) then
        call lo_stop_gracefully(['The points needs to be labelled to read a path from file'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
    endif

    if ( verbosity .gt. 1 ) write(lo_iou,*) '... getting q-path from file'

    !@TODO Fix read only on one rank.
    u=open_file('in','infile.qpoints_dispersion')
        read(u,*) dum
        select case(trim(dum))
        case('CUSTOM')
            ! A path that can go wherever
            read(u,*) bs%n_point_per_path
            read(u,*) bs%n_path
            allocate(bs%symb_q_start(bs%n_path))
            allocate(bs%symb_q_end(bs%n_path))
            allocate(bs%segment(bs%n_path))
            do i=1,bs%n_path
                read(u,*) v0,v1,bs%symb_q_start(i),bs%symb_q_end(i)
                bs%segment(i)%r1=p%fractional_to_cartesian(v0,reciprocal=.true.)
                bs%segment(i)%r2=p%fractional_to_cartesian(v1,reciprocal=.true.)
            enddo
        case default
            ! just a normal specification, such as FCC or whatever
            read(u,*) bs%n_point_per_path
            read(u,*) bs%n_path
            allocate(bs%symb_q_start(bs%n_path))
            allocate(bs%symb_q_end(bs%n_path))
            allocate(bs%segment(bs%n_path))
            v0=[1000.0_r8,0.0_r8,0.0_r8]
            do i=1,bs%n_path
                read(u,*) bs%symb_q_start(i),bs%symb_q_end(i)
                bs%segment(i)%r1=p%coordinate_from_high_symmetry_point_label(bs%symb_q_start(i),previous=v0)
                v0=bs%segment(i)%r1
                bs%segment(i)%r2=p%coordinate_from_high_symmetry_point_label(bs%symb_q_end(i),previous=v0)
                v0=bs%segment(i)%r2
            enddo
        end select
    close(u)

    call interpolate_kpoint_path(bs,p,mw,verbosity)
end subroutine

!> Tabulated paths for most Bravais lattices in the brillouin zone
module subroutine standardpath(bs,p,mw,verbosity)
    !> bandstructure
    class(lo_bandstructure), intent(inout) :: bs
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> talk a lot?
    integer, intent(in) :: verbosity

    character(len=8), dimension(:,:), allocatable :: points
    character(len=8), dimension(:), allocatable :: allpts,unpts
    real(r8), dimension(3) :: v0
    integer, dimension(:), allocatable :: dum
    integer :: i,j,k,l

    ! Doublecheck that the structure has what is needed
    if ( p%info%havewedge .eqv. .false. ) then
        call lo_stop_gracefully(['Irreducuble wedge is needed for the standard path'],lo_exitcode_param,__FILE__,__LINE__)
    endif
    if ( p%info%pointslabelled .eqv. .false. ) then
        call lo_stop_gracefully(['The points needs to be labelled for the standard path'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    if ( verbosity .ge. 1 ) write(*,*) '... getting standard q-path'
    select case(trim(p%info%bravaislattice))
        case('CUB')
            allocate(points(5,2))
            points(1,:)=(/'GM','X '/)
            points(2,:)=(/'X ','M '/)
            points(3,:)=(/'M ','GM'/)
            points(4,:)=(/'GM','R '/)
            points(5,:)=(/'R ','M '/)
        case('FCC')
            allocate(points(4,2))
            points(1,:)=(/'GM','X '/)
            points(2,:)=(/'X ','U '/)
            points(3,:)=(/'K ','GM'/)
            points(4,:)=(/'GM','L '/)
        case('BCC')
            allocate(points(5,2))
            points(1,:)=(/'GM','H '/)
            points(2,:)=(/'H ','N '/)
            points(3,:)=(/'N ','GM'/)
            points(4,:)=(/'GM','P '/)
            points(5,:)=(/'P ','H '/)
        case('HEX')
            allocate(points(5,2))
            points(1,:)=(/'GM','M '/)
            points(2,:)=(/'M ','K '/)
            points(3,:)=(/'K ','GM'/)
            points(4,:)=(/'GM','A '/)
            points(5,:)=(/'A ','L '/)
        case('HEX2')
            allocate(points(5,2))
            points(1,:)=(/'GM','M '/)
            points(2,:)=(/'M ','K '/)
            points(3,:)=(/'K ','GM'/)
            points(4,:)=(/'GM','A '/)
            points(5,:)=(/'A ','L '/)
        case('RHL1')
            allocate(points(5,2))
            points(1,:)=(/'GM','Z '/)
            points(2,:)=(/'Z ','F '/)
            points(3,:)=(/'F ','GM'/)
            points(4,:)=(/'GM','X '/)
            points(5,:)=(/'X ','L '/)
        case('TRI')
            allocate(points(7,2))
            points(1,:)=(/'X ','GM'/)
            points(2,:)=(/'GM','Y '/)
            points(3,:)=(/'L ','GM'/)
            points(4,:)=(/'GM','Z '/)
            points(5,:)=(/'N ','GM'/)
            points(6,:)=(/'GM','M '/)
            points(7,:)=(/'R ','GM'/)
        case('ORC')
            allocate(points(7,2))
            points(1,:)=(/'GM','X '/)
            points(2,:)=(/'X ','S '/)
            points(3,:)=(/'S ','Y '/)
            points(4,:)=(/'Y ','GM'/)
            points(5,:)=(/'GM','Z '/)
            points(6,:)=(/'Z ','U '/)
            points(7,:)=(/'U ','R '/)
        case('ORCC1')
            allocate(points(8,2))
            points(1,:)=(/'GM','X '/)
            points(2,:)=(/'X ','S '/)
            points(3,:)=(/'S ','R '/)
            points(4,:)=(/'R ','A '/)
            points(5,:)=(/'A ','Z '/)
            points(6,:)=(/'Z ','GM'/)
            points(7,:)=(/'GM','Y '/)
            points(8,:)=(/'Y ','X1'/)
        case('ORCC2')
            allocate(points(8,2))
            points(1,:)=(/'GM','X '/)
            points(2,:)=(/'X ','S '/)
            points(3,:)=(/'S ','R '/)
            points(4,:)=(/'R ','A '/)
            points(5,:)=(/'A ','Z '/)
            points(6,:)=(/'Z ','GM'/)
            points(7,:)=(/'GM','Y '/)
            points(8,:)=(/'Y ','X1'/)
        case('TET')
            allocate(points(9,2))
            points(1,:)=(/'GM','X '/)
            points(2,:)=(/'X ','M '/)
            points(3,:)=(/'M ','GM'/)
            points(4,:)=(/'GM','Z '/)
            points(5,:)=(/'Z ','R '/)
            points(6,:)=(/'R ','A '/)
            points(7,:)=(/'A ','Z '/)
            points(8,:)=(/'X ','R '/)
            points(9,:)=(/'M ','A '/)
        case('MCL')
            allocate(points(7,2))
            points(1,:)=(/'GM','Y '/)
            points(2,:)=(/'Y ','H '/)
            points(3,:)=(/'H ','C '/)
            points(4,:)=(/'C ','GM'/)
            points(5,:)=(/'GM','X '/)
            points(6,:)=(/'X ','H2'/)
            points(7,:)=(/'H2','Y '/)
        case('MCLC1')
            allocate(points(5,2))
            points(1,:)=(/'GM','L '/)
            points(2,:)=(/'L ','F '/)
            points(3,:)=(/'F ','GM'/)
            points(4,:)=(/'GM','X '/)
            points(5,:)=(/'X ','M '/)
        case('MCLC3')
            allocate(points(5,2))
            points(1,:)=(/'GM','L '/)
            points(2,:)=(/'L ','Z '/)
            points(3,:)=(/'Z ','GM'/)
            points(4,:)=(/'GM','X '/)
            points(5,:)=(/'X ','M '/)
        case('MCLC5')
            allocate(points(5,2))
            points(1,:)=(/'GM','N '/)
            points(2,:)=(/'N ','F '/)
            points(3,:)=(/'F ','GM'/)
            points(4,:)=(/'GM','X '/)
            points(5,:)=(/'X ','M '/)
        case('ORCF2')
            allocate(points(5,2))
            points(1,:)=(/'GM','L '/)
            points(2,:)=(/'L ','H '/)
            points(3,:)=(/'H ','GM'/)
            points(4,:)=(/'GM','X '/)
            points(5,:)=(/'X ','D '/)
        case('TRI1a')
            allocate(points(5,2))
            points(1,:)=(/'GM','L '/)
            points(2,:)=(/'L ','R '/)
            points(3,:)=(/'R ','GM'/)
            points(4,:)=(/'GM','Y '/)
            points(5,:)=(/'Y ','M '/)
        case('TRI1b')
            allocate(points(5,2))
            points(1,:)=(/'GM','L '/)
            points(2,:)=(/'L ','R '/)
            points(3,:)=(/'R ','GM'/)
            points(4,:)=(/'GM','Y '/)
            points(5,:)=(/'Y ','M '/)
        case('TRI2a')
            allocate(points(5,2))
            points(1,:)=(/'GM','L '/)
            points(2,:)=(/'L ','R '/)
            points(3,:)=(/'R ','GM'/)
            points(4,:)=(/'GM','Y '/)
            points(5,:)=(/'Y ','M '/)
        case('TRI2b')
            allocate(points(5,2))
            points(1,:)=(/'GM','L '/)
            points(2,:)=(/'L ','R '/)
            points(3,:)=(/'R ','GM'/)
            points(4,:)=(/'GM','Y '/)
            points(5,:)=(/'Y ','M '/)
        case('ORCI')
            allocate(points(6,2))
            points(1,:)=['GM','S ']
            points(2,:)=['S ','Y ']
            points(3,:)=['Y ','W ']
            points(4,:)=['W ','GM']
            points(5,:)=['GM','X ']
            points(6,:)=['X ','T ']
        case('BCT1')
            allocate(points(5,2))
            points(1,:)=['GM','X ']
            points(2,:)=['X ','M ']
            points(3,:)=['M ','GM']
            points(4,:)=['GM','Z ']
            points(5,:)=['Z ','P ']
        case('BCT2')
            allocate(points(5,2))
            points(1,:)=['GM','X ']
            points(2,:)=['X ','N ']
            points(3,:)=['N ','GM']
            points(4,:)=['GM','R ']
            points(5,:)=['R ','P ']
        case default
            write(*,*) 'could not find a standard path for this Bravais lattice. Nag on me and I might fix it.'
            write(*,*) 'or if you don not want to wait, add it in "type_qpointmesh_bandstructure.f90"'
            stop
    end select

    ! To make life easier, I want to make sure that these points actually
    ! exist. If they are missing, I pick a random one instead, annying if it
    ! crashes here. Not really important, it's just some random German dude
    ! who decided what things are called, not that it has any actual importance.

    ! Flat list of all points
    allocate(allpts(size(points,1)*size(points,2)))
    l=0
    do i=1,size(points,1)
    do j=1,2
        l=l+1
        allpts(l)=points(i,j)
    enddo
    enddo
    ! Get the unique
    call lo_return_unique(allpts,unpts)

    ! Test if these exist:
    allocate(dum(size(unpts,1)))
    dum=0
    l=0
    unloop: do i=1,size(unpts,1)
        do j=1,p%bz%nhighsymmetrypoints
            if ( trim(unpts(i)) .eq. trim(p%bz%label(j)) ) then
                ! found it, ok
                cycle unloop
            endif
        enddo
        ! did not find it. Find some other point instead, that does exist
        do j=1,p%bz%nhighsymmetrypoints
        do k=1,p%bz%nhighsymmetrypoints
            if ( trim(p%bz%label(j)) .eq. 'NP'//tochar(k) ) then
            if ( k .gt. l ) then
                dum(i)=k
                l=k
                ! got it!
                cycle unloop
            endif
            endif
        enddo
        enddo
    enddo unloop

    ! Relabel points, if there are any missing
    do i=1,size(points,1)
    do j=1,2
        do l=1,size(unpts,1)
            if ( trim(unpts(l)) .eq. trim(points(i,j)) ) then
            if ( dum(l) .ne. 0 ) then
                ! Relabel this point!
                points(i,j)='NP'//tochar(dum(l))
            endif
            endif
        enddo
    enddo
    enddo

    ! Store stuff in the bs structure
    bs%n_path=size(points,1)
    allocate(bs%symb_q_start(bs%n_path))
    allocate(bs%symb_q_end(bs%n_path))
    allocate(bs%segment(bs%n_path))
    v0=[100.0_r8,2.0_r8,1.0_r8] ! to bias the starting point
    j=0
    do i=1,bs%n_path
        bs%symb_q_start(i)=points(i,1)
        bs%symb_q_end(i)=points(i,2)
        bs%segment(i)%r1=p%coordinate_from_high_symmetry_point_label(bs%symb_q_start(i),previous=v0)
        v0=bs%segment(i)%r1
        bs%segment(i)%r2=p%coordinate_from_high_symmetry_point_label(bs%symb_q_end(i),previous=v0)
        v0=bs%segment(i)%r2
    enddo

    ! get all the q-points
    call interpolate_kpoint_path(bs,p,mw,verbosity)

    if ( allocated(points) ) deallocate(points)
    if ( allocated(allpts) ) deallocate(allpts)
    if ( allocated(unpts ) ) deallocate(unpts )
    if ( allocated(dum   ) ) deallocate(dum   )
end subroutine

end submodule
