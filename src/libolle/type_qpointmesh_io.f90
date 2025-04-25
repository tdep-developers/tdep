submodule (type_qpointmesh) type_qpointmesh_io
use hdf5_wrappers, only: lo_hdf5_helper
implicit none
contains

!> dump the q-point mesh to file
module subroutine write_to_file(qp,p,filename,mem,verbosity,input_id)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> filename
    character(len=*), intent(in) :: filename
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> write to a specific unit
    integer(HID_T), intent(in), optional :: input_id

    type(lo_hdf5_helper) :: h5
    real(r8) :: timer
    real(r8), dimension(:,:), allocatable :: r2d
    real(r8), dimension(:), allocatable :: r1d
    integer, dimension(:,:), allocatable :: i2d
    integer, dimension(:), allocatable :: i1d
    integer :: i,j,maxngp

    if ( verbosity .gt. 0 ) then
        timer=walltime()
        write(lo_iou,*) ' '
        write(lo_iou,*) 'Writing q-mesh to file'
    endif

    ! Open the file
    if ( present(input_id) ) then
        h5%file_id=input_id
    else
        call h5%init(__FILE__,__LINE__)
        call h5%open_file('write',trim(filename))
    endif

    ! some metadata first:
    select type(qp)
        type is(lo_monkhorst_pack_mesh)
            i=1
        type is(lo_fft_mesh)
            i=2
        type is(lo_wedge_mesh)
            i=3
    end select
    call h5%store_attribute(i,h5%file_id,'mesh_type')

    ! dimensions of grid, number of points, number of tetrahedrons
    call h5%store_attribute(qp%n_full_point,h5%file_id,'number_of_qpoints')
    call h5%store_attribute(qp%n_irr_point,h5%file_id,'number_of_irreducible_qpoints')
    call h5%store_attribute(qp%n_full_tet,h5%file_id,'number_of_tetrahedrons')
    call h5%store_attribute(qp%n_irr_tet,h5%file_id,'number_of_irreducible_tetrahedrons')
    call h5%store_attribute(qp%timereversal,h5%file_id,'time_reversal_symmetry')
    select type(qp)
    type is(lo_fft_mesh)
        call h5%store_attribute(qp%griddensity(1),h5%file_id,'number_of_gridpoints_dimension_1')
        call h5%store_attribute(qp%griddensity(2),h5%file_id,'number_of_gridpoints_dimension_2')
        call h5%store_attribute(qp%griddensity(3),h5%file_id,'number_of_gridpoints_dimension_3')
    type is(lo_monkhorst_pack_mesh)
        call h5%store_attribute(qp%griddensity(1),h5%file_id,'number_of_gridpoints_dimension_1')
        call h5%store_attribute(qp%griddensity(2),h5%file_id,'number_of_gridpoints_dimension_2')
        call h5%store_attribute(qp%griddensity(3),h5%file_id,'number_of_gridpoints_dimension_3')
    end select

    ! Will store the reciprocal basis, does not hurt
    call h5%store_data(p%reciprocal_latticevectors/lo_bohr_to_A,h5%file_id,'reciprocal_latticevectors',enhet='1/A')
    ! Also store the BZ and that stuff, good for plotting
    call h5%open_group('write','brillouin_zone')
    call p%write_bz_to_hdf5('null',h5%group_id)
    call h5%close_group()

    if ( verbosity .gt. 0 ) write(lo_iou,*) '... stored metadata'

    ! Get the max number of gridpoints the irreducible can transform to
    maxngp=0
    do i=1,qp%n_irr_point
        maxngp=max(maxngp,qp%ip(i)%n_full_point)
    enddo

    ! Store the points
    call mem%allocate(r2d,[5,qp%n_irr_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(i1d,qp%n_irr_point,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(i2d,[2*maxngp+1,qp%n_irr_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    r2d=0.0_r8
    i1d=0
    i2d=0
    do i=1,qp%n_irr_point
        r2d(1:3,i)=matmul(p%inv_reciprocal_latticevectors,qp%ip(i)%r)
        r2d(4,i)=qp%ip(i)%integration_weight
        r2d(5,i)=qp%ip(i)%radius
        i1d(i)=qp%ip(i)%full_index
        j=qp%ip(i)%n_full_point
        i2d(1,i)=j
        i2d(2:j+1,i)=qp%ip(i)%operation_full_point
        i2d(maxngp+2:maxngp+1+j,i)=qp%ip(i)%index_full_point
    enddo
    call h5%store_data(r2d,h5%file_id,'irreducible_points')
    call h5%store_data(i1d,h5%file_id,'irreducible_points_gridind')
    call h5%store_data(i2d,h5%file_id,'irreducible_points_operations')
    call mem%deallocate(r2d,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(i1d,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(i2d,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    call mem%allocate(r2d,[5,qp%n_full_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(i2d,[2,qp%n_full_point],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    r2d=0.0_r8
    i2d=0
    do i=1,qp%n_full_point
        r2d(1:3,i)=matmul(p%inv_reciprocal_latticevectors,qp%ap(i)%r)
        r2d(4,i)=qp%ap(i)%integration_weight
        r2d(5,i)=qp%ap(i)%radius
        i2d(1,i)=qp%ap(i)%irreducible_index
        i2d(2,i)=qp%ap(i)%operation_from_irreducible
    enddo
    call h5%store_data(r2d,h5%file_id,'points')
    call h5%store_data(i2d,h5%file_id,'points_irrind_operation')
    call mem%deallocate(r2d,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(i2d,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    call mem%allocate(r1d,qp%n_irr_tet,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(i2d,[8,qp%n_irr_tet],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    r1d=0.0_r8
    i2d=0
    do i=1,qp%n_irr_tet
        i2d(1:4,i)=qp%it(i)%full_index
        i2d(5:8,i)=qp%it(i)%irreducible_index
        r1d(i)=qp%it(i)%integration_weight
    enddo
    call h5%store_data(i2d,h5%file_id,'irreducible_tetrahedrons')
    call h5%store_data(r1d,h5%file_id,'irreducible_tetrahedron_weights')
    call mem%deallocate(r1d,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(i2d,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    call mem%allocate(r1d,qp%n_full_tet,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(i2d,[8,qp%n_full_tet],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    r1d=0.0_r8
    i2d=0
    do i=1,qp%n_full_tet
        i2d(1:4,i)=qp%at(i)%full_index
        i2d(5:8,i)=qp%at(i)%irreducible_index
        r1d(i)=qp%at(i)%integration_weight
    enddo
    call h5%store_data(i2d,h5%file_id,'tetrahedrons')
    call h5%store_data(r1d,h5%file_id,'tetrahedron_weights')
    call mem%deallocate(r1d,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(i2d,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

    if ( present(input_id) ) then
        ! Do nothing
    else
        call h5%close_file()
        call h5%destroy()
    endif

    if ( verbosity .gt. 0 ) write(lo_iou,*) '... stored points and tetrahedrons'
    if ( verbosity .gt. 0 ) write(lo_iou,*) 'Finished writing mesh to file (',tochar(walltime()-timer),'s)'
end subroutine

!> read q-mesh from file, and initiate the type.
module subroutine lo_read_qmesh_from_file(qp,p,filename,mem,verbosity,input_id)
    !> the q-point mesh
    class(lo_qpoint_mesh), intent(out), allocatable :: qp
    !> the crystal structure
    type(lo_crystalstructure), intent(inout) :: p
    !> the filename
    character(len=*), intent(in) :: filename
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> how much to talk
    integer, intent(in) :: verbosity
    !> write to a specific unit
    integer(HID_T), intent(in), optional :: input_id

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:,:), allocatable :: r2d
    real(r8), dimension(:), allocatable :: r1d
    real(r8), dimension(3,4) :: nodes
    real(r8), dimension(3) :: v0,v1,v2
    real(r8) :: sqrc,f0,f1,t0
    integer, dimension(:,:), allocatable :: i2d
    integer, dimension(:), allocatable :: i1d
    integer :: i,j,k,l,maxngp

    ! start timer
    t0=walltime()
    call mem%tick()

    if ( present(input_id) ) then
        h5%file_id=input_id
    else
        call h5%init(__FILE__,__LINE__)
        call h5%open_file('read',trim(filename))
    endif

    ! Decide on meshtype
    call h5%read_attribute(i,h5%file_id,'mesh_type')
    select case(i)
        case(1)
            allocate(lo_monkhorst_pack_mesh::qp)
        case(2)
            allocate(lo_fft_mesh::qp)
        case(3)
            allocate(lo_wedge_mesh::qp)
        case default
            call lo_stop_gracefully(['ERROR: strange kind of mesh in "infile.qgrid"'],lo_exitcode_io,__FILE__,__LINE__)
    end select

    ! Now the meshtype is not abstract. Figure out the time-reversal stuff
    call h5%read_attribute(qp%timereversal,h5%file_id,'time_reversal_symmetry')
    if ( p%info%havespacegroup .eqv. .false. ) then
        call p%classify('wedge',timereversal=qp%timereversal)
    else
        if ( p%sym%timereversal .neqv. qp%timereversal ) then
            call lo_stop_gracefully(['Conflicting orders regarding time reversal symmetry'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
    endif

    if ( verbosity .gt. 0 ) then
        write(lo_iou,*) ''
        select type(qp)
        type is(lo_monkhorst_pack_mesh)
            write(lo_iou,*) 'Reading Monkhorst-Pack q-mesh from file'
        type is(lo_fft_mesh)
            write(lo_iou,*) 'Reading FFT q-mesh from file'
        type is(lo_wedge_mesh)
            write(lo_iou,*) 'Reading wedge q-mesh from file'
        end select
    endif

    ! Read size of grid
    call h5%read_attribute(qp%n_full_point,h5%file_id,'number_of_qpoints')
    call h5%read_attribute(qp%n_irr_point,h5%file_id,'number_of_irreducible_qpoints')
    call h5%read_attribute(qp%n_full_tet,h5%file_id,'number_of_tetrahedrons')
    call h5%read_attribute(qp%n_irr_tet,h5%file_id,'number_of_irreducible_tetrahedrons')
    allocate(qp%ap(qp%n_full_point))
    allocate(qp%ip(qp%n_irr_point))
    allocate(qp%it(qp%n_irr_tet))
    allocate(qp%at(qp%n_full_tet))
    select type(qp)
        type is(lo_monkhorst_pack_mesh)
            call h5%read_attribute(qp%griddensity(1),h5%file_id,'number_of_gridpoints_dimension_1')
            call h5%read_attribute(qp%griddensity(2),h5%file_id,'number_of_gridpoints_dimension_2')
            call h5%read_attribute(qp%griddensity(3),h5%file_id,'number_of_gridpoints_dimension_3')
        type is(lo_fft_mesh)
            call h5%read_attribute(qp%griddensity(1),h5%file_id,'number_of_gridpoints_dimension_1')
            call h5%read_attribute(qp%griddensity(2),h5%file_id,'number_of_gridpoints_dimension_2')
            call h5%read_attribute(qp%griddensity(3),h5%file_id,'number_of_gridpoints_dimension_3')
    end select

    if ( verbosity .gt. 0 ) then
        write(lo_iou,*) '       number of irreducible points: ',tochar(qp%n_irr_point)
        write(lo_iou,*) '             total number of points: ',tochar(qp%n_full_point)
        write(lo_iou,*) ' irreducible number of tetrahedrons: ',tochar(qp%n_irr_tet)
        write(lo_iou,*) '       total number of tetrahedrons: ',tochar(qp%n_full_tet)
    endif

    ! Get the irreducible points
    sqrc=p%bz%rmin**2
    call h5%read_data(r2d,h5%file_id,'irreducible_points')
    call h5%read_data(i1d,h5%file_id,'irreducible_points_gridind')
    call h5%read_data(i2d,h5%file_id,'irreducible_points_operations')
    maxngp=maxval(i2d(1,:))
    do i=1,qp%n_irr_point
        v0=r2d(1:3,i)
        qp%ip(i)%integration_weight=r2d(4,i)
        qp%ip(i)%radius=r2d(5,i)
        qp%ip(i)%full_index=i1d(i)
        select type(qp)
        type is(lo_monkhorst_pack_mesh)
            if ( qp%is_point_on_grid(v0) ) then
                v1=matmul(p%reciprocal_latticevectors,v0)
                if ( lo_sqnorm(v1) .lt. sqrc ) then
                    v2=0.0_r8
                else
                    v2=p%bz%gshift(v1+lo_degenvector)
                endif
                qp%ip(i)%r=lo_chop(v1-v2,1E-13_r8)
            else
                call lo_stop_gracefully(['bad irreducible qpoint'],lo_exitcode_io,__FILE__,__LINE__)
            endif
        type is(lo_fft_mesh)
            if ( qp%is_point_on_grid(v0) ) then
                v1=matmul(p%reciprocal_latticevectors,v0)
                if ( lo_sqnorm(v1) .lt. sqrc ) then
                    v2=0.0_r8
                else
                    v2=p%bz%gshift(v1+lo_degenvector)
                endif
                qp%ip(i)%r=lo_chop(v1-v2,1E-13_r8)
            else
                call lo_stop_gracefully(['bad irreducible qpoint'],lo_exitcode_io,__FILE__,__LINE__)
            endif
        type is(lo_wedge_mesh)
            v1=matmul(p%reciprocal_latticevectors,v0)
            qp%ip(i)%r=lo_chop(v1,1E-13_r8)
        end select
        ! and the small group
        call lo_get_small_group_of_qpoint(qp%ip(i),p)
        ! and the outfolding things
        j=i2d(1,i)
        qp%ip(i)%n_full_point=j
        allocate(qp%ip(i)%operation_full_point( j ))
        allocate(qp%ip(i)%index_full_point( j ))
        qp%ip(i)%operation_full_point=i2d(2:j+1,i)
        qp%ip(i)%index_full_point=i2d(maxngp+2:maxngp+1+j,i)
    enddo
    deallocate(r2d)
    deallocate(i1d)
    deallocate(i2d)

    ! Get the all the points
    call h5%read_data(r2d,h5%file_id,'points')
    call h5%read_data(i2d,h5%file_id,'points_irrind_operation')
    do i=1,qp%n_full_point
        v0=r2d(1:3,i)
        qp%ap(i)%irreducible_index=i2d(1,i)
        qp%ap(i)%operation_from_irreducible=i2d(2,i)
        qp%ap(i)%integration_weight=r2d(4,i)
        qp%ap(i)%radius=r2d(5,i)
        select type(qp)
        type is(lo_monkhorst_pack_mesh)
            if ( qp%is_point_on_grid(v0) ) then
                v1=matmul(p%reciprocal_latticevectors,v0)
                if ( lo_sqnorm(v1) .lt. sqrc ) then
                    v2=0.0_r8
                else
                    v2=p%bz%gshift(v1+lo_degenvector)
                endif
                qp%ap(i)%r=lo_chop(v1-v2,1E-13_r8)
            else
                call lo_stop_gracefully(['bad qpoint'],lo_exitcode_io,__FILE__,__LINE__)
            endif
        type is(lo_fft_mesh)
            if ( qp%is_point_on_grid(v0) ) then
                v1=matmul(p%reciprocal_latticevectors,v0)
                if ( lo_sqnorm(v1) .lt. sqrc ) then
                    v2=0.0_r8
                else
                    v2=p%bz%gshift(v1+lo_degenvector)
                endif
                qp%ap(i)%r=lo_chop(v1-v2,1E-13_r8)
            else
                call lo_stop_gracefully(['bad qpoint'],lo_exitcode_io,__FILE__,__LINE__)
            endif
        type is(lo_wedge_mesh)
            v1=matmul(p%reciprocal_latticevectors,v0)
            qp%ap(i)%r=lo_chop(v1,1E-13_r8)
        end select
        ! and the small group
        call lo_get_small_group_of_qpoint(qp%ap(i),p)
    enddo
    deallocate(r2d)
    deallocate(i2d)

    if ( verbosity .gt. 0 ) write(lo_iou,*) '... read points'

    ! Sanity check the points and the symmetryoperations so that everything is consistent
    do i=1,qp%n_full_point
        j=qp%ap(i)%irreducible_index
        k=qp%ap(i)%operation_from_irreducible
        if ( k .gt. 0 ) then
            v0=lo_operate_on_vector(p%sym%op(k),qp%ap(i)%r,reciprocal=.true.,inverse=.true.)-qp%ip(j)%r
        else
            v0=lo_operate_on_vector(p%sym%op(abs(k)),qp%ap(i)%r,reciprocal=.true.,inverse=.true.)+qp%ip(j)%r
        endif
        v0=matmul(p%inv_reciprocal_latticevectors,v0)
        v0=lo_clean_fractional_coordinates(v0+0.5_r8)-0.5_r8
        if ( lo_sqnorm(v0) .gt. lo_sqtol ) then
            call lo_stop_gracefully(['bad symmetryoperation in mesh'],lo_exitcode_io,__FILE__,__LINE__)
        endif
    enddo

    ! Read tetrahedrons, first the irreducible
    call h5%read_data(i2d,h5%file_id,'irreducible_tetrahedrons')
    call h5%read_data(r1d,h5%file_id,'irreducible_tetrahedron_weights')
    do i=1,qp%n_irr_tet
        qp%it(i)%full_index=i2d(1:4,i)
        qp%it(i)%irreducible_index=i2d(5:8,i)
        qp%it(i)%integration_weight=r1d(i)
    enddo
    deallocate(r1d)
    deallocate(i2d)
    if ( verbosity .gt. 0 ) write(lo_iou,*) '... read irreducible tetrahedrons'
    ! Then the full set
    call h5%read_data(i2d,h5%file_id,'tetrahedrons')
    call h5%read_data(r1d,h5%file_id,'tetrahedron_weights')
    do i=1,qp%n_full_tet
        qp%at(i)%full_index=i2d(1:4,i)
        qp%at(i)%irreducible_index=i2d(5:8,i)
        qp%at(i)%integration_weight=r1d(i)
    enddo
    deallocate(r1d)
    deallocate(i2d)
    if ( verbosity .gt. 0 ) write(lo_iou,*) '... read all tetrahedrons'

    ! Get tetrahedrons and aux stuff.
    select type(qp)
    type is(lo_wedge_mesh)
        ! Recalculate weights to be on the safe side?
        f0=0.0_r8
        do i=1,qp%n_irr_tet
            do j=1,4
                nodes(:,j)=qp%ip( qp%it(i)%irreducible_index(j) )%r
            enddo
            qp%it(i)%integration_weight=lo_unsigned_tetrahedron_volume(nodes)*p%sym%n*p%volume
            f0=f0+qp%it(i)%integration_weight
        enddo
        f1=0.0_r8
        do i=1,qp%n_full_tet
            do j=1,4
                nodes(:,j)=qp%ap( qp%at(i)%full_index(j) )%r
            enddo
            qp%at(i)%integration_weight=lo_unsigned_tetrahedron_volume(nodes)*p%volume
            f1=f1+qp%at(i)%integration_weight
        enddo
        do i=1,qp%n_irr_tet
            qp%it(i)%integration_weight=qp%it(i)%integration_weight/f0
        enddo
        do i=1,qp%n_full_tet
            qp%at(i)%integration_weight=qp%at(i)%integration_weight/f1
        enddo
        if ( verbosity .gt. 0 ) write(*,*) '... recalculated weights',f0,f1

        ! and the scaled basis for adaptive gaussians
        f0=(1.0_r8*qp%n_full_point)**(1.0_r8/3.0_r8) ! points per distance, sort of
        v0(1)=norm2(p%reciprocal_latticevectors(:,1))
        v0(2)=norm2(p%reciprocal_latticevectors(:,2))
        v0(3)=norm2(p%reciprocal_latticevectors(:,3))
        v0=v0*f0/sum(v0) ! get it per axis or somthing
        ! and to normal units
        qp%scaledrecbasis=p%reciprocal_latticevectors*lo_twopi
        qp%scaledrecbasis(:,1)=qp%scaledrecbasis(:,1)/v0(1)
        qp%scaledrecbasis(:,2)=qp%scaledrecbasis(:,2)/v0(2)
        qp%scaledrecbasis(:,3)=qp%scaledrecbasis(:,3)/v0(3)
    type is(lo_monkhorst_pack_mesh)
        allocate(qp%ind2gridind(3,qp%n_full_point))
        allocate(qp%gridind2ind(qp%griddensity(1),qp%griddensity(2),qp%griddensity(3)))
        l=0
        do i=1,qp%griddensity(1)
        do j=1,qp%griddensity(2)
        do k=1,qp%griddensity(3)
            l=l+1
            qp%ind2gridind(:,l)=[i,j,k]
            qp%gridind2ind(i,j,k)=l
        enddo
        enddo
        enddo
        ! and the scaled basis for adaptive gaussians
        qp%scaledrecbasis=p%reciprocal_latticevectors*lo_twopi
        qp%scaledrecbasis(:,1)=qp%scaledrecbasis(:,1)/qp%griddensity(1)
        qp%scaledrecbasis(:,2)=qp%scaledrecbasis(:,2)/qp%griddensity(2)
        qp%scaledrecbasis(:,3)=qp%scaledrecbasis(:,3)/qp%griddensity(3)
    type is(lo_fft_mesh)
        allocate(qp%ind2gridind(3,qp%n_full_point))
        allocate(qp%gridind2ind(qp%griddensity(1),qp%griddensity(2),qp%griddensity(3)))
        l=0
        do i=1,qp%griddensity(1)
        do j=1,qp%griddensity(2)
        do k=1,qp%griddensity(3)
            l=l+1
            qp%ind2gridind(:,l)=[i,j,k]
            qp%gridind2ind(i,j,k)=l
        enddo
        enddo
        enddo
        ! and the scaled basis for adaptive gaussians
        qp%scaledrecbasis=p%reciprocal_latticevectors*lo_twopi
        qp%scaledrecbasis(:,1)=qp%scaledrecbasis(:,1)/qp%griddensity(1)
        qp%scaledrecbasis(:,2)=qp%scaledrecbasis(:,2)/qp%griddensity(2)
        qp%scaledrecbasis(:,3)=qp%scaledrecbasis(:,3)/qp%griddensity(3)
    end select
    ! And everything should be done!
    if ( present(input_id) ) then
        ! Do nothing
    else
        call h5%close_file()
        call h5%destroy()
    endif
    ! Check that I did not waste any memory.
    call mem%tock(__FILE__,__LINE__)
    if ( verbosity .gt. 0 ) write(*,*) 'Finished reading q-mesh from file (',tochar(walltime()-t0),'s)'
end subroutine

!> dump minimal about the mesh to an hdf5 handle
module subroutine write_metadata_to_hdf5(qp,filename,input_id)
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> filename
    character(len=*), intent(in), optional :: filename
    !> in case I want to write it as a part of another file
    integer(HID_T), intent(in), optional :: input_id

    type(lo_hdf5_helper) :: h5

    ! Write in some other file, or create a new one
    if ( present(input_id) ) then
        ! now I assume hdf is open, and a file is open
        !@todo insert sanity check that this is really the case
        h5%file_id=input_id
    elseif ( present(filename) ) then
        ! open a new file
        call h5%init(__FILE__,__LINE__)
        call h5%open_file('write',trim(filename))
    else
        call lo_stop_gracefully(['Provide filename or input id, but only ony of them.'],lo_exitcode_param,__FILE__,__LINE__)
    endif

    ! Will probably depend a little
    select type(qp)
    type is(lo_monkhorst_pack_mesh)
        call h5%store_attribute('Monkhorst-Pack',h5%file_id,'q_mesh_type')
        call h5%store_data(qp%griddensity,h5%file_id,'grid_dimensions')
    type is(lo_fft_mesh)
        call h5%store_attribute('FFT mesh',h5%file_id,'q_mesh_type')
        call h5%store_data(qp%griddensity,h5%file_id,'grid_dimensions')
    type is(lo_wedge_mesh)
        call h5%store_attribute('Wedge mesh',h5%file_id,'q_mesh_type')
    end select
    ! Store size of mesh?
    call h5%store_attribute(qp%n_full_point,h5%file_id,'number_of_qpoints')
    call h5%store_attribute(qp%n_irr_point,h5%file_id,'number_of_irreducible_qpoints')

    if ( present(input_id) .eqv. .false. ) then
        call h5%close_file()
        call h5%destroy()
    endif

    ! Add more as needed. Symmetry operations? Entire mesh, with tetrahedrons?
end subroutine

end submodule
