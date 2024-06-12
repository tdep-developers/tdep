#include "precompilerdefinitions"
submodule (type_qpointmesh) type_qpointmesh_gridgeneration
use mpi
implicit none
contains

!> Generate a q-point grid using fancy-pants polymorphism
module subroutine build_grid(qp,p,mw,mem,verbosity,nosymmetry)
    !> q-point grid
    class(lo_qpoint_grid), intent(inout) :: qp
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: p
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity
    !> Skip point-group symmetries?
    logical, intent(in) :: nosymmetry

    real(r8), dimension(:,:), allocatable :: points
    real(r8) :: t0,t1
    integer, dimension(:,:), allocatable :: connections
    integer :: np,nsym

    ! Start the timer
    t0=walltime()
    t1=t0

    ! Set the basic things first
    init: block
        real(r8), dimension(:,:), allocatable :: dr0,dr1
        real(r8) :: f0,f1
        integer :: i,j,isym,ctr

        ! Say what kind of mesh we are building?
        if ( verbosity .gt. 0 ) then
            select type(qp)
            type is(lo_monkhorst_pack_mesh)
                write(lo_iou,'(1X,A)') '... building '//tochar(qp%griddensity(1))//' X '//tochar(qp%griddensity(2))//' X '//tochar(qp%griddensity(3))//' Monkhorst-Pack k-point mesh using '//tochar(mw%n)//' ranks.'
            type is(lo_fft_mesh)
                write(lo_iou,'(1X,A)') '... building '//tochar(qp%griddensity(1))//' X '//tochar(qp%griddensity(2))//' X '//tochar(qp%griddensity(3))//' DFT k-point mesh using '//tochar(mw%n)//' ranks.'
            end select
            if ( nosymmetry ) then
                write(lo_iou,*) '... not enforcing symmetry'
            endif
            if ( qp%timereversal ) then
                write(lo_iou,*) '... enforcing time reversal symmetry'
            endif
        endif

        ! Decide on the number of symmetry operations.
        if ( nosymmetry ) then
            nsym=1
        else
            nsym=p%sym%n
        endif

        ! Also, now I know the total number of points.
        np=product(qp%griddensity)

        ! Decide on which operations are ok to use. This can catch some really odd corner cases.
        allocate(qp%operationok(nsym))
        qp%operationok=.false.
        call mem%allocate(dr0,[3,p%bz%nnodes],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(dr1,[3,p%bz%nnodes],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        dr0=p%bz%r
        dr1=0.0_r8
        ctr=0
        do isym=1,nsym
            ! Rotate the BZ
            do i=1,p%bz%nnodes
                dr1(:,i)=lo_operate_on_vector(p%sym%op(isym),dr0(:,i),reciprocal=.true.,fractional=.false.)
            enddo
            ! Check that it's still the same
            f0=0.0_r8
            do i=1,p%bz%nnodes
                f1=lo_huge
                do j=1,p%bz%nnodes
                    f1=min(f1,lo_sqnorm(dr0(:,i)-dr1(:,j)))
                enddo
                f0=f0+abs(f1)
            enddo
            if ( f0 .lt. lo_sqtol ) then
                ctr=ctr+1
                qp%operationok(isym)=.true.
            endif
        enddo
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) '... found ',tochar(ctr),' operations out of ',tochar(nsym),' that leave the BZ invariant.'
        endif

        call mem%deallocate(dr0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(dr1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block init

    ! Then we build all the points. Do this only on the head rank, for now.
    buildpoints: block
        real(r8), dimension(3) :: v0,v1
        real(r8) :: sqrc
        integer :: i,j,k,l

        ! Space in the structure for indexing
        allocate(qp%ind2gridind(3,np))
        allocate(qp%gridind2ind( qp%griddensity(1) , qp%griddensity(2) , qp%griddensity(3) ))
        qp%ind2gridind=0
        qp%gridind2ind=0

        ! Space for points
        call mem%allocate(points,[3,np],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        points=0.0_r8
        sqrc=p%bz%rmin**2
        l=0
        do i=1,qp%griddensity(1)
        do j=1,qp%griddensity(2)
        do k=1,qp%griddensity(3)
            l=l+1
            if ( mod(l,mw%n) .ne. mw%r ) cycle
            ! this gives points in fractional coordinates
            v0=qp%coordinate_from_index([i,j,k])
            qp%ind2gridind(:,l)=[i,j,k]
            qp%gridind2ind(i,j,k)=l
            ! Directly move it to the first BZ
            if ( lo_sqnorm(matmul(p%reciprocal_latticevectors,v0)) .gt. sqrc ) then
                ! just to make it predictable on what side of the BZ things end up I
                ! add a small number to the q-vectors when moving them to the first BZ.
                v1=matmul(p%reciprocal_latticevectors,v0)+lo_degenvector
                v1=anint(matmul(p%inv_reciprocal_latticevectors,p%bz%gshift(v1)))
                points(:,l)=v0-v1
            else
                points(:,l)=v0
            endif
        enddo
        enddo
        enddo
        ! Sync
        call mw%allreduce('sum',points)
        call mw%allreduce('sum',qp%ind2gridind)
        call mw%allreduce('sum',qp%gridind2ind)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... built raw mesh ('//tochar(t1-t0,6)//'s)'
            t0=t1
        endif
    end block buildpoints

    ! Now start with the symmetry reduction thingy. Also headrank only, for now.
    buildconnections: block
        real(r8), dimension(3) :: v0,v1,v2
        integer, dimension(3) :: gi,gif
        integer :: i,j,l,ii,jj,ctr

        ! Figure out the unique ones by check which points transform to which points under what
        ! operation. I build a list, such that
        !   op * qi = qj
        ! This section is a little more careful than it needs to be, but it is fast enough anyway
        ! and great for catching strange corner cases.
        if ( qp%timereversal ) then
            call mem%allocate(connections,[np,nsym*2],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        else
            call mem%allocate(connections,[np,nsym],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        endif
        connections=0

        if ( verbosity .gt. 0 ) call lo_progressbar_init()

        ! Count number of points on this rank, for the progressbar.
        ctr=0
        do i=1,np
            if ( mod(i,mw%n) .eq. mw%r ) ctr=ctr+1
        enddo

        l=0
        do i=1,np
            ! make it parallel
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            ! the q-point to be tested
            v0=points(:,i)
            ! check vs all symmetry operations
            syml1: do j=1,nsym
                if ( qp%operationok(j) .eqv. .false. ) cycle syml1
                ! rotate the vector.
                v1=lo_operate_on_vector(p%sym%op(j),v0,reciprocal=.true.,fractional=.true.,inverse=.false.)
                ! is it on the grid?
                if ( qp%is_point_on_grid(v1) .eqv. .false. ) cycle syml1
                ! it's on the grid, get what point it transforms to
                gi=qp%index_from_coordinate(v1)
                ii=qp%gridind2ind(gi(1),gi(2),gi(3))
                ! rotate it backwards
                v2=lo_operate_on_vector(p%sym%op(j),points(:,ii),reciprocal=.true.,fractional=.true.,inverse=.true.)
                ! check if the backwards point is on the grid
                if ( qp%is_point_on_grid(v2) .eqv. .false. ) cycle syml1
                ! this is starting to look good. But I'm not convinced yet.
                gif=qp%index_from_coordinate(v2)
                jj=qp%gridind2ind(gif(1),gif(2),gif(3))
                ! it has to return to the same point
                if ( i .ne. jj ) cycle syml1
                ! Yep, is all of this passes, we have a connection!
                connections(i,j)=ii
            enddo syml1

            if ( qp%timereversal ) then
                ! The exact same thing again, but testing includes time-reversal symmetry
                ! which just mean I test if
                !   op * qi = -qj
                ! is valid as well.
                syml2: do j=1,nsym
                    if ( qp%operationok(j) .eqv. .false. ) cycle syml2
                    v1=-lo_operate_on_vector(p%sym%op(j),v0,reciprocal=.true.,fractional=.true.,inverse=.false.)
                    ! is it on the grid?
                    if ( qp%is_point_on_grid(v1) .eqv. .false. ) cycle syml2
                    ! it's on the grid, get what point it transforms to
                    gi=qp%index_from_coordinate(v1)
                    ii=qp%gridind2ind(gi(1),gi(2),gi(3))
                    ! rotate it backwards
                    v2=-lo_operate_on_vector(p%sym%op(j),points(:,ii),reciprocal=.true.,fractional=.true.,inverse=.true.)
                    ! check if the backwards point is on the grid
                    if ( qp%is_point_on_grid(v2) .eqv. .false. ) cycle syml2
                    ! this is starting to look good. But I'm not convinced yet.
                    gif=qp%index_from_coordinate(v2)
                    jj=qp%gridind2ind(gif(1),gif(2),gif(3))
                    ! it has to return to the same point
                    if ( i .ne. jj ) cycle syml2
                    ! Yep, is all of this passes, we have a connection!
                    connections(i,j+nsym)=-ii
                enddo syml2
            endif

            if ( verbosity .gt. 0 ) then
                l=l+1
                if ( lo_trueNtimes(l,127,ctr) ) call lo_progressbar(' ... building connections',l,ctr,walltime()-t0)
            endif
        enddo
        ! Sync across ranks
        call mw%allreduce('sum',connections,filename=__FILE__,linenumber=__LINE__)
        ! Make final report
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            call lo_progressbar(' ... building connections',ctr,ctr,t1-t0)
            t0=t1
        endif



        ! ! This is an agressive test to make sure all the symmetry operations
        ! ! I identified actually hold. This has not triggered in many years,
        ! ! but cost almost nothing so I'm keeping it.
        ! ii=0
        ! do i=1,np
        !     if ( mod(i,mw%n) .ne. mw%r ) cycle
        !     v0=points(:,i)
        !     do j=1,p%sym%n
        !         k=connections(i,j)
        !         if ( k .eq. 0 ) cycle
        !         ! ok, first connection is op(j)*v(i)=v(k)
        !         v1=points(:,k)
        !         v2=lo_operate_on_vector(p%sym%op(j),v0,reciprocal=.true.,fractional=.true.,inverse=.false.)
        !         if ( compare_qpoints_with_pbc(v1,v2,lo_sqtol) .eqv. .false. ) then
        !             ii=ii+1
        !         endif
        !         ! should also be checked that v(i)=op(j)^{-1}*v(k)
        !         v3=lo_operate_on_vector(p%sym%op(j),v1,reciprocal=.true.,fractional=.true.,inverse=.true.)
        !         if ( compare_qpoints_with_pbc(v0,v3,lo_sqtol) .eqv. .false. ) then
        !             ii=ii+1
        !         endif
        !     enddo
        !     ! Skip the next test if we don't have timereversal symmetry
        !     if ( qp%timereversal .eqv. .false. ) cycle
        !     do j=1,nsym
        !         k=connections(i,p%sym%n+j)
        !         if ( k .eq. 0 ) cycle
        !         ! ok, first connection is op(j)*v(i)=-v(k)
        !         v1=-points(:,abs(k))
        !         v2=lo_operate_on_vector(p%sym%op(j),v0,reciprocal=.true.,fractional=.true.,inverse=.false.)
        !         if ( compare_qpoints_with_pbc(v1,v2,lo_sqtol) .eqv. .false. ) then
        !             ii=ii+1
        !         endif
        !         ! should also be checked that v(i)=-op(j)^{-1}*v(k)
        !         v3=lo_operate_on_vector(p%sym%op(j),v1,reciprocal=.true.,fractional=.true.,inverse=.true.)
        !         if ( compare_qpoints_with_pbc(v0,v3,lo_sqtol) .eqv. .false. ) then
        !             ii=ii+1
        !         endif
        !     enddo
        ! enddo
        ! call mw%allreduce('sum',ii)
        ! if ( ii .ne. 0 ) then
        !     call lo_stop_gracefully(['Bad connection in connectionlist. Should never happen.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        ! endif
        ! if ( verbosity .gt. 0 ) then
        !     t1=walltime()
        !     write(lo_iou,*) '... built and tested connections ('//tochar(t1-t0,6)//'s)'
        !     t0=t1
        ! endif
    end block buildconnections

    ! Use the connection table to figure out the unique points. Also serially, for now.
    findunique: block
        integer, dimension(:), allocatable :: ind,listofirr
        integer :: i,j,k,l

        ! This little thingy actually gets the list of unique q-points. Neat.
        call mem%allocate(ind,np,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ind=1
        do i=1,np
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            do j=1,size(connections,2)
                k=abs(connections(i,j))
                ! j*i=k, k>i then we remove k?
                ! some clarification, j is an operation, i and k are q-points. If k can be constructed from
                ! i, and k>i, we don't need k.
                if ( k .gt. i ) ind(k)=0
            enddo
        enddo
        call mw%allreduce('min',ind)

        ! Build the list of irreducible points
        qp%n_irr_point=sum(ind)
        call mem%allocate(listofirr,qp%n_irr_point,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        listofirr=0
        l=0
        do i=1,np
            if ( ind(i) .eq. 1 ) then
                l=l+1
                listofirr(l)=i
            endif
        enddo

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... reduced the number of points from '//tochar(np)//' to '//tochar(qp%n_irr_point)//' ('//tochar(t1-t0,6)//'s)'
            t0=t1
        endif

        ! Store the list of irreducible points for future reference.
        lo_allocate(qp%ip(qp%n_irr_point))
        do i=1,qp%n_irr_point
            ! store position in BZ
            qp%ip(i)%r=points(:,listofirr(i))
            ! the index in the large array
            qp%ip(i)%full_index=listofirr(i)
            ! initialize the weight to 0
            qp%ip(i)%integration_weight=0.0_r8
        enddo

        ! And list of all points
        qp%n_full_point=np
        lo_allocate(qp%ap(qp%n_full_point))
        do i=1,np
            qp%ap(i)%r=points(:,i)
            qp%ap(i)%irreducible_index=0
            qp%ap(i)%operation_from_irreducible=0
            qp%ap(i)%integration_weight=0.0_r8
        enddo
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... stored points ('//tochar(t1-t0,6)//'s)'
            t0=t1
        endif

        ! Try the mapping, again. Hmm. Think long and hard. Think this is the
        ! fastest and safest way of doing it. Could perhaps transpose some array
        ! but this is fast enough anyway.
        do k=1,qp%n_irr_point
            do j=1,nsym
                i=connections( listofirr(k),j )
                if ( i .eq. 0 ) cycle
                if ( qp%ap(i)%irreducible_index .ne. 0 ) cycle
                qp%ap(i)%irreducible_index=k
                qp%ap(i)%operation_from_irreducible=j
            enddo
            if ( qp%timereversal ) then
                do j=1,nsym
                    i=-connections( listofirr(k),j+nsym )
                    if ( i .eq. 0 ) cycle
                    if ( qp%ap(i)%irreducible_index .ne. 0 ) cycle
                    qp%ap(i)%irreducible_index=k
                    qp%ap(i)%operation_from_irreducible=-j
                enddo
            endif
        enddo

        ! ! See if I got it right
        ! do i=1,qp%n_full_point
        !     ! First make sure the point got classified
        !     j=qp%ap(i)%irrind
        !     k=qp%ap(i)%operation
        !     if ( j .eq. 0 ) then
        !         call lo_stop_gracefully(['Failed miserably mapping irreducible grid to full'],lo_exitcode_symmetry,__FILE__,__LINE__)
        !     endif
        !     if ( k .eq. 0 ) then
        !         call lo_stop_gracefully(['Failed miserably mapping irreducible grid to full'],lo_exitcode_symmetry,__FILE__,__LINE__)
        !     endif
        !     ! Then make sure it was classified properly!
        !     if ( k .gt. 0 ) then
        !         v1=lo_operate_on_vector(p%sym%op(k),qp%ip(j)%v,reciprocal=.true.,fractional=.true.,inverse=.false.)
        !     else
        !         v1=-lo_operate_on_vector(p%sym%op(abs(k)),qp%ip(j)%v,reciprocal=.true.,fractional=.true.,inverse=.false.)
        !     endif
        !     v0=lo_clean_fractional_coordinates(v1-qp%ap(i)%v+0.5_r8)-0.5_r8
        !     if ( lo_sqnorm(v0) .gt. lo_sqtol ) then
        !         call lo_stop_gracefully(['Failed miserably mapping irreducible grid to full'],lo_exitcode_symmetry,__FILE__,__LINE__)
        !     endif
        ! enddo

        ! A little cleanup
        call mem%deallocate(ind,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(listofirr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(points,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block findunique

    ! Now that I figured out the unique and all that, get the weights and some other aux stuff. Also serial.
    finalize: block
        real(r8), dimension(3) :: v0,v1
        real(r8) :: f0
        integer, dimension(:), allocatable :: dum1,dum2
        integer, dimension(3) :: gi
        integer :: i,j,k,l,ii,jj,kk

        ! And the weight of each of the irreducible points
        do i=1,qp%n_full_point
            j=qp%ap(i)%irreducible_index
            qp%ip(j)%integration_weight=qp%ip(j)%integration_weight+1.0_r8
        enddo

        ! Adjust the weights so that they integrate to the correct thing
        do i=1,qp%n_irr_point
            qp%ip(i)%integration_weight=qp%ip(i)%integration_weight/real(qp%n_full_point,r8)
        enddo
        do i=1,qp%n_full_point
            qp%ap(i)%integration_weight=1.0_r8/real(qp%n_full_point,r8)
        enddo

        ! Get the q-radius per point, simplest possible way
        f0=( 3.0_r8/p%volume/real(qp%n_full_point,r8)/4.0_r8/lo_pi )**(1.0_r8/3.0_r8)
        do i=1,qp%n_irr_point
            qp%ip(i)%radius=f0
        enddo
        do i=1,qp%n_full_point
            qp%ap(i)%radius=f0
        enddo

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... weights and mapping ('//tochar(t1-t0,6)//'s)'
            t0=t1
        endif

        ! Now it's about time to convert it to Cartesian coordinates
        do i=1,qp%n_irr_point
            qp%ip(i)%r=lo_chop( p%fractional_to_cartesian(qp%ip(i)%r,reciprocal=.true.),1E-14_r8)
        enddo
        do i=1,qp%n_full_point
            qp%ap(i)%r=lo_chop( p%fractional_to_cartesian(qp%ap(i)%r,reciprocal=.true.),1E-14_r8)
        enddo

        ! I need the small group of each q-vector to accurately calculate group velocities and so on.
        do i=1,qp%n_irr_point
            call lo_get_small_group_of_qpoint(qp%ip(i),p)
        enddo
        do i=1,qp%n_full_point
            call lo_get_small_group_of_qpoint(qp%ap(i),p)
        enddo
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... got small groups ('//tochar(t1-t0,6)//'s)'
            t0=t1
        endif

        ! And the backwards list to unfold the irreducible to the full grid
        call mem%allocate(dum1,size(connections,2),persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(dum2,size(connections,2),persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        dum1=0
        dum2=0
        do i=1,qp%n_irr_point
            ! Reset the counters
            ii=qp%ip(i)%full_index
            v0=matmul(p%inv_reciprocal_latticevectors,qp%ap(ii)%r)
            kk=0
            l=0
            dum1=0
            dum2=0
            ! check vs all operations
            sml1: do k=1,nsym
                if ( qp%operationok(k) .eqv. .false. ) cycle
                ! first check the normal
                v1=lo_operate_on_vector(p%sym%op(k),v0,reciprocal=.true.,fractional=.true.)
                if ( qp%is_point_on_grid(v1) ) then
                    gi=qp%index_from_coordinate(v1)
                    jj=qp%gridind2ind(gi(1),gi(2),gi(3))
                    kk=kk+1
                    ! check if it's a new point
                    do j=1,l
                        if ( dum1(j) .eq. jj ) cycle sml1
                    enddo
                    ! if we made it here, it's a new point
                    l=l+1
                    dum1(l)=jj
                    dum2(l)=k
                endif
            enddo sml1

            ! check with time-reversal symmetry, or not
            if ( qp%timereversal ) then
                sml2: do k=1,nsym
                    if ( qp%operationok(k) .eqv. .false. ) cycle
                    v1=-lo_operate_on_vector(p%sym%op(k),v0,reciprocal=.true.,fractional=.true.)
                    if ( qp%is_point_on_grid(v1) ) then
                        gi=qp%index_from_coordinate(v1)
                        jj=qp%gridind2ind(gi(1),gi(2),gi(3))
                        kk=kk+1
                        ! check if it's a new point
                        do j=1,l
                            if ( dum1(j) .eq. jj ) cycle sml2
                        enddo
                        ! if we made it here, it's a new point
                        l=l+1
                        dum1(l)=jj
                        dum2(l)=-k
                    endif
                enddo sml2
            endif
            ! sanity check
            if ( abs(qp%ip(i)%integration_weight*qp%n_full_point-l*1.0_r8) .gt. lo_tol ) then
                call lo_stop_gracefully(['Failed mapping operation backwards'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
            ! store
            qp%ip(i)%n_full_point=l
            lo_allocate(qp%ip(i)%index_full_point(l))
            lo_allocate(qp%ip(i)%operation_full_point(l))
            qp%ip(i)%index_full_point=dum1(1:l)
            qp%ip(i)%operation_full_point=dum2(1:l)
        enddo

        call mem%deallocate(dum1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(dum2,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(connections,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        ! And that's it! For now, at least.
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) '... mapped symmetry backwards ('//tochar(walltime()-t0,6)//'s)'
        endif
    end block finalize
end subroutine

!> Get the coordinate from indices for a Monkhorst-Pack mesh
module pure function mp_coordinate_from_index(qp,i) result(r)
    !> q-mesh
    class(lo_monkhorst_pack_mesh), intent(in) :: qp
    !> indices
    integer, dimension(3), intent(in) :: i
    !> resulting coordinates
    real(r8), dimension(3) :: r

    integer :: j

    do j=1,3
        r(j)=(2.0_r8*i(j)-qp%griddensity(j)*1.0_r8-1.0_r8)/(2.0_r8*qp%griddensity(j))
    enddo
end function

!> Get the indices from coordinates for a Monkhorst-Pack mesh
module pure function mp_index_from_coordinate(qp,r) result(i)
    !> q-mesh
    class(lo_monkhorst_pack_mesh), intent(in) :: qp
    !> resulting coordinates
    real(r8), dimension(3), intent(in) :: r
    !> indices
    integer, dimension(3) :: i

    real(r8), dimension(3) :: v1
    integer :: j

    v1=r
    do j=1,3
        ! get indices
        v1(j)=(v1(j)*2.0_r8*qp%griddensity(j)+qp%griddensity(j)+1.0_r8)*0.5_r8
        ! adjust for pbc
        v1(j)=lo_index_in_periodic_array(int(anint(v1(j))),qp%griddensity(j))
    enddo
    ! convert to integer
    i=int(anint(v1))
end function

!> Checks if a vector is on the specified MP-grid
module pure function mp_is_point_on_grid(qp,r) result(ongrid)
    !> the mesh
    class(lo_monkhorst_pack_mesh), intent(in) :: qp
    !> the coordinates
    real(r8), dimension(3), intent(in) :: r
    !> is it on the grid?
    logical :: ongrid

    real(r8), dimension(3) :: v1
    integer :: j

    v1=r
    do j=1,3
        ! get indices
        v1(j)=(v1(j)*2.0_r8*qp%griddensity(j)+qp%griddensity(j)+1.0_r8)*0.5_r8
    enddo
    ! check if the generating indices are integers
    if ( sum(abs(v1-anint(v1))) .lt. lo_tol ) then
        ongrid=.true.
    else
        ongrid=.false.
    endif
end function

!> Checks if a vector is on the specified FFT-grid
module pure function fft_is_point_on_grid(qp,r) result(ongrid)
    !> the mesh
    class(lo_fft_mesh), intent(in) :: qp
    !> the resulting coordinates
    real(r8), dimension(3), intent(in) :: r
    !> is it on the grid?
    logical :: ongrid

    integer :: j
    real(r8), dimension(3) :: v1

    v1=r
    do j=1,3
        ! get indices
        v1(j)=v1(j)*qp%griddensity(j)+1.0_r8
    enddo
    ! check if the generating indices are integers
    if ( sum(abs(v1-anint(v1))) .lt. lo_tol ) then
        ongrid=.true.
    else
        ongrid=.false.
    endif
end function

!> Gets the indices from the coordinate:
!> $$ i = xn+1  $$
!> where \(x\) the coordinate, \(i\) the index and \(n\) the number of
!> gridpoints in this dimension.
module pure function fft_index_from_coordinate(qp,r) result(i)
    !> the mesh
    class(lo_fft_mesh), intent(in) :: qp
    !> the coordinates, in fractional units
    real(r8), dimension(3), intent(in) :: r
    !> the indices
    integer, dimension(3) :: i

    integer :: j
    real(r8), dimension(3) :: v1

    ! Get it to 0-1 coordinates, inclusive zero not inclusive 1
    v1=lo_clean_fractional_coordinates(r)
    do j=1,3
        ! get indices
        ! x=(i-1)/n
        ! x*n+1=i
        v1(j)=v1(j)*qp%griddensity(j)+1.0_r8
        !v1(j)=lo_index_in_periodic_array(int(anint(v1(j))),qp%griddensity(j))
    enddo
    ! convert to integer
    i=int(anint(v1))
end function

!> Get coordinates for q-points from indices.
!> $$ x = \frac{i-1}{n}  $$
!> where \(x\) the coordinate, \(i\) the index and \(n\) the number of
!> gridpoints in this dimension.
module pure function fft_coordinate_from_index(qp,i) result(r)
    !> the mesh
    class(lo_fft_mesh), intent(in) :: qp
    !> the indices
    integer, dimension(3), intent(in) :: i
    !> the resulting coordinates in fractional coordinates
    real(r8), dimension(3) :: r
    !
    integer :: j
    do j=1,3
        r(j)=( i(j)-1.0_r8)/(qp%griddensity(j)*1.0_r8)
    enddo
end function

!> Chop a q-point grid into tetrahedrons
module subroutine tesselate_qgrid(qp,p,mw,mem,verbosity)
    !> The q-point grid
    class(lo_qpoint_grid), intent(inout) :: qp
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:,:), allocatable :: r
    real(r8) :: t0,t1
    integer, dimension(:,:), allocatable :: shell_member
    integer, dimension(:), allocatable :: shell_ctr,prototype_shell
    integer :: np

    ! Start timer
    t0=walltime()
    t1=t0

    ! Set some basic things
    init: block
        real(r8), dimension(3) :: v0
        integer :: i,ix,iy,iz

        if ( verbosity .gt. 0 ) write(lo_iou,*) '... tesselating grid'

        ! How many centers of tetrahedrons are there?
        np=qp%n_full_point
        ! Make some space for tetrahedron centers.
        call mem%allocate(r,[3,np],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        r=0.0_r8
        ! Store coordinates of tetrahedron centers.
        do i=1,qp%n_full_point
            ! Make it parallel
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            ix=qp%ind2gridind(1,i)
            iy=qp%ind2gridind(2,i)
            iz=qp%ind2gridind(3,i)
            select type(qp)
            type is(lo_monkhorst_pack_mesh)
                v0(1)=( 2*(ix + 0.5_r8) - qp%griddensity(1) - 1 )/( 2*qp%griddensity(1) )
                v0(2)=( 2*(iy + 0.5_r8) - qp%griddensity(2) - 1 )/( 2*qp%griddensity(2) )
                v0(3)=( 2*(iz + 0.5_r8) - qp%griddensity(3) - 1 )/( 2*qp%griddensity(3) )
            type is(lo_fft_mesh)
                v0(1)=( ix-0.5_r8 )/real(qp%griddensity(1),r8)
                v0(2)=( iy-0.5_r8 )/real(qp%griddensity(2),r8)
                v0(3)=( iz-0.5_r8 )/real(qp%griddensity(3),r8)
            end select
            v0=lo_clean_fractional_coordinates(lo_chop(v0,1E-11_r8))
            v0=lo_clean_fractional_coordinates(lo_chop(v0,1E-11_r8))
            r(:,i)=v0
        enddo
        ! And sync
        call mw%allreduce('sum',r)
    end block init

    ! Reduce the points by symmetry by using a connectionlist instead.
    clsred: block
        real(r8), dimension(3) :: v0,v1,w0
        integer, dimension(:,:), allocatable :: connections
        integer, dimension(:), allocatable :: ind
        integer :: nsym,ix,iy,iz
        integer :: ipt,jpt,iop,i,l,ctr

        ! Number of operations to consider (ugly hack to cover the case of switched-off symmetries)
        nsym=size(qp%operationok)

        ! Space for connections
        if ( qp%timereversal ) then
            call mem%allocate(connections,[np,nsym*2],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        else
            call mem%allocate(connections,[np,nsym],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        endif
        connections=0

        if ( verbosity .gt. 0 ) call lo_progressbar_init()
        ! Count number of points on this rank, for the progressbar.
        ctr=0
        do i=1,np
            if ( mod(i,mw%n) .eq. mw%r ) ctr=ctr+1
        enddo

        l=0
        do ipt=1,np
            ! Make it parallel
            if ( mod(ipt,mw%n) .ne. mw%r ) cycle
            ! Point to be tested.
            v0=r(:,ipt)
            ! Check with all symmetry operations.
            do iop=1,nsym
                if ( qp%operationok(iop) .eqv. .false. ) cycle
                ! rotate the vector.
                v1=lo_operate_on_vector(p%sym%op(iop),v0,reciprocal=.true.,fractional=.true.,inverse=.false.)
                ! Check if it is on the grid?
                select type(qp)
                type is(lo_monkhorst_pack_mesh)
                    w0(1)=qp%griddensity(1)*(1+2*v1(1))*0.5_r8
                    w0(2)=qp%griddensity(2)*(1+2*v1(2))*0.5_r8
                    w0(3)=qp%griddensity(3)*(1+2*v1(3))*0.5_r8
                type is(lo_fft_mesh)
                    w0(1)=(1+2*qp%griddensity(1)*v1(1))*0.5_r8
                    w0(2)=(1+2*qp%griddensity(2)*v1(2))*0.5_r8
                    w0(3)=(1+2*qp%griddensity(3)*v1(3))*0.5_r8
                end select
                ! Skip if not on grid.
                if ( sum(abs(w0-anint(w0))) .gt. lo_sqtol ) cycle
                ! This is on the grid, get the indices
                ix=lo_index_in_periodic_array( int(anint(w0(1))), qp%griddensity(1) )
                iy=lo_index_in_periodic_array( int(anint(w0(2))), qp%griddensity(2) )
                iz=lo_index_in_periodic_array( int(anint(w0(3))), qp%griddensity(3) )
                connections(ipt,iop)=qp%gridind2ind(ix,iy,iz)
            enddo

            ! And maybe again, in case we time reversal symmetry
            if ( qp%timereversal ) then
                do iop=1,nsym
                    if ( qp%operationok(iop) .eqv. .false. ) cycle
                    ! rotate the vector.
                    v1=-lo_operate_on_vector(p%sym%op(iop),v0,reciprocal=.true.,fractional=.true.,inverse=.false.)
                    ! Check if it is on the grid?
                    select type(qp)
                    type is(lo_monkhorst_pack_mesh)
                        w0(1)=qp%griddensity(1)*(1+2*v1(1))*0.5_r8
                        w0(2)=qp%griddensity(2)*(1+2*v1(2))*0.5_r8
                        w0(3)=qp%griddensity(3)*(1+2*v1(3))*0.5_r8
                    type is(lo_fft_mesh)
                        w0(1)=(1+2*qp%griddensity(1)*v1(1))*0.5_r8
                        w0(2)=(1+2*qp%griddensity(2)*v1(2))*0.5_r8
                        w0(3)=(1+2*qp%griddensity(3)*v1(3))*0.5_r8
                    end select
                    ! Skip if not on grid.
                    if ( sum(abs(w0-anint(w0))) .gt. lo_sqtol ) cycle
                    ! This is on the grid, get the indices
                    ix=lo_index_in_periodic_array( int(anint(w0(1))), qp%griddensity(1) )
                    iy=lo_index_in_periodic_array( int(anint(w0(2))), qp%griddensity(2) )
                    iz=lo_index_in_periodic_array( int(anint(w0(3))), qp%griddensity(3) )
                    connections(ipt,nsym+iop)=qp%gridind2ind(ix,iy,iz)
                enddo
            endif

            if ( verbosity .gt. 0 ) then
                l=l+1
                if ( lo_trueNtimes(l,127,ctr) ) call lo_progressbar(' ... building connections',l,ctr,walltime()-t0)
            endif
        enddo
        ! Sync the connectionlist
        call mw%allreduce('sum',connections)
        ! Make final report
        if ( verbosity .gt. 0 ) then
            t1=walltime()
            call lo_progressbar(' ... building connections',ctr,ctr,t1-t0)
            t0=t1
        endif

        ! This little thingy actually gets the list of unique q-points. Neat.
        call mem%allocate(ind,np,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ind=1
        do ipt=1,np
            if ( mod(ipt,mw%n) .ne. mw%r ) cycle
            do iop=1,size(connections,2)
                jpt=abs(connections(ipt,iop))
                ! j*i=k, k>i then we remove k?
                ! some clarification, j is an operation, i and k are q-points. If k can be constructed from
                ! i, and k>i, we don't need k.
                if ( jpt .gt. ipt ) ind(jpt)=0
            enddo
        enddo
        call mw%allreduce('min',ind)

        ! Get the list of unique points
        l=sum(ind)
        call mem%allocate(prototype_shell,l,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(shell_ctr,l,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        prototype_shell=0
        shell_ctr=0
        l=0
        do i=1,np
            if ( ind(i) .eq. 1 ) then
                l=l+1
                prototype_shell(l)=i
            endif
        enddo

        call mem%deallocate(ind,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(ind,size(connections,2),persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        ind=0

        ! Count members per shell
        shell_ctr=0
        do ipt=1,size(prototype_shell)
            if ( mod(ipt,mw%n) .ne. mw%r ) cycle
            ind=0
            l=0
            sl1: do iop=1,nsym
                jpt=connections( prototype_shell(ipt),iop )

                if ( jpt .eq. 0 ) cycle sl1
                do i=1,l
                    if ( jpt .eq. ind(i) ) cycle sl1
                enddo
                l=l+1
                ind(l)=jpt
            enddo sl1

            if ( qp%timereversal ) then
                sl2: do iop=1,nsym
                    jpt=connections( prototype_shell(ipt),iop+nsym )
                    if ( jpt .eq. 0 ) cycle sl2
                    do i=1,l
                        if ( jpt .eq. ind(i) ) cycle sl2
                    enddo
                    l=l+1
                    ind(l)=jpt
                enddo sl2
            endif
            shell_ctr(ipt)=l
        enddo
        call mw%allreduce('sum',shell_ctr)

        ! Store members per shell
        l=maxval(shell_ctr)
        call mem%allocate(shell_member,[l,size(prototype_shell)],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        shell_member=0
        do ipt=1,size(prototype_shell)
            if ( mod(ipt,mw%n) .ne. mw%r ) cycle
            ind=0
            l=0
            sl3: do iop=1,nsym
                jpt=connections( prototype_shell(ipt),iop )

                if ( jpt .eq. 0 ) cycle sl3
                do i=1,l
                    if ( jpt .eq. ind(i) ) cycle sl3
                enddo
                l=l+1
                ind(l)=jpt
            enddo sl3

            if ( qp%timereversal ) then
                sl4: do iop=1,nsym
                    jpt=connections( prototype_shell(ipt),iop+nsym )
                    if ( jpt .eq. 0 ) cycle sl4
                    do i=1,l
                        if ( jpt .eq. ind(i) ) cycle sl4
                    enddo
                    l=l+1
                    ind(l)=jpt
                enddo sl4
            endif
            shell_member(1:l,ipt)=ind(1:l)
        enddo
        call mw%allreduce('sum',shell_member)

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... defined shells (',tochar(t1-t0),'s)'
            t0=t1
        endif

        call mem%deallocate(ind,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(connections,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block clsred

    ! Here I make sure that the prototype shells are all bunched together as closely as possible.
    !shuffleshells: block
    !end block shuffleshells

    ! Build tetrahedrons somehow.
    buildtet: block
        real(r8), dimension(:,:), allocatable :: center_it,center_at
        integer, dimension(:,:,:), allocatable :: shift_it,shift_at
        integer, dimension(:,:), allocatable :: ind_it,ind_at
        integer, dimension(:), allocatable :: offset_shell

        real(r8), dimension(3,4,6) :: tc0
        real(r8), dimension(3,6) :: tc1
        real(r8), dimension(3,8) :: bc0
        real(r8), dimension(3) :: v0,v1,v2,w0,w1
        integer, dimension(4,6) :: tetind,tetjnd
        integer, dimension(8) :: boxind
        integer, dimension(3) :: gi
        integer :: i,j,ii,jj
        integer :: ibox,jbox,iop,itet,jtet
        integer :: ctr0,ctr1,ctr2,ctr3
        integer :: nt_irr,nt_tot

        ! Make a little space for tetrahedron indices:
        nt_irr=size(prototype_shell)*6
        nt_tot=qp%n_full_point*6
        call mem%allocate(ind_it,[4,nt_irr],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(ind_at,[4,nt_tot],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(center_it,[3,nt_irr],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(center_at,[3,nt_tot],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(shift_it,[3,4,nt_irr],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(shift_at,[3,4,nt_tot],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        ind_it=0
        ind_at=0
        center_it=0.0_r8
        center_at=0.0_r8
        shift_it=0
        shift_at=0

        ! And a dummy counter thing
        call mem%allocate(offset_shell,size(prototype_shell),persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        offset_shell=0
        j=0
        do i=1,size(prototype_shell)
            offset_shell(i)=j
            j=j+shell_ctr(i)*6
        enddo

        ! Now I have a little bit of information I think.
        ctr0=0
        ptloop: do ii=1,size(prototype_shell)
            ! Make it parallel
            if ( mod(ii,mw%n) .ne. mw%r ) cycle

            ibox=prototype_shell(ii)
            v0=r(:,ibox)
            call neighbour_corners_from_center(qp,v0,boxind,bc0)
            call chop_box_into_tetrahedrons(p,bc0,tetind,tc0,tc1)

            ! Store this irreducible tetrahedron
            do i=1,6
                itet=(ii-1)*6+i
                ind_it(:,itet)=boxind(tetind(:,i))
                center_it(:,itet)=tc1(:,i)
                do j=1,4
                    ! absolute coordinate of point
                    w0=qp%ap( ind_it(j,itet) )%r
                    ! to fractional
                    w0=matmul(p%inv_reciprocal_latticevectors,w0)
                    ! distance to box corner
                    w0=tc0(:,j,i)-w0
                    ! Store as lattice vector shift
                    shift_it(:,j,itet)=anint(w0)
                enddo
            enddo

            ctr1=0
            ! Now go over all the points in this shell.
            shloop: do jj=1,shell_ctr(ii)
                jbox=shell_member(jj,ii)
                v1=r(:,jbox)
                ctr2=0

                ! Try to find a neat operation
                opl1: do iop=1,p%sym%n
                    v2=matmul(p%sym%op(iop)%rfm,v0)-v1
                    v2=lo_clean_fractional_coordinates(v2+0.5_r8)-0.5_r8
                    if ( lo_sqnorm(v2) .gt. lo_sqtol ) cycle

                    ! Ok, at least the center point matches. Probably a good idea
                    ! to rotate each tetrahedron on its own?
                    tetjnd=-1
                    ctr3=0
                    do i=1,6
                    do j=1,4
                        w0=lo_clean_fractional_coordinates( matmul(p%sym%op(iop)%rfm,tc0(:,j,i)) )
                        if ( qp%is_point_on_grid(w0) ) then
                            ! Note that we have a match
                            ctr3=ctr3+1
                            ! Fetch the index on the grid
                            gi=qp%index_from_coordinate(w0)
                            tetjnd(j,i)=qp%gridind2ind(gi(1),gi(2),gi(3))
                        endif
                    enddo
                    enddo
                    ! If we found them all, things are good!
                    if ( ctr3 .eq. 6*4 ) then
                        ctr2=ctr2+1
                        ctr1=ctr1+1
                        ! Make sure to store this tetrahedron
                        do i=1,6
                            jtet=offset_shell(ii)+(jj-1)*6+i
                            ind_at(:,jtet)=tetjnd(:,i)
                            center_at(:,jtet)=lo_clean_fractional_coordinates( matmul(p%sym%op(iop)%rfm,tc1(:,i)) )
                            do j=1,4
                                w0=matmul(p%sym%op(iop)%rfm,tc0(:,j,i))
                                shift_at(:,j,jtet)=int(anint(v0-lo_clean_fractional_coordinates(v0)))
                            enddo
                            do j=1,4
                                w1=matmul(p%sym%op(iop)%rfm,tc0(:,j,i))
                                ! absolute coordinate of point
                                w0=qp%ap( ind_at(j,jtet) )%r
                                ! to fractional
                                w0=matmul(p%inv_reciprocal_latticevectors,w0)
                                ! distance to box corner
                                w0=w1-w0
                                ! Store as lattice vector shift
                                shift_at(:,j,jtet)=anint(w0)
                            enddo
                        enddo
                        exit opl1
                    endif
                enddo opl1

                ! If we found it, all's well.
                if ( ctr2 .gt. 0 ) then
                    cycle shloop
                endif

                ! New attempt, with inversion added.
                opl2: do iop=1,p%sym%n
                    v2=-matmul(p%sym%op(iop)%rfm,v0)-v1
                    v2=lo_clean_fractional_coordinates(v2+0.5_r8)-0.5_r8
                    if ( lo_sqnorm(v2) .gt. lo_sqtol ) cycle

                    ! Ok, at least the center point matches. Probably a good idea
                    ! to rotate each tetrahedron on its own?
                    tetjnd=-1
                    ctr3=0
                    do i=1,6
                    do j=1,4
                        w0=lo_clean_fractional_coordinates( -matmul(p%sym%op(iop)%rfm,tc0(:,j,i)) )
                        if ( qp%is_point_on_grid(w0) ) then
                            ! Note that we have a match
                            ctr3=ctr3+1
                            ! Fetch the index on the grid
                            gi=qp%index_from_coordinate(w0)
                            tetjnd(j,i)=qp%gridind2ind(gi(1),gi(2),gi(3))
                        endif
                    enddo
                    enddo
                    ! If we found them all, things are good!
                    if ( ctr3 .eq. 6*4 ) then
                        ctr2=ctr2+1
                        ctr1=ctr1+1
                        ! Make sure to store this tetrahedron
                        do i=1,6
                            jtet=offset_shell(ii)+(jj-1)*6+i
                            ind_at(:,jtet)=tetjnd(:,i)
                            center_at(:,jtet)=lo_clean_fractional_coordinates( -matmul(p%sym%op(iop)%rfm,tc1(:,i)) )
                            do j=1,4
                                w1=-matmul(p%sym%op(iop)%rfm,tc0(:,j,i))
                                ! absolute coordinate of point
                                w0=qp%ap( ind_at(j,jtet) )%r
                                ! to fractional
                                w0=matmul(p%inv_reciprocal_latticevectors,w0)
                                ! distance to box corner
                                w0=w1-w0
                                ! Store as lattice vector shift
                                shift_at(:,j,jtet)=anint(w0)
                            enddo
                        enddo
                        exit opl2
                    endif
                enddo opl2
            enddo shloop

            ! Now a sanity test:
            if ( ctr1 .ne. shell_ctr(ii) ) then
                ! This did not work out.
                ctr0=-1
                exit ptloop
            else
                ! Went fine
                ctr0=ctr0+1
            endif
        enddo ptloop

        ! Sync, and sanity checks
        call mw%allreduce('sum',ind_it)
        call mw%allreduce('sum',ind_at)
        call mw%allreduce('sum',center_it)
        call mw%allreduce('sum',center_at)
        call mw%allreduce('sum',shift_it)
        call mw%allreduce('sum',shift_at)
        call mw%allreduce('sum',ctr0)
        if ( minval(ind_it) .le. 0               ) ctr0=ctr0+1
        if ( maxval(ind_it) .gt. qp%n_full_point ) ctr0=ctr0+1
        if ( minval(ind_at) .le. 0               ) ctr0=ctr0+1
        if ( maxval(ind_at) .gt. qp%n_full_point ) ctr0=ctr0+1

        if ( verbosity .gt. 0 ) then
            t1=walltime()
            write(lo_iou,*) '... rotated tetrahedrons (',tochar(t1-t0),'s)'
            t0=t1
        endif

        ! Now decide on what to do going forward:
        if ( ctr0 .ne. size(prototype_shell) ) then
            ! This did not go well. Skip symmetrization of tetrahedrons.
            if ( verbosity .gt. 0 ) then
                write(lo_iou,*) '... mesh not quite commensurate with symmetry'
            endif

            ! Just regenerate all the tetrahedrons the stupid way.
            ind_at=0
            shift_at=0
            center_at=0
            do ii=1,qp%n_full_point
                if ( mod(ii,mw%n) .ne. mw%r ) cycle
                v0=r(:,ii)
                call neighbour_corners_from_center(qp,v0,boxind,bc0)
                call chop_box_into_tetrahedrons(p,bc0,tetind,tc0,tc1)
                ! Store this tetrahedron
                do i=1,6
                    itet=(ii-1)*6+i
                    ind_at(:,itet)=boxind(tetind(:,i))
                    center_at(:,itet)=tc1(:,i)
                    do j=1,4
                        ! absolute coordinate of point
                        w0=qp%ap( ind_at(j,itet) )%r
                        ! to fractional
                        w0=matmul(p%inv_reciprocal_latticevectors,w0)
                        ! distance to box corner
                        w0=tc0(:,j,i)-w0
                        ! Store as lattice vector shift
                        shift_at(:,j,itet)=anint(w0)
                    enddo
                enddo
            enddo
            call mw%allreduce('sum',ind_at)
            call mw%allreduce('sum',center_at)
            call mw%allreduce('sum',shift_at)

            ! Make space for tetrahedrons:
            qp%n_irr_tet=6*qp%n_full_point
            qp%n_full_tet=6*qp%n_full_point
            allocate(qp%it(qp%n_irr_tet))
            allocate(qp%at(qp%n_full_tet))
            ! Store tetrahedrons
            do itet=1,qp%n_full_tet
                qp%at(itet)%full_index=ind_at(:,itet)
                qp%it(itet)%full_index=ind_at(:,itet)
                qp%at(itet)%integration_weight=1.0_r8
                qp%it(itet)%integration_weight=1.0_r8
                qp%at(itet)%center_of_mass=center_at(:,itet)
                qp%it(itet)%center_of_mass=center_at(:,itet)
                qp%at(itet)%lattice_vector_shift=shift_at(:,:,itet)
                qp%it(itet)%lattice_vector_shift=shift_at(:,:,itet)
            enddo
        else
            ! This went well!
            if ( verbosity .gt. 0 ) then
                write(lo_iou,*) '... mesh commensurate with symmetry'
            endif
            qp%n_irr_tet=nt_irr
            qp%n_full_tet=nt_tot
            allocate(qp%it(qp%n_irr_tet))
            allocate(qp%at(qp%n_full_tet))

            do ii=1,size(prototype_shell)
                do i=1,6
                    itet=(ii-1)*6+i
                    qp%it(itet)%full_index=ind_it(:,itet)
                    qp%it(itet)%integration_weight=shell_ctr(ii)
                    qp%it(itet)%center_of_mass=center_it(:,itet)
                    qp%it(itet)%lattice_vector_shift=shift_it(:,:,itet)
                enddo
            enddo
            do i=1,qp%n_full_tet
                qp%at(i)%full_index=ind_at(:,i)
                qp%at(i)%integration_weight=1.0_r8
                qp%at(i)%center_of_mass=center_at(:,i)
                qp%at(i)%lattice_vector_shift=shift_at(:,:,i)
            enddo
        endif

        ! Now we can clean a lot of stuff!
        call mem%deallocate(offset_shell,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ind_it,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(ind_at,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(center_it,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(center_at,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(shift_it,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(shift_at,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(r,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(shell_ctr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(shell_member,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(prototype_shell,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block buildtet

    ! And now I can can set the missing indices as well as the integration weights.
    tetwts: block
        real(r8), dimension(:), allocatable :: wt
        real(r8), dimension(3,4) :: tet
        real(r8), dimension(3) :: v0,v1
        real(r8) :: f0
        integer :: itet,icrn,ikp,jkp

        ! The appropriate weight for one tetrahedron:
        f0=1.0_r8/qp%n_full_tet
        ! Sort out the weights? First the irreducible I think.
        call mem%allocate(wt,qp%n_irr_tet,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        wt=0.0_r8
        do itet=1,qp%n_irr_tet
            ! Make it parallel
            if ( mod(itet,mw%n) .ne. mw%r ) cycle
            ! Get unwrapped coordinates of tetrahedron
            v0=qp%it(itet)%center_of_mass
            do icrn=1,4
                v1=qp%ap( qp%it(itet)%full_index(icrn) )%r
                v1=matmul(p%inv_reciprocal_latticevectors,v1)-v0
                v1=v1+qp%it(itet)%lattice_vector_shift(:,icrn)
                tet(:,icrn)=matmul(p%reciprocal_latticevectors,v1)
            enddo
            if ( abs(lo_unsigned_tetrahedron_volume(tet)*p%volume-f0) .gt. 1E-10_r8 ) then
                call lo_stop_gracefully(['Unexpected tetrahedron weight'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            else
                wt(itet)=f0*qp%it(itet)%integration_weight
            endif
        enddo
        call mw%allreduce('sum',wt)

        ! Sanity check, do it again for the full set
        do itet=1,qp%n_full_tet
            ! Make it parallel
            if ( mod(itet,mw%n) .ne. mw%r ) cycle
            ! Get unwrapped coordinates of tetrahedron
            v0=qp%at(itet)%center_of_mass
            do icrn=1,4
                v1=qp%ap( qp%at(itet)%full_index(icrn) )%r
                v1=matmul(p%inv_reciprocal_latticevectors,v1)-v0
                v1=v1+qp%at(itet)%lattice_vector_shift(:,icrn)
                tet(:,icrn)=matmul(p%reciprocal_latticevectors,v1)
            enddo
            if ( abs(lo_unsigned_tetrahedron_volume(tet)*p%volume-f0) .gt. 1E-10_r8 ) then
                call lo_stop_gracefully(['Unexpected tetrahedron weight'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
            endif
        enddo



        ! Now store the proper weights, and sort out the irreducible indices
        do itet=1,qp%n_irr_tet
            qp%it(itet)%integration_weight=wt(itet)
            do icrn=1,4
                ikp=qp%it(itet)%full_index(icrn)
                jkp=qp%ap(ikp)%irreducible_index
                qp%it(itet)%irreducible_index(icrn)=jkp
            enddo
        enddo
        ! And then the same thing for all the tetrahedrons
        do itet=1,qp%n_full_tet
            qp%at(itet)%integration_weight=f0
            do icrn=1,4
                ikp=qp%at(itet)%full_index(icrn)
                jkp=qp%ap(ikp)%irreducible_index
                qp%at(itet)%irreducible_index(icrn)=jkp
            enddo
        enddo

        call mem%deallocate(wt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block tetwts

    if ( verbosity .gt. 0 ) then
        t1=walltime()
        write(lo_iou,*) '... finalized tetrahedron construction (',tochar(t1-t0),'s)'
        t0=t1
    endif
end subroutine

!> get the 'box' of points that surround a box center.
subroutine neighbour_corners_from_center(qp,r,ind,br)
    !> q-point grid
    class(lo_qpoint_grid), intent(in) :: qp
    !> point in question
    real(r8), dimension(3), intent(in) :: r
    !> indices to corners
    integer, dimension(8), intent(out) :: ind
    !> actual coordinates
    real(r8), dimension(3,8), intent(out) :: br

    real(r8), dimension(3) :: v
    real(r8), dimension(2) :: dlx,dly,dlz
    integer, dimension(3) :: gi
    integer :: i,j,k,l

    ! The spacing of the grid:
    dlx(1)=-0.5_r8/qp%griddensity(1)
    dlx(2)= 0.5_r8/qp%griddensity(1)
    dly(1)=-0.5_r8/qp%griddensity(2)
    dly(2)= 0.5_r8/qp%griddensity(2)
    dlz(1)=-0.5_r8/qp%griddensity(3)
    dlz(2)= 0.5_r8/qp%griddensity(3)

    ! This is a sanity check to make sure I know how to index, basically.
    l=0
    do i=1,2
    do j=1,2
    do k=1,2
        v=r+[dlx(i),dly(j),dlz(k)]
        l=l+1
        if ( qp%is_point_on_grid(v) ) then
            ! it's on the grid, get what point it transforms to
            gi=qp%index_from_coordinate(v)
            ind(l)=qp%gridind2ind(gi(1),gi(2),gi(3))
        else
            call lo_stop_gracefully(['Did not think correctly when it comes to tetrahedrons.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
    enddo
    enddo
    enddo
    ! Then it's better to return the box in relative coordinates, since I actually know that it should be.
    l=0
    do i=1,2
    do j=1,2
    do k=1,2
        v=r+[dlx(i),dly(j),dlz(k)]
        l=l+1
        br(:,l)=v
    enddo
    enddo
    enddo

    ! ! And eventually fetch the coordinates?
    ! do i=1,8
    !     gi=qp%ind2gridind(:,ind(i))
    !     v=qp%coordinate_from_index(gi)
    !     v=lo_clean_fractional_coordinates(lo_chop(v,1E-11_r8))
    !     v=lo_clean_fractional_coordinates(lo_chop(v,1E-11_r8))
    !     br(:,i)=v
    ! enddo
end subroutine

!> slice a box into tetrahedrons
subroutine chop_box_into_tetrahedrons(p,box,tetind,tetcoord,tetctr)
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> coordinates of box (in fractional coordinates)
    real(r8), dimension(3,8), intent(in) :: box
    !> indices to tetrahedrons
    integer, dimension(4,6), intent(out) :: tetind
    !> coordinates of tetrahedrons, fractional coordinates.
    real(r8), dimension(3,4,6), intent(out) :: tetcoord
    !> center of mass of the tetrahedrons
    real(r8), dimension(3,6), intent(out) :: tetctr

    type(lo_plane) :: plane
    real(r8), dimension(3,8) :: r
    real(r8), dimension(3,6) :: u
    real(r8), dimension(3) :: v0,v1,com
    real(r8) :: f0,f1,f2
    integer, dimension(6) :: otherpoints,angleind
    integer, dimension(2) :: diagonal
    integer :: i,j,k,l

    ! First we have to unwrap the box and get the center of mass.
    v0=box(:,1) !+[1.0_r8,2.0_r8,3.0_r8]*1E-5
    com=0.0_r8
    r=0.0_r8
    do i=1,8
        v1=box(:,i)-v0
        !v1=lo_clean_fractional_coordinates(v1+0.5_r8)-0.5_r8
        r(:,i)=matmul(p%reciprocal_latticevectors,v1)
        com=com+r(:,i)/8.0_r8
    enddo

    ! Then I have to find the shortest diagonal.
    f2=lo_huge
    diagonal=-1
    do i=1,8
    do j=i+1,8
        v0=r(:,i)-com
        v1=r(:,j)-com
        f0=0.0_r8
        f0=f0+abs(dot_product(v0,v1)+dot_product(v0,v0))
        f0=f0+abs(dot_product(v0,v1)+dot_product(v0,v0))
        if ( f0 .lt. 1E-10_r8 ) then
            f1=norm2(r(:,j)-r(:,i))
            if ( f1 .lt. f2 ) then
                f2=f1
                diagonal=[i,j]
            endif
        endif
    enddo
    enddo

    ! Also keep the points not on the diagonal.
    l=0
    do i=1,8
        if ( i .ne. diagonal(1) .and. i .ne. diagonal(2) ) then
            l=l+1
            otherpoints(l)=i
            u(:,l)=r(:,i)
        endif
    enddo

    ! Get a plane with the normal aligned with the shortest diagonal.
    v0=r(:,diagonal(1))-r(:,diagonal(2))
    v0=v0/norm2(v0)
    call plane%generate(normal=v0,point=com)
    call plane%anglesort(u,angleind)
    ! get the tetrahedron mapping
    do i=1,5
        tetind(1,i)=diagonal(1)
        tetind(2,i)=diagonal(2)
        tetind(3,i)=otherpoints(angleind(i))
        tetind(4,i)=otherpoints(angleind(i+1))
    enddo
    tetind(1,6)=diagonal(1)
    tetind(2,6)=diagonal(2)
    tetind(3,6)=otherpoints(angleind(6))
    tetind(4,6)=otherpoints(angleind(1))

    ! Get the center of the tetrahedron, think that makes some sense
    ! at least.
    tetctr=0.0_r8
    do i=1,6
        v0=0.0_r8
        do j=1,4
            k=tetind(j,i)
            v0=v0+box(:,k)*0.25_r8
        enddo
        tetctr(:,i)=lo_clean_fractional_coordinates(v0)
    enddo

    ! And return the coordinates of the tetrahedrons.
    do i=1,6
    do j=1,4
        k=tetind(j,i)
        tetcoord(:,j,i)=box(:,k)
    enddo
    enddo
end subroutine

!> check if two points are on top of each other
pure function compare_qpoints_with_pbc(u,v,tol) result(match)
    ! the first q-point
    real(r8), dimension(3), intent(in) :: u
    ! the second q-point
    real(r8), dimension(3), intent(in) :: v
    ! the tolerance
    real(r8), intent(in) :: tol
    ! is it a match
    logical :: match
    !
    real(r8), dimension(3) :: r
    r=lo_clean_fractional_coordinates(u-v+0.5_r8)-0.5_r8
    if ( lo_sqnorm(r) .lt. tol ) then
        match=.true.
    else
        match=.false.
    endif
end function

end submodule
