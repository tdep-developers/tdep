module lo_distribute_and_classify_things
!!
!! Helper routines that distribute and sort/classify things in various ways.
!!
use konstanter, only: r8,lo_iou,lo_huge,lo_hugeint,lo_exitcode_symmetry,lo_exitcode_param
use gottochblandat, only: qsort,walltime,tochar,lo_symmetric_eigensystem_3x3matrix
use lo_memtracker, only: lo_mem_helper
use lo_sorting, only: lo_return_unique,lo_return_unique_indices,lo_qsort
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully

implicit none
private
public :: lo_minmax_division_of_points
public :: lo_coordination_shells_from_permutation_list
public :: lo_distribute_weighted

!> helper type to assist with partitioning
type lo_set_of_points
    !> how many points
    integer :: n=-lo_hugeint
    !> Which points are in this set
    integer, dimension(:), allocatable :: ind
end type

contains

!> distribute a set of weighted things across ranks (that are already distributed, but badly)
subroutine lo_distribute_weighted(mw,weight,n_elem,elem,mem)
    !> MPI communicator to distribute across
    type(lo_mpi_helper), intent(in) :: mw
    !> weight of things I am to distribute
    integer, dimension(:), intent(in) :: weight
    !> how many elements does this rank get?
    integer, intent(out) :: n_elem
    !> which elements does this rank get
    integer, dimension(:), allocatable, intent(out) :: elem
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    integer, dimension(:), allocatable :: buf_wt,buf_ind
    integer :: n,i,j,ihi,ilo,iter,niter

    ! The purpose of this routine is to distribute things over MPI ranks. Each of the things has a certain
    ! weight, or cost, associated with it. We want to ensure that we get a distribution of sums of weights
    ! over the ranks that is about as constant as possible. Pretty sure it is an NP-hard problem in general,
    ! e.g. box-packing, so I'm just going for something that is not horrible. If load-balancing becomes an
    ! issue eventually, I would have to revise.

    ! Number of things we have to distribute.
    n=size(weight)
    ! We need a little temporary space
    call mem%allocate(buf_wt,n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(buf_ind,n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    buf_wt=0
    buf_ind=0

    ! Start by sorting the things by weight
    buf_wt=weight
    call lo_qsort(buf_wt,buf_ind)

    ! Now count things per rank.
    n_elem=0
    do i=1,n
        if ( mod(i,mw%n) .eq. mw%r ) n_elem=n_elem+1
    enddo
    if ( n_elem .gt. 0 ) then
        call mem%allocate(elem,n_elem,persistent=.true.,scalable=.true.,file=__FILE__,line=__LINE__)
        elem=0
    else
        ! Dummy allocate, and return early in case there are no elements on this rank.
        call mem%allocate(elem,1,persistent=.true.,scalable=.true.,file=__FILE__,line=__LINE__)
        elem=-lo_hugeint
        return
    endif

    ! A little convoluted way of looping over values, but should help with load balancing.
    niter=0
    do
        if ( niter*mw%n .gt. n ) then
            exit
        else
            niter=niter+1
        endif
    enddo

    ! So I go through the sorted values, and alternate direction every nrank steps.
    ! that should neatly fill up the ranks in a somewhat evenly distributed manner.
    j=0
    n_elem=0
    do iter=1,niter
        ilo=(iter-1)*mw%n+1
        ihi=min(iter*mw%n,n)
        if ( mod(iter,2) .eq. 0 ) then
            do i=ihi,ilo,-1
                j=j+1
                if ( mod(j,mw%n) .eq. mw%r ) then
                    n_elem=n_elem+1
                    elem(n_elem)=buf_ind(i)
                endif
            enddo
        else
            do i=ilo,ihi
                j=j+1
                if ( mod(j,mw%n) .eq. mw%r ) then
                    n_elem=n_elem+1
                    elem(n_elem)=buf_ind(i)
                endif
            enddo
        endif
    enddo

    ! Cleanup
    call mem%deallocate(buf_wt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(buf_ind,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> build coordination shells from a list of transformation rules
subroutine lo_coordination_shells_from_permutation_list(perm,shell_ctr,shell_member,shell_index,mem,prototype_shell,mw)
    !> list of transformation rules
    integer, dimension(:,:), intent(in) :: perm
    !> counter for the number of elements in each shell
    integer, dimension(:), allocatable, intent(out) :: shell_ctr
    !> the members of each shell
    integer, dimension(:,:), allocatable, intent(out) :: shell_member
    !> which shell does each point belong to
    integer, dimension(:), allocatable, intent(out) :: shell_index
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> which are the unique points
    integer, dimension(:), allocatable, intent(out), optional :: prototype_shell
    !> mpi helper
    type(lo_mpi_helper), intent(inout), optional :: mw

    integer, dimension(:,:), allocatable :: pt,pu
    integer, dimension(:), allocatable :: di,dk
    integer :: i,j,np,no,n_shell

    ! This is a fairly general routine to aid in symmetry reduction. The things
    ! that are to be reduced are given as a list of integers. The action of a
    ! symmetry operation is to permute these indices. From the list of all
    ! permutaitons, in array perm(n_things,n_operations), we can deduce what the
    ! irreducible things are, and how they should be divided into irreducible groups.

    ! Some temporary space
    np=size(perm,1) ! number of elements
    no=size(perm,2) ! number of operations
    call mem%allocate(pt,[no,np],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(pu,[no,np],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    pt=0
    pu=transpose(perm)

    ! It can either be done in parallel or serially, so two version are given here.
    if ( present(mw) ) then
        ! So, sort each column.
        do i=1,size(pt,2)
            ! parallel over particles
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            pt(:,i)=pu(:,i)
            call lo_qsort(pt(:,i))
        enddo
        call mw%allreduce('sum',pt)

        ! The unique columns will be the coordination shells, I think.
        ! Should really make a parallel version of this somehow.
        call mem%allocate(shell_index,np,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        shell_index=0
        call lo_return_unique_indices(pt,di,mem,redind=shell_index)

        ! return which the unique are
        if ( present(prototype_shell) ) then
            call mem%allocate(prototype_shell,size(di),persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
            prototype_shell=di
        endif

        ! Store things as shells
        n_shell=size(di)
        call mem%allocate(shell_member,[no,n_shell],persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(shell_ctr,n_shell,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        shell_member=0
        shell_ctr=0
        do i=1,n_shell
            ! parallel over shells
            if ( mod(i,mw%n) .ne. mw%r ) cycle
            j=di(i)
            call lo_return_unique(pt(:,j),dk,mem)
            shell_ctr(i)=size(dk)
            shell_member(1:shell_ctr(i),i)=dk
            call mem%deallocate(dk,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        enddo
        ! Reduce over MPI.
        call mw%allreduce('sum',shell_ctr)
        call mw%allreduce('sum',shell_member)

        ! Small sanity test that I think always must hold
        if ( sum(shell_ctr) .ne. np ) then
            call lo_stop_gracefully(['Lost particle when dividing into shells.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
        ! Then I could think of one more test
        call mem%deallocate(di,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        !call mem%deallocate(di,persistent=.true.,scalable=.false.)

        call mem%allocate(di,np,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        di=-1
        do i=1,n_shell
        do j=1,shell_ctr(i)
            di( shell_member(j,i) )=di( shell_member(j,i) )+1
        enddo
        enddo
        if ( sum(abs(di)) .ne. 0 ) then
            call lo_stop_gracefully(['Cluster does not divide cleanly into shells.'],lo_exitcode_symmetry,__FILE__,__LINE__,mw%comm)
        endif
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    else
        ! Same thing as above, but not parallel.
        do i=1,size(pt,2)
            pt(:,i)=pu(:,i)
            call lo_qsort(pt(:,i))
        enddo

        ! The unique columns will be the coordination shells, I think.
        call mem%allocate(shell_index,np,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        shell_index=0
        call lo_return_unique_indices(pt,di,mem,redind=shell_index)

        ! return which the unique are
        if ( present(prototype_shell) ) then
            call mem%allocate(prototype_shell,size(di),persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
            prototype_shell=di
        endif

        ! Store things as shells
        n_shell=size(di)
        call mem%allocate(shell_member,[no,n_shell],persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(shell_ctr,n_shell,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        shell_member=0
        shell_ctr=0
        do i=1,n_shell
            j=di(i)
            call lo_return_unique(pt(:,j),dk,mem)
            shell_ctr(i)=size(dk)
            shell_member(1:shell_ctr(i),i)=dk
            call mem%deallocate(dk,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        enddo

        ! Small sanity test that I think always must hold
        if ( sum(shell_ctr) .ne. np ) then
            call lo_stop_gracefully(['Cluster does not divide cleanly into shells.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        ! Then I could think of one more test
        call mem%deallocate(di,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)

        call mem%allocate(di,np,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        di=-1
        do i=1,n_shell
        do j=1,shell_ctr(i)
            di( shell_member(j,i) )=di( shell_member(j,i) )+1
        enddo
        enddo
        if ( sum(abs(di)) .ne. 0 ) then
            call lo_stop_gracefully(['Cluster does not divide cleanly into shells.'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    endif

    ! And cleanup
    call mem%deallocate(pt,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(pu,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> Routine that min-max divides a list of points into evenly sized sets that are spatially localized.
subroutine lo_minmax_division_of_points(r,ri,nset,mem,verbosity)
    !> list of points to divide
    real(r8), dimension(:,:), intent(in) :: r
    !> resulting set for each point
    integer, dimension(:), intent(out) :: ri
    !> how many sets should they be divided in
    integer, intent(in) :: nset
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    integer, dimension(:), allocatable :: divisors
    integer, dimension(:,:), allocatable :: setctr
    real(r8) :: timer,t0,t1

    timer=walltime()
    t0=timer
    t1=timer

    ! Decide on the divisors
    init: block
        integer, dimension(:), allocatable :: di,dj,dk
        integer :: i,j,k,ll,nmult,nval
        character(len=1000) :: opf

        ! Some sanity tests never hurt anyone
        if ( nset .le. 0 ) then
            call lo_stop_gracefully(['The number of sets has to be a positive number'],lo_exitcode_param,__FILE__,__LINE__)
        endif
        if ( nset .gt. size(r,2) ) then
            call lo_stop_gracefully(['More sets than points to divide, makes no sense.'],lo_exitcode_param,__FILE__,__LINE__)
        endif

        ! It could be very very simple, i.e. just one set. Catch that right away.
        if ( nset .eq. 1 ) then
            ri=1
            return
        endif

        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) '... sorting '//tochar(size(r,2))//' points into '//tochar(nset)// ' sets'
        endif

        ! Prime-factorize the number of sets, that gives my how many divisions I should do.
        call lo_find_prime_factors(nset, divisors, mem)

        ! Now I need to know how things are to be divided in each step. I also
        ! want an even number of points per rank, in the end.
        ll=size(divisors)+1
        call mem%allocate(setctr,[nset,ll],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(di,nset,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(dj,ll,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(dk,ll,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        ! Count points per set
        di=0
        do i=1,size(r,2)
            j=mod(i,nset)+1
            di(j)=di(j)+1
        enddo

        ! Init the things I know. I know what the final distribution should be
        setctr=0
        setctr(:,ll)=di
        ! dummy thingy for the divisors
        dj(1)=1
        dj(2:ll)=divisors
        ! dummy thing with the number of sets per iteration
        dk=0
        j=1
        do i=1,ll
            j=j*dj(i)
            dk(i)=j
        enddo
        ! figure out the number of points per set backwards, from the final distribution.
        do i=ll-1,1,-1
            nmult=dj(i+1)
            nval=dk(i)
            do j=1,nval
            do k=1,nmult
                setctr(j,i)=setctr(j,i)+setctr(j+(k-1)*nval,i+1)
            enddo
            enddo
        enddo

        ! Report a little
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) '... decided on division:'
            k=10
            do i=1,ll

                if ( dk(i) .lt. k ) then
                    j=dk(i)
                    opf="(1X,I4,':',"//tochar(j)//"(1X,I9))"
                    write(lo_iou,opf) dk(i),setctr(1:j,i)
                else
                    opf="(1X,I4,':',"//tochar(k)//"(1X,I9))"
                    write(lo_iou,opf) dk(i),setctr(1:k,i)
                endif
            enddo
        endif

        ! And cleanup the dummys
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(dj,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(dk,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    end block init

    ! Now start dividing them according to plan.
    splitstuff: block
        type(lo_set_of_points), dimension(:), allocatable :: set
        integer, dimension(:), allocatable :: splitcount
        integer :: i,j,is,iter,nst,mult

        ! Start by assigning all the points to the first set.
        ri=1
        ! We only have one set
        nst=1
        ! Populate this set with all the points.
        allocate(set(nset))
        set(1)%n=size(r,2)
        call mem%allocate(set(1)%ind,set(1)%n,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        do i=1,set(1)%n
            set(1)%ind(i)=i
        enddo

        ! Start splitting
        iterloop: do iter=1,size(divisors)
            ! How many times should each set be split?
            mult=divisors(iter)
            call mem%allocate(splitcount,mult,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

            do is=1,nst
                ! Fetch how it should be split
                do i=1,mult
                    j=is+(i-1)*nst
                    splitcount(i)=setctr(j,iter+1)
                enddo
                ! Split it!
                call split_set_of_points(set,is,nst,r,splitcount,mem)
            enddo
            call mem%deallocate(splitcount,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

            ! Update the number of sets!
            nst=nst*mult
        enddo iterloop

        ! A small sanity test perhaps. Never hurts.
        do is=1,nset
            if ( set(is)%n .ne. setctr(is,size(divisors)+1) ) then
                call lo_stop_gracefully(['Could not divide the points the way I thought I could'],lo_exitcode_symmetry,__FILE__,__LINE__)
            endif
        enddo

        ! Now we are done, actually! Store the division!
        ri=0
        do is=1,nset
            do i=1,set(is)%n
                j=set(is)%ind(i)
                ri(j)=is
            enddo
        enddo

        ! Some small sanity tests to make sure I did not mess up.
        if ( minval(ri) .le. 0 ) then
            call lo_stop_gracefully(['Could not divide the points the way I thought I could'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif
        if ( maxval(ri) .gt. nset ) then
            call lo_stop_gracefully(['Could not divide the points the way I thought I could'],lo_exitcode_symmetry,__FILE__,__LINE__)
        endif

        ! And cleanup
        call mem%deallocate(divisors,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(setctr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        do is=1,nset
            call mem%deallocate(set(is)%ind,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        enddo
        deallocate(set)

        ! And say it went fine!
        if ( verbosity .gt. 0 ) then
            write(lo_iou,*) '... done splitting (',tochar(walltime()-timer),'s)'
        endif

    end block splitstuff
end subroutine

!> Split a set of points
subroutine split_set_of_points(set,is,nset,r,splitcount,mem)
    !> sets of points
    type(lo_set_of_points), dimension(:), intent(inout) :: set
    !> which set are we splitting
    integer, intent(in) :: is
    !> current number of sets
    integer, intent(in) :: nset
    !> coordinates of all points
    real(r8), dimension(:,:), intent(in) :: r
    !> How should the set be split?
    integer, dimension(:), intent(in) :: splitcount
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:,:), allocatable :: dr
    real(r8), dimension(:), allocatable :: dv0,dv1
    real(r8), dimension(3,3) :: m0,m1
    real(r8), dimension(3) :: com,v0,normal
    real(r8) :: f0
    integer, dimension(:), allocatable :: di,dj
    integer :: i,j,l,js,im

    ! Some dummy space needed
    call mem%allocate(di,set(is)%n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(dj,set(is)%n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(dr,[3,set(is)%n],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(dv0,set(is)%n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(dv1,set(is)%n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    di=0
    dj=0
    dr=0.0_r8
    dv0=0.0_r8
    dv1=0.0_r8

    ! Calculate the center of mass
    com=0.0_r8
    do i=1,set(is)%n
        com=com+r(:, set(is)%ind(i) )
    enddo
    com=com/real(set(is)%n,r8)

    ! Get the normal of the plane that passes through the center of mass.
    ! That will split things as evenly as possible. First get the vectors
    ! that point from the center of mass to the grid.
    do i=1,set(is)%n
        dr(:,i)=r(:, set(is)%ind(i) ) -com
    enddo
    ! calculate m0 = dr0*dr0^T
    m0=0.0_r8
    call dsyrk('U','N',3,set(is)%n,1.0_r8,dr,3,0.0_r8,m0,3)
    ! Fill out the missing pieces of the matrix, don't think
    ! it's necessary but does not hurt.
    m0(2,1)=m0(1,2)
    m0(3,1)=m0(1,3)
    m0(3,2)=m0(2,3)
    ! get eigenvalues
    call lo_symmetric_eigensystem_3x3matrix(m0,v0,m1)
    ! set the normal to the eigenvector corresponding to largest eigenvalue
    f0=-lo_huge
    do i=1,3
        if ( v0(i) .gt. f0 ) then
            f0=v0(i)
            normal=m1(:,i)/norm2(m1(:,i))
        endif
    enddo

    ! calculate distances to this plane
    f0=-dot_product(normal,com)
    do i=1,set(is)%n
        v0=r(:, set(is)%ind(i) )
        dv0(i)=dot_product(normal,v0)+f0
    enddo

    ! sort the distances
    dv1=dv0
    call qsort(dv1,dj)

    ! Store indices temporarily
    di=set(is)%ind

    ! Destroy current indices
    call mem%deallocate(set(is)%ind,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)

    ! Store the new sets.
    l=0
    do im=1,size(splitcount)
        ! What is the index of the new set?
        js=is+(im-1)*nset
        set(js)%n=splitcount(im)
        call mem%allocate(set(js)%ind,splitcount(im),persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        do i=1,splitcount(im)
            j=dj(l+i)
            set(js)%ind(i)=di(j)
        enddo
        l=l+splitcount(im)
    enddo

    ! And a little cleanup
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(dj,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(dr,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(dv0,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%deallocate(dv1,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> Prime factorization of numbers that are not particularly large.
subroutine lo_find_prime_factors(n,d,mem)
    !> number to factorize
    integer, intent(in) :: n
    !> array with factors such that product(d)=n, sorted largest to smallest.
    integer, dimension(:), allocatable, intent(out) :: d
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    integer :: div, next, rest
    integer :: i,ctr

    ! Quick catch for the simplest
    if ( n .eq. 1 ) then
        call mem%allocate(d,1,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        d=1
        return
    endif

    ctr=0
    i = 1
    div=2
    next=3
    rest=n
    do while ( rest .ne. 1 )
        do while ( mod(rest, div) == 0 )
            ctr=ctr+1
            i=i+1
            rest=rest/div
        enddo
        div=next
        next=next+2
    end do

    call mem%allocate(d,ctr,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
    d=0
    i = 1
    div=2
    next=3
    rest=n
    do while ( rest .ne. 1 )
        do while ( mod(rest, div) == 0 )
            d(i) = div
            i=i+1
            rest=rest/div
        enddo
        div=next
        next=next+2
    end do

    ! And sort it largest to smallest.
    d=-d
    call qsort(d)
    d=-d
end subroutine

end module
