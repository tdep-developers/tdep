#include "precompilerdefinitions"
!> wrap ScaLAPACK functions
module scalapack_wrapper
    use constants
    use gottochblandat
    implicit none

    private
    public :: lo_blacs_helper
    public :: lo_scalapack_matrix
    public :: lo_pdgels

    !> BLACS parameters, essentially the same as the MPI parameters.
    type lo_blacs_helper
        !> context, suppose it's like an MPI communicator
        integer :: icontxt=-lo_hugeint
        !> is it initialized?
        logical :: initialized=.false.
        !> number of processors
        integer :: n=-lo_hugeint
        !> whoamI?
        integer :: r=-lo_hugeint
        !> dimensions of processor grid
        integer :: nrow=-lo_hugeint,ncol=-lo_hugeint
        !> where in the grid am I?
        integer :: row=-lo_hugeint,col=-lo_hugeint
        contains
            procedure :: init=>init_blacs
            procedure :: destroy=>destroy_blacs
    end type

    !> A distributed matrix, the way ScaLAPACK likes it
    type lo_scalapack_matrix
        !> size of local matrix
        integer :: nc=-lo_hugeint,nr=-lo_hugeint
        !> size of global matrix
        integer :: ncol=-lo_hugeint,nrow=-lo_hugeint
        !> local matrix
        real(flyt), dimension(:,:), allocatable :: buf
        !> temporary buffer, for solutions and stuff
        real(flyt), dimension(:,:), allocatable :: tmpbuf
        !> blocksize in each direction
        integer :: rowblocksize=-lo_hugeint
        integer :: colblocksize=-lo_hugeint
        !> where the processors start indexing. Should always be zero, I think
        integer :: proc_col_zero=-lo_hugeint,proc_row_zero=-lo_hugeint
        !> Scalapack descriptors that I don't understand
        integer, dimension(9) :: desc=-lo_hugeint
        contains
            !> setup basic things for this matrix
            procedure :: init=>init_scw
            !> take local slices and rearrange it block-cyclic
            procedure :: reshuffle_matrix
    end type

contains

!> Solve normal linear system: arg min ||Ax-B||
subroutine lo_pdgels(mA,mB,solution,bw,mw,nodestroy)
    !> ScaLAPACK-friendly matrix A
    type(lo_scalapack_matrix), intent(inout) :: mA
    !> ScaLAPACK-friendly matrix B
    type(lo_scalapack_matrix), intent(inout) :: mB
    !> Least squares solution
    real(flyt), dimension(:), intent(out) :: solution
    !> BLACS helper
    type(lo_blacs_helper), intent(inout) :: bw
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> do I need the matrices afterwards?
    logical, intent(in), optional :: nodestroy
    
    logical :: keepmatrix
    real(flyt), dimension(:), allocatable :: solbuf,work
    real(flyt), dimension(1) :: dwrk
    character(len=1) :: trans
    integer :: lwork,i,j,ii,jj,nsol
    
    ! Some simple heuristics
    if ( present(nodestroy) ) then
        keepmatrix=nodestroy
    else
        keepmatrix=.false.
    endif
    ! Some other things I assume
    if ( mB%ncol .ne. 1 ) then
        write(*,*) 'I need to update the routine to handle several RHS'
        call bw%destroy()
        call mw%destroy()
        stop
    endif
    if ( mA%nrow .ne. mB%nrow ) then
        write(*,*) 'Dimensions of distributed matrices do not match'
        call bw%destroy()
        call mw%destroy()
        stop
    endif

    ! Make space for the solution buffer
    nsol=size(solution,1)
    lo_allocate(solbuf(nsol))

    ! Actually solve stuff
    if ( keepmatrix ) then
        ! solve the system, and keep a copy of the matrices
        lo_allocate(mA%tmpbuf(size(mA%buf,1),size(mA%buf,2)))
        lo_allocate(mB%tmpbuf(size(mB%buf,1),size(mB%buf,2)))
        mA%tmpbuf=mA%buf
        mB%tmpbuf=mB%buf
        ! query for workspace
        lwork=-1
        dwrk=0.0_flyt
        call pdgels('N',mA%nrow,mA%ncol,mB%ncol,mA%tmpbuf,1,1,mA%desc,mB%tmpbuf,1,1,mB%desc, dwrk, lwork, lo_status)
        lwork=int(anint(dwrk(1)))
        lo_allocate(work(lwork))
        work=0.0_flyt
        ! actual solution
        call pdgels('N',mA%nrow,mA%ncol,mB%ncol,mA%tmpbuf,1,1,mA%desc,mB%tmpbuf,1,1,mB%desc, work, lwork, lo_status)
        ! store partial solution in buffer
        solbuf=0.0_flyt
        do i=1,mB%nr
        do j=1,mB%nc
            ii=global_from_local(i,bw%row,mB%nrow,bw%nrow,mB%rowblocksize)
            jj=global_from_local(j,bw%col,mB%ncol,bw%ncol,mB%colblocksize)
            if ( ii .le. nsol ) solbuf(ii)=mB%tmpbuf(i,j)
        enddo
        enddo
        lo_deallocate(mA%tmpbuf)
        lo_deallocate(mB%tmpbuf)
    else
        ! query for workspace
        lwork=-1
        dwrk=0.0_flyt
        call pdgels('N',mA%nrow,mA%ncol,mB%ncol,mA%buf,1,1,mA%desc,mB%buf,1,1,mB%desc, dwrk, lwork, lo_status)
        lwork=int(anint(dwrk(1))) 
        lo_allocate(work(lwork))
        work=0.0_flyt
        ! actual solution
        call pdgels('N',mA%nrow,mA%ncol,mB%ncol,mA%buf,1,1,mA%desc,mB%buf,1,1,mB%desc, work, lwork, lo_status)
        ! store partial solution in buffer
        solbuf=0.0_flyt
        do i=1,mB%nr
        do j=1,mB%nc
            ii=global_from_local(i,bw%row,mB%nrow,bw%nrow,mB%rowblocksize)
            jj=global_from_local(j,bw%col,mB%ncol,bw%ncol,mB%colblocksize)
            if ( ii .le. nsol ) solbuf(ii)=mB%buf(i,j)
        enddo
        enddo
    endif

    ! Add up the solution over ranks
    solution=0.0_flyt
    call mpi_allreduce(solbuf,solution,nsol,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    
    ! Some cleanup
    if ( allocated(mA%tmpbuf) ) lo_deallocate(mA%tmpbuf)
    if ( allocated(mb%tmpbuf) ) lo_deallocate(mB%tmpbuf)
    if ( allocated(solbuf) ) lo_deallocate(solbuf)
    if ( allocated(work) ) lo_deallocate(work)
end subroutine

!> initialize BLACS
subroutine init_blacs(bw)
    !> blacs helper
    class(lo_blacs_helper), intent(out) :: bw

    integer, dimension(2) :: dims
    ! basic info
    call blacs_pinfo(bw%r,bw%n)
    ! get neat dimensions of processor grid
    dims=0
    call MPI_Dims_create(bw%n,2,dims,lo_status)
    if ( lo_status .ne. 0 ) then
        write(*,*) 'MPI_Dims_create status:',lo_status
    endif
    bw%nrow=dims(1)
    bw%ncol=dims(2)
    call blacs_get(0,0,bw%icontxt,lo_status)
    call blacs_gridinit(bw%icontxt,'R',bw%nrow,bw%ncol,lo_status)
    call blacs_gridinfo(bw%icontxt,bw%nrow,bw%ncol,bw%row,bw%col,lo_status)
!    write(*,*) 'rank ',tochar(bw%r),' nproc ',tochar(bw%n),' coord ',tochar([bw%row,bw%col]),bw%row*bw%ncol+bw%col
    bw%initialized=.true.
end subroutine

!> destroy BLACS
subroutine destroy_blacs(bw)
    !> blacs helper
    class(lo_blacs_helper), intent(out) :: bw
    ! the "1" means that I should be able to do normal MPI after this.
    call blacs_exit(1)
end subroutine

!> set up scalapack thingies for a large matrix
subroutine init_scw(scw,nrow,ncol,bw)
    !> ScaLAPACK matrix
    class(lo_scalapack_matrix), intent(out) :: scw
    !> number of columns in the global matrix
    integer, intent(in) :: ncol
    !> number of rows in the global matrix
    integer, intent(in) :: nrow
    !> BLACS helper
    type(lo_blacs_helper), intent(in) :: bw

    integer, parameter :: defblock=10 !20 
    integer :: i,j
    ! Store some basic info
    scw%nrow=nrow
    scw%ncol=ncol
    ! start all indexing at zero
    scw%proc_col_zero=0 
    scw%proc_row_zero=0 
    ! set a blocksize
    scw%colblocksize=min(defblock,ncol)
    scw%rowblocksize=min(defblock,nrow)
    ! size of the local buffer
    i=numrc(scw%nrow,scw%rowblocksize,bw%row,scw%proc_row_zero,bw%nrow)
    j=numrc(scw%ncol,scw%colblocksize,bw%col,scw%proc_col_zero,bw%ncol)
    ! fix the descriptor thingy
    call descinit(scw%desc,scw%nrow,scw%ncol,scw%rowblocksize,scw%colblocksize,&
                  scw%proc_row_zero,scw%proc_col_zero,bw%icontxt,i,lo_status)
end subroutine

!> Given a slice of matrix A, rearrange it so that ScaLAPACK likes it
subroutine reshuffle_matrix(scw,buf,rows,mw,bw,verbosity)
    !> ScaLAPACK matrix
    class(lo_scalapack_matrix), intent(inout) :: scw
    !> the input slice
    real(flyt), dimension(:,:), intent(in) :: buf
    !> which columns did I get?
    integer, dimension(:), intent(in) :: rows
    !> mpi helper
    type(lo_mpi_helper), intent(in) :: mw
    !> blacs helper
    type(lo_blacs_helper), intent(in) :: bw 
    !> how much to talk
    integer, intent(in), optional :: verbosity

    real(flyt), dimension(:), allocatable :: dsendbuf,drecvbuf
    real(flyt) :: t0,t1,t2
    integer, dimension(:,:), allocatable :: isendbuf,irecvbuf
    integer, dimension(mw%n,mw%n) :: local_ctr,global_ctr
    integer, dimension(mw%n) :: rankctr,rankoffset
    integer :: i,j,k,l,nr
    integer :: li,lj,lr,pcol,prow
    integer :: destrank,recvrank,nrecv
    logical :: talk

    talk=.false.
    if ( present(verbosity) ) then
        if ( verbosity .gt. 0 ) talk=.true.
    endif

    nr=size(rows,1)                               ! number of rows on this rank
    if ( allocated(scw%buf) ) deallocate(scw%buf) ! reset the buffer, if necessary

    ! Start by counting how much data each rank should recieve from each other rank?
    local_ctr=0
    do lr=1,nr
        i=rows(lr)
        do j=1,scw%ncol
            ! Where should this data go?
            prow=proc_from_global(i,bw%nrow,scw%rowblocksize)
            pcol=proc_from_global(j,bw%ncol,scw%colblocksize)
            destrank=prow*bw%ncol+pcol+1 
            local_ctr(destrank,mw%r+1)=local_ctr(destrank,mw%r+1)+1
        enddo
    enddo
    ! Sum up the counters
    global_ctr=0
    call mpi_allreduce(local_ctr,global_ctr,mw%n*mw%n,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)

    if ( talk ) write(*,*) '... estimated data to be sent'

    ! Now do a series of gatherv to get everything to the right place.
    do recvrank=1,mw%n
        t0=walltime()
        ! how many things should this rank get from all others?
        nrecv=sum(global_ctr(recvrank,:))
        ! Count number of things to send per rank, and their offset.
        rankctr=0
        rankoffset=0
        do i=1,mw%n
            rankctr(i)=global_ctr(recvrank,i)
        enddo
        j=0
        do i=1,mw%n-1
            j=j+rankctr(i)
            rankoffset(i+1)=j
        enddo
        ! Fetch what to send
        if ( rankctr(mw%r+1) .gt. 0 ) then
            lo_allocate(isendbuf(2,rankctr(mw%r+1)))
            lo_allocate(dsendbuf(rankctr(mw%r+1)))
            isendbuf=0
            dsendbuf=0.0_flyt
            l=0
                do j=1,scw%ncol
            do lr=1,nr
                i=rows(lr)
                    prow=proc_from_global(i,bw%nrow,scw%rowblocksize)
                    pcol=proc_from_global(j,bw%ncol,scw%colblocksize)
                    destrank=prow*bw%ncol+pcol+1
                    if ( destrank .eq. recvrank ) then
                        l=l+1
                        isendbuf(:,l)=[i,j]
                        dsendbuf(l)=buf(lr,j)
                    endif
                enddo
            enddo
        else
            lo_allocate(isendbuf(1,1))
            lo_allocate(dsendbuf(1))
            isendbuf=0
            dsendbuf=0
        endif
        ! some space to recieve things
        lo_allocate(drecvbuf(nrecv))
        lo_allocate(irecvbuf(2,nrecv))
        drecvbuf=0.0_flyt
        irecvbuf=0        
        ! Send the indices
        if ( talk ) write(*,*) '... packed data for rank ',tochar(recvrank-1),' (',tochar(walltime()-t0),'s)'
        t0=walltime()
        call mpi_gatherv(isendbuf,rankctr(mw%r+1)*2,MPI_INTEGER,irecvbuf,rankctr*2,rankoffset*2,&
                         MPI_INTEGER,recvrank-1,mw%comm,mw%error)
        ! Send the values
        call mpi_gatherv(dsendbuf,rankctr(mw%r+1),MPI_DOUBLE_PRECISION,drecvbuf,rankctr,rankoffset,&
                         MPI_DOUBLE_PRECISION,recvrank-1,mw%comm,mw%error)
        if ( talk ) write(*,*) '... sent data to rank ',tochar(recvrank-1),' (',tochar(walltime()-t0),'s)'
        ! see what I got, and store it in the local buffer
        if ( mw%r+1 .eq. recvrank ) then
            scw%nc=0
            scw%nr=0
            do lr=1,size(irecvbuf,2)
                li=local_from_global(irecvbuf(1,lr),scw%nrow,bw%nrow,scw%rowblocksize)
                lj=local_from_global(irecvbuf(2,lr),scw%ncol,bw%ncol,scw%colblocksize)
                scw%nr=max(scw%nr,li)
                scw%nc=max(scw%nc,lj)
            enddo
            if ( scw%nr .gt. 0 .and. scw%nc .gt. 0 ) then
                lo_allocate(scw%buf(scw%nr,scw%nc))
                scw%buf=0.0_flyt
            endif
            do lr=1,size(irecvbuf,2)
                li=local_from_global(irecvbuf(1,lr),scw%nrow,bw%nrow,scw%rowblocksize)
                lj=local_from_global(irecvbuf(2,lr),scw%ncol,bw%ncol,scw%colblocksize)
                scw%buf(li,lj)=drecvbuf(lr)
            enddo
        endif
        ! cleanup
        if ( allocated(isendbuf) ) deallocate(isendbuf)
        if ( allocated(dsendbuf) ) deallocate(dsendbuf)
        if ( allocated(drecvbuf) ) deallocate(drecvbuf)
        if ( allocated(irecvbuf) ) deallocate(irecvbuf)
    enddo

    ! maybe this rank has an empty buffer, try to fix that
    if ( allocated(scw%buf) .eqv. .false. ) then
        lo_allocate(scw%buf(1,1))
        scw%buf=0.0_flyt
        scw%nr=0
        scw%nc=0
    endif
    ! And it should be ok!
end subroutine

!> convert global index to local index in block-cyclic distribution
elemental function local_from_global(i,n,np,nb) result(il)
   integer, intent(in)  :: i    ! global array index, input
   integer, intent(in)  :: n    ! global array dimension, input
   integer, intent(in)  :: np   ! processor array dimension, input
   integer, intent(in)  :: nb   ! block size, input
   integer :: il                ! local array index, output

   integer :: im1 
   im1 = i-1
   il  = (im1/(np*nb))*nb + mod(im1,nb) + 1
end function

!> convert local index to global index in block-cyclic distribution
elemental function global_from_local(il,p,n,np,nb) result(i)
   integer, intent(in) :: il   ! local array index, input
   integer, intent(in) :: p    ! processor array index, input
   integer, intent(in) :: n    ! global array dimension, input
   integer, intent(in) :: np   ! processor array dimension, input
   integer, intent(in) :: nb   ! block size, input
   integer :: i                ! global array index, output

   integer :: ilm1   
   ilm1 = il-1
   i    = (((ilm1/nb) * np) + p)*nb + mod(ilm1,nb) + 1
end function

!> processor index from global
elemental function proc_from_global(i,np,nb) result(p)
   integer, intent(in)  :: i    ! global array index, input
   integer, intent(in)  :: np   ! processor array dimension, input
   integer, intent(in)  :: nb   ! block size, input
   integer  :: p                ! processor array index, output

   integer :: im1 
   im1 = i-1
   p   = mod((im1/nb),np)
end function

!> stolen from the ScaLAPACK reference
integer function numrc( n, nb, iproc, isrcproc, nprocs )
    integer, intent(in) :: iproc, isrcproc, n, nb, nprocs
    ! Purpose
    ! =======
    ! NUMROC computes the NUMber of Rows Or Columns of a distributed
    ! matrix owned by the process indicated by IPROC.
    ! Arguments
    ! =========
    ! N         (global input) INTEGER
    !           The number of rows/columns in distributed matrix.
    ! NB        (global input) INTEGER
    !           Block size, size of the blocks the distributed matrix is
    !           split into.
    ! IPROC     (local input) INTEGER
    !           The coordinate of the process whose local array row or
    !           column is to be determined.
    ! ISRCPROC  (global input) INTEGER
    !           The coordinate of the process that possesses the first
    !           row or column of the distributed matrix.
    ! NPROCS    (global input) INTEGER
    !           The total number processes over which the matrix is
    !           distributed.
    ! =====================================================================
    integer :: extrablks, mydist, nblocks
    ! Figure PROC's distance from source process
    mydist = mod( nprocs+iproc-isrcproc, nprocs )
    ! Figure the total number of whole NB blocks N is split up into
    nblocks = n / nb
    ! Figure the minimum number of rows/cols a process can have
    numrc = (nblocks/nprocs) * nb
    ! See if there are any extra blocks
    extrablks = mod( nblocks, nprocs )
    ! If I have an extra block
    if( mydist .lt. extrablks ) then
        numrc = numrc + nb
        ! If I have last block, it may be a partial block
    elseif( mydist .eq. extrablks ) then
        numrc = numrc + mod( n, nb )
    end if
end function

end module

