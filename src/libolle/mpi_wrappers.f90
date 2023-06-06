#include "precompilerdefinitions"
module mpi_wrappers
!!
!! I have wrappers for all the most common MPI procedures. Mostly for catching errors
!! and making the code less verbose. Also I never remember the correct syntax.
!!
use konstanter, only: r8,i8,lo_status,lo_hugeint,lo_exitcode_mpi
use gottochblandat, only: tochar
use mpi

implicit none
private

public :: lo_mpi_helper
public :: lo_stop_gracefully
! Expose some MPI constants, in case someone needs them
public :: MPI_CHARACTER,MPI_DOUBLE_COMPLEX,MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_LOGICAL
public :: MPI_MAX,MPI_MIN,MPI_SUM,MPI_IN_PLACE,MPI_ANY_TAG

!> Helper that keeps track of all things MPI.
type lo_mpi_helper
    !> which communicator
    integer :: comm=-lo_hugeint
    !> current rank
    integer :: r=-lo_hugeint
    !> total number of ranks
    integer :: n=-lo_hugeint
    !> exit status
    integer :: error=-lo_hugeint
    !> is this rank allowed to write to stdout?
    logical :: talk=.false.
    !> current status
    integer, dimension(MPI_STATUS_SIZE) :: stat
    !> counter for size per rank
    integer, dimension(:), allocatable :: ctr_size_per_rank
    !> counter for offset per rank
    integer, dimension(:), allocatable :: ctr_offset_per_rank
    contains
        !> initialize
        procedure :: init=>init_mpi
        !> close
        procedure :: destroy=>kill_mpi
        !> split communicator
        procedure :: split=>split_communicator
        !> barrier
        procedure :: barrier
        !> free a communicator
        procedure :: free
        !> get size and offset of distributed array
        procedure :: size_and_offset
        !> wait for nonblocking communication
        procedure :: wait

        !> mpi allreduce
        generic :: allreduce=>allreduce_int,allreduce_real,allreduce_1d_real,allreduce_2d_real,allreduce_3d_real,&
                              allreduce_4d_real,allreduce_5d_real,allreduce_6d_real,allreduce_1d_int,allreduce_2d_int,allreduce_3d_int,&
                              allreduce_2d_complex,allreduce_3d_complex,allreduce_4d_complex,allreduce_5d_complex,allreduce_1d_logical,&
                              allreduce_int8
        procedure, private :: allreduce_int
        procedure, private :: allreduce_1d_int
        procedure, private :: allreduce_2d_int
        procedure, private :: allreduce_3d_int
        procedure, private :: allreduce_real
        procedure, private :: allreduce_1d_real
        procedure, private :: allreduce_2d_real
        procedure, private :: allreduce_3d_real
        procedure, private :: allreduce_4d_real
        procedure, private :: allreduce_5d_real
        procedure, private :: allreduce_6d_real
        procedure, private :: allreduce_2d_complex
        procedure, private :: allreduce_3d_complex
        procedure, private :: allreduce_4d_complex
        procedure, private :: allreduce_5d_complex
        procedure, private :: allreduce_1d_logical
        procedure, private :: allreduce_int8
        !> mpi broadcast
        generic :: bcast=>broadcast_int,broadcast_1d_int,broadcast_2d_int,broadcast_real,broadcast_1d_real,&
                          broadcast_2d_real,broadcast_3d_real,broadcast_4d_real,broadcast_4d_complex,&
                          broadcast_character,broadcast_logical
        procedure, private :: broadcast_int
        procedure, private :: broadcast_1d_int
        procedure, private :: broadcast_2d_int
        procedure, private :: broadcast_real
        procedure, private :: broadcast_1d_real
        procedure, private :: broadcast_2d_real
        procedure, private :: broadcast_3d_real
        procedure, private :: broadcast_4d_real
        procedure, private :: broadcast_4d_complex
        procedure, private :: broadcast_character
        procedure, private :: broadcast_logical

        !>@TODO add tolerance for check and sync
        !> check and synchronize variables
        generic :: check_and_sync=>check_and_sync_int,check_and_sync_1d_int,check_and_sync_real,check_and_sync_1d_real,&
                                   check_and_sync_2d_real
        procedure, private :: check_and_sync_int
        procedure, private :: check_and_sync_1d_int
        procedure, private :: check_and_sync_real
        procedure, private :: check_and_sync_1d_real
        procedure, private :: check_and_sync_2d_real

        !> mpi reduce
        generic :: reduce=>reduce_2d_real,reduce_3d_real,reduce_4d_real,reduce_5d_real,reduce_4d_complex
        procedure, private :: reduce_2d_real
        procedure, private :: reduce_3d_real
        procedure, private :: reduce_4d_real
        procedure, private :: reduce_5d_real
        procedure, private :: reduce_4d_complex

        !> mpi pack
        generic :: pack=>pack_integer,pack_1d_integer,pack_2d_integer,&
        pack_real,pack_1d_real,pack_2d_real,pack_3d_real,pack_4d_real,pack_5d_real,&
        pack_complex,pack_1d_complex,pack_2d_complex,pack_3d_complex,&
        pack_1d_logical
        procedure, private :: pack_integer
        procedure, private :: pack_1d_integer
        procedure, private :: pack_2d_integer
        procedure, private :: pack_real
        procedure, private :: pack_1d_real
        procedure, private :: pack_2d_real
        procedure, private :: pack_3d_real
        procedure, private :: pack_4d_real
        procedure, private :: pack_5d_real
        procedure, private :: pack_complex
        procedure, private :: pack_1d_complex
        procedure, private :: pack_2d_complex
        procedure, private :: pack_3d_complex
        procedure, private :: pack_1d_logical

        !> mpi unpack
        generic :: unpack=>unpack_integer,unpack_1d_integer,unpack_2d_integer,&
        unpack_real,unpack_1d_real,unpack_2d_real,unpack_3d_real,unpack_4d_real,unpack_5d_real,&
        unpack_complex,unpack_1d_complex,unpack_2d_complex,unpack_3d_complex,&
        unpack_1d_logical
        procedure, private :: unpack_integer
        procedure, private :: unpack_1d_integer
        procedure, private :: unpack_2d_integer
        procedure, private :: unpack_real
        procedure, private :: unpack_1d_real
        procedure, private :: unpack_2d_real
        procedure, private :: unpack_3d_real
        procedure, private :: unpack_4d_real
        procedure, private :: unpack_5d_real
        procedure, private :: unpack_complex
        procedure, private :: unpack_1d_complex
        procedure, private :: unpack_2d_complex
        procedure, private :: unpack_3d_complex
        procedure, private :: unpack_1d_logical

        !> mpi allgatherv
        generic :: gatherv=>gatherv_1d_integer
        procedure, private :: gatherv_1d_integer

        !> mpi allgatherv
        generic :: allgatherv=>allgatherv_character,allgatherv_1d_integer
        procedure, private :: allgatherv_character,allgatherv_1d_integer

        !> mpi alltoallv
        generic :: alltoallv=>alltoallv_character
        procedure, private :: alltoallv_character

        !> mpi scatterv
        generic :: scatterv=>scatterv_character
        procedure, private :: scatterv_character

        !> mpi send
        generic :: send=>send_integer,send_1d_integer,send_2d_integer,send_1d_real,send_2d_real,send_character
        procedure, private :: send_integer
        procedure, private :: send_1d_integer
        procedure, private :: send_2d_integer
        procedure, private :: send_1d_real
        procedure, private :: send_2d_real
        procedure, private :: send_character

        !> mpi recv
        generic :: recv=>recv_integer,recv_1d_integer,recv_2d_integer,recv_1d_real,recv_2d_real,recv_character
        procedure, private :: recv_integer
        procedure, private :: recv_1d_integer
        procedure, private :: recv_2d_integer
        procedure, private :: recv_1d_real
        procedure, private :: recv_2d_real
        procedure, private :: recv_character

        !> mpi isend
        generic :: isend=>isend_character
        procedure, private :: isend_character

        !> mpi_irecv
        generic :: irecv=>irecv_character
        procedure, private :: irecv_character
end type

contains

!> initialize MPI
subroutine init_mpi(mw,communicator,notalk)
    !> communicator helper
    class(lo_mpi_helper), intent(out) :: mw
    !> specific communicator, if not world
    integer, intent(in), optional :: communicator
    !> are all ranks forbidden to talk?
    logical, intent(in), optional :: notalk

    logical :: init

    ! start MPI, if not already done so
    call mpi_initialized(init,mw%error)
    if ( init .eqv. .false. ) call mpi_init(mw%error)

    ! communicator
    if ( present(communicator) ) then
        mw%comm=communicator
    else
        mw%comm=MPI_COMM_WORLD
    endif

    ! current rank
    call mpi_comm_rank(mw%comm,mw%r,mw%error)
    ! number of ranks
    call mpi_comm_size(mw%comm,mw%n,mw%error)
    ! which rank can talk? Always rank 0.
    mw%talk=.false.
    if ( mw%r .eq. 0 ) mw%talk=.true.
    if ( present(notalk) ) then
        if ( notalk ) then
            mw%talk=.false.
        endif
    endif

    ! space for counters, might as well keep it here
    allocate(mw%ctr_size_per_rank(mw%n))
    allocate(mw%ctr_offset_per_rank(mw%n))
    mw%ctr_size_per_rank=0
    mw%ctr_offset_per_rank=0

    ! status
    mw%stat=0
    mw%error=0
end subroutine

!> get amount of data per rank, and the offset if it is to be collected
subroutine size_and_offset(mw,size_on_self,total_size)
    !> communicator helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what is the size on this rank
    integer, intent(in) :: size_on_self
    !> what is the total size
    integer, intent(out) :: total_size

    integer :: i,j

    mw%ctr_size_per_rank=0
    mw%ctr_offset_per_rank=0
    mw%ctr_size_per_rank(mw%r+1)=size_on_self

    call mpi_allreduce(MPI_IN_PLACE,mw%ctr_size_per_rank,mw%n,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
    if ( mw%error .ne. 0 ) then
        call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
    endif

    j=0
    do i=1,mw%n-1
        j=j+mw%ctr_size_per_rank(i)
        mw%ctr_offset_per_rank(i+1)=j
    enddo
    total_size=sum(mw%ctr_size_per_rank)
end subroutine

!> free communicator
subroutine free(mw,filename,linenumber)
    !> communicator helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> file name from where this is called
    character(len=*), intent(in), optional :: filename
    !> line number where this was called
    integer, intent(in), optional :: linenumber

    call mpi_comm_free(mw%comm,mw%error)
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_comm_free exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_comm_free exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> barrier
subroutine barrier(mw,filename,linenumber)
    !> communicator helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> file name from where this is called
    character(len=*), intent(in), optional :: filename
    !> line number where this was called
    integer, intent(in), optional :: linenumber

    call mpi_barrier(mw%comm,mw%error)
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_barrier exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_barrier exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> kill MPI
subroutine kill_mpi(mw)
    !> communicator helper
    class(lo_mpi_helper), intent(inout) :: mw

    ! just make sure everything is destroyed
    call mpi_finalize(mw%error)
    mw%comm=-1
    mw%r=-1
    mw%n=-1
    mw%error=-1
    mw%stat=-1
end subroutine

!> split an MPI communicator into chunks
subroutine split_communicator(mw,mn,color,filename,linenumber)
    !> communicator to split
    class(lo_mpi_helper), intent(inout) :: mw
    !> split communicator
    type(lo_mpi_helper), intent(out) :: mn
    !> per-rank color, ranks with equal color end up in the same communicator
    integer, intent(in) :: color
    !> file name from where this is called
    character(len=*), intent(in), optional :: filename
    !> line number where this was called
    integer, intent(in), optional :: linenumber

    logical :: dlog

    ! Check that MPI is running
    call mpi_initialized(dlog,mw%error)
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['MPI error when checking initialization'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['MPI error when checking initialization'],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
    if ( dlog .eqv. .false. ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['No MPI initialized when splitting.'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['No MPI initialized when splitting.'],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif

    ! No error checks are done on the color. Not sure if good.
    call MPI_comm_split(mw%comm,color,mw%r,mn%comm,mw%error)
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['MPI error when splitting.'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['MPI error when splitting.'],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif

    ! current rank
    call mpi_comm_rank(mn%comm,mn%r,mn%error)
    ! number of ranks
    call mpi_comm_size(mn%comm,mn%n,mn%error)
    ! which rank can talk? None of them.
    mn%talk=.false.

    ! space for counters, might as well keep it here
    allocate(mn%ctr_size_per_rank(mn%n))
    allocate(mn%ctr_offset_per_rank(mn%n))
    mn%ctr_size_per_rank=0
    mn%ctr_offset_per_rank=0

    ! status
    mn%stat=0
    mn%error=0
end subroutine

!> mpi wait for nonblocking communication
subroutine wait(mw,request,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> request
    integer, intent(inout) :: request
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    call mpi_wait(request,mw%stat,mw%error)
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_wait exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_wait exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> Will check that a variable is synced across ranks
subroutine check_and_sync_int(mw,x,ref,bigdeal,vname,filename,linenumber)
    !> mpi helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> int
    integer, intent(inout) :: x
    !> which rank to use as reference
    integer, intent(in) :: ref
    !> how big of an issue?
    integer, intent(in), optional :: bigdeal
    !> name of variable for message
    character(len=*), intent(in), optional :: vname
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    integer :: j
    integer :: loudness

    if ( present(bigdeal) ) then
        loudness=bigdeal
    else
        loudness=2 ! Default to angry error
    endif

    j=0
    call mpi_allreduce(x,j,1,MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
    j=j/mw%n
    if ( abs(x-j) .gt. 0 ) then
        select case(loudness)
        case(0)
            ! Just sync and move on with life
            call MPI_Bcast(x,1,MPI_INTEGER,ref,mw%comm,mw%error)
        case(1)
            ! Sync and throw warning
            call MPI_Bcast(x,1,MPI_INTEGER,ref,mw%comm,mw%error)
            if ( mw%talk ) then
                if ( present(vname) ) then
                    write(*,*) 'WARNING: "'//trim(vname)//'" out of sync:'
                else
                    write(*,*) 'WARNING: out of sync:'
                endif
                write(*,*) '    average:',j
                write(*,*) '     rank 0:',x
            endif
        case(2)
            ! KILL
            if ( present(vname) .and. present(filename) .and. present(linenumber) ) then
                call lo_stop_gracefully(['"'//trim(vname)//'" out of sync'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            else
                call lo_stop_gracefully(['out of sync'],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
            endif
        end select
    endif
end subroutine
subroutine check_and_sync_1d_int(mw,x,ref,bigdeal,vname,filename,linenumber)
    !> mpi helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> integer to check
    integer, dimension(:), intent(inout) :: x
    !> which rank to use as reference
    integer, intent(in) :: ref
    !> how big of an issue?
    integer, intent(in), optional :: bigdeal
    !> name of variable for message
    character(len=*), intent(in), optional :: vname
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    integer, dimension(:), allocatable :: dr
    integer :: loudness

    if ( present(bigdeal) ) then
        loudness=bigdeal
    else
        loudness=2
    endif
    lo_allocate(dr(size(x)))
    dr=0
    call mpi_allreduce(x,dr,size(x),MPI_INTEGER,MPI_SUM,mw%comm,mw%error)
    dr=dr/mw%n
    if ( sum(abs(x-dr)) .ne. 0 ) then
        select case(loudness)
        case(0)
            ! Just sync and move on with life
            call MPI_Bcast(x,size(x),MPI_INTEGER,ref,mw%comm,mw%error)
        case(1)
            ! Sync and throw warning
            call MPI_Bcast(x,size(x),MPI_INTEGER,ref,mw%comm,mw%error)
            if ( mw%talk ) then
                if ( present(vname) ) then
                    write(*,*) 'WARNING: "'//trim(vname)//'" out of sync:'
                else
                    write(*,*) 'WARNING: out of sync:'
                endif
            endif
        case(2)
            ! KILL
            if ( present(vname) ) then
                call lo_stop_gracefully(['"'//trim(vname)//'" out of sync'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            else
                call lo_stop_gracefully(['out of sync'],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
            endif
        end select
    endif
end subroutine
subroutine check_and_sync_real(mw,x,ref,bigdeal,vname,filename,linenumber)
    !> mpi helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> real to check
    real(r8), intent(inout) :: x
    !> which rank to use as reference
    integer, intent(in) :: ref
    !> how big of an issue?
    integer, intent(in), optional :: bigdeal
    !> name of variable for message
    character(len=*), intent(in), optional :: vname
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    real(r8) :: dr
    integer :: loudness

    if ( present(bigdeal) ) then
        loudness=bigdeal
    else
        loudness=2 ! Default to angry error
    endif

    dr=0.0_r8
    call mpi_allreduce(x,dr,1,MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    dr=dr/mw%n
    if ( abs(x-dr) .gt. 1E-10_r8 ) then
        select case(loudness)
        case(0)
            ! Just sync and move on with life
            call MPI_Bcast(x,1,MPI_DOUBLE_PRECISION,ref,mw%comm,mw%error)
        case(1)
            ! Sync and throw warning
            call MPI_Bcast(x,1,MPI_DOUBLE_PRECISION,ref,mw%comm,mw%error)
            if ( mw%talk ) then
                if ( present(vname) ) then
                    write(*,*) 'WARNING: "'//trim(vname)//'" out of sync:'
                else
                    write(*,*) 'WARNING: out of sync:'
                endif
                write(*,*) '    average:',dr
                write(*,*) '     rank 0:',x
            endif
        case(2)
            ! KILL
            if ( present(vname) .and. present(filename) .and. present(linenumber) ) then
                call lo_stop_gracefully(['"'//trim(vname)//'" out of sync'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            else
                call lo_stop_gracefully(['out of sync'],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
            endif
        end select
    endif
end subroutine
subroutine check_and_sync_1d_real(mw,x,tol,ref,bigdeal,vname,filename,linenumber)
    !> mpi helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> real to check
    real(r8), dimension(:), intent(inout) :: x
    !> tolerance for what is ok
    real(r8), intent(in) :: tol
    !> which rank to use as reference
    integer, intent(in) :: ref
    !> how big of an issue?
    integer, intent(in), optional :: bigdeal
    !> name of variable for message
    character(len=*), intent(in), optional :: vname
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    real(r8), dimension(:), allocatable :: dr
    integer :: loudness

    if ( present(bigdeal) ) then
        loudness=bigdeal
    else
        loudness=2
    endif
    lo_allocate(dr(size(x)))
    dr=0.0_r8
    call mpi_allreduce(x,dr,size(x),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    dr=dr/mw%n
    if ( sum(abs(x-dr))/size(x) .gt. tol ) then
        select case(loudness)
        case(0)
            ! Just sync and move on with life
            call MPI_Bcast(x,size(x),MPI_DOUBLE_PRECISION,ref,mw%comm,mw%error)
        case(1)
            ! Sync and throw warning
            call MPI_Bcast(x,size(x),MPI_DOUBLE_PRECISION,ref,mw%comm,mw%error)
            if ( mw%talk ) then
                if ( present(vname) ) then
                    write(*,*) 'WARNING: "'//trim(vname)//'" out of sync:'
                else
                    write(*,*) 'WARNING: out of sync:'
                endif
            endif
        case(2)
            ! KILL
            if ( present(vname) ) then
                call lo_stop_gracefully(['"'//trim(vname)//'" out of sync'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            else
                call lo_stop_gracefully(['out of sync'],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
            endif
        end select
    endif
end subroutine
subroutine check_and_sync_2d_real(mw,x,ref,bigdeal,vname,filename,linenumber)
    !> mpi helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> real to check
    real(r8), dimension(:,:), intent(inout) :: x
    !> which rank to use as reference
    integer, intent(in) :: ref
    !> how big of an issue?
    integer, intent(in), optional :: bigdeal
    !> name of variable for message
    character(len=*), intent(in), optional :: vname
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    real(r8), dimension(:,:), allocatable :: dr
    integer :: loudness

    if ( present(bigdeal) ) then
        loudness=bigdeal
    else
        loudness=2
    endif
    lo_allocate(dr(size(x,1),size(x,2)))
    dr=0.0_r8
    call mpi_allreduce(x,dr,size(x),MPI_DOUBLE_PRECISION,MPI_SUM,mw%comm,mw%error)
    dr=dr/mw%n
    if ( sum(abs(x-dr))/size(x) .gt. 1E-10_r8 ) then
        select case(loudness)
        case(0)
            ! Just sync and move on with life
            call MPI_Bcast(x,size(x),MPI_DOUBLE_PRECISION,ref,mw%comm,mw%error)
        case(1)
            ! Sync and throw warning
            call MPI_Bcast(x,size(x),MPI_DOUBLE_PRECISION,ref,mw%comm,mw%error)
            if ( mw%talk ) then
                if ( present(vname) ) then
                    write(*,*) 'WARNING: "'//trim(vname)//'" out of sync:'
                else
                    write(*,*) 'WARNING: out of sync.'
                endif
            endif
        case(2)
            ! KILL
            if ( present(vname) ) then
                call lo_stop_gracefully(['"'//trim(vname)//'" out of sync'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            else
                call lo_stop_gracefully(['out of sync'],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
            endif
        end select
    endif
end subroutine

!> Wrappers for mpi_bcast
subroutine broadcast_int(mw,i,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Integer to broadcast
    integer, intent(inout) :: i
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    call MPI_Bcast(i,1,MPI_INTEGER,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_1d_int(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Array to broadcast
    integer, dimension(:), intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber
    ! Broadcast it
    call MPI_Bcast(d,size(d),MPI_INTEGER,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_2d_int(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Array to broadcast
    integer, dimension(:,:), intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber
    ! Broadcast it
    call MPI_Bcast(d,size(d),MPI_INTEGER,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_real(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Integer to broadcast
    real(r8), intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    call MPI_Bcast(d,1,MPI_DOUBLE_PRECISION,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_1d_real(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Array to broadcast
    real(r8), dimension(:), intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber
    ! Broadcast it
    call MPI_Bcast(d,size(d),MPI_DOUBLE_PRECISION,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_2d_real(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Array to broadcast
    real(r8), dimension(:,:), intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber
    ! Broadcast it
    call MPI_Bcast(d,size(d),MPI_DOUBLE_PRECISION,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_3d_real(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Array to broadcast
    real(r8), dimension(:,:,:), intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber
    ! Broadcast it
    call MPI_Bcast(d,size(d),MPI_DOUBLE_PRECISION,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_4d_real(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Array to broadcast
    real(r8), dimension(:,:,:,:), intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber
    ! Broadcast it
    call MPI_Bcast(d,size(d),MPI_DOUBLE_PRECISION,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_4d_complex(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Array to broadcast
    complex(r8), dimension(:,:,:,:), intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber
    ! Broadcast it
    call MPI_Bcast(d,size(d),MPI_DOUBLE_PRECISION,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_character(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Array to broadcast
    character, dimension(:), intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber
    ! Broadcast it
    call MPI_Bcast(d,size(d),MPI_CHARACTER,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine broadcast_logical(mw,d,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> Array to broadcast
    logical, intent(inout) :: d
    !> Where to broadcast from
    integer, intent(in) :: from
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber
    ! Broadcast it
    call MPI_Bcast(d,1,MPI_LOGICAL,from,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_bcast exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> Wrappers for mpi_allreduce
subroutine allreduce_int(mw,operation,i,j,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> integer to allreduce
    integer, intent(inout) :: i
    !> destination to allgather to, if omitted default to in-place
    integer, intent(inout), optional :: j
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(j) ) then
        j=0
        call mpi_allreduce(i,j,1,MPI_INTEGER,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,i,1,MPI_INTEGER,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_1d_int(mw,operation,i,j,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> integer to allreduce
    integer, intent(inout), dimension(:) :: i
    !> destination to allgather to, if omitted default to in-place
    integer, intent(inout), dimension(:), optional :: j
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(j) ) then
        j=0
        call mpi_allreduce(i,j,size(i),MPI_INTEGER,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,i,size(i),MPI_INTEGER,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_2d_int(mw,operation,i,j,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> integer to allreduce
    integer, intent(inout), dimension(:,:) :: i
    !> destination to allgather to, if omitted default to in-place
    integer, intent(inout), dimension(:,:), optional :: j
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(j) ) then
        j=0
        call mpi_allreduce(i,j,size(i),MPI_INTEGER,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,i,size(i),MPI_INTEGER,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_3d_int(mw,operation,i,j,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> integer to allreduce
    integer, intent(inout), dimension(:,:,:) :: i
    !> destination to allgather to, if omitted default to in-place
    integer, intent(inout), dimension(:,:,:), optional :: j
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(j) ) then
        j=0
        call mpi_allreduce(i,j,size(i),MPI_INTEGER,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,i,size(i),MPI_INTEGER,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_real(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> number to allreduce
    real(r8), intent(inout) :: x
    !> destination to allreduce to, if omitted default to in-place
    real(r8), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        y=0.0_r8
        call mpi_allreduce(x,y,1,MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,1,MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_1d_real(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    real(r8), dimension(:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    real(r8), dimension(:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_2d_real(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    real(r8), dimension(:,:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    real(r8), dimension(:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_3d_real(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    real(r8), dimension(:,:,:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    real(r8), dimension(:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_4d_real(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    real(r8), dimension(:,:,:,:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    real(r8), dimension(:,:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_5d_real(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    real(r8), dimension(:,:,:,:,:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    real(r8), dimension(:,:,:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_6d_real(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    real(r8), dimension(:,:,:,:,:,:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    real(r8), dimension(:,:,:,:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_2d_complex(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    complex(r8), dimension(:,:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    complex(r8), dimension(:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_3d_complex(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    complex(r8), dimension(:,:,:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    complex(r8), dimension(:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_4d_complex(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    complex(r8), dimension(:,:,:,:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    complex(r8), dimension(:,:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_5d_complex(mw,operation,x,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to allreduce
    complex(r8), dimension(:,:,:,:,:), intent(inout) :: x
    !> array to allreduce to, if omitted default to in-place
    complex(r8), dimension(:,:,:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( size(x) .ne. size(y) ) then
            call lo_stop_gracefully(['mpi_allreduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
        endif
        y=0.0_r8
        call mpi_allreduce(x,y,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_COMPLEX,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_1d_logical(mw,operation,i,j,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> logical to allreduce
    logical, intent(inout), dimension(:) :: i
    !> destination to allgather to, if omitted default to in-place
    logical, intent(inout), dimension(:), optional :: j
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(j) ) then
        j=.false.
        call mpi_allreduce(i,j,size(i),MPI_LOGICAL,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,i,size(i),MPI_LOGICAL,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allreduce_int8(mw,operation,i,j,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> integer to allreduce
    integer(i8), intent(inout) :: i
    !> destination to allgather to, if omitted default to in-place
    integer(i8), intent(inout), optional :: j
    !> filename we call from for debugging
    character(len=*), intent(in), optional :: filename
    !> line number we call from for debugging
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(j) ) then
        j=0
        call mpi_allreduce(i,j,1,MPI_LONG,mpiop,mw%comm,mw%error)
    else
        call mpi_allreduce(MPI_IN_PLACE,i,1,MPI_LONG,mpiop,mw%comm,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allreduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> Wrappers for mpi_reduce
subroutine reduce_2d_real(mw,operation,x,recvrank,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to reduce
    real(r8), dimension(:,:), intent(inout) :: x
    !> rank to reduce to
    integer, intent(in) :: recvrank
    !> array to reduce to, if omitted default to in-place
    real(r8), dimension(:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( mw%r .eq. recvrank ) then
            if ( size(x) .ne. size(y) ) then
                call lo_stop_gracefully(['mpi_reduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            endif
            y=0.0_r8
            call mpi_reduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    else
        if ( mw%r .eq. recvrank ) then
            call mpi_reduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine reduce_3d_real(mw,operation,x,recvrank,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to reduce
    real(r8), dimension(:,:,:), intent(inout) :: x
    !> rank to reduce to
    integer, intent(in) :: recvrank
    !> array to reduce to, if omitted default to in-place
    real(r8), dimension(:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( mw%r .eq. recvrank ) then
            if ( size(x) .ne. size(y) ) then
                call lo_stop_gracefully(['mpi_reduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            endif
            y=0.0_r8
            call mpi_reduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    else
        if ( mw%r .eq. recvrank ) then
            call mpi_reduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine reduce_4d_real(mw,operation,x,recvrank,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to reduce
    real(r8), dimension(:,:,:,:), intent(inout) :: x
    !> rank to reduce to
    integer, intent(in) :: recvrank
    !> array to reduce to, if omitted default to in-place
    real(r8), dimension(:,:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( mw%r .eq. recvrank ) then
            if ( size(x) .ne. size(y) ) then
                call lo_stop_gracefully(['mpi_reduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            endif
            y=0.0_r8
            call mpi_reduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    else
        if ( mw%r .eq. recvrank ) then
            call mpi_reduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine reduce_5d_real(mw,operation,x,recvrank,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to reduce
    real(r8), dimension(:,:,:,:,:), intent(inout) :: x
    !> rank to reduce to
    integer, intent(in) :: recvrank
    !> array to reduce to, if omitted default to in-place
    real(r8), dimension(:,:,:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( mw%r .eq. recvrank ) then
            if ( size(x) .ne. size(y) ) then
                call lo_stop_gracefully(['mpi_reduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            endif
            y=0.0_r8
            call mpi_reduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    else
        if ( mw%r .eq. recvrank ) then
            call mpi_reduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine reduce_4d_complex(mw,operation,x,recvrank,y,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> what operation to do
    character(len=*), intent(in) :: operation
    !> array to reduce
    complex(r8), dimension(:,:,:,:), intent(inout) :: x
    !> rank to reduce to
    integer, intent(in) :: recvrank
    !> array to reduce to, if omitted default to in-place
    complex(r8), dimension(:,:,:,:), intent(out), optional :: y
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    integer :: mpiop
    ! Fetch proper operation code
    mpiop=mpi_operation_code(operation)
    ! Do the actual communication
    if ( present(y) ) then
        if ( mw%r .eq. recvrank ) then
            if ( size(x) .ne. size(y) ) then
                call lo_stop_gracefully(['mpi_reduce inconsistent array sizes'],lo_exitcode_mpi,filename,linenumber,mw%comm)
            endif
            y=0.0_r8
            call mpi_reduce(x,y,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    else
        if ( mw%r .eq. recvrank ) then
            call mpi_reduce(MPI_IN_PLACE,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        else
            call mpi_reduce(x,x,size(x),MPI_DOUBLE_PRECISION,mpiop,recvrank,mw%comm,mw%error)
        endif
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_reduce exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> Wrappers for mpi_pack
subroutine pack_integer(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    integer, intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,1,MPI_INTEGER,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_1d_integer(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    integer, dimension(:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_INTEGER,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_2d_integer(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    integer, dimension(:,:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_INTEGER,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,1,MPI_DOUBLE_PRECISION,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_1d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_DOUBLE_PRECISION,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_2d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:,:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_DOUBLE_PRECISION,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_3d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:,:,:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_DOUBLE_PRECISION,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_4d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:,:,:,:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_DOUBLE_PRECISION,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_5d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:,:,:,:,:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_DOUBLE_PRECISION,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_complex(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    complex(r8), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,1,MPI_DOUBLE_COMPLEX,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_1d_complex(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    complex(r8), dimension(:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_DOUBLE_COMPLEX,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_2d_complex(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    complex(r8), dimension(:,:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_DOUBLE_COMPLEX,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_3d_complex(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    complex(r8), dimension(:,:,:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_DOUBLE_COMPLEX,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine pack_1d_logical(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    logical, dimension(:), intent(in) :: x
    !> buffer to pack to
    character, dimension(:), intent(inout) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! pack it
    call mpi_pack(x,size(x),MPI_LOGICAL,buf,size(buf),pos,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_pack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> Wrappers for mpi_unpack
subroutine unpack_integer(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    integer, intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,1,MPI_INTEGER,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_1d_integer(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    integer, dimension(:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_INTEGER,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_2d_integer(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    integer, dimension(:,:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_INTEGER,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,1,MPI_DOUBLE_PRECISION,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_1d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_DOUBLE_PRECISION,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_2d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:,:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_DOUBLE_PRECISION,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_3d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:,:,:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_DOUBLE_PRECISION,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_4d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:,:,:,:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_DOUBLE_PRECISION,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_5d_real(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    real(r8), dimension(:,:,:,:,:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_DOUBLE_PRECISION,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_complex(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    complex(r8), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,1,MPI_DOUBLE_COMPLEX,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_1d_complex(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    complex(r8), dimension(:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_DOUBLE_COMPLEX,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_2d_complex(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    complex(r8), dimension(:,:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_DOUBLE_COMPLEX,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_3d_complex(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    complex(r8), dimension(:,:,:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_DOUBLE_COMPLEX,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine unpack_1d_logical(mw,x,buf,pos,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> thing to pack
    logical, dimension(:), intent(out) :: x
    !> buffer to unpack from
    character, dimension(:), intent(in) :: buf
    !> current position in buffer
    integer, intent(inout) :: pos
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! unpack it
    call mpi_unpack(buf,size(buf),pos,x,size(x),MPI_LOGICAL,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_unpack exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> wrappers for mpi allgatherv
subroutine gatherv_1d_integer(mw,sbuf,rbuf,ctr,offset,recvrank,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    integer, dimension(:), intent(in) :: sbuf
    !> buffer to recieve in
    integer, dimension(:), intent(out) :: rbuf
    !> number of things to receive from each rank
    integer, dimension(:), intent(in) :: ctr
    !> offset
    integer, dimension(:), intent(in) :: offset
    !> rank to gather to
    integer, intent(in) :: recvrank
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_gatherv(sbuf,ctr(mw%r+1),MPI_INTEGER,rbuf,ctr,offset,MPI_CHARACTER,recvrank,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_gatherv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_gatherv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> wrappers for mpi allgatherv
subroutine allgatherv_character(mw,sbuf,rbuf,ctr,offset,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    character, dimension(:), intent(in) :: sbuf
    !> buffer to recieve in
    character, dimension(:), intent(out) :: rbuf
    !> number of things to receive from each rank
    integer, dimension(:), intent(in) :: ctr
    !> offset
    integer, dimension(:), intent(in) :: offset
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_allgatherv(sbuf,ctr(mw%r+1),MPI_CHARACTER,rbuf,ctr,offset,MPI_CHARACTER,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allgatherv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allgatherv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine allgatherv_1d_integer(mw,sbuf,rbuf,ctr,offset,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    integer, dimension(:), intent(in) :: sbuf
    !> buffer to recieve in
    integer, dimension(:), intent(out) :: rbuf
    !> number of things to receive from each rank
    integer, dimension(:), intent(in) :: ctr
    !> offset
    integer, dimension(:), intent(in) :: offset
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_allgatherv(sbuf,ctr(mw%r+1),MPI_INTEGER,rbuf,ctr,offset,MPI_INTEGER,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allgatherv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allgatherv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> wrappers for mpi alltoallv
subroutine alltoallv_character(mw,sbuf,rbuf,sendcount,sendoffset,recvcount,recvoffset,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    character, dimension(:), intent(in) :: sbuf
    !> buffer to recieve in
    character, dimension(:), intent(out) :: rbuf
    !> number of things to send to each rank
    integer, dimension(:), intent(in) :: sendcount
    !> offset of send things
    integer, dimension(:), intent(in) :: sendoffset
    !> number of things to recieve from each rank
    integer, dimension(:), intent(in) :: recvcount
    !> offset of recieve things
    integer, dimension(:), intent(in) :: recvoffset
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_alltoallv(sbuf, sendcount, sendoffset, MPI_CHARACTER, rbuf, recvcount, recvoffset, MPI_CHARACTER, mw%comm, mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_alltoallv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_alltoallv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> wrappers for mpi allgatherv
subroutine scatterv_character(mw,sbuf,rbuf,ctr,offset,from,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    character, dimension(:), intent(in) :: sbuf
    !> buffer to recieve in
    character, dimension(:), intent(out) :: rbuf
    !> number of things to receive from each rank
    integer, dimension(:), intent(in) :: ctr
    !> offset
    integer, dimension(:), intent(in) :: offset
    !> which rank to scatter from
    integer, intent(in) :: from
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_scatterv(sbuf,ctr,offset,MPI_CHARACTER,rbuf,ctr(mw%r+1),MPI_CHARACTER,from,mw%comm,mw%error)

    !int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
                     !MPI_Datatype sendtype, void *recvbuf, int recvcount,
                     !MPI_Datatype recvtype, int root, MPI_Comm comm)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_allgatherv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_allgatherv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> wrappers for mpi send
subroutine send_integer(mw,x,destrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    integer, intent(in) :: x
    !> destination rank
    integer, intent(in) :: destrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_send(x,1,MPI_INTEGER,destrank,tag,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine send_1d_integer(mw,x,destrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    integer, dimension(:), intent(in) :: x
    !> destination rank
    integer, intent(in) :: destrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_send(x,size(x),MPI_INTEGER,destrank,tag,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine send_2d_integer(mw,x,destrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    integer, dimension(:,:), intent(in) :: x
    !> destination rank
    integer, intent(in) :: destrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_send(x,size(x),MPI_INTEGER,destrank,tag,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine send_1d_real(mw,x,destrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    real(r8), dimension(:), intent(in) :: x
    !> destination rank
    integer, intent(in) :: destrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_send(x,size(x),MPI_DOUBLE_PRECISION,destrank,tag,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine send_2d_real(mw,x,destrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    real(r8), dimension(:,:), intent(in) :: x
    !> destination rank
    integer, intent(in) :: destrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_send(x,size(x),MPI_DOUBLE_PRECISION,destrank,tag,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine send_character(mw,x,destrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    character, dimension(:), intent(in) :: x
    !> destination rank
    integer, intent(in) :: destrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_send(x,size(x),MPI_CHARACTER,destrank,tag,mw%comm,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> wrappers for mpi recv
subroutine recv_integer(mw,x,fromrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    integer, intent(out) :: x
    !> destination rank
    integer, intent(in) :: fromrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_recv(x,1,MPI_INTEGER,fromrank,tag,mw%comm,mw%stat,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine recv_1d_integer(mw,x,fromrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    integer, dimension(:), intent(out) :: x
    !> destination rank
    integer, intent(in) :: fromrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_recv(x,size(x),MPI_INTEGER,fromrank,tag,mw%comm,mw%stat,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine recv_2d_integer(mw,x,fromrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    integer, dimension(:,:), intent(out) :: x
    !> destination rank
    integer, intent(in) :: fromrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_recv(x,size(x),MPI_INTEGER,fromrank,tag,mw%comm,mw%stat,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine recv_1d_real(mw,x,fromrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    real(r8), dimension(:), intent(out) :: x
    !> destination rank
    integer, intent(in) :: fromrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_recv(x,size(x),MPI_DOUBLE_PRECISION,fromrank,tag,mw%comm,mw%stat,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine recv_2d_real(mw,x,fromrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    real(r8), dimension(:,:), intent(out) :: x
    !> destination rank
    integer, intent(in) :: fromrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_recv(x,size(x),MPI_DOUBLE_PRECISION,fromrank,tag,mw%comm,mw%stat,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine
subroutine recv_character(mw,x,fromrank,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    character, dimension(:), intent(out) :: x
    !> destination rank
    integer, intent(in) :: fromrank
    !> message tag
    integer, intent(in) :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    call mpi_recv(x,size(x),MPI_CHARACTER,fromrank,tag,mw%comm,mw%stat,mw%error)
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> wrappers for mpi isend
subroutine isend_character(mw,x,destrank,request,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    character, dimension(:), intent(in) :: x
    !> destination rank
    integer, intent(in) :: destrank
    !> request
    integer, intent(out) :: request
    !> message tag, if we want it
    integer, intent(in), optional :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    if ( present(tag) ) then
        call mpi_isend(x,size(x),MPI_CHARACTER,destrank,tag,mw%comm,request,mw%error)
    else
        call mpi_isend(x,size(x),MPI_CHARACTER,destrank,0,mw%comm,request,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_send exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> wrappers for mpi irecv
subroutine irecv_character(mw,x,fromrank,request,tag,filename,linenumber)
    !> MPI helper
    class(lo_mpi_helper), intent(inout) :: mw
    !> buffer to send
    character, dimension(:), intent(out) :: x
    !> destination rank
    integer, intent(in) :: fromrank
    !> request
    integer, intent(out) :: request
    !> message tag, if we want it
    integer, intent(in), optional :: tag
    !> filename we call from
    character(len=*), intent(in), optional :: filename
    !> line number we call from
    integer, intent(in), optional :: linenumber

    ! Communicate!
    if ( present(tag) ) then
        call mpi_irecv(x,size(x),MPI_CHARACTER,fromrank,tag,mw%comm,request,mw%error)
    else
        call mpi_irecv(x,size(x),MPI_CHARACTER,fromrank,MPI_ANY_TAG,mw%comm,request,mw%error)
    endif
    ! Check that things went ok
    if ( mw%error .ne. 0 ) then
        if ( present(filename) .and. present(linenumber) ) then
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,filename,linenumber,mw%comm)
        else
            call lo_stop_gracefully(['mpi_recv exit code '//tochar(mw%error)],lo_exitcode_mpi,__FILE__,__LINE__,mw%comm)
        endif
    endif
end subroutine

!> return the proper MPI operation from my string thingy.
function mpi_operation_code(operation) result(code)
    !> string describing the operation
    character(len=*), intent(in) :: operation
    !> proper MPI code
    integer :: code

    ! Decide what operation to use
    select case(trim(adjustl(operation)))
    case('sum')
        code=MPI_SUM
    case('max')
        code=MPI_MAX
    case('min')
        code=MPI_MIN
    case('and')
        code=MPI_LAND
    case('or')
        code=MPI_LOR
    case default
        call lo_stop_gracefully(['Unknown MPI operation: '//operation],lo_exitcode_mpi,__FILE__,__LINE__)
    end select
end function

!> Stop gracefully, but this version handles and kills MPI as well. Had to repeat it, only way to not get circular dependencies I think.
subroutine lo_stop_gracefully(msg,exitcode,filename,line,communicator)
    !> message
    character(len=*), dimension(:), intent(in) :: msg
    !> what exit code to give?
    integer, intent(in) :: exitcode
    !> which file?
    character(len=*), intent(in), optional :: filename
    !> which line?
    integer, intent(in), optional :: line
    !> if it's an MPI application, send in the communicator
    integer, intent(in), optional :: communicator

    integer :: rank,nrank,i

    write(*,*) ''
    write(*,*) 'ERROR'
    select case(exitcode)
    case(1)
        write(*,*) 'exit code 1: unexpected dimensions'
    case(2)
        write(*,*) 'exit code 1: blas/lapack returned nonzero exitcode'
    case(3)
        write(*,*) 'exit code 3: unphysical value detected'
    case(4)
        write(*,*) 'exit code 4: symmetry error'
    case(5)
        write(*,*) 'exit code 5: bad parameters sent to routine'
    case(6)
        write(*,*) 'exit code 6: I/O error'
    case(7)
        write(*,*) 'exit code 7: MPI error'
    end select
    write(*,*) ''
    do i=1,size(msg)
        write(*,*) trim(msg(i))
    enddo
    write(*,*) ''
    if ( present(filename) ) write(*,*) '    occurs in file: ',filename
    if ( present(line) )     write(*,*) '    occurs on line: ',tochar(line)

    if ( present(communicator) ) then
        call mpi_comm_rank(communicator,rank,lo_status)
        call mpi_comm_size(communicator,nrank,lo_status)
        write(*,*) '       active rank: ',tochar(rank),' out of ',tochar(nrank)
        call mpi_finalize(lo_status)
    endif

    ! Seems I can not use a variable as the errorcode with Ifort.
    ! No, worries, got a brute force solution:
    select case(exitcode)
    case(1)
        error stop 1
    case(2)
        error stop 2
    case(3)
        error stop 3
    case(4)
        error stop 4
    case(5)
        error stop 5
    case(6)
        error stop 6
    case(7)
        error stop 7
    case default
        error stop
    end select
end subroutine

end module
