module lo_timetracker
!!
!! Barebones timer, just to make code less verbose
!!
use konstanter, only: r8,lo_iou,lo_hugeint,lo_huge,lo_exitcode_param
use mpi_wrappers, only: lo_mpi_helper,lo_stop_gracefully
use gottochblandat, only: walltime,lo_mean,lo_stddev
implicit none
private
public :: lo_timer

! Will increase this as needed
integer, parameter :: max_n_things_to_time=10
integer, parameter :: max_string_length=200

!> keep track of the number of allocations/deallocations and memory used.
type lo_timer
    !> current number of things I have timed
    integer :: n=-lo_hugeint
    !> runtime of things
    real(r8), dimension(:), allocatable :: mellantid
    !> what have I been timing
    character(len=max_string_length), dimension(:), allocatable :: whatItimed
    !> time when we started
    real(r8), private :: time_zero=-lo_huge
    !> time when I last timed something
    real(r8), private :: time_tick=-lo_huge
    !> time when we stopped
    real(r8), private :: time_end=-lo_huge
    !> is the timer running?
    logical, private :: running=.false.
    contains
        !> initialize the tracker
        procedure :: start
        !> stop timing
        procedure :: stop
        !> dump info to stdout
        procedure :: dump
        !> start timing around checkpoint
        procedure :: tick
        !> end timing around checkpoint
        procedure :: tock
end type

contains

!> dump info from timer
subroutine dump(tmr,mw,msg)
    !> timer
    class(lo_timer), intent(inout) :: tmr
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> header message
    character(len=*), intent(in) :: msg

    real(r8), dimension(:,:), allocatable :: ri,rj
    character(len=100) :: opf
    integer :: i
    ! Do something when the timer is running? Not sure.
    ! Anyway, try to dump things

    if ( tmr%n .gt. 0 ) then
        allocate(ri(mw%n,tmr%n+2))
        allocate(rj(4,tmr%n+2))
        ri=0.0_r8
        rj=0.0_r8
    else
        ! We don't have anything timed. Likely user error.
        call lo_stop_gracefully(['Nothing has been timed.'],lo_exitcode_param,__FILE__,__LINE__,mw%comm)
    endif

    ! individual timings
    do i=1,tmr%n
        ri(mw%r+1,i)=tmr%mellantid(i)
    enddo
    ! total
    ri(mw%r+1,tmr%n+2)=tmr%time_end-tmr%time_zero
    ! idle time
    ri(mw%r+1,tmr%n+1)=ri(mw%r+1,tmr%n+2)-sum(ri(mw%r+1,1:tmr%n))
    ! sync across ranks?
    call mw%allreduce('sum',ri)

    ! Massage a little
    do i=1,tmr%n+2
        rj(1,i)=lo_mean(ri(:,i))
        rj(2,i)=lo_stddev(ri(:,i))
        rj(3,i)=maxval(ri(:,i))
    enddo
    rj(4,1:tmr%n+1)=100.0_r8*rj(1,1:tmr%n+1)/sum(rj(1,1:tmr%n+1))
    rj(4,tmr%n+2)=100.0_r8

    if ( mw%talk ) then

        opf='(1X,A40,3(1X,F12.3),1X,F6.2,"%")'

        write(lo_iou,*) trim(adjustl(msg))

        opf='(41X,3(1X,A12))'
        write(lo_iou,opf) 'avg (s)','dev (s)','max (s)'
        opf='(1X,A40,3(1X,F12.3),1X,F6.2,"%")'
        do i=1,tmr%n
            write(lo_iou,opf) trim(adjustl(tmr%whatItimed(i))),rj(1:4,i)
        enddo
        write(lo_iou,opf) 'idle',rj(1:4,tmr%n+1)
        write(lo_iou,opf) 'total',rj(1:4,tmr%n+2)
    endif

    deallocate(ri)
    deallocate(rj)
end subroutine

!> start the timer
subroutine start(tmr)
    !> timer
    class(lo_timer), intent(out) :: tmr

    tmr%time_zero=walltime()
    tmr%n=0
    allocate(tmr%mellantid(max_n_things_to_time))
    allocate(tmr%whatItimed(max_n_things_to_time))
    tmr%mellantid=0.0_r8
    tmr%whatItimed='nothing'
    tmr%time_tick=tmr%time_zero
    tmr%time_end=tmr%time_zero ! not sure what to set this to.
    tmr%running=.true.
end subroutine

!> make a 'tick' of the current time
subroutine tick(tmr)
    !> timer
    class(lo_timer), intent(inout) :: tmr

    tmr%time_tick=walltime()
end subroutine

!> add a measurement since last 'tick'
subroutine tock(tmr,what)
    !> timer
    class(lo_timer), intent(inout) :: tmr
    !> what I timed
    character(len=*), intent(in) :: what

    real(r8), dimension(:), allocatable :: dr
    real(r8) :: t1
    integer :: n,i
    character(len=max_string_length), dimension(:), allocatable :: ds

    ! Make a note of the time
    t1=walltime()
    ! First check if the field already exists, and if so, add
    ! the timing information there.
    do i=1,tmr%n
        if ( trim(adjustl(tmr%whatItimed(i))) .eq. trim(adjustl(what)) ) then
            tmr%mellantid(i)=tmr%mellantid(i)+t1-tmr%time_tick
            tmr%time_tick=t1
            return
        endif
    enddo

    ! Ok, if we make it here it did not exist. Check if I have to extend arrays?
    if ( tmr%n+1 .gt. size(tmr%mellantid) ) then
        n=size(tmr%mellantid)+max_n_things_to_time
        allocate(dr(n))
        allocate(ds(n))
        dr=0.0_r8
        ds='nothing'
        dr(1:tmr%n)=tmr%mellantid
        ds(1:tmr%n)=tmr%whatItimed
        call move_alloc(dr,tmr%mellantid)
        call move_alloc(ds,tmr%whatItimed)
    endif

    ! Store this as a new entry
    tmr%n=tmr%n+1
    tmr%mellantid(tmr%n)=t1-tmr%time_tick
    tmr%whatItimed(tmr%n)=trim(adjustl(what))
    tmr%time_tick=t1
end subroutine

subroutine stop(tmr)
    !> timer
    class(lo_timer), intent(inout) :: tmr

    tmr%time_end=walltime()
    tmr%time_tick=tmr%time_end
    tmr%running=.false.
end subroutine


end module
