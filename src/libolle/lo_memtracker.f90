module lo_memtracker
!!
!! Barebones memory tracker.
!!
use konstanter, only: i4,i8,r8,lo_iou,lo_hugeint
use mpi
implicit none
private
public :: lo_mem_helper

!> keep track of the number of allocations/deallocations and memory used.
type lo_mem_helper
    !> number of calls to allocate
    integer(i8) :: n_allocate=-lo_hugeint
    !> number of calls to deallocate
    integer(i8) :: n_deallocate=-lo_hugeint
    !> current amount of memory allocated
    integer(i8) :: persistent_scalable=-lo_hugeint
    integer(i8) :: persistent_nonscalable=-lo_hugeint
    integer(i8) :: temporary_scalable=-lo_hugeint
    integer(i8) :: temporary_nonscalable=-lo_hugeint
    !> keep track of peak memory of each kind?
    integer(i8) :: peak_persistent_scalable=-lo_hugeint
    integer(i8) :: peak_persistent_nonscalable=-lo_hugeint
    integer(i8) :: peak_temporary_scalable=-lo_hugeint
    integer(i8) :: peak_temporary_nonscalable=-lo_hugeint
    !> space for catching error codes
    integer, private :: ierr=0
    !> checkpoint memory usage
    integer(i8), private :: chkpoint_persistent_scalable=-lo_hugeint
    integer(i8), private :: chkpoint_persistent_nonscalable=-lo_hugeint
    integer(i8), private :: chkpoint_temporary_scalable=-lo_hugeint
    integer(i8), private :: chkpoint_temporary_nonscalable=-lo_hugeint
    contains
        !> initialize the tracker
        procedure :: init
        !> dump info to stdout
        procedure :: dump
        !> check that there is no temporary memory allocated
        procedure :: assertzero
        !> insert checkpoint
        procedure :: tick
        !> test checkpoint checkpoint (tock)
        procedure :: tock

        !> allocate memory
        generic :: allocate=>&
            allocate_1d_i4,&
            allocate_2d_i4,&
            allocate_3d_i4,&
            allocate_1d_r8,&
            allocate_2d_r8,&
            allocate_3d_r8,&
            allocate_4d_r8,&
            allocate_5d_r8,&
            allocate_1d_c8,&
            allocate_2d_c8,&
            allocate_3d_c8,&
            allocate_4d_c8,&
            allocate_5d_c8,&
            allocate_1d_char,&
            allocate_1d_logical,&
            allocate_2d_logical

        !> free memory
        generic :: deallocate=>&
            deallocate_1d_i4,&
            deallocate_2d_i4,&
            deallocate_3d_i4,&
            deallocate_1d_r8,&
            deallocate_2d_r8,&
            deallocate_3d_r8,&
            deallocate_4d_r8,&
            deallocate_5d_r8,&
            deallocate_1d_c8,&
            deallocate_2d_c8,&
            deallocate_3d_c8,&
            deallocate_4d_c8,&
            deallocate_5d_c8,&
            deallocate_1d_char,&
            deallocate_1d_logical,&
            deallocate_2d_logical

        procedure, private :: allocate_1d_i4
        procedure, private :: allocate_2d_i4
        procedure, private :: allocate_3d_i4

        procedure, private :: allocate_1d_r8
        procedure, private :: allocate_2d_r8
        procedure, private :: allocate_3d_r8
        procedure, private :: allocate_4d_r8
        procedure, private :: allocate_5d_r8

        procedure, private :: allocate_1d_c8
        procedure, private :: allocate_2d_c8
        procedure, private :: allocate_3d_c8
        procedure, private :: allocate_4d_c8
        procedure, private :: allocate_5d_c8

        procedure, private :: allocate_1d_char

        procedure, private :: allocate_1d_logical
        procedure, private :: allocate_2d_logical

        procedure, private :: deallocate_1d_i4
        procedure, private :: deallocate_2d_i4
        procedure, private :: deallocate_3d_i4

        procedure, private :: deallocate_1d_r8
        procedure, private :: deallocate_2d_r8
        procedure, private :: deallocate_3d_r8
        procedure, private :: deallocate_4d_r8
        procedure, private :: deallocate_5d_r8

        procedure, private :: deallocate_1d_c8
        procedure, private :: deallocate_2d_c8
        procedure, private :: deallocate_3d_c8
        procedure, private :: deallocate_4d_c8
        procedure, private :: deallocate_5d_c8

        procedure, private :: deallocate_1d_char

        procedure, private :: deallocate_1d_logical
        procedure, private :: deallocate_2d_logical
end type

! Should probably throw an error if I try to allocate something very large.
! This can be supressed, but is good for catching errors when you allocate
! something you think is small but it turns out to be large!
integer, parameter :: size_for_throwing_error=1000000

contains

!> create the memory tracker
subroutine init(mem)
    !> memory tracker
    class(lo_mem_helper), intent(out) :: mem

    ! Just set everything to zero
    mem%n_allocate=0
    mem%n_deallocate=0
    mem%persistent_scalable=0
    mem%persistent_nonscalable=0
    mem%temporary_scalable=0
    mem%temporary_nonscalable=0
    mem%peak_persistent_scalable=0
    mem%peak_persistent_nonscalable=0
    mem%peak_temporary_scalable=0
    mem%peak_temporary_nonscalable=0
    mem%ierr=0
end subroutine

!> Make a snapshot of the current usage.
subroutine tick(mem)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    ! Make a note what the current usage is.
    mem%chkpoint_persistent_scalable    = mem%persistent_scalable
    mem%chkpoint_persistent_nonscalable = mem%persistent_nonscalable
    mem%chkpoint_temporary_scalable     = mem%temporary_scalable
    mem%chkpoint_temporary_nonscalable  = mem%temporary_nonscalable
end subroutine

!> Assert that no memory usage has change since last 'tick'.
subroutine tock(mem,file,line,comm)
    !> memory tracker
    class(lo_mem_helper), intent(in) :: mem
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> MPI communicator, if we are doing things in parallel.
    integer, intent(in), optional :: comm

    integer(i8) :: i
    integer :: j,k,rank

    i=0_i8
    i=i+abs(mem%persistent_scalable    - mem%chkpoint_persistent_scalable    )
    i=i+abs(mem%persistent_nonscalable - mem%chkpoint_persistent_nonscalable )
    i=i+abs(mem%temporary_scalable     - mem%chkpoint_temporary_scalable     )
    i=i+abs(mem%temporary_nonscalable  - mem%chkpoint_temporary_nonscalable  )
    if ( present(comm) ) then
        ! Add up over ranks
        j=int(abs(i),i4)
        call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_SUM,comm,k)
        i=j
    endif
    if ( i .eq. 0_i8 ) then
        ! All is well
        return
    else
        ! All is not well.
        if ( present(comm) ) then
            ! Get the rank thing.
            call mpi_comm_rank(comm,rank,k)
        else
            rank=0
        endif

        if ( rank == 0 ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'ERROR: memory changed since last "tick".'
            write(lo_iou,*) '                   file:',trim(file)
            write(lo_iou,*) '                   line:',line
        endif

        if ( present(comm) ) then
            j=int(abs(mem%persistent_scalable),i4)
            call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_MAX,comm,k)
            if ( rank == 0 ) write(lo_iou,*) '    persistent_scalable:',j

            j=int(abs(mem%persistent_nonscalable),i4)
            call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_MAX,comm,k)
            if ( rank == 0 ) write(lo_iou,*) ' persistent_nonscalable:',j

            j=int(abs(mem%temporary_scalable),i4)
            call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_MAX,comm,k)
            if ( rank == 0 ) write(lo_iou,*) '     temporary_scalable:',j

            j=int(abs(mem%temporary_nonscalable),i4)
            call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_MAX,comm,k)
            if ( rank == 0 ) write(lo_iou,*) '  temporary_nonscalable:',j
        else
            write(lo_iou,*) '    persistent_scalable:',mem%persistent_scalable
            write(lo_iou,*) ' persistent_nonscalable:',mem%persistent_nonscalable
            write(lo_iou,*) '     temporary_scalable:',mem%temporary_scalable
            write(lo_iou,*) '  temporary_nonscalable:',mem%temporary_nonscalable
        endif

        ! Kill MPI, if applicable
        if ( present(comm) ) then
            call mpi_finalize(k)
        endif
        ! And kill
        stop
    endif
end subroutine

!> make sure we do not have anything allocated. If called with communicator, it also acts as barrier.
subroutine assertzero(mem,file,line,comm)
    !> memory tracker
    class(lo_mem_helper), intent(in) :: mem
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> MPI communicator, if we are doing things in parallel.
    integer, intent(in), optional :: comm

    integer(i8) :: i
    integer :: j,k,rank

    i=0_i8
    i=i+mem%persistent_scalable
    i=i+mem%persistent_nonscalable
    i=i+mem%temporary_scalable
    i=i+mem%temporary_nonscalable
    if ( present(comm) ) then
        ! Add up over ranks
        j=int(abs(i),i4)
        call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_SUM,comm,k)
        i=j
    endif
    if ( i .eq. 0_i8 ) then
        ! All is well
        return
    else
        ! All is not well.
        if ( present(comm) ) then
            ! Get the rank thing.
            call mpi_comm_rank(comm,rank,k)
        else
            rank=0
        endif

        if ( rank == 0 ) then
            write(lo_iou,*) ''
            write(lo_iou,*) 'ERROR: memory not cleared.'
            write(lo_iou,*) '                   file:',trim(file)
            write(lo_iou,*) '                   line:',line
        endif

        if ( present(comm) ) then
            j=int(abs(mem%persistent_scalable),i4)
            call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_MAX,comm,k)
            if ( rank == 0 ) write(lo_iou,*) '    persistent_scalable:',j

            j=int(abs(mem%persistent_nonscalable),i4)
            call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_MAX,comm,k)
            if ( rank == 0 ) write(lo_iou,*) ' persistent_nonscalable:',j

            j=int(abs(mem%temporary_scalable),i4)
            call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_MAX,comm,k)
            if ( rank == 0 ) write(lo_iou,*) '     temporary_scalable:',j

            j=int(abs(mem%temporary_nonscalable),i4)
            call mpi_allreduce(MPI_IN_PLACE,j,1,MPI_INTEGER,MPI_MAX,comm,k)
            if ( rank == 0 ) write(lo_iou,*) '  temporary_nonscalable:',j
        else
            write(lo_iou,*) '    persistent_scalable:',mem%persistent_scalable
            write(lo_iou,*) ' persistent_nonscalable:',mem%persistent_nonscalable
            write(lo_iou,*) '     temporary_scalable:',mem%temporary_scalable
            write(lo_iou,*) '  temporary_nonscalable:',mem%temporary_nonscalable
        endif

        ! Kill MPI, if applicable
        if ( present(comm) ) then
            call mpi_finalize(k)
        endif
        ! And kill
        stop
    endif
end subroutine

!> dump snapshot
subroutine dump(mem,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(in) :: mem
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    write(lo_iou,*) ''
    write(lo_iou,*) 'TEMPORARY MEMORY DUMP'
    write(lo_iou,*) '                     file:',trim(file)
    write(lo_iou,*) '                     line:',line
    write(lo_iou,*) '      persistent scalable:',mem%persistent_scalable
    write(lo_iou,*) '   persistent nonscalable:',mem%persistent_nonscalable
    write(lo_iou,*) '       temporary scalable:',mem%temporary_scalable
    write(lo_iou,*) '    temporary nonscalable:',mem%temporary_nonscalable
end subroutine

subroutine allocate_1d_i4(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    integer(i4), dimension(:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( n .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( n .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_2d_i4(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    integer(i4), dimension(:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(2), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_3d_i4(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    integer(i4), dimension(:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(3), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_1d_r8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    real(r8), dimension(:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( n .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( n .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_2d_r8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    real(r8), dimension(:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(2), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_3d_r8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    real(r8), dimension(:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(3), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_4d_r8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    real(r8), dimension(:,:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(4), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3),n(4)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_5d_r8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    real(r8), dimension(:,:,:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(5), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3),n(4),n(5)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_1d_c8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    complex(r8), dimension(:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( n .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_2d_c8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    complex(r8), dimension(:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(2), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_3d_c8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    complex(r8), dimension(:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(3), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_4d_c8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    complex(r8), dimension(:,:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(4), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3),n(4)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_5d_c8(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    complex(r8), dimension(:,:,:,:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(5), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2),n(3),n(4),n(5)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_1d_char(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    character, dimension(:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( n .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_1d_logical(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    logical, dimension(:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size
    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( n .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine
subroutine allocate_2d_logical(mem,x,n,persistent,scalable,file,line,supress_error)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to allocate
    logical, dimension(:,:), allocatable, intent(inout) :: x
    !> size to allocate
    integer, dimension(2), intent(in) :: n
    !> should this be a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> is this a scalable or non-scalable array?
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line
    !> suppress error for too large array
    logical, intent(in), optional :: supress_error

    logical :: check_size

    ! maybe do not check for large arrays?
    if ( present(supress_error) ) then
        check_size=.not.supress_error
    else
        check_size=.true.
    endif
    ! allocate if everything seems sensible
    if ( minval(n) .le. 0 ) then
        ! Throw error because I can not allocate empty thingy
        write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
        write(lo_iou,*) '       that could not have been intended.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    ! elseif ( product(n) .gt. size_for_throwing_error .and. check_size ) then
    !     ! Throw error because the array will be too large
    !     write(lo_iou,*) 'ERROR: trying to allocate an array of size ',n
    !     write(lo_iou,*) '       that could not have been intended.'
    !     stop
    else
        ! Things seem sensible, go ahead with the allocation
        allocate(x(n(1),n(2)),stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            ! This went badly
            write(lo_iou,*) 'ERROR: failed allocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_allocate=mem%n_allocate+1
            if ( persistent ) then
                if ( scalable ) then
                    mem%persistent_scalable=mem%persistent_scalable+size(x)*storage_size(x)
                    mem%peak_persistent_scalable=max(mem%peak_persistent_scalable,mem%persistent_scalable)
                else
                    mem%persistent_nonscalable=mem%persistent_nonscalable+size(x)*storage_size(x)
                    mem%peak_persistent_nonscalable=max(mem%peak_persistent_nonscalable,mem%persistent_nonscalable)
                endif
            else
                if ( scalable ) then
                    mem%temporary_scalable=mem%temporary_scalable+size(x)*storage_size(x)
                    mem%peak_temporary_scalable=max(mem%peak_temporary_scalable,mem%temporary_scalable)
                else
                    mem%temporary_nonscalable=mem%temporary_nonscalable+size(x)*storage_size(x)
                    mem%peak_temporary_nonscalable=max(mem%peak_temporary_nonscalable,mem%temporary_nonscalable)
                endif
            endif
        endif
    endif
end subroutine

subroutine deallocate_1d_i4(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    integer(i4), dimension(:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_2d_i4(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    integer(i4), dimension(:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_3d_i4(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    integer(i4), dimension(:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_1d_r8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    real(r8), dimension(:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_2d_r8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    real(r8), dimension(:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_3d_r8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    real(r8), dimension(:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_4d_r8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    real(r8), dimension(:,:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_5d_r8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    real(r8), dimension(:,:,:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_1d_c8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    complex(r8), dimension(:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_2d_c8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    complex(r8), dimension(:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_3d_c8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    complex(r8), dimension(:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_4d_c8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    complex(r8), dimension(:,:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_5d_c8(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    complex(r8), dimension(:,:,:,:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_1d_char(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    character, dimension(:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_1d_logical(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    logical, dimension(:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine
subroutine deallocate_2d_logical(mem,x,persistent,scalable,file,line)
    !> memory tracker
    class(lo_mem_helper), intent(inout) :: mem
    !> array to deallocate
    logical, dimension(:,:), allocatable, intent(inout) :: x
    !> was this a persistent array or a temporary array
    logical, intent(in) :: persistent
    !> was this scalable or non-scalable memory
    logical, intent(in) :: scalable
    !> filename
    character(len=*), intent(in) :: file
    !> line number
    integer, intent(in) :: line

    integer :: n

    if ( .not.allocated(x) ) then
        write(lo_iou,*) 'ERROR: trying to deallocate array that is already deallocated.'
        write(lo_iou,*) '   file: ',trim(file)
        write(lo_iou,*) '   line: ',line
        stop
    else
        n=size(x)*storage_size(x)
    endif
    if ( persistent ) then
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%persistent_scalable=mem%persistent_scalable-n
            else
                mem%persistent_nonscalable=mem%persistent_nonscalable-n
            endif
        endif
    else
        deallocate(x,stat=mem%ierr)
        if ( mem%ierr .ne. 0 ) then
            write(lo_iou,*) 'ERROR: failed deallocation.'
            write(lo_iou,*) '   file: ',trim(file)
            write(lo_iou,*) '   line: ',line
            stop
        else
            mem%n_deallocate=mem%n_deallocate+1
            if ( scalable ) then
                mem%temporary_scalable=mem%temporary_scalable-n
            else
                mem%temporary_nonscalable=mem%temporary_nonscalable-n
            endif
        endif
    endif
end subroutine

end module
