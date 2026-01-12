module lo_sorting
!!
!! Sort things in a sensible way.
!!
use konstanter, only: r8,lo_hugeint,lo_huge,lo_tol
use lo_memtracker, only: lo_mem_helper

implicit none
private

public :: lo_qsort
public :: lo_return_unique
public :: lo_return_unique_indices

! Quicksort things
interface lo_qsort
   module procedure quick_sort_real
   module procedure quick_sort_int
   !module procedure quick_sortv
   !module procedure quick_sortiv
   !module procedure sortchar
end interface lo_qsort

! Return the unique things from a larger set of things
interface lo_return_unique
    module procedure lo_return_unique_characters
    module procedure lo_return_unique_integers
    module procedure lo_return_unique_integer_columns
    module procedure lo_return_unique_double_columns
    module procedure lo_return_unique_double_matrices
end interface lo_return_unique

! Return unique things from many things as an index array
interface lo_return_unique_indices
    module procedure lo_return_unique_indices_real_vectors
    module procedure lo_return_unique_indices_integers
    module procedure lo_return_unique_indices_integer_vectors
end interface lo_return_unique_indices

contains

!> Return a unique list of characters
subroutine lo_return_unique_characters(a,u,ind)
    !> array of strings
    character(len=*), dimension(:), intent(in) :: a
    !> unique array of strings
    character(len=len(a(1))), dimension(:), allocatable, intent(out) :: u
    !> auxiliary index saying which unique each of the original correspond to
    integer, dimension(:), intent(out), optional :: ind

    character(len=len(a(1))), dimension(:), allocatable :: dum
    integer :: i,j,k

    allocate(dum(size(a)))
    dum=''
    k=0
    vl: do i=1,size(a)
        do j=1,k
            if ( trim(a(i)) .eq. trim(dum(j)) ) then
                if ( present(ind) ) ind(i)=j
                cycle vl
            endif
        enddo
        k=k+1
        dum(k)=trim(a(i))
        if ( present(ind) ) ind(i)=k
    enddo vl

    allocate(u(k))
    u=dum(1:k)
    deallocate(dum)
end subroutine

!> From an integer array, return only the unique values. Two algorithms available, depending on what one can expect from the data.
subroutine lo_return_unique_integers(a,u,mem,ind,sort)
    !> list of integers
    integer, dimension(:), intent(in) :: a
    !> unique integers
    integer, dimension(:), allocatable, intent(out) :: u
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> auxiliary index saying which unique each of the original correspond to
    integer, dimension(:), intent(out), optional :: ind
    !> sorting
    logical, intent(in), optional :: sort

    integer, dimension(:), allocatable :: di,dj,dk
    integer :: i,j,ctr,n
    logical :: sortvals

    ! Decide wich algorithm to use
    if ( present(sort) ) then
        sortvals=sort
    else
        sortvals=.false.
    endif
    n=size(a)

    if ( present(ind) ) then
        ind=0
    endif

    if ( sortvals ) then
        ! This version scales as quicksort. With a large number of unique, it is
        ! faster. Get a copy and a buffer array, sort the copy.
        call mem%allocate(di,n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%allocate(dj,n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        di=a
        dj=0
        if ( present(ind) ) then
            ! Count and store the unique, also keep track of which unique each of the original becomes.
            call mem%allocate(dk,n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
            dk=0
            call lo_qsort(di,dk)
            ctr=1           ! at least one unique
            dj(ctr)=di(1)   ! it's the first one
            ind(1)=ctr
            do i=2,n
                if ( di(i) .ne. di(i-1) ) then
                    ctr=ctr+1
                    dj(ctr)=di(i)
                endif
                ind(i)=ctr
            enddo
            ! now, ind(i) points to the right thing, but is in the
            ! wrong order. this should sort it out, possibly.
            di=ind
            do i=1,n
                ind(dk(i))=di(i)
            enddo
            call mem%deallocate(dk,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        else
            ! count and store unique in buffer
            call lo_qsort(di)
            ctr=1
            dj(ctr)=di(1)
            do i=2,n
                if ( di(i) .ne. di(i-1) ) then
                    ctr=ctr+1
                    dj(ctr)=di(i)
                endif
            enddo
        endif

        ! return the unique.
        call mem%allocate(u,ctr,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        u=dj(1:ctr)
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        call mem%deallocate(dj,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    else
        ! This version is faster if we have a constant small number of
        ! unique, or if the list of not particularly long, say sub 1000.

        call mem%allocate(di,n,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)

        di=-lo_hugeint
        ctr=0
        ! Start checking list
        vl: do i=1,n
            ! Check if this lines is in the list of unique lines
            do j=1,ctr
                if ( a(i)-di(j) .eq. 0 ) then
                    if ( present(ind) ) ind(i)=j
                    cycle vl
                endif
            enddo
            ! If we made it here it's a new entry
            ctr=ctr+1
            di(ctr)=a(i)
            if ( present(ind) ) ind(i)=ctr
        enddo vl
        ! And return the list of unique
        call mem%allocate(u,ctr,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
        u=di(1:ctr)
        call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    endif
end subroutine

!> From an integer array, return only the unique columns @todo add sort to make it fast
subroutine lo_return_unique_integer_columns(a,u,mem,ind)
    !> list of integer vectors
    integer, dimension(:,:), intent(in) :: a
    !> unique integer vectors
    integer, dimension(:,:), allocatable, intent(out) :: u
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> auxiliary index saying which unique each of the original correspond to
    integer, dimension(:), intent(out), optional :: ind

    integer, dimension(:,:), allocatable :: dum
    integer :: i,j,k

    if ( present(ind) ) then
        ind=0
    endif

    call mem%allocate(dum,[size(a,1),size(a,2)],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dum=0
    k=0
    ! Start checking list
    vl: do i=1,size(a,2)
        do j=1,k
            if ( sum(abs(a(:,i)-dum(:,j))) .eq. 0 ) then
                if ( present(ind) ) ind(i)=j
                cycle vl
            endif
        enddo
        k=k+1
        dum(:,k)=a(:,i)
        if ( present(ind) ) ind(i)=k
    enddo vl

    call mem%allocate(u,[size(a,1),k],persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
    u=dum(:,1:k)
    call mem%deallocate(dum,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> Return unique columns. @todo Should be updated to NlogN via sorting, but can't get it completely stable.
subroutine lo_return_unique_double_columns(a,u,mem,tol,ind)
    !> list of vectors
    real(r8), dimension(:,:), intent(in) :: a
    !> unique vectors
    real(r8), dimension(:,:), allocatable, intent(out) :: u
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> tolerance
    real(r8), intent(in), optional :: tol
    !> auxiliary index saying which unique each of the original correspond to
    integer, dimension(:), intent(out), optional :: ind

    real(r8), dimension(:,:), allocatable :: dum
    real(r8) :: tolerance
    integer :: i,j,k

    if ( present(tol) ) then
        tolerance=tol
    else
        tolerance=lo_tol
    endif
    if ( present(ind) ) then
        ind=0
    endif

    call mem%allocate(dum,[size(a,1),size(a,2)],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dum=lo_huge
    k=0
    vl: do i=1,size(a,2)
        do j=1,k
            if ( sum(abs(a(:,i)-dum(:,j))) .lt. tolerance ) then
                if ( present(ind) ) ind(i)=j
                cycle vl
            endif
        enddo
        k=k+1
        dum(:,k)=a(:,i)
        if ( present(ind) ) ind(i)=k
    enddo vl

    call mem%allocate(u,[size(a,1),k],persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
    u=dum(:,1:k)
    call mem%deallocate(dum,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> Return unique matrices, within a tolerance
subroutine lo_return_unique_double_matrices(a,u,mem,tol,ind)
    !> list of matrices
    real(r8), dimension(:,:,:), intent(in) :: a
    !> unique matrices
    real(r8), dimension(:,:,:), allocatable, intent(out) :: u
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> tolerance
    real(r8), intent(in), optional :: tol
    !> auxiliary index saying which unique each of the original correspond to
    integer, dimension(:), intent(out), optional :: ind

    real(r8), dimension(:,:,:), allocatable :: dum
    real(r8) :: tolerance
    integer :: i,j,k

    if ( present(tol) ) then
        tolerance=tol
    else
        tolerance=lo_tol
    endif

    if ( present(ind) ) then
        ind=0
    endif

    ! Initialize counter
    call mem%allocate(dum,[size(a,1),size(a,2),size(a,3)],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dum=lo_huge
    k=0
    ! Start checking list
    vl: do i=1,size(a,3)
        do j=1,k
            if ( sum(abs(a(:,:,i)-dum(:,:,j))) .lt. tolerance ) then
                if ( present(ind) ) ind(i)=j
                cycle vl
            endif
        enddo
        k=k+1
        dum(:,:,k)=a(:,:,i)
        if ( present(ind) ) ind(i)=k
    enddo vl

    call mem%allocate(u,[size(a,1),size(a,2),k],persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
    u=dum(:,:,1:k)
    call mem%deallocate(dum,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> return unique values as an index array.
subroutine lo_return_unique_indices_real_vectors(a,ind,mem,redind,tol)
    !> list of vectors
    real(r8), dimension(:,:), intent(in) :: a
    !> indices to unique vectors
    integer, dimension(:), allocatable, intent(out) :: ind
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> optional index array that says which unique each element in the original array correspond to
    integer, dimension(:), intent(out), optional :: redind
    !> tolerance
    real(r8), intent(in), optional :: tol

    real(r8), dimension(:,:), allocatable :: dum
    real(r8) :: tolerance
    integer, dimension(:), allocatable :: di
    integer :: i,j,k

    if ( present(tol) ) then
        tolerance=tol
    else
        tolerance=lo_tol
    endif
    if ( present(redind) ) then
        redind=-lo_hugeint
    endif

    call mem%allocate(dum,[size(a,1),size(a,2)],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(di,size(a,2),persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dum=lo_huge
    di=-1
    k=0
    ! Start checking list
    vl: do i=1,size(a,2)
        do j=1,k
            if ( sum(abs(a(:,i)-dum(:,j))) .lt. tolerance ) then
                if ( present(redind) ) redind(i)=j
                cycle vl
            endif
        enddo
        k=k+1
        dum(:,k)=a(:,i)
        di(k)=i
        if ( present(redind) ) redind(i)=k
    enddo vl

    ! And return the unique
    call mem%deallocate(dum,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(ind,k,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
    ind=di(1:k)
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> return unique values as an index array.
subroutine lo_return_unique_indices_integers(a,ind,mem,redind)
    !> list of integers
    integer, dimension(:), intent(in) :: a
    !> indices to unique integers
    integer, dimension(:), allocatable, intent(out) :: ind
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> optional index array that says which unique each element in the original array correspond to
    integer, dimension(:), intent(out), optional :: redind

    integer, dimension(:), allocatable :: dum,di
    integer :: i,j,k

    if ( present(redind) ) then
        redind=-lo_hugeint
    endif

    call mem%allocate(dum,size(a),persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(di,size(a),persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dum=lo_hugeint
    di=-1
    k=0
    ! Start checking list
    vl: do i=1,size(a)
        do j=1,k
            if ( a(i) == dum(j) ) then
                if ( present(redind) ) redind(i)=j
                cycle vl
            endif
        enddo
        k=k+1
        dum(k)=a(i)
        di(k)=i
        if ( present(redind) ) redind(i)=k
    enddo vl

    ! And return the unique
    call mem%deallocate(dum,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(ind,k,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
    ind=di(1:k)
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> return unique values as an index array.
subroutine lo_return_unique_indices_integer_vectors(a,ind,mem,redind)
    !> list of integers
    integer, dimension(:,:), intent(in) :: a
    !> indices to unique integers
    integer, dimension(:), allocatable, intent(out) :: ind
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> optional index array that says which unique each element in the original array correspond to
    integer, dimension(:), intent(out), optional :: redind

    integer, dimension(:,:), allocatable :: dum
    integer, dimension(:), allocatable :: di
    integer :: i,j,k

    if ( present(redind) ) then
        redind=-1
    endif

    call mem%allocate(dum,[size(a,1),size(a,2)],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(di,size(a,2),persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    dum=lo_hugeint
    di=-1
    k=0
    ! Start checking list
    vl: do i=1,size(a,2)
        do j=1,k
            if ( sum(abs(a(:,i)-dum(:,j))) .eq. 0 ) then
                if ( present(redind) ) redind(i)=j
                cycle vl
            endif
        enddo
        k=k+1
        dum(:,k)=a(:,i)
        di(k)=i
        if ( present(redind) ) redind(i)=k
    enddo vl

    ! And return the unique
    call mem%deallocate(dum,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    call mem%allocate(ind,k,persistent=.true.,scalable=.false.,file=__FILE__,line=__LINE__)
    ind=di(1:k)
    call mem%deallocate(di,persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
end subroutine

!> Quick sort routine from:  Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to  Fortran 90", McGraw-Hill ISBN 0-07-000248-7, pages 149-150.
pure subroutine quick_sort_int(list,order)
    integer, dimension(:), intent(inout) :: list
    integer, dimension(:), intent(inout), optional :: order

    integer, parameter  :: max_simple_sort_size=6
    integer :: k
    ! initialize the order thingy
    if ( present(order) ) then
        do k=1,size(list)
            order(k) = k
        enddo
    endif

    if ( present(order) ) then
        call qsint(list,1,size(list),order)
    else
        call qsint(list,1,size(list))
    endif

    contains
    pure recursive subroutine qsint(buf,left_end,right_end,ind)
        integer, dimension(:), intent(inout) :: buf
        integer, dimension(:), intent(inout), optional :: ind
        integer, intent(in) :: left_end, right_end

        integer :: i,j,reference,tmp

        if (right_end < left_end + max_simple_sort_size) then
            ! Use interchange sort for small lists
            if ( present(ind) ) then
                call interchange_sort_int(buf,left_end, right_end,ind)
            else
                call interchange_sort_int(buf,left_end, right_end)
            endif
        else
            ! Use partition quicksort
            reference = buf((left_end+right_end)/2)
            i=left_end-1
            j=right_end+1
            do
                ! Scan list from left end until element >= reference is found
                do
                    i = i + 1
                    if ( buf(i) >= reference) exit
                enddo
                ! Scan list from right end until element <= reference is found
                do
                    j = j - 1
                    if ( buf(j) <= reference) exit
                enddo
                if (i < j) then
                    ! Swap two out-of-order elements
                    tmp=buf(i)
                    buf(i)=buf(j)
                    buf(j)=tmp
                    if ( present(order) ) then
                        tmp=ind(i)
                        ind(i)=ind(j)
                        ind(j)=tmp
                    endif
                elseif (i == j) then
                    i = i + 1
                    exit
                else
                    exit
                endif
            enddo
            if ( present(ind) ) then
                if (left_end < j) call qsint(buf,left_end,j,ind)
                if (i < right_end) call qsint(buf,i,right_end,ind)
            else
                if (left_end < j) call qsint(buf,left_end,j)
                if (i < right_end) call qsint(buf,i,right_end)
            endif
        endif
    end subroutine
end subroutine

!> Interchange sort for small integer arrays
pure subroutine interchange_sort_int(list,left_end,right_end,order)
    integer, dimension(:), intent(inout) :: list
    integer, dimension(:), intent(inout), optional :: order
    integer, intent(in) :: left_end, right_end

    integer :: i,j,tmp

    if ( present(order) ) then
        do i=left_end,right_end-1
        do j=i+1,right_end
            if ( list(i) > list(j) ) then
                tmp=list(i)
                list(i)=list(j)
                list(j)=tmp
                tmp=order(i)
                order(i)=order(j)
                order(j)=tmp
            endif
        enddo
        enddo
    else
        do i=left_end,right_end-1
        do j=i+1,right_end
            if ( list(i) > list(j) ) then
                tmp=list(i)
                list(i)=list(j)
                list(j)=tmp
            endif
        enddo
        enddo
    endif
end subroutine

!> qsort for double precision
pure subroutine quick_sort_real(list,order)
    !> array to sort
    real(r8), dimension(:), intent(inout)  :: list
    !> optionally, an index array
    integer, dimension(:), intent(inout), optional :: order

    integer :: k

    if ( present(order) ) then
        do k=1,size(list)
            order(k) = k
        enddo
    endif

    if ( present(order) ) then
        call qsreal(1,size(list),list,order)
    else
        call qsreal(1,size(list),list)
    endif
    contains

    pure recursive subroutine qsreal(left_end,right_end,buf,ind)
        integer, intent(in) :: left_end,right_end
        real(r8), dimension(:), intent(inout) :: buf
        integer, dimension(:), intent(inout), optional :: ind

        integer, parameter :: max_simple_sort_size=6
        real(r8) :: reference, temp
        integer :: i,j,itemp

        if (right_end < left_end + max_simple_sort_size) then
            ! Use interchange sort for small lists
            if ( present(order) ) then
                call interchange_sort_real(left_end,right_end,buf,ind)
            else
                call interchange_sort_real(left_end,right_end,buf)
            endif
        else
            ! Use partition ("quick") sort
            reference = buf((left_end + right_end)/2)
            i = left_end - 1; j = right_end + 1
            do
                ! Scan list from left end until element >= reference is
                ! found
                do
                    i = i + 1
                    if (buf(i) >= reference) exit
                enddo
                ! Scan list from right end until element <= reference is found
                do
                    j = j - 1
                    if (buf(j) <= reference) exit
                enddo

                if (i < j) then
                    ! Swap two out-of-order elements
                    temp = buf(i); buf(i) = buf(j); buf(j) = temp
                    if ( present(order) ) then
                        itemp = ind(i); ind(i) = ind(j); ind(j) = itemp
                    endif
                elseif (i == j) then
                    i = i + 1
                    exit
                else
                    exit
                endif
            enddo

            if ( present(ind) ) then
                if (left_end < j) call qsreal(left_end, j,buf,ind)
                if (i < right_end) call qsreal(i, right_end,buf,ind)
            else
                if (left_end < j) call qsreal(left_end, j,buf)
                if (i < right_end) call qsreal(i, right_end,buf)
            endif
        endif
    end subroutine
end subroutine

!> Interchange sort to be used for small real arrays
pure subroutine interchange_sort_real(left_end,right_end,list,order)
    integer, intent(in) :: left_end, right_end
    real(r8), dimension(:), intent(inout) :: list
    integer, dimension(:), intent(inout), optional :: order

    real(r8) :: temp
    integer :: i,j,itemp

    if ( present(order) ) then
        do i = left_end, right_end-1
        do j = i+1, right_end
            if (list(i) > list(j)) then
                temp = list(i); list(i) = list(j); list(j) = temp
                itemp = order(i); order(i) = order(j); order(j) = itemp
            endif
        enddo
        enddo
    else
        do i = left_end, right_end-1
        do j = i+1, right_end
            if (list(i) > list(j)) then
                temp = list(i)
                list(i) = list(j)
                list(j) = temp
            endif
        enddo
        enddo
    endif
end subroutine

end module
