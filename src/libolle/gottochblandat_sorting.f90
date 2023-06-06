#include "precompilerdefinitions"
submodule (gottochblandat) gottochblandat_sorting
implicit none
contains

!> Quick sort routine from:  Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to  Fortran 90", McGraw-Hill ISBN 0-07-000248-7, pages 149-150.
module pure subroutine quick_sort_int(list,order)
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

!> quicksort, but with vectors instead.
module pure subroutine quick_sortv(A,ind,tol)
    !> list of vectors
    real(flyt), dimension(:,:), intent(inout) :: A
    !> index that sorts A
    integer, dimension(:), intent(out), optional :: ind
    !> tolerance
    real(flyt), intent(in) :: tol

    real(flyt) :: diff
    real(flyt), dimension(:), allocatable :: dr
    integer, dimension(:), allocatable :: di1
    integer :: i,n,m,col,ctr,istart

    ! size of the input
    n=size(A,2)
    m=size(A,1)
    if ( n .eq. 1 ) then
        ! might get weird, and also, we have to do nothing
        if ( present(ind) ) ind=1
        return
    endif
    ! some temporary things
    allocate(dr(n))
    allocate(di1(n))

    dr=A(1,:) ! all the x-components
    call quick_sort_real(dr,di1)
    ! rearrange A so that the first element is in descending order
    A=A(:,di1)
    if ( present(ind) ) then
        do i=1,n
            ind(i)=i
        enddo
        ind=ind(di1)
    endif

    ! Now rearrange stuff per each column, in smaller and smaller groups
    dr=0.0_flyt
    di1=0
    cloop: do col=2,m
        istart=0
        ctr=1
        dr=0.0_flyt
        dr(ctr)=A(col,1)
        do i=2,n+1
            if ( i .le. n ) then
                diff=sum(abs(A(1:col-1,i)-A(1:col-1,i-1)))
            else
                diff=lo_huge
            endif
            if ( diff .lt. tol .and. i .le. n ) then
                ! one more of the same element
                ctr=ctr+1
                dr(ctr)=A(col,i)
            else
                ! now we are finished with this group, sort it
                ! we only have to do that if there is more than one row in that group.
                if ( ctr .gt. 1 ) then
                    call quick_sort_real(dr(1:ctr),di1(1:ctr))
                    ! Rearrange the array
                    A(:,istart+1:istart+ctr)=A(:,di1(1:ctr)+istart)
                    if ( present(ind) ) ind(istart+1:istart+ctr)=ind(istart+di1(1:ctr))
                endif
                ! reset counters and stuff
                if ( i .le. n ) then
                    ctr=1
                    istart=i-1
                    di1=0
                    dr=0.0_flyt
                    dr(1)=A(col,i)
                endif
            endif
        enddo
    enddo cloop
    deallocate(dr)
    deallocate(di1)
end subroutine

!> quicksort, but with integer vector
module pure subroutine quick_sortiv(A,ind)
    !> list of integer vectors
    integer, dimension(:,:), intent(inout) :: A
    !> index that sorts A
    integer, dimension(:), intent(out), optional :: ind

    integer :: i,j,nval,ncol
    integer, dimension(:,:), allocatable :: dj
    integer, dimension(:), allocatable :: di
    real(flyt), dimension(:), allocatable :: hash
    real(flyt), dimension(size(A,1)) :: mult
    real(flyt) :: imax

    nval=size(A,2)
    ncol=size(A,1)

    ! Get some prefactors
    imax=maxval(abs(A))
    do i=1,ncol
        j=ncol-i
        mult(i)=(imax)*10.0_flyt**j
    enddo

    ! Hash the data
    allocate(hash(nval))
    hash=0.0_flyt
    do i=1,nval
        hash(i)=sum(A(:,i)*mult)
    enddo

    ! Sort it
    if ( present(ind) ) then
        call quick_sort_real(hash,ind)
        allocate(dj(ncol,nval))
        dj=A
        A=dj(:,ind)
        deallocate(dj)
    else
        allocate(di(nval))
        call quick_sort_real(hash,di)
        allocate(dj(ncol,nval))
        dj=A
        A=dj(:,di)
        deallocate(dj)
        deallocate(di)
    endif
end subroutine

!> qsort for double precision
pure module subroutine quick_sort_real(list,order)
    !> array to sort
    real(flyt), dimension(:), intent(inout)  :: list
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
        real(flyt), dimension(:), intent(inout) :: buf
        integer, dimension(:), intent(inout), optional :: ind

        integer, parameter :: max_simple_sort_size=6
        real(flyt) :: reference, temp
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
    real(flyt), dimension(:), intent(inout) :: list
    integer, dimension(:), intent(inout), optional :: order

    real(flyt) :: temp
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

!> Interchange sort for small integer arrays
pure subroutine interchange_sort_int(list,left_end,right_end,order)
    integer, dimension(:), intent(inout) :: list
    integer, dimension(:), intent(inout), optional :: order
    integer, intent(in) :: left_end, right_end

    integer :: i,j,tmp

    do i=left_end,right_end-1
    do j=i+1,right_end
        if ( list(i) > list(j) ) then
            tmp=list(i)
            list(i)=list(j)
            list(j)=tmp
            if ( present(order) ) then
                tmp=order(i)
                order(i)=order(j)
                order(j)=tmp
            endif
        endif
    enddo
    enddo
end subroutine

!> qsort for character arrays
module subroutine sortchar(StringArray, indexarray, CaseInsensitive)
    ! Found this on the internet somewhere.
    character (len = *), dimension (:), intent (inout) :: StringArray
    integer, dimension(:), intent(out), optional :: indexarray
    logical, intent (in), optional :: CaseInsensitive
    !
    integer, dimension(:), allocatable :: indarr
    integer :: low, high, k !,ios
    logical :: CaseSensitive,ReturnIndex

    if (present(CaseInsensitive)) then
      CaseSensitive = .not. CaseInsensitive
    else
      CaseSensitive = .true.
    end if
    !
    if ( present(indexarray) ) then
        ReturnIndex=.true.
    else
        ReturnIndex=.false.
    endif

    low = 1
    high = size(StringArray,1)

    lo_allocate(indarr(high))

    ! fancy implied loop. Should maybe start using these
    indarr = (/ (k, k = low, high) /)

    ! Get the index array
    call quicksort_string(StringArray, indarr, low, high, CaseSensitive)
    StringArray = StringArray(indarr)
    !
    if ( ReturnIndex ) then
        indexarray=indarr
    endif
    ! The local routines for this kind of sort
    contains

    recursive subroutine quicksort_string(StringArray, indarr, low, high, CaseSensitive)
        character (len = *), dimension (:), intent (inout) :: StringArray
        integer, dimension(:), intent(inout) :: indarr
        integer, intent (in) :: low, high
        logical, intent(in) :: CaseSensitive
        !
        integer :: pivotlocation
        if (low < high) then
            call partition(StringArray, indarr, low, high, pivotlocation, CaseSensitive)
            call quicksort_string(StringArray, indarr, low, pivotlocation - 1, CaseSensitive)
            call quicksort_string(StringArray, indarr, pivotlocation + 1, high, CaseSensitive)
        end if
    end subroutine

    subroutine partition(StringArray, indarr, low, high, pivotlocation, CaseSensitive)
        character(len=*), dimension(:), intent(inout) :: StringArray
        integer, dimension(:), intent(inout) :: indarr
        integer, intent(in) :: low, high
        integer, intent(out) :: pivotlocation
        logical, intent(in) :: CaseSensitive
        !
        integer :: k, lastsmall
        call swap(indarr(low), indarr((low + high)/2))
        lastsmall = low
        do k=low+1, high
            if ( stringComp(StringArray(indarr(k)),StringArray(indarr(low)),CaseSensitive) ) then
                lastsmall = lastsmall + 1
                call swap(indarr(lastsmall), indarr(k))
            end if
        end do
        call swap(indarr(low), indarr(lastsmall))
        pivotlocation = lastsmall
    end subroutine

    subroutine swap(m,n)
        integer, intent(inout) :: m, n
        integer :: temp
        temp = m
        m = n
        n = temp
    end subroutine

    function stringComp(p,q,CaseSensitive) result(lexicalLess)
        character (len = *), intent (in) :: p, q
        logical, intent(in) :: CaseSensitive
        logical :: lexicalLess
        integer :: kq, k
        !
        if (CaseSensitive) then
            lexicalLess = p < q
        else
            kq = 1
            do k = 1, max(len_trim(p), len_trim(q))
                if (UpperCase(p(k:k)) == UpperCase(q(k:k)) ) then
                    cycle
                else
                    kq = k
                    exit
                end if
            end do
            lexicalLess = UpperCase(p(kq:kq)) < UpperCase(q(kq:kq))
        end if
    end function

    function UpperCase(letter) result(L)
        character (len=*), intent (in) :: letter
        character (len=1) :: L
        character (len=26), parameter :: Lower = "abcdefghijklmnopqrstuvwxyz"
        character (len=26), parameter :: Upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        !
        integer :: k
        ! Where in the lowercase alphabet is the character
        k = index(Lower,letter)
        if (k > 0) then
            ! Change to uppercase
            L = Upper(k:k)
        else
            ! It is already uppercase
            L = letter
        end if
    end function UpperCase
end subroutine

!> Return a unique list of characters
module pure subroutine lo_return_unique_characters(a,u)
    character(len=*), dimension(:), intent(in) :: a
    character(len=len(a(1))), dimension(:), allocatable, intent(out) :: u
    !
    character(len=len(a(1))), dimension(:), allocatable :: dum
    integer :: i,j,k
    logical :: newval

    ! Initialize counter and dummy arrays
    allocate(dum(size(a,1)))
    dum=''
    k=0
    ! Start checking list
    vl: do i=1,size(a,1)
        newval=.true.
        ! Check if this lines is in the list of unique lines
        do j=1,k
            if ( trim(a(i)) .eq. trim(dum(j)) ) then
                newval=.false.
                cycle vl
            endif
        enddo
        ! if l is 0, then this means it's a new entry, add it
        if ( newval ) then
            k=k+1
            dum(k)=trim(a(i))
        endif
    enddo vl
    ! return the unique
    allocate(u(k))
    u=dum(1:k)
    deallocate(dum)
end subroutine

!> From an integer array, return only the unique values
module pure subroutine lo_return_unique_integers(a,u,sort)
    !> list of integers
    integer, dimension(:), intent(in) :: a
    !> unique integers
    integer, dimension(:), allocatable, intent(out) :: u
    !> sorting
    logical, intent(in), optional :: sort

    integer, dimension(:), allocatable :: di,dj
    integer :: i,j,n,ctr
    logical :: sortvals

    ! Decide wich algorithm to use
    if ( present(sort) ) then
        sortvals=sort
    else
        sortvals=.false.
    endif
    n=size(a)

    if ( sortvals ) then
        ! This version scales as quicksort. With a large number of unique, it is
        ! faster. Get a copy and a buffer array, sort the copy.
        lo_allocate(di(n))
        lo_allocate(dj(n))
        di=a
        dj=0
        call qsort(di)

        ! count and store unique in buffer
        ctr=1
        dj(ctr)=di(1)
        do i=2,n
            if ( di(i) .ne. di(i-1) ) then
                ctr=ctr+1
                dj(ctr)=di(i)
            endif
        enddo

        ! return the unique
        allocate(u(ctr))
        u=dj(1:ctr)

        deallocate(di)
        deallocate(dj)
    else
        ! This version is faster if we have a constant small number of unique, then it's strictly O(N).
        allocate(di(n))
        di=-lo_hugeint
        ctr=0
        ! Start checking list
        vl: do i=1,n
            ! Check if this lines is in the list of unique lines
            do j=1,ctr
                if ( a(i)-di(j) .eq. 0 ) cycle vl
            enddo
            ! If we made it here it's a new entry
            ctr=ctr+1
            di(ctr)=a(i)
        enddo vl
        ! And return the list of unique
        allocate(u(ctr))
        u=di(1:ctr)
        deallocate(di)
    endif
end subroutine

!> From an integer array, return only the unique columns @todo add sort to make it fast
module pure subroutine lo_return_unique_integer_columns(a,u)
    !> list of integer vectors
    integer, dimension(:,:), intent(in) :: a
    !> unique integer vectors
    integer, dimension(:,:), allocatable, intent(out) :: u

    integer, dimension(:,:), allocatable :: dum
    integer :: i,j,k
    logical :: newval
    !
    ! Initialize counters
    allocate(dum(size(a,1),size(a,2)))
    dum=0
    k=0
    ! Start checking list
    vl: do i=1,size(a,2)
        newval=.true.
        ! Check if this line is in the list of unique lines
        do j=1,k
            if ( sum(abs(a(:,i)-dum(:,j))) .eq. 0 ) then
                newval=.false.
                cycle vl
            endif
        enddo
        ! this means it's a new entry, add it
        if ( newval ) then
            k=k+1
            dum(:,k)=a(:,i)
        endif
    enddo vl

    allocate(u(size(a,1),k))
    u=dum(:,1:k)
    deallocate(dum)
end subroutine

!> Return unique columns. @todo Should be update to NlogN by sorting, but can't get it completely stable.
module pure subroutine lo_return_unique_double_columns(a,u,tol,ind)
    !> list of vectors
    real(flyt), dimension(:,:), intent(in) :: a
    !> unique vectors
    real(flyt), dimension(:,:), allocatable, intent(out) :: u
    !> tolerance
    real(flyt), intent(in), optional :: tol
    !> auxiliary index saying which unique each of the original correspond to
    integer, dimension(:), intent(out), optional :: ind

    real(flyt) :: tolerance
    real(flyt), dimension(:,:), allocatable :: dum
    integer :: i,j,k !,l
    logical :: newval

    allocate(dum(size(a,1),size(a,2)))
    if ( present(tol) ) then
        tolerance=tol
    else
        tolerance=lo_tol
    endif

    if ( present(ind) ) then
        ind=0
    endif

    ! Initialize counter
    dum=lo_huge
    k=0
    ! Start checking list
    vl: do i=1,size(a,2)
        newval=.true.
        ! Check if this line is in the list of unique lines
        do j=1,k
            if ( sum(abs(a(:,i)-dum(:,j))) .lt. tolerance ) then
                newval=.false.
                if ( present(ind) ) ind(i)=j
                cycle vl
            endif
        enddo
        ! this means it's a new entry, add it
        if ( newval ) then
            k=k+1
            dum(:,k)=a(:,i)
            if ( present(ind) ) ind(i)=k
        endif
    enddo vl

    allocate(u(size(a,1),k))
    u=dum(:,1:k)
    deallocate(dum)
end subroutine

!> Return unique matrices, within a tolerance
module pure subroutine lo_return_unique_double_matrices(a,u,tol)
    !> list of matrices
    real(flyt), dimension(:,:,:), intent(in) :: a
    !> unique matrices
    real(flyt), dimension(:,:,:), allocatable, intent(out) :: u
    !> tolerance
    real(flyt), intent(in), optional :: tol

    real(flyt) :: tolerance
    real(flyt), dimension(:,:,:), allocatable :: dum
    integer :: i,j,k
    logical :: newval

    allocate(dum(size(a,1),size(a,2),size(a,3)))
    if ( present(tol) ) then
        tolerance=tol
    else
        tolerance=lo_tol
    endif
    ! Initialize counter
    dum=lo_huge
    k=0
    ! Start checking list
    vl: do i=1,size(a,3)
        newval=.true.
        ! Check if this line is in the list of unique lines
        do j=1,k
            if ( sum(abs(a(:,:,i)-dum(:,:,j))) .lt. tolerance ) then
                newval=.false.
                cycle vl
            endif
        enddo
        ! this means it's a new entry, add it
        if ( newval ) then
            k=k+1
            dum(:,:,k)=a(:,:,i)
        endif
    enddo vl

    allocate(u(size(a,1),size(a,2),k))
    u=dum(:,:,1:k)
    deallocate(dum)
end subroutine

!> Return unique scalars, within a tolerance @TODO add sort to make it NlogN
module pure subroutine lo_return_unique_doubles(a,u,tol)
    !> list of numbers
    real(flyt), dimension(:), intent(in) :: a
    !> list of unique
    real(flyt), dimension(:), allocatable, intent(out) :: u
    !> tolerance
    real(flyt), intent(in), optional :: tol
    !
    real(flyt), dimension(:), allocatable :: dum
    real(flyt) :: tolerance
    integer :: i,j,k
    logical :: newval

    if ( present(tol) ) then
        tolerance=tol
    else
        tolerance=lo_tol
    endif

    ! Initialize counter
    allocate(dum(size(a,1)))
    dum=lo_huge
    k=0
    ! Start checking list
    vl: do i=1,size(a,1)
        newval=.true.
        ! Check if this lines is in the list of unique lines
        do j=1,k
            if ( abs(a(i)-dum(j)) .lt. tolerance ) then
                newval=.false.
                cycle vl
            endif
        enddo
        ! this means it's a new entry, add it
        if ( newval ) then
            k=k+1
            dum(k)=a(i)
        endif
    enddo vl

    allocate(u(k))
    u=dum(1:k)
    deallocate(dum)
end subroutine

!> Return all permuations of n integers.
module subroutine lo_permutations(p,n)
    !> all permutations
    integer, dimension(:,:), allocatable, intent(out) :: p
    !> will permute the numbers from 1 to n
    integer, intent(in) :: n
    !
    integer, dimension(n) :: permutation
    integer :: val_min,val_max,pos_min,pos_max,i,j

    ! number of permutations:
    j=1
    do i=1,n
        j=j*i
    enddo
    lo_allocate(p(n,j))
    val_min=1
    val_max=n
    pos_min=1
    pos_max=n
    i=0
    call generate(pos_min,p,i)

    contains
    recursive subroutine generate(pos,pm,k)
        integer, intent(in) :: pos
        integer, dimension(:,:), intent(inout) :: pm
        integer, intent(inout) :: k
        integer :: val

        if (pos > pos_max) then
            ! store this one
            k=k+1
            pm(:,i)=permutation
        else
            do val = val_min,val_max
                if (.not. any (permutation(:pos-1) == val)) then
                    permutation(pos) = val
                    call generate(pos+1,pm,k)
                end if
            end do
        end if
    end subroutine
end subroutine

!> return unique values as an index array.
module subroutine lo_return_unique_indices_real_vectors(a,ind,redind,tol)
    !> list of vectors
    real(flyt), dimension(:,:), intent(in) :: a
    !> indices to unique vectors
    integer, dimension(:), allocatable, intent(out) :: ind
    !> optional index array that says which unique each element in the original array correspond to
    integer, dimension(:), intent(out), optional :: redind
    !> tolerance
    real(flyt), intent(in), optional :: tol
    !
    real(flyt) :: tolerance
    real(flyt), dimension(:,:), allocatable :: dum
    integer, dimension(:), allocatable :: di
    integer :: i,j,k !,l
    logical :: newval

    allocate(dum(size(a,1),size(a,2)))
    allocate(di(size(a,2)))
    if ( present(tol) ) then
        tolerance=tol
    else
        tolerance=lo_tol
    endif
    if ( present(redind) ) then
        redind=-lo_hugeint
    endif

    ! Initialize counter
    dum=lo_huge
    di=-1
    k=0
    ! Start checking list
    vl: do i=1,size(a,2)
        newval=.true.
        ! Check if this line is in the list of unique
        do j=1,k
            if ( sum(abs(a(:,i)-dum(:,j))) .lt. tolerance ) then
                newval=.false.
                if ( present(redind) ) redind(i)=j
                cycle vl
            endif
        enddo
        ! this means it's a new entry, add it
        if ( newval ) then
            k=k+1
            dum(:,k)=a(:,i)
            di(k)=i
            if ( present(redind) ) redind(i)=k
        endif
    enddo vl

    ! And return the unique
    allocate(ind(k))
    ind=di(1:k)
    deallocate(dum,di)
end subroutine

!> return unique values as an index array.
module subroutine lo_return_unique_indices_integers(a,ind,redind)
    !> list of integers
    integer, dimension(:), intent(in) :: a
    !> indices to unique integers
    integer, dimension(:), allocatable, intent(out) :: ind
    !> optional index array that says which unique each element in the original array correspond to
    integer, dimension(:), intent(out), optional :: redind

    integer, dimension(:), allocatable :: dum,di
    integer :: i,j,k
    logical :: newval

    allocate(dum(size(a)))
    allocate(di(size(a)))
    if ( present(redind) ) then
        redind=-lo_hugeint
    endif

    ! Initialize counter
    dum=lo_hugeint
    di=-1
    k=0
    ! Start checking list
    vl: do i=1,size(a)
        newval=.true.
        ! Check if this line is in the list of unique
        do j=1,k
            if ( abs(a(i)-dum(j)) .eq. 0 ) then
                newval=.false.
                if ( present(redind) ) redind(i)=j
                cycle vl
            endif
        enddo
        ! this means it's a new entry, add it
        if ( newval ) then
            k=k+1
            dum(k)=a(i)
            di(k)=i
            if ( present(redind) ) redind(i)=k
        endif
    enddo vl

    ! And return the unique
    allocate(ind(k))
    ind(1:k)=di(1:k)
    deallocate(dum,di)
end subroutine

!> return unique values as an index array.
module subroutine lo_return_unique_indices_integer_vectors(a,ind,redind)
    !> list of integers
    integer, dimension(:,:), intent(in) :: a
    !> indices to unique integers
    integer, dimension(:), allocatable, intent(out) :: ind
    !> optional index array that says which unique each element in the original array correspond to
    integer, dimension(:), intent(out), optional :: redind

    integer, dimension(:,:), allocatable :: dum
    integer, dimension(:), allocatable :: di
    integer :: i,j,k
    logical :: newval

    allocate(dum(size(a,1),size(a,2)))
    allocate(di(size(a,2)))
    if ( present(redind) ) then
        redind=-1
    endif

    ! Initialize counter
    dum=minval(a)-1
    di=-1
    k=0
    ! Start checking list
    vl: do i=1,size(a,2)
        newval=.true.
        ! Check if this line is in the list of unique
        do j=1,k
            if ( sum(abs(a(:,i)-dum(:,j))) .eq. 0 ) then
                newval=.false.
                if ( present(redind) ) redind(i)=j
                cycle vl
            endif
        enddo
        ! this means it's a new entry, add it
        if ( newval ) then
            k=k+1
            dum(:,k)=a(:,i)
            di(k)=i
            if ( present(redind) ) redind(i)=k
        endif
    enddo vl

    ! And return the unique
    allocate(ind(k))
    ind(1:k)=di(1:k)
    deallocate(dum,di)
end subroutine



end submodule
