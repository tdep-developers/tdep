#include "precompilerdefinitions"
submodule (gottochblandat) gottochblandat_tensors
implicit none

contains

!> generate all tensor transposition things, as simple arrays as well as matrices that operate on flattened tensors
module subroutine lo_return_tensor_transpositions(trm_pair,trm_triplet,trm_quartet,prm_pair,prm_triplet,prm_quartet)
    real(flyt), dimension(9,9,2), intent(out), optional :: trm_pair
    real(flyt), dimension(27,27,6), intent(out), optional :: trm_triplet
    real(flyt), dimension(81,81,24), intent(out), optional :: trm_quartet
    integer, dimension(2,2), intent(out), optional :: prm_pair
    integer, dimension(3,6), intent(out), optional :: prm_triplet
    integer, dimension(4,24), intent(out), optional :: prm_quartet

    integer, dimension(:,:), allocatable :: pp,pt,pq
    integer :: d2(2),dd2(2),d3(3),dd3(3),d4(4),dd4(4)
    integer :: i1,i2,i3,i4,j1,j2,j3,j4,l,prm

    ! Start by generating the permutations
    call lo_permutations(pp,2)
    call lo_permutations(pt,3)
    call lo_permutations(pq,4)
    !prm_pair=pp
    !prm_triplet=pt
    !prm_quartet=pq

    if ( present(prm_pair) ) prm_pair=pp
    if ( present(prm_triplet) ) prm_triplet=pt
    if ( present(prm_quartet) ) prm_quartet=pq

    if ( present(trm_pair) ) then
        trm_pair=0.0_flyt
        l=0
        do i1=1,3
        do i2=1,3
            l=l+1
            d2=[i1,i2]
            do prm=1,size(pp,2)
                dd2=d2(pp(:,prm))
                j1=dd2(2)
                j2=dd2(1)
                trm_pair(j1 + (j2-1)*3 , l, prm)=1.0_flyt
            enddo
        enddo
        enddo
    endif
    if ( present(trm_triplet) ) then
        trm_triplet=0.0_flyt
        l=0
        do i1=1,3
        do i2=1,3
        do i3=1,3
            l=l+1
            d3=[i1,i2,i3]
            do prm=1,size(pt,2)
                dd3=d3(pt(:,prm))
                j1=dd3(3)
                j2=dd3(2)
                j3=dd3(1)
                trm_triplet(j1 + (j2-1)*3 + (j3-1)*9 , l, prm)=1.0_flyt
            enddo
        enddo
        enddo
        enddo
    endif
    if ( present(trm_quartet) ) then
        trm_quartet=0.0_flyt
        l=0
        do i1=1,3
        do i2=1,3
        do i3=1,3
        do i4=1,3
            l=l+1
            d4=[i1,i2,i3,i4]
            do prm=1,size(pq,2)
                dd4=d4(pq(:,prm))
                j1=dd4(4)
                j2=dd4(3)
                j3=dd4(2)
                j4=dd4(1)
                trm_quartet(j1 + (j2-1)*3 + (j3-1)*9 + (j4-1)*27, l, prm)=1.0_flyt
            enddo
        enddo
        enddo
        enddo
        enddo
    endif
end subroutine

!> routines that permute tensors
module pure function lo_transpose_2tensor(m,perm) result(pm)
    real(flyt), dimension(3,3), intent(in) :: m
    integer, dimension(2), intent(in) :: perm
    real(flyt), dimension(3,3) :: pm

    integer :: i,j
    integer, dimension(2) :: d
    do i=1,3
    do j=1,3
        d=[i,j]
        d=d(perm)
        pm(d(1),d(2))=m(i,j)
    enddo
    enddo
end function
module pure function lo_transpose_3tensor(m,perm) result(pm)
    real(flyt), dimension(3,3,3), intent(in) :: m
    integer, dimension(3), intent(in) :: perm
    real(flyt), dimension(3,3,3) :: pm

    integer :: i,j,k
    integer, dimension(3) :: d
    do i=1,3
    do j=1,3
    do k=1,3
        d=[i,j,k]
        d=d(perm)
        pm(d(1),d(2),d(3))=m(i,j,k)
    enddo
    enddo
    enddo
end function
module pure function lo_transpose_4tensor(m,perm) result(pm)
    real(flyt), dimension(3,3,3,3), intent(in) :: m
    integer, dimension(4), intent(in) :: perm
    real(flyt), dimension(3,3,3,3) :: pm

    integer :: i,j,k,l
    integer, dimension(4) :: d
    do i=1,3
    do j=1,3
    do k=1,3
    do l=1,3
        d=[i,j,k,l]
        d=d(perm)
        pm(d(1),d(2),d(3),d(4))=m(i,j,k,l)
    enddo
    enddo
    enddo
    enddo
end function

!> routines to flatten tensors
module pure function lo_flatten_2tensor(m) result(fm)
    real(flyt), dimension(3,3), intent(in) :: m
    real(flyt), dimension(9) :: fm

    integer :: i,j,l
    l=0
    do i=1,3
    do j=1,3
        l=l+1
        fm(l)=m(i,j)
    enddo
    enddo
end function
module pure function lo_flatten_3tensor(m) result(fm)
    real(flyt), dimension(3,3,3), intent(in) :: m
    real(flyt), dimension(27) :: fm

    integer :: i,j,k,l
    l=0
    do i=1,3
    do j=1,3
    do k=1,3
        l=l+1
        fm(l)=m(i,j,k)
    enddo
    enddo
    enddo
end function
module pure function lo_flatten_4tensor(m) result(fm)
    real(flyt), dimension(3,3,3,3), intent(in) :: m
    real(flyt), dimension(81) :: fm

    integer :: i,j,k,l,ll
    ll=0
    do i=1,3
    do j=1,3
    do k=1,3
    do l=1,3
        ll=ll+1
        fm(ll)=m(i,j,k,l)
    enddo
    enddo
    enddo
    enddo
end function
!> and unflatten them
module pure function lo_unflatten_2tensor(fm) result(m)
    real(flyt), dimension(9), intent(in) :: fm
    real(flyt), dimension(3,3) :: m

    integer :: i,j,l
    l=0
    do i=1,3
    do j=1,3
        l=l+1
        m(i,j)=fm(l)
    enddo
    enddo
end function
module pure function lo_unflatten_3tensor(fm) result(m)
    real(flyt), dimension(27), intent(in) :: fm
    real(flyt), dimension(3,3,3) :: m

    integer :: i,j,k,l
    l=0
    do i=1,3
    do j=1,3
    do k=1,3
        l=l+1
        m(i,j,k)=fm(l)
    enddo
    enddo
    enddo
end function
module pure function lo_unflatten_4tensor(fm) result(m)
    real(flyt), dimension(81), intent(in) :: fm
    real(flyt), dimension(3,3,3,3) :: m

    integer :: i,j,k,l,ll
    ll=0
    do i=1,3
    do j=1,3
    do k=1,3
    do l=1,3
        ll=ll+1
        m(i,j,k,l)=fm(ll)
    enddo
    enddo
    enddo
    enddo
end function

end submodule
