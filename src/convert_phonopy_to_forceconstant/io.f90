#include "precompilerdefinitions"
module io
! Read files and do slight massaging
use konstanter, only: flyt, lo_forceconstant_2nd_eVA_to_HartreeBohr, lo_sqtol, lo_exitcode_io, lo_huge
use gottochblandat, only: open_file, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use type_symmetryoperation, only: lo_operate_on_secondorder_tensor
implicit none

private
public :: read_new_phonopy_forceconstant
public :: read_phonopy_forceconstant
public :: read_phonopy_born_file

contains

!> Read the phonopy forceconstant from file
subroutine read_new_phonopy_forceconstant(ss, filename, phonopy_fc, ssind)
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> filename
    character(len=*), intent(in) :: filename
    !> forceconstant
    real(flyt), dimension(:, :, :, :), allocatable, intent(out) :: phonopy_fc
    !> which supercell indices do I have forceconstants for?
    integer, dimension(:), allocatable, intent(out) :: ssind

    real(flyt), dimension(3) :: v0
    integer :: u, i, j, k, ii, jj, kk
    integer :: nuc

    nuc = maxval(ss%info%index_in_unitcell)
    lo_allocate(phonopy_fc(3, 3, ss%na, nuc))
    lo_allocate(ssind(nuc))
    ssind = -1
    phonopy_fc = lo_huge
    ! Now parse magic file
    u = open_file('in', trim(filename))
    read (u, *) ii, jj
    if (ii .ne. nuc) then
        call lo_stop_gracefully(['Thinking is hard'], lo_exitcode_io, __FILE__, __LINE__)
    end if
    if (jj .ne. ss%na) then
        call lo_stop_gracefully(['Thinking is hard'], lo_exitcode_io, __FILE__, __LINE__)
    end if
    do i = 1, nuc
    do j = 1, ss%na
        read (u, *) ii, jj
        kk = ss%info%index_in_unitcell(ii)
        ! Store index
        ssind(kk) = ii
        ! Store forceconstant
        do k = 1, 3
            read (u, *) v0
            phonopy_fc(:, k, jj, kk) = v0*lo_forceconstant_2nd_eVA_to_HartreeBohr
        end do
    end do
    end do
    close (u)
end subroutine

!> Read the phonopy forceconstant from file
subroutine read_phonopy_forceconstant(ss, filename, phonopy_fc)
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> filename
    character(len=*), intent(in) :: filename
    !> forceconstant
    real(flyt), dimension(:, :, :, :), allocatable, intent(out) :: phonopy_fc

    real(flyt), dimension(3) :: v0
    integer :: u, i, j, k, ii, jj, n1, n2

    lo_allocate(phonopy_fc(3, 3, ss%na, ss%na))
    phonopy_fc = 0.0_flyt
    u = open_file('in', trim(filename))
    read (u, *) n1, n2
    do i = 1, ss%na
    do j = 1, ss%na
        read (u, *) ii, jj
        do k = 1, 3
            read (u, *) v0
            phonopy_fc(:, k, jj, ii) = v0*lo_forceconstant_2nd_eVA_to_HartreeBohr
        end do
    end do
    end do
    close (u)
end subroutine

!> Read phonopy Born effective charges from file
subroutine read_phonopy_born_file(uc, filename, born_charges, dielectric_tensor)
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> filename
    character(len=*), intent(in) :: filename
    !> forceconstant
    real(flyt), dimension(:, :, :), allocatable, intent(out) :: born_charges
    !> dielectric tensor
    real(flyt), dimension(3, 3), intent(out) :: dielectric_tensor

    character(len=1000) :: ds1
    character(len=10), dimension(6) :: ds2
    integer, dimension(:), allocatable :: di
    integer :: u, i, l, a1, a2, o

    lo_allocate(born_charges(3, 3, uc%na))
    born_charges = 0.0_flyt

    ! Read raw born charges from file
    u = open_file('in', trim(filename))
    ! First grab the list which atoms the Born-charges are defined for
    read (u, '(A)') ds1
    i = lo_count_words_in_string(trim(adjustl(ds1))) - 6
    lo_allocate(di(i))
    read (ds1, *) ds2, di
    ! read dielectric tensor
    read (u, *) dielectric_tensor
    do i = 1, size(di)
        read (u, *) born_charges(:, :, di(i))
    end do
    close (u)

    !@TODO Figure out if I have to transpose them.

    ! Fill in the missing born charges
    do a1 = 1, uc%na
        if (sum(abs(born_charges(:, :, a1))) .gt. lo_sqtol) cycle
        ! ok found an empty
        l = 0
        a2 = 0
        ol1: do o = 1, uc%sym%n
        do i = 1, uc%na
            if (uc%sym%op(o)%fmap(i) .ne. a1) cycle
            if (sum(abs(born_charges(:, :, i))) .lt. lo_sqtol) cycle
            l = o
            a2 = i
            exit ol1
        end do
        end do ol1
        born_charges(:, :, a1) = lo_operate_on_secondorder_tensor(uc%sym%op(l), born_charges(:, :, a2))
    end do
end subroutine

!> Count the number of words in a string. Not smart or safe.
function lo_count_words_in_string(s) result(nw)
    !> string
    character(len=*), intent(in) :: s
    !> number of words
    integer :: nw

    integer, parameter :: n_max_words = 1000
    integer :: i, j, k

    i = 0
    j = 0
    l1: do
        if (j .gt. n_max_words) exit l1
        j = j + 1
        if (s(j:j) .ne. ' ') then
            i = i + 1
            if (j .eq. len(s)) exit l1
            do k = 1, n_max_words
                if (s(j + 1:j + 1) .eq. ' ') exit
                if (s(j + 1:j + 1) .ne. ' ') j = j + 1
                if (j .eq. len(s)) exit l1
            end do
        end if
    end do l1
    nw = i
end function

end module
