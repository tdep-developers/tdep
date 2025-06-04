#include "precompilerdefinitions"
program remap_forceconstant
!!{!src/remap_forceconstant/manual.md!}
use konstanter, only: flyt, lo_sqtol, lo_freqtol
use gottochblandat, only: lo_invert3x3matrix, open_file, lo_does_file_exist, lo_frobnorm
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use options, only: lo_opts
use lo_memtracker, only: lo_mem_helper


implicit none

type(lo_opts) :: opts
type(lo_crystalstructure) :: uc, nc
type(lo_forceconstant_secondorder) :: fc, fcn
type(lo_forceconstant_thirdorder) :: fct, fctn
type(lo_forceconstant_fourthorder) :: fcf !,fcfn
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
integer :: i,j,k,l,u,ii
real(flyt), dimension(3) :: v0
real(flyt), dimension(3, 3) :: tm
real(flyt), dimension(:, :, :), allocatable :: dfc1, dfc2
!
call opts%parse()
call mem%init()
call mw%init()
! Start memory tracker
call mem%init()

! read files
call uc%readfromfile('infile.ucposcar')
call nc%readfromfile('infile.newposcar')
write (*, *) '... read structures'

if (lo_does_file_exist('infile.forceconstant')) then
    call fc%readfromfile(uc, 'infile.forceconstant', mem, 1)
    write (*, *) '... read secondorder forceconstant'
end if
!
if (lo_does_file_exist('infile.forceconstant_thirdorder')) then
    call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
    write (*, *) '... read thirdorder forceconstant'
end if
!
if (lo_does_file_exist('infile.forceconstant_fourthorder')) then
    call fcf%readfromfile(uc, 'infile.forceconstant_fourthorder')
    write (*, *) '... read fourthorder forceconstant'
end if

! check for a transformation matrix
if (lo_does_file_exist('infile.transformation')) then
    write (*, *) 'Trouble: I should fix this'
    stop
    u = open_file('in', 'infile.transformation')
    do i = 1, 3
        read (u, *) tm(i, :)
    end do
    write (*, *) '... read a transformation matrix:'
    do i = 1, 3
        write (*, "(3(2X,F11.8))") tm(i, :)
    end do
else
    ! use the identity transformation
    tm(1, :) = [1.0_flyt, 0.0_flyt, 0.0_flyt]
    tm(2, :) = [0.0_flyt, 1.0_flyt, 0.0_flyt]
    tm(3, :) = [0.0_flyt, 0.0_flyt, 1.0_flyt]
end if

! Transform the original cell
tm = transpose(tm)
uc%latticevectors = matmul(uc%latticevectors, tm)
uc%inv_latticevectors = lo_invert3x3matrix(uc%latticevectors)
! need to transpose here, otherwise it does not work
tm = transpose(tm)

! Transform and remap the forceconstants
if (lo_does_file_exist('infile.forceconstant')) then
    ! transform
    do i = 1, fc%na
        do j = 1, fc%atom(i)%n
            fc%atom(i)%pair(j)%lv1 = matmul(tm, fc%atom(i)%pair(j)%lv1)
            fc%atom(i)%pair(j)%lv2 = matmul(tm, fc%atom(i)%pair(j)%lv2)
            fc%atom(i)%pair(j)%r = matmul(tm, fc%atom(i)%pair(j)%r)
            fc%atom(i)%pair(j)%m = rotate(tm, fc%atom(i)%pair(j)%m)
        end do
    end do
    ! remap
    write (*, *) 'Remapping, cutoff:', fc%cutoff
    call nc%classify('supercell', uc)
    call fc%remap(uc, nc, fcn)
    ! print forceconstant file
    call fcn%writetofile(nc, 'outfile.forceconstant_remapped')
    write (*, *) 'remapped the second order force constants'
    ! print this in LAMMPS-compliant format.
    if (opts%lammps) then
        ! print the output in a lammps-friendly way,
        ! first we write the equilibrium positions:
        u = open_file('out', 'outfile.lammps_eqpos')
        do i = 1, nc%na
            !v0=nc%r(:,i)
            !do j=1,3
            !    if ( abs(v0(j)) .lt. lo_sqtol ) v0(j)=1.0_flyt
            !enddo
            !write(u,*) matmul(nc%latticevectors,v0)
            write (u, *) nc%rcart(:, i) !-matmul(nc%latticevectors,[1_flyt,1_flyt,1_flyt])
        end do
        close (u)
        ! second is a list of the forceconstants, only the small cell ones
        l = 0
        do i = 1, fc%na
            l = l + fc%atom(i)%n
        end do
        lo_allocate(dfc1(3, 3, l))

        l = 0
        do i = 1, fc%na
        do j = 1, fc%atom(i)%n
            l = l + 1
            dfc1(:, :, l) = fc%atom(i)%pair(j)%m
        end do
        end do

        u = open_file('out', 'outfile.lammps_pairfc')
        ! number of distinct forceconstants
        write (u, *) l
        ! all the forceconstants
        do i = 1, l
        do j = 1, 3
            write (u, *) dfc1(j, :, i)
        end do
        end do
        close (u)

        ! next print a map of all pairs.
        u = open_file('out', 'outfile.lammps_pairmap')
        ! how many atoms, cutoff
        write (u, *) fcn%na, fcn%cutoff
        ! the number of neighbours for each atom
        do i = 1, fcn%na
            write (u, *) i, fcn%atom(i)%n
        end do
        ! now massive loop over all pairs
        do i = 1, fcn%na
        do j = 1, fcn%atom(i)%n
            ! index to atom 2
            k = fcn%atom(i)%pair(j)%i2
            ! vector to atom 2
            v0 = fcn%atom(i)%pair(j)%r
            ! which forceconstant?
            ii = -1
            do l = 1, size(dfc1, 3)
                if (lo_frobnorm(dfc1(:, :, l) - fcn%atom(i)%pair(j)%m) .lt. lo_sqtol) then
                    ii = l
                    exit
                end if
            end do
            write (u, "(4(1X,I5),3(1X,F18.10))") i, j, k, ii, v0
        end do
        end do
        close (u)

        ! get the harmonic stuff
        call fcn%initialize_cell(nc, uc, fc, 0.0_flyt, .true., .false., lo_freqtol, mw=mw)
        u = open_file('out', 'outfile.commensurate_mode_frequencies')
        do i = 1, fcn%na*3
            write (u, *) fcn%omega(i)
        end do
        close (u)
    end if
end if

if (lo_does_file_exist('infile.forceconstant_thirdorder')) then
    ! transform
    do i = 1, fct%na
        do j = 1, fct%atom(i)%n
            fct%atom(i)%triplet(j)%lv1 = matmul(tm, fct%atom(i)%triplet(j)%lv1)
            fct%atom(i)%triplet(j)%lv2 = matmul(tm, fct%atom(i)%triplet(j)%lv2)
            fct%atom(i)%triplet(j)%lv3 = matmul(tm, fct%atom(i)%triplet(j)%lv3)
            fct%atom(i)%triplet(j)%rv1 = matmul(tm, fct%atom(i)%triplet(j)%rv1)
            fct%atom(i)%triplet(j)%rv2 = matmul(tm, fct%atom(i)%triplet(j)%rv2)
            fct%atom(i)%triplet(j)%rv3 = matmul(tm, fct%atom(i)%triplet(j)%rv3)
            !
            fct%atom(i)%triplet(j)%m = rotate3(tm, fct%atom(i)%triplet(j)%m)
        end do
    end do
    ! remap
    call fct%remap(uc, nc, fctn)
    call fctn%writetofile(nc, 'outfile.forceconstant_thirdorder_remapped')
    write (*, *) 'remapped the third order force constants'
end if

!if ( lo_does_file_exist('infile.forceconstant_fourthorder') ) then
!    ! transform
!    do i=1,fcf%na
!        do j=1,fcf%atom(i)%n
!            fcf%atom(i)%quartet(j)%lv1=matmul(tm,fcf%atom(i)%quartet(j)%lv1)
!            fcf%atom(i)%quartet(j)%lv2=matmul(tm,fcf%atom(i)%quartet(j)%lv2)
!            fcf%atom(i)%quartet(j)%lv3=matmul(tm,fcf%atom(i)%quartet(j)%lv3)
!            fcf%atom(i)%quartet(j)%lv4=matmul(tm,fcf%atom(i)%quartet(j)%lv4)
!            fcf%atom(i)%quartet(j)%rv1=matmul(tm,fcf%atom(i)%quartet(j)%rv1)
!            fcf%atom(i)%quartet(j)%rv2=matmul(tm,fcf%atom(i)%quartet(j)%rv2)
!            fcf%atom(i)%quartet(j)%rv3=matmul(tm,fcf%atom(i)%quartet(j)%rv3)
!            fcf%atom(i)%quartet(j)%rv4=matmul(tm,fcf%atom(i)%quartet(j)%rv4)
!            !
!            fcf%atom(i)%quartet(j)%m=rotate4(tm,fcf%atom(i)%quartet(j)%m)
!            !
!        enddo
!    enddo
!    ! remap
!    call fcf%remap(uc,nc,fcfn)
!    call fcfn%writetofile(nc,'outfile.forceconstant_fourthorder_remapped')
!    write(*,*) 'remapped the fourth order force constants'
!endif

contains
!
function rotate(op, m) result(n)
    real(flyt), dimension(3, 3), intent(in) :: op, m
    real(flyt), dimension(3, 3) :: n
    !
    integer :: i, j, ii, jj
    !
    n = 0.0_flyt
    do j = 1, 3
    do i = 1, 3
        do jj = 1, 3
        do ii = 1, 3
            n(i, j) = n(i, j) + m(ii, jj)*op(i, ii)*op(j, jj)
        end do
        end do
    end do
    end do
    !
end function
!
function rotate3(op, m) result(n)
    real(flyt), dimension(3, 3, 3), intent(in) :: m
    real(flyt), dimension(3, 3), intent(in) :: op
    real(flyt), dimension(3, 3, 3) :: n
    !
    integer :: i, j, k, ii, jj, kk
    !
    n = 0.0_flyt
    do i = 1, 3
    do j = 1, 3
    do k = 1, 3
        do ii = 1, 3
        do jj = 1, 3
        do kk = 1, 3
            n(i, j, k) = n(i, j, k) + m(ii, jj, kk)*op(i, ii)*op(j, jj)*op(k, kk)
        end do
        end do
        end do
    end do
    end do
    end do
    !
end function
!
function rotate4(op, m) result(n)
    real(flyt), dimension(3, 3, 3, 3), intent(in) :: m
    real(flyt), dimension(3, 3), intent(in) :: op
    real(flyt), dimension(3, 3, 3, 3) :: n
    !
    integer :: i, j, k, l, ii, jj, kk, ll
    !
    n = 0.0_flyt
    do i = 1, 3
    do j = 1, 3
    do k = 1, 3
    do l = 1, 3
        do ii = 1, 3
        do jj = 1, 3
        do kk = 1, 3
        do ll = 1, 3
            n(i, j, k, l) = n(i, j, k, l) + m(ii, jj, kk, ll)*op(i, ii)*op(j, jj)*op(k, kk)*op(l, ll)
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do
    !
end function
!
end program
