#include "precompilerdefinitions"
module vectordist
use konstanter, only: flyt, lo_huge, lo_hugeint, lo_status
use gottochblandat, only: lo_sqnorm, walltime, lo_progressbar_init, lo_progressbar, tochar, lo_linspace
use pairmapping, only: lo_pairmapping
use type_crystalstructure, only: lo_crystalstructure
use type_mdsim, only: lo_mdsim
use hdf5_wrappers, only: lo_h5_store_data, lo_h5_store_attribute, HID_T, H5F_ACC_TRUNC_F, &
                         h5close_f, h5open_f, h5fclose_f, h5fopen_f, h5fcreate_f, h5dclose_f, h5dopen_f
implicit none
private
public :: lo_vectordist

!> representation of all the atoms in a single coordination shell
type lo_vd_atom_shell
    !> what are the limits of the histogram?
    real(flyt), dimension(3) :: maxcoord = -lo_huge
    real(flyt), dimension(3) :: mincoord = -lo_huge
    !> and the same in all Cartesian directions
    real(flyt) :: maxv = -lo_huge
    !> the bin centers
    real(flyt), dimension(:), allocatable :: bincenters
    !> the actual histogram
    real(flyt), dimension(:, :, :), allocatable :: histogram
end type

!> representation of the pairs originating from this atom
type lo_vd_atom
    !> each atom will have a certain number of coordination shells
    integer :: nshell = -lo_hugeint
    !> but it will also have a certain number of neighbours
    integer :: npair = -lo_hugeint
    !> and the actual shells
    type(lo_vd_atom_shell), dimension(:), allocatable :: shell
end type

!> vector distribution
type lo_vectordist
    ! number of bins in the histogram
    integer :: nh = -lo_hugeint
    ! axis that things are rotated towards
    real(flyt) :: axis = -lo_huge
    ! largest distance in +- something
    real(flyt) :: hmax = -lo_huge
    ! number of unique atoms
    integer :: na = -lo_hugeint
    ! one distribution per unique atom
    type(lo_vd_atom), dimension(:), allocatable :: atom
contains
    !> Bin the timesteps
    procedure :: generate
    !> dump it
    procedure :: write_to_hdf5
end type

contains

!> Dump everything to file
subroutine write_to_hdf5(vd, pm)
    class(lo_vectordist), intent(in) :: vd
    class(lo_pairmapping), intent(in) :: pm

    integer(HID_T) :: file_id, dset_id
    integer :: i, a1, sh
    character(len=1000) :: filename, dname

    filename = 'outfile.vector_distribution.hdf5'

    call h5open_f(lo_status)
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, lo_status)

    ! Some general info
    call lo_h5_store_attribute(vd%na, file_id, 'number_unique_atoms', lo_status)
    do i = 1, vd%na
        call lo_h5_store_attribute(vd%atom(i)%nshell, file_id, 'number_shells_atom_'//tochar(i), lo_status)
    end do

    do a1 = 1, vd%na
    do sh = 1, vd%atom(a1)%nshell
        ! Store the histogram
        dname = 'distribution_atom_'//tochar(a1)//'_shell_'//tochar(sh)
        call lo_h5_store_data(vd%atom(a1)%shell(sh)%histogram, file_id, trim(dname))
        ! Set some metadata on the histogram
        call h5dopen_f(file_id, trim(dname), dset_id, lo_status)
        !
        call lo_h5_store_attribute(vd%atom(a1)%shell(sh)%bincenters, dset_id, 'bincenters', lo_status)
        call lo_h5_store_attribute(pm%atom(a1)%shell(sh)%prototype_vector, dset_id, 'bond_vector', lo_status)
        call lo_h5_store_attribute(pm%atom(a1)%shell(sh)%coordinate_transformation, dset_id, 'rotation_matrix', lo_status)
        call lo_h5_store_attribute('Angstrom', dset_id, 'unit', lo_status)
        !
        call h5dclose_f(dset_id, lo_status)
    end do
    end do

!    ! store the total
!    call lo_h5_store_data(pdf%x,file_id,'r_axis')
!        call h5dopen_f(file_id, 'r_axis', dset_id, lo_status)
!        call lo_h5_store_attribute('Angstrom',dset_id,'unit',lo_status)
!        call h5dclose_f(dset_id,lo_status)
!    call lo_h5_store_data(pdf%y,file_id,'radial_pair_distribution_function')
!
!    ! store the projected
!    call lo_h5_store_attribute(pdf%na,file_id,'number_unique_atoms',lo_status)
!    do a1=1,pdf%na
!        ! pack the data
!        if ( pdf%diffusion ) then
!            nb=size(pdf%atom(a1)%shell(1)%x,1)
!            lo_allocate(dumx(pdf%atom(a1)%nshell,nb))
!            lo_allocate(dumy(pdf%atom(a1)%nshell,nb))
!            do sh=1,pdf%atom(a1)%nshell
!                dumx(sh,:)=pdf%atom(a1)%shell(sh)%x
!                dumy(sh,:)=pdf%atom(a1)%shell(sh)%y
!            enddo
!        else
!            nb=size(pdf%atom(a1)%shell(1)%x,1)
!            lo_allocate(dumx(pdf%atom(a1)%nshell-1,nb))
!            lo_allocate(dumy(pdf%atom(a1)%nshell-1,nb))
!            do sh=2,pdf%atom(a1)%nshell
!                dumx(sh-1,:)=pdf%atom(a1)%shell(sh)%x
!                dumy(sh-1,:)=pdf%atom(a1)%shell(sh)%y
!            enddo
!        endif
!
!        dname='projected_r_axis_'//trim(int2char(a1))
!        call lo_h5_store_data(dumx,file_id,trim(dname))
!        dname='projected_pair_distribution_function_atom_'//trim(int2char(a1))
!        call lo_h5_store_data(dumy,file_id,trim(dname))
!
!        lo_deallocate(dumx)
!        lo_deallocate(dumy)
!    enddo

    call h5fclose_f(file_id, lo_status)
    call h5close_f(lo_status)
end subroutine

!> build the actual histograms
subroutine generate(vd, pm, sim, bintype, transform, nh)
    !> the histogram stuff
    class(lo_vectordist), intent(inout) :: vd
    !> pair mapping
    type(lo_pairmapping), intent(in) :: pm
    !> the simulation
    type(lo_mdsim), intent(in) :: sim
    !> how to bin it
    integer, intent(in) :: bintype
    !> should the coordinates be transformed
    logical, intent(in) :: transform
    !> how many bins?
    integer, intent(in) :: nh

    integer :: i, j, k
    integer :: a1, a2, sh, una, unsh, t, stride
    integer :: ii, jj, kk
    real(flyt) :: t0
    real(flyt), dimension(3) :: v0
    real(flyt), dimension(3, 3, 3) :: kernel

    ! Make some space
    vd%na = pm%na
    vd%nh = nh
    allocate (vd%atom(pm%na))
    do a1 = 1, pm%na
        vd%atom(a1)%nshell = pm%atom(a1)%nshell
        allocate (vd%atom(a1)%shell(vd%atom(a1)%nshell))
        do sh = 1, vd%atom(a1)%nshell
            lo_allocate(vd%atom(a1)%shell(sh)%bincenters(vd%nh))
            lo_allocate(vd%atom(a1)%shell(sh)%histogram(vd%nh, vd%nh, vd%nh))
            vd%atom(a1)%shell(sh)%histogram = 0.0_flyt
            ! initialize the limits
            vd%atom(a1)%shell(sh)%mincoord = lo_huge
            vd%atom(a1)%shell(sh)%maxcoord = -lo_huge
        end do
    end do

    ! Get the approximate max and min coordinates for the histogram. Check 100 random timesteps, should be ok.
    call lo_progressbar_init()
    t0 = walltime()
    stride = max(sim%nt/100, 1)
    do t = 1, sim%nt, stride
        ! loop over all atoms
        do a1 = 1, pm%nass
            ! unique atom
            una = pm%ssatom(a1)%unique_atom
            do i = 1, pm%ssatom(a1)%npair
                ! second atom in the pair
                a2 = pm%ssatom(a1)%pair(i)%i2
                ! I only care about displacements from equilibrium.
                if (a1 .ne. a2) then
                    v0 = sim%u(:, a2, t) - sim%u(:, a1, t)
                else
                    v0 = sim%u(:, a1, t)
                end if
                ! Which shell is it
                unsh = pm%ssatom(a1)%pair(i)%shell
                ! Transform it properly
                v0 = matmul(pm%ssatom(a1)%pair(i)%m, v0)
                if (transform) then
                    v0 = matmul(pm%atom(una)%shell(unsh)%coordinate_transformation, v0)
                end if
                ! Check limits
                do k = 1, 3
                    vd%atom(una)%shell(unsh)%mincoord(k) = min(vd%atom(una)%shell(unsh)%mincoord(k), v0(k))
                    vd%atom(una)%shell(unsh)%maxcoord(k) = max(vd%atom(una)%shell(unsh)%maxcoord(k), v0(k))
                end do
            end do
        end do
        call lo_progressbar(' ... checking limits of histograms', t, sim%nt)
    end do
    call lo_progressbar(' ... checking limits of histograms', sim%nt, sim%nt)

    ! Print the limits for troubleshooting
!    write(*,*) '... got limits',walltime()-t0
!    do a1=1,vd%na
!        do sh=1,vd%atom(a1)%nshell
!            write(*,*) 'Atom '//tochar(a1),' shell '//tochar(sh)
!            write(*,"(' ... lower bounds ',3(2X,F11.6))") vd%atom(a1)%shell(sh)%mincoord
!            write(*,"(' ... upper bounds ',3(2X,F11.6))") vd%atom(a1)%shell(sh)%maxcoord
!            vd%atom(a1)%shell(sh)%maxv=max( maxval(vd%atom(a1)%shell(sh)%maxcoord), -minval(vd%atom(a1)%shell(sh)%mincoord) )*1.0_flyt
!            write(*,"(' ... range +-',1(1X,F9.6))") vd%atom(a1)%shell(sh)%maxv
!        enddo
!    enddo

    ! Get the ranges for the histogram
    do a1 = 1, vd%na
        do sh = 1, vd%atom(a1)%nshell
            call lo_linspace(-vd%atom(a1)%shell(sh)%maxv, vd%atom(a1)%shell(sh)%maxv, vd%atom(a1)%shell(sh)%bincenters)
        end do
    end do

    ! Magic kernel
    do i = 1, 3
    do j = 1, 3
    do k = 1, 3
        v0 = [i - 2, j - 2, k - 2]*1.0_flyt
        kernel(i, j, k) = exp(-lo_sqnorm(v0*0.5_flyt))
    end do
    end do
    end do

    ! Actually bin everything
    t0 = walltime()
    call lo_progressbar_init()
    do t = 1, sim%nt
        ! loop over all atoms
        do a1 = 1, pm%nass
            ! unique atom
            una = pm%ssatom(a1)%unique_atom
            do i = 1, pm%ssatom(a1)%npair
                ! second atom in the pair
                a2 = pm%ssatom(a1)%pair(i)%i2
                ! I only care about displacements from equilibrium.
                if (a1 .ne. a2) then
                    v0 = sim%u(:, a2, t) - sim%u(:, a1, t)
                else
                    v0 = sim%u(:, a1, t)
                end if
                ! Which shell is it
                unsh = pm%ssatom(a1)%pair(i)%shell
                ! Transform it properly
                v0 = matmul(pm%ssatom(a1)%pair(i)%m, v0)
                ! and maybe a bit more
                if (transform) then
                    v0 = matmul(pm%atom(una)%shell(unsh)%coordinate_transformation, v0)
                end if

                ! Actually bin it
                select case (bintype)
                case (1)
                    ! Straight binning
                    ii = binindex(v0(1), vd%nh, vd%atom(una)%shell(unsh)%maxv)
                    jj = binindex(v0(2), vd%nh, vd%atom(una)%shell(unsh)%maxv)
                    kk = binindex(v0(3), vd%nh, vd%atom(una)%shell(unsh)%maxv)
                    if (ii .gt. 0 .and. ii .lt. vd%nh) then
                    if (jj .gt. 0 .and. jj .lt. vd%nh) then
                    if (kk .gt. 0 .and. kk .lt. vd%nh) then
                        vd%atom(una)%shell(unsh)%histogram(ii, jj, kk) = vd%atom(una)%shell(unsh)%histogram(ii, jj, kk) + 1.0_flyt
                    end if
                    end if
                    end if
                case (2)
                    ! add it like a gaussian
                    call add_to_histogram_with_proper_gaussian(v0, vd%atom(una)%shell(unsh)%histogram, &
                                                               vd%atom(una)%shell(unsh)%bincenters, vd%atom(una)%shell(unsh)%maxv, 4)
                case (3)
                    ! add like a gaussian but stupid
                    ii = binindex(v0(1), vd%nh, vd%atom(una)%shell(unsh)%maxv)
                    jj = binindex(v0(2), vd%nh, vd%atom(una)%shell(unsh)%maxv)
                    kk = binindex(v0(3), vd%nh, vd%atom(una)%shell(unsh)%maxv)
                    if (ii .gt. 1 .and. ii .lt. vd%nh - 1) then
                    if (jj .gt. 1 .and. jj .lt. vd%nh - 1) then
                    if (kk .gt. 1 .and. kk .lt. vd%nh - 1) then
                        vd%atom(una)%shell(unsh)%histogram(ii - 1:ii + 1, jj - 1:jj + 1, kk - 1:kk + 1) = &
                            vd%atom(una)%shell(unsh)%histogram(ii - 1:ii + 1, jj - 1:jj + 1, kk - 1:kk + 1) + kernel
                    end if
                    end if
                    end if
                end select
            end do
        end do
        call lo_progressbar(' ... binning steps', t, sim%nt, walltime() - t0)
    end do

    ! Normalise it to something. I just take the largest values as 1 for now.
    do a1 = 1, vd%na
        do sh = 1, vd%atom(a1)%nshell
            vd%atom(a1)%shell(sh)%histogram = vd%atom(a1)%shell(sh)%histogram/maxval(vd%atom(a1)%shell(sh)%histogram)
        end do
    end do

contains

    !> add things to the histogram with a slight gaussian smearing, to get subpixel resolution
    pure subroutine add_to_histogram_with_proper_gaussian(r, h, x, maxv, sigma)
        !> the coordinates
        real(flyt), dimension(3), intent(in) :: r
        !> the histogram
        real(flyt), dimension(:, :, :), intent(inout) :: h
        !> the bincenters
        real(flyt), dimension(:), intent(in) :: x
        !> the max coordinate
        real(flyt), intent(in) :: maxv
        !> the smearing, in some kind of integer units
        integer, intent(in) :: sigma
        !
        real(flyt), dimension(3) :: v0
        real(flyt) :: smear, f0
        integer :: n
        integer :: i, j, k, ii, jj, kk
        ! number of bins
        n = size(x, 1)
        smear = x(2) - x(1)
        smear = 2.0_flyt/smear/smear
        !
        i = binindex(r(1), n, maxv)
        j = binindex(r(2), n, maxv)
        k = binindex(r(3), n, maxv)
        !
        do ii = max(i - sigma, 1), min(i + sigma, n)
        do jj = max(j - sigma, 1), min(j + sigma, n)
        do kk = max(k - sigma, 1), min(k + sigma, n)
            v0 = [x(ii), x(jj), x(kk)] - r
            f0 = exp(-lo_sqnorm(v0)*smear)
            h(ii, jj, kk) = h(ii, jj, kk) + f0
        end do
        end do
        end do
    end subroutine

    !> convert coordinates to a binindex
    pure function binindex(r, n, maxv) result(i)
        !> the coordinate
        real(flyt), intent(in) :: r
        !> how many bins
        integer, intent(in) :: n
        !> max range in +-
        real(flyt), intent(in) :: maxv
        !> the index
        integer :: i
        !
        i = floor((r + maxv)*n*0.5_flyt/maxv) + 1
    end function

end subroutine

end module
