module mean_square_displacement

use konstanter, only: r8, lo_status, lo_bohr_to_A, lo_time_au_to_fs
use dump_data, only: lo_dump_gnuplot_2d_real
use gottochblandat, only: lo_sqnorm, tochar
use type_crystalstructure, only: lo_crystalstructure
use type_mdsim, only: lo_mdsim
use mpi_wrappers, only: lo_mpi_helper
use hdf5_wrappers, only: lo_hdf5_helper, lo_h5_store_data, lo_h5_store_attribute
use pairmapping, only: lo_pairmapping

implicit none
private

public :: lo_mean_square_displacement

!> symmetry decomposed pair distribution function
type lo_mean_square_displacement
    !> number of unique atoms
    integer :: na
    !> number of timesteps
    integer :: nt
    !> number of delta-times
    integer :: ndt
    !> time-axis for raw msd
    real(r8), dimension(:), allocatable :: x1, x2
    !> msd per atom
    real(r8), dimension(:, :), allocatable :: y1, y2
contains
    !> Calculate the mean square displacement
    procedure :: generate
    !> Write it to file
    procedure :: write_to_hdf5
    procedure :: write_to_plaintext
end type

contains

!> calculate the mean square displacement
subroutine generate(msd, uc, ss, sim, mw)
    !> mean square displacement
    class(lo_mean_square_displacement), intent(out) :: msd
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> simulation data
    type(lo_mdsim), intent(in) :: sim
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw

    real(r8), dimension(:, :), allocatable :: msd1
    real(r8), dimension(:, :), allocatable :: r0
    real(r8), dimension(3) :: v0
    integer, dimension(:), allocatable :: dum
    integer :: i, j, t

    ! ! Get the multiplicity of each atom
    ! allocate(dum(pm%na))
    ! dum=0
    ! do i=1,pm%nass
    !     j=pm%ssatom(i)%unique_atom
    !     dum(j)=dum(j)+1
    ! enddo

    ! Reference position per atom
    allocate (r0(3, ss%na))
    r0 = ss%r
    do i = 1, ss%na
        v0 = sim%r_npbc(:, i, 1) - r0(:, i)
        do j = 1, 3
            if (v0(j) .gt. 0.5_r8) r0(j, i) = r0(j, i) + 1.0_r8
            if (v0(j) .lt. -0.5_r8) r0(j, i) = r0(j, i) - 1.0_r8
        end do
    end do

    ! Total mean square displacement for each atom, just straight calculation
    allocate (msd1(sim%nt, sim%na))
    msd1 = 0.0_r8
    do t = 1, sim%nt
        if (mod(t, mw%n) .ne. mw%r) cycle
        do i = 1, ss%na
            msd1(t, i) = lo_sqnorm(ss%fractional_to_cartesian(sim%r_npbc(:, i, t) - r0(:, i)))
        end do
    end do
    call mw%allreduce('sum', msd1)

    ! Group this together by symmetry
    msd%na = uc%sym%n_irreducible_atom
    msd%nt = sim%nt
    allocate (dum(msd%na))
    allocate (msd%x1(msd%nt))
    allocate (msd%y1(msd%nt, msd%na))
    msd%x1 = 0.0_r8
    msd%y1 = 0.0_r8
    dum = 0
    do i = 1, ss%na
        j = ss%info%index_in_unitcell(i)
        j = uc%sym%all_to_irr(j)
        dum(j) = dum(j) + 1
    end do

    do i = 1, ss%na
        j = ss%info%index_in_unitcell(i)
        j = uc%sym%all_to_irr(j)
        msd%y1(:, j) = msd%y1(:, j) + msd1(:, i)/dum(j)
    end do
    do i = 1, sim%nt
        msd%x1(i) = (i - 1)*sim%timestep
    end do
end subroutine

!> Dump everything to file
subroutine write_to_hdf5(msd)
    !> mean square displacement
    class(lo_mean_square_displacement), intent(in) :: msd

    type(lo_hdf5_helper) :: h5
    real(r8), dimension(:), allocatable :: dr
    integer :: a1
    character(len=1000) :: filename, dname

    filename = 'outfile.mean_square_displacement.hdf5'
    allocate (dr(msd%nt))

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', 'outfile.mean_square_displacement.hdf5')

    call lo_h5_store_attribute(msd%na, h5%file_id, 'number_unique_atoms', lo_status)
    ! store the time axis
    dr = msd%x1*lo_time_au_to_fs
    call lo_h5_store_data(dr, h5%file_id, 'time_axis', enhet='fs')
    ! store the total MSD
    dr = 0.0_r8
    do a1 = 1, msd%na
        dr = dr + msd%y1(:, a1)
    end do
    dr = dr*lo_bohr_to_A**2
    call lo_h5_store_data(dr, h5%file_id, 'mean_square_displacement', enhet='Angstrom^2')
    ! store the projected MSD per atom
    if (msd%na .gt. 1) then
        do a1 = 1, msd%na
            dr = msd%y1(:, a1)
            dname = 'mean_square_displacement_atom_'//tochar(a1)
            call lo_h5_store_data(dr, h5%file_id, trim(dname), enhet='Angstrom^2')
        end do
    end if
    deallocate (dr)

    call h5%close_file()
    call h5%destroy()
end subroutine

!> Dump it to a plaintext file
subroutine write_to_plaintext(msd, uc)
    !> mean square displacement
    class(lo_mean_square_displacement), intent(in) :: msd
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc

    integer :: a1, i, j
    character(len=1000), dimension(:), allocatable :: legend
    real(r8), dimension(:, :), allocatable :: dum

    ! Only one kind of atom
    if (msd%na .eq. 1) then
        ! Only one kind of atom
        allocate (dum(2, msd%nt))
        dum(1, :) = msd%x1*lo_time_au_to_fs
        dum(2, :) = msd%y1(:, 1)*(lo_bohr_to_A**2)
        call lo_dump_gnuplot_2d_real(dum, 'outfile.mean_square_displacement', xlabel='Simulation time (fs)', ylabel='Mean square displacement (Å^2)')
        deallocate (dum)
    else
        ! More than one kind of atom
        allocate (dum(2 + msd%na, msd%nt))
        ! allocate(legend(msd%na+1))
        ! legend(1)='Total'
        dum = 0.0_r8
        dum(1, :) = msd%x1*lo_time_au_to_fs
        do a1 = 1, msd%na
            dum(2, :) = dum(2, :) + msd%y1(:, a1)*(lo_bohr_to_A**2)
            dum(2 + a1, :) = msd%y1(:, a1)*(lo_bohr_to_A**2)
            ! get the annoying legend reasonable
            !j=pm%atom(a1)%unitcell_indices(1)
            ! legend(a1+1)=trim(uc%atomic_symbol(uc%species(j)))//':'
            ! do i=1,pm%atom(a1)%n_uc_atoms
            !     j=pm%atom(a1)%unitcell_indices(i)
            !     legend(a1+1)=trim(legend(a1+1))//' atom '//tochar(j)
            ! enddo
        end do
        !call lo_dump_gnuplot_2d_real(dum,'outfile.mean_square_displacement',xlabel='Simulation time (fs)',ylabel='Mean square displacement (Å^2)',legend=legend)
        call lo_dump_gnuplot_2d_real(dum, 'outfile.mean_square_displacement', xlabel='Simulation time (fs)', ylabel='Mean square displacement (Å^2)')
        deallocate (dum)
    end if
end subroutine

end module
