#include "precompilerdefinitions"
module type_forceconstant_fourthorder
use konstanter, only: flyt, lo_huge, lo_hugeint, lo_tol, lo_sqtol, lo_status, lo_bohr_to_A, lo_forceconstant_4th_HartreeBohr_to_eVA, &
                      lo_forceconstant_4th_eVa_to_HartreeBohr, lo_exitcode_symmetry
use gottochblandat, only: open_file, lo_sqnorm, lo_chop, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use type_distancetable, only: lo_distancetable
implicit none

private
public :: lo_forceconstant_fourthorder

type dumegv
    complex(flyt), dimension(3, 3, 3, 3) :: m = lo_huge
end type

!> one quartet
type lo_fc4_quartet
    ! index in the unit cell, for atom 1 2 3 4 in the quartet
    integer :: i1 = -lo_hugeint, i2 = -lo_hugeint, i3 = -lo_hugeint, i4 = -lo_hugeint
    ! lattice vector pointing to the unit cell of the atom
    real(flyt), dimension(3) :: lv1 = lo_huge, lv2 = lo_huge, lv3 = lo_huge, lv4 = lo_huge
    ! absolute vector between the atoms
    real(flyt), dimension(3) :: rv1 = lo_huge, rv2 = lo_huge, rv3 = lo_huge, rv4 = lo_huge
    ! force constant tensor
    real(flyt), dimension(3, 3, 3, 3) :: m = lo_huge
    ! force constant tensor with masses
    real(flyt), dimension(3, 3, 3, 3) :: mwm = lo_huge
end type

!> list of quartets
type lo_fc4_atom
    !> how many quartets
    integer :: n = -lo_hugeint
    !> quartets
    type(lo_fc4_quartet), dimension(:), allocatable :: quartet
end type

!> fourth-order forceconstant
type lo_forceconstant_fourthorder
    !> Need to know number of atoms/unitcell
    integer :: na = -lo_hugeint
    !> Cutoff radius, not really that necessary
    real(flyt) :: cutoff = -lo_huge
    !> Atoms, these contain the quartets
    type(lo_fc4_atom), allocatable, dimension(:) :: atom
    !> How much to talk
    integer :: verbosity = -lo_hugeint
contains
    !> read from file
    procedure :: readfromfile
    !> write to file
    procedure :: writetofile
    !> remap to another cell
    procedure :: remap
    !> potential energy
    procedure :: potential_energy
    !> forces
    procedure :: forces
    !> four-phonon scattering amplitude
    procedure :: scatteringamplitude
    !> destroy it
    procedure :: unallocate
#ifdef AGRESSIVE_SANITY
    !> size in memory, in bytes
    procedure :: size_in_mem => fc4_size_in_mem
#endif
end type

contains

#ifdef AGRESSIVE_SANITY
!> size in memory, in bytes
function fc4_size_in_mem(f) result(mem)
    !> forceconstant
    class(lo_forceconstant_fourthorder), intent(in) :: f
    !> size in memory, in bytes
    integer :: mem

    integer :: i, j

    mem = 0
    mem = mem + storage_size(f)
    if (allocated(f%atom)) then
        do i = 1, size(f%atom)
            mem = mem + storage_size(f%atom(i))
            if (allocated(f%atom(i)%quartet)) then
                do j = 1, size(f%atom(i)%quartet)
                    mem = mem + storage_size(f%atom(i)%quartet(j))
                end do
            end if
        end do
    end if
    mem = mem/8
end function
#endif

!> destroy the forceconstant
subroutine unallocate(fc)
    !> forceconstant
    class(lo_forceconstant_fourthorder), intent(inout) :: fc

    integer :: i

    if (allocated(fc%atom)) then
        do i = 1, size(fc%atom, 1)
            if (allocated(fc%atom(i)%quartet)) lo_deallocate(fc%atom(i)%quartet)
        end do
        lo_deallocate(fc%atom)
    end if
end subroutine

!> get the potential energy from a set of displacements
real(flyt) function potential_energy(fc, u)
    !> forceconstant
    class(lo_forceconstant_fourthorder), intent(in) :: fc
    !> displacements
    real(flyt), dimension(:, :), intent(in) :: u
    !
    real(flyt) :: e
    integer :: a1, a2, a3, a4, t, i, j, k, l
    !
    if (fc%na .ne. size(u, 2)) then
        write (*, *) 'Size mismatch between fourth order force constant and'
        write (*, *) 'displacement array'
        stop
    end if
    !
    e = 0.0_flyt
    do a1 = 1, fc%na
        do t = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%quartet(t)%i2
            a3 = fc%atom(a1)%quartet(t)%i3
            a4 = fc%atom(a1)%quartet(t)%i4
            do i = 1, 3
            do j = 1, 3
            do k = 1, 3
            do l = 1, 3
                e = e + fc%atom(a1)%quartet(t)%m(i, j, k, l)*u(i, a1)*u(j, a2)*u(k, a3)*u(l, a4)/24.0_flyt
            end do
            end do
            end do
            end do
        end do
    end do
    potential_energy = e
end function

subroutine forces(fc, u, f)
    class(lo_forceconstant_fourthorder), intent(in) :: fc
    real(flyt), dimension(:, :), intent(in) :: u
    real(flyt), dimension(:, :), intent(out) :: f
    !
    integer :: a1, a2, a3, a4, t, i, j, k, l
    !
    f = 0.0_flyt
    do a1 = 1, fc%na
        do t = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%quartet(t)%i2
            a3 = fc%atom(a1)%quartet(t)%i3
            a4 = fc%atom(a1)%quartet(t)%i4
            do i = 1, 3
            do j = 1, 3
            do k = 1, 3
            do l = 1, 3
                f(i, a1) = f(i, a1) - fc%atom(a1)%quartet(t)%m(i, j, k, l)*u(j, a2)*u(k, a3)*u(l, a4)/6.0_flyt
            end do
            end do
            end do
            end do
        end do
    end do
end subroutine

!> Map the fourth order force constant to a supercell
subroutine remap(fc, uc, ss, fcss)
    class(lo_forceconstant_fourthorder), intent(in) :: fc
    type(lo_crystalstructure), intent(in) :: ss, uc
    type(lo_forceconstant_fourthorder), intent(out) :: fcss
    !
    type(lo_distancetable) :: dt
    integer :: i, j, k, l
    integer :: uca, ii, kk, a1
    real(flyt) :: sqrc, rc

    ! Get the actual cutoff
    rc = 0.0_flyt
    do a1 = 1, uc%na
    do i = 1, fc%atom(a1)%n
        rc = max(rc, norm2(fc%atom(a1)%quartet(i)%rv2))
        rc = max(rc, norm2(fc%atom(a1)%quartet(i)%rv3))
        rc = max(rc, norm2(fc%atom(a1)%quartet(i)%rv4))
    end do
    end do

    ! Build a supercell neighbourlist
    call dt%generate(ss%r, ss%latticevectors, rc + 10*lo_tol, verbosity=0)

    ! Build the skeleton of the forceconstant
    fcss%na = ss%na
    fcss%cutoff = fc%cutoff
    sqrc = rc**2 + 10*lo_sqtol
    allocate (fcss%atom(fcss%na))
    do a1 = 1, fcss%na
        uca = ss%info%index_in_unitcell(a1)
        fcss%atom(a1)%n = fc%atom(uca)%n
        allocate (fcss%atom(a1)%quartet(fcss%atom(a1)%n))

        l = 0
        do i = 1, dt%particle(a1)%n
        do j = 1, dt%particle(a1)%n
        do k = 1, dt%particle(a1)%n
            if (lo_sqnorm(dt%particle(a1)%v(:, i) - dt%particle(a1)%v(:, j)) .lt. sqrc) then
            if (lo_sqnorm(dt%particle(a1)%v(:, j) - dt%particle(a1)%v(:, k)) .lt. sqrc) then
            if (lo_sqnorm(dt%particle(a1)%v(:, i) - dt%particle(a1)%v(:, k)) .lt. sqrc) then
                l = l + 1
                fcss%atom(a1)%quartet(l)%i1 = a1
                fcss%atom(a1)%quartet(l)%i2 = dt%particle(a1)%ind(i)
                fcss%atom(a1)%quartet(l)%i3 = dt%particle(a1)%ind(j)
                fcss%atom(a1)%quartet(l)%i4 = dt%particle(a1)%ind(k)
                fcss%atom(a1)%quartet(l)%lv1 = 0.0_flyt
                fcss%atom(a1)%quartet(l)%lv2 = dt%particle(a1)%lv(:, i)
                fcss%atom(a1)%quartet(l)%lv3 = dt%particle(a1)%lv(:, j)
                fcss%atom(a1)%quartet(l)%lv4 = dt%particle(a1)%lv(:, k)
                fcss%atom(a1)%quartet(l)%rv1 = 0.0_flyt
                fcss%atom(a1)%quartet(l)%rv2 = dt%particle(a1)%v(:, i)
                fcss%atom(a1)%quartet(l)%rv3 = dt%particle(a1)%v(:, j)
                fcss%atom(a1)%quartet(l)%rv4 = dt%particle(a1)%v(:, k)
                fcss%atom(a1)%quartet(l)%m = 0.0_flyt
                fcss%atom(a1)%quartet(l)%mwm = 0.0_flyt
                kk = 0
                do ii = 1, fc%atom(uca)%n
                    if (lo_sqnorm(fc%atom(uca)%quartet(ii)%rv2 - fcss%atom(a1)%quartet(l)%rv2) .gt. lo_sqtol) cycle
                    if (lo_sqnorm(fc%atom(uca)%quartet(ii)%rv3 - fcss%atom(a1)%quartet(l)%rv3) .gt. lo_sqtol) cycle
                    if (lo_sqnorm(fc%atom(uca)%quartet(ii)%rv4 - fcss%atom(a1)%quartet(l)%rv4) .gt. lo_sqtol) cycle
                    kk = ii
                    exit
                end do
                if (kk .eq. 0) then
                    call lo_stop_gracefully(['Could not locate quartet, should be impossible.'], lo_exitcode_symmetry, __FILE__, __LINE__)
                else
                    fcss%atom(a1)%quartet(l)%m = fc%atom(uca)%quartet(kk)%m
                    fcss%atom(a1)%quartet(l)%mwm = fc%atom(uca)%quartet(kk)%mwm
                end if
            end if
            end if
            end if
        end do
        end do
        end do

        if (l .ne. fc%atom(uca)%n) then
            call lo_stop_gracefully(['Inconsistent number of quartets of atom, should be impossible.'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
    end do
end subroutine

!> read the forceconstant from file
subroutine readfromfile(fc, p, fn)
    !
    class(lo_forceconstant_fourthorder), intent(out) :: fc
    type(lo_crystalstructure), intent(in) :: p
    character(len=*), intent(in) :: fn
    !
    integer :: u, a1, a2, a3, a4, i, ii, jj, kk
    real(flyt), dimension(3) :: v0, v1, v2, v3, v4, lv1, lv2, lv3, lv4, ucv1, ucv2, ucv3, ucv4

    ! get stuff from file
    u = open_file('in', trim(fn))
    read (u, *) fc%na
    read (u, *) fc%cutoff
    lo_allocate(fc%atom(fc%na))
    do a1 = 1, fc%na
        read (u, *) fc%atom(a1)%n
        lo_allocate(fc%atom(a1)%quartet(fc%atom(a1)%n))
        do i = 1, fc%atom(a1)%n
            read (u, *) fc%atom(a1)%quartet(i)%i1
            read (u, *) fc%atom(a1)%quartet(i)%i2
            read (u, *) fc%atom(a1)%quartet(i)%i3
            read (u, *) fc%atom(a1)%quartet(i)%i4
            read (u, *) fc%atom(a1)%quartet(i)%lv1
            read (u, *) fc%atom(a1)%quartet(i)%lv2
            read (u, *) fc%atom(a1)%quartet(i)%lv3
            read (u, *) fc%atom(a1)%quartet(i)%lv4
            do ii = 1, 3
            do jj = 1, 3
            do kk = 1, 3
                read (u, *) v0
                fc%atom(a1)%quartet(i)%m(ii, jj, kk, :) = v0*lo_forceconstant_4th_eVA_to_HartreeBohr
            end do
            end do
            end do
            a2 = fc%atom(a1)%quartet(i)%i2
            a3 = fc%atom(a1)%quartet(i)%i3
            a4 = fc%atom(a1)%quartet(i)%i4
            fc%atom(a1)%quartet(i)%mwm = fc%atom(a1)%quartet(i)%m*p%invsqrtmass(a1)*p%invsqrtmass(a2)*p%invsqrtmass(a3)*p%invsqrtmass(a4)
        end do
    end do
    close (u)

    ! convert stuff to Cartesian coordinates, and find the actual cutoff
    fc%cutoff = 0.0_flyt
    do a1 = 1, fc%na
        do i = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%quartet(i)%i2
            a3 = fc%atom(a1)%quartet(i)%i3
            a4 = fc%atom(a1)%quartet(i)%i4
            ! fix all the vectors
            ucv1 = p%r(:, a1)
            ucv2 = p%r(:, a2)
            ucv3 = p%r(:, a3)
            ucv4 = p%r(:, a4)
            lv1 = fc%atom(a1)%quartet(i)%lv1
            lv2 = fc%atom(a1)%quartet(i)%lv2
            lv3 = fc%atom(a1)%quartet(i)%lv3
            lv4 = fc%atom(a1)%quartet(i)%lv4
            v1 = p%fractional_to_cartesian(lv1 + ucv1)
            v2 = p%fractional_to_cartesian(lv2 + ucv2)
            v3 = p%fractional_to_cartesian(lv3 + ucv3)
            v4 = p%fractional_to_cartesian(lv4 + ucv4)
            !
            fc%atom(a1)%quartet(i)%lv1 = lo_chop(p%fractional_to_cartesian(lv1), lo_sqtol)
            fc%atom(a1)%quartet(i)%lv2 = lo_chop(p%fractional_to_cartesian(lv2), lo_sqtol)
            fc%atom(a1)%quartet(i)%lv3 = lo_chop(p%fractional_to_cartesian(lv3), lo_sqtol)
            fc%atom(a1)%quartet(i)%lv4 = lo_chop(p%fractional_to_cartesian(lv4), lo_sqtol)
            fc%atom(a1)%quartet(i)%rv1 = 0.0_flyt
            fc%atom(a1)%quartet(i)%rv2 = lo_chop(v2 - v1, lo_sqtol)
            fc%atom(a1)%quartet(i)%rv3 = lo_chop(v3 - v1, lo_sqtol)
            fc%atom(a1)%quartet(i)%rv4 = lo_chop(v4 - v1, lo_sqtol)
            !
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%quartet(i)%rv2))
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%quartet(i)%rv3))
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%quartet(i)%rv4))
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%quartet(i)%rv2 - fc%atom(a1)%quartet(i)%rv3))
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%quartet(i)%rv2 - fc%atom(a1)%quartet(i)%rv4))
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%quartet(i)%rv3 - fc%atom(a1)%quartet(i)%rv4))
        end do
    end do
end subroutine

!> write it to file
subroutine writetofile(fc, p, fn)
    !> forceconstant
    class(lo_forceconstant_fourthorder), intent(in) :: fc
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: p
    !> filename
    character(len=*), intent(in) :: fn

    integer :: u, a1, a2, a3, a4, i, ii, jj, kk
    real(flyt) :: f0
    real(flyt), dimension(3) :: v0, v1, v2, v3, v4

    ! find the cutoff
    f0 = 0.0_flyt
    do a1 = 1, fc%na
        do i = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%quartet(i)%i2
            a3 = fc%atom(a1)%quartet(i)%i3
            a4 = fc%atom(a1)%quartet(i)%i3
            v2 = p%r(:, a2) - p%r(:, a1) + fc%atom(a1)%quartet(i)%lv2
            v3 = p%r(:, a3) - p%r(:, a1) + fc%atom(a1)%quartet(i)%lv3
            v4 = p%r(:, a4) - p%r(:, a1) + fc%atom(a1)%quartet(i)%lv4
            f0 = max(f0, norm2(v2))
            f0 = max(f0, norm2(v3))
            f0 = max(f0, norm2(v4))
            f0 = max(f0, norm2(v2 - v3))
            f0 = max(f0, norm2(v2 - v4))
            f0 = max(f0, norm2(v3 - v4))
        end do
    end do

    ! Dump it
    u = open_file('out', trim(fn))
    write (u, *) fc%na
    write (u, *) (f0 + 1E-5_flyt)*lo_bohr_to_A
    do a1 = 1, fc%na
        write (u, *) fc%atom(a1)%n
        do i = 1, fc%atom(a1)%n
            write (u, *) fc%atom(a1)%quartet(i)%i1
            write (u, *) fc%atom(a1)%quartet(i)%i2
            write (u, *) fc%atom(a1)%quartet(i)%i3
            write (u, *) fc%atom(a1)%quartet(i)%i4
            v1 = anint(matmul(p%inv_latticevectors, fc%atom(a1)%quartet(i)%lv1))*1.0_flyt
            v2 = anint(matmul(p%inv_latticevectors, fc%atom(a1)%quartet(i)%lv2))*1.0_flyt
            v3 = anint(matmul(p%inv_latticevectors, fc%atom(a1)%quartet(i)%lv3))*1.0_flyt
            v4 = anint(matmul(p%inv_latticevectors, fc%atom(a1)%quartet(i)%lv4))*1.0_flyt
            write (u, *) v1
            write (u, *) v2
            write (u, *) v3
            write (u, *) v4
            do ii = 1, 3
            do jj = 1, 3
            do kk = 1, 3
                v0 = fc%atom(a1)%quartet(i)%m(ii, jj, kk, :)
                write (u, *) v0*lo_forceconstant_4th_HartreeBohr_to_eVA
            end do
            end do
            end do
        end do
    end do
    close (u)
end subroutine

!> four-phonon scattering amplitude
complex(flyt) function scatteringamplitude(fc, omega, egv, q2, q3, q4)
    !> forceconstant
    class(lo_forceconstant_fourthorder), intent(in) :: fc
    !> four frequencies
    real(flyt), dimension(4), intent(in) :: omega
    !> four eigenvectors
    complex(flyt), dimension(fc%na*3, 4), intent(in) :: egv
    !> three q-points
    real(flyt), dimension(3), intent(in) :: q2, q3, q4
    !
    type(dumegv), dimension(:, :, :, :), allocatable :: evp
    complex(flyt) :: expiqr, c0
    integer :: atom1, atom2, atom3, atom4, i, j, k, l, k1, k2, k3, k4, t
    real(flyt), dimension(3) :: rv2, rv3, rv4
    real(flyt) :: iqr, omegaprod

    ! Eigenvector product
    allocate (evp(fc%na, fc%na, fc%na, fc%na))
    do atom1 = 1, fc%na
    do atom2 = 1, fc%na
    do atom3 = 1, fc%na
    do atom4 = 1, fc%na
        do i = 1, 3
        do j = 1, 3
        do k = 1, 3
        do l = 1, 3
            k1 = (atom1 - 1)*3 + i
            k2 = (atom2 - 1)*3 + j
            k3 = (atom3 - 1)*3 + k
            k4 = (atom4 - 1)*3 + l
            evp(atom1, atom2, atom3, atom4)%m(i, j, k, l) = egv(k1, 1)*egv(k2, 2)*egv(k3, 3)*egv(k4, 4)
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do
    ! Frequency product
    omegaprod = omega(1)*omega(2)*omega(3)*omega(4)
    omegaprod = sqrt(omegaprod)
    ! The actual matrix element
    c0 = 0.0_flyt
    do atom1 = 1, fc%na
    do t = 1, fc%atom(atom1)%n
        !
        atom2 = fc%atom(atom1)%quartet(t)%i2
        atom3 = fc%atom(atom1)%quartet(t)%i3
        atom4 = fc%atom(atom1)%quartet(t)%i4

        rv2 = fc%atom(atom1)%quartet(t)%lv2
        rv3 = fc%atom(atom1)%quartet(t)%lv3
        rv4 = fc%atom(atom1)%quartet(t)%lv4

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3) + dot_product(q4, rv4)
        expiqr = cmplx(cos(iqr), sin(iqr), flyt)
        c0 = c0 + sum(fc%atom(atom1)%quartet(t)%mwm*evp(atom1, atom2, atom3, atom4)%m)*expiqr
    end do
    end do
    scatteringamplitude = c0/omegaprod
end function

end module
