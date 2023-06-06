module type_forceconstant_thirdorder
use konstanter, only: r8, i8, lo_huge, lo_hugeint, lo_tol, lo_sqtol, lo_freqtol, lo_status, lo_twopi, lo_exitcode_physical, &
                      lo_bohr_to_A, lo_forceconstant_3rd_HartreeBohr_to_eVA, lo_forceconstant_3rd_eVA_to_HartreeBohr, &
                      lo_exitcode_symmetry, lo_exitcode_param
use gottochblandat, only: open_file, lo_chop, lo_sqnorm, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use type_distancetable, only: lo_distancetable
use type_qpointmesh, only: lo_qpoint
implicit none

private
public :: lo_forceconstant_thirdorder

!> all information for one triplet
type lo_fc3_triplet
    ! index in the unit cell, for atom 1 2 3 in the triplet
    integer :: i1 = -lo_hugeint, i2 = -lo_hugeint, i3 = -lo_hugeint
    ! lattice vector pointing to the unit cell of the atom
    real(r8), dimension(3) :: lv1 = lo_huge, lv2 = lo_huge, lv3 = lo_huge
    ! absolute vector between the atoms
    real(r8), dimension(3) :: rv1 = lo_huge, rv2 = lo_huge, rv3 = lo_huge
    ! the force constant matrix
    real(r8), dimension(3, 3, 3) :: m = lo_huge
    ! the force constant matrix with premultiplied masses
    real(r8), dimension(3, 3, 3) :: mwm = lo_huge
end type

!> list of triplets that contain this atom
type lo_fc3_atom
    integer :: n = -lo_hugeint
    type(lo_fc3_triplet), dimension(:), allocatable :: triplet
end type

!> third order forceconstant
type lo_forceconstant_thirdorder
    !> number of atoms
    integer :: na = -lo_hugeint
    !> cutoff radius
    real(r8) :: cutoff = -lo_huge
    !> list of triplets per atom
    type(lo_fc3_atom), allocatable, dimension(:) :: atom
    !> how much to talk
    integer :: verbosity = -lo_hugeint
contains
    !> read from file
    procedure :: readfromfile
    !> write to file
    procedure :: writetofile
    !> remap to another unitcell
    procedure :: remap
    !> the potential energy
    procedure :: potential_energy
    !> forces
    procedure :: forces
    !> three-phonon scattering amplitude
    procedure :: scatteringamplitude
    !> pre-transform to q-space
    procedure :: pretransform
    !> deallocate everything
    procedure :: unallocate
    !> mode gruneisen parameter
    procedure :: mode_gruneisen_parameter
    !> size in memory, in bytes
    procedure :: size_in_mem => fc3_size_in_mem
end type

contains

!> size in memory, in bytes
function fc3_size_in_mem(f) result(mem)
    !> forceconstant
    class(lo_forceconstant_thirdorder), intent(in) :: f
    !> size in memory, in bytes
    integer(i8) :: mem

    integer :: i, j

    mem = 0
    mem = mem + storage_size(f)
    if (allocated(f%atom)) then
        do i = 1, size(f%atom)
            mem = mem + storage_size(f%atom(i))
            if (allocated(f%atom(i)%triplet)) then
                do j = 1, size(f%atom(i)%triplet)
                    mem = mem + storage_size(f%atom(i)%triplet(j))
                end do
            end if
        end do
    end if
    mem = mem/8
end function

!> almost a destructor
subroutine unallocate(fc)
    class(lo_forceconstant_thirdorder), intent(inout) :: fc
    integer :: i
    ! destroy all information
    if (allocated(fc%atom)) then
        do i = 1, size(fc%atom, 1)
            if (allocated(fc%atom(i)%triplet)) deallocate (fc%atom(i)%triplet)
        end do
        deallocate (fc%atom)
    end if
    fc%na = -1
    fc%cutoff = -1.0_r8
end subroutine

!> get the potential energy in eV from a set of displacements
real(r8) function potential_energy(fc, u)
    !> forceconstant
    class(lo_forceconstant_thirdorder), intent(in) :: fc
    !> displacements in Cartesian coordinates
    real(r8), dimension(:, :), intent(in) :: u

    real(r8) :: e
    integer :: a1, t, i, j, k, atom2, atom3

    e = 0.0_r8
    do a1 = 1, fc%na
        do t = 1, fc%atom(a1)%n
            atom2 = fc%atom(a1)%triplet(t)%i2
            atom3 = fc%atom(a1)%triplet(t)%i3
            do i = 1, 3
            do j = 1, 3
            do k = 1, 3
                e = e + fc%atom(a1)%triplet(t)%m(i, j, k)*u(i, a1)*u(j, atom2)*u(k, atom3)/6.0_r8
            end do
            end do
            end do
        end do
    end do
    potential_energy = e
end function

!> forces from displacements
subroutine forces(fc, u, f)
    !> forceconstant
    class(lo_forceconstant_thirdorder), intent(in) :: fc
    !> displacements in Cartesian coordinates
    real(r8), dimension(:, :), intent(in) :: u
    !> forces
    real(r8), dimension(:, :), intent(out) :: f

    integer :: a1, t, i, j, k, atom2, atom3

    f = 0.0_r8
    do a1 = 1, fc%na
        do t = 1, fc%atom(a1)%n
            atom2 = fc%atom(a1)%triplet(t)%i2
            atom3 = fc%atom(a1)%triplet(t)%i3
            do i = 1, 3
            do j = 1, 3
            do k = 1, 3
                f(i, a1) = f(i, a1) - fc%atom(a1)%triplet(t)%m(i, j, k)*u(j, atom2)*u(k, atom3)*0.5_r8
            end do
            end do
            end do
        end do
    end do
end subroutine

!> Map the third order force constant to a supercell
subroutine remap(fc, uc, ss, fcss)
    class(lo_forceconstant_thirdorder), intent(in) :: fc
    type(lo_crystalstructure), intent(in) :: ss, uc
    type(lo_forceconstant_thirdorder), intent(out) :: fcss
    !
    type(lo_distancetable) :: dt
    integer :: i, j, k, l, a1, uca, ii
    real(r8) :: sqrc, rc

    if (fc%verbosity .gt. 0) write (*, *) 'Remapping thirdorder forceconstants'

    ! Just a small sanity check
    if (ss%info%supercell .eqv. .false.) then
        call lo_stop_gracefully(['The supercell needs to be related to the unitcell'], lo_exitcode_param, __FILE__, __LINE__)
    end if

    ! Make sure we have the proper cutoff.
    rc = 0.0_r8
    do a1 = 1, uc%na
    do i = 1, fc%atom(a1)%n
        rc = max(rc, norm2(fc%atom(a1)%triplet(i)%rv2))
        rc = max(rc, norm2(fc%atom(a1)%triplet(i)%rv3))
    end do
    end do
    rc = rc + lo_tol

    ! Build a supercell neighbourlist
    call dt%generate(ss%r, ss%latticevectors, rc, verbosity=0)

    ! Build the skeleton of the forceconstant
    fcss%na = ss%na
    fcss%cutoff = fc%cutoff
    sqrc = rc**2
    allocate (fcss%atom(fcss%na))
    do a1 = 1, fcss%na
        uca = ss%info%index_in_unitcell(a1)
        fcss%atom(a1)%n = fc%atom(uca)%n
        allocate (fcss%atom(a1)%triplet(fcss%atom(a1)%n))
        l = 0
        do i = 1, dt%particle(a1)%n
        do j = 1, dt%particle(a1)%n
            if (lo_sqnorm(dt%particle(a1)%v(:, i) - dt%particle(a1)%v(:, j)) .lt. sqrc) then
                l = l + 1
                fcss%atom(a1)%triplet(l)%i1 = a1
                fcss%atom(a1)%triplet(l)%i2 = dt%particle(a1)%ind(i)
                fcss%atom(a1)%triplet(l)%i3 = dt%particle(a1)%ind(j)
                fcss%atom(a1)%triplet(l)%lv1 = 0.0_r8
                fcss%atom(a1)%triplet(l)%lv2 = dt%particle(a1)%lv(:, i)
                fcss%atom(a1)%triplet(l)%lv3 = dt%particle(a1)%lv(:, j)
                fcss%atom(a1)%triplet(l)%rv1 = 0.0_r8
                fcss%atom(a1)%triplet(l)%rv2 = dt%particle(a1)%v(:, i)
                fcss%atom(a1)%triplet(l)%rv3 = dt%particle(a1)%v(:, j)
                fcss%atom(a1)%triplet(l)%m = 0.0_r8
                fcss%atom(a1)%triplet(l)%mwm = 0.0_r8
                k = 0
                do ii = 1, fc%atom(uca)%n
                    if (lo_sqnorm(fc%atom(uca)%triplet(ii)%rv2 - fcss%atom(a1)%triplet(l)%rv2) .lt. lo_sqtol) then
                    if (lo_sqnorm(fc%atom(uca)%triplet(ii)%rv3 - fcss%atom(a1)%triplet(l)%rv3) .lt. lo_sqtol) then
                        k = ii
                        exit
                    end if
                    end if
                end do
                if (k .eq. 0) then
                    call lo_stop_gracefully(['Failed remapping forceconstants, cells probably do not match'], lo_exitcode_symmetry, __FILE__, __LINE__)
                else
                    fcss%atom(a1)%triplet(l)%m = fc%atom(uca)%triplet(k)%m
                    fcss%atom(a1)%triplet(l)%mwm = fc%atom(uca)%triplet(k)%mwm
                end if
            end if
        end do
        end do

        if (l .ne. fcss%atom(a1)%n) then
            call lo_stop_gracefully(['Failed remapping forceconstants, cells probably do not match'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
    end do
end subroutine

!> read the forceconstant from file
subroutine readfromfile(fc, p, fn)
    !> forceconstant
    class(lo_forceconstant_thirdorder), intent(out) :: fc
    !> crystalstructure
    type(lo_crystalstructure), intent(in) :: p
    !> filename
    character(len=*), intent(in) :: fn

    integer :: u, a1, a2, a3, i, ii, jj
    real(r8), dimension(3) :: v1, v2, v3, lv1, lv2, lv3, ucv1, ucv2, ucv3

    ! Get stuff from file
    u = open_file('in', trim(fn))
    read (u, *) fc%na
    read (u, *) fc%cutoff
    allocate (fc%atom(fc%na))
    do a1 = 1, fc%na
        read (u, *) fc%atom(a1)%n
        allocate (fc%atom(a1)%triplet(fc%atom(a1)%n))
        do i = 1, fc%atom(a1)%n
            read (u, *) fc%atom(a1)%triplet(i)%i1
            read (u, *) fc%atom(a1)%triplet(i)%i2
            read (u, *) fc%atom(a1)%triplet(i)%i3
            read (u, *) fc%atom(a1)%triplet(i)%lv1
            read (u, *) fc%atom(a1)%triplet(i)%lv2
            read (u, *) fc%atom(a1)%triplet(i)%lv3
            do ii = 1, 3
            do jj = 1, 3
                read (u, *) v1
                fc%atom(a1)%triplet(i)%m(ii, jj, :) = v1*lo_forceconstant_3rd_eVA_to_HartreeBohr
            end do
            end do
            a2 = fc%atom(a1)%triplet(i)%i2
            a3 = fc%atom(a1)%triplet(i)%i3
            fc%atom(a1)%triplet(i)%mwm = fc%atom(a1)%triplet(i)%m*p%invsqrtmass(a1)*p%invsqrtmass(a2)*p%invsqrtmass(a3)
        end do
    end do
    close (u)

    ! Convert the coordinates from fractional to Cartesian
    fc%cutoff = 0.0_r8
    do a1 = 1, fc%na
        do i = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%triplet(i)%i2
            a3 = fc%atom(a1)%triplet(i)%i3
            ! fix all the vectors
            ucv1 = p%r(:, a1)
            ucv2 = p%r(:, a2)
            ucv3 = p%r(:, a3)
            lv1 = fc%atom(a1)%triplet(i)%lv1
            lv2 = fc%atom(a1)%triplet(i)%lv2
            lv3 = fc%atom(a1)%triplet(i)%lv3
            v1 = p%fractional_to_cartesian(lv1 + ucv1)
            v2 = p%fractional_to_cartesian(lv2 + ucv2)
            v3 = p%fractional_to_cartesian(lv3 + ucv3)
            !
            fc%atom(a1)%triplet(i)%lv1 = lo_chop(p%fractional_to_cartesian(lv1), lo_sqtol)
            fc%atom(a1)%triplet(i)%lv2 = lo_chop(p%fractional_to_cartesian(lv2), lo_sqtol)
            fc%atom(a1)%triplet(i)%lv3 = lo_chop(p%fractional_to_cartesian(lv3), lo_sqtol)
            fc%atom(a1)%triplet(i)%rv1 = 0.0_r8
            fc%atom(a1)%triplet(i)%rv2 = lo_chop(v2 - v1, lo_sqtol)
            fc%atom(a1)%triplet(i)%rv3 = lo_chop(v3 - v1, lo_sqtol)
            ! And get the actual cutoff
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%triplet(i)%rv2))
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%triplet(i)%rv3))
            fc%cutoff = max(fc%cutoff, norm2(fc%atom(a1)%triplet(i)%rv3 - fc%atom(a1)%triplet(i)%rv2))
        end do
    end do
end subroutine

!> write the forceconstant to file
subroutine writetofile(fc, p, fn)
    !> forceconstant
    class(lo_forceconstant_thirdorder), intent(in) :: fc
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> filename
    character(len=*), intent(in) :: fn
    !
    integer :: u, a1, a2, a3, i, ii, jj
    real(r8) :: f0
    real(r8), dimension(3) :: v0, v1, v2, v3

    ! find the cutoff
    f0 = 0.0_r8
    do a1 = 1, fc%na
        do i = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%triplet(i)%i2
            a3 = fc%atom(a1)%triplet(i)%i3
            v2 = p%r(:, a2) - p%r(:, a1) + fc%atom(a1)%triplet(i)%lv2
            v3 = p%r(:, a3) - p%r(:, a1) + fc%atom(a1)%triplet(i)%lv3
            f0 = max(f0, norm2(v2))
            f0 = max(f0, norm2(v3))
            f0 = max(f0, norm2(v2 - v3))
        end do
    end do

    u = open_file('out', trim(fn))
    write (u, "(1X,I10,15X,'How many atoms per unit cell')") fc%na
    write (u, "(1X,F20.15,5X,'Realspace cutoff (A)')") (f0 + 1E-5_r8)*lo_bohr_to_A
    do a1 = 1, fc%na
        write (u, "(1X,I10,15X,'How many triplets are atom ',I3,' part of')") fc%atom(a1)%n, a1
        do i = 1, fc%atom(a1)%n
            write (u, "(1X,I10,15X,' atom 1 in triplet',I3)") fc%atom(a1)%triplet(i)%i1, i
            write (u, "(1X,I10,15X,' atom 2 in triplet',I3)") fc%atom(a1)%triplet(i)%i2, i
            write (u, "(1X,I10,15X,' atom 3 in triplet',I3)") fc%atom(a1)%triplet(i)%i3, i
            v1 = matmul(p%inv_latticevectors, fc%atom(a1)%triplet(i)%lv1)
            v2 = matmul(p%inv_latticevectors, fc%atom(a1)%triplet(i)%lv2)
            v3 = matmul(p%inv_latticevectors, fc%atom(a1)%triplet(i)%lv3)
            v1 = anint(v1)*1.0_r8
            v2 = anint(v2)*1.0_r8
            v3 = anint(v3)*1.0_r8
            write (u, *) v1
            write (u, *) v2
            write (u, *) v3
            do ii = 1, 3
            do jj = 1, 3
                v0 = fc%atom(a1)%triplet(i)%m(ii, jj, :)
                write (u, *) v0*lo_forceconstant_3rd_HartreeBohr_to_eVA
            end do
            end do
        end do
    end do
    close (u)
end subroutine

!> The three-phonon matrix element
#ifdef AGRESSIVE_SANITY
complex(r8) function scatteringamplitude(fc, omega, egv, q2, q3)
#else
    pure complex(r8) function scatteringamplitude(fc, omega, egv, q2, q3)
#endif
        !> the third order force constant
        class(lo_forceconstant_thirdorder), intent(in) :: fc
        !> the frequencies in question
        real(r8), dimension(3), intent(in) :: omega
        !> the eigenvectors
        complex(r8), dimension(:, :), intent(in) :: egv
        !> the two q-vectors that matter
        real(r8), dimension(3), intent(in) :: q2, q3

        complex(r8), dimension(:, :, :, :, :, :), allocatable :: egvprod
        complex(r8) :: expiqr, c0
        real(r8), dimension(3) :: rv2, rv3
        real(r8) :: iqr, omegaprod
        integer :: atom1, atom2, atom3, i, j, k, k1, k2, k3, t

        ! Not enough room on the stack for large unitcells, have to allocate space
        allocate (egvprod(3, 3, 3, fc%na, fc%na, fc%na))
        ! Eigenvector product
        do atom3 = 1, fc%na
        do atom2 = 1, fc%na
        do atom1 = 1, fc%na
            do k = 1, 3
            do j = 1, 3
            do i = 1, 3
                k1 = (atom1 - 1)*3 + i
                k2 = (atom2 - 1)*3 + j
                k3 = (atom3 - 1)*3 + k
                egvprod(i, j, k, atom1, atom2, atom3) = egv(k1, 1)*egv(k2, 2)*egv(k3, 3)
            end do
            end do
            end do
        end do
        end do
        end do
#ifdef AGRESSIVE_SANITY
        if (minval(omega) .lt. lo_freqtol) then
            call lo_stop_gracefully(['Trying to calculate three-phonon matrix element with zero/imaginary phonons'], lo_exitcode_physical, __FILE__, __LINE__)
        end if
#endif
        ! Frequency product
        omegaprod = omega(1)*omega(2)*omega(3)
        omegaprod = sqrt(omegaprod)
        ! The actual matrix element
        c0 = 0.0_r8
        do atom1 = 1, fc%na
            do t = 1, fc%atom(atom1)%n
                atom2 = fc%atom(atom1)%triplet(t)%i2
                atom3 = fc%atom(atom1)%triplet(t)%i3

                rv2 = fc%atom(atom1)%triplet(t)%lv2
                rv3 = fc%atom(atom1)%triplet(t)%lv3

                iqr = dot_product(q2, rv2) + dot_product(q3, rv3)
                expiqr = cmplx(cos(iqr), sin(iqr), r8)
                ! Note that the masses are premultiplied into the forceconstants, much faster.
                c0 = c0 + sum(fc%atom(atom1)%triplet(t)%mwm*egvprod(:, :, :, atom1, atom2, atom3))*expiqr
            end do
        end do
        scatteringamplitude = c0/omegaprod
! fix the formatter
#ifdef AGRESSIVE_SANITY
    end function
#else
end function
#endif

!> do half the transform, just to q-space not to normalmode space, for fast evaluation of matrix elements
pure subroutine pretransform(fct, q2, q3, ptf)
    !> third order forceconstant
    class(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-vectors (not including two pi)
    real(r8), dimension(3), intent(in) :: q2, q3
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l

    real(r8), parameter :: prefactor = 1.0_r8/6.0_r8
    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3
    real(r8) :: iqr
    integer :: a1, a2, a3, ia, ib, ic, t, nb

    nb = fct%na*3
    ptf = 0.0_r8
    do a1 = 1, fct%na
    do t = 1, fct%atom(a1)%n
        a2 = fct%atom(a1)%triplet(t)%i2
        a3 = fct%atom(a1)%triplet(t)%i3

        rv2 = fct%atom(a1)%triplet(t)%lv2
        rv3 = fct%atom(a1)%triplet(t)%lv3

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do k = 1, 3
        do j = 1, 3
        do i = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            ! Now for the grand flattening scheme, consistent with the zger operations later
            l = (ia - 1)*nb*nb + (ib - 1)*nb + ic
            ptf(l) = ptf(l) + fct%atom(a1)%triplet(t)%m(i, j, k)*expiqr
        end do
        end do
        end do
    end do
    end do
    ! and add the prefactor 1/3!
    ptf = ptf*prefactor
end subroutine

!> Get the mode Gruneisen parameter from the third order force constants
subroutine mode_gruneisen_parameter(fct, p, qpoint, omega, eigenvectors, g)
    !> the third order force constants
    class(lo_forceconstant_thirdorder), intent(in) :: fct
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> the q-point
    class(lo_qpoint), intent(in) :: qpoint
    !> the frequencies
    real(r8), dimension(:), intent(in) :: omega
    !> the eigenvectors
    complex(r8), dimension(:, :), intent(in) :: eigenvectors
    !> the mode gruneisen parameters
    real(r8), dimension(:), intent(out) :: g
    !
    real(r8), dimension(3) :: qv, rv2, rv3
    real(r8), dimension(3, 3, 3) :: m
    real(r8) :: effmass, iqr
    complex(r8), dimension(:, :, :), allocatable :: egv
    complex(r8) :: c0, evprod, expiqr
    integer :: i, j, k, s, t
    integer :: nb, atom1, atom2, atom3

    nb = fct%na*3
    qv = qpoint%r*lo_twopi

    ! Sort eigenvectors to be per atom
    allocate (egv(3, fct%na, nb))
    do s = 1, nb
        k = 0
        do i = 1, fct%na
            do j = 1, 3
                k = k + 1
                ! sort by mode,atom,component
                egv(j, i, s) = eigenvectors(k, s)
            end do
        end do
    end do

    g = 0.0_r8
    ! Get the actual gruneisen parameters
    do s = 1, nb
        ! A parameter for each band
        c0 = 0.0_r8
        do atom1 = 1, fct%na
            do t = 1, fct%atom(atom1)%n
                atom2 = fct%atom(atom1)%triplet(t)%i2
                atom3 = fct%atom(atom1)%triplet(t)%i3
                rv2 = fct%atom(atom1)%triplet(t)%lv2
                rv3 = fct%atom(atom1)%triplet(t)%rv3
                effmass = p%invsqrtmass(atom1)*p%invsqrtmass(atom2)
                ! If I don't have minus here, everything breaks. No idea why.
                iqr = -dot_product(qv, rv2)
                expiqr = cmplx(cos(iqr), sin(iqr), r8)
                m = fct%atom(atom1)%triplet(t)%m
                do i = 1, 3
                do j = 1, 3
                    evprod = conjg(egv(i, atom1, s))*egv(j, atom2, s)
                    do k = 1, 3
                        c0 = c0 + m(i, j, k)*evprod*expiqr*rv3(k)*effmass
                    end do
                end do
                end do
            end do
        end do
        ! I don't want infinity at gamma
        if (omega(s) .lt. lo_freqtol) then
            c0 = 0.0_r8
        else
            c0 = c0/(6.0_r8*(omega(s)**2))
        end if
#ifdef AGRESSIVE_SANITY
        if (abs(aimag(c0)) .gt. lo_sqtol) then
            call lo_stop_gracefully(['Imaginary component of gruneisen parameters nonzero, not correct'], lo_exitcode_symmetry, __FILE__, __LINE__)
        end if
#endif
        g(s) = -real(c0)
    end do
    deallocate (egv)
end subroutine

end module
