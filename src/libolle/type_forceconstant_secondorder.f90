module type_forceconstant_secondorder
!!
!! Container for the second order forceconstants. Calculates dynamical matrices, frequencies, eigenvectors,
!! but those are low-level routines and there are better ones to use (in phonon_dispersions) that are more
!! user-friendly. Here it's a bunch of low-level stuff and io, nothing interesting.
!!
use konstanter, only: r8, i8, lo_iou, lo_tol, lo_sqtol, lo_pi, lo_twopi, lo_huge, lo_hugeint, lo_exitcode_param, lo_exitcode_symmetry, &
                      lo_freqtol, lo_kb_hartree, lo_frequency_Hartree_to_THz
use gottochblandat, only: tochar, walltime, lo_clean_fractional_coordinates, lo_chop, lo_sqnorm
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use lo_randomnumbers, only: lo_mersennetwister
use lo_longrange_electrostatics, only: lo_ewald_parameters
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint
implicit none

private
public :: lo_forceconstant_secondorder

!> information for a single pair
type lo_fc2_pair
    !> index in the unit cell to atom 1 and atom 2
    integer :: i1 = -lo_hugeint, i2 = -lo_hugeint
    !> lattice vectors positioning the unit cell
    real(r8), dimension(3) :: lv1 = lo_huge, lv2 = lo_huge
    !> vector between atom 1 and 2
    real(r8), dimension(3) :: r = lo_huge
    !> forceconstant
    real(r8), dimension(3, 3) :: m = lo_huge
    !> dynamically adjust forceconstant
    real(r8), dimension(3, 3) :: adjusted_m = lo_huge
    !> weight @todo Check if I actually use this anymore
    real(r8) :: weight = lo_huge
    !> which irreducible shell does it belong to?
    integer :: irreducible_shell = -lo_hugeint
    !> what operation takes the shell to here
    integer :: irreducible_operation = -lo_hugeint
end type

!> list of pairs per atom
type lo_fc2_atom
    !> how many pairs
    integer :: n = -lo_hugeint
    !> information about each pair
    type(lo_fc2_pair), dimension(:), allocatable :: pair
contains
    procedure :: size_in_mem => fc2_atom_size_in_mem
end type

!> coordination shell
type lo_fc2_shell
    !> how many pairs in this shell?
    integer :: n = -lo_hugeint
    !> for each pair, what is the pair vector
    real(r8), dimension(:, :), allocatable :: vec
    !> for each pair, which the atom at the origin
    integer, dimension(:), allocatable :: atind
    !> for each pair, which the index to the pair of that atom
    integer, dimension(:), allocatable :: pairind
    !> for each pair, which operation takes the prototype to this
    integer, dimension(:), allocatable :: opind
    !> radius of this shell
    real(r8) :: rad = lo_huge
    !> norm of the forceconstant for this shell
    real(r8) :: norm = lo_huge
    !> average weight of this shell
    real(r8) :: weight = lo_huge
    !> coefficient matrix
    real(r8), dimension(:, :), allocatable :: coeffM
    !> matrix that keeps this IFC invariant
    real(r8), dimension(9, 9) :: invarM = lo_huge
    !> indices of the IFCs, globally
    integer, dimension(:), allocatable :: ifcind
    !> how many irreducible forceconstants
    integer :: nfc = -lo_hugeint
contains
    procedure :: size_in_mem => fc2_shell_size_in_mem
end type

!> all the information needed to add polar corrections
type lo_forceconstant_secondorder_polarstuff
    !> what kind of correction is intended
    integer :: correctiontype = -lo_hugeint
    !> dielectric tensor
    real(r8), dimension(3, 3) :: eps = lo_huge
    !> Born effective charges
    real(r8), dimension(:, :, :), allocatable :: born_effective_charges
    !> On-site polar corrections
    real(r8), dimension(:, :, :), allocatable :: born_onsite_correction
    !> How many irreducible Born effective charges
    integer :: nx_Z = -lo_hugeint
    !> Irreducible representation of the Born effective charges
    real(r8), dimension(:), allocatable :: x_Z
    !> Coefficient matrix that generate Z from the irreducible
    real(r8), dimension(:, :), allocatable :: coeff_Z
    !> Largest q where polaritons are relevant
    real(r8) :: polariton_qmax = -lo_huge
contains
    procedure :: size_in_mem => fc2_polar_size_in_mem
end type

!> Secondorder forceconstant
type lo_forceconstant_secondorder
    !> Number of atoms in the cell
    integer :: na = -lo_hugeint
    !> Length of the longest pair + a tiny margin
    real(r8) :: cutoff = lo_huge
    !> Is it a polar material?
    logical :: polar = .false.
    !> Born charges, dielectric tensor and related things
    type(lo_forceconstant_secondorder_polarstuff) :: loto
    !> list of atoms
    type(lo_fc2_atom), allocatable, dimension(:) :: atom
    !> settings for ewald summations
    type(lo_ewald_parameters), private :: ew

    !> elastic constants, could be fun or something
    real(r8), dimension(6, 6) :: elastic_constants_voigt = lo_huge
    real(r8), dimension(3, 3, 3, 3) :: elastic_constants_tensor = lo_huge

    ! Can be useful with the commensurate modes:

    !> commensurate eigenvectors
    real(r8), dimension(:, :), allocatable :: eigenvectors
    !> commensurate frequencies
    real(r8), dimension(:), allocatable :: omega
    !> commensurate normal mode amplitudes
    real(r8), dimension(:), allocatable :: amplitudes

    ! It's also quite useful to have all the coordination shell info at hand, for some reason.
    integer :: npairshells = -lo_hugeint, npairop = -lo_hugeint, nifc = -lo_hugeint, nconstraints = -lo_hugeint
    type(lo_fc2_shell), dimension(:), allocatable :: pairshell
    real(r8), dimension(:, :, :), allocatable :: pairop
    real(r8), dimension(:, :, :), allocatable :: pairop3
    real(r8), dimension(:, :), allocatable :: linear_constraints

contains
    !> calculate the potential energy for a set of displacements
    procedure :: potential_energy
    !> the forces due to the given displacements
    procedure :: forces
    !> read from file
    procedure :: readfromfile
    !> write to file
    procedure :: writetofile
    !> set the sum to zero, numerically enforce the acoustic sum rules
    procedure :: setsumtozero
    !> remap it to another force constant
    procedure :: remap
    !> calculate the elastic constants from the force constant
    procedure :: get_elastic_constants
    !> initialise a cell
    procedure :: initialize_cell
    !> create a fake forceconstant from a debye temperature
    procedure :: fake_forceconstant
    !> dump the forceconstant in Abinit anaddb format
    procedure :: write_to_anaddb
    !> Set parameters for Ewald summation and enforce hermiticity of Born effective charges
    procedure :: set_ewald_and_enforce_borncharge_hermiticity
    !> Return a dynamical matrix
    procedure :: dynamicalmatrix
    !> Calculate the longrange dynamical matrix
    procedure :: longrange_dynamical_matrix
    !> Calculate the non-analytical longrange dynamical matrix
    procedure :: nonanalytical_dynamical_matrix
    !> Calculate frequencies
    procedure :: frequencies_eigenvectors_groupvelocities
    !> Return the longrange forceconstant for a supercell (dynamical matrix at Gamma, without masses)
    procedure :: supercell_longrange_dynamical_matrix_at_gamma
    !> Get the commensurate modes for a supercell
    procedure :: commensurate_modes
    !> Size in memory, in bytes
    procedure :: size_in_mem => fc2_size_in_mem
end type

! Interfaces to type_forceconstant_secondorder_aux
interface
    module subroutine commensurate_modes(fcss, ss, uc, fc, mw)
        class(lo_forceconstant_secondorder), intent(inout) :: fcss
        type(lo_crystalstructure), intent(in) :: ss
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_forceconstant_secondorder), intent(in) :: fc
        type(lo_mpi_helper), intent(inout) :: mw
    end subroutine
    module subroutine fake_forceconstant(fc, uc, ss, debye_temperature, maximum_frequency, verbosity)
        class(lo_forceconstant_secondorder), intent(out) :: fc
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), intent(in), optional :: debye_temperature
        real(r8), intent(in), optional :: maximum_frequency
        integer, intent(in), optional :: verbosity
    end subroutine
    module subroutine initialize_cell(fcss, ss, uc, fc, temperature, quantum, exact, closest_distance, mw, nosync, imode, invert, tw)
        class(lo_forceconstant_secondorder), intent(inout) :: fcss
        type(lo_crystalstructure), intent(inout) :: ss
        type(lo_crystalstructure), intent(in) :: uc
        type(lo_forceconstant_secondorder), intent(in) :: fc
        real(r8), intent(in) :: temperature
        logical, intent(in) :: quantum
        logical, intent(in) :: exact
        real(r8), intent(in) :: closest_distance
        type(lo_mpi_helper), intent(inout) :: mw
        logical, intent(in), optional :: nosync
        integer, intent(in), optional :: imode
        logical, intent(in), optional :: invert
        type(lo_mersennetwister), intent(inout), optional :: tw
    end subroutine
    module subroutine setsumtozero(fc)
        class(lo_forceconstant_secondorder), intent(inout) :: fc
    end subroutine
end interface
! Interfaces to type_forceconstant_secondorder_dynamicalmatrix
interface
    module subroutine dynamicalmatrix(fc, p, qpoint, dynamical_matrix, mem, dynamical_matrix_gradient, dynamical_matrix_hessian, qdirection, skipnonanalytical)
        class(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_crystalstructure), intent(in) :: p
        class(lo_qpoint), intent(in) :: qpoint
        complex(r8), dimension(:, :), intent(out) :: dynamical_matrix
        type(lo_mem_helper), intent(inout) :: mem
        complex(r8), dimension(:, :, :), intent(out), optional :: dynamical_matrix_gradient
        complex(r8), dimension(:, :, :), intent(out), optional :: dynamical_matrix_hessian
        real(r8), dimension(3), intent(in), optional :: qdirection
        logical, intent(in), optional :: skipnonanalytical
    end subroutine
    module subroutine frequencies_eigenvectors_groupvelocities(fc, dynamical_matrix, omega, mem, dynamical_matrix_gradient, dynamical_matrix_hessian, eigenvectors, groupvelocities, grouphessian, qpoint)
        class(lo_forceconstant_secondorder), intent(in) :: fc
        complex(r8), dimension(:, :), intent(in) :: dynamical_matrix
        real(r8), dimension(:), intent(out) :: omega
        type(lo_mem_helper), intent(inout) :: mem
        complex(r8), dimension(:, :, :), intent(in), optional :: dynamical_matrix_gradient
        complex(r8), dimension(:, :, :), intent(in), optional :: dynamical_matrix_hessian
        complex(r8), dimension(:, :), intent(out), optional :: eigenvectors
        real(r8), dimension(:, :), intent(out), optional :: groupvelocities
        real(r8), dimension(:, :, :), intent(out), optional :: grouphessian
        class(lo_qpoint), intent(in), optional :: qpoint
    end subroutine
end interface
! Interfaces to type_forceconstant_secondorder_loto
interface
    module subroutine longrange_dynamical_matrix(fc, D, p, q, Dx, Dy, Dz)
        class(lo_forceconstant_secondorder), intent(in) :: fc
        complex(r8), dimension(:, :, :, :), intent(out) :: D
        type(lo_crystalstructure), intent(in) :: p
        real(r8), dimension(3), intent(in) :: q
        complex(r8), dimension(:, :, :, :), intent(out), optional :: Dx, Dy, Dz
    end subroutine
    module subroutine nonanalytical_dynamical_matrix(fc, p, qdir, D)
        class(lo_forceconstant_secondorder), intent(in) :: fc
        type(lo_crystalstructure), intent(in) :: p
        real(r8), dimension(3), intent(in) :: qdir
        real(r8), dimension(:, :, :, :), intent(out) :: D
    end subroutine
    module subroutine supercell_longrange_dynamical_matrix_at_gamma(fc, ss, dynmat, thres)
        class(lo_forceconstant_secondorder), intent(in) :: fc
        type(lo_crystalstructure), intent(in) :: ss
        real(r8), dimension(:, :, :, :), intent(out) :: dynmat
        real(r8), intent(in) :: thres
    end subroutine
    module subroutine set_ewald_and_enforce_borncharge_hermiticity(fc, p, mem, verbosity, fixlambda)
        class(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_crystalstructure), intent(in) :: p
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
        real(r8), intent(in), optional :: fixlambda
    end subroutine
end interface
! Interfaces to type_forceconstant_secondorder_io
interface
    module subroutine writetofile(fc, p, fn)
        class(lo_forceconstant_secondorder), intent(in) :: fc
        type(lo_crystalstructure), intent(in) :: p
        character(len=*), intent(in) :: fn
    end subroutine
    module subroutine readfromfile(fc, p, fn, mem, verbosity)
        class(lo_forceconstant_secondorder), intent(out) :: fc
        type(lo_crystalstructure), intent(inout) :: p
        character(len=*), intent(in) :: fn
        type(lo_mem_helper), intent(inout) :: mem
        integer, intent(in) :: verbosity
    end subroutine
    module subroutine write_to_anaddb(fc, uc, qgrid, mw, mem)
        class(lo_forceconstant_secondorder), intent(inout) :: fc
        type(lo_crystalstructure), intent(inout) :: uc
        integer, dimension(3), intent(in) :: qgrid
        type(lo_mpi_helper), intent(inout) :: mw
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
end interface

contains

!> The potential energy for a harmonic force constant from a set of displecements
function potential_energy(fc, u) result(energy)
    !> forceconstant
    class(lo_forceconstant_secondorder), intent(in) :: fc
    !> displacements
    real(r8), dimension(:, :), intent(in) :: u
    !> energy
    real(r8) :: energy
    !
    integer :: a1, a2, k
    real(r8), dimension(3) :: v
    real(r8) :: ep

    ! sanity check
    if (fc%na .ne. size(u, 2)) then
        write (*, *) 'Displacements and forceconstants do not seem to match, perhaps they'
        write (*, *) 'are defined for different cells?'
        stop
    end if

    ep = 0.0_r8
    do a1 = 1, fc%na
        do k = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%pair(k)%i2
            v = matmul(fc%atom(a1)%pair(k)%m(:, :), u(:, a2))
            ep = ep + 0.5_r8*dot_product(u(:, a1), v)
        end do
    end do
    energy = ep
end function

!> The force from a harmonic force constant from a set of displecements
subroutine forces(fc, u, f)
    !> forceconstant
    class(lo_forceconstant_secondorder), intent(in) :: fc
    !> displacements
    real(r8), dimension(:, :), intent(in) :: u
    !> forces
    real(r8), dimension(:, :), intent(out) :: f
    !
    integer :: a1, a2, k
    real(r8), dimension(3) :: v

    f = 0.0_r8
    do a1 = 1, fc%na
        do k = 1, fc%atom(a1)%n
            a2 = fc%atom(a1)%pair(k)%i2
            v = matmul(fc%atom(a1)%pair(k)%m(:, :), u(:, a2))
            f(:, a1) = f(:, a1) - v
        end do
    end do
end subroutine

!> Map the force constant to a larger cell
subroutine remap(fc, uc, ss, fcss, ind)
    !> unitcell forceconstant
    class(lo_forceconstant_secondorder), intent(in) :: fc
    !> unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> supercell forceconstant
    type(lo_forceconstant_secondorder), intent(out) :: fcss
    !> mapping, in case I need to do it many times
    integer, dimension(:, :, :), intent(out), optional :: ind

    real(r8), dimension(3) :: r0, r1, r2
    integer :: uca, a1, a2, i

    ! Just a small sanity check
    if (ss%info%supercell .eqv. .false.) then
        call lo_stop_gracefully(['The supercell needs to be related to the unitcell'], lo_exitcode_param, __FILE__, __LINE__)
    end if
    if (uc%na .ne. fc%na) then
        call lo_stop_gracefully(['The number of atoms in the cell and in the forceconstant do not match'], lo_exitcode_param, __FILE__, __LINE__)
    end if
    !if ( fc%polar ) then
    !    call lo_stop_gracefully(['Revisit remapping for polar systems. Problem for future Olle.'],lo_exitcode_param,__FILE__,__LINE__)
    !endif

    ! Build empty forceconstants, set everything to nothing.
    fcss%na = -1
    fcss%cutoff = -1
    fcss%elastic_constants_voigt = -1
    fcss%elastic_constants_tensor = -1
    fcss%npairshells = -1
    fcss%npairop = -1
    fcss%nifc = -1
    fcss%nconstraints = -1
    fcss%polar = .false.

    fcss%loto%correctiontype = 0
    fcss%loto%eps = -1
    fcss%loto%nx_Z = -1
    ! Now set the real things that I can specify
    fcss%na = ss%na
    fcss%cutoff = fc%cutoff
    fcss%polar = fc%polar
    if (fcss%polar) then
        fcss%loto%eps = fc%loto%eps
        allocate (fcss%loto%born_effective_charges(3, 3, ss%na))
        fcss%loto%born_effective_charges = 0.0_r8
        fcss%loto%nx_Z = fc%loto%nx_Z
        allocate (fcss%loto%coeff_Z(fcss%na*9, fcss%loto%nx_Z))
        fcss%loto%coeff_Z = 0.0_r8
        allocate (fcss%loto%x_Z(fcss%loto%nx_Z))
        fcss%loto%x_Z = 0.0_r8
        ! copy irreducible Born charges
        fcss%loto%x_Z = fcss%loto%x_Z
        ! and copy Ewald parameter
        fcss%ew%lambda = fc%ew%lambda
    end if

    if (present(ind)) ind = 0
    allocate (fcss%atom(fcss%na))
    do a1 = 1, fcss%na
        r0 = ss%rcart(:, a1)
        uca = ss%info%index_in_unitcell(a1)
        if (fcss%polar) then
            fcss%loto%born_effective_charges(:, :, a1) = fc%loto%born_effective_charges(:, :, uca)
            fcss%loto%coeff_Z(3*(a1 - 1) + 1:3*a1, :) = fc%loto%coeff_Z(3*(uca - 1) + 1:3*uca, :)
        end if
        fcss%atom(a1)%n = fc%atom(uca)%n
        allocate (fcss%atom(a1)%pair(fcss%atom(a1)%n))
        do i = 1, fcss%atom(a1)%n
            fcss%atom(a1)%pair(i)%m = fc%atom(uca)%pair(i)%m
            fcss%atom(a1)%pair(i)%r = fc%atom(uca)%pair(i)%r
            fcss%atom(a1)%pair(i)%lv1 = 0.0_r8
            fcss%atom(a1)%pair(i)%lv2 = 0.0_r8
            fcss%atom(a1)%pair(i)%i1 = a1
            fcss%atom(a1)%pair(i)%i2 = 0
            ! Now, the only thing missing are the indices. Try and fix that.
            r1 = fcss%atom(a1)%pair(i)%r + r0
            r1 = ss%cartesian_to_fractional(r1)
            r1 = lo_clean_fractional_coordinates(r1, 1E-10_r8)
            do a2 = 1, ss%na
                if (sum(abs(ss%r(:, a2) - r1)) .lt. lo_tol) then
                    fcss%atom(a1)%pair(i)%i2 = a2
                    ! Store the mapping
                    if (present(ind)) then
                        ind(1, 1, a1) = uca
                        ind(2, 1, a1) = fcss%atom(a1)%n
                        ind(3, i, a1) = a2
                    end if
                    exit
                end if
            end do
            ! sanity check
            if (fcss%atom(a1)%pair(i)%i2 .eq. 0) then
                call lo_stop_gracefully(['Failed mapping secondorder forceconstants for atom '//tochar(a1)], lo_exitcode_symmetry, __FILE__, __LINE__)
            end if
            ! fix the a lattice vector?
            ! r = lv2 + r2 - r1
            ! lv2 = r - r2 + r1
            r2 = fcss%atom(a1)%pair(i)%r - ss%rcart(:, fcss%atom(a1)%pair(i)%i2) + ss%rcart(:, fcss%atom(a1)%pair(i)%i1)
            r2 = anint(matmul(ss%inv_latticevectors, r2))
            fcss%atom(a1)%pair(i)%lv2 = matmul(ss%latticevectors, r2)
        end do
    end do
end subroutine

!> Returns all 81 elastic constants, provided I got the Huang invariances right.
subroutine get_elastic_constants(fc, uc)
    !> the forceconstant
    class(lo_forceconstant_secondorder), intent(inout) :: fc
    !> the crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !
    integer :: i, j
    integer :: a1, pair
    integer :: al, be, gm, la
    real(r8), dimension(6, 6) :: elc2
    real(r8), dimension(3, 3, 3, 3) :: elc, bracket
    real(r8), dimension(3, 3) :: m
    real(r8), dimension(3) :: r

    bracket = 0.0_r8
    do a1 = 1, fc%na
    do pair = 1, fc%atom(a1)%n
        m = fc%atom(a1)%pair(pair)%m
        r = fc%atom(a1)%pair(pair)%r
        do la = 1, 3
        do gm = 1, 3
        do be = 1, 3
        do al = 1, 3
            bracket(al, be, gm, la) = bracket(al, be, gm, la) + m(al, be)*r(gm)*r(la)
        end do
        end do
        end do
        end do
    end do
    end do
    bracket = -bracket*0.5_r8/uc%volume

    do la = 1, 3
    do gm = 1, 3
    do be = 1, 3
    do al = 1, 3
        elc(al, be, gm, la) = bracket(al, gm, be, la) + bracket(be, gm, al, la) - bracket(al, be, gm, la)
    end do
    end do
    end do
    end do

    elc2 = 0.0_r8
    do la = 1, 3
    do gm = 1, 3
    do be = 1, 3
    do al = 1, 3
        i = contract_elastic_constant_indices(al, be)
        j = contract_elastic_constant_indices(gm, la)
        elc2(i, j) = elc(al, be, gm, la)
    end do
    end do
    end do
    end do

    fc%elastic_constants_voigt = elc2
    fc%elastic_constants_tensor = elc
contains
    !> Voigt-notation
    function contract_elastic_constant_indices(mu, nu) result(ind)
        integer, intent(in) :: mu, nu
        integer, dimension(2) :: d
        integer :: ind

        ind = 0
        if (mu < nu) then
            d = [mu, nu]
        else
            d = [nu, mu]
        end if
        if (d(1) .eq. d(2)) ind = d(1)
        if (d(1) .eq. 1 .and. d(2) .eq. 2) ind = 6
        if (d(1) .eq. 1 .and. d(2) .eq. 3) ind = 5
        if (d(1) .eq. 2 .and. d(2) .eq. 3) ind = 4
    end function
end subroutine

!> size in memory, in bytes
function fc2_atom_size_in_mem(a) result(mem)
    !> atom
    class(lo_fc2_atom), intent(in) :: a
    !> size in memory, bytes
    integer(i8) :: mem

    mem = 0
    mem = mem + storage_size(a)
    if (allocated(a%pair)) mem = mem + storage_size(a%pair)*size(a%pair)
    mem = mem/8
end function

!> size in memory, in bytes
function fc2_shell_size_in_mem(s) result(mem)
    !> atom
    class(lo_fc2_shell), intent(in) :: s
    !> size in memory, bytes
    integer(i8) :: mem

    mem = 0
    mem = mem + storage_size(s)
    if (allocated(s%vec)) mem = mem + storage_size(s%vec)*size(s%vec)
    if (allocated(s%atind)) mem = mem + storage_size(s%atind)*size(s%atind)
    if (allocated(s%pairind)) mem = mem + storage_size(s%pairind)*size(s%pairind)
    if (allocated(s%opind)) mem = mem + storage_size(s%opind)*size(s%opind)
    if (allocated(s%coeffM)) mem = mem + storage_size(s%coeffM)*size(s%coeffM)
    if (allocated(s%ifcind)) mem = mem + storage_size(s%ifcind)*size(s%ifcind)
    mem = mem/8
end function

!> size in memory, in bytes
function fc2_polar_size_in_mem(p) result(mem)
    !> atom
    class(lo_forceconstant_secondorder_polarstuff), intent(in) :: p
    !> size in memory, bytes
    integer(i8) :: mem

    mem = 0
    mem = mem + storage_size(p)
    if (allocated(p%born_effective_charges)) mem = mem + storage_size(p%born_effective_charges)*size(p%born_effective_charges)
    if (allocated(p%born_onsite_correction)) mem = mem + storage_size(p%born_onsite_correction)*size(p%born_onsite_correction)
    if (allocated(p%x_Z)) mem = mem + storage_size(p%x_Z)*size(p%x_Z)
    if (allocated(p%coeff_Z)) mem = mem + storage_size(p%coeff_Z)*size(p%coeff_Z)
    mem = mem/8
end function

!> size in memory, in bytes
function fc2_size_in_mem(f) result(mem)
    !> atom
    class(lo_forceconstant_secondorder), intent(in) :: f
    !> size in memory, bytes
    integer(i8) :: mem

    integer :: i
    mem = 0
    mem = mem + storage_size(f)
    if (allocated(f%eigenvectors)) mem = mem + storage_size(f%eigenvectors)*size(f%eigenvectors)
    if (allocated(f%omega)) mem = mem + storage_size(f%omega)*size(f%omega)
    if (allocated(f%amplitudes)) mem = mem + storage_size(f%amplitudes)*size(f%amplitudes)
    if (allocated(f%pairop)) mem = mem + storage_size(f%pairop)*size(f%pairop)
    if (allocated(f%pairop3)) mem = mem + storage_size(f%pairop3)*size(f%pairop3)
    if (allocated(f%linear_constraints)) mem = mem + storage_size(f%linear_constraints)*size(f%linear_constraints)
    mem = mem/8

    mem = mem + f%loto%size_in_mem()
    mem = mem + f%ew%size_in_mem()
    if (allocated(f%pairshell)) then
        do i = 1, size(f%pairshell)
            mem = mem + f%pairshell(i)%size_in_mem()
        end do
    end if
    if (allocated(f%atom)) then
        do i = 1, size(f%atom)
            mem = mem + f%atom(i)%size_in_mem()
        end do
    end if
end function

end module
